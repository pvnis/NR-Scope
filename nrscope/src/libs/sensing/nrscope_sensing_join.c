/*
 * SPDX-License-Identifier: LicenseRef-CSSL-1.0
 */
/*
 * One snapshot per stream per slot from every grant of the slot. See
 * nrscope_sensing_join.h for what is joined and when.
 */
#include "nrscope/hdr/sensing/nrscope_sensing_join.h"

#include <complex.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "nrscope/hdr/sensing/nrscope_sensing.h"
#include "srsran/phy/utils/vector.h"

#define NSYMB 14
#define JOIN_MAX_KEYS 16
#define JOIN_MAX_SEGS 128
#define JOIN_MAX_RNTI 32
#define HIST_BINS 19 // |phase step| in 10 degree bins, the last one 180

nrscope_sensing_join_args_t nrscope_sensing_join_args = {
    .enable = true, .max_gap_prb = 0, .check = true, .max_step = 0.3f, .learn = true};

/* chi-square, 2 degrees of freedom, 99%: a junction of one channel is cut by noise once in
   a hundred */
#define JOIN_CHI2_99 9.21f
/* a fit is informative when its noise is 13 dB under the channel: a random precoder step
   is then seen for what it is */
#define JOIN_INFORMATIVE 0.05
/* the cell's verdict, from the informative junctions between UEs (decayed counts) */
#define JOIN_LEARN_MIN 50.0
#define JOIN_LEARN_DECAY 0.995
#define JOIN_LEARN_PRECODED 0.6 // under this share passing: the cell precodes per UE
#define JOIN_LEARN_COHERENT 0.85

/* Pilots one grant wrote in a row, as a run of consecutive lattice positions. A grant
   is several of them when its PRBs are not contiguous (FDRA type 0) or another grant
   already held some of its pilots. */
typedef struct {
  uint16_t rnti;
  int      lo, hi; // first and last subcarrier written
} seg_t;

typedef struct {
  bool     used;
  uint32_t aarx;
  uint16_t ports;
  uint8_t  layer, n_layers;
  float    beta;
  uint16_t mask; // rows written this slot
  uint32_t sfn, slot_idx;
  uint64_t slot_abs;
  int      n_grants;
  int      n_rnti;
  uint16_t rnti[JOIN_MAX_RNTI];
  int      overlap; // pilots another grant already held
  bool     seg_full;
  int      n_seg[NSYMB];
  seg_t    seg[NSYMB][JOIN_MAX_SEGS];
  cf_t*    H; // [NSYMB][n_sc]
  bool*    V;
} jkey_t;

struct nrscope_sensing_join_s {
  pthread_mutex_t           lock;
  uint32_t                  crb0, n_sc;
  bool                      open;
  uint32_t                  sfn, slot_idx;
  struct nrscope_sensing_s* s; // the context the grants came with
  int                       ssb_lo, ssb_hi; // the SSB's subcarriers on the grid, -1 when unknown
  jkey_t                    key[JOIN_MAX_KEYS];
  cf_t*                     oH; // flush output, [NSYMB][n_sc]
  bool*                     oV;
};

/* Totals, every worker */
static struct {
  pthread_mutex_t lock;
  uint64_t        slots, rows, kept, kept_single, kept_joined, dropped, dropped_single_ok;
  uint64_t        used_pos, best_pos, cover_pos; // lattice positions, summed over kept rows
  uint64_t        junctions, junction_cut, gaps, gaps_wide, gaps_filled, gaps_cut, overlap, seg_full, key_full,
      fallback;
  uint64_t        hist_junc[HIST_BINS];
  double          prb_deg[3], prb_db[3];
  uint64_t        prb_n[3];
  uint64_t        hist_prb[3][HIST_BINS];
  FILE*           csv;
  bool            csv_tried;
  double          learn_n, learn_pass; // informative junctions between UEs, and how many passed
  int             cell;                // 0 undecided, 1 coherent (joined), 2 precoded per UE (not joined)
  uint64_t        cell_changes;
  uint64_t        ssb_rows, ssb_cut; // rows over the SSB's PRBs, and those its block was cut from
} st = {.lock = PTHREAD_MUTEX_INITIALIZER};

static const char* cell_name(int c)
{
  return c == 1 ? "coherent across UEs" : c == 2 ? "precoded per UE" : "undecided";
}

nrscope_sensing_join_t* nrscope_sensing_join_alloc(uint32_t grid_crb0, uint32_t n_sc_grid)
{
  nrscope_sensing_join_t* j = calloc(1, sizeof(*j));
  if (j == NULL) {
    return NULL;
  }
  pthread_mutex_init(&j->lock, NULL);
  j->ssb_lo = j->ssb_hi = -1;
  j->crb0 = grid_crb0;
  j->n_sc = n_sc_grid;
  bool ok = true;
  for (int i = 0; i < JOIN_MAX_KEYS; i++) {
    j->key[i].H = srsran_vec_cf_malloc(NSYMB * n_sc_grid);
    j->key[i].V = calloc((size_t)NSYMB * n_sc_grid, sizeof(bool));
    ok          = ok && j->key[i].H != NULL && j->key[i].V != NULL;
  }
  j->oH = srsran_vec_cf_malloc(NSYMB * n_sc_grid);
  j->oV = calloc((size_t)NSYMB * n_sc_grid, sizeof(bool));
  if (!ok || j->oH == NULL || j->oV == NULL) {
    nrscope_sensing_join_free(j);
    return NULL;
  }
  return j;
}

void nrscope_sensing_join_free(nrscope_sensing_join_t* j)
{
  if (j == NULL) {
    return;
  }
  for (int i = 0; i < JOIN_MAX_KEYS; i++) {
    free(j->key[i].H);
    free(j->key[i].V);
  }
  free(j->oH);
  free(j->oV);
  pthread_mutex_destroy(&j->lock);
  free(j);
}

void nrscope_sensing_join_set_ssb(nrscope_sensing_join_t* j, int k_lo, int k_hi)
{
  if (j != NULL && k_hi >= k_lo && k_hi >= 0 && k_lo < (int)j->n_sc) {
    j->ssb_lo = k_lo < 0 ? 0 : k_lo;
    j->ssb_hi = k_hi >= (int)j->n_sc ? (int)j->n_sc - 1 : k_hi;
  }
}

void nrscope_sensing_join_begin(nrscope_sensing_join_t* j, uint32_t sfn, uint32_t slot_idx)
{
  if (j == NULL) {
    return;
  }
  pthread_mutex_lock(&j->lock);
  for (int i = 0; i < JOIN_MAX_KEYS; i++) {
    j->key[i].used = false;
  }
  j->open     = true;
  j->sfn      = sfn;
  j->slot_idx = slot_idx;
  pthread_mutex_unlock(&j->lock);
}

static void seg_close(jkey_t* k, int l, uint16_t rnti, int lo, int hi)
{
  if (lo < 0) {
    return;
  }
  if (k->n_seg[l] >= JOIN_MAX_SEGS) {
    k->seg_full = true;
    return;
  }
  k->seg[l][k->n_seg[l]++] = (seg_t){.rnti = rnti, .lo = lo, .hi = hi};
}

int nrscope_sensing_join_add(nrscope_sensing_join_t* j,
                             uint32_t                aarx,
                             uint16_t                ports,
                             int                     layer,
                             int                     n_layers,
                             float                   beta,
                             uint16_t                rnti,
                             uint16_t                dmrs_mask,
                             const cf_t*             H,
                             const bool*             valid,
                             uint32_t                n_sc_grid,
                             uint32_t                grid_crb0,
                             uint32_t                sfn,
                             uint32_t                slot_idx,
                             uint64_t                slot_abs)
{
  if (j == NULL || H == NULL || valid == NULL || dmrs_mask == 0) {
    return 0;
  }
  if (!isnormal(beta)) {
    beta = 1.0f;
  }
  const int shift = ((int)grid_crb0 - (int)j->crb0) * 12;
  pthread_mutex_lock(&j->lock);
  if (!j->open) {
    pthread_mutex_unlock(&j->lock);
    return 0;
  }
  jkey_t* k = NULL;
  for (int i = 0; i < JOIN_MAX_KEYS && k == NULL; i++) {
    jkey_t* c = &j->key[i];
    if (c->used && c->aarx == aarx && c->ports == ports && c->layer == (uint8_t)layer && fabsf(c->beta - beta) < 1e-3f) {
      k = c;
    }
  }
  for (int i = 0; i < JOIN_MAX_KEYS && k == NULL; i++) {
    jkey_t* c = &j->key[i];
    if (!c->used) {
      k           = c;
      k->used     = true;
      k->aarx     = aarx;
      k->ports    = ports;
      k->layer    = (uint8_t)layer;
      k->n_layers = (uint8_t)n_layers;
      k->beta     = beta;
      k->mask     = 0;
      k->n_grants = 0;
      k->n_rnti   = 0;
      k->overlap  = 0;
      k->seg_full = false;
      k->sfn      = sfn;
      k->slot_idx = slot_idx;
      k->slot_abs = slot_abs;
    }
  }
  if (k == NULL) {
    pthread_mutex_unlock(&j->lock);
    pthread_mutex_lock(&st.lock);
    st.key_full++;
    pthread_mutex_unlock(&st.lock);
    return 0;
  }
  if (k->slot_abs != slot_abs) {
    /* a grant of another slot reached this worker's open slot: not ours to join */
    pthread_mutex_unlock(&j->lock);
    return 0;
  }
  k->n_grants++;
  bool seen = false;
  for (int i = 0; i < k->n_rnti; i++) {
    seen = seen || k->rnti[i] == rnti;
  }
  if (!seen && k->n_rnti < JOIN_MAX_RNTI) {
    k->rnti[k->n_rnti++] = rnti;
  }

  for (int l = 0; l < NSYMB; l++) {
    if (!(dmrs_mask & (1u << l))) {
      continue;
    }
    const cf_t* hr = H + (size_t)l * n_sc_grid;
    const bool* vr = valid + (size_t)l * n_sc_grid;
    cf_t*       kh = k->H + (size_t)l * j->n_sc;
    bool*       kv = k->V + (size_t)l * j->n_sc;
    if (!(k->mask & (1u << l))) {
      memset(kv, 0, j->n_sc * sizeof(bool));
      k->n_seg[l] = 0;
      k->mask |= (uint16_t)(1u << l);
    }
    int lo = -1, hi = -1, step = 0;
    for (int kk = 0; kk < (int)n_sc_grid; kk++) {
      if (!vr[kk]) {
        continue;
      }
      const int kj = kk + shift;
      if (kj < 0 || kj >= (int)j->n_sc) {
        continue;
      }
      if (kv[kj]) { // the first grant keeps it
        k->overlap++;
        seg_close(k, l, rnti, lo, hi);
        lo = hi = -1;
        step    = 0;
        continue;
      }
      kh[kj] = hr[kk];
      kv[kj] = true;
      if (lo >= 0) {
        const int d = kj - hi;
        if (step == 0 && d <= 6) {
          step = d;
        } else if (d != step) {
          seg_close(k, l, rnti, lo, hi);
          lo   = kj;
          step = 0;
        }
        hi = kj;
      } else {
        lo = hi = kj;
      }
    }
    seg_close(k, l, rnti, lo, hi);
  }
  pthread_mutex_unlock(&j->lock);
  return 1;
}

/* Complex line a + b x through up to m pilots from k0 walking in dir (+1 / -1), on the
   row's common slope ph (radians per lattice step) removed, x = lattice position - x0.
   *var: variance of a from the fit's own residuals. Returns the pilots used. */
static int join_fit(const cf_t* H, const bool* V, int n_sc, int k0, int dir, int k_step, int m, float ph, double x0,
                    cf_t* a, cf_t* b, double* var)
{
  double x[32];
  cf_t   y[32];
  int    n = 0;
  if (m > 32) {
    m = 32;
  }
  for (int k = k0; k >= 0 && k < n_sc && n < m && V[k]; k += dir * k_step) {
    const double q = (double)k / k_step;
    x[n]           = q - x0;
    y[n]           = H[k] * cexpf(-I * ph * (float)q);
    n++;
  }
  if (n == 0) {
    return 0;
  }
  double s1 = 0, s2 = 0;
  cf_t   sy = 0, sxy = 0;
  for (int i = 0; i < n; i++) {
    s1 += x[i], s2 += x[i] * x[i], sy += y[i], sxy += (float)x[i] * y[i];
  }
  const double det = n * s2 - s1 * s1;
  if (n < 3 || det <= 0) {
    *a   = sy / (float)n;
    *b   = 0;
    *var = crealf(*a * conjf(*a)); // no noise estimate: never informative
    return n;
  }
  *a          = ((float)s2 * sy - (float)s1 * sxy) / (float)det;
  *b          = ((float)n * sxy - (float)s1 * sy) / (float)det;
  double res2 = 0;
  for (int i = 0; i < n; i++) {
    const cf_t r = y[i] - *a - *b * (float)x[i];
    res2 += crealf(r * conjf(r));
  }
  *var = res2 / (n - 2) * s2 / det;
  return n;
}

typedef struct {
  float  db, deg;
  bool   pass;        // the step is within noise and join_max_step
  bool   informative; // the fits are clean enough to tell a precoder step from noise
  cf_t   la, lb, ra, rb;
  double x0;
} jm_t;

/* The junction between the pilots ending at kl and the ones starting at kr: each side's
   line extrapolated to the middle, compared within their noise plus join_max_step of the
   channel. */
static bool join_measure(const cf_t* H, const bool* V, int n_sc, int kl, int kr, int k_step, int m, float ph,
                         jm_t* r)
{
  r->x0     = 0.5 * ((double)kl + (double)kr) / k_step;
  double va = 0, vb = 0;
  if (join_fit(H, V, n_sc, kl, -1, k_step, m, ph, r->x0, &r->la, &r->lb, &va) == 0 ||
      join_fit(H, V, n_sc, kr, +1, k_step, m, ph, r->x0, &r->ra, &r->rb, &vb) == 0) {
    return false;
  }
  const double ma = cabsf(r->la), mb = cabsf(r->ra);
  if (!(ma > 0) || !(mb > 0)) {
    return false;
  }
  r->db              = 20.0f * log10f((float)(mb / ma));
  r->deg             = cargf(r->ra * conjf(r->la)) * (float)(180.0 / M_PI);
  const cf_t   d     = r->ra - r->la;
  const double h2    = 0.5 * (ma * ma + mb * mb);
  const double kappa = nrscope_sensing_join_args.max_step;
  r->pass            = crealf(d * conjf(d)) <= JOIN_CHI2_99 * (va + vb) + kappa * kappa * h2;
  r->informative     = (va + vb) < JOIN_INFORMATIVE * h2;
  return true;
}

/* Phase slope of a row per lattice step, from the pairs of neighbouring pilots inside
   each run (never across a junction): the receiver's timing offset, mostly. */
static float row_slope(const cf_t* H, const bool* V, const seg_t* sg, int ns, int k_step)
{
  cf_t acc = 0;
  for (int i = 0; i < ns; i++) {
    for (int k = sg[i].lo + k_step; k <= sg[i].hi; k += k_step) {
      if (V[k] && V[k - k_step]) {
        acc += H[k] * conjf(H[k - k_step]);
      }
    }
  }
  return cabsf(acc) > 0 ? cargf(acc) : 0.0f;
}

int nrscope_sensing_join_test(const cf_t* H,
                              const bool* valid,
                              int         n_sc,
                              int         kl,
                              int         kr,
                              int         k_step,
                              int         m,
                              float       ph,
                              float*      db,
                              float*      deg,
                              bool*       informative)
{
  jm_t r;
  if (!join_measure(H, valid, n_sc, kl, kr, k_step, m, ph, &r)) {
    return -1;
  }
  *db          = r.db;
  *deg         = r.deg;
  *informative = r.informative;
  return r.pass ? 1 : 0;
}

static int hist_bin(float deg)
{
  int b = (int)(fabsf(deg) / 10.0f);
  return b >= HIST_BINS ? HIST_BINS - 1 : b;
}

static void csv_open(void)
{
  st.csv_tried = true;
  if (nrscope_sensing_args.dump[0] == 0) {
    return;
  }
  char path[300];
  snprintf(path, sizeof(path), "%s.join.csv", nrscope_sensing_args.dump);
  st.csv = fopen(path, "w");
  if (st.csv == NULL) {
    return;
  }
  setvbuf(st.csv, NULL, _IOFBF, 4 << 20);
  fprintf(st.csv, "sfn,slot,ports,layer,symbol,n_grants,n_rnti,cover_prb,best_single_prb,used_prb,kept,junctions,cuts,"
                  "gaps,gap_prb,filled,max_step_deg,max_step_db,overlap,segments\n");
}

/* One row: walk the grants in frequency, join what may be joined, keep the widest run.
   oh/ov receive the row (zero outside the kept run). Returns the run's lattice
   positions, 0 when the row is dropped. */
static int join_row(const nrscope_sensing_join_t* j, jkey_t* k, int l, cf_t* oh, bool* ov, bool stats_row)
{
  const int   n_sc = (int)j->n_sc;
  const cf_t* kh   = k->H + (size_t)l * n_sc;
  const bool* kv   = k->V + (size_t)l * n_sc;
  seg_t*      sg   = k->seg[l];
  const int   ns   = k->n_seg[l];
  memset(ov, 0, n_sc * sizeof(bool));
  if (ns == 0) {
    return 0;
  }
  // sort the runs by frequency (few of them)
  for (int i = 1; i < ns; i++) {
    seg_t x = sg[i];
    int   q = i - 1;
    while (q >= 0 && sg[q].lo > x.lo) {
      sg[q + 1] = sg[q];
      q--;
    }
    sg[q + 1] = x;
  }
  // the comb: the spacing of the first run with two pilots, 2 when none has
  int step = 0;
  for (int i = 0; i < ns && step == 0; i++) {
    if (sg[i].hi > sg[i].lo) {
      for (int q = sg[i].lo + 1; q <= sg[i].hi; q++) {
        if (kv[q]) {
          step = q - sg[i].lo;
          break;
        }
      }
    }
  }
  if (step == 0) {
    step = 2;
  }
  const int offset = sg[0].lo % step;
  const int n_lat  = (n_sc - offset + step - 1) / step;
  const int   m      = 12 / step; // one PRB of pilots on each side of a junction
  const float ph     = row_slope(kh, kv, sg, ns, step);

  // copy every run, then join
  int cover = 0, best_single = 0;
  for (int i = 0; i < ns; i++) {
    for (int q = sg[i].lo; q <= sg[i].hi; q += step) {
      oh[q] = kh[q];
      ov[q] = kv[q];
    }
    const int w = (sg[i].hi - sg[i].lo) / step + 1;
    cover += w;
  }
  // best single grant: the widest stretch of one RNTI's touching runs
  {
    int i = 0;
    while (i < ns) {
      int w = (sg[i].hi - sg[i].lo) / step + 1, q = i + 1;
      while (q < ns && sg[q].rnti == sg[i].rnti && sg[q].lo - sg[q - 1].hi == step) {
        w += (sg[q].hi - sg[q].lo) / step + 1;
        q++;
      }
      best_single = w > best_single ? w : best_single;
      i           = q;
    }
  }

  int   best_lo = sg[0].lo, best_hi = sg[0].hi;
  int   run_lo = sg[0].lo, run_hi = sg[0].hi;
  int   junctions = 0, cuts = 0, gaps = 0, gap_pos = 0, filled = 0;
  float max_deg = 0, max_db = 0;
  for (int i = 1; i <= ns; i++) {
    bool join = false;
    if (i < ns) {
      const seg_t* nx = &sg[i];
      const int    d  = nx->lo - run_hi;
      const bool   same_grant = nx->rnti == sg[i - 1].rnti;
      if (d > 0 && d % step == 0) {
        const int missing = d / step - 1;
        if (missing == 0 && same_grant) {
          join = true; // one grant, split where another one held pilots
        } else {
          if (missing > 0) {
            gaps++;
            gap_pos += missing;
          }
          if (!same_grant) {
            junctions++;
          }
          const bool gap_ok = missing * step <= nrscope_sensing_join_args.max_gap_prb * 12;
          jm_t       r      = {0};
          const bool have   = join_measure(oh, ov, n_sc, run_hi, nx->lo, step, m, ph, &r);
          if (have) {
            max_deg = fmaxf(max_deg, fabsf(r.deg));
            max_db  = fmaxf(max_db, fabsf(r.db));
          }
          /* between two UEs, the cell's verdict first: one that precodes every UE is never
            joined across them, whatever a junction looks like at low SNR */
          const bool cell_no = !same_grant && nrscope_sensing_join_args.learn && st.cell == 2;
          join = gap_ok && !cell_no && (!nrscope_sensing_join_args.check || (have && r.pass));
          if (stats_row) {
            pthread_mutex_lock(&st.lock);
            if (!same_grant && have) {
              st.hist_junc[hist_bin(r.deg)]++;
              if (r.informative && nrscope_sensing_join_args.learn) {
                st.learn_n    = st.learn_n * JOIN_LEARN_DECAY + 1.0;
                st.learn_pass = st.learn_pass * JOIN_LEARN_DECAY + (r.pass ? 1.0 : 0.0);
                if (st.learn_n >= JOIN_LEARN_MIN) {
                  const double f    = st.learn_pass / st.learn_n;
                  const int    cell = f < JOIN_LEARN_PRECODED ? 2 : f > JOIN_LEARN_COHERENT ? 1 : st.cell;
                  if (cell != st.cell) {
                    printf("sensing join: the cell looks %s (%.0f%% of %.0f clean junctions between UEs continuous)%s\n",
                           cell_name(cell), 100.0 * f, st.learn_n,
                           cell == 2 ? ": UEs no longer joined" : ": UEs joined");
                    st.cell = cell;
                    st.cell_changes++;
                  }
                }
              }
            }
            if (!same_grant) {
              st.junctions++;
              st.junction_cut += !join;
            }
            if (missing > 0) {
              st.gaps++;
              st.gaps_wide += !gap_ok;
              st.gaps_cut += gap_ok && !join;
              st.gaps_filled += join;
            }
            pthread_mutex_unlock(&st.lock);
          }
          if (join && missing > 0) {
            // across the gap: each side's line, blended from one edge to the other, on the slope
            for (int q = run_hi + step; q < nx->lo; q += step) {
              const float  t  = (float)(q - run_hi) / (float)(nx->lo - run_hi);
              const double x  = (double)q / step - r.x0;
              const cf_t   vl = r.la + r.lb * (float)x, vr = r.ra + r.rb * (float)x;
              oh[q]           = (vl * (1.0f - t) + vr * t) * cexpf(I * ph * (float)((double)q / step));
              ov[q]           = true;
            }
            filled += missing;
          }
          cuts += !join;
        }
      }
    }
    if (join) {
      run_hi = sg[i].hi;
      continue;
    }
    if ((run_hi - run_lo) > (best_hi - best_lo)) {
      best_lo = run_lo, best_hi = run_hi;
    }
    if (i < ns) {
      run_lo = sg[i].lo, run_hi = sg[i].hi;
    }
  }
  const int width = (best_hi - best_lo) / step + 1;
  const bool keep = width * 2 >= n_lat; // what nr_ue_sensing_slot_profile() asks too
  for (int q = 0; q < n_sc; q++) {
    if (q < best_lo || q > best_hi || !keep) {
      ov[q] = false;
    }
  }

  if (stats_row) {
    pthread_mutex_lock(&st.lock);
    st.rows++;
    if (keep) {
      st.kept++;
      st.used_pos += width;
      st.best_pos += best_single;
      st.cover_pos += cover;
      if (best_single * 2 >= n_lat) {
        st.kept_single++;
      } else {
        st.kept_joined++;
      }
    } else {
      st.dropped++;
    }
    st.overlap += k->overlap;
    st.seg_full += k->seg_full;
    if (!st.csv_tried) {
      csv_open();
    }
    if (st.csv != NULL) {
      const float prb = (float)step / 12.0f;
      fprintf(st.csv, "%u,%u,%u,%u,%d,%d,%d,%.1f,%.1f,%.1f,%d,%d,%d,%d,%.1f,%d,%.1f,%.2f,%d,", k->sfn, k->slot_idx,
              k->ports, k->layer, l, k->n_grants, k->n_rnti, cover * prb, best_single * prb, keep ? width * prb : 0.0f,
              keep, junctions, cuts, gaps, gap_pos * prb, filled, max_deg, max_db, k->overlap);
      for (int i = 0; i < ns && i < 48; i++) {
        fprintf(st.csv, "%s%x:%d-%d", i ? ";" : "", sg[i].rnti, sg[i].lo / 12, sg[i].hi / 12);
      }
      fputc('\n', st.csv);
    }
    pthread_mutex_unlock(&st.lock);
  }
  return keep ? width : 0;
}

/* Phase and amplitude steps at the PRB boundaries inside each grant, by CRB class:
   a precoder that changes per PRG of 2 or 4 PRBs steps at CRB multiples of 2 or 4
   (PRGs are aligned to CRB 0, 38.214 5.1.2.3), a wideband or no precoder nowhere. */
static void prg_stats(const nrscope_sensing_join_t* j, const jkey_t* k, int l)
{
  const int   n_sc = (int)j->n_sc;
  const cf_t* kh   = k->H + (size_t)l * n_sc;
  const bool* kv   = k->V + (size_t)l * n_sc;
  double      sdeg[3] = {0}, sdb[3] = {0};
  uint64_t    n[3] = {0}, hist[3][HIST_BINS];
  memset(hist, 0, sizeof(hist));
  for (int i = 0; i < k->n_seg[l]; i++) {
    const seg_t* s = &k->seg[l][i];
    // the comb inside this run
    int step = 0;
    for (int q = s->lo + 1; q <= s->hi && step == 0; q++) {
      step = kv[q] ? q - s->lo : 0;
    }
    if (step == 0) {
      continue;
    }
    const int   m  = 12 / step;
    const float ph = row_slope(kh, kv, s, 1, step);
    for (int kb = (s->lo / 12 + 1) * 12; kb <= s->hi; kb += 12) {
      int kl = kb - 1;
      while (kl > s->lo && !kv[kl]) {
        kl--;
      }
      int kr = kb;
      while (kr < s->hi && !kv[kr]) {
        kr++;
      }
      if (kl < kb - 12 || kr >= kb + 12) {
        continue;
      }
      jm_t r;
      if (!join_measure(kh, kv, n_sc, kl, kr, step, m, ph, &r)) {
        continue;
      }
      const float db = r.db, deg = r.deg;
      const uint32_t crb = j->crb0 + (uint32_t)(kb / 12);
      const int      c   = (crb % 4 == 0) ? 0 : (crb % 2 == 0) ? 1 : 2;
      sdeg[c] += fabsf(deg);
      sdb[c] += fabsf(db);
      n[c]++;
      hist[c][hist_bin(deg)]++;
    }
  }
  pthread_mutex_lock(&st.lock);
  for (int c = 0; c < 3; c++) {
    st.prb_deg[c] += sdeg[c];
    st.prb_db[c] += sdb[c];
    st.prb_n[c] += n[c];
    for (int b = 0; b < HIST_BINS; b++) {
      st.hist_prb[c][b] += hist[c][b];
    }
  }
  pthread_mutex_unlock(&st.lock);
}

/* DM-RS coherence of the pilots in [a, b]: pairs 4 subcarriers apart (the same cover
   weight on comb 0, neighbours once despread). Returns the pairs used. */
static int coh_range(const cf_t* H, const bool* V, int a, int b, int n_sc, double* c)
{
  cf_t   acc = 0;
  double den = 0;
  int    n   = 0;
  for (int k = a < 0 ? 0 : a; k + 4 <= b && k + 4 < n_sc; k++) {
    if (V[k] && V[k + 4]) {
      acc += H[k + 4] * conjf(H[k]);
      den += 0.5 * (crealf(H[k] * conjf(H[k])) + crealf(H[k + 4] * conjf(H[k + 4])));
      n++;
    }
  }
  *c = den > 0 ? cabsf(acc) / den : 0;
  return n;
}

/* An SSB in this symbol: its PRBs are not available to the PDSCH (38.214 5.1.4), so the
   grants carry no DM-RS there and the estimates are the SSB's own REs. Seen on air rather
   than predicted from the burst, half frame and periodicity: the pilots over the SSB's
   PRBs lose their coherence while the rest of the row keeps it. They are taken out, which
   splits the runs there. Returns true when cut. */
static bool ssb_cut(nrscope_sensing_join_t* j, jkey_t* k, int l)
{
  if (j->ssb_lo < 0 || k->n_seg[l] == 0) {
    return false;
  }
  const int n_sc = (int)j->n_sc;
  cf_t*     kh   = k->H + (size_t)l * n_sc;
  bool*     kv   = k->V + (size_t)l * n_sc;
  const int a    = (j->ssb_lo / 12) * 12, b = (j->ssb_hi / 12) * 12 + 11; // whole PRBs
  double    c_in = 0, c_lo = 0, c_hi = 0;
  const int n_in = coh_range(kh, kv, a, b, n_sc, &c_in);
  if (n_in < 8) {
    return false; // no grant over the SSB
  }
  const int    n_lo = coh_range(kh, kv, 0, a - 13, n_sc, &c_lo);
  const int    n_hi = coh_range(kh, kv, b + 13, n_sc - 1, n_sc, &c_hi);
  const double c_out = (n_lo * c_lo + n_hi * c_hi) / (n_lo + n_hi > 0 ? n_lo + n_hi : 1);
  bool         cut   = n_lo + n_hi >= 24 && c_out > 0.6 && c_in < 0.5 * c_out;
  if (k->aarx == 0) {
    pthread_mutex_lock(&st.lock);
    st.ssb_rows++;
    st.ssb_cut += cut;
    pthread_mutex_unlock(&st.lock);
  }
  if (!cut) {
    return false;
  }
  for (int q = a; q <= b && q < n_sc; q++) {
    kv[q] = false;
  }
  // split the runs at the hole
  seg_t     keep[JOIN_MAX_SEGS];
  int       n = 0;
  const int ns = k->n_seg[l];
  for (int i = 0; i < ns && n < JOIN_MAX_SEGS; i++) {
    seg_t sg = k->seg[l][i];
    if (sg.hi < a || sg.lo > b) {
      keep[n++] = sg;
      continue;
    }
    if (sg.lo < a) {
      int hi = a - 1;
      while (hi > sg.lo && !kv[hi]) {
        hi--;
      }
      if (kv[hi]) {
        keep[n++] = (seg_t){sg.rnti, sg.lo, hi};
      }
    }
    if (sg.hi > b && n < JOIN_MAX_SEGS) {
      int lo = b + 1;
      while (lo < sg.hi && !kv[lo]) {
        lo++;
      }
      if (kv[lo]) {
        keep[n++] = (seg_t){sg.rnti, lo, sg.hi};
      }
    }
  }
  memcpy(k->seg[l], keep, n * sizeof(seg_t));
  k->n_seg[l] = n;
  (void)kh;
  return true;
}

typedef void (*emit_t)(void* ctx, const jkey_t* k, uint16_t mask, const cf_t* H, const bool* V, uint32_t n_sc);

static int flush_impl(nrscope_sensing_join_t* j, emit_t emit, void* ctx)
{
  pthread_mutex_lock(&j->lock);
  j->open = false;
  // chain 0 first, then by stream, which is the order the live path pushed in
  int order[JOIN_MAX_KEYS], n = 0;
  for (int i = 0; i < JOIN_MAX_KEYS; i++) {
    if (j->key[i].used) {
      order[n++] = i;
    }
  }
  for (int i = 1; i < n; i++) {
    const int x = order[i];
    int       q = i - 1;
    while (q >= 0) {
      const jkey_t* a = &j->key[order[q]];
      const jkey_t* b = &j->key[x];
      const bool    after =
          a->aarx > b->aarx ||
          (a->aarx == b->aarx && (a->ports > b->ports || (a->ports == b->ports && a->layer > b->layer)));
      if (!after) {
        break;
      }
      order[q + 1] = order[q];
      q--;
    }
    order[q + 1] = x;
  }
  int pushed = 0;
  for (int i = 0; i < n; i++) {
    jkey_t*  k        = &j->key[order[i]];
    uint16_t out_mask = 0;
    bool     first    = true;
    for (int l = 0; l < NSYMB; l++) {
      if (!(k->mask & (1u << l))) {
        continue;
      }
      const bool stats_row = k->aarx == 0;
      ssb_cut(j, k, l);
      if (stats_row && k->layer == 0 && first) {
        prg_stats(j, k, l);
      }
      first = false;
      if (join_row(j, k, l, j->oH + (size_t)l * j->n_sc, j->oV + (size_t)l * j->n_sc, stats_row) > 0) {
        out_mask |= (uint16_t)(1u << l);
      }
    }
    if (out_mask != 0 && emit != NULL) {
      emit(ctx, k, out_mask, j->oH, j->oV, j->n_sc);
      pushed++;
    }
  }
  if (n > 0) {
    pthread_mutex_lock(&st.lock);
    st.slots++;
    pthread_mutex_unlock(&st.lock);
  }
  pthread_mutex_unlock(&j->lock);
  return pushed;
}

typedef struct {
  nrscope_sensing_join_sink_t sink;
  void*                       ctx;
} user_sink_t;

static void emit_user(void* ctx, const jkey_t* k, uint16_t mask, const cf_t* H, const bool* V, uint32_t n_sc)
{
  const user_sink_t* u = (const user_sink_t*)ctx;
  u->sink(u->ctx, k->aarx, k->ports, k->layer, mask, H, V, n_sc);
}

int nrscope_sensing_join_flush_to(nrscope_sensing_join_t* j, nrscope_sensing_join_sink_t sink, void* ctx)
{
  if (j == NULL) {
    return 0;
  }
  user_sink_t u = {sink, ctx};
  return flush_impl(j, sink != NULL ? emit_user : NULL, &u);
}

typedef struct {
  struct nrscope_sensing_s* s;
  int                       pushed;
} live_sink_t;

static void emit_live(void* ctx, const jkey_t* k, uint16_t mask, const cf_t* H, const bool* V, uint32_t n_sc)
{
  live_sink_t* c = (live_sink_t*)ctx;
  c->pushed += nrscope_sensing_push_recorded(c->s, k->aarx, k->ports, k->layer, k->n_layers, k->sfn, k->slot_idx,
                                             k->slot_abs, mask, H, V, n_sc);
}

int nrscope_sensing_join_flush(nrscope_sensing_join_t* j, struct nrscope_sensing_s* s)
{
  if (j == NULL) {
    return 0;
  }
  live_sink_t c = {.s = s != NULL ? s : j->s, .pushed = 0};
  if (c.s == NULL) { // nothing was added: no context seen yet
    flush_impl(j, NULL, NULL);
    return 0;
  }
  flush_impl(j, emit_live, &c);
  return c.pushed;
}

void nrscope_sensing_join_set_context(nrscope_sensing_join_t* j, struct nrscope_sensing_s* s)
{
  if (j != NULL && s != NULL) {
    j->s = s;
  }
}

static float hist_quantile(const uint64_t* h, double q)
{
  uint64_t tot = 0;
  for (int b = 0; b < HIST_BINS; b++) {
    tot += h[b];
  }
  if (tot == 0) {
    return NAN;
  }
  uint64_t acc = 0;
  for (int b = 0; b < HIST_BINS; b++) {
    acc += h[b];
    if ((double)acc >= q * (double)tot) {
      return (float)(b * 10 + 10);
    }
  }
  return 180.0f;
}

int nrscope_sensing_join_cell(void)
{
  pthread_mutex_lock(&st.lock);
  const int c = st.cell;
  pthread_mutex_unlock(&st.lock);
  return c;
}

void nrscope_sensing_join_reset(void)
{
  pthread_mutex_lock(&st.lock);
  FILE* f = st.csv;
  memset((char*)&st + sizeof(st.lock), 0, sizeof(st) - sizeof(st.lock));
  st.csv       = f;
  st.csv_tried = f != NULL;
  pthread_mutex_unlock(&st.lock);
}

void nrscope_sensing_join_report(void)
{
  pthread_mutex_lock(&st.lock);
  if (st.csv != NULL) {
    fclose(st.csv);
    st.csv = NULL;
  }
  if (st.rows == 0) {
    pthread_mutex_unlock(&st.lock);
    return;
  }
  printf("sensing join: %lu slots, %lu rows (stream x DM-RS symbol, chain 0)\n", (unsigned long)st.slots,
         (unsigned long)st.rows);
  printf("  kept %lu (%.1f%%): %lu wide enough with one grant, %lu only by joining; %lu dropped (under half the "
         "carrier after joining)\n",
         (unsigned long)st.kept, 100.0 * st.kept / st.rows, (unsigned long)st.kept_single,
         (unsigned long)st.kept_joined, (unsigned long)st.dropped);
  if (st.kept > 0) {
    printf("  kept rows: mean %.1f%% of the pilots used, the widest single grant held %.1f%% of them\n",
           100.0 * st.used_pos / st.cover_pos, 100.0 * st.best_pos / st.cover_pos);
  }
  printf("  junctions between UEs %lu: %lu cut; |phase step| median <%.0f deg, 90%% <%.0f deg; cell %s (%.0f%% of "
         "%.0f clean junctions continuous)\n",
         (unsigned long)st.junctions, (unsigned long)st.junction_cut, hist_quantile(st.hist_junc, 0.5),
         hist_quantile(st.hist_junc, 0.9), cell_name(st.cell), st.learn_n > 0 ? 100.0 * st.learn_pass / st.learn_n : 0.0,
         st.learn_n);
  printf("  gaps %lu: %lu wider than join_max_gap_prb=%d, %lu cut, %lu filled; %lu pilots overlapping (first kept)%s%s\n",
         (unsigned long)st.gaps, (unsigned long)st.gaps_wide, nrscope_sensing_join_args.max_gap_prb,
         (unsigned long)st.gaps_cut, (unsigned long)st.gaps_filled, (unsigned long)st.overlap,
         st.seg_full ? ", segment table full" : "", st.key_full ? ", key table full" : "");
  if (st.ssb_rows > 0) {
    printf("  rows over the SSB's PRBs %lu: its block cut from %lu (no DM-RS there: an SSB in that symbol)\n",
           (unsigned long)st.ssb_rows, (unsigned long)st.ssb_cut);
  }
  static const char* cls[3] = {"CRB%4==0", "CRB%4==2", "odd CRB"};
  printf("  steps at PRB boundaries inside one grant (precoder granularity):");
  for (int c = 0; c < 3; c++) {
    if (st.prb_n[c] > 0) {
      printf(" %s %.1f deg %.2f dB (90%% <%.0f deg, n %lu);", cls[c], st.prb_deg[c] / st.prb_n[c],
             st.prb_db[c] / st.prb_n[c], hist_quantile(st.hist_prb[c], 0.9), (unsigned long)st.prb_n[c]);
    }
  }
  printf("\n");
  pthread_mutex_unlock(&st.lock);
}
