/*
 * SPDX-License-Identifier: LicenseRef-CSSL-1.0
 */
/*
 * sensing_offline: replay a sensing.record_estimates file through the sensing pipeline.
 *
 *   sensing_offline <config.yaml> <estimates.bin> [--dump <map2d.csv>] [--inject <spec>]...
 *
 * The file holds every grant's DM-RS channel estimates as the live receiver produced
 * them (nrscope_sensing_record.h), so everything upstream of the estimates -- the air,
 * the gNB, the radio, the sniffer -- is the recorded run's, impairments included.
 * Each record is handed to nrscope_sensing_push_estimates(), the same call the live
 * receiver makes, in the order it arrived live. From there on it is the live code:
 * random drop, alignment, history, map trigger, clutter removal, Doppler, detector,
 * AoA, and the same dumps, with the settings of the yaml's sensing: block. Change the
 * pipeline or the yaml, rebuild, and rerun on identical data.
 *
 * Only the sensing: block of the yaml is read; no radio is opened. --dump overrides
 * sensing.dump so a replay does not overwrite the live run's maps.
 *
 * INJECTED TARGETS (--inject, up to MAX_INJECT of them)
 *
 * A target with a known trajectory is added to the recorded channel, so what the
 * pipeline finds can be scored against the truth on the real run's residue.
 *
 * It is added as a delayed, Doppler-shifted, attenuated copy of the measured channel
 * itself, not as a clean synthetic path:
 *
 *     H'[k] = H[k] * (1 + g * exp(-j 2 pi (k - K/2) df tau(t)) * exp(-j 2 pi R(t) / lambda) * exp(j psi_a))
 *
 * An echo is the transmitted signal arriving later, so it shares whatever happened to
 * the direct path at that instant: the gNB's power and precoder, the timing jumps,
 * the residual frequency offset, each chain's LO phase. The copy carries all of that,
 * and the alignment then treats the target exactly as it would a real one. A clean path
 * would carry none of it and would instead receive the inverse of the scene's
 * corrections from the alignment.
 *
 * The copy also brings a copy of the static scene, g down and tau later, and of the
 * noise, g down. At the levels worth testing (-20 dB and below) both are small, and
 * the static part is removed by the clutter stage like any static path.
 *
 * The trajectory is the excess path length over the direct path, the map's own range
 * axis, and its rate sets the Doppler, so range drift and Doppler always agree:
 *     R(t) = range - 2 * speed * t' + swing * sin(2 pi t' / period),  t' = t - start
 * speed is in the map's convention (f = 2 speed / lambda, positive approaching).
 * Spec, comma-separated key=value, every key optional except range:
 *     range=<m>     excess path length at start (map range axis, from the direct path)
 *     speed=<m/s>   radial speed, map convention                         (default 0)
 *     level=<dB>    g, relative to the direct path                       (default -25)
 *     swing=<m>     back-and-forth amplitude of R, with period           (default 0)
 *     period=<s>    period of the swing                                  (default 0: none)
 *     angle=<deg>   direction relative to the direct path's, half-wavelength ULA:
 *                   psi_a = 2 pi * NR_AOA_ELEMENT_SPACING * a * sin(angle)  (default 0).
 *                   The copy keeps the direct path's per-chain phases (its angle and
 *                   the chains' LO offsets, as a real echo would), so the AoA reports
 *                   asin(sin(angle of the direct path) + sin(angle)), not angle itself
 *     start=<s>, stop=<s>  when it is present, from the recording's first slot
 *
 * The truth goes to <dump>.truth.csv, one row per slot and target, for
 * scripts/sensing/score_injection.py.
 */
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include <yaml-cpp/yaml.h>

#include "nrscope/hdr/load_config.h"
#include "nrscope/hdr/sensing/nrscope_sensing.h"
#include "nrscope/hdr/sensing/nrscope_sensing_record.h"

// nr_ue_sensing.h is C with VLA prototypes, so only the one global is declared here
extern "C" int nr_sensing_tdd_period_slots;

#define NSYMB SRSRAN_NSYMB_PER_SLOT_NR
#define MAX_INJECT 4
#define C_LIGHT 299792458.0
// receive array spacing in wavelengths, NR_AOA_ELEMENT_SPACING in nr_ue_aoa.h
#define ULA_SPACING 0.5

struct inject_t {
  double range_m  = NAN;
  double speed_ms = 0.0;
  double level_db = -25.0;
  double swing_m  = 0.0;
  double period_s = 0.0;
  double angle_deg = 0.0;
  double start_s  = 0.0;
  double stop_s   = INFINITY;

  bool active(double t) const { return t >= start_s && t < stop_s; }
  double range_at(double t) const
  {
    const double u = t - start_s;
    double       r = range_m - 2.0 * speed_ms * u;
    if (period_s > 0.0)
      r += swing_m * sin(2.0 * M_PI * u / period_s);
    return r;
  }
  // map convention: speed = -dR/dt / 2
  double speed_at(double t) const
  {
    const double u = t - start_s;
    double       v = speed_ms;
    if (period_s > 0.0)
      v -= swing_m * M_PI / period_s * cos(2.0 * M_PI * u / period_s);
    return v;
  }
};

static bool parse_inject(const char* spec, inject_t* out)
{
  inject_t    t;
  std::string str(spec);
  size_t      pos = 0;
  while (pos <= str.size()) {
    size_t            end = str.find(',', pos);
    const std::string kv  = str.substr(pos, end == std::string::npos ? std::string::npos : end - pos);
    pos                   = end == std::string::npos ? str.size() + 1 : end + 1;
    const size_t eq       = kv.find('=');
    if (eq == std::string::npos) {
      fprintf(stderr, "--inject: '%s' is not key=value\n", kv.c_str());
      return false;
    }
    const std::string k = kv.substr(0, eq);
    char*             e = nullptr;
    const double      v = strtod(kv.c_str() + eq + 1, &e);
    if (e == kv.c_str() + eq + 1 || *e != 0) {
      fprintf(stderr, "--inject: '%s' has no number\n", kv.c_str());
      return false;
    }
    if (k == "range")
      t.range_m = v;
    else if (k == "speed")
      t.speed_ms = v;
    else if (k == "level")
      t.level_db = v;
    else if (k == "swing")
      t.swing_m = v;
    else if (k == "period")
      t.period_s = v;
    else if (k == "angle")
      t.angle_deg = v;
    else if (k == "start")
      t.start_s = v;
    else if (k == "stop")
      t.stop_s = v;
    else {
      fprintf(stderr, "--inject: unknown key '%s'\n", k.c_str());
      return false;
    }
  }
  if (std::isnan(t.range_m)) {
    fprintf(stderr, "--inject: range= is required\n");
    return false;
  }
  *out = t;
  return true;
}

/* Add the targets to one record's estimates, in place. t_slot is the slot's start
   from the recording's first slot, in seconds. */
static void inject_record(const std::vector<inject_t>& tg, const nrsr_rec_hdr_t& rh, const nrsr_file_hdr_t& fh,
                          double t_slot, double t_symbol, cf_t* H, const bool* valid)
{
  const uint32_t n_sc   = fh.n_sc_grid;
  const double   lambda = C_LIGHT / (double)fh.carrier_hz;
  for (int l = 0; l < (int)NSYMB; l++) {
    if (!(rh.dmrs_mask & (1u << l)))
      continue;
    const double t = t_slot + l * t_symbol;
    // per target: the gain at subcarrier 0 and the rotation from one subcarrier to the next
    std::complex<double> g0[MAX_INJECT], step[MAX_INJECT];
    int                  n = 0;
    for (const inject_t& x : tg) {
      if (!x.active(t))
        continue;
      const double R   = x.range_at(t);
      const double tau = R / C_LIGHT;
      const double ph  = -2.0 * M_PI * R / lambda + 2.0 * M_PI * ULA_SPACING * rh.aarx * sin(x.angle_deg * M_PI / 180.0)
                        - 2.0 * M_PI * (0.0 - n_sc / 2.0) * fh.scs_hz * tau;
      g0[n]   = std::polar(pow(10.0, x.level_db / 20.0), ph);
      step[n] = std::polar(1.0, -2.0 * M_PI * fh.scs_hz * tau);
      n++;
    }
    if (n == 0)
      continue;
    cf_t*       h = &H[(size_t)l * n_sc];
    const bool* v = &valid[(size_t)l * n_sc];
    for (uint32_t k = 0; k < n_sc; k++) {
      std::complex<double> m = 1.0;
      for (int i = 0; i < n; i++) {
        m += g0[i];
        g0[i] *= step[i];
      }
      if (v[k]) {
        // cf_t is GCC's __complex__ float here, not std::complex
        const std::complex<double> y = std::complex<double>(__real__ h[k], __imag__ h[k]) * m;
        __real__ h[k]                = (float)y.real();
        __imag__ h[k]                = (float)y.imag();
      }
    }
  }
}

static void usage(const char* argv0)
{
  fprintf(stderr,
          "usage: %s <config.yaml> <estimates.bin> [--dump <map2d.csv>] [--inject <spec>]...\n"
          "  spec: range=<m>[,speed=<m/s>][,level=<dB>][,swing=<m>,period=<s>][,angle=<deg>][,start=<s>][,stop=<s>]\n",
          argv0);
}

int main(int argc, char** argv)
{
  setvbuf(stdout, nullptr, _IOLBF, 0);
  if (argc < 3) {
    usage(argv[0]);
    return EXIT_FAILURE;
  }
  const std::string cfg_path = argv[1];
  const std::string rec_path = argv[2];
  std::string           dump;
  std::vector<inject_t> targets;
  for (int i = 3; i < argc; i++) {
    if (strcmp(argv[i], "--dump") == 0 && i + 1 < argc) {
      dump = argv[++i];
    } else if (strcmp(argv[i], "--inject") == 0 && i + 1 < argc) {
      inject_t t;
      if (!parse_inject(argv[++i], &t))
        return EXIT_FAILURE;
      if (targets.size() == MAX_INJECT) {
        fprintf(stderr, "--inject: at most %d targets\n", MAX_INJECT);
        return EXIT_FAILURE;
      }
      targets.push_back(t);
    } else {
      usage(argv[0]);
      return EXIT_FAILURE;
    }
  }

  load_sensing_config(YAML::LoadFile(cfg_path));
  nrscope_sensing_args.enable              = true;
  nrscope_sensing_args.record_estimates[0] = 0; // never record the replay into itself
  if (!dump.empty()) {
    snprintf(nrscope_sensing_args.dump, sizeof(nrscope_sensing_args.dump), "%s", dump.c_str());
  }

  FILE* f = fopen(rec_path.c_str(), "rb");
  if (f == nullptr) {
    perror(rec_path.c_str());
    return EXIT_FAILURE;
  }
  setvbuf(f, nullptr, _IOFBF, 8 << 20);
  nrsr_file_hdr_t fh;
  if (fread(&fh, sizeof(fh), 1, f) != 1 || memcmp(fh.magic, NRSR_MAGIC, sizeof(NRSR_MAGIC)) != 0) {
    fprintf(stderr, "%s: not a sensing estimates file\n", rec_path.c_str());
    return EXIT_FAILURE;
  }
  if (fh.version != NRSR_VERSION) {
    fprintf(stderr, "%s: format version %u, this build reads %u\n", rec_path.c_str(), fh.version, NRSR_VERSION);
    return EXIT_FAILURE;
  }
  printf("%s: %u chain(s), carrier %.3f MHz, %.2f Msps, FFT %u, %u kHz SCS, %u subcarriers\n", rec_path.c_str(),
         fh.nof_antennas, fh.carrier_hz / 1e6, fh.srate_hz / 1e6, fh.ofdm_size, fh.scs_hz / 1000, fh.n_sc_grid);
  printf("dump: %s\n", nrscope_sensing_args.dump);

  nrscope_sensing_t* s =
      nrscope_sensing_get(fh.carrier_hz, fh.srate_hz, fh.ofdm_size, fh.scs_hz, fh.nof_antennas);
  if (s == nullptr) {
    fprintf(stderr, "could not create the sensing context\n");
    return EXIT_FAILURE;
  }
  nrscope_sensing_set_offline(s, true);

  const double t_slot_s   = 1e-3 / (fh.scs_hz / 15000);
  FILE*        truth      = nullptr;
  uint64_t     truth_slot = UINT64_MAX;
  if (!targets.empty()) {
    const std::string tp = std::string(nrscope_sensing_args.dump) + ".truth.csv";
    truth                = fopen(tp.c_str(), "w");
    if (truth == nullptr) {
      perror(tp.c_str());
      return EXIT_FAILURE;
    }
    fprintf(truth, "t_s,frame,slot,slot_abs,target,active,range_m,speed_ms,level_db,angle_deg\n");
    for (size_t i = 0; i < targets.size(); i++) {
      const inject_t& x = targets[i];
      printf("inject %zu: range %.1f m, speed %+.2f m/s, %.1f dB, swing %.1f m / %.1f s, angle %.0f deg, %.1f..%.1f s\n",
             i, x.range_m, x.speed_ms, x.level_db, x.swing_m, x.period_s, x.angle_deg, x.start_s, x.stop_s);
    }
    printf("truth: %s\n", tp.c_str());
  }

  const uint32_t    n_sc     = fh.n_sc_grid;
  const size_t      n_bitmap = (n_sc + 7) / 8;
  std::vector<cf_t> H((size_t)NSYMB * n_sc);
  // std::vector<bool> is a bitset, and push_estimates needs a bool array
  bool*                valid = static_cast<bool*>(calloc((size_t)NSYMB * n_sc, sizeof(bool)));
  std::vector<uint8_t> bitmap(n_bitmap);
  std::vector<cf_t>    vals(n_sc);

  const auto t_start   = std::chrono::steady_clock::now();
  uint64_t   n_records = 0, n_pushed = 0, slot_first = 0, slot_last = 0;
  bool       truncated = false;
  for (;;) {
    nrsr_rec_hdr_t rh;
    if (fread(&rh, sizeof(rh), 1, f) != 1) {
      break;
    }
    if (rh.magic != NRSR_REC_MAGIC) {
      fprintf(stderr, "record %lu: bad magic, file corrupt; stopping here\n", (unsigned long)n_records);
      break;
    }
    memset(valid, 0, (size_t)NSYMB * n_sc * sizeof(bool));
    bool ok = true;
    for (int l = 0; l < (int)NSYMB && ok; l++) {
      if (!(rh.dmrs_mask & (1u << l))) {
        continue;
      }
      uint32_t n_valid = 0;
      ok = fread(&n_valid, sizeof(n_valid), 1, f) == 1 && n_valid <= n_sc &&
           fread(bitmap.data(), 1, n_bitmap, f) == n_bitmap &&
           fread(vals.data(), sizeof(cf_t), n_valid, f) == n_valid;
      uint32_t j = 0;
      for (uint32_t k = 0; ok && k < n_sc; k++) {
        if (bitmap[k / 8] & (1u << (k % 8))) {
          ok = j < n_valid;
          if (ok) {
            H[(size_t)l * n_sc + k]     = vals[j++];
            valid[(size_t)l * n_sc + k] = true;
          }
        }
      }
      ok = ok && j == n_valid;
    }
    if (!ok) {
      truncated = true;
      break;
    }

    if (n_records == 0) {
      slot_first = rh.slot_abs;
    }
    if (!targets.empty()) {
      const double t_slot = (double)(int64_t)(rh.slot_abs - slot_first) * t_slot_s;
      inject_record(targets, rh, fh, t_slot, t_slot_s / NSYMB, H.data(), valid);
      if (rh.slot_abs != truth_slot) {
        truth_slot = rh.slot_abs;
        for (size_t i = 0; i < targets.size(); i++) {
          const inject_t& x = targets[i];
          fprintf(truth, "%.6f,%u,%u,%lu,%zu,%d,%.3f,%.4f,%.1f,%.1f\n", t_slot, rh.sfn, rh.slot_idx,
                  (unsigned long)rh.slot_abs, i, (int)x.active(t_slot), x.range_at(t_slot), x.speed_at(t_slot),
                  x.level_db, x.angle_deg);
        }
      }
    }

    // set from SIB1 live, so it is whatever it was when this grant arrived
    if (rh.tdd_period_slots > 0) {
      nr_sensing_tdd_period_slots = rh.tdd_period_slots;
    }
    n_pushed += (uint64_t)nrscope_sensing_push_estimates(s, rh.aarx, rh.ports, rh.layer, rh.n_layers, rh.sfn,
                                                         rh.slot_idx, rh.slot_abs, rh.dmrs_mask, H.data(), valid,
                                                         n_sc);
    slot_last = rh.slot_abs;
    n_records++;
  }
  fclose(f);
  if (truncated) {
    fprintf(stderr, "record %lu is incomplete (the recording was cut short); replayed up to it\n",
            (unsigned long)n_records);
  }

  nrscope_sensing_wait_maps(s);
  if (truth != nullptr)
    fclose(truth);
  const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_start).count();
  printf("replayed %lu records, %lu snapshots pushed, %.2f s of air (slots %lu..%lu) in %.1f s\n",
         (unsigned long)n_records, (unsigned long)n_pushed, (slot_last - slot_first) * 1e-3 / (fh.scs_hz / 15000),
         (unsigned long)slot_first, (unsigned long)slot_last, secs);
  free(valid);
  return EXIT_SUCCESS;
}
