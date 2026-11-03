/*
 * SPDX-License-Identifier: LicenseRef-CSSL-1.0
 */
/*
 * \brief GoldenSniffer reads a tracked UE's DM-RS structure on air (gs_dmrs_shape /
 * gs_dmrs_decide): its antenna-ports field and RRC are ciphered. Grids are built from the
 * PCI sequence with srsRAN's mapper for port 1000, port 1001 as its covered copy, comb 1
 * empty, carrying PDSCH, or carrying ports 1002/1003, and extra DM-RS symbols, through a
 * multipath channel at several SNRs. Prints the share decided right per case.
 */
#include <complex>
#include <cstdio>
#include <cstring>
#include <random>

#include "nrscope/hdr/gs_fdra.h"
extern "C" {
#include "srsran/phy/ch_estimation/dmrs_sch.h"
}

static int failures = 0;

int main()
{
  const uint32_t nprb = 273, nsc = nprb * 12, pci = 1;
  srsran_carrier_nr_t carrier = SRSRAN_DEFAULT_CARRIER_NR;
  carrier.pci = pci, carrier.nof_prb = nprb, carrier.scs = srsran_subcarrier_spacing_30kHz;
  srsran_dmrs_sch_t dmrs = {};
  srsran_dmrs_sch_init(&dmrs, false);
  srsran_dmrs_sch_set_carrier(&dmrs, &carrier);
  static cf_t                        tx[273 * 12 * 14], rx[273 * 12 * 14];
  std::mt19937                       gen(7);
  std::normal_distribution<float>    nd(0.0f, (float)M_SQRT1_2);
  std::uniform_real_distribution<float> ud(0.0f, 1.0f);
  /* case: ports on comb 0 and what comb 1 holds */
  struct Case {
    const char* name;
    uint16_t    ports;
    int         comb1; // 0 empty (CDM groups 2), 1 PDSCH (CDM groups 1), 2 DM-RS of ports 1002
  } cases[] = {{"port 1000, comb 1 PDSCH", 0x1, 1},       {"port 1000, comb 1 empty", 0x1, 0},
               {"port 1001, comb 1 PDSCH", 0x2, 1},       {"ports 1000+1001, comb 1 PDSCH", 0x3, 1},
               {"ports 1000+1001, comb 1 empty", 0x3, 0}, {"ports 1000+1001, 1002+ on comb 1", 0x3, 2}};
  const float snrs[] = {0, 5, 10, 20};
  printf("%-34s %9s %9s %9s %9s   (safe/exact, 100 grants of 4..60 PRB, 1 us of delay spread)\n", "case",
         "0 dB", "5 dB", "10 dB", "20 dB");
  double right_10[6] = {};
  for (int ci = 0; ci < 6; ci++) {
    printf("%-34s", cases[ci].name);
    for (int si = 0; si < 4; si++) {
      int ok = 0, exact = 0;
      for (int t = 0; t < 100; t++) {
        const uint32_t n_rb = 4 + (uint32_t)(ud(gen) * 56), rb0 = (uint32_t)(ud(gen) * (nprb - n_rb));
        srsran_sch_cfg_nr_t cfg = {};
        cfg.dmrs.type = srsran_dmrs_sch_type_1, cfg.dmrs.typeA_pos = srsran_dmrs_sch_typeA_pos_2;
        cfg.dmrs.additional_pos = srsran_dmrs_sch_add_pos_0, cfg.dmrs.length = srsran_dmrs_sch_len_1;
        cfg.grant.mapping = srsran_sch_mapping_type_A, cfg.grant.S = 1, cfg.grant.L = 13;
        cfg.grant.nof_layers = 1, cfg.grant.nof_dmrs_cdm_groups_without_data = 1, cfg.grant.beta_dmrs = 1.0f;
        bool prb[SRSRAN_MAX_PRB_NR] = {};
        for (uint32_t p = rb0; p < rb0 + n_rb; p++) {
          cfg.grant.prb_idx[p] = prb[p] = true;
        }
        cfg.grant.nof_prb       = n_rb;
        const uint32_t    slot_idx = (uint32_t)(ud(gen) * 20);
        srsran_slot_cfg_t slot     = {.idx = slot_idx};
        memset(tx, 0, sizeof(tx));
        srsran_dmrs_sch_put_sf(&dmrs, &slot, &cfg, &cfg.grant, tx);
        // two layers' channels, 6 paths up to 1 us, and the comb-1 layers'
        std::complex<float> g[4][6];
        float               tau[6];
        for (int p = 0; p < 6; p++) {
          tau[p] = p == 0 ? 0.2e-6f : 0.2e-6f + ud(gen) * 1e-6f;
          for (int a = 0; a < 4; a++) {
            g[a][p] = std::complex<float>(nd(gen), nd(gen)) * (p == 0 ? 1.0f : expf(-2 * ud(gen))) * 0.45f;
          }
        }
        /* layers of one transmission: a precoder scales them within 6 dB of each other */
        float lp[4];
        for (int a = 0; a < 4; a++) {
          float e = 0;
          for (int p = 0; p < 6; p++) {
            e += std::norm(g[a][p]);
          }
          lp[a] = 1.0f / sqrtf(e) * (a == 0 ? 1.0f : powf(10.0f, -0.3f * ud(gen)));
        }
        for (int a = 0; a < 4; a++) {
          for (int p = 0; p < 6; p++) {
            g[a][p] *= lp[a];
          }
        }
        const float sn = powf(10.0f, -snrs[si] / 20.0f);
        const uint32_t l = 2;
        for (uint32_t k = 0; k < nsc; k++) {
          std::complex<float> h[4] = {};
          for (int p = 0; p < 6; p++) {
            const std::complex<float> e = std::polar(1.0f, (float)(-2 * M_PI * 30e3 * k * tau[p]));
            for (int a = 0; a < 4; a++) {
              h[a] += g[a][p] * e;
            }
          }
          const bool in = prb[k / 12];
          std::complex<float> y = 0;
          if (in && (k % 2) == 0) { // comb 0: port 1000 and its covered copy 1001
            const std::complex<float> x  = std::complex<float>(tx[l * nsc + k]);
            const float               w1 = ((k / 2) % 2 == 1) ? -1.0f : 1.0f;
            y = x * ((cases[ci].ports & 1 ? h[0] : 0.0f) + (cases[ci].ports & 2 ? w1 * h[1] : 0.0f));
          } else if (in) { // comb 1
            if (cases[ci].comb1 == 1) { // PDSCH: the DM-RS layers carrying QPSK
              auto qpsk = [&]() {
                return std::complex<float>(ud(gen) < 0.5f ? 0.707f : -0.707f, ud(gen) < 0.5f ? 0.707f : -0.707f);
              };
              y = (cases[ci].ports & 1 ? h[0] * qpsk() : 0.0f) + (cases[ci].ports & 2 ? h[1] * qpsk() : 0.0f);
              y *= (float)M_SQRT1_2 * (__builtin_popcount(cases[ci].ports) == 2 ? (float)M_SQRT1_2 : 1.0f) * 1.41421356f;
            } else if (cases[ci].comb1 == 2) { // port 1002: the same sequence one subcarrier up
              y = h[2] * std::complex<float>(tx[l * nsc + k - 1]);
            }
          }
          const std::complex<float> v = y * 1.41421356f + sn * std::complex<float>(nd(gen), nd(gen)); // |x| = 1/sqrt2
          __real__ rx[l * nsc + k] = v.real();
          __imag__ rx[l * nsc + k] = v.imag();
        }
        GsDmrsShape sh;
        gs_dmrs_shape(rx + l * nsc, 0, prb, nprb, slot_idx, l, pci, &sh);
        uint16_t ports;
        bool     cdm2;
        int      dp;
        gs_dmrs_decide(sh, &ports, &cdm2, &dp);
        // safe: no layer that is on air left out, and the DM-RS power right
        ok += (ports & cases[ci].ports) == cases[ci].ports && cdm2 == (cases[ci].comb1 != 1);
        exact += ports == cases[ci].ports && cdm2 == (cases[ci].comb1 != 1);
      }
      printf(" %3d%%/%3d%%", ok, exact);
      if (si == 2) {
        right_10[ci] = ok / 100.0;
      }
    }
    printf("\n");
  }
  for (int ci = 0; ci < 6; ci++) {
    if (right_10[ci] < 0.95) {
      printf("FAIL: %s decided safely %.0f%% at 10 dB\n", cases[ci].name, 100 * right_10[ci]);
      failures++;
    }
  }
  srsran_dmrs_sch_free(&dmrs);
  printf("%s: %d failure(s)\n", failures ? "FAILED" : "PASSED", failures);
  return failures ? 1 : 0;
}
