#include "nrscope/hdr/sibs_decoder.h"
#include "nrscope/hdr/run_recorder.h"

#include <atomic>
#include <mutex>
#include <time.h>
#include <stdlib.h>

/* WHY THIS SUMMARY EXISTS

A SIB1 blind search that finds nothing used to return silently (the "No DCI found" line
is behind NRSCOPE_TRACE_PER_SLOT, which at 2000 slots/s cannot be left on). The whole
sniffer then sits there decoding nothing, because without SIB1 there is no RACH
configuration, so no TC-RNTI, no Msg4 and no UE DCIs -- and the only outward sign is
"0 RNTI(s)" in the status line, which looks exactly like a cell with no traffic.

One line a second, from whichever worker gets there first, with the counts and the best
PDCCH measurement of the window. It separates the cases that need completely different
answers: no candidate measured at all (nothing in CORESET#0, so the grid or the
frequency mapping is wrong), candidates with EPRE but a poor correlation (the samples
reaching the decoder are not the cell's), or a good correlation failing CRC (the
candidate is there and the payload is being corrupted). */
static struct {
  std::mutex mutex;
  uint64_t   slots;
  uint64_t   candidates;
  uint64_t   crc_fail;
  double     best_epre;
  double     best_corr;
  /* Aggregation level and CCE of the best-correlating candidate. The cell's real SIB1
  PDCCH sits at one fixed (L, ncce) -- L=2, ncce=0 on the Benetel cell -- so a best
  candidate anywhere else says the search is locking onto something that is not SIB1,
  which is a different fault from finding the right candidate and failing its CRC. */
  uint32_t   best_L;
  uint32_t   best_ncce;
  uint64_t   last_report;
} sib_search_stats = {{}, 0, 0, 0, -INFINITY, -INFINITY, 0, 0, 0};

SIBsDecoder::SIBsDecoder(){
  /* NR-Scope fills only some fields of the higher-layer PDSCH config and
    leaves the rest (the NZP/ZP CSI-RS sets among them) as "not configured".
    Without this they held whatever was on the heap, and the grant conversion
    walked garbage CSI-RS sets: "csi_rs.c: Unhandled configuration row=invalid",
    then "Error in resource mapping" and a grant cut short. */
  pdsch_hl_cfg = {};
  data_pdcch = srsran_vec_u8_malloc(SRSRAN_SLOT_MAX_NOF_BITS_NR);
  if (data_pdcch == NULL) {
    ERROR("Error malloc");
  }
}

SIBsDecoder::~SIBsDecoder(){
    
}

int SIBsDecoder::SIBDecoderandReceptionInit(WorkState* state,
                                            cf_t* input[SRSRAN_MAX_PORTS]){  
  std::cout << "Initializing SIB decoder" << std::endl;
  memcpy(&coreset0_t, &state->coreset0_t, sizeof(srsran_coreset_t));

  dci_cfg.bwp_dl_initial_bw   = 275;
  dci_cfg.bwp_ul_initial_bw   = 275;
  dci_cfg.bwp_dl_active_bw    = 275;
  dci_cfg.bwp_ul_active_bw    = 275;
  dci_cfg.monitor_common_0_0  = true;
  dci_cfg.monitor_0_0_and_1_0 = true;
  dci_cfg.monitor_0_1_and_1_1 = true;
  // set coreset0 bandwidth
  dci_cfg.coreset0_bw = srsran_coreset_get_bw(&coreset0_t);

  pdcch_cfg.coreset_present[0] = true;
  search_space = &pdcch_cfg.search_space[0];
  pdcch_cfg.search_space_present[0]   = true;
  search_space->id                    = 0;
  search_space->coreset_id            = 0;
  search_space->type                  = srsran_search_space_type_common_0;
  search_space->formats[0]            = srsran_dci_format_nr_1_0;
  search_space->nof_formats           = 1;
  for (uint32_t L = 0; L < SRSRAN_SEARCH_SPACE_NOF_AGGREGATION_LEVELS_NR; L++) {
    search_space->nof_candidates[L] = 
      srsran_pdcch_nr_max_candidates_coreset(&coreset0_t, L);
  }
  pdcch_cfg.coreset[0] = coreset0_t; 

  arg_scs = state->arg_scs;
  memcpy(&base_carrier, &state->args_t.base_carrier, sizeof(srsran_carrier_nr_t));
  cell = state->cell;
  pdsch_hl_cfg.typeA_pos = cell.mib.dmrs_typeA_pos;

  ue_dl_args.nof_rx_antennas               = 1;
  ue_dl_args.pdsch.sch.disable_simd        = false;
  ue_dl_args.pdsch.sch.decoder_use_flooded = false;
  ue_dl_args.pdsch.measure_evm             = true;
  ue_dl_args.pdcch.disable_simd            = false;
  ue_dl_args.pdcch.measure_evm             = true;
  ue_dl_args.nof_max_prb                   = 275;

  ue_dl_args.pdcch_dmrs_corr_thr           = 0.05;

  if (srsran_ue_dl_nr_init_nrscope(&ue_dl_sibs, input, &ue_dl_args, arg_scs)) {
    ERROR("Error UE DL");
    return SRSRAN_ERROR;
  }

  if (srsran_ue_dl_nr_set_carrier_nrscope(&ue_dl_sibs, &base_carrier, arg_scs)){
    ERROR("Error setting SCH NR carrier");
    return SRSRAN_ERROR;
  }

  if (srsran_ue_dl_nr_set_pdcch_config(&ue_dl_sibs, &pdcch_cfg, &dci_cfg)) {
    ERROR("Error setting CORESET");
    return SRSRAN_ERROR;
  }

  if (srsran_softbuffer_rx_init_guru(&softbuffer, 
      SRSRAN_SCH_NR_MAX_NOF_CB_LDPC, SRSRAN_LDPC_MAX_LEN_ENCODED_CB) <
      SRSRAN_SUCCESS) {
    ERROR("Error init soft-buffer");
    return SRSRAN_ERROR;
  }

  // task_scheduler_nrscope->sib1_inited = true;
  /* Everything that decides WHERE in the grid this decoder looks for CORESET#0. The
  PDCCH search reports strong EPRE with noise-level correlation both when the samples
  are wrong and when the mapping is wrong, and these are the mapping: if they are
  identical across two runs that behave differently, the mapping is exonerated. */
  std::cout << "SIB Decoder Initialized.. srate=" << arg_scs.srate << " scs=" << (int)arg_scs.scs
            << " coreset_offset_scs=" << arg_scs.coreset_offset_scs << " coreset_slot=" << arg_scs.coreset_slot
            << " coreset0_bw=" << dci_cfg.coreset0_bw << " carrier_prb=" << base_carrier.nof_prb
            << " abs_pdcch_scs=" << cell.abs_pdcch_scs << " pci=" << base_carrier.pci << std::endl;

  return SRSRAN_SUCCESS;
}

int SIBsDecoder::DecodeandParseSIB1fromSlot(srsran_slot_cfg_t* slot,
                                            WorkState* state,
                                            SlotResult* result) {
  if (state->all_sibs_found) {
    result->sib_result = false;
    return SRSRAN_SUCCESS;
  }
  
  /* reset the result's elements */
  result->sib_result = true;
  result->found_sib1 = false;
  result->found_sib.clear();
  result->sibs.clear();

  // uint8_t hardcoded_sib1[111] = {0x6e, 0x88, 0x18, 0x0a, 0x09, 0x88, 0x49, 
  //   0x81, 0x5b, 0x98, 0x00, 0x51, 0x43, 0x68, 0x00, 0x68, 0x26, 0x23, 0x49, 
  //   0x05, 0x6e, 0x60, 0x01, 0x45, 0x0d, 0xa0, 0x01, 0x00, 0x0c, 0x50, 0x01, 
  //   0x08, 0x30, 0xc0, 0x4a, 0x08, 0x46, 0x01, 0x40, 0x00, 0x00, 0x4e, 0x65, 
  //   0x3c, 0xa1, 0x0f, 0x9b, 0x82, 0x01, 0x00, 0x00, 0x00, 0x84, 0x00, 0x02, 
  //   0x08, 0x34, 0x20, 0x21, 0x06, 0x10, 0x01, 0x9c, 0x21, 0x18, 0x20, 0x64, 
  //   0x60, 0x00, 0x00, 0x4e, 0xd5, 0xca, 0x79, 0x42, 0xd2, 0x11, 0x00, 0x78, 
  //   0xc0, 0xd4, 0x41, 0x99, 0x00, 0x5f, 0xcb, 0xce, 0x08, 0xdc, 0x21, 0xb8, 
  //   0x63, 0x71, 0x06, 0xfb, 0x65, 0xf1, 0x21, 0x16, 0x3c, 0x40, 0x54, 0x55, 
  //   0xee, 0x43, 0x14, 0x41, 0x0e, 0x00, 0x00, 0x00};
  
  memset(&dci_sibs, 0, sizeof(srsran_dci_dl_nr_t));

  // FILE *fp;
  // fp = fopen("/home/wanhr/Documents/codes/cpp/srsRAN_4G/build/srsue/src/SIB_debug.txt", "r");
  // fseek(fp, file_position * sizeof(cf_t), SEEK_SET);
  // uint32_t a = fread(ue_dl.fft[0].cfg.in_buffer, sizeof(cf_t), ue_dl.fft[0].sf_sz, fp);
  // uint32_t a = fread(ue_dl_sibs.fft[0].cfg.in_buffer, sizeof(cf_t), ue_dl_sibs.fft[0].sf_sz, fp);
  // file_position += ue_dl_sibs.fft[0].sf_sz;
  
  // Check the fft plan and how does it manipulate the buffer
  /* Diagnostic: the exact time-domain slot this decoder is about to demodulate, as the
  worker handed it over. Gated on NRSCOPE_DUMP_SIB_SLOTS=<count>, off when unset, so it
  costs nothing in a normal run.

  It exists because a SIB1 that never decodes cannot be told apart, from the outside,
  from a cell that is not transmitting one: the PDCCH search reports energy with no
  correlation either way. Dumping what reaches the FFT settles whether the samples are
  the cell's and properly aligned (fault downstream, in the mapping or the descrambling)
  or not (fault upstream, in the chain-0 path through the slot queue).

  File: /tmp/sib_slots.bin, one record per slot: sfn, slot_idx, n_samples as uint32,
  then n_samples interleaved float re/im. Read it with
  scripts/benchmark/sib_slot_check.py. */
  {
    /* Read once, thread-safely (C++11 guarantees function-local static init), so the
      count cannot be reset by a race between workers. */
    static const int dump_total = [] {
      const char* e = getenv("NRSCOPE_DUMP_SIB_SLOTS");
      return e != NULL ? atoi(e) : 0;
    }();
    static const char* dump_path = [] {
      const char* p = getenv("NRSCOPE_DUMP_SIB_PATH");
      return p != NULL ? p : "/tmp/sib_slots.bin";
    }();
    static std::atomic<int> dumped{0};
    static std::mutex       dump_mutex;
    if (dump_total > 0) {
      const int mine = dumped.fetch_add(1, std::memory_order_relaxed);
      if (mine < dump_total) {
        std::lock_guard<std::mutex> lock(dump_mutex);
        FILE*                       f = fopen(dump_path, mine == 0 ? "wb" : "ab");
        if (f == NULL) {
          if (mine == 0)
            printf("SIBDecoder: cannot open %s for the slot dump\n", dump_path);
        } else {
          /* The FFT's own input pointer and length, so this is exactly the buffer it
            reads and exactly how much of it, not a guess at either. */
          const cf_t*    in     = ue_dl_sibs.fft[0].cfg.in_buffer;
          const uint32_t n      = ue_dl_sibs.fft[0].sf_sz;
          const uint32_t hdr[3] = {state->sfn, slot->idx, n};
          fwrite(hdr, sizeof(uint32_t), 3, f);
          fwrite(in, sizeof(cf_t), n, f);
          fclose(f);
          if (mine == 0)
            printf("SIBDecoder: dumping %d slot(s) of %u samples to %s\n", dump_total, n, dump_path);
          if (mine == dump_total - 1)
            printf("SIBDecoder: slot dump complete (%d slots) in %s\n", dump_total, dump_path);
        }
      }
    }
  }
  srsran_ue_dl_nr_estimate_fft_nrscope(&ue_dl_sibs, slot, arg_scs);
  // Blind search
  int nof_found_dci = srsran_ue_dl_nr_find_dl_dci(&ue_dl_sibs, slot, 0xFFFF, 
    srsran_rnti_type_si, &dci_sibs, 1);
  if (nof_found_dci < SRSRAN_SUCCESS){
    ERROR("SIBDecoder -- Error in blind search");
    return SRSRAN_ERROR;
  }
  /* Print PDCCH blind search candidates */
  for (uint32_t pdcch_idx = 0; NRSCOPE_TRACE_PER_SLOT && pdcch_idx < ue_dl_sibs.pdcch_info_count; pdcch_idx++) {
    const srsran_ue_dl_nr_pdcch_info_t* info = &(ue_dl_sibs.pdcch_info[pdcch_idx]);
    printf("PDCCH: %s-rnti=0x%x, crst_id=%d, ss_type=%s, ncce=%d, al=%d, EPRE=%+.2f, RSRP=%+.2f, corr=%.3f; "
    "nof_bits=%d; crc=%s;\n",
    srsran_rnti_type_str_short(info->dci_ctx.rnti_type),
    info->dci_ctx.rnti,
    info->dci_ctx.coreset_id,
    srsran_ss_type_str(info->dci_ctx.ss_type),
    info->dci_ctx.location.ncce,
    info->dci_ctx.location.L,
    info->measure.epre_dBfs,
    info->measure.rsrp_dBfs,
    info->measure.norm_corr,
    info->nof_bits,
    info->result.crc ? "OK" : "KO");
  }
  if (nof_found_dci < 1) {
    NRSCOPE_SLOT_TRACE("SIBDecoder -- No DCI found :'(\n");
    /* Summarise the failure once a second rather than per slot; see sib_search_stats. */
    bool     report = false;
    uint64_t slots = 0, candidates = 0, crc_fail = 0;
    double   best_epre = 0, best_corr = 0;
    uint32_t best_L = 0, best_ncce = 0;
    {
      std::lock_guard<std::mutex> lock(sib_search_stats.mutex);
      sib_search_stats.slots++;
      for (uint32_t i = 0; i < ue_dl_sibs.pdcch_info_count; i++) {
        const srsran_dmrs_pdcch_measure_t& m = ue_dl_sibs.pdcch_info[i].measure;
        sib_search_stats.candidates++;
        sib_search_stats.crc_fail += ue_dl_sibs.pdcch_info[i].result.crc ? 0 : 1;
        if (isnormal(m.norm_corr)) {
          sib_search_stats.best_epre = SRSRAN_MAX(sib_search_stats.best_epre, m.epre_dBfs);
          if (m.norm_corr > sib_search_stats.best_corr) {
            sib_search_stats.best_corr = m.norm_corr;
            sib_search_stats.best_L    = ue_dl_sibs.pdcch_info[i].dci_ctx.location.L;
            sib_search_stats.best_ncce = ue_dl_sibs.pdcch_info[i].dci_ctx.location.ncce;
          }
        }
      }
      const uint64_t now = (uint64_t)time(NULL);
      if (now != sib_search_stats.last_report) {
        sib_search_stats.last_report = now;
        report = true, slots = sib_search_stats.slots, candidates = sib_search_stats.candidates;
        crc_fail = sib_search_stats.crc_fail;
        best_epre = sib_search_stats.best_epre, best_corr = sib_search_stats.best_corr;
        best_L = sib_search_stats.best_L, best_ncce = sib_search_stats.best_ncce;
        sib_search_stats.slots = sib_search_stats.candidates = sib_search_stats.crc_fail = 0;
        sib_search_stats.best_epre = sib_search_stats.best_corr = -INFINITY;
        sib_search_stats.best_L = sib_search_stats.best_ncce = 0;
      }
    }
    if (report) {
      printf("SIBDecoder: no SIB1 in %lu slot(s), %lu PDCCH candidate(s) measured (%lu CRC fail), "
             "best EPRE %+.1f dBfs, best corr %.3f at L=%u ncce=%u (threshold %.3f)\n",
             (unsigned long)slots,
             (unsigned long)candidates,
             (unsigned long)crc_fail,
             best_epre,
             best_corr,
             best_L,
             best_ncce,
             ue_dl_sibs.pdcch_dmrs_corr_thr);
    }
    return SRSRAN_ERROR;
  }

  char str[1024] = {};
  srsran_dci_dl_nr_to_str(&(ue_dl_sibs.dci), &dci_sibs, str, 
    (uint32_t)sizeof(str));
  printf("SIBDecoder -- Found DCI: %s\n", str);

  pdsch_cfg = {};
  pdsch_cfg.dmrs.typeA_pos = state->cell.mib.dmrs_typeA_pos;

  if (srsran_ra_dl_dci_to_grant_nr(&base_carrier, slot, &pdsch_hl_cfg, 
      &dci_sibs, &pdsch_cfg, &pdsch_cfg.grant) < SRSRAN_SUCCESS) {
    ERROR("SIBDecoder -- Error decoding PDSCH search");
    return SRSRAN_ERROR;
  }

  srsran_sch_cfg_nr_info(&pdsch_cfg, str, (uint32_t)sizeof(str));
  if (RunRecorder::verbose()) printf("PDSCH_cfg:\n%s", str);

  if (srsran_softbuffer_rx_init_guru(&softbuffer, SRSRAN_SCH_NR_MAX_NOF_CB_LDPC, 
      SRSRAN_LDPC_MAX_LEN_ENCODED_CB) < SRSRAN_SUCCESS) {
    ERROR("SIBDecoder -- Error init soft-buffer");
    return SRSRAN_ERROR;
  }

  // Reset the data_pdcch to zeros
  srsran_vec_u8_zero(data_pdcch, SRSRAN_SLOT_MAX_NOF_BITS_NR);
  
  pdsch_cfg.grant.tb[0].softbuffer.rx = &softbuffer; // Set softbuffer
  pdsch_res = {}; // Prepare PDSCH result
  pdsch_res.tb[0].payload = data_pdcch;

  // Decode PDSCH
  if (srsran_ue_dl_nr_decode_pdsch(&ue_dl_sibs, slot, &pdsch_cfg, &pdsch_res) < 
      SRSRAN_SUCCESS) {
    printf("SIBDecoder -- Error decoding PDSCH search\n");
    return SRSRAN_ERROR;
  }
  if (!pdsch_res.tb[0].crc) {
    printf("SIBDecoder -- Error decoding PDSCH (CRC)\n");
    return SRSRAN_ERROR;
  }
  // printf("Decoded PDSCH (%d B)\n", pdsch_cfg.grant.tb[0].tbs / 8);
  // srsran_vec_fprint_byte(stdout, pdsch_res.tb[0].payload, 
  // pdsch_cfg.grant.tb[0].tbs / 8);

   // check payload is not all null
  bool all_zero = true;
  for (int i = 0; i < pdsch_cfg.grant.tb[0].tbs / 8; ++i) {
    if (pdsch_res.tb[0].payload[i] != 0x0) {
      all_zero = false;
      break;
    }
  }
  if (all_zero) {
    ERROR("PDSCH payload is all zeros");
    return SRSRAN_ERROR;
  }
  std::cout << "Try to decode SIBs..." << std::endl;
  asn1::rrc_nr::bcch_dl_sch_msg_s dlsch_msg;
  asn1::cbit_ref dlsch_bref(pdsch_res.tb[0].payload, 
    pdsch_cfg.grant.tb[0].tbs / 8);
  asn1::SRSASN_CODE err = dlsch_msg.unpack(dlsch_bref);

  /* Try to decode the SIB, we need a better way to provide the result */
  if(srsran_unlikely(asn1::rrc_nr::bcch_dl_sch_msg_type_c::c1_c_::types_opts::
      sib_type1 != dlsch_msg.msg.c1().type())){
    // Try to decode other SIBs
    // Get the sib_id, sib_id is uint8_t and is 2 for sib2, 3 for sib 3, etc...
    auto sib_id = dlsch_msg.msg.c1().sys_info().crit_exts.sys_info().
      sib_type_and_info[0].type().to_number();
    auto decoded_sib = dlsch_msg.msg.c1().sys_info();
    result->sibs.emplace_back(decoded_sib);
    result->found_sib.emplace_back(sib_id);
    // sibs[sib_id - 2] = dlsch_msg.msg.c1().sys_info();
    // found_sib[sib_id - 2] = 1;

    // If we collect all the SIBs, we can skip the thread.
    // long unsigned int found_result = 0;
    // for(long unsigned int i=0; i<found_sib.size(); i++){
    //   found_result += found_sib[i];
    // }
    // if(found_result >= found_sib.size()){
    //   *all_sibs_found = true;
    // }

    std::cout << "SIB " << (int)sib_id << " Decoded." << std::endl;
    /* Uncomment to print the decode SIBs. */
    asn1::json_writer js_sibs;
    (decoded_sib).to_json(js_sibs);
    if (RunRecorder::verbose()) printf("Decoded SIBs: %s\n", js_sibs.to_string().c_str());
  }else if(srsran_unlikely(asn1::rrc_nr::bcch_dl_sch_msg_type_c::c1_c_::
      types_opts::sys_info != dlsch_msg.msg.c1().type())){
    result->found_sib1 = true;
    result->sib1 = dlsch_msg.msg.c1().sib_type1();
    std::cout << "SIB 1 Decoded." << std::endl;

    /* Record where this SIB1 PDSCH landed, so the cell map can draw it. The
      grant's prb_idx is relative to the CORESET#0 grid used for SI decoding. */
    uint32_t first_prb = SRSRAN_MAX_PRB_NR;
    for (uint32_t i = 0; i < SRSRAN_MAX_PRB_NR && first_prb == SRSRAN_MAX_PRB_NR; i++) {
      if (pdsch_cfg.grant.prb_idx[i]) {
        first_prb = i;
      }
    }
    result->sib1_prb_start    = (first_prb == SRSRAN_MAX_PRB_NR) ? 0 : first_prb;
    result->sib1_nof_prb      = pdsch_cfg.grant.nof_prb;
    result->sib1_symbol_start = pdsch_cfg.grant.S;
    result->sib1_nof_symbols  = pdsch_cfg.grant.L;
    result->sib1_slot_idx     = slot->idx;
    result->sib1_pdcch_L      = dci_sibs.ctx.location.L;
    result->sib1_pdcch_ncce   = dci_sibs.ctx.location.ncce;

    /* Uncomment to print the decode SIB1. */
    asn1::json_writer js_sib1;
    (result->sib1).to_json(js_sib1);
    if (RunRecorder::verbose()) printf("Decoded SIB1: %s\n", js_sib1.to_string().c_str());
  }

  srsran_softbuffer_rx_free(&softbuffer);

  return SRSRAN_SUCCESS;
}
