#include "nrscope/hdr/nrscope_worker.h"
#include <chrono>
#include <atomic>
#include <ctime>
#include <mutex>
#include <semaphore>

namespace NRScopeTask {

std::vector<SlotResult> global_slot_results;
std::mutex              queue_lock;
std::mutex              task_scheduler_lock;
std::mutex              slot_data_lock;
std::mutex              worker_locks[128];

sem_t smph_data; // counts ready slots
sem_t smph_idle; // counts idle workers

NRScopeWorker::NRScopeWorker() : rf_buffer_t(1), rach_decoder(), sibs_decoder()
{
  worker_state.sib1_inited = false;
  worker_state.rach_inited = false;
  worker_state.dci_inited  = false;

  worker_state.sib1_found = false;
  worker_state.rach_found = false;

  initializing = false;

  worker_state.nof_known_rntis = 0;
  worker_state.known_rntis.resize(worker_state.nof_known_rntis);
  sem_init(&smph_has_job, 0, 0);
  busy = false;
}

NRScopeWorker::~NRScopeWorker() {}

int NRScopeWorker::InitWorker(WorkState task_scheduler_state, int worker_id_)
{
  worker_id = worker_id_;
  /* Copy initial values */
  worker_state.nof_threads            = task_scheduler_state.nof_threads;
  worker_state.nof_rnti_worker_groups = task_scheduler_state.nof_rnti_worker_groups;
  worker_state.nof_bwps               = task_scheduler_state.nof_bwps;
  worker_state.args_t                 = task_scheduler_state.args_t;
  worker_state.slot_sz                = task_scheduler_state.slot_sz;
  worker_state.nof_antennas           = task_scheduler_state.nof_antennas;
  worker_state.cpu_affinity           = task_scheduler_state.cpu_affinity;
  worker_state.rrc_recfg_user         = task_scheduler_state.rrc_recfg_user;

  /* Size of one subframe, per receive chain */
  const uint32_t buf_samples = SRSRAN_NOF_SLOTS_PER_SF_NR(worker_state.args_t.ssb_scs) * worker_state.slot_sz;
  for (uint32_t a = 0; a < NRSCOPE_MAX_RX_ANTENNAS; a++) {
    rx_buffer[a] = (a < worker_state.nof_antennas) ? srsran_vec_cf_malloc(buf_samples) : nullptr;
  }

  /* A wrapper for the synchronised chain. Only channel 0 is populated: the
    decoders below are all initialised with nof_rx_antennas = 1 and read
    input[0], so handing them the other chains would be misleading rather than
    useful. The sensing path reads rx_buffer[] directly instead. */
  rf_buffer_t = srsran::rf_buffer_t(rx_buffer[0], buf_samples);
  /* Start the worker thread */
  // std::cout << "Starting the worker..." << std::endl;
  StartWorker();
  return SRSRAN_SUCCESS;
}

void NRScopeWorker::StartWorker()
{
  // std::cout << "Creating the thread. " << std::endl;
  if (worker_state.cpu_affinity) {
    cpu_set_t cpu_set_worker;
    CPU_ZERO(&cpu_set_worker);
    CPU_SET(worker_id * (3 + worker_state.nof_threads), &cpu_set_worker);
    worker_thread = std::thread{&NRScopeWorker::Run, this};
    assert(pthread_setaffinity_np(worker_thread.native_handle(), sizeof(cpu_set_t), &cpu_set_worker) == 0);
  } else {
    worker_thread = std::thread{&NRScopeWorker::Run, this};
  }

  worker_thread.detach();
}

void NRScopeWorker::CopySlotandBuffer(uint64_t                    sf_round_,
                                      srsran_slot_cfg_t           slot_,
                                      srsran_ue_sync_nr_outcome_t outcome_,
                                      cf_t* const*                rx_buffer_)
{
  sf_round = sf_round_;
  slot     = slot_;
  outcome  = outcome_;
  for (uint32_t a = 0; a < worker_state.nof_antennas; a++) {
    srsran_vec_cf_copy(rx_buffer[a], rx_buffer_[a], worker_state.slot_sz);
  }
}

int NRScopeWorker::PrewarmDecoders(WorkState* task_scheduler_state)
{
  if (SyncState(task_scheduler_state) < SRSRAN_SUCCESS) {
    return SRSRAN_ERROR;
  }
  if (!worker_state.sib1_inited) {
    if (InitSIBDecoder() < SRSRAN_SUCCESS) {
      return SRSRAN_ERROR;
    }
    worker_state.sib1_inited = true;
  }
  return SRSRAN_SUCCESS;
}

int NRScopeWorker::InitSIBDecoder()
{
  /* Will always be called before any tasks */
  initializing = true;
  if (sibs_decoder.SIBDecoderandReceptionInit(&worker_state, rf_buffer_t.to_cf_t()) < SRSASN_SUCCESS) {
    return SRSRAN_ERROR;
  }
  std::cout << "SIBs decoder initialized..." << std::endl;
  initializing = false;
  return SRSRAN_SUCCESS;
}

int NRScopeWorker::InitRACHDecoder()
{
  // std::thread rach_init_thread {&RachDecoder::rach_decoder_init,
  //  &rach_decoder, task_scheduler_nrscope.sib1, args_t.base_carrier};
  initializing = true;
  rach_decoder.RACHDecoderInit(worker_state);
  if (rach_decoder.RACHReceptionInit(&worker_state, rf_buffer_t.to_cf_t()) < SRSASN_SUCCESS) {
    ERROR("RACHDecoder Init Error");
    return SRSRAN_ERROR;
  }
  std::cout << "RACH decoder initialized.." << std::endl;
  initializing = false;
  return SRSRAN_SUCCESS;
}

int NRScopeWorker::InitDCIDecoders()
{
  initializing = true;
  sharded_results.resize(worker_state.nof_threads);
  nof_sharded_rntis.resize(worker_state.nof_threads);
  sharded_rntis.resize(worker_state.nof_threads);
  results.resize(worker_state.nof_bwps);
  for (uint32_t i = 0; i < worker_state.nof_rnti_worker_groups; i++) {
    // for each rnti worker group, for each bwp, spawn a decoder
    for (uint8_t j = 0; j < worker_state.nof_bwps; j++) {
      DCIDecoder* decoder = new DCIDecoder(100);
      if (decoder->DCIDecoderandReceptionInit(&worker_state, j, rf_buffer_t.to_cf_t()) < SRSASN_SUCCESS) {
        ERROR("DCIDecoder Init Error");
        return SRSRAN_ERROR;
      }
      decoder->dci_decoder_id       = i * worker_state.nof_bwps + j;
      decoder->rnti_worker_group_id = i;
      dci_decoders.push_back(std::unique_ptr<DCIDecoder>(decoder));
    }
  }
  initializing = false;
  return SRSRAN_SUCCESS;
}

int NRScopeWorker::SyncState(WorkState* task_scheduler_state)
{
  worker_state.args_t = task_scheduler_state->args_t;
  worker_state.cell   = task_scheduler_state->cell;
  worker_state.cs_ret = task_scheduler_state->cs_ret;
  memcpy(&worker_state.srsran_searcher_cfg_t,
         &task_scheduler_state->srsran_searcher_cfg_t,
         sizeof(srsue::nr::cell_search::cfg_t));
  worker_state.coreset0_args_t = task_scheduler_state->coreset0_args_t;
  memcpy(&worker_state.coreset0_t, &task_scheduler_state->coreset0_t, sizeof(srsran_coreset_t));

  worker_state.arg_scs = task_scheduler_state->arg_scs;

  worker_state.sib1 = task_scheduler_state->sib1;
  worker_state.sibs.resize(task_scheduler_state->sibs.size());
  for (long unsigned int i = 0; i < task_scheduler_state->sibs.size(); i++) {
    worker_state.sibs[i] = task_scheduler_state->sibs[i];
  }
  worker_state.found_sib.resize(task_scheduler_state->found_sib.size());
  for (long unsigned int i = 0; i < task_scheduler_state->found_sib.size(); i++) {
    worker_state.found_sib[i] = task_scheduler_state->found_sib[i];
  }
  worker_state.sibs_to_be_found.resize(task_scheduler_state->sibs_to_be_found.size());
  for (long unsigned int i = 0; i < task_scheduler_state->sibs_to_be_found.size(); i++) {
    worker_state.sibs_to_be_found[i] = task_scheduler_state->sibs_to_be_found[i];
  }

  worker_state.rrc_setup         = task_scheduler_state->rrc_setup;
  worker_state.master_cell_group = task_scheduler_state->master_cell_group;
  worker_state.rrc_recfg         = task_scheduler_state->rrc_recfg;

  /* SIB 1 decoded, we can start the RACH thread */
  worker_state.sib1_found = task_scheduler_state->sib1_found;
  worker_state.rach_found = task_scheduler_state->rach_found;

  worker_state.all_sibs_found = task_scheduler_state->all_sibs_found;

  worker_state.nof_known_rntis = task_scheduler_state->nof_known_rntis;
  worker_state.known_rntis.resize(worker_state.nof_known_rntis);
  for (long unsigned int i = 0; i < worker_state.nof_known_rntis; i++) {
    worker_state.known_rntis[i] = task_scheduler_state->known_rntis[i];
  }

  return SRSRAN_SUCCESS;
}

int NRScopeWorker::MergeResults()
{
  for (uint8_t b = 0; b < worker_state.nof_bwps; b++) {
    DCIFeedback new_result;
    results[b] = new_result;
    results[b].dl_grants.resize(worker_state.nof_known_rntis);
    results[b].ul_grants.resize(worker_state.nof_known_rntis);
    results[b].spare_dl_prbs.resize(worker_state.nof_known_rntis);
    results[b].spare_dl_tbs.resize(worker_state.nof_known_rntis);
    results[b].spare_dl_bits.resize(worker_state.nof_known_rntis);
    results[b].spare_ul_prbs.resize(worker_state.nof_known_rntis);
    results[b].spare_ul_tbs.resize(worker_state.nof_known_rntis);
    results[b].spare_ul_bits.resize(worker_state.nof_known_rntis);
    results[b].dl_dcis.resize(worker_state.nof_known_rntis);
    results[b].ul_dcis.resize(worker_state.nof_known_rntis);

    uint32_t rnti_s = 0;
    uint32_t rnti_e = 0;
    for (uint32_t i = 0; i < worker_state.nof_rnti_worker_groups; i++) {
      if (rnti_s >= worker_state.nof_known_rntis) {
        continue;
      }
      uint32_t n_rntis = nof_sharded_rntis[i * worker_state.nof_bwps];
      rnti_e           = rnti_s + n_rntis;
      if (rnti_e > worker_state.nof_known_rntis) {
        rnti_e = worker_state.nof_known_rntis;
      }

      uint32_t thread_id = i * worker_state.nof_bwps + b;
      results[b].nof_dl_used_prbs += sharded_results[thread_id].nof_dl_used_prbs;
      results[b].nof_ul_used_prbs += sharded_results[thread_id].nof_ul_used_prbs;

      /* Concatenated rather than indexed by RNTI position: this is every grant
      the slot carried, and a UE may appear in it more than once. Each entry names
      its own RNTI, so the shard offset does not apply. */
      results[b].all_dl_grants.insert(results[b].all_dl_grants.end(),
                                      sharded_results[thread_id].all_dl_grants.begin(),
                                      sharded_results[thread_id].all_dl_grants.end());

      for (uint32_t k = 0; k < n_rntis; k++) {
        results[b].dl_dcis[k + rnti_s]   = sharded_results[thread_id].dl_dcis[k];
        results[b].ul_dcis[k + rnti_s]   = sharded_results[thread_id].ul_dcis[k];
        results[b].dl_grants[k + rnti_s] = sharded_results[thread_id].dl_grants[k];
        results[b].ul_grants[k + rnti_s] = sharded_results[thread_id].ul_grants[k];
      }
      rnti_s = rnti_e;
    }

    /* TO-DISCUSS: to obtain even more precise result,
      here maybe we should total user payload prb in that bwp - used prb */
    results[b].nof_dl_spare_prbs = worker_state.args_t.base_carrier.nof_prb * (14 - 2) - results[b].nof_dl_used_prbs;
    for (uint32_t idx = 0; idx < worker_state.nof_known_rntis; idx++) {
      results[b].spare_dl_prbs[idx] = results[b].nof_dl_spare_prbs / worker_state.nof_known_rntis;
      if (abs(results[b].spare_dl_prbs[idx]) > worker_state.args_t.base_carrier.nof_prb * (14 - 2)) {
        results[b].spare_dl_prbs[idx] = 0;
      }
      results[b].spare_dl_tbs[idx]  = (int)((float)results[b].spare_dl_prbs[idx] * dl_prb_rate[idx]);
      results[b].spare_dl_bits[idx] = (int)((float)results[b].spare_dl_prbs[idx] * dl_prb_bits_rate[idx]);
    }

    results[b].nof_ul_spare_prbs = worker_state.args_t.base_carrier.nof_prb * (14 - 2) - results[b].nof_ul_used_prbs;
    for (uint32_t idx = 0; idx < worker_state.nof_known_rntis; idx++) {
      results[b].spare_ul_prbs[idx] = results[b].nof_ul_spare_prbs / worker_state.nof_known_rntis;
      if (abs(results[b].spare_ul_prbs[idx]) > worker_state.args_t.base_carrier.nof_prb * (14 - 2)) {
        results[b].spare_ul_prbs[idx] = 0;
      }
      results[b].spare_ul_tbs[idx]  = (int)((float)results[b].spare_ul_prbs[idx] * ul_prb_rate[idx]);
      results[b].spare_ul_bits[idx] = (int)((float)results[b].spare_ul_prbs[idx] * ul_prb_bits_rate[idx]);
    }
  }

  return SRSRAN_SUCCESS;
}

void NRScopeWorker::Run()
{
  while (true) {
    /* When there is a job, the semaphore is set and buffer is copied */
    sem_wait(&smph_has_job);
    // worker_locks[worker_id].lock();
    // busy = true;
    busy.store(true, std::memory_order_release);
    // worker_locks[worker_id].unlock();
    // (per-slot timing removed with the sub-threads it measured)

    if (NRSCOPE_TRACE_PER_SLOT)
      std::cout << "Processing sf_round: " << sf_round << ", sfn: " << outcome.sfn << ", slot.idx: " << slot.idx
              << std::endl;

    SlotResult slot_result = {};
    /* Set the all the results to be false, will be set inside the decoder
    threads */
    slot_result.sib_result  = false;
    slot_result.rach_result = false;
    slot_result.dci_result  = false;
    slot_result.slot        = slot;
    slot_result.outcome     = outcome;
    slot_result.sf_round    = sf_round;

    /* Decoder construction, at most one worker at a time.

    Building a decoder allocates srsRAN objects -- transforms, softbuffers,
    tables -- and takes long enough that the worker stops consuming slots while
    it happens. Every worker discovers it needs the same decoder in the same
    slot, so all sixteen used to stall together and the slot queue overflowed
    each time: once for SIB, once for RACH, and once for the DCI decoders when
    the first C-RNTI appeared. Those three moments accounted for every dropped
    slot in a seventy second run; the steady state never dropped one.

    SIB is normally built before capture starts, by PrewarmWorkers(). The other
    two cannot be, because they need a SIB1 and an RNTI that only exist once the
    capture is running. So they are serialised instead: a worker that cannot take
    the lock leaves the decoder for a later slot and processes this one without
    it, which is what it did before being initialised anyway. One worker stalls
    rather than sixteen. */
    static std::mutex init_lock;

    if (!worker_state.sib1_inited || (!worker_state.rach_inited && worker_state.sib1_found)
        || (!worker_state.dci_inited && worker_state.rach_found)) {
      std::unique_lock<std::mutex> lk(init_lock, std::try_to_lock);
      if (lk.owns_lock()) {
        if (!worker_state.sib1_inited) {
          InitSIBDecoder();
          worker_state.sib1_inited = true;
        }
        if (!worker_state.rach_inited && worker_state.sib1_found) {
          InitRACHDecoder();
          worker_state.rach_inited = true;
        }
        if (!worker_state.dci_inited && worker_state.rach_found) {
          InitDCIDecoders();
          worker_state.dci_inited = true;
        }
      }
    }

    /* SIB, RACH and DCI decoding run inline on this worker thread.
    
    They used to be spawned as std::threads per slot and joined immediately,
    which cost a create and a join for each -- measured at 15.6 us a pair on this
    machine, about 4000 a second across the pool. The concurrency bought nothing:
    the worker pool already runs one slot per worker, so sixteen slots are in
    flight at once, and splitting a single slot three ways on top of that only
    added syscalls and scheduler churn. Running them in order here also makes the
    worker a plain function again, which matters because the sensing path will
    add per-slot work to it. */

    /* If SIB1 is not found we decode it; once found this is skipped for good. */
    if (worker_state.sib1_inited && !worker_state.sib1_found) {
      sibs_decoder.DecodeandParseSIB1fromSlot(&slot, &worker_state, &slot_result);
    }

    if (worker_state.rach_inited) {
      rach_decoder.DecodeandParseMS4fromSlot(&slot, &worker_state, &slot_result);
    }

    if (worker_state.dci_inited) {
      slot_result.dci_result = true;

      dl_prb_rate.resize(worker_state.nof_known_rntis);
      ul_prb_rate.resize(worker_state.nof_known_rntis);
      dl_prb_bits_rate.resize(worker_state.nof_known_rntis);
      ul_prb_bits_rate.resize(worker_state.nof_known_rntis);

      for (uint32_t i = 0; i < worker_state.nof_threads; i++) {
        dci_decoders[i]->DecodeandParseDCIfromSlot(&slot,
                                                   &worker_state,
                                                   sharded_results,
                                                   sharded_rntis,
                                                   nof_sharded_rntis,
                                                   dl_prb_rate,
                                                   dl_prb_bits_rate,
                                                   ul_prb_rate,
                                                   ul_prb_bits_rate);
      }
    }

    if (worker_state.dci_inited) {
      MergeResults();
      slot_result.dci_feedback_results = results;

      /* How much DM-RS this slot actually offers, reported once a second.
      
      Grants, not UEs: each grant places its own pilots, so this is the number of
      independent channel estimates a slot can yield, and the PRB total is how
      much of the band they cover. Rate limited because it is a running health
      figure rather than an event. */
      {
        static std::atomic<uint64_t> n_slots{0}, n_grants{0}, n_prbs{0}, last_report{0};
        uint64_t                     g = 0, prb = 0;
        for (const auto& r : results) {
          g += r.all_dl_grants.size();
          for (const auto& rec : r.all_dl_grants) {
            prb += rec.grant.grant.nof_prb;
          }
        }
        n_slots.fetch_add(1, std::memory_order_relaxed);
        n_grants.fetch_add(g, std::memory_order_relaxed);
        n_prbs.fetch_add(prb, std::memory_order_relaxed);

        const uint64_t now  = (uint64_t)time(NULL);
        uint64_t       prev = last_report.load(std::memory_order_relaxed);
        if (now != prev && last_report.compare_exchange_strong(prev, now)) {
          const uint64_t sl = n_slots.exchange(0, std::memory_order_relaxed);
          const uint64_t gr = n_grants.exchange(0, std::memory_order_relaxed);
          const uint64_t pb = n_prbs.exchange(0, std::memory_order_relaxed);
          if (gr > 0) {
            printf("DM-RS available: %lu grant(s) over %lu slot(s), %.1f PRB per grant\n",
                   (unsigned long)gr,
                   (unsigned long)sl,
                   (double)pb / (double)gr);
          }
        }
      }
    }

    // std::cout << "After processing sf_round: " << sf_round << ", sfn: "
    //   << outcome.sfn << ", slot.idx: " << slot.idx << std::endl;
    // std::cout << "slot_result sf_round: " << slot_result.sf_round << ", sfn: "
    //   << slot_result.outcome.sfn << ", slot.idx: " << slot_result.slot.idx
    //   << std::endl;

    /* Post the result into the result queue*/
    queue_lock.lock();
    global_slot_results.push_back(slot_result);
    queue_lock.unlock();
    // worker_locks[worker_id].lock();
    busy.store(false, std::memory_order_release);
    sem_post(&smph_idle);
    // worker_locks[worker_id].unlock();
  }
}
} // namespace NRScopeTask