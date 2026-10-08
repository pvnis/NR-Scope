#include "nrscope/hdr/run_recorder.h"
#include <iostream>
#include <string>
#include <unistd.h>
#include <csignal>
#include <sched.h>
#include <thread>

#include "nrscope/hdr/nrscope_def.h"
#include "nrscope/hdr/load_config.h"
#include "nrscope/hdr/sensing/nrscope_sensing.h"

#include "srsran/common/band_helper.h"
#include "srsran/phy/common/phy_common_nr.h"

int main(int argc, char** argv){

  /* Line-buffer stdout: piped into tee it is block-buffered, and a run stopped by
    timeout or Ctrl-C lost its last few KB of output, often the very line that
    said why it stopped. */
  setvbuf(stdout, nullptr, _IOLBF, 0);

  // Initialise logging infrastructure
  srslog::init();

  // Usage: nrscope [config file], defaults to config.yaml in the working directory
  std::string file_name = (argc > 1) ? argv[1] : "config.yaml";
  if (access(file_name.c_str(), R_OK) != 0) {
    std::cout << "Cannot read config file: " << file_name << std::endl;
    return NR_FAILURE;
  }
  std::cout << "Using config: " << file_name << std::endl;

  int nof_usrp = get_nof_usrp(file_name);
  std::vector<Radio> radios(nof_usrp);

  /* Ctrl-C: stop the radio streams, then leave. The receive loops stop on nrscope_stop;
    once they have had time to leave UHD, the streams are stopped and flushed, so the
    X410 is idle for the next run, and the process exits without running the destructors
    of objects other threads still use. See my_sig_handler. */
  sem_init(&nrscope_stop_sem, 0, 0);
  {
    struct sigaction sa = {};
    sa.sa_handler = my_sig_handler;
    sigemptyset(&sa.sa_mask);
    sigaction(SIGINT, &sa, nullptr);
    sigaction(SIGTERM, &sa, nullptr);
    /* Run as "nrscope ... | tee log", Ctrl-C also ends tee, and the next print would kill
      the process with SIGPIPE before the streams are stopped. */
    signal(SIGPIPE, SIG_IGN);
  }
  std::thread([&radios]() {
    while (sem_wait(&nrscope_stop_sem) != 0) {
    }
    usleep(300000);
    for (auto& r : radios) {
      r.StopStreams();
    }
    // the estimates still queued for sensing.record_estimates, before _exit drops them
    nrscope_sensing_record_close();
    printf("Radio streams stopped, exiting\n");
    fflush(nullptr);
    _exit(0);
  }).detach();

  // TODO: Add a USRP as cell searcher -- always searching for the cell 
  if(load_config(radios, file_name) == NR_FAILURE){
    std::cout << "Load config fail." << std::endl;
    return NR_FAILURE;
  }

  /* Keep every thread that does not pin itself off the serial capture stages' CPUs.
    UHD starts its own threads unpinned during the radio init, inheriting this thread's
    mask; its control endpoint thread (uhd_ctrl_ep*, ~50% of a core with four chains)
    was seen on the fetch and consumer CPUs, taking time the capture needs. The stages
    and the workers pin themselves explicitly, so they are unaffected. */
  {
    cpu_set_t mask;
    if (sched_getaffinity(0, sizeof(mask), &mask) == 0) {
      int n_cut = 0;
      for (auto& r : radios) {
        if (!r.cpu_affinity)
          continue;
        for (int c : {r.fetch_cpu, r.consumer_cpu, r.dispatcher_cpu}) {
          if (c >= 0 && CPU_ISSET(c, &mask) && CPU_COUNT(&mask) > 1) {
            CPU_CLR(c, &mask);
            n_cut++;
          }
        }
      }
      if (n_cut > 0 && sched_setaffinity(0, sizeof(mask), &mask) != 0)
        perror("sched_setaffinity");
    }
  }

  // All the radios have the same setting for local log or push to google
  if(radios[0].local_log){
    std::vector<std::string> log_names(nof_usrp);
    for(int i = 0; i < nof_usrp; i++){
      log_names[i] = radios[i].log_name;
    }
    NRScopeLog::init_logger(log_names);
  }

  if(radios[0].to_google){
    ToGoogle::init_to_google(radios[0].google_credential, radios[0].google_dataset_id, nof_usrp);
  }

  std::vector<std::thread> radio_threads;

  for (auto& my_radio : radios) {
    radio_threads.emplace_back(&Radio::RadioThread, &my_radio);
  }

  for (auto& t : radio_threads) {
    if(t.joinable()){
      t.join();
    }
  }

  /* Stop the writer threads before returning. Without this the process aborts
    on the way out and whatever the radio was reporting is lost behind
    "terminate called without an active exception". */
  if(radios[0].to_google){
    ToGoogle::exit_to_google();
  }
  if(radios[0].local_log){
    NRScopeLog::exit_logger();
  }
  RunRecorder::close();

  return NR_SUCCESS;
}