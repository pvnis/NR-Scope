#include "nrscope/hdr/run_recorder.h"
#include <iostream>
#include <string>
#include <unistd.h>

#include "nrscope/hdr/nrscope_def.h"
#include "nrscope/hdr/load_config.h"

#include "srsran/common/band_helper.h"
#include "srsran/phy/common/phy_common_nr.h"

int main(int argc, char** argv){

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

  // TODO: Add a USRP as cell searcher -- always searching for the cell 
  if(load_config(radios, file_name) == NR_FAILURE){
    std::cout << "Load config fail." << std::endl;
    return NR_FAILURE;
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