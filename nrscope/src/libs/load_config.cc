#include <stdlib.h>
#include <string.h>

#include "nrscope/hdr/sensing/nrscope_sensing.h"
#include "srsran/common/band_helper.h"
#include "srsran/common/crash_handler.h"
#include "srsran/common/string_helpers.h"
#include "srsran/phy/common/phy_common_nr.h"
#include "srsran/radio/radio.h"
#include "srsran/srsran.h"
#include "srsue/hdr/phy/phy_nr_sa.h"
#include "srsue/hdr/stack/ue_stack_nr.h"
#include "test/phy/dummy_ue_stack.h"
#include <boost/program_options.hpp>
#include <boost/program_options/parsers.hpp>

#include "nrscope/hdr/load_config.h"
#include "nrscope/hdr/run_recorder.h"
#include "nrscope/hdr/nrscope_def.h"
#include "nrscope/hdr/radio_nr.h"
#include "nrscope/hdr/rrc_recfg_parse.h"

#include "srsran/phy/ue/srsgui_plot.h"

using namespace std;

int get_nof_usrp(std::string file_name)
{
  YAML::Node config_yaml = YAML::LoadFile(file_name);
  if (config_yaml["nof_usrp_dev"]) {
    return config_yaml["nof_usrp_dev"].as<int>();
  } else {
    return -1;
  }
}

int load_config(std::vector<Radio>& radios, std::string file_name)
{
  /*To get the config from file, including the number of usrp devices,
  cell central frequency and etc.*/
  YAML::Node config_yaml = YAML::LoadFile(file_name);

  std::cout << "Reading configs." << std::endl;

  int nof_usrp = config_yaml["nof_usrp_dev"].as<int>();

  for (int i = 0; i < nof_usrp; i++) {
    radios[i].rf_index       = i;
    std::string setting_name = "usrp_setting_" + to_string(i);
    std::cout << "USRP Device: " << i << std::endl;
    if (config_yaml[setting_name]) {
      srsran::rf_args_t rf_args = {};
      if (config_yaml[setting_name]["rf_args"]) {
        std::string rf_args_config    = config_yaml[setting_name]["rf_args"].as<string>();
        radios[i].rf_args.device_args = rf_args_config;
      }
      std::cout << "    rf_args: " << radios[i].rf_args.device_args << std::endl;

      if (config_yaml[setting_name]["device_name"]) {
        std::string rf_args_devicename = config_yaml[setting_name]["device_name"].as<string>();
        radios[i].rf_args.device_name  = rf_args_devicename;
      }
      std::cout << "    device_name: " << radios[i].rf_args.device_name << endl;

      if (config_yaml[setting_name]["log_level"]) {
        std::string rf_args_loglevel = config_yaml[setting_name]["log_level"].as<string>();
        radios[i].rf_args.log_level  = rf_args_loglevel;
      }

      if (config_yaml[setting_name]["srsran_srate_hz"]) {
        radios[i].rf_args.srsran_srate_hz = config_yaml[setting_name]["srsran_srate_hz"].as<double>();
      }
      std::cout << "    srsran_srate_hz: " << radios[i].rf_args.srsran_srate_hz / 1e6 << " MHz" << std::endl;

      if (config_yaml[setting_name]["srate_hz"]) {
        radios[i].rf_args.srate_hz = config_yaml[setting_name]["srate_hz"].as<double>();
      }
      std::cout << "    srate_hz: " << radios[i].rf_args.srate_hz / 1e6 << " MHz" << std::endl;

      if (config_yaml[setting_name]["rx_gain"]) {
        radios[i].rf_args.rx_gain = config_yaml[setting_name]["rx_gain"].as<float>();
      }
      std::cout << "    rx_gain: " << radios[i].rf_args.rx_gain << std::endl;

      if (config_yaml[setting_name]["min_rx_gain"]) {
        radios[i].min_rx_gain = config_yaml[setting_name]["min_rx_gain"].as<float>();
      }
      std::cout << "    AGC min_rx_gain: " << radios[i].min_rx_gain << std::endl;

      if (config_yaml[setting_name]["band_list"]) {
        radios[i].band_list = config_yaml[setting_name]["band_list"].as<std::vector<uint16_t> >();
        std::cout << "    band_list: ";
        for (const auto& b : radios[i].band_list)
          std::cout << b << " ";
        std::cout << std::endl;
      }

      if (config_yaml[setting_name]["pci"]) {
        radios[i].pci = config_yaml[setting_name]["pci"].as<uint16_t>();
        if (radios[i].pci > 1007) {
          ERROR("PCI ranges from 0 to 1007, are you doing it right?");
        }
        std::cout << "    pci: " << radios[i].pci << std::endl;
      } else {
        radios[i].pci = 9999; // PCI ranges from 0 to 1007, thus 9999 is invalid
        // std::cout << "    pci: " << radios[i].pci << std::endl;
      }

      if (config_yaml[setting_name]["max_rx_gain"]) {
        radios[i].max_rx_gain = config_yaml[setting_name]["max_rx_gain"].as<float>();
      }
      std::cout << "    AGC max_rx_gain: " << radios[i].max_rx_gain << std::endl;

      if (config_yaml[setting_name]["nof_carriers"]) {
        radios[i].rf_args.nof_carriers = config_yaml[setting_name]["nof_carriers"].as<int>();
      }
      std::cout << "    nof_carriers: " << radios[i].rf_args.nof_carriers << std::endl;

      if (config_yaml[setting_name]["nof_antennas"]) {
        radios[i].rf_args.nof_antennas = config_yaml[setting_name]["nof_antennas"].as<int>();
      }
      std::cout << "    nof_antennas: " << radios[i].rf_args.nof_antennas << std::endl;

      if (config_yaml[setting_name]["freq_offset"]) {
        radios[i].rf_args.freq_offset = config_yaml[setting_name]["freq_offset"].as<float>();
      }
      std::cout << "    freq_offset: " << radios[i].rf_args.freq_offset << std::endl;

      if (config_yaml[setting_name]["scs_index"]) {
        radios[i].ssb_scs = (srsran_subcarrier_spacing_t)config_yaml[setting_name]["scs_index"].as<int>();
      }
      std::cout << "    scs: " << radios[i].ssb_scs << std::endl;

      /* ssb_freq is where the SS/PBCH block is; rx_center_freq is where the radio
        tunes. They default to the same value, which only captures the whole
        carrier when the SSB sits near its centre. dl_center_frequency_hz holds
        the tuning frequency from here on: every offset into the FFT grid, and
        the per-symbol phase compensation, is relative to it. */
      if (config_yaml[setting_name]["ssb_freq"]) {
        radios[i].args_t.base_carrier.ssb_center_freq_hz     = config_yaml[setting_name]["ssb_freq"].as<double>();
        radios[i].args_t.base_carrier.dl_center_frequency_hz = radios[i].args_t.base_carrier.ssb_center_freq_hz;
      }
      if (config_yaml[setting_name]["rx_center_freq"]) {
        radios[i].args_t.base_carrier.dl_center_frequency_hz = config_yaml[setting_name]["rx_center_freq"].as<double>();
      }
      std::cout << "    ssb_freq: " << radios[i].args_t.base_carrier.ssb_center_freq_hz / 1e6
                << " MHz, rx_center_freq: " << radios[i].args_t.base_carrier.dl_center_frequency_hz / 1e6 << " MHz"
                << std::endl;

      if (config_yaml[setting_name]["rf_log_level"]) {
        radios[i].rf_args.log_level = config_yaml[setting_name]["rf_log_level"].as<string>();
      } else {
        radios[i].rf_args.log_level = "info";
      }

      if (config_yaml[setting_name]["log_name"]) {
        radios[i].log_name = config_yaml[setting_name]["log_name"].as<string>();
      }

      if (config_yaml[setting_name]["google_dataset_id"]) {
        radios[i].google_dataset_id = config_yaml[setting_name]["google_dataset_id"].as<string>();
      }

      if (config_yaml[setting_name]["nof_rnti_worker_groups"]) {
        radios[i].nof_rnti_worker_groups = config_yaml[setting_name]["nof_rnti_worker_groups"].as<int>();
      } else {
        radios[i].nof_rnti_worker_groups = 1;
      }

      radios[i].nof_threads = radios[i].nof_rnti_worker_groups;

      if (config_yaml[setting_name]["nof_bwps"]) {
        radios[i].nof_bwps = config_yaml[setting_name]["nof_bwps"].as<int>();
      } else {
        radios[i].nof_bwps = 1;
      }

      if (config_yaml[setting_name]["cpu_affinity"]) {
        radios[i].cpu_affinity = config_yaml[setting_name]["cpu_affinity"].as<bool>();
      } else {
        radios[i].cpu_affinity = false;
      }
      for (auto& [key, cpu] : {std::pair<const char*, int*>{"fetch_cpu", &radios[i].fetch_cpu},
                               {"consumer_cpu", &radios[i].consumer_cpu},
                               {"dispatcher_cpu", &radios[i].dispatcher_cpu},
                               {"results_cpu", &radios[i].results_cpu}}) {
        if (config_yaml[setting_name][key]) {
          *cpu = config_yaml[setting_name][key].as<int>();
        }
      }
      if (config_yaml[setting_name]["worker_cpus"]) {
        radios[i].worker_cpus = config_yaml[setting_name]["worker_cpus"].as<std::vector<int> >();
      }
      if (radios[i].cpu_affinity) {
        std::cout << "    cpu layout: fetch " << radios[i].fetch_cpu << ", consumer " << radios[i].consumer_cpu
                  << ", dispatcher " << radios[i].dispatcher_cpu << ", results " << radios[i].results_cpu << ", "
                  << radios[i].worker_cpus.size() << " worker CPUs" << std::endl;
      }

      if (config_yaml[setting_name]["disable_cfo"]) {
        radios[i].disable_cfo = config_yaml[setting_name]["disable_cfo"].as<bool>();
      } else {
        radios[i].disable_cfo = false;
      }
      std::cout << "    disable_cfo: " << (radios[i].disable_cfo ? "true" : "false") << std::endl;

      if (config_yaml[setting_name]["rnti_expiry_s"]) {
        radios[i].rnti_expiry_s = config_yaml[setting_name]["rnti_expiry_s"].as<double>();
      }
      std::cout << "    rnti_expiry_s: " << radios[i].rnti_expiry_s << std::endl;

      if (config_yaml[setting_name]["agc"]) {
        radios[i].agc = config_yaml[setting_name]["agc"].as<bool>();
      } else {
        radios[i].agc = true;
      }
      std::cout << "    agc: " << (radios[i].agc ? "true" : "false") << std::endl;

      if (config_yaml[setting_name]["nof_workers"]) {
        radios[i].nof_workers = config_yaml[setting_name]["nof_workers"].as<int>();
        if (radios[i].nof_workers > 128) {
          ERROR("Worker number shouldn't be > 128");
          return SRSRAN_ERROR;
        }
      } else {
        radios[i].nof_workers = 1;
      }

      if (config_yaml[setting_name]["rrc_recfg_config"]) {
        parse_rrc_recfg(radios[i], config_yaml[setting_name]["rrc_recfg_config"].as<string>());
      }

      radios[i].nof_threads = radios[i].nof_threads * radios[i].nof_bwps;

      // std::cout << "    nof_thread: " << radios[i].nof_thread << std::endl;
    } else {
      std::cout << "Please set the usrp_setting_" << i << " in config.yaml properly." << std::endl;
      return NR_FAILURE;
    }
  }

  /* Check that every pinned CPU exists. Workers share one pool of CPUs, and the
    threads each one spawns per slot inherit it, so there is no per-thread core
    count to satisfy any more. */
  const int nof_cores = (int)std::thread::hardware_concurrency();
  for (int i = 0; i < nof_usrp; i++) {
    if (!radios[i].cpu_affinity) {
      continue;
    }
    std::vector<int> all = radios[i].worker_cpus;
    all.insert(all.end(), {radios[i].fetch_cpu, radios[i].consumer_cpu, radios[i].dispatcher_cpu, radios[i].results_cpu});
    for (int cpu : all) {
      if (cpu >= nof_cores) {
        ERROR("CPU affinity set, but CPU %d is not on this machine (%d CPUs)", cpu, nof_cores);
        return NR_FAILURE;
      }
    }
  }

  std::string setting_name = "log_config";
  if (config_yaml[setting_name]["local_log"]) {
    for (int i = 0; i < nof_usrp; i++) {
      radios[i].local_log = config_yaml[setting_name]["local_log"].as<bool>();
    }

    if (config_yaml[setting_name]["enable_gui"]) {
      if (config_yaml[setting_name]["enable_gui"].as<bool>()) {
        init_plots();
      }
    }
  } else {
    for (int i = 0; i < nof_usrp; i++) {
      radios[i].local_log = false;
    }
  }

  /* Record RRCSetups and UE DCIs to per-run CSVs at the project root, and keep
    the terminal to the "Found DCI" and "hooray" lines. See run_recorder.h. */
  RunRecorder::init(config_yaml[setting_name]["recording_mode"] &&
                    config_yaml[setting_name]["recording_mode"].as<bool>());
  /* Keep the terminal clean without writing CSVs: recording_mode's quiet terminal for a
  run that does not want its files (e.g. a sensing-only or sample-recording run). */
  RunRecorder::set_quiet(config_yaml[setting_name]["quiet"] &&
                         config_yaml[setting_name]["quiet"].as<bool>());
  RunRecorder::enable_pdcch_candidates(config_yaml[setting_name]["record_pdcch_candidates"] &&
                                       config_yaml[setting_name]["record_pdcch_candidates"].as<bool>());

  /* sensing: block, OAI's --sensing-* options under the same names without the
  prefix (dashes as underscores). Off unless enabled; see nrscope_sensing.h. */
  nrscope_sensing_default_args(&nrscope_sensing_args);
  snprintf(nrscope_sensing_args.dump, sizeof(nrscope_sensing_args.dump), "%s/sensing/map2d.csv", NRSCOPE_ROOT_DIR);
  if (config_yaml["sensing"]) {
    const YAML::Node sn = config_yaml["sensing"];
    auto flag = [&sn](const char* key, bool* v) {
      if (sn[key]) {
        *v = sn[key].as<bool>();
      }
    };
    flag("enable", &nrscope_sensing_args.enable);
    if (sn["symbols"]) {
      nrscope_sensing_args.symbols = sn["symbols"].as<int>();
      // the upper end, NR_SENSING_HISTORY_DEPTH, is checked where the histories are made
      if (nrscope_sensing_args.symbols < 32) { // NR_SENSING_MIN_SNAPSHOTS
        std::cerr << "sensing: symbols must be at least 32, got " << nrscope_sensing_args.symbols << std::endl;
        exit(EXIT_FAILURE);
      }
    }
    if (sn["max_speed"]) {
      nrscope_sensing_args.max_speed_ms = sn["max_speed"].as<double>();
    }
    if (sn["dump"]) {
      snprintf(nrscope_sensing_args.dump, sizeof(nrscope_sensing_args.dump), "%s", sn["dump"].as<std::string>().c_str());
    }
    flag("clutter_removal", &nrscope_sensing_args.clutter_removal);
    flag("clutter_kernel", &nrscope_sensing_args.clutter_kernel);
    flag("clutter_compare", &nrscope_sensing_args.clutter_compare);
    flag("antenna_avg", &nrscope_sensing_args.antenna_avg);
    flag("layer_avg", &nrscope_sensing_args.layer_avg);
    flag("random_drop", &nrscope_sensing_args.random_drop);
    flag("tdd_detect", &nrscope_sensing_args.tdd_detect);
    flag("music", &nrscope_sensing_args.music);
    flag("compensate_window_shifts", &nrscope_sensing_args.compensate_window_shifts);
    flag("spatial_null", &nrscope_sensing_args.spatial_null);
    flag("mirror_reject", &nrscope_sensing_args.mirror_reject);

    /* Parameter sweep (mode b). sweep.param names the clutter field to vary, lo/hi/step
    its range; empty or absent disables it. See nrscope_sensing_args_t.sweep_param. */
    if (sn["sweep"]) {
      const YAML::Node sw = sn["sweep"];
      if (sw["param"]) {
        snprintf(nrscope_sensing_args.sweep_param, sizeof(nrscope_sensing_args.sweep_param), "%s",
                 sw["param"].as<std::string>().c_str());
      }
      if (sw["lo"])   nrscope_sensing_args.sweep_lo   = sw["lo"].as<double>();
      if (sw["hi"])   nrscope_sensing_args.sweep_hi   = sw["hi"].as<double>();
      if (sw["step"]) nrscope_sensing_args.sweep_step = sw["step"].as<double>();
    }

    flag("test_record_samples", &nrscope_sensing_args.test_record_samples);
    if (sn["record_max_files"])
      nrscope_sensing_args.record_max_files = sn["record_max_files"].as<int>();
    if (sn["avg_maps"]) {
      nrscope_sensing_args.avg_maps = sn["avg_maps"].as<int>();
      if (nrscope_sensing_args.avg_maps < 1 || nrscope_sensing_args.avg_maps > NR_SENSING_AVG_MAX) {
        std::cerr << "sensing: avg_maps must be 1.." << NR_SENSING_AVG_MAX << ", got "
                  << nrscope_sensing_args.avg_maps << std::endl;
        exit(EXIT_FAILURE);
      }
    }
    flag("avg_vcomp", &nrscope_sensing_args.avg_vcomp);
    if (sn["map_max_gap_ms"])
      nrscope_sensing_args.map_max_gap_ms = sn["map_max_gap_ms"].as<double>();
  }

  if (config_yaml[setting_name]["push_to_google"]) {
    for (int i = 0; i < nof_usrp; i++) {
      radios[i].to_google = config_yaml[setting_name]["push_to_google"].as<bool>();
      if (config_yaml[setting_name]["google_service_account_credential"]) {
        radios[i].google_credential = config_yaml[setting_name]["google_service_account_credential"].as<string>();
      }
    }
    // if(config_yaml[setting_name]["google_project_id"]){
    //   radios[i].google_project_id =
    //      config_yaml[setting_name]["google_project_id"].as<string>();
    // }
  }

  return NR_SUCCESS;
}
