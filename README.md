NR-Scope
========

**NR-Scope** is a 5G wireless network telemetry tool for 5G Standalone network measurement and performance optimization.  NR-Scope decodes DCI, SIB and RACH information from a 5G Standalone (New Radio) base station.  It is distinct from [NG-Scope](https://github.com/PrincetonUniversity/NG-Scope), our prior 4G cellular network telemetry tool.

NR-Scope is [published](https://doi.org/10.1145/3680121.3697808) in ACM CoNEXT 2024: if you use this tool, we ask that you cite our paper and reach out to the developers at the email address `ng-scope-dev` *at* `lists.cs.princeton.edu`

```
@inproceedings{10.1145/3680121.3697808,
  author = {Wan, Haoran and Cao, Xuyang and Marder, Alexander and Jamieson, Kyle},
  title = {NR-Scope: A Practical 5G Standalone Telemetry Tool},
  year = {2024},
  isbn = {9798400711084},
  publisher = {Association for Computing Machinery},
  address = {New York, NY, USA},
  url = {https://doi.org/10.1145/3680121.3697808},
  doi = {10.1145/3680121.3697808},
  booktitle = {Proceedings of the 20th International Conference on Emerging Networking EXperiments and Technologies},
  pages = {73–80},
  numpages = {8},
  keywords = {5G network, network measurement, telemetry, wireless network},
  location = {Los Angeles, CA, USA},
  series = {CoNEXT '24}
}
```

We welcome contributions.  If you use NR-Scope and would like to contribute improvements, modifications, bug fixes, and the like, please initiate a pull request on this repository and one of our team will review your contribution.

NR-Scope and NG-Scope are open source products of NSF Award CNS-2223556 (IMR: MT).  Any opinions, findings, and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Science Foundation.

# Features

The main features are as follows:

1. Supports 5G SA cells with 15kHz and 30kHz SCS, tested with 10MHz, 15 MHz, 20MHz, 40MHz and 100MHz cell bandwidth, achieving 1.4% DCI missing rate under the 100MHz T-Mobile TDD cell!
2. Carrier aggregation decoding, NR-Scope would use one SDR to try to decode both carrier aggregated and not carrier aggregated UEs within the RAN. Tested in commercial T-Mobile 15MHz FDD and 100MHz TDD cells.
3. 5G SA cell search and MIB decoding.
4. SIB 1 and all other SIBs decoding.
5. RRCSetup (MSG 4) decoding for continuously UE attatch detection.
6. Threaded DCI decoding for detected UEs.
7. SIB, RRCSetup and DCI decoding are threaded for performance and independent processing.
8. Can concurrently decode up to 64 UEs (max number we can create with groundtruth) in the same base station. In T-Mobile cell, we can decode up to hundreds of UEs.
9. Accurate DCIs, PRB and bit rate estimation for each UEs (please stay tuned for our paper).
10. Local logging functions in `.csv` file and remote google BigQuery uploading function in BigQuery table.
11. Support using multiple USRPs to decode multiple base stations independently, each has its own log.
12. Cell search for all possible frequencies in a separate file and the results are stored in `.csv` log file.
13. Non-initial (BWP id > 0) and plaintext (configured in SIB 1 or MSG 4) BWPs decoding in separate DCI decoder threads.
14. Threaded resampling function for better-fidelity TwinRX USRP X310 daughterboard (better signal quality but doesn't support 5G sampling rate without resampling).
15. Worker pool function is implemented. Each worker works on the slot data asynchronously. Now the processing doesn't need to keep up with the short slot time, and the throughput is increase through more workers.
16. Band list can be configured for cell scan functions.
17. A circular buffer is used for slot data dispatch among workers, now the CPU/worker number requirement is largely reduced.
18. PCI can be configured in the config.yaml. If configured, NR-Scope will proceed only if the detected cell's PCI matches what is configured.
19. Enable real-time constellation graph for PDCCH symbols, by setting `enable_gui: true` in the config file, now this function only support single USRP setting.
20. User provided RRCReconfiguration is supported, please check `rrc_recfg_config` folder for more information. We provide several examples of RRCReconfiguration for different scenarios (T-Mobile FDD and TDD SA cells).
21. GoldenSniffer blind DCI detection (`golden_sniffer:` in the config): every DCI of every UE, without having seen its RACH and without knowing its RNTI, from SIB1 on. A C port of [5G-GoldenSniffer](https://github.com/5G-GoldenSniffer/5G-GoldenSniffer) (MobiCom '26); see [GoldenSniffer](#goldensniffer) below.
22. Stay tuned... 😄

Please refer to the [wiki page](https://github.com/PrincetonUniversity/NG-Scope-5G/wiki) for more feature description and documentation.

## Notice about the base stations

If you are working with srsRAN base stations, please check out the `srsran` branch. The only difference is that during DCI decoding, the 64 QAM MCS table is selected, while the `main` branch uses 256 QAM MCS table, which is widely used in the commercial base stations. This MCS table information could be encrypted during the communication, thus it's set manually.

The `main` branch will decode the DCIs based on the RRCSetup decoded during the RACH, and will try to decode DCIs both with and without carrier aggregation bits. According to our experience with T-Mobile cells, it can decode the DCIs for most of the UEs in the RAN. In the `fixed_rnti` branch, we add another worker group to manually include the RRCReconfiguration information in their DCI decoders, because the RRCReconfiguration is encrypted and can't be decoded by NR-Scope. In this case, NR-Scope can decode nearly all of the UEs in the cells we tested with (~ 99%). Also, you can set a known RNTI in the config file in `fixed_rnti` branch, then this RNTI will be assigned to the worker group using the manually configuered RRCReconfiguration.

For a commercial cell, if you see a lot of RACH information but can't find further DCIs, this is because that the UEs are using updated configuration in the RRCReconfiguration. Please raise an issue and discuss with us, and if you could provide the RRCReconfiguration (from your phone), we can help you config NR-Scope properly.

## Requirements

We tested this system on Ubuntu 22.04 system and it may support other version of Ubuntu. To build this project and make it run properly, the following libraries are needed.

[UHD libraries](https://files.ettus.com/manual/page_install.html):

```
sudo add-apt-repository ppa:ettusresearch/uhd
sudo apt-get update
sudo apt-get install libuhd-dev uhd-host
```

[srsRAN_4G&#39;s requirements](https://docs.srsran.com/projects/4g/en/latest/general/source/1_installation.html):

```
sudo apt-get install build-essential cmake libfftw3-dev libmbedtls-dev libboost-program-options-dev libconfig++-dev libsctp-dev
```

[srsGUI](https://github.com/srsran/srsGUI):
```
git clone https://github.com/srsLTE/srsGUI.git
cd srsgui
mkdir build
cd build
cmake ../
make 
sudo make install
```

[yaml-cpp](https://github.com/jbeder/yaml-cpp):

```
# In a different directory
git clone https://github.com/jbeder/yaml-cpp.git
cd yaml-cpp
mkdir build
cd build
cmake ..
make
sudo make install
```

We need [liquid-dsp](https://github.com/jgaeddert/liquid-dsp) for resampling if better-fidelity TwinRX USRP X310 daughterboard is used:

```
# In a different directory
sudo apt-get install automake autoconf
# download source codes
git clone https://github.com/jgaeddert/liquid-dsp.git
cd liquid-dsp
# Building and installing the main library
./bootstrap.sh
./configure
make
sudo make install
sudo ldconfig
# to double check, libs should appear at /usr/local/lib and header liquid.h should appear at /usr/local/include/liquid/
```

For different USRP daughterboard, different `config.yaml` should be used. Please refer to the sample `config.yaml` in `./nrscope/config/config.yaml` and the explanatory comment in it.

CBX:

```
......
rf_args: "clock=external,type=x300,sampling_rate=23040000" #"type=x300" #"clock=external"
rx_gain: 30 # for x310, max rx gain is 31.5, for b210, it's around 80
srate_hz: 23040000 #11520000 #11520000 #23040000
srsran_srate_hz: 23040000
......
```

TwinRX (note TwinRX has a significantly higher rx gain limit):

```
......
rf_args: "clock=external,type=x300,master_clock_rate=200000000,sampling_rate=25000000" #"type=x300" #"clock=external"
rx_gain: 90 # for x310, max rx gain is 31.5, for b210, it's around 80
srate_hz: 25000000 #11520000 #11520000 #23040000
srsran_srate_hz: 23040000
......
```

Also, you should turn on the performance mode:

```
# This script is copied from srsRAN_Project
cd scripts
sudo ./srsran_performance
```

Push DCI logs to google storage BigQuery table (optional):
Using google client c++ library requires c++ 14.0, where in ubuntu 22.04, the default c++ version is 11.0. We thought that nobody wants to mess with the system compiling environments, so we implement the function of pushing data to google cloud storage with python. Here is a step-by-step instruction of how to push the DCI log to our google cloud storage:

```
1. sudo pip install google-cloud-storage geocoder
2. sudo gcloud auth application-default login --impersonate-service-account bigquery-writer@tutorial-explore.iam.gserviceaccount.com
3. # There will be a link from goole, open that link and login with any of your google account.
4. sudo gcloud init
5. # In the config file (./nrscope/config.yaml), set push_to_google: true.
6. # In the config file (./nrscope/config.yaml), set google_service_account_credential: "/home/wanhr/Downloads/nsf-2223556-222187-b5d2ea50f5d1.json" with the google service account credential file (provided by us)'s location on your file system.
7. # In the config file (./nrscope/config.yaml), set google_dataset_id: "ngscope5g_dci_log" with the google cloud Bigquery dataset name, everyone has his/her own dataset. The code will create one if the dataset with this dataset id is not existed.
8. sudo gcloud auth application-default set-quota-project <Your Google Storage Project Name>
```

## Files and functions

```
Entry: /nrscope/src/main.cc
Load config: /nrscope/src/libs/load_config.cc
Radio thread (cell search, mib decoding, coreset decoding, etc.): /nrscope/src/libs/radio_nr.cc
Config file: /nrscope/config/config.yaml
Cell scan entry (scan all GSCN/SSB points): /nrscope/src/scan_main.cc
```

## Usage

```
mkdir build
cd build
cmake ../
cmake -DCMAKE_CXX_FLAGS="-DBOOST_TIMER_ENABLE_DEPRECATED" .. 		# If you are using Ubuntu 24.04, constraint by the srsGUI library.
make all -j ${nof_proc}
cd nrscope/src/
sudo ./nrscope | tee ./$(date +"%Y-%m-%d_%H:%M:%S").txt | grep Found 	# This command can save the trace and show if cell or DCIs are found.

# or to scan all 5G SA cells (in nrscope/src/)
sudo ./nrscan
```

## GoldenSniffer

The regular DCI decoders only search the RNTIs learned from an RRCSetup (Msg4), with
the CORESET it configured. GoldenSniffer finds DCIs blind, as soon as SIB1 is decoded:

1. the busy subcarriers of the first PDCCH symbols are clustered into runs, and PDCCH
   candidates are cut from them (non-interleaved, and interleaved pairs/triples);
2. the QPSK modulation is removed by phase unwrapping, and the PDCCH DM-RS c_init, hence
   the scrambling ID N_ID, is read by inverting the Gold sequence (no ID search);
3. with the scrambling known, polar code + CRC + RNTI mask are affine over GF(2): one
   linear solve per DCI size gives the payload and the RNTI (no RNTI search).

Code: `nrscope/src/libs/golden_sniffer/` (C core, `nrscope/hdr/golden_sniffer/gs_core.h`,
MIT, ported from the MATLAB of the 5G-GoldenSniffer authors) and
`nrscope/src/libs/golden_sniffer.cc` (the worker thread, CSV, report). The encrypted
`PDCCH_DMRS_solve_underdet.p` of the MATLAB is replaced by `gs_dmrs_solve_underdet()`.
Differences from the MATLAB are listed at the top of `gs_search.c`.

```
golden_sniffer:
  enable: true
  max_duration: 2       # PDCCH symbols scanned (3 for cells with 3-symbol CORESETs)
  css_sizes: auto       # sizes allowed for n_RNTI = 0 DCIs that may be in CORESET#0
  min_parity: 16        # bits checked beyond the unknowns (false-detection guard)
  feed_rntis: false     # hand the C-RNTIs found to the regular DCI decoders
  csv: /tmp/golden_sniffer.csv
```

Every DCI goes to the CSV (SFN, slot, RNTI, N_ID, AL, size, payload); a per-RNTI count is
printed every `report_s` seconds and at exit.

Offline, on a recorded capture (srsRAN's `file` radio; NR-Scope then drops no slot and
stops at the end of the file), e.g. the GoldenSniffer dataset
(https://doi.org/10.5281/zenodo.21678295):

```
nrscope/scripts/sc16_to_fc32.py iq_1980M_23.04M_conf3_0.sc16 /tmp/conf3.fc32
cd build/nrscope/src && ./nrscope -c config_gs_file.yaml          # rx_file, ssb_freq inside
../../../nrscope/scripts/gs_compare.py gnb_conf3.log /tmp/golden_sniffer.csv --rnti 0x4601
```

`gs_compare.py` scores the CSV against the gNB log's `PDCCH:` lines. On the dataset:
conf3 (FDD, 15 kHz, 1 UE) 10939/10939 DCIs, conf17 (TDD, 30 kHz, 2 UEs, 2-symbol
CORESETs) 8290/8290 and 5223/5223, no false DCI. Self-test (no radio):
`build/nrscope/src/tests/gs_selftest` (also `ctest -R gs`).

### Fallback: UEs whose RRC config went stale

`mode: fallback` (default) runs the regular decoders as before and GoldenSniffer only
where they cannot follow a UE (its RRC reconfiguration is ciphered, or it attached before
NR-Scope started):

- **RRC**: decoded by the regular decoders from its RRCSetup.
- **ARMED**: a new RRCSetup or an RRC PDU on SRB1; the blind search runs for `arm_s`.
- **STALE**: seen blind (outside CORESET#0) where the regular decoders found nothing;
  its PDCCH is learned: N_ID, scrambling, duration, sizes, ALs, DL and UL.
- **CONVERGED**: decoded directly in the workers (`gs_decode_tracked_multi`): all
  converged UEs at once, over the whole band on their CCE grid, every AL; the DM-RS is
  correlated once per RB, each candidate is solved once per size and the solve gives the
  RNTI; candidates the hard decision cannot solve get a soft decision (`gs_ec.c`,
  5GDescrambler's error correction; `soft`). Other C-RNTIs found there are sightings.
- **DISCOVERED**: a C-RNTI found blind with no RRCSetup (startup window, sweeps, or the
  tracked decode's sightings), learned the same way.

A converged UE re-learns in a shadow when a DCI disagrees with its track or an RRC PDU
goes to it; a sweep starts when it goes quiet or one of its directions stops; it is kept
`keep_s` after its last DCI. In a replay the timers run on air time.

```
golden_sniffer:
  mode: fallback          # | always (search every slot: recorded captures without RACH)
  startup_s: 10           # discovery window after SIB1
  scan_period_s: 30       # periodic sweep, scan_window_s long
  converge_hits: 10       # DCIs learned before decoding a UE directly...
  converge_s: 3           # ...with DL and UL both seen, or this long
  keep_s: 300
  size_min: 20
  size_max: 140
  soft: true
  searcher_idle: true     # searchers at SCHED_IDLE: never in the pipeline's way
```

The final report gives the execution time of every stage (worker FFT, tracked decode,
search of empty and busy slots, queue wait, slot-to-DCI latency).

### Evaluation

- `build/nrscope/src/tests/gs_selftest [seed] [test]`: DM-RS inversion, GF(2), blind
  search, tracked decode (one UE, many UEs, vs SNR with and without soft decisions), soft
  error correction (detection, false alarms, timing).
- `nrscope/scripts/gs_eval_dataset.py --conf 1,3,6,...`: GoldenSniffer dataset captures
  (srsRAN gNB, FDD 15 kHz and TDD 30 kHz), blind search and fallback, scored per gNB log.
- `nrscope/scripts/gs_cellsim_run.sh run_dir nrscope/config/cellsim/c.ini`: synthetic
  cells with ground truth (`gs_cellsim`): many UEs, interleaved CORESETs, 3-symbol
  CORESETs, ScramblingID, idle periods, size changes, re-attaches, multipath;
  `gs_compare.py --truth` scores them.
- Live: `gs_compare.py gnb.log golden_sniffer.csv` against OCUDU's MAC "Slot decisions".

### Field test on commercial networks (sens6 + X410)

```
nrscope/scripts/field_prep.sh                       # X410 management IP, RT settings, checks
sudo nrscope/scripts/field_scan.py                  # cells around: SSB frequency and PCI
nrscope/scripts/field_run.sh --ssb 3624480000 --pci 123 --secs 120 --tag site1
#   -> prints "GoldenSniffer: carrier ... set rx_center_freq: F" if the capture does not
#      hold the whole carrier: run again with --center F
nrscope/scripts/field_run.sh --ssb 3624480000 --pci 123 --center F --secs 600 --tag site1
```

Each run directory (`~/field/<date>_<tag>`) keeps the config, the log, the DCI CSV, the
SRB PDUs and a summary. There is no gNB log to score against outdoors: a C-RNTI seen
many times with stable sizes is a UE; one-off RNTIs are false decodes (rare: 3 in 46 000
DCIs on the 100-UE synthetic cell). Keep the `[status]` RX peak below about -10 dBFS
(`--gain`). The X410's reference is about -4 ppm off (`freq_offset_ppm: -4.02`).

Please refer to this [wiki page](https://github.com/PrincetonUniversity/NR-Scope/wiki/1.-Introduction#output-format) for the description of log output and [FAQ](https://github.com/PrincetonUniversity/NR-Scope/wiki/7.-Usage-Recommendations-and-FAQ).
