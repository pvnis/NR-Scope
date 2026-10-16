#ifndef NRSCOPE_CALIBRATION_H
#define NRSCOPE_CALIBRATION_H

#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>

#include "srsran/interfaces/radio_interfaces.h"

/* Tone calibration of the receive chains: the calibration: block, or --calibrate.
  Once the radio is on its final frequency, a bladeRF sends a CW tone and nrscope
  logs, for every chain, the tone's power and its phase against chain 0, next to
  what the bladeRF antenna's distance to each antenna predicts. Measured minus
  predicted is the chain's own gain and phase offset, what NR_AOA_CAL_DEG holds.
  The result is printed and appended to a CSV; nothing uses it yet.

  The tone must sit inside the capture but outside the cell's carrier, e.g. in the
  guard band between the carrier's edge and the edge of the capture. The phases
  only repeat from run to run because the X410 tunes on whole seconds after a PPS
  time reset (rf_uhd_tune_rx_timed); any retune after the calibration starts a new
  set of offsets. */
struct CalibrationArgs {
  bool                enable          = false;
  double              tone_freq_hz    = 0;  // required
  std::string         bladerf_device  = ""; // libbladeRF device identifier, e.g. "*:serial=8a2"; empty: the first found
  int                 bladerf_tx_gain = 40; // dB
  int                 duration_ms     = 200;
  std::vector<double> distances_mm;         // bladeRF antenna to chain 0, 1, ...; empty: no prediction
  std::string         output          = "calibration.csv";
};
extern CalibrationArgs nrscope_calibration_args;

/* The calibration: block, into nrscope_calibration_args; load_config() calls it. */
void load_calibration_config(const YAML::Node& config_yaml);

/* Runs the calibration on nof_chains chains of radio, which is tuned to rx_center_hz,
  samples at srate_hz and is read slot_sz samples at a time. Returns false when it
  could not run (no bladeRF, tone outside the capture, ...); nrscope goes on anyway. */
bool run_tone_calibration(srsran::radio_interface_phy& radio,
                          uint32_t                     nof_chains,
                          double                       rx_center_hz,
                          double                       srate_hz,
                          uint32_t                     slot_sz);

/* Takes the tone off the air if a calibration is sending it. main's shutdown thread
  calls it, so a Ctrl-C during the calibration cannot leave the bladeRF transmitting. */
void stop_calibration_tone();

#endif
