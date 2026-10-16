#ifndef NRSCOPE_CALIBRATION_H
#define NRSCOPE_CALIBRATION_H

#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>

#include "srsran/interfaces/radio_interfaces.h"

/* Tone calibration of the receive chains.

  nrscope --calibrate tunes the radio, has a bladeRF send a CW tone, measures every
  chain's tone power and its phase against chain 0, saves the result and exits. Every
  normal run after that loads the saved calibration and the AoA removes its phase
  offsets, until the next --calibrate replaces it. This holds because the X410 tunes on
  whole seconds after a PPS time reset (rf_uhd_tune_rx_timed), so its chains come back
  to the same relative phases on every launch at the same frequency.

  The measurement: each chain's tone in every 0.5 ms slot as one complex number a_c.
  Its phase against chain 0, averaged over the slots, minus what the bladeRF antenna's
  distances to the antennas (distances_mm) predict, is the chain's own phase offset.
  Coherence, |sum a_c conj(a_0)| / sqrt(sum |a_c|^2 sum |a_0|^2), is how steady that
  phase stays from slot to slot: 1 for a clean tone, about 0.99 at 16 dB SNR, lower with
  interference, a moving path or a loose connector. The calibration passes when every
  chain reaches min_snr_db and min_coherence, and only a passing one is saved.

  Files, in <repo>/calibration/: calibration.yaml, the calibration in use, and
  history.csv, every attempt, passed or not. */
struct CalibrationArgs {
  double              tone_freq_hz    = 0;    // required for --calibrate
  std::string         bladerf_device  = "";   // libbladeRF device identifier, e.g. "*:serial=8a2"; empty: the first found
  int                 bladerf_tx_gain = 40;   // dB
  int                 duration_ms     = 200;
  std::vector<double> distances_mm;           // bladeRF antenna to chain 0, 1, ...; required to save
  double              min_coherence   = 0.95; // every chain, to pass
  double              min_snr_db      = 10;   // every chain, to pass
};
extern CalibrationArgs nrscope_calibration_args;

/* --calibrate: measure, save if it passes, and end the run. */
extern bool nrscope_calibrate_only;

/* The calibration: block, into nrscope_calibration_args; load_config() calls it. */
void load_calibration_config(const YAML::Node& config_yaml);

/* --calibrate: measures nof_chains chains of radio, which is tuned to rx_center_hz at
  rx_gain_db, samples at srate_hz and is read slot_sz samples at a time; saves the
  result if it passes. Returns whether it passed. */
bool run_tone_calibration(srsran::radio_interface_phy& radio,
                          uint32_t                     nof_chains,
                          double                       rx_center_hz,
                          double                       rx_gain_db,
                          double                       srate_hz,
                          uint32_t                     slot_sz);

/* A normal run: gives the AoA the saved calibration's phase offsets if it was taken at
  rx_center_hz, and says which offsets are in use. */
void apply_saved_calibration(double rx_center_hz, double rx_gain_db);

/* Takes the tone off the air if a calibration is sending it. main's shutdown thread
  calls it, so a Ctrl-C during the calibration cannot leave the bladeRF transmitting. */
void stop_calibration_tone();

#endif
