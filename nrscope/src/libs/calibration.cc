#include "nrscope/hdr/calibration.h"
#include "nrscope/hdr/nrscope_def.h"

#include "srsran/phy/utils/vector.h"
#include "srsran/radio/rf_buffer.h"
#include "srsran/radio/rf_timestamp.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <complex>
#include <ctime>
#include <fstream>
#include <iostream>
#include <mutex>
#include <thread>
#include <unistd.h>

#ifdef NRSCOPE_HAVE_BLADERF
#include <libbladeRF.h>
#endif

CalibrationArgs nrscope_calibration_args;

void load_calibration_config(const YAML::Node& config_yaml)
{
  CalibrationArgs& a = nrscope_calibration_args;
  if (!config_yaml["calibration"]) {
    return;
  }
  const YAML::Node c = config_yaml["calibration"];
  if (c["enable"])
    a.enable = c["enable"].as<bool>();
  if (c["tone_freq_hz"])
    a.tone_freq_hz = c["tone_freq_hz"].as<double>();
  if (c["bladerf_device"])
    a.bladerf_device = c["bladerf_device"].as<std::string>();
  if (c["bladerf_tx_gain"])
    a.bladerf_tx_gain = c["bladerf_tx_gain"].as<int>();
  if (c["duration_ms"])
    a.duration_ms = c["duration_ms"].as<int>();
  if (c["distances_mm"])
    a.distances_mm = c["distances_mm"].as<std::vector<double> >();
  if (c["output"])
    a.output = c["output"].as<std::string>();
}

/* The tone: a constant sample on bladeRF TX1, i.e. a carrier exactly at its LO. */
#ifdef NRSCOPE_HAVE_BLADERF
static std::mutex        tone_mutex;
static struct bladerf*   tone_dev = nullptr;
static std::atomic<bool> tone_running{false};
static std::thread       tone_thread;

static bool start_tone(const CalibrationArgs& a)
{
  std::lock_guard<std::mutex> lock(tone_mutex);
  const bladerf_channel       ch  = BLADERF_CHANNEL_TX(0);
  const char*                 dev = a.bladerf_device.empty() ? nullptr : a.bladerf_device.c_str();
  int                         err = bladerf_open(&tone_dev, dev);
  const char*                 step = "open";
  if (err == 0 && (step = "sample rate", err = bladerf_set_sample_rate(tone_dev, ch, 1000000, nullptr)) == 0 &&
      (step = "bandwidth", err = bladerf_set_bandwidth(tone_dev, ch, 1500000, nullptr)) == 0 &&
      (step = "frequency", err = bladerf_set_frequency(tone_dev, ch, (bladerf_frequency)a.tone_freq_hz)) == 0 &&
      (step = "gain", err = bladerf_set_gain(tone_dev, ch, a.bladerf_tx_gain)) == 0 &&
      (step = "stream", err = bladerf_sync_config(tone_dev, BLADERF_TX_X1, BLADERF_FORMAT_SC16_Q11, 16, 8192, 8, 1000)) ==
          0 &&
      (step = "enable", err = bladerf_enable_module(tone_dev, ch, true)) == 0) {
    tone_running = true;
    tone_thread  = std::thread([]() {
      // half of full scale (2048) on I, so the DAC and the PA stay linear
      std::vector<int16_t> buf(2 * 8192);
      for (size_t i = 0; i < buf.size(); i += 2) {
        buf[i] = 1024;
      }
      while (tone_running && bladerf_sync_tx(tone_dev, buf.data(), 8192, nullptr, 1000) == 0) {
      }
    });
    return true;
  }
  std::cout << "Calibration: bladeRF " << step << " failed: " << bladerf_strerror(err) << std::endl;
  if (tone_dev != nullptr) {
    bladerf_close(tone_dev);
    tone_dev = nullptr;
  }
  return false;
}

void stop_calibration_tone()
{
  std::lock_guard<std::mutex> lock(tone_mutex);
  if (tone_dev == nullptr) {
    return;
  }
  tone_running = false;
  if (tone_thread.joinable()) {
    tone_thread.join();
  }
  bladerf_enable_module(tone_dev, BLADERF_CHANNEL_TX(0), false);
  bladerf_close(tone_dev);
  tone_dev = nullptr;
}
#else
static bool start_tone(const CalibrationArgs&)
{
  std::cout << "Calibration: nrscope was built without libbladeRF" << std::endl;
  return false;
}
void stop_calibration_tone() {}
#endif

static double wrap_deg(double d)
{
  d = std::fmod(d + 180.0, 360.0);
  return d < 0 ? d + 180.0 : d - 180.0;
}

bool run_tone_calibration(srsran::radio_interface_phy& radio,
                          uint32_t                     nof_chains,
                          double                       rx_center_hz,
                          double                       srate_hz,
                          uint32_t                     slot_sz)
{
  typedef std::complex<double> cd;
  const CalibrationArgs&       a     = nrscope_calibration_args;
  const double                 f_off = a.tone_freq_hz - rx_center_hz;
  if (a.tone_freq_hz <= 0 || std::fabs(f_off) > 0.45 * srate_hz) {
    std::cout << "Calibration: tone_freq_hz " << a.tone_freq_hz / 1e6 << " MHz is not inside the capture ("
              << (rx_center_hz - 0.45 * srate_hz) / 1e6 << " - " << (rx_center_hz + 0.45 * srate_hz) / 1e6
              << " MHz); skipped" << std::endl;
    return false;
  }

  std::vector<cf_t*>  bufs(nof_chains);
  srsran::rf_buffer_t rf_buffer = {};
  rf_buffer.set_nof_samples(slot_sz);
  for (uint32_t c = 0; c < nof_chains; c++) {
    bufs[c] = srsran_vec_cf_malloc(slot_sz);
    rf_buffer.set(c, bufs[c]);
  }
  auto receive = [&]() {
    srsran::rf_timestamp_t ts;
    return radio.rx_now(rf_buffer, ts);
  };
  auto chain = [&](uint32_t c) { return reinterpret_cast<const std::complex<float>*>(bufs[c]); };

  /* Mixing table for a tone f Hz from the centre, over one slot. Every slot is
    projected on its own, from its own first sample, and only ratios between chains
    of the same slot are kept, so a gap between two receives does not matter. */
  std::vector<cd> mix(slot_sz);
  auto            set_mix = [&](double f) {
    for (uint32_t n = 0; n < slot_sz; n++) {
      mix[n] = std::polar(1.0, -2.0 * M_PI * f * n / srate_hz);
    }
  };
  auto project = [&](uint32_t c) {
    const std::complex<float>* x = chain(c);
    cd                         s = 0;
    for (uint32_t n = 0; n < slot_sz; n++) {
      s += cd(x[n].real(), x[n].imag()) * mix[n];
    }
    return s / (double)slot_sz;
  };

  bool ok = start_tone(a);
  if (ok) {
    std::cout << "Calibration: bladeRF tone at " << a.tone_freq_hz / 1e6 << " MHz, " << f_off / 1e6
              << " MHz from the centre" << std::endl;

    // let the bladeRF settle, then find the tone near f_off: its oscillator is a few kHz off
    for (int i = 0; i < 200 && !nrscope_stop; i++) {
      receive();
    }
    double              best_df = 0, best_p = 0;
    std::vector<double> grid_p;
    if (receive()) {
      const uint32_t D = 64, K = slot_sz / D;
      set_mix(f_off);
      std::vector<std::vector<cd> > y(nof_chains, std::vector<cd>(K));
      for (uint32_t c = 0; c < nof_chains; c++) {
        const std::complex<float>* x = chain(c);
        for (uint32_t k = 0; k < K; k++) {
          for (uint32_t n = k * D; n < (k + 1) * D; n++) {
            y[c][k] += cd(x[n].real(), x[n].imag()) * mix[n];
          }
        }
      }
      for (double df = -20e3; df <= 20e3; df += 50) {
        double p = 0;
        for (uint32_t c = 0; c < nof_chains; c++) {
          cd s = 0;
          for (uint32_t k = 0; k < K; k++) {
            s += y[c][k] * std::polar(1.0, -2.0 * M_PI * df * (k * D + D / 2.0) / srate_hz);
          }
          p += std::norm(s);
        }
        grid_p.push_back(p);
        if (p > best_p) {
          best_p  = p;
          best_df = df;
        }
      }
    }
    set_mix(f_off + best_df);
    // a tone stands well clear of the median over the 40 kHz searched; noise alone reaches a few times it
    std::nth_element(grid_p.begin(), grid_p.begin() + grid_p.size() / 2, grid_p.end());
    if (grid_p.empty() || best_p < 10 * grid_p[grid_p.size() / 2]) {
      std::cout << "Calibration: no tone found within 20 kHz of " << a.tone_freq_hz / 1e6
                << " MHz; check the bladeRF antenna and bladerf_tx_gain" << std::endl;
      ok = false;
    }

    // the tone on every chain, then the background at the same frequency once it is off
    std::vector<double> tone_p(nof_chains, 0), back_p(nof_chains, 0);
    std::vector<cd>     cross(nof_chains, 0);
    const int           nof_slots = std::max(1, a.duration_ms * 2);
    int                 nof_tone = 0, nof_back = 0;
    for (int i = 0; ok && i < nof_slots && !nrscope_stop; i++) {
      if (!receive()) {
        continue;
      }
      const cd ref = project(0);
      for (uint32_t c = 0; c < nof_chains; c++) {
        const cd v = c == 0 ? ref : project(c);
        tone_p[c] += std::norm(v);
        cross[c] += v * std::conj(ref);
      }
      nof_tone++;
    }
    stop_calibration_tone();
    for (int i = 0; ok && i < 100 && !nrscope_stop; i++) { // 50 ms for the tone to die away
      receive();
    }
    for (int i = 0; ok && i < nof_slots && !nrscope_stop; i++) {
      if (!receive()) {
        continue;
      }
      for (uint32_t c = 0; c < nof_chains; c++) {
        back_p[c] += std::norm(project(c));
      }
      nof_back++;
    }
    ok = ok && nof_tone > 0 && nof_back > 0 && !nrscope_stop;

    if (ok) {
      const double tone_hz = a.tone_freq_hz + best_df;
      const double lambda  = 299792458.0 / tone_hz * 1e3; // mm
      const bool   geom    = a.distances_mm.size() >= nof_chains;
      if (!geom && !a.distances_mm.empty()) {
        std::cout << "Calibration: distances_mm has " << a.distances_mm.size() << " entries for " << nof_chains
                  << " chains; no prediction" << std::endl;
      }

      const bool    new_file = access(a.output.c_str(), F_OK) != 0;
      std::ofstream csv(a.output, std::ios::app);
      if (new_file) {
        csv << "time,tone_hz,rx_center_hz,chain,tone_dbfs,background_dbfs,snr_db,rel_db,expected_rel_db,"
               "gain_offset_db,phase_deg,expected_phase_deg,phase_offset_deg,coherence\n";
      }
      const std::time_t now = std::time(nullptr);
      char              when[32];
      std::strftime(when, sizeof(when), "%Y-%m-%dT%H:%M:%S", std::gmtime(&now));

      printf("Calibration at %.4f MHz (bladeRF %+.0f Hz off), %d slots, relative to chain 0:\n",
             tone_hz / 1e6,
             best_df,
             nof_tone);
      printf("  chain  tone dBFS  SNR dB   rel dB  expect  offset dB   phase  expect  offset deg  coherence\n");
      std::string offsets;
      const double p0 = 10 * std::log10(tone_p[0] / nof_tone);
      for (uint32_t c = 0; c < nof_chains; c++) {
        const double p     = 10 * std::log10(tone_p[c] / nof_tone);
        const double b     = 10 * std::log10(back_p[c] / nof_back);
        const double ph    = std::arg(cross[c]) * 180.0 / M_PI;
        const double coh   = std::abs(cross[c]) / std::sqrt(tone_p[c] * tone_p[0]);
        const double e_db  = geom ? 20 * std::log10(a.distances_mm[0] / a.distances_mm[c]) : NAN;
        const double e_deg = geom ? wrap_deg(-360.0 * (a.distances_mm[c] - a.distances_mm[0]) / lambda) : NAN;
        const double o_db  = p - p0 - e_db;
        const double o_deg = geom ? wrap_deg(ph - e_deg) : NAN;
        printf("  %5u  %9.1f  %6.1f  %7.2f  %6.2f  %9.2f  %6.1f  %6.1f  %10.1f  %9.3f\n",
               c, p, p - b, p - p0, e_db, o_db, ph, e_deg, o_deg, coh);
        csv << when << "," << (long long)tone_hz << "," << (long long)rx_center_hz << "," << c << "," << p << "," << b
            << "," << p - b << "," << p - p0 << "," << e_db << "," << o_db << "," << ph << "," << e_deg << ","
            << o_deg << "," << coh << "\n";
        char one[24];
        snprintf(one, sizeof(one), "%s%.1f", c ? ", " : "", geom ? o_deg : ph);
        offsets += one;
      }
      printf("  phase offsets (deg, %s): {%s}, appended to %s\n",
             geom ? "measured minus expected" : "measured, no distances",
             offsets.c_str(),
             a.output.c_str());
    }
  }
  stop_calibration_tone();

  for (uint32_t c = 0; c < nof_chains; c++) {
    free(bufs[c]);
  }
  return ok;
}
