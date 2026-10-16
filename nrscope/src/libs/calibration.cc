#include "nrscope/hdr/calibration.h"
#include "nrscope/hdr/nrscope_def.h"
#include "nrscope/hdr/sensing/nrscope_sensing.h"

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
#include <sys/stat.h>
#include <thread>
#include <unistd.h>

#ifdef NRSCOPE_HAVE_BLADERF
#include <libbladeRF.h>
#endif

CalibrationArgs nrscope_calibration_args;
bool            nrscope_calibrate_only = false;

static const std::string cal_dir      = std::string(NRSCOPE_ROOT_DIR) + "/calibration";
static const std::string cal_file     = cal_dir + "/calibration.yaml";
static const std::string history_file = cal_dir + "/history.csv";

void load_calibration_config(const YAML::Node& config_yaml)
{
  CalibrationArgs& a = nrscope_calibration_args;
  if (!config_yaml["calibration"]) {
    return;
  }
  const YAML::Node c = config_yaml["calibration"];
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
  if (c["min_coherence"])
    a.min_coherence = c["min_coherence"].as<double>();
  if (c["min_snr_db"])
    a.min_snr_db = c["min_snr_db"].as<double>();
}

/* nrscope runs under sudo: hand what it writes back to the user who ran it. */
static void give_to_invoking_user(const std::string& path)
{
  const char* uid = getenv("SUDO_UID");
  const char* gid = getenv("SUDO_GID");
  if (uid != nullptr && gid != nullptr && chown(path.c_str(), (uid_t)atoi(uid), (gid_t)atoi(gid)) != 0) {
    // not fatal: the file stays root's
  }
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
  const bladerf_channel       ch   = BLADERF_CHANNEL_TX(0);
  const char*                 dev  = a.bladerf_device.empty() ? nullptr : a.bladerf_device.c_str();
  int                         err  = bladerf_open(&tone_dev, dev);
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

/* One chain's result, all relative to chain 0. */
struct ChainCal {
  double tone_dbfs, snr_db, rel_db, expected_rel_db, phase_deg, expected_phase_deg, offset_deg, coherence;
  bool   ok;
};

static std::string yaml_list(const std::vector<ChainCal>& r, double ChainCal::*v, const char* fmt)
{
  std::string s = "[";
  for (size_t c = 0; c < r.size(); c++) {
    char one[32];
    snprintf(one, sizeof(one), fmt, r[c].*v);
    s += (c ? ", " : "") + std::string(one);
  }
  return s + "]";
}

bool run_tone_calibration(srsran::radio_interface_phy& radio,
                          uint32_t                     nof_chains,
                          double                       rx_center_hz,
                          double                       rx_gain_db,
                          double                       srate_hz,
                          uint32_t                     slot_sz)
{
  typedef std::complex<double> cd;
  const CalibrationArgs&       a     = nrscope_calibration_args;
  const double                 f_off = a.tone_freq_hz - rx_center_hz;
  if (a.tone_freq_hz <= 0 || std::fabs(f_off) > 0.45 * srate_hz) {
    std::cout << "Calibration: calibration.tone_freq_hz " << a.tone_freq_hz / 1e6
              << " MHz is not inside the capture (" << (rx_center_hz - 0.45 * srate_hz) / 1e6 << " - "
              << (rx_center_hz + 0.45 * srate_hz) / 1e6 << " MHz)" << std::endl;
    return false;
  }
  if (a.distances_mm.size() < nof_chains) {
    std::cout << "Calibration: calibration.distances_mm needs one distance per chain, " << nof_chains << std::endl;
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

  std::vector<double> tone_p(nof_chains, 0), back_p(nof_chains, 0);
  std::vector<cd>     cross(nof_chains, 0);
  int                 nof_tone = 0, nof_back = 0;
  double              best_df  = 0;
  bool                measured = start_tone(a);
  if (measured) {
    std::cout << "Calibration: bladeRF tone at " << a.tone_freq_hz / 1e6 << " MHz" << std::endl;

    // let the bladeRF settle, then find the tone near f_off: its oscillator is a few kHz off
    for (int i = 0; i < 200 && !nrscope_stop; i++) {
      receive();
    }
    double              best_p = 0;
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
      measured = false;
    }

    // the tone on every chain, then the background at the same frequency once it is off
    const int nof_slots = std::max(1, a.duration_ms * 2);
    for (int i = 0; measured && i < nof_slots && !nrscope_stop; i++) {
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
    for (int i = 0; measured && i < 100 && !nrscope_stop; i++) { // 50 ms for the tone to die away
      receive();
    }
    for (int i = 0; measured && i < nof_slots && !nrscope_stop; i++) {
      if (!receive()) {
        continue;
      }
      for (uint32_t c = 0; c < nof_chains; c++) {
        back_p[c] += std::norm(project(c));
      }
      nof_back++;
    }
    measured = measured && nof_tone > 0 && nof_back > 0 && !nrscope_stop;
  }
  stop_calibration_tone();
  for (uint32_t c = 0; c < nof_chains; c++) {
    free(bufs[c]);
  }
  if (!measured) {
    std::cout << "Calibration: not measured; " << cal_file << " is unchanged" << std::endl;
    return false;
  }

  const double          tone_hz = a.tone_freq_hz + best_df;
  const double          lambda  = 299792458.0 / tone_hz * 1e3; // mm
  const double          p0      = 10 * std::log10(tone_p[0] / nof_tone);
  std::vector<ChainCal> r(nof_chains);
  std::string           failed;
  for (uint32_t c = 0; c < nof_chains; c++) {
    ChainCal& k          = r[c];
    k.tone_dbfs          = 10 * std::log10(tone_p[c] / nof_tone);
    k.snr_db             = k.tone_dbfs - 10 * std::log10(back_p[c] / nof_back);
    k.rel_db             = k.tone_dbfs - p0;
    k.expected_rel_db    = 20 * std::log10(a.distances_mm[0] / a.distances_mm[c]);
    k.phase_deg          = std::arg(cross[c]) * 180.0 / M_PI;
    k.expected_phase_deg = wrap_deg(-360.0 * (a.distances_mm[c] - a.distances_mm[0]) / lambda);
    k.offset_deg         = wrap_deg(k.phase_deg - k.expected_phase_deg);
    k.coherence          = std::abs(cross[c]) / std::sqrt(tone_p[c] * tone_p[0]);
    k.ok                 = k.snr_db >= a.min_snr_db && k.coherence >= a.min_coherence;
    if (!k.ok) {
      char why[80];
      snprintf(why, sizeof(why), "%schain %u: SNR %.1f dB, coherence %.3f", failed.empty() ? "" : "; ", c, k.snr_db,
               k.coherence);
      failed += why;
    }
  }
  const bool passed = failed.empty();

  const std::time_t now = std::time(nullptr);
  char              when[32];
  std::strftime(when, sizeof(when), "%Y-%m-%dT%H:%M:%SZ", std::gmtime(&now));
  printf("Calibration at %.4f MHz (bladeRF %+.0f Hz off), %d slots, relative to chain 0:\n", tone_hz / 1e6, best_df,
         nof_tone);
  printf("  chain  tone dBFS  SNR dB   rel dB  expect   phase  expect  offset deg  coherence  check\n");
  for (uint32_t c = 0; c < nof_chains; c++) {
    const ChainCal& k = r[c];
    printf("  %5u  %9.1f  %6.1f  %7.2f  %6.2f  %6.1f  %6.1f  %10.1f  %9.3f  %s\n", c, k.tone_dbfs, k.snr_db, k.rel_db,
           k.expected_rel_db, k.phase_deg, k.expected_phase_deg, k.offset_deg, k.coherence, k.ok ? "ok" : "FAIL");
  }

  mkdir(cal_dir.c_str(), 0755); // already there is fine
  give_to_invoking_user(cal_dir);
  const bool new_history = access(history_file.c_str(), F_OK) != 0;
  {
    std::ofstream csv(history_file, std::ios::app);
    if (new_history) {
      csv << "time,passed,tone_hz,rx_center_hz,rx_gain_db,chain,tone_dbfs,snr_db,rel_db,expected_rel_db,phase_deg,"
             "expected_phase_deg,offset_deg,coherence\n";
    }
    for (uint32_t c = 0; c < nof_chains; c++) {
      const ChainCal& k = r[c];
      csv << when << "," << passed << "," << (long long)tone_hz << "," << (long long)rx_center_hz << "," << rx_gain_db
          << "," << c << "," << k.tone_dbfs << "," << k.snr_db << "," << k.rel_db << "," << k.expected_rel_db << ","
          << k.phase_deg << "," << k.expected_phase_deg << "," << k.offset_deg << "," << k.coherence << "\n";
    }
  }
  give_to_invoking_user(history_file);

  if (!passed) {
    printf("Calibration FAILED (needs SNR >= %.1f dB and coherence >= %.3f on every chain): %s.\n"
           "  %s is unchanged; the attempt is in %s\n",
           a.min_snr_db, a.min_coherence, failed.c_str(), cal_file.c_str(), history_file.c_str());
    return false;
  }
  {
    std::ofstream y(cal_file);
    y << "# written by nrscope --calibrate; every normal run at rx_center_hz uses aoa_cal_deg\n"
      << "time: " << when << "\n"
      << "unix_time: " << (long long)now << "\n"
      << "rx_center_hz: " << (long long)rx_center_hz << "\n"
      << "rx_gain_db: " << rx_gain_db << "\n"
      << "tone_hz: " << (long long)tone_hz << "\n"
      << "aoa_cal_deg: " << yaml_list(r, &ChainCal::offset_deg, "%.2f") << "\n"
      << "coherence: " << yaml_list(r, &ChainCal::coherence, "%.4f") << "\n"
      << "snr_db: " << yaml_list(r, &ChainCal::snr_db, "%.1f") << "\n"
      << "rel_db: " << yaml_list(r, &ChainCal::rel_db, "%.2f") << "\n";
  }
  give_to_invoking_user(cal_file);
  printf("Calibration PASSED: saved to %s; the next runs use aoa_cal_deg %s\n", cal_file.c_str(),
         yaml_list(r, &ChainCal::offset_deg, "%.1f").c_str());
  return true;
}

void apply_saved_calibration(double rx_center_hz, double rx_gain_db)
{
  if (access(cal_file.c_str(), F_OK) != 0) {
    printf("No saved calibration (%s): the AoA uses sensing.aoa_cal_deg. Run nrscope --calibrate to make one.\n",
           cal_file.c_str());
    return;
  }
  YAML::Node          y;
  double              center;
  std::vector<double> cal;
  long long           age_s;
  try {
    y      = YAML::LoadFile(cal_file);
    center = y["rx_center_hz"].as<double>();
    cal    = y["aoa_cal_deg"].as<std::vector<double> >();
    age_s  = (long long)std::time(nullptr) - y["unix_time"].as<long long>();
  } catch (const std::exception& e) {
    printf("Cannot read %s (%s): the AoA uses sensing.aoa_cal_deg\n", cal_file.c_str(), e.what());
    return;
  }
  if (std::fabs(center - rx_center_hz) > 1e3) {
    printf("Saved calibration (%s) was taken at %.3f MHz, not %.3f: the AoA uses sensing.aoa_cal_deg. "
           "Run nrscope --calibrate at this frequency.\n",
           y["time"].as<std::string>().c_str(), center / 1e6, rx_center_hz / 1e6);
    return;
  }
  nrscope_sensing_set_aoa_cal_deg(cal.data(), (int)cal.size());
  std::string deg;
  for (size_t c = 0; c < cal.size(); c++) {
    char one[16];
    snprintf(one, sizeof(one), "%s%.1f", c ? ", " : "", cal[c]);
    deg += one;
  }
  printf("AoA phase offsets [%s] deg from the calibration of %s (%lldh%02lldm ago)\n", deg.c_str(),
         y["time"].as<std::string>().c_str(), age_s / 3600, age_s / 60 % 60);
  if (std::fabs(y["rx_gain_db"].as<double>() - rx_gain_db) > 0.5) {
    // measured on the X410: 20 to 60 dB of gain moved the offsets by 2 to 5 degrees
    printf("  note: it was taken at %.1f dB of Rx gain, this run uses %.1f dB\n", y["rx_gain_db"].as<double>(),
           rx_gain_db);
  }
}
