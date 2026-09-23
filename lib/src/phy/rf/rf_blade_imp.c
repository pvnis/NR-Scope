/**
 * Copyright 2013-2023 Software Radio Systems Limited
 *
 * This file is part of srsRAN.
 *
 * srsRAN is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * srsRAN is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * A copy of the GNU Affero General Public License can be found in
 * the LICENSE file in the top-level directory of this distribution
 * and at http://www.gnu.org/licenses/.
 *
 */

#include <libbladeRF.h>
#include <string.h>
#include <unistd.h>

#include "rf_blade_imp.h"
#include "rf_plugin.h"
#include "srsran/phy/common/timestamp.h"
#include "srsran/phy/utils/debug.h"
#include "srsran/phy/utils/vector.h"

#define UNUSED __attribute__((unused))
#define CONVERT_BUFFER_SIZE (240 * 1024)

/* Receive channels this driver can stream at once. The bladeRF 2.0 micro (xA4,
 * xA9) carries an AD9361 with two RX chains sharing one LO and one sample
 * clock, which is what makes them phase coherent and therefore usable for
 * direction finding. The original bladeRF x40/x115 has a single chain and must
 * stay at 1; asking for 2 there fails when the stream is configured.
 *
 * Note that BLADERF_RX_X2 == 2 == BLADERF_CHANNEL_RX(1). The layout enum and
 * the channel macro overlap numerically, so a call meant for "both channels"
 * that is handed the layout silently configures channel 1 alone. Every
 * per-channel call below therefore uses BLADERF_CHANNEL_RX(i) explicitly, and
 * the layout appears only in bladerf_sync_config(). */
#define BLADE_MAX_RX_CHANNELS 2

typedef struct {
  struct bladerf*     dev;
  bladerf_sample_rate rx_rate;
  bladerf_sample_rate tx_rate;
  /* In a 2 channel layout bladerf_sync_rx() returns the chains interleaved
   * sample by sample, so this holds nof_rx_channels * nsamples of them. */
  int16_t          rx_buffer[CONVERT_BUFFER_SIZE];
  int16_t          tx_buffer[CONVERT_BUFFER_SIZE];
  bool             rx_stream_enabled;
  bool             tx_stream_enabled;
  uint32_t         nof_rx_channels;
  srsran_rf_info_t info;
} rf_blade_handler_t;

/// Layout to configure the RX stream with, given the number of chains in use
static inline bladerf_channel_layout blade_rx_layout(uint32_t nof_rx_channels)
{
  return (nof_rx_channels > 1) ? BLADERF_RX_X2 : BLADERF_RX_X1;
}

static srsran_rf_error_handler_t blade_error_handler     = NULL;
static void*                     blade_error_handler_arg = NULL;

void rf_blade_suppress_stdout(UNUSED void* h)
{
  bladerf_log_set_verbosity(BLADERF_LOG_LEVEL_SILENT);
}

void rf_blade_register_error_handler(UNUSED void* ptr, srsran_rf_error_handler_t new_handler, void* arg)
{
  blade_error_handler     = new_handler;
  blade_error_handler_arg = arg;
}

const unsigned int num_buffers       = 256;
const unsigned int ms_buffer_size_rx = 1024;
const unsigned int buffer_size_tx    = 1024;
const unsigned int num_transfers     = 32;
const unsigned int timeout_ms        = 4000;

const char* rf_blade_devname(UNUSED void* h)
{
  return DEVNAME;
}

int rf_blade_start_tx_stream(void* h)
{
  int                 status;
  rf_blade_handler_t* handler = (rf_blade_handler_t*)h;

  status = bladerf_sync_config(handler->dev,
                               BLADERF_TX_X1,
                               BLADERF_FORMAT_SC16_Q11_META,
                               num_buffers,
                               buffer_size_tx,
                               num_transfers,
                               timeout_ms);
  if (status != 0) {
    ERROR("Failed to configure TX sync interface: %s", bladerf_strerror(status));
    return status;
  }
  status = bladerf_enable_module(handler->dev, BLADERF_TX_X1, true);
  if (status != 0) {
    ERROR("Failed to enable TX module: %s", bladerf_strerror(status));
    return status;
  }
  handler->tx_stream_enabled = true;
  return 0;
}

int rf_blade_start_rx_stream(void* h, UNUSED bool now)
{
  int                 status;
  rf_blade_handler_t* handler = (rf_blade_handler_t*)h;

  /* Configure the device's RX module for use with the sync interface.
   * SC16 Q11 samples *with* metadata are used.
   *
   * buffer_size counts samples across every chain, not per chain: libbladeRF
   * documents 2048 samples on two channels as generating 4096 total. So it is
   * scaled by the chain count, otherwise a 2 chain stream would hold half the
   * time depth of a 1 chain one and the host would get half the slack before
   * the device FIFO overruns. The result stays a multiple of 1024, which
   * bladerf_sync_config() requires, because the single chain value already is. */
  uint32_t buffer_size_rx = ms_buffer_size_rx * (handler->rx_rate / 1000 / 1024) * handler->nof_rx_channels;

  status = bladerf_sync_config(handler->dev,
                               blade_rx_layout(handler->nof_rx_channels),
                               BLADERF_FORMAT_SC16_Q11_META,
                               num_buffers,
                               buffer_size_rx,
                               num_transfers,
                               timeout_ms);
  if (status != 0) {
    ERROR("Failed to configure RX sync interface for %d channel(s): %s",
          handler->nof_rx_channels,
          bladerf_strerror(status));
    return status;
  }
  status = bladerf_sync_config(handler->dev,
                               BLADERF_TX_X1,
                               BLADERF_FORMAT_SC16_Q11_META,
                               num_buffers,
                               buffer_size_tx,
                               num_transfers,
                               timeout_ms);
  if (status != 0) {
    ERROR("Failed to configure TX sync interface: %s", bladerf_strerror(status));
    return status;
  }
  /* Enabled one chain at a time: bladerf_enable_module() takes a channel, not a
  layout, so enabling "X2" would only ever turn on channel 1. */
  for (uint32_t i = 0; i < handler->nof_rx_channels; i++) {
    status = bladerf_enable_module(handler->dev, BLADERF_CHANNEL_RX(i), true);
    if (status != 0) {
      ERROR("Failed to enable RX channel %d: %s", i, bladerf_strerror(status));
      return status;
    }
  }
  status = bladerf_enable_module(handler->dev, BLADERF_TX_X1, true);
  if (status != 0) {
    ERROR("Failed to enable TX module: %s", bladerf_strerror(status));
    return status;
  }
  handler->rx_stream_enabled = true;
  return 0;
}

int rf_blade_stop_rx_stream(void* h)
{
  rf_blade_handler_t* handler = (rf_blade_handler_t*)h;
  int                 status  = 0;
  for (uint32_t i = 0; i < handler->nof_rx_channels; i++) {
    int s = bladerf_enable_module(handler->dev, BLADERF_CHANNEL_RX(i), false);
    if (s != 0) {
      status = s;
    }
  }
  if (status != 0) {
    ERROR("Failed to disable RX channels: %s", bladerf_strerror(status));
    return status;
  }
  status = bladerf_enable_module(handler->dev, BLADERF_TX_X1, false);
  if (status != 0) {
    ERROR("Failed to enable TX module: %s", bladerf_strerror(status));
    return status;
  }
  handler->rx_stream_enabled = false;
  handler->tx_stream_enabled = false;
  return 0;
}

void rf_blade_flush_buffer(UNUSED void* h) {}

bool rf_blade_has_rssi(UNUSED void* h)
{
  return false;
}

float rf_blade_get_rssi(UNUSED void* h)
{
  return 0;
}

int rf_blade_open_multi(char* args, void** h, uint32_t nof_channels)
{
  return rf_blade_open_nof_rx(args, h, nof_channels);
}

int rf_blade_open(char* args, void** h)
{
  return rf_blade_open_nof_rx(args, h, 1);
}

int rf_blade_open_nof_rx(char* args, void** h, uint32_t nof_rx_channels)
{
  const struct bladerf_range* range_tx = NULL;
  const struct bladerf_range* range_rx = NULL;
  *h                                   = NULL;

  if (nof_rx_channels < 1) {
    nof_rx_channels = 1;
  }
  if (nof_rx_channels > BLADE_MAX_RX_CHANNELS) {
    ERROR("bladeRF supports at most %d RX channels, %d requested", BLADE_MAX_RX_CHANNELS, nof_rx_channels);
    return -1;
  }

  rf_blade_handler_t* handler = (rf_blade_handler_t*)malloc(sizeof(rf_blade_handler_t));
  if (!handler) {
    perror("malloc");
    return -1;
  }
  *h                        = handler;
  handler->nof_rx_channels  = nof_rx_channels;

  printf("Opening bladeRF with %d RX channel(s)...\n", nof_rx_channels);
  int status = bladerf_open(&handler->dev, args);
  if (status) {
    ERROR("Unable to open device: %s", bladerf_strerror(status));
    goto clean_exit;
  }

  /* A second chain only exists on a bladeRF 2.0 micro. Asking an x40/x115 for
  one here gives a clearer failure than letting bladerf_sync_config() reject the
  X2 layout later, by which point the cause is several frames away. */
  if (nof_rx_channels > 1) {
    size_t nof_dev_rx = bladerf_get_channel_count(handler->dev, BLADERF_RX);
    if (nof_dev_rx < nof_rx_channels) {
      ERROR("Device has %zu RX channel(s), %d requested. A second RX chain needs a bladeRF 2.0 micro (xA4/xA9)",
            nof_dev_rx,
            nof_rx_channels);
      status = -1;
      goto clean_exit;
    }
  }

  for (uint32_t i = 0; i < nof_rx_channels; i++) {
    status = bladerf_set_gain_mode(handler->dev, BLADERF_CHANNEL_RX(i), BLADERF_GAIN_MGC);
    if (status) {
      ERROR("Failed to set manual gain mode on RX channel %d: %s", i, bladerf_strerror(status));
      goto clean_exit;
    }
  }

  // bladerf_log_set_verbosity(BLADERF_LOG_LEVEL_VERBOSE);

  /* Get Gain ranges and set Rx to maximum */
  status = bladerf_get_gain_range(handler->dev, BLADERF_CHANNEL_RX(0), &range_rx);
  if ((status != 0) || (range_rx == NULL)) {
    ERROR("Failed to get RX gain range: %s", bladerf_strerror(status));
    goto clean_exit;
  }

  status = bladerf_get_gain_range(handler->dev, BLADERF_TX_X1, &range_tx);
  if ((status != 0) || (range_tx == NULL)) {
    ERROR("Failed to get TX gain range: %s", bladerf_strerror(status));
    goto clean_exit;
  }

  /* Same gain on every chain. An AoA estimate reads the phase differences
  between chains, so a per-chain gain difference would show up as an amplitude
  taper across the array and bias the spatial spectrum. */
  for (uint32_t i = 0; i < nof_rx_channels; i++) {
    status = bladerf_set_gain(handler->dev, BLADERF_CHANNEL_RX(i), (bladerf_gain)range_rx->max);
    if (status != 0) {
      ERROR("Failed to set RX LNA gain on channel %d: %s", i, bladerf_strerror(status));
      goto clean_exit;
    }
  }
  handler->rx_stream_enabled = false;
  handler->tx_stream_enabled = false;

  /* Set default sampling rates */
  rf_blade_set_tx_srate(handler, 1.92e6);
  rf_blade_set_rx_srate(handler, 1.92e6);

  /* Set info structure */
  handler->info.min_tx_gain = range_tx->min;
  handler->info.max_tx_gain = range_tx->max;
  handler->info.min_rx_gain = range_rx->min;
  handler->info.max_rx_gain = range_rx->max;

  return SRSRAN_SUCCESS;

clean_exit:
  free(handler);
  return status;
}

int rf_blade_close(void* h)
{
  rf_blade_handler_t* handler = (rf_blade_handler_t*)h;
  bladerf_close(handler->dev);
  return 0;
}

double rf_blade_set_rx_srate(void* h, double freq)
{
  uint32_t            bw      = 0;
  rf_blade_handler_t* handler = (rf_blade_handler_t*)h;

  /* Every chain is configured identically. On the AD9361 the two RX chains
  share a sample clock, so this is really asserting that shared setting once per
  chain rather than setting two independent rates. */
  for (uint32_t i = 0; i < handler->nof_rx_channels; i++) {
    int status = bladerf_set_sample_rate(handler->dev, BLADERF_CHANNEL_RX(i), (uint32_t)freq, &handler->rx_rate);
    if (status != 0) {
      ERROR("Failed to set samplerate = %u on channel %d: %s", (uint32_t)freq, i, bladerf_strerror(status));
      return -1;
    }

    const bladerf_bandwidth want_bw =
        (handler->rx_rate < 2000000) ? handler->rx_rate : (bladerf_bandwidth)(handler->rx_rate * 0.8);
    status = bladerf_set_bandwidth(handler->dev, BLADERF_CHANNEL_RX(i), want_bw, &bw);
    if (status != 0) {
      ERROR("Failed to set bandwidth = %u on channel %d: %s", handler->rx_rate, i, bladerf_strerror(status));
      return -1;
    }
  }

  printf("Set RX sampling rate %.2f Mhz, filter BW: %.2f Mhz, %d channel(s)\n",
         (float)handler->rx_rate / 1e6,
         (float)bw / 1e6,
         handler->nof_rx_channels);
  return (double)handler->rx_rate;
}

double rf_blade_set_tx_srate(void* h, double freq)
{
  uint32_t            bw;
  rf_blade_handler_t* handler = (rf_blade_handler_t*)h;
  int                 status  = bladerf_set_sample_rate(handler->dev, BLADERF_TX_X1, (uint32_t)freq, &handler->tx_rate);
  if (status != 0) {
    ERROR("Failed to set samplerate = %u: %s", (uint32_t)freq, bladerf_strerror(status));
    return -1;
  }
  status = bladerf_set_bandwidth(handler->dev, BLADERF_TX_X1, handler->tx_rate, &bw);
  if (status != 0) {
    ERROR("Failed to set bandwidth = %u: %s", handler->tx_rate, bladerf_strerror(status));
    return -1;
  }
  return (double)handler->tx_rate;
}

int rf_blade_set_rx_gain(void* h, double gain)
{
  int                 status  = 0;
  rf_blade_handler_t* handler = (rf_blade_handler_t*)h;
  /* Applied to every chain, so the array stays amplitude balanced. AGC drives
  this, and a gain step landing on only one chain would look like the scene
  changing on that chain alone. */
  for (uint32_t i = 0; i < handler->nof_rx_channels && status == 0; i++) {
    status = bladerf_set_gain(handler->dev, BLADERF_CHANNEL_RX(i), (bladerf_gain)gain);
  }
  if (status != 0) {
    ERROR("Failed to set RX gain: %s", bladerf_strerror(status));
    return SRSRAN_ERROR;
  }
  return SRSRAN_SUCCESS;
}

int rf_blade_set_rx_gain_ch(void* h, uint32_t ch, double gain)
{
  return rf_blade_set_rx_gain(h, gain);
}

int rf_blade_set_tx_gain(void* h, double gain)
{
  int                 status;
  rf_blade_handler_t* handler = (rf_blade_handler_t*)h;
  status                      = bladerf_set_gain(handler->dev, BLADERF_TX_X1, (bladerf_gain)gain);
  if (status != 0) {
    ERROR("Failed to set TX gain: %s", bladerf_strerror(status));
    return SRSRAN_ERROR;
  }
  return SRSRAN_SUCCESS;
}

int rf_blade_set_tx_gain_ch(void* h, uint32_t ch, double gain)
{
  return rf_blade_set_tx_gain(h, gain);
}

double rf_blade_get_rx_gain(void* h)
{
  int                 status;
  bladerf_gain        gain    = 0;
  rf_blade_handler_t* handler = (rf_blade_handler_t*)h;
  status                      = bladerf_get_gain(handler->dev, BLADERF_RX_X1, &gain);
  if (status != 0) {
    ERROR("Failed to get RX gain: %s", bladerf_strerror(status));
    return -1;
  }
  return gain;
}

double rf_blade_get_tx_gain(void* h)
{
  int                 status;
  bladerf_gain        gain    = 0;
  rf_blade_handler_t* handler = (rf_blade_handler_t*)h;
  status                      = bladerf_get_gain(handler->dev, BLADERF_TX_X1, &gain);
  if (status != 0) {
    ERROR("Failed to get TX gain: %s", bladerf_strerror(status));
    return -1;
  }
  return gain;
}

srsran_rf_info_t* rf_blade_get_info(void* h)
{
  srsran_rf_info_t* info = NULL;

  if (h) {
    rf_blade_handler_t* handler = (rf_blade_handler_t*)h;

    info = &handler->info;
  }
  return info;
}

double rf_blade_set_rx_freq(void* h, uint32_t ch, double freq)
{
  rf_blade_handler_t* handler = (rf_blade_handler_t*)h;
  bladerf_frequency   f_int   = (uint32_t)round(freq);

  /* ch is honoured rather than ignored, because srsran::radio::set_rx_freq()
  already walks the antennas of a carrier and calls this once per chain. Tuning
  channel 0 every time, as this did before, left a second chain on whatever
  frequency it powered up with. */
  if (ch >= handler->nof_rx_channels) {
    ERROR("set_rx_freq: channel %d out of range, %d configured", ch, handler->nof_rx_channels);
    return -1;
  }

  int status = bladerf_set_frequency(handler->dev, BLADERF_CHANNEL_RX(ch), f_int);
  if (status != 0) {
    ERROR("Failed to set RX frequency = %u on channel %d: %s", (uint32_t)freq, ch, bladerf_strerror(status));
    return -1;
  }
  f_int = 0;
  bladerf_get_frequency(handler->dev, BLADERF_CHANNEL_RX(ch), &f_int);
  printf("set RX channel %d frequency to %lu\n", ch, f_int);

  return freq;
}

double rf_blade_set_tx_freq(void* h, UNUSED uint32_t ch, double freq)
{
  rf_blade_handler_t* handler = (rf_blade_handler_t*)h;
  bladerf_frequency   f_int   = (uint32_t)round(freq);
  int                 status  = bladerf_set_frequency(handler->dev, BLADERF_TX_X1, f_int);
  if (status != 0) {
    ERROR("Failed to set samplerate = %u: %s", (uint32_t)freq, bladerf_strerror(status));
    return -1;
  }

  f_int = 0;
  bladerf_get_frequency(handler->dev, BLADERF_TX_X1, &f_int);
  printf("set TX frequency to %lu\n", f_int);
  return freq;
}

static void timestamp_to_secs(uint32_t rate, uint64_t timestamp, time_t* secs, double* frac_secs)
{
  double totalsecs = (double)timestamp / rate;
  time_t secs_i    = (time_t)totalsecs;
  if (secs) {
    *secs = secs_i;
  }
  if (frac_secs) {
    *frac_secs = totalsecs - secs_i;
  }
}

void rf_blade_get_time(void* h, time_t* secs, double* frac_secs)
{
  rf_blade_handler_t*     handler = (rf_blade_handler_t*)h;
  struct bladerf_metadata meta;

  int status = bladerf_get_timestamp(handler->dev, BLADERF_RX, &meta.timestamp);
  if (status != 0) {
    ERROR("Failed to get current RX timestamp: %s", bladerf_strerror(status));
  }
  timestamp_to_secs(handler->rx_rate, meta.timestamp, secs, frac_secs);
}

int rf_blade_recv_with_time_multi(void*    h,
                                  void**   data,
                                  uint32_t nsamples,
                                  bool     blocking,
                                  time_t*  secs,
                                  double*  frac_secs)
{
  rf_blade_handler_t* handler = (rf_blade_handler_t*)h;

  /* One chain is the old path unchanged, including its conversion call. */
  if (handler->nof_rx_channels < 2) {
    return rf_blade_recv_with_time(h, *data, nsamples, blocking, secs, frac_secs);
  }

  struct bladerf_metadata meta;
  int                     status;
  const uint32_t          nof_ch = handler->nof_rx_channels;

  memset(&meta, 0, sizeof(meta));
  meta.flags = BLADERF_META_FLAG_RX_NOW;

  /* bladerf_sync_rx() counts samples across all chains in a multi-channel
  layout, and each is an I/Q pair of int16. */
  if (2 * nsamples * nof_ch > CONVERT_BUFFER_SIZE) {
    ERROR("RX failed: nsamples exceeds buffer size (%d*%d>%d)", nsamples, nof_ch, CONVERT_BUFFER_SIZE);
    return -1;
  }

  status = bladerf_sync_rx(handler->dev, handler->rx_buffer, nsamples * nof_ch, &meta, 2000);
  if (status) {
    ERROR("RX failed: %s; nsamples=%d; channels=%d", bladerf_strerror(status), nsamples, nof_ch);
    return -1;
  } else if (meta.status & BLADERF_META_STATUS_OVERRUN) {
    if (blade_error_handler) {
      srsran_rf_error_t error;
      error.opt  = meta.actual_count;
      error.type = SRSRAN_RF_ERROR_OVERFLOW;
      blade_error_handler(blade_error_handler_arg, error);
    }
  }

  timestamp_to_secs(handler->rx_rate, meta.timestamp, secs, frac_secs);

  /* The stream arrives interleaved, one I/Q pair per chain per sample instant:
  ch0_I, ch0_Q, ch1_I, ch1_Q, ch0_I, ... Deinterleaving it by hand with a stride
  is a scalar, cache-hostile pass over 61 M samples/s at 2x30.72 Msps, and it
  sits directly on the receive thread, which is what makes the device FIFO
  overrun.

  So the reordering is left to libbladeRF, which does it in place and knows the
  format, and the conversion is then two contiguous SIMD passes.
  The 1/2048 matches the scale the single channel path passes to
  srsran_vec_convert_if(), so both reach the same full scale. */
  /* Plain SC16_Q11 and not the _META variant the stream is configured with:
  bladerf_sync_rx() reports the metadata through its own argument and hands back
  a buffer of samples alone, whereas passing _META here would skip the first 16
  bytes as a header that is not there and shift every chain by two samples. */
  status = bladerf_deinterleave_stream_buffer(
      blade_rx_layout(nof_ch), BLADERF_FORMAT_SC16_Q11, nsamples * nof_ch, handler->rx_buffer);
  if (status != 0) {
    ERROR("RX deinterleave failed: %s", bladerf_strerror(status));
    return -1;
  }

  /* Now laid out as nof_ch contiguous blocks of nsamples, so chain c starts at
  int16 offset 2 * c * nsamples. */
  for (uint32_t c = 0; c < nof_ch; c++) {
    if (data[c] == NULL) {
      continue;
    }
    srsran_vec_convert_if(&handler->rx_buffer[2 * c * nsamples], 2048, (float*)data[c], 2 * nsamples);
  }

  return nsamples;
}

int rf_blade_recv_with_time(void*       h,
                            void*       data,
                            uint32_t    nsamples,
                            UNUSED bool blocking,
                            time_t*     secs,
                            double*     frac_secs)
{
  rf_blade_handler_t*     handler = (rf_blade_handler_t*)h;
  struct bladerf_metadata meta;
  int                     status;

  memset(&meta, 0, sizeof(meta));
  meta.flags = BLADERF_META_FLAG_RX_NOW;

  if (2 * nsamples > CONVERT_BUFFER_SIZE) {
    ERROR("RX failed: nsamples exceeds buffer size (%d>%d)", nsamples, CONVERT_BUFFER_SIZE);
    return -1;
  }
  status = bladerf_sync_rx(handler->dev, handler->rx_buffer, nsamples, &meta, 2000);
  if (status) {
    ERROR("RX failed: %s; nsamples=%d;", bladerf_strerror(status), nsamples);
    return -1;
  } else if (meta.status & BLADERF_META_STATUS_OVERRUN) {
    if (blade_error_handler) {
      srsran_rf_error_t error;
      error.opt  = meta.actual_count;
      error.type = SRSRAN_RF_ERROR_OVERFLOW;
      blade_error_handler(blade_error_handler_arg, error);
    } else {
      /*ERROR("Overrun detected in scheduled RX. "
            "%u valid samples were read.", meta.actual_count);*/
    }
  }

  timestamp_to_secs(handler->rx_rate, meta.timestamp, secs, frac_secs);
  srsran_vec_convert_if(handler->rx_buffer, 2048, data, 2 * nsamples);

  return nsamples;
}

int rf_blade_send_timed_multi(void*  h,
                              void*  data[4],
                              int    nsamples,
                              time_t secs,
                              double frac_secs,
                              bool   has_time_spec,
                              bool   blocking,
                              bool   is_start_of_burst,
                              bool   is_end_of_burst)
{
  return rf_blade_send_timed(
      h, data[0], nsamples, secs, frac_secs, has_time_spec, blocking, is_start_of_burst, is_end_of_burst);
}

int rf_blade_send_timed(void*       h,
                        void*       data,
                        int         nsamples,
                        time_t      secs,
                        double      frac_secs,
                        bool        has_time_spec,
                        UNUSED bool blocking,
                        bool        is_start_of_burst,
                        bool        is_end_of_burst)
{
  rf_blade_handler_t*     handler = (rf_blade_handler_t*)h;
  struct bladerf_metadata meta;
  int                     status;

  if (!handler->tx_stream_enabled) {
    rf_blade_start_tx_stream(h);
  }

  if (2 * nsamples > CONVERT_BUFFER_SIZE) {
    ERROR("TX failed: nsamples exceeds buffer size (%d>%d)", nsamples, CONVERT_BUFFER_SIZE);
    return -1;
  }

  srsran_vec_convert_fi(data, 2048, handler->tx_buffer, 2 * nsamples);

  memset(&meta, 0, sizeof(meta));
  if (is_start_of_burst) {
    if (has_time_spec) {
      // Convert time to ticks
      srsran_timestamp_t ts = {.full_secs = secs, .frac_secs = frac_secs};
      meta.timestamp        = srsran_timestamp_uint64(&ts, handler->tx_rate);
    } else {
      meta.flags |= BLADERF_META_FLAG_TX_NOW;
    }
    meta.flags |= BLADERF_META_FLAG_TX_BURST_START;
  }
  if (is_end_of_burst) {
    meta.flags |= BLADERF_META_FLAG_TX_BURST_END;
  }
  srsran_rf_error_t error;
  bzero(&error, sizeof(srsran_rf_error_t));

  status = bladerf_sync_tx(handler->dev, handler->tx_buffer, nsamples, &meta, 2000);
  if (status == BLADERF_ERR_TIME_PAST) {
    if (blade_error_handler) {
      error.type = SRSRAN_RF_ERROR_LATE;
      blade_error_handler(blade_error_handler_arg, error);
    } else {
      ERROR("TX failed: %s", bladerf_strerror(status));
    }
  } else if (status) {
    ERROR("TX failed: %s", bladerf_strerror(status));
    return status;
  } else if (meta.status == BLADERF_META_STATUS_UNDERRUN) {
    if (blade_error_handler) {
      error.type = SRSRAN_RF_ERROR_UNDERFLOW;
      blade_error_handler(blade_error_handler_arg, error);
    } else {
      ERROR("TX warning: underflow detected.");
    }
  }

  return nsamples;
}

rf_dev_t srsran_rf_dev_blade = {"bladeRF",
                                rf_blade_devname,
                                rf_blade_start_rx_stream,
                                rf_blade_stop_rx_stream,
                                rf_blade_flush_buffer,
                                rf_blade_has_rssi,
                                rf_blade_get_rssi,
                                rf_blade_suppress_stdout,
                                rf_blade_register_error_handler,
                                rf_blade_open,
                                .srsran_rf_open_multi = rf_blade_open_multi,
                                rf_blade_close,
                                rf_blade_set_rx_srate,
                                rf_blade_set_rx_gain,
                                rf_blade_set_rx_gain_ch,
                                rf_blade_set_tx_gain,
                                rf_blade_set_tx_gain_ch,
                                rf_blade_get_rx_gain,
                                rf_blade_get_tx_gain,
                                rf_blade_get_info,
                                rf_blade_set_rx_freq,
                                rf_blade_set_tx_srate,
                                rf_blade_set_tx_freq,
                                rf_blade_get_time,
                                NULL,
                                rf_blade_recv_with_time,
                                rf_blade_recv_with_time_multi,
                                rf_blade_send_timed,
                                .srsran_rf_send_timed_multi = rf_blade_send_timed_multi};

#ifdef ENABLE_RF_PLUGINS
int register_plugin(rf_dev_t** rf_api)
{
  if (rf_api == NULL) {
    return SRSRAN_ERROR;
  }
  *rf_api = &srsran_rf_dev_blade;
  return SRSRAN_SUCCESS;
}
#endif /* ENABLE_RF_PLUGINS */
