// rec_direct: continuous multi-channel IQ recording without dropped samples.
//
// One UHD rx streamer (all channels sample-aligned, timed start) receives into a large RAM ring
// allocated and touched up front. One writer thread per channel writes the ring to
// <dir>/iq_ch<c>.sc16 with O_DIRECT, so the receive loop never waits for the disk or the page cache.
// Gaps (radio overflow, ring full) are found from the radio timestamps, logged to gaps.csv and filled
// with zeros, so sample n of every file is always at time t0 + n/rate.
//
// Build: g++ -O2 -std=c++17 -o rec_direct rec_direct.cpp -luhd -lpthread
// Run (DPDK needs root): sudo ./rec_direct --dir D --duration 120 [--args ..] [--rate 122.88e6]
//   [--freq 3750e6] [--gain 38] [--ant TX/RX0] [--channels 0,1,2,3] [--ram-gb 40] [--setup 1]
#include <uhd/usrp/multi_usrp.hpp>
#include <uhd/utils/thread.hpp>
#include <atomic>
#include <chrono>
#include <cmath>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <string>
#include <sys/stat.h>
#include <thread>
#include <unistd.h>
#include <vector>

static std::atomic<bool> g_stop{false};
static void on_sig(int) { g_stop = true; }

static const size_t BLK_SAMPS = size_t(1) << 22;  // 4 Mi samples = 16 MiB per channel per block
static const size_t SAMP_BYTES = 4;               // sc16 I/Q

int main(int argc, char** argv)
{
    std::string args = "type=x4xx,master_clock_rate=245760000,addr=192.168.20.2,mgmt_addr=192.168.1.190,use_dpdk=1234";
    std::string dir, ant = "TX/RX0", chs = "0,1,2,3";
    double rate = 122.88e6, freq = 3750e6, gain = 38, duration = 10, ram_gb = 40, setup = 1;
    for (int i = 1; i + 1 < argc; i += 2) {
        std::string k = argv[i], v = argv[i + 1];
        if (k == "--args") args = v;
        else if (k == "--dir") dir = v;
        else if (k == "--ant") ant = v;
        else if (k == "--channels") chs = v;
        else if (k == "--rate") rate = atof(v.c_str());
        else if (k == "--freq") freq = atof(v.c_str());
        else if (k == "--gain") gain = atof(v.c_str());
        else if (k == "--duration") duration = atof(v.c_str());
        else if (k == "--ram-gb") ram_gb = atof(v.c_str());
        else if (k == "--setup") setup = atof(v.c_str());
        else { fprintf(stderr, "unknown option %s\n", k.c_str()); return 2; }
    }
    if (dir.empty()) { fprintf(stderr, "--dir is required\n"); return 2; }
    std::vector<size_t> chans;
    for (size_t p = 0; p < chs.size();) {
        size_t q = chs.find(',', p); if (q == std::string::npos) q = chs.size();
        chans.push_back(std::stoul(chs.substr(p, q - p))); p = q + 1;
    }
    const size_t nch = chans.size();
    const uint64_t total = uint64_t(llround(duration * rate));  // samples per channel to record
    mkdir(dir.c_str(), 0755);

    // RAM ring: nblk blocks x nch channels x 16 MiB, aligned for O_DIRECT, touched now.
    const size_t blk_bytes = BLK_SAMPS * SAMP_BYTES;
    const size_t nblk = std::max<size_t>(4, size_t(ram_gb * 1e9 / double(blk_bytes * nch)));
    std::vector<std::vector<char*>> ring(nblk, std::vector<char*>(nch));
    for (auto& b : ring) for (auto& p : b) {
        if (posix_memalign((void**)&p, 4096, blk_bytes)) { fprintf(stderr, "ring allocation failed\n"); return 1; }
        memset(p, 0, blk_bytes);
    }
    printf("ring: %zu blocks x %zu ch x %zu MiB = %.1f GB (%.1f s of recording)\n", nblk, nch, blk_bytes >> 20,
        nblk * nch * blk_bytes / 1e9, nblk * double(BLK_SAMPS) / rate);

    std::vector<int> fds(nch);
    for (size_t c = 0; c < nch; c++) {
        std::string f = dir + "/iq_ch" + std::to_string(c) + ".sc16";
        fds[c] = open(f.c_str(), O_WRONLY | O_CREAT | O_TRUNC | O_DIRECT, 0644);
        if (fds[c] < 0) { perror(f.c_str()); return 1; }
    }
    FILE* gaps = fopen((dir + "/gaps.csv").c_str(), "w");
    fprintf(gaps, "sample_index,gap_samples,reason\n");

    // produced = blocks handed to the writers; done[c] = blocks channel c has written.
    // Writers always write whole blocks; the files are trimmed to the exact length at the end.
    std::atomic<uint64_t> produced{0};
    std::vector<std::atomic<uint64_t>> done(nch);
    for (auto& d : done) d = 0;
    std::atomic<bool> finished{false};
    std::atomic<int> write_err{0};

    std::vector<std::thread> writers;
    for (size_t c = 0; c < nch; c++) writers.emplace_back([&, c] {
        uint64_t k = 0;
        while (true) {
            uint64_t p = produced.load(std::memory_order_acquire);
            if (k == p) { if (finished.load()) break; std::this_thread::sleep_for(std::chrono::milliseconds(1)); continue; }
            const size_t nw = blk_bytes;
            const char* src = ring[k % nblk][c];
            off_t off = off_t(k) * off_t(blk_bytes);
            for (size_t w = 0; w < nw;) {
                ssize_t r = pwrite(fds[c], src + w, nw - w, off + off_t(w));
                if (r <= 0) { perror("pwrite"); write_err = 1; g_stop = true; break; }
                w += size_t(r);
            }
            done[c].store(++k, std::memory_order_release);
        }
    });

    signal(SIGINT, on_sig); signal(SIGTERM, on_sig);
    auto usrp = uhd::usrp::multi_usrp::make(args);
    usrp->set_rx_rate(rate);
    for (size_t i = 0; i < nch; i++) {
        usrp->set_rx_freq(uhd::tune_request_t(freq), chans[i]);
        usrp->set_rx_gain(gain, chans[i]);
        usrp->set_rx_antenna(ant, chans[i]);
    }
    printf("rate %.3f Msps, freq %.3f MHz, gain %.1f dB, antenna %s, channels %s\n",
        usrp->get_rx_rate() / 1e6, usrp->get_rx_freq(chans[0]) / 1e6, usrp->get_rx_gain(chans[0]), ant.c_str(), chs.c_str());
    uhd::stream_args_t sa("sc16", "sc16");
    sa.channels = chans;
    auto rx = usrp->get_rx_stream(sa);
    usrp->set_time_now(uhd::time_spec_t(0.0));
    std::this_thread::sleep_for(std::chrono::milliseconds(200));
    uhd::set_thread_priority_safe(1.0, true);

    uhd::stream_cmd_t cmd(uhd::stream_cmd_t::STREAM_MODE_START_CONTINUOUS);
    cmd.stream_now = false;
    cmd.time_spec = usrp->get_time_now() + uhd::time_spec_t(setup);
    rx->issue_stream_cmd(cmd);

    const size_t spp = rx->get_max_num_samps();
    std::vector<std::vector<char>> scratch(nch, std::vector<char>(spp * 64 * SAMP_BYTES));
    std::vector<void*> bp(nch);
    uhd::rx_metadata_t md;
    uint64_t idx = 0;        // samples per channel placed in the ring (data + zero fill)
    uint64_t blk = 0;        // block being filled
    size_t fill = 0;         // samples in that block
    bool have_t0 = false;
    uhd::time_spec_t t0;
    uint64_t n_ovf = 0, n_seq = 0, n_gap = 0, gap_samps = 0, ring_full_samps = 0;
    double timeout = setup + 2.0;
    auto t_last = std::chrono::steady_clock::now();
    int rc = 0;

    auto ring_has_room = [&](uint64_t b) {
        uint64_t mn = UINT64_MAX;
        for (auto& d : done) mn = std::min<uint64_t>(mn, d.load(std::memory_order_acquire));
        return b - mn < nblk;
    };
    auto close_block = [&]() { produced.store(++blk, std::memory_order_release); fill = 0; };
    auto wait_room = [&]() { while (!ring_has_room(blk) && !g_stop) std::this_thread::sleep_for(std::chrono::milliseconds(1)); };
    // zero-fill n samples into the ring (gap filler); waits for the writers if the ring is full.
    auto zero_fill = [&](uint64_t n) {
        while (n > 0 && idx < total && !g_stop) {
            wait_room();
            size_t k = std::min<uint64_t>({n, BLK_SAMPS - fill, total - idx});
            for (size_t c = 0; c < nch; c++) memset(ring[blk % nblk][c] + fill * SAMP_BYTES, 0, k * SAMP_BYTES);
            fill += k; idx += k; n -= k;
            if (fill == BLK_SAMPS) close_block();
        }
    };
    bool dropped_since = false;  // ring-full samples dropped since the last stored buffer

    while (!g_stop && idx < total) {
        bool room = ring_has_room(blk);
        size_t want;
        if (room) {
            want = std::min<uint64_t>(BLK_SAMPS - fill, total - idx);
            for (size_t c = 0; c < nch; c++) bp[c] = ring[blk % nblk][c] + fill * SAMP_BYTES;
        } else {  // ring full: keep the radio drained, the timestamp jump fills the hole later
            want = spp * 64;
            for (size_t c = 0; c < nch; c++) bp[c] = scratch[c].data();
        }
        size_t n = rx->recv(bp, want, md, timeout, false);
        timeout = 0.5;

        if (md.error_code == uhd::rx_metadata_t::ERROR_CODE_TIMEOUT) { fprintf(stderr, "receive timeout\n"); rc = 1; break; }
        if (md.error_code == uhd::rx_metadata_t::ERROR_CODE_OVERFLOW) { if (md.out_of_sequence) n_seq++; else n_ovf++; }
        else if (md.error_code != uhd::rx_metadata_t::ERROR_CODE_NONE) { fprintf(stderr, "receive error: %s\n", md.strerror().c_str()); }
        if (n == 0) continue;

        if (!have_t0) { t0 = md.time_spec; have_t0 = true; }
        int64_t at = (md.time_spec - t0).to_ticks(rate);  // sample index of the first returned sample
        if (!room) { ring_full_samps += n; dropped_since = true; continue; }
        if (at > int64_t(idx)) {  // samples missing before this buffer
            uint64_t g = uint64_t(at) - idx;
            fprintf(gaps, "%llu,%llu,%s\n", (unsigned long long)idx, (unsigned long long)g, dropped_since ? "ring_full" : "radio");
            n_gap++; gap_samps += g;
            // move the just-received data behind the hole: copy it out, fill zeros, put it back
            std::vector<std::vector<char>> tmp(nch);
            for (size_t c = 0; c < nch; c++) tmp[c].assign((char*)bp[c], (char*)bp[c] + n * SAMP_BYTES);
            zero_fill(g);
            size_t off = 0;
            while (off < n && idx < total && !g_stop) {
                wait_room();
                size_t k = std::min<uint64_t>({n - off, BLK_SAMPS - fill, total - idx});
                for (size_t c = 0; c < nch; c++) memcpy(ring[blk % nblk][c] + fill * SAMP_BYTES, tmp[c].data() + off * SAMP_BYTES, k * SAMP_BYTES);
                fill += k; idx += k; off += k;
                if (fill == BLK_SAMPS) close_block();
            }
        } else if (at < int64_t(idx)) {
            fprintf(gaps, "%llu,%lld,time_went_back\n", (unsigned long long)idx, (long long)(at - int64_t(idx)));
            fill += n; idx += n;
            if (fill == BLK_SAMPS) close_block();
        } else {
            fill += n; idx += n;
            if (fill == BLK_SAMPS) close_block();
        }

        dropped_since = false;
        auto now = std::chrono::steady_clock::now();
        if (now - t_last > std::chrono::seconds(1)) {
            t_last = now;
            uint64_t mn = UINT64_MAX;
            for (auto& d : done) mn = std::min<uint64_t>(mn, d.load());
            printf("t %7.1f s  ring %5.1f%%  overflows %llu  seq_err %llu  gaps %llu (%llu samples)\n", idx / rate,
                100.0 * double(blk - mn) / double(nblk), (unsigned long long)n_ovf, (unsigned long long)n_seq,
                (unsigned long long)n_gap, (unsigned long long)gap_samps);
            fflush(stdout); fflush(gaps);
        }
    }
    rx->issue_stream_cmd(uhd::stream_cmd_t(uhd::stream_cmd_t::STREAM_MODE_STOP_CONTINUOUS));
    for (size_t c = 0; c < nch; c++) bp[c] = scratch[c].data();
    while (rx->recv(bp, spp * 64, md, 0.2, false) > 0 || md.error_code == uhd::rx_metadata_t::ERROR_CODE_OVERFLOW) {}

    // hand over the last partial block, wait for the writers, trim the O_DIRECT padding
    if (fill > 0) produced.store(++blk, std::memory_order_release);
    finished = true;
    for (auto& t : writers) t.join();
    for (size_t c = 0; c < nch; c++) {
        if (ftruncate(fds[c], off_t(idx) * off_t(SAMP_BYTES))) perror("ftruncate");
        close(fds[c]);
    }
    fclose(gaps);

    const char* su = getenv("SUDO_UID"); const char* sg = getenv("SUDO_GID");
    FILE* m = fopen((dir + "/meta.txt").c_str(), "w");
    fprintf(m, "center_hz %.0f\nrate_sps %.0f\ngain_db %.1f\nantenna %s\nchannels %s\nformat sc16 interleaved I/Q, one file per channel iq_ch<c>.sc16, one streamer (sample-aligned)\n"
               "samples_per_channel %llu\nrequested %llu\nt0_radio_s %.9f\noverflows %llu\nseq_errors %llu\ngaps %llu\ngap_samples_zero_filled %llu\n"
               "ring_full_samples_dropped %llu\nwrite_error %d\n",
        freq, rate, gain, ant.c_str(), chs.c_str(), (unsigned long long)idx, (unsigned long long)total, t0.get_real_secs(),
        (unsigned long long)n_ovf, (unsigned long long)n_seq, (unsigned long long)n_gap, (unsigned long long)gap_samps,
        (unsigned long long)ring_full_samps, write_err.load());
    fclose(m);
    if (su && sg) {
        std::string cmdl = "chown -R " + std::string(su) + ":" + std::string(sg) + " '" + dir + "'";
        if (system(cmdl.c_str())) {}
    }
    printf("done: %llu samples per channel (%.3f s), overflows %llu, seq_err %llu, gaps %llu (%llu samples zero-filled)\n",
        (unsigned long long)idx, idx / rate, (unsigned long long)n_ovf, (unsigned long long)n_seq,
        (unsigned long long)n_gap, (unsigned long long)gap_samps);
    return rc || write_err ? 1 : 0;
}
