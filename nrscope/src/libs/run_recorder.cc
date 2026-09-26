#include "nrscope/hdr/run_recorder.h"

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <mutex>
#include <sys/stat.h>
#include <sys/time.h>
#include <unistd.h>

#include "srsran/phy/common/phy_common.h"

#ifndef NRSCOPE_ROOT_DIR
#define NRSCOPE_ROOT_DIR "."
#endif

namespace RunRecorder {
namespace {

/* DCI rows can arrive thousands of times a second from several decoder
  threads, so they are written through a large buffer and flushed at most this
  often, rather than one write per row. Runs usually end by being killed, which
  loses at most this much of the tail. */
const double DCI_FLUSH_PERIOD_S = 0.25;

bool        is_enabled = false;
std::string run_stamp; // when the run started, shared by both files of the run

std::mutex msg4_mtx;
FILE*      msg4_file   = nullptr;
bool       msg4_failed = false;

std::mutex dci_mtx;
FILE*      dci_file       = nullptr;
bool       dci_failed     = false;
double     dci_last_flush = 0;

const char* MSG4_HEADER = "timestamp,pci,sfn,slot,tc_rnti,c_rnti,rrc_transaction_id,rrc_offset,nof_bytes,dci,"
                          "msg4_bytes,master_cell_group\n";

const char* DCI_HEADER =
    "timestamp,pci,sfn,slot,direction,rnti,rnti_type,dci_format,ss_type,coreset_id,aggregation_level,cce,"
    "freq_alloc,time_alloc,dci_mcs,dci_ndi,dci_rv,harq_id,tpc,ports,dmrs_id,srs_request,"
    "k,mapping,time_start,time_length,prbs,nof_prb,nof_layers,"
    "dmrs_type,dmrs_add_pos,dmrs_len,dmrs_typeA_pos,nof_dmrs_cdm_groups,n_scid,beta_dmrs,"
    "modulation,mcs,tbs,code_rate,rv,ndi,nof_re,nof_bits,mcs_table,xoverhead,dci\n";

double now_s()
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return (double)tv.tv_sec + (double)tv.tv_usec / 1e6;
}

/* nrscope runs under sudo: hand what it creates back to the invoking user, so
  the recordings can be moved or deleted without root. */
void give_to_invoking_user(const std::string& path)
{
  const char* uid = getenv("SUDO_UID");
  const char* gid = getenv("SUDO_GID");
  if (uid != nullptr && gid != nullptr) {
    if (chown(path.c_str(), (uid_t)atoi(uid), (gid_t)atoi(gid)) != 0) {
      // Not fatal: the file is still written, just owned by root
    }
  }
}

FILE* open_run_file(const char* dir, const char* prefix, uint32_t pci, const char* header, size_t buf_sz)
{
  std::string dir_path = std::string(NRSCOPE_ROOT_DIR) + "/" + dir;
  mkdir(dir_path.c_str(), 0755); // already there is fine
  give_to_invoking_user(dir_path);

  std::string path = dir_path + "/" + prefix + "_" + run_stamp + "_pci" + std::to_string(pci) + ".csv";
  FILE*       f    = fopen(path.c_str(), "w");
  if (f == nullptr) {
    fprintf(stderr, "RunRecorder: cannot create %s, not recording it\n", path.c_str());
    return nullptr;
  }
  if (buf_sz > 0) {
    setvbuf(f, nullptr, _IOFBF, buf_sz);
  }
  fputs(header, f);
  fflush(f);
  give_to_invoking_user(path);
  printf("Recording to %s\n", path.c_str());
  return f;
}

/* A CSV field in double quotes, with inner quotes doubled. */
void append_quoted(std::string& out, const char* s)
{
  out += '"';
  for (; *s != '\0'; s++) {
    if (*s == '"') {
      out += '"';
    }
    out += *s;
  }
  out += '"';
}

/* The DCI line as printed, without its trailing space. */
std::string trimmed(const char* s)
{
  std::string t(s);
  while (!t.empty() && (t.back() == ' ' || t.back() == '\n')) {
    t.pop_back();
  }
  return t;
}

/* Allocated PRBs as ranges, e.g. "8-128" or "0-3;10-20", which is what the
  DM-RS lattice of the grant spans. */
std::string prb_ranges(const bool* prb_idx)
{
  std::string out;
  int         start = -1;
  for (int i = 0; i <= SRSRAN_MAX_PRB_NR; i++) {
    bool on = i < SRSRAN_MAX_PRB_NR && prb_idx[i];
    if (on && start < 0) {
      start = i;
    } else if (!on && start >= 0) {
      if (!out.empty()) {
        out += ';';
      }
      out += std::to_string(start) + "-" + std::to_string(i - 1);
      start = -1;
    }
  }
  return out;
}

uint32_t dmrs_add_pos_value(srsran_dmrs_sch_add_pos_t p)
{
  switch (p) {
    case srsran_dmrs_sch_add_pos_0:
      return 0;
    case srsran_dmrs_sch_add_pos_1:
      return 1;
    case srsran_dmrs_sch_add_pos_2:
      return 2;
    case srsran_dmrs_sch_add_pos_3:
      return 3;
  }
  return 0;
}

uint32_t xoverhead_value(srsran_xoverhead_t x)
{
  switch (x) {
    case srsran_xoverhead_0:
      return 0;
    case srsran_xoverhead_6:
      return 6;
    case srsran_xoverhead_12:
      return 12;
    case srsran_xoverhead_18:
      return 18;
  }
  return 0;
}

/* The JSON srsRAN prints spans hundreds of lines; one CSV row wants it on one. */
std::string compact_json(const std::string& js)
{
  std::string out;
  out.reserve(js.size());
  bool skip_indent = false;
  for (char c : js) {
    if (c == '\n') {
      skip_indent = true;
      continue;
    }
    if (skip_indent && c == ' ') {
      continue;
    }
    skip_indent = false;
    out += c;
  }
  return out;
}

} // namespace

void init(bool enable)
{
  is_enabled = enable;
  if (!is_enabled) {
    return;
  }
  time_t    t = time(NULL);
  struct tm tm_now;
  localtime_r(&t, &tm_now);
  char stamp[32];
  strftime(stamp, sizeof(stamp), "%Y-%m-%d_%H-%M-%S", &tm_now);
  run_stamp = stamp;
  printf("Recording mode: RRCSetups to %s/msg4/, UE DCIs to %s/DCIs/, terminal kept quiet\n",
         NRSCOPE_ROOT_DIR,
         NRSCOPE_ROOT_DIR);
}

bool enabled()
{
  return is_enabled;
}

void close()
{
  {
    std::lock_guard<std::mutex> lock(msg4_mtx);
    if (msg4_file != nullptr) {
      fclose(msg4_file);
      msg4_file = nullptr;
    }
  }
  {
    std::lock_guard<std::mutex> lock(dci_mtx);
    if (dci_file != nullptr) {
      fclose(dci_file);
      dci_file = nullptr;
    }
  }
}

void record_rrc_setup(uint32_t           pci,
                      uint32_t           sfn,
                      uint32_t           slot_idx,
                      uint16_t           tc_rnti,
                      uint16_t           c_rnti,
                      uint32_t           rrc_transaction_id,
                      const char*        dci_str,
                      uint32_t           rrc_offset,
                      const uint8_t*     msg4_bytes,
                      uint32_t           nof_bytes,
                      const std::string& master_cell_group_json)
{
  if (!is_enabled) {
    return;
  }

  char head[160];
  snprintf(head,
           sizeof(head),
           "%.6f,%u,%u,%u,%u,%u,%u,%u,%u,",
           now_s(),
           pci,
           sfn,
           slot_idx,
           tc_rnti,
           c_rnti,
           rrc_transaction_id,
           rrc_offset,
           nof_bytes);
  std::string row(head);
  append_quoted(row, trimmed(dci_str).c_str());
  row += ',';
  static const char hex[] = "0123456789abcdef";
  for (uint32_t i = 0; i < nof_bytes; i++) {
    row += hex[msg4_bytes[i] >> 4];
    row += hex[msg4_bytes[i] & 0xf];
  }
  row += ',';
  append_quoted(row, compact_json(master_cell_group_json).c_str());
  row += '\n';

  std::lock_guard<std::mutex> lock(msg4_mtx);
  if (msg4_file == nullptr && !msg4_failed) {
    // RRCSetups are rare: default buffering, flushed on every row
    msg4_file   = open_run_file("msg4", "msg4", pci, MSG4_HEADER, 0);
    msg4_failed = msg4_file == nullptr;
  }
  if (msg4_file != nullptr) {
    fputs(row.c_str(), msg4_file);
    fflush(msg4_file);
  }
}

void record_dci(uint32_t                   pci,
                uint32_t                   sfn,
                uint32_t                   slot_idx,
                bool                       downlink,
                const srsran_dci_ctx_t&    ctx,
                uint32_t                   freq_alloc,
                uint32_t                   time_alloc,
                uint32_t                   mcs,
                uint32_t                   ndi,
                uint32_t                   rv,
                uint32_t                   harq_id,
                uint32_t                   tpc,
                uint32_t                   ports,
                uint32_t                   dmrs_id,
                uint32_t                   srs_request,
                const srsran_sch_cfg_nr_t* sch_cfg,
                const char*                dci_str)
{
  if (!is_enabled) {
    return;
  }

  const double t = now_s();
  char         buf[1024];
  int          n = snprintf(buf,
                   sizeof(buf),
                   "%.6f,%u,%u,%u,%s,%u,%s,%s,%s,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,",
                   t,
                   pci,
                   sfn,
                   slot_idx,
                   downlink ? "DL" : "UL",
                   ctx.rnti,
                   srsran_rnti_type_str_short(ctx.rnti_type),
                   srsran_dci_format_nr_string(ctx.format),
                   srsran_ss_type_str(ctx.ss_type),
                   ctx.coreset_id,
                   1U << ctx.location.L,
                   ctx.location.ncce,
                   freq_alloc,
                   time_alloc,
                   mcs,
                   ndi,
                   rv,
                   harq_id,
                   tpc,
                   ports,
                   dmrs_id,
                   srs_request);
  std::string row(buf, n > 0 ? (size_t)n : 0);

  if (sch_cfg != nullptr) {
    const srsran_sch_grant_nr_t& g  = sch_cfg->grant;
    const srsran_sch_tb_t&       tb = g.tb[0];
    n = snprintf(buf,
                 sizeof(buf),
                 "%u,%s,%u,%u,%s,%u,%u,%u,%u,%s,%u,%u,%u,%.4f,%s,%u,%d,%.4f,%d,%d,%u,%u,%s,%u,",
                 g.k,
                 srsran_sch_mapping_type_to_str(g.mapping),
                 g.S,
                 g.L,
                 prb_ranges(g.prb_idx).c_str(),
                 g.nof_prb,
                 g.nof_layers,
                 sch_cfg->dmrs.type == srsran_dmrs_sch_type_1 ? 1U : 2U,
                 dmrs_add_pos_value(sch_cfg->dmrs.additional_pos),
                 sch_cfg->dmrs.length == srsran_dmrs_sch_len_1 ? "single" : "double",
                 sch_cfg->dmrs.typeA_pos == srsran_dmrs_sch_typeA_pos_2 ? 2U : 3U,
                 g.nof_dmrs_cdm_groups_without_data,
                 g.n_scid ? 1U : 0U,
                 g.beta_dmrs,
                 srsran_mod_string(tb.mod),
                 tb.mcs,
                 tb.tbs,
                 tb.R,
                 tb.rv,
                 tb.ndi,
                 tb.nof_re,
                 tb.nof_bits,
                 srsran_mcs_table_to_str(sch_cfg->sch_cfg.mcs_table),
                 xoverhead_value(sch_cfg->sch_cfg.xoverhead));
    row.append(buf, n > 0 ? (size_t)n : 0);
  } else {
    row += ",,,,,,,,,,,,,,,,,,,,,,,,"; // no grant derived: 24 empty columns
  }
  append_quoted(row, trimmed(dci_str).c_str());
  row += '\n';

  std::lock_guard<std::mutex> lock(dci_mtx);
  if (dci_file == nullptr && !dci_failed) {
    dci_file   = open_run_file("DCIs", "dci", pci, DCI_HEADER, 1 << 20);
    dci_failed = dci_file == nullptr;
  }
  if (dci_file != nullptr) {
    fputs(row.c_str(), dci_file);
    if (t - dci_last_flush >= DCI_FLUSH_PERIOD_S) {
      fflush(dci_file);
      dci_last_flush = t;
    }
  }
}

} // namespace RunRecorder
