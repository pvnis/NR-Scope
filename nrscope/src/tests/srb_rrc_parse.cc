/*
 * Parses the SRB PDUs the RACH decoder logged (msg4_bytes.log, outcome=mac:srb) as RRC
 * DL-DCCH messages. They are integrity protected and, unless the network ciphers RRC
 * (NEA1-3), readable: OCUDU's default nea_pref_list starts with nea0, null ciphering.
 * Prints each message's type and, for an RRCReconfiguration, its masterCellGroup as JSON.
 *
 * usage: srb_rrc_parse [msg4_bytes.log] [rnti]
 */
#include "srsran/asn1/rrc_nr.h"
#include "srsran/mac/mac_sch_pdu_nr.h"
#include <cstdio>
#include <cctype>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

static std::vector<uint8_t> unhex(const std::string& h)
{
  std::vector<uint8_t> v;
  for (size_t i = 0; i + 1 < h.size() && isxdigit((unsigned char)h[i]) && isxdigit((unsigned char)h[i + 1]); i += 2) {
    v.push_back((uint8_t)std::stoul(h.substr(i, 2), nullptr, 16));
  }
  return v;
}

int main(int argc, char** argv)
{
  const char* path   = argc > 1 ? argv[1] : "msg4_bytes.log";
  long        want   = argc > 2 ? strtol(argv[2], nullptr, 16) : -1;
  std::ifstream in(path);
  std::string   line;
  int           n = 0, ok = 0;
  while (std::getline(in, line)) {
    if (line.find("outcome=mac:srb") == std::string::npos) {
      continue;
    }
    size_t r = line.find("rnti=0x"), b = line.rfind(" bytes=");
    if (r == std::string::npos || b == std::string::npos) {
      continue;
    }
    uint16_t rnti = (uint16_t)strtol(line.c_str() + r + 7, nullptr, 16);
    if (want >= 0 && rnti != want) {
      continue;
    }
    std::vector<uint8_t>  tb = unhex(line.substr(b + 7));
    srsran::mac_sch_pdu_nr mac(false);
    if (mac.unpack(tb.data(), tb.size()) != SRSRAN_SUCCESS) {
      continue;
    }
    for (uint32_t i = 0; i < mac.get_num_subpdus(); i++) {
      srsran::mac_sch_subpdu_nr& sp = mac.get_subpdu(i);
      if (!sp.is_sdu() || sp.get_lcid() < 1 || sp.get_lcid() > 3) {
        continue;
      }
      const uint8_t* sdu = sp.get_sdu();
      uint32_t       len = sp.get_sdu_length();
      uint32_t       si  = (sdu[0] >> 4) & 0x3;
      if (si != 0) {
        printf("0x%04x SRB%u: RLC segment (SI %u), %u bytes: not parsed\n", rnti, sp.get_lcid(), si, len);
        continue;
      }
      /* RLC AM 12-bit SN (2 bytes), PDCP SRB header (2 bytes), RRC, MAC-I (4 bytes) */
      if (len < 2 + 2 + 4 + 1) {
        continue;
      }
      const uint8_t* rrc     = sdu + 4;
      uint32_t       rrc_len = len - 4 - 4;
      uint32_t       pdcp_sn = ((sdu[2] & 0x0f) << 8) | sdu[3];
      n++;
      asn1::rrc_nr::dl_dcch_msg_s msg;
      asn1::cbit_ref              bref(rrc, rrc_len);
      if (msg.unpack(bref) != asn1::SRSASN_SUCCESS ||
          msg.msg.type().value != asn1::rrc_nr::dl_dcch_msg_type_c::types_opts::c1) {
        printf("0x%04x SRB%u PDCP SN %u, %u bytes: does not parse as DL-DCCH (ciphered?)\n", rnti, sp.get_lcid(),
               pdcp_sn, len);
        continue;
      }
      ok++;
      printf("0x%04x SRB%u PDCP SN %u, %u bytes: %s\n", rnti, sp.get_lcid(), pdcp_sn, len,
             msg.msg.c1().type().to_string());
      if (msg.msg.c1().type().value == asn1::rrc_nr::dl_dcch_msg_type_c::c1_c_::types_opts::rrc_recfg) {
        const auto& recfg = msg.msg.c1().rrc_recfg();
        if (recfg.crit_exts.type().value == asn1::rrc_nr::rrc_recfg_s::crit_exts_c_::types_opts::rrc_recfg &&
            recfg.crit_exts.rrc_recfg().non_crit_ext_present &&
            recfg.crit_exts.rrc_recfg().non_crit_ext.master_cell_group.size() > 0) {
          const auto&                    mcg = recfg.crit_exts.rrc_recfg().non_crit_ext.master_cell_group;
          asn1::rrc_nr::cell_group_cfg_s cg;
          asn1::cbit_ref                 cbr(mcg.data(), mcg.size());
          if (cg.unpack(cbr) == asn1::SRSASN_SUCCESS) {
            asn1::json_writer js;
            cg.to_json(js);
            printf("  masterCellGroup: %s\n", js.to_string().c_str());
          }
        }
      }
    }
  }
  printf("%d SRB PDUs, %d parsed as RRC\n", n, ok);
  return 0;
}
