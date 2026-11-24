// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "delphipp/io/fz_zebra_parser.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cctype>
#include <cstring>
#include <optional>
#include <set>
#include <unordered_set>

#include "delphipp/io/fz_exchange.hpp"

namespace delphipp::io {

namespace {

std::uint32_t rd_u32_be(const unsigned char* p) noexcept {
  return (static_cast<std::uint32_t>(p[0]) << 24) |
         (static_cast<std::uint32_t>(p[1]) << 16) |
         (static_cast<std::uint32_t>(p[2]) << 8) |
         (static_cast<std::uint32_t>(p[3]));
}

std::uint32_t rd_u32_le(const unsigned char* p) noexcept {
  return (static_cast<std::uint32_t>(p[0])) |
         (static_cast<std::uint32_t>(p[1]) << 8) |
         (static_cast<std::uint32_t>(p[2]) << 16) |
         (static_cast<std::uint32_t>(p[3]) << 24);
}

bool match_sync(const unsigned char* p, bool big_endian) noexcept {
  auto rd = big_endian ? rd_u32_be : rd_u32_le;
  for (std::size_t i = 0; i < 4; ++i) {
    if (rd(p + i * 4) != kFzSyncWords[i]) {
      return false;
    }
  }
  return true;
}

bool is_valid_hollerith_chars(const std::array<char, 4>& c) noexcept {
  if (c[0] < 'A' || c[0] > 'Z') {
    return false;
  }
  bool seen_space = false;
  int alnum = 0;
  for (int i = 0; i < 4; ++i) {
    const unsigned char uc = static_cast<unsigned char>(c[i]);
    const bool is_upper = (uc >= 'A' && uc <= 'Z');
    const bool is_digit = (uc >= '0' && uc <= '9');
    const bool is_sym = (uc == '_' || uc == '-' || uc == '+');
    if (uc == ' ') {
      seen_space = true;
    } else if (is_upper || is_digit || is_sym) {
      if (seen_space) {
        return false;  // Disallow internal spaces such as "A B "
      }
      if (is_upper || is_digit) {
        ++alnum;
      }
    } else {
      return false;
    }
  }
  return alnum >= 2;
}

const std::unordered_set<std::string>& known_bank_priority_set() {
  static const std::unordered_set<std::string> kKnown = {
      "HEAD", "DFIX", "LDTO", "PV  ", "PA  ", "TK  ", "TS  ", "TE  ",
      "TD  ", "EGRA", "EQRA", "EDRA", "EMRA", "EFRA", "MAIN", "EMNC",
      "MUID", "ELID", "HAID", "BTAG", "SIMU", "EVHD", "RUNH", "TX  ",
      "TV  ", "VD  ", "TPC ", "ID  ", "OD  ", "HPC ", "HAC ", "MUF ",
      "MUB ", "MUS ", "EMF ", "TOF ", "SAT ", "STC ", "RIC ", "FCA ",
      "FCB ", "RAW ", "DST ",
  };
  return kKnown;
}

// Decode a 32-bit word into canonical 4-char Hollerith bank ID string if valid.
std::optional<std::string> decode_hollerith_word(std::uint32_t w, bool big_endian) {
  std::array<char, 4> c_fwd{};
  std::array<char, 4> c_rev{};
  if (big_endian) {
    c_fwd = {static_cast<char>((w >> 24) & 0xFFu),
             static_cast<char>((w >> 16) & 0xFFu),
             static_cast<char>((w >> 8) & 0xFFu),
             static_cast<char>(w & 0xFFu)};
    c_rev = {c_fwd[3], c_fwd[2], c_fwd[1], c_fwd[0]};
  } else {
    c_fwd = {static_cast<char>(w & 0xFFu),
             static_cast<char>((w >> 8) & 0xFFu),
             static_cast<char>((w >> 16) & 0xFFu),
             static_cast<char>((w >> 24) & 0xFFu)};
    c_rev = {c_fwd[3], c_fwd[2], c_fwd[1], c_fwd[0]};
  }

  const bool fwd_ok = is_valid_hollerith_chars(c_fwd);
  const bool rev_ok = is_valid_hollerith_chars(c_rev);
  if (!fwd_ok && !rev_ok) {
    return std::nullopt;
  }
  const std::string s_fwd(c_fwd.begin(), c_fwd.end());
  const std::string s_rev(c_rev.begin(), c_rev.end());
  if (fwd_ok && rev_ok) {
    if (known_bank_priority_set().contains(s_rev) &&
        !known_bank_priority_set().contains(s_fwd)) {
      return s_rev;
    }
    return s_fwd;
  }
  return fwd_ok ? s_fwd : s_rev;
}

ParsedZebraBank make_bank_at(const std::vector<std::int32_t>& words,
                             std::size_t p,
                             const std::string& idh) {
  ParsedZebraBank b;
  b.idn = words[p];
  b.idh_raw = words[p + 1];
  b.idh = idh;
  b.nl = words[p + 2];
  b.ns = words[p + 3];
  b.nd = words[p + 4];
  b.status = words[p + 5];

  // Preceding NL links (skipping LQ(L)=next, LQ(L+1)=up, LQ(L+2)=orig -> 3 words)
  if (b.nl > 0 && p >= static_cast<std::size_t>(b.nl + 3)) {
    b.links.reserve(static_cast<std::size_t>(b.nl));
    const std::size_t link_base = p - 3 - static_cast<std::size_t>(b.nl);
    for (std::int32_t j = 0; j < b.nl; ++j) {
      b.links.push_back(words[link_base + static_cast<std::size_t>(j)]);
    }
  }

  b.iwords.reserve(static_cast<std::size_t>(b.nd));
  b.fwords.reserve(static_cast<std::size_t>(b.nd));
  for (std::int32_t k = 0; k < b.nd; ++k) {
    const std::int32_t iw = words[p + 6 + static_cast<std::size_t>(k)];
    b.iwords.push_back(iw);
    b.fwords.push_back(std::bit_cast<float>(static_cast<std::uint32_t>(iw)));
  }
  return b;
}

void update_event_header_from_head(const ParsedZebraBank& head,
                                   ParsedFzEvent& ev) {
  // DELPHI HEAD bank layout:
  //   IHEAD(1): bank version/NWRDS, IHEAD(2): date (YYMMDD), IHEAD(3): time,
  //   IHEAD(4): NRUN (run number), IHEAD(5): IEVT (event number).
  if (head.iwords.size() >= 5 && head.iwords[3] > 0 && head.iwords[4] >= 0) {
    ev.run_number = head.iwords[3];
    ev.event_number = head.iwords[4];
  } else if (head.iwords.size() >= 2 && head.iwords[0] > 0 && head.iwords[1] >= 0) {
    ev.run_number = head.iwords[0];
    ev.event_number = head.iwords[1];
  }
}

void scan_banks_in_words(const std::vector<std::int32_t>& words,
                         bool big_endian,
                         std::size_t max_events,
                         FzParseResult& result) {
  std::optional<ParsedFzEvent> cur_event;
  std::int32_t seq_event = 1;

  auto flush_cur = [&]() {
    if (cur_event.has_value()) {
      result.events.push_back(std::move(*cur_event));
      cur_event.reset();
    }
  };

  for (std::size_t p = 0; p + 6 <= words.size();) {
    const std::uint32_t idh_u = static_cast<std::uint32_t>(words[p + 1]);
    const std::int32_t nl = words[p + 2];
    const std::int32_t ns = words[p + 3];
    const std::int32_t nd = words[p + 4];

    if (nl >= 0 && nl <= 256 && ns >= 0 && ns <= nl && nd >= 0 && nd <= 262144 &&
        p + 6 + static_cast<std::size_t>(nd) <= words.size()) {
      if (const auto idh = decode_hollerith_word(idh_u, big_endian)) {
        // Require either a known bank ID or at least 1 word/link so plain integers don't false-match
        if (known_bank_priority_set().contains(*idh) || nd > 0 || nl > 0) {
          ParsedZebraBank bank = make_bank_at(words, p, *idh);
          ++result.total_banks;
          ++result.bank_counts[*idh];

          if (*idh == "HEAD") {
            flush_cur();
            if (max_events > 0 && result.events.size() >= max_events) {
              return;
            }
            cur_event = ParsedFzEvent{};
            cur_event->event_number = seq_event++;
            update_event_header_from_head(bank, *cur_event);
            cur_event->banks.push_back(std::move(bank));
          } else if (cur_event.has_value()) {
            cur_event->banks.push_back(std::move(bank));
          } else {
            result.preamble_banks.push_back(std::move(bank));
          }

          p += 6 + static_cast<std::size_t>(nd);
          continue;
        }
      }
    }
    ++p;
  }

  flush_cur();
  if (result.events.empty() && !result.preamble_banks.empty()) {
    ParsedFzEvent fallback_ev;
    fallback_ev.event_number = 1;
    fallback_ev.banks = std::move(result.preamble_banks);
    result.preamble_banks.clear();
    result.events.push_back(std::move(fallback_ev));
  }
  if (max_events > 0 && result.events.size() > max_events) {
    result.events.resize(max_events);
  }
}

}  // namespace

const ParsedZebraBank* ParsedFzEvent::find_first(const std::string& idh) const noexcept {
  for (const auto& b : banks) {
    if (b.idh == idh) {
      return &b;
    }
  }
  return nullptr;
}

std::vector<const ParsedZebraBank*> ParsedFzEvent::find_all(const std::string& idh) const {
  std::vector<const ParsedZebraBank*> out;
  for (const auto& b : banks) {
    if (b.idh == idh) {
      out.push_back(&b);
    }
  }
  return out;
}

std::map<std::string, std::int32_t> ParsedFzEvent::populate_zebra_store(
    ZebraStoreView& store, std::int32_t start_lq_offset) const {
  std::map<std::string, std::int32_t> heads;
  std::map<std::string, std::int32_t> tails;
  std::int32_t cursor = std::max<std::int32_t>(32, start_lq_offset);

  for (const auto& b : banks) {
    const std::int32_t nl_alloc = std::max<std::int32_t>(b.nl, 0);
    const std::int32_t nd_alloc = std::max<std::int32_t>(b.nd, 0);
    const std::int32_t l = cursor + nl_alloc + 10;
    const std::int32_t needed_lq_end = l + 8 + nd_alloc + 4;
    if (static_cast<std::size_t>(needed_lq_end) >= store.size_words()) {
      break;
    }

    store.set_lq(l, 0);
    store.set_lq(l + 1, 0);
    store.set_lq(l + 2, 0);
    store.set_iq(l - 5, b.idn);
    store.set_iq(l - 4, b.idh_raw);
    store.set_iq(l - 3, nl_alloc);
    store.set_iq(l - 2, std::min(b.ns, nl_alloc));
    store.set_iq(l - 1, nd_alloc);
    store.set_iq(l, b.status);

    for (std::int32_t j = 0; j < nl_alloc; ++j) {
      const std::int32_t lv = (static_cast<std::size_t>(j) < b.links.size()) ? b.links[j] : 0;
      store.set_lq(l - (j + 1), lv);
    }
    for (std::int32_t k = 0; k < nd_alloc; ++k) {
      store.set_iq(l + 1 + k, b.iwords[static_cast<std::size_t>(k)]);
    }

    if (auto it = tails.find(b.idh); it != tails.end()) {
      store.set_lq(it->second, l);
      it->second = l;
    } else {
      heads[b.idh] = l;
      tails[b.idh] = l;
    }
    cursor = needed_lq_end + 8;
  }
  return heads;
}

FzParseResult parse_fz_stream(std::istream& in, std::size_t max_events) {
  FzParseResult res;
  std::vector<unsigned char> buf(
      (std::istreambuf_iterator<char>(in)),
      std::istreambuf_iterator<char>());
  res.file_bytes = buf.size();
  if (buf.size() < 32) {
    return res;
  }

  bool big_endian = true;
  std::size_t first_off = buf.size();
  for (std::size_t off = 0; off + 32 <= buf.size(); off += 4) {
    if (match_sync(buf.data() + off, true)) {
      big_endian = true;
      first_off = off;
      break;
    }
    if (match_sync(buf.data() + off, false)) {
      big_endian = false;
      first_off = off;
      break;
    }
  }
  if (first_off >= buf.size()) {
    return res;
  }

  res.detected_zebra_fz = true;
  res.is_big_endian = big_endian;
  auto rd = big_endian ? rd_u32_be : rd_u32_le;

  const std::uint32_t raw_len = rd(buf.data() + first_off + 16);
  const std::uint32_t phys_words_u = (raw_len & 0x00FFFFFFu);
  res.physical_record_words =
      (phys_words_u >= 16 && phys_words_u <= 65536) ? phys_words_u : 900u;
  const std::size_t phys_words = static_cast<std::size_t>(res.physical_record_words);
  const std::size_t phys_bytes = phys_words * 4u;

  // Unwrap steering headers and collect logical word stream
  std::vector<std::int32_t> logical_words;
  logical_words.reserve(buf.size() / 4u);

  for (std::size_t off = first_off; off + 32 <= buf.size(); off += phys_bytes) {
    const std::size_t avail_words =
        std::min(phys_words, (buf.size() - off) / 4u);
    std::size_t start_word = 0;
    if (avail_words >= kSteeringHeaderWords &&
        match_sync(buf.data() + off, big_endian)) {
      ++res.steering_blocks;
      start_word = kSteeringHeaderWords;
    } else {
      ++res.fast_blocks;
    }
    for (std::size_t w = start_word; w < avail_words; ++w) {
      logical_words.push_back(static_cast<std::int32_t>(rd(buf.data() + off + w * 4u)));
    }
  }

  // First try structured DELPHI++ FzExchangeReader format if block 0 matches its exact pilot layout:
  // [total_len (= 2 + nuh + npay), type (1..6), nuh, npay]
  if (logical_words.size() >= 4 &&
      static_cast<std::uint32_t>(logical_words[1]) >= 1u &&
      static_cast<std::uint32_t>(logical_words[1]) <= 6u &&
      static_cast<std::uint32_t>(logical_words[2]) <= 1024u &&
      static_cast<std::uint32_t>(logical_words[3]) < logical_words.size() &&
      static_cast<std::uint32_t>(logical_words[0]) ==
          2u + static_cast<std::uint32_t>(logical_words[2]) +
               static_cast<std::uint32_t>(logical_words[3])) {
    in.clear();
    in.seekg(static_cast<std::streamoff>(first_off));
    FzExchangeReader reader(in);
    while (auto rec = reader.next_record()) {
      if (rec->type != FzRecordType::EventNormal &&
          rec->type != FzRecordType::EventContinue) {
        continue;
      }
      // Check if payload encodes pack_zebra_bank: [idn, idh, nl, ns, nd, status, links..., data...]
      if (rec->payload.size() >= 6) {
        const std::int32_t nl = rec->payload[2];
        const std::int32_t ns = rec->payload[3];
        const std::int32_t nd = rec->payload[4];
        if (nl >= 0 && ns >= 0 && ns <= nl && nd >= 0 &&
            rec->payload.size() == static_cast<std::size_t>(6 + nl + nd)) {
          if (const auto idh = decode_hollerith_word(
                  static_cast<std::uint32_t>(rec->payload[1]), true)) {
            ParsedZebraBank bank;
            bank.idn = rec->payload[0];
            bank.idh_raw = rec->payload[1];
            bank.idh = *idh;
            bank.nl = nl;
            bank.ns = ns;
            bank.nd = nd;
            bank.status = rec->payload[5];
            std::size_t idx = 6;
            for (std::int32_t j = 0; j < nl; ++j) {
              bank.links.push_back(rec->payload[idx++]);
            }
            for (std::int32_t k = 0; k < nd; ++k) {
              const std::int32_t iw = rec->payload[idx++];
              bank.iwords.push_back(iw);
              bank.fwords.push_back(std::bit_cast<float>(static_cast<std::uint32_t>(iw)));
            }

            ParsedFzEvent ev;
            ev.user_header = rec->user_header;
            if (rec->user_header.size() >= 2) {
              ev.run_number = rec->user_header[0];
              ev.event_number = rec->user_header[1];
            }
            if (*idh == "HEAD") {
              update_event_header_from_head(bank, ev);
            }
            ev.banks.push_back(std::move(bank));
            ++res.total_banks;
            ++res.bank_counts[*idh];
            res.events.push_back(std::move(ev));
            if (max_events > 0 && res.events.size() >= max_events) {
              return res;
            }
            continue;
          }
        }
      }
      // Fallback: scan payload words directly for embedded ZEBRA banks
      scan_banks_in_words(rec->payload, big_endian, max_events, res);
      if (max_events > 0 && res.events.size() >= max_events) {
        return res;
      }
    }
    if (!res.events.empty()) {
      return res;
    }
  }

  // Standard CERNLIB ZEBRA FZ Exchange stream: scan unwrapped physical blocks
  scan_banks_in_words(logical_words, big_endian, max_events, res);
  return res;
}

}  // namespace delphipp::io
