// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <sstream>
#include <vector>

#include "delphipp/db/kapack_reader.hpp"
#include "delphipp/helix.hpp"
#include "delphipp/io/raw_blocklet.hpp"
#include "delphipp/magnetic_field.hpp"
#include "delphipp/reco/kalman_fitter.hpp"
#include "delphipp/zebra_store.hpp"

namespace {

void write_le32(std::vector<unsigned char>& buf, std::size_t one_based_word,
                std::uint32_t v) {
  const std::size_t off = (one_based_word - 1) * 4;
  buf[off + 0] = static_cast<unsigned char>(v & 0xFFu);
  buf[off + 1] = static_cast<unsigned char>((v >> 8) & 0xFFu);
  buf[off + 2] = static_cast<unsigned char>((v >> 16) & 0xFFu);
  buf[off + 3] = static_cast<unsigned char>((v >> 24) & 0xFFu);
}

}  // namespace

int main() {
  using namespace delphipp;

  // 1. KAPACK 2-level B-tree traversal (root index block 1 -> leaf block 2)
  {
    constexpr std::size_t kW = 64;
    std::vector<unsigned char> img(2 * kW * 4, 0);
    // Block 1: level=2 index block, first_entry=25, 1 entry of length 4 at word 25 -> child block 2
    write_le32(img, 1, 1u);
    write_le32(img, 4, 1u);   // root_block = 1
    write_le32(img, 5, 2u);   // tree_level = 2 (index)
    write_le32(img, 6, 25u);  // first_entry
    write_le32(img, 7, 29u);  // next_free_word
    write_le32(img, 13, 2u);  // KAVRSN = 2
    write_le32(img, 14, 2u);  // num_blocks = 2
    write_le32(img, 15, kW);  // block_words
    write_le32(img, 16, 4u);  // max_pack
    write_le32(img, 24, 0xCAFEBABEu);
    // Index entry at word 25: len=4, key_len=2, key_word=9401, child_block=2
    write_le32(img, 25, 4u);
    write_le32(img, 26, 2u);
    write_le32(img, 27, 9401u);
    write_le32(img, 28, 2u);

    // Block 2: level=1 leaf block, first_entry=10, 1 entry of length 5 -> payload {110208, 672}
    write_le32(img, kW + 1, 2u);
    write_le32(img, kW + 5, 1u);   // leaf
    write_le32(img, kW + 6, 10u);  // first_entry
    write_le32(img, kW + 7, 15u);  // next_free_word
    write_le32(img, kW + 10, 5u);  // entry_len = 5
    write_le32(img, kW + 11, 2u);  // key_len = 2 (len + 1 key word)
    write_le32(img, kW + 12, 14u); // key
    write_le32(img, kW + 13, 110208u);
    write_le32(img, kW + 14, 672u);

    std::string raw(reinterpret_cast<const char*>(img.data()), img.size());
    std::istringstream iss(raw, std::ios::in | std::ios::binary);
    db::KapackReader kr(iss);
    const auto leaves = kr.list_leaf_entries();
    assert(leaves.size() == 1);
    assert(leaves[0].key_words.size() == 1 && leaves[0].key_words[0] == 14);
    assert(leaves[0].payload_words.size() == 2 && leaves[0].payload_words[0] == 110208);
  }

  // 2. MagneticField barrel grid interpolation + -z radial sign flip
  {
    MagneticField field;
    field.set_barrel_layers({
        FieldLayerTable{
            .coord_min_cm = 180.0f,
            .coord_max_cm = 220.0f,
            .bz_grid_cm = {0.0f, 100.0f, 200.0f},
            .bz_tesla = {1.20f, 1.18f, 1.10f},
            .br_grid_cm = {0.0f, 100.0f, 200.0f},
            .br_tesla = {0.00f, 0.04f, 0.12f},
        },
    });
    const auto b_pos = field.evaluate({200.0f, 0.0f, 100.0f});
    const auto b_neg = field.evaluate({200.0f, 0.0f, -100.0f});
    assert(std::fabs(b_pos.bz_tesla() - 1.18f) < 1e-3f);
    assert(b_pos.bx > 0.0f && b_neg.bx < 0.0f);
    assert(std::fabs(b_pos.bx + b_neg.bx) < 1e-6f);
  }

  // 3. 5-parameter Kalman track refit across VD + ID + TPC + OD layers
  {
    MagneticField field;
    const float bz = field.evaluate({0.0f, 0.0f, 0.0f}).bz;
    const PerigeeTrack truth{
        .eps = 0.02f,
        .z0 = 1.5f,
        .theta = 1.35f,
        .phi0 = 0.6f,
        .inv_r = bz / 20.0f,  // 20 GeV/c track
    };
    const std::vector<float> radii = {6.3f, 9.0f, 10.9f, 16.0f, 40.0f, 65.0f, 90.0f, 110.0f, 200.0f};
    std::vector<reco::tracking::CylinderMeasurement> meas;
    for (const float r : radii) {
      const auto pt = propagate_to_cylinder(truth, r);
      assert(pt.has_value());
      meas.push_back(reco::tracking::CylinderMeasurement{
          .r_cm = r,
          .phi_rad = std::atan2(pt->xyz[1], pt->xyz[0]),
          .z_cm = pt->xyz[2],
          .sigma_rphi_cm = (r < 12.0f) ? 0.002f : 0.02f,
          .sigma_z_cm = (r < 12.0f) ? 0.005f : 0.10f,
          .has_z = true,
          .material_x0 = 0.005f,
      });
    }

    PerigeeTrack perturbed = truth;
    perturbed.eps += 0.05f;
    perturbed.z0 -= 0.4f;
    perturbed.phi0 += 0.02f;
    perturbed.inv_r *= 1.15f;

    const auto fit = reco::tracking::kalman_refit(perturbed, bz, meas);
    assert(fit.hits_used == static_cast<std::int32_t>(radii.size()));
    assert(std::fabs(fit.track.eps - truth.eps) < 0.02f);
    assert(std::fabs(fit.track.z0 - truth.z0) < 0.05f);
    assert(std::fabs(fit.track.phi0 - truth.phi0) < 0.005f);
    assert(std::fabs(fit.track.inv_r - truth.inv_r) / std::fabs(truth.inv_r) < 0.05f);
  }

  // 4. RawBlocklet chain unpacking
  {
    OwnedZebraStore store(512);
    auto v = store.view();
    v.set_lq(100, 200);
    v.set_iq(100 - 1, 3);
    v.set_iq(100 + 2, (29 << 8) | (1 << 16));
    v.set_lq(200, 0);
    v.set_iq(200 - 1, 3);
    v.set_iq(200 + 2, (14 << 8) | (2 << 16));

    const auto blks = io::read_raw_blocklet_chain(v, 100);
    assert(blks.size() == 2);
    assert(blks[0].detector_id == 29 && blks[0].subdet_id == 1);
    assert(blks[1].detector_id == 14 && blks[1].subdet_id == 2);
  }

  std::cout << "test_kalman_and_plumbing: PASS\n";
  return 0;
}
