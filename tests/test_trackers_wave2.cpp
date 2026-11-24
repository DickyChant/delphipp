// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include <cassert>
#include <cmath>
#include <iostream>
#include <numbers>
#include <vector>

#include "delphipp/reco/fwd_chambers.hpp"
#include "delphipp/reco/id.hpp"
#include "delphipp/reco/od.hpp"
#include "delphipp/reco/tpc.hpp"
#include "delphipp/reco/vd.hpp"
#include "delphipp/tanagra_store.hpp"

int main() {
  using namespace delphipp;

  TanagraStore store;

  // 1. VD silicon strip clustering
  std::vector<reco::vd::VdStripHit> vd_strips = {
      {.layer = 1, .module = 3, .strip = 100, .pulse_height = 3.0f, .noise = 1.0f},
      {.layer = 1, .module = 3, .strip = 101, .pulse_height = 9.0f, .noise = 1.0f},
      {.layer = 1, .module = 3, .strip = 102, .pulse_height = 3.0f, .noise = 1.0f},
  };
  const auto vd_clus = reco::vd::find_clusters(vd_strips);
  assert(vd_clus.size() == 1);
  assert(std::fabs(vd_clus[0].centroid_strip - 101.0f) < 1e-4f);
  assert(reco::vd::write_ter_banks(store, vd_clus).size() == 1);

  // 2. TPC z-reconstruction verifying exact 94D ABST0 shift (~232 um = 0.0232 cm)
  reco::tpc::TpcCalibration cal_94d;
  cal_94d.abst0_bins = reco::tpc::kAbst0_94D;
  reco::tpc::TpcCalibration cal_sim;
  cal_sim.abst0_bins = reco::tpc::kAbst0_Sim;

  const float z_94d = reco::tpc::reconstruct_z(0, 100.0f, cal_94d);
  const float z_sim = reco::tpc::reconstruct_z(0, 100.0f, cal_sim);
  // 94D ABST0 (1.10208) is smaller than sim (1.147) -> longer drift -> hit closer to z=0 by ~222-232 um
  assert(z_sim - z_94d > 0.020f && z_sim - z_94d < 0.025f);

  std::vector<float> dedx_samples = {1.0f, 1.1f, 1.2f, 1.3f, 8.0f};  // 8.0 is Landau tail
  const float dedx = reco::tpc::truncated_mean_dedx(dedx_samples, 0.80f);
  assert(std::fabs(dedx - 1.15f) < 1e-4f);  // mean of {1.0, 1.1, 1.2, 1.3}

  std::vector<reco::tpc::TpcPadHit> tpc_hits = {
      {.endplate = 0, .sector = 2, .pad_row = 5, .r_cm = 60.0f, .r_phi_cm = 12.0f, .time_bin = 100.0f, .charge = 50.0f},
  };
  const auto tpc_pts = reco::tpc::reconstruct_space_points(tpc_hits, cal_94d);
  assert(reco::tpc::write_ter_banks(store, tpc_pts).size() == 1);

  // 3. ID jet chamber segment
  std::vector<reco::id::IdJetHit> id_hits = {
      {.sector = 1, .wire = 1, .r_cm = 13.0f, .drift_time_ns = 100.0f, .lr_sign = +1},
      {.sector = 1, .wire = 12, .r_cm = 18.0f, .drift_time_ns = 100.0f, .lr_sign = +1},
  };
  const auto id_segs = reco::id::reconstruct_segments(id_hits);
  assert(id_segs.size() == 1);
  assert(reco::id::write_ter_banks(store, id_segs).size() == 1);

  // 4. OD plank reconstruction
  std::vector<reco::od::OdTubeHit> od_hits = {
      {.plank = 4, .layer = 1, .r_cm = 198.0f, .r_phi_cm = 40.0f, .t_end_a_ns = 10.0f, .t_end_c_ns = 14.0f},
      {.plank = 4, .layer = 3, .r_cm = 202.0f, .r_phi_cm = 40.2f, .t_end_a_ns = 10.0f, .t_end_c_ns = 14.0f},
  };
  const auto od_els = reco::od::reconstruct_planks(od_hits);
  assert(od_els.size() == 1);
  assert(std::fabs(od_els[0].z_cm - 40.0f) < 1e-3f);
  assert(reco::od::write_ter_banks(store, od_els).size() == 1);

  // 5. FCA / FCB multi-plane least-squares intersection for true (x=30, y=-40)
  const float deg60 = std::numbers::pi_v<float> / 3.0f;
  const float x_true = 30.0f, y_true = -40.0f;
  std::vector<reco::fwd::ForwardWireHit> fca_hits = {
      {.endcap = 1, .z_cm = 160.0f, .alpha_rad = 0.0f, .u_cm = x_true},
      {.endcap = 1, .z_cm = 161.0f, .alpha_rad = deg60, .u_cm = x_true * std::cos(deg60) + y_true * std::sin(deg60)},
      {.endcap = 1, .z_cm = 162.0f, .alpha_rad = -deg60, .u_cm = x_true * std::cos(-deg60) + y_true * std::sin(-deg60)},
  };
  const auto fca_seg = reco::fwd::fit_segment(fca_hits);
  assert(fca_seg.has_value());
  assert(std::fabs(fca_seg->x_cm - x_true) < 1e-3f);
  assert(std::fabs(fca_seg->y_cm - y_true) < 1e-3f);
  std::vector<reco::fwd::ForwardSegment> segs = {*fca_seg};
  assert(reco::fwd::write_ter_banks(store, segs, reco::fwd::kDetectorIdFca).size() == 1);
  assert(reco::fwd::write_ter_banks(store, segs, reco::fwd::kDetectorIdFcb).size() == 1);

  assert(store.ter_count() == 6);
  std::cout << "test_trackers_wave2: PASS\n";
  return 0;
}
