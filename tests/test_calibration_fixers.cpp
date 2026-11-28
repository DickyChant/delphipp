// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <vector>

#include "delphipp/db/condition_db.hpp"
#include "delphipp/magnetic_field.hpp"
#include "delphipp/reco/calibration/calo_muon_fixers.hpp"
#include "delphipp/reco/calibration/global_alignment.hpp"
#include "delphipp/reco/calibration/vd_alignment.hpp"

namespace {

void require(bool cond, const char* expr) {
  if (!cond) {
    std::cerr << "Check failed: " << expr << '\n';
    std::abort();
  }
}

}  // namespace

#define CHECK(expr) require((expr), #expr)

int main() {
  using namespace delphipp;
  using namespace delphipp::reco::calibration;

  // 1. 6-DOF VD alignment round-trip via delphipp::dtrans
  {
    const std::vector<VdLadderSurvey> surveys = {
        VdLadderSurvey{
            .layer = 1,
            .module = 105,
            .translation_cm = {6.312f, -0.045f, 0.180f},
            .rotation_rad = {0.012f, -0.008f, 0.2618f},
            .humidity_bow_cm = 0.0025f,
        },
    };
    const Vec3 local_pt{0.42f, 0.01f, -4.75f};
    const Vec3 global_pt = align_vd_point(local_pt, 1, 105, surveys, true);
    const Vec3 roundtrip = align_vd_point(global_pt, 1, 105, surveys, false);
    CHECK(std::fabs(roundtrip[0] - local_pt[0]) < 1e-5f);
    CHECK(std::fabs(roundtrip[1] - local_pt[1]) < 1e-5f);
    CHECK(std::fabs(roundtrip[2] - local_pt[2]) < 1e-5f);

    // Unmatched module returns input unchanged
    const Vec3 unmatched = align_vd_point(local_pt, 2, 201, surveys, true);
    CHECK(unmatched == local_pt);
  }

  // 2. Iterative beamspot fit rejecting outlier secondary vertices
  {
    std::vector<Vec3> vertices;
    // 25 primary vertices tightly clustered around (0.12, -0.05, -0.70) cm
    for (int i = 0; i < 25; ++i) {
      const float dx = static_cast<float>((i % 5) - 2) * 0.006f;
      const float dy = static_cast<float>(((i / 5) % 5) - 2) * 0.001f;
      const float dz = static_cast<float>((i % 7) - 3) * 0.15f;
      vertices.push_back({0.12f + dx, -0.05f + dy, -0.70f + dz});
    }
    // Displaced secondary vertices (heavy-flavor / hadronic interaction outliers)
    vertices.push_back({0.95f, 0.60f, 1.40f});
    vertices.push_back({-0.80f, -0.95f, -3.50f});
    vertices.push_back({0.45f, -0.45f, 0.25f});
    vertices.push_back({0.12f, 0.35f, -0.70f});

    const auto fit = fit_fill_beamspot(vertices, 4701);
    CHECK(fit.fill_number == 4701);
    CHECK(fit.input_count == 29U);
    CHECK(fit.accepted_count == 25U);
    CHECK(fit.converged);
    CHECK(fit.quality_flag == 0);
    CHECK(std::fabs(fit.position_cm[0] - 0.12f) < 0.002f);
    CHECK(std::fabs(fit.position_cm[1] - (-0.05f)) < 0.001f);
    CHECK(std::fabs(fit.position_cm[2] - (-0.70f)) < 0.05f);
    CHECK(fit.beam_size_cm[0] < 0.02f);
    CHECK(fit.beam_size_cm[1] < 0.005f);
    CHECK(fit.covariance_cm2[0] > 0.0f);
  }

  // 3. VD pedestal/gain calibration & hot/dead channel masking
  {
    const VdChannelCalibration good_ch{
        .layer = 1,
        .module = 3,
        .strip = 100,
        .pedestal = 20.0f,
        .relative_gain = 1.25f,
        .noise = 2.0f,
        .status = VdChannelStatus::Good,
    };
    const VdChannelCalibration dead_ch{
        .layer = 1,
        .module = 3,
        .strip = 101,
        .pedestal = 20.0f,
        .relative_gain = 1.0f,
        .noise = 2.0f,
        .status = VdChannelStatus::Dead,
    };
    const VdChannelCalibration hot_ch{
        .layer = 1,
        .module = 3,
        .strip = 102,
        .pedestal = 20.0f,
        .relative_gain = 1.0f,
        .noise = 8.0f,
        .status = VdChannelStatus::HotNoisy,
    };
    const auto g_sig = apply_vd_pedestal_gain(45.0f, good_ch);
    CHECK(g_sig.has_value());
    CHECK(std::fabs(*g_sig - 20.0f) < 1e-4f);
    CHECK(!apply_vd_pedestal_gain(45.0f, dead_ch).has_value());
    CHECK(!apply_vd_pedestal_gain(45.0f, hot_ch).has_value());

    const std::vector<VdRawStripInput> raw = {
        {.layer = 1, .module = 3, .strip = 100, .raw_adc = 45.0f, .is_z_side = false},
        {.layer = 1, .module = 3, .strip = 101, .raw_adc = 99.0f, .is_z_side = false},
        {.layer = 1, .module = 3, .strip = 102, .raw_adc = 99.0f, .is_z_side = false},
    };
    const std::vector<VdChannelCalibration> table = {good_ch, dead_ch, hot_ch};
    const auto cal_strips = apply_vd_pedestal_gain(raw, table);
    CHECK(cal_strips.size() == 1U);
    CHECK(cal_strips[0].strip == 100);
    CHECK(std::fabs(cal_strips[0].sn_ratio - 10.0f) < 1e-4f);
  }

  // 4. Global detector alignment round-trip (glbali)
  {
    GlobalDetectorAlignment glb;
    glb.set_transform(DetectorId::Tpc, {0.02f, -0.03f, 0.15f}, {0.001f, -0.002f, 0.003f});
    glb.set_transform(static_cast<std::int32_t>(DetectorId::Fcb),
                      {0.10f, 0.05f, 275.2f},
                      {12.0f, 0.15f, -12.0f},
                      true);
    const Vec3 pt{60.0f, 25.0f, 45.0f};
    const Vec3 pt_glb = apply_global_transform(pt, static_cast<std::int32_t>(DetectorId::Tpc), glb);
    const Vec3 pt_back = glb.apply_global_transform(pt_glb, DetectorId::Tpc, true);
    CHECK(std::fabs(pt_back[0] - pt[0]) < 1e-4f);
    CHECK(std::fabs(pt_back[1] - pt[1]) < 1e-4f);
    CHECK(std::fabs(pt_back[2] - pt[2]) < 1e-4f);

    const Vec3 fcb_pt{30.0f, -40.0f, 0.0f};
    const Vec3 fcb_glb = glb.apply_global_transform(fcb_pt, DetectorId::Fcb, false);
    const Vec3 fcb_back = glb.apply_global_transform(fcb_glb, DetectorId::Fcb, true);
    CHECK(std::fabs(fcb_back[0] - fcb_pt[0]) < 1e-4f);
    CHECK(std::fabs(fcb_back[1] - fcb_pt[1]) < 1e-4f);
    CHECK(std::fabs(fcb_back[2] - fcb_pt[2]) < 1e-4f);
  }

  // 5. 94D TPC ABST0 shift (-232 um difference between 1.147 and 1.10208) + ExB + B-map
  {
    db::ConditionDatabase db;
    std::istringstream asc(R"(
*DBRECORD
CALB   GLOB$TPC*.GLOB$ENP0.W
  940101,000001,990101,120000
*DBFIELD SLOW
!94C reprocessing: absolute T0 1.147 - 0.04492 time bins (TPC z +-232 um)
 3,SLOW,14., 1.10208
*DBRECORD
CALB   GLOB$TPC*.GLOB$ENP1.W
  940101,000001,990101,120000
*DBFIELD SLOW
 3,SLOW,14., 1.10208
*DBRECORD
CALB   GLOB$TPC*.GLOB$ENP0.W
  950101,000001,990101,120000
*DBFIELD SLOW
 3,SLOW,14., 1.147
*DBRECORD
CALB   GLOB$TPC*.GLOB$ENP1.W
  950101,000001,990101,120000
*DBFIELD SLOW
 3,SLOW,14., 1.147
)");
    CHECK(db.load_ascii(asc) == 4U);

    MagneticField bfield;
    bfield.set_barrel_layers({
        FieldLayerTable{
            .coord_min_cm = 20.0f,
            .coord_max_cm = 120.0f,
            .bz_grid_cm = {0.0f, 150.0f},
            .bz_tesla = {1.20f, 1.195f},
            .br_grid_cm = {0.0f, 150.0f},
            .br_tesla = {0.0f, 0.006f},
        },
    });

    const Vec3 tpc_pos_p{60.0f, 0.0f, 50.0f};   // ENP0 (+z)
    const Vec3 tpc_pos_n{60.0f, 0.0f, -50.0f};  // ENP1 (-z)

    TpcCorrectionConfig cfg_94d;
    cfg_94d.query_time = {940526, 120000};
    TpcCorrectionConfig cfg_95;
    cfg_95.query_time = {950201, 0};

    const auto res_94d_p = correct_tpc_coordinates_detailed(tpc_pos_p, cfg_94d, &db, &bfield);
    const auto res_95_p = correct_tpc_coordinates_detailed(tpc_pos_p, cfg_95, &db, &bfield);
    const auto res_94d_n = correct_tpc_coordinates_detailed(tpc_pos_n, cfg_94d, &db, &bfield);

    CHECK(std::fabs(res_94d_p.applied_abst0_bins - 1.10208f) < 1e-5f);
    CHECK(std::fabs(res_95_p.applied_abst0_bins - 1.147f) < 1e-5f);

    // Exact -232 um (-0.0232 cm) shift on +z endplate between 1.147 and 1.10208
    const float dz_um_p = (res_94d_p.corrected_xyz_cm[2] - res_95_p.corrected_xyz_cm[2]) * 10000.0f;
    CHECK(std::fabs(dz_um_p - (-232.0f)) < 1.0f);
    // Symmetric +232 um shift on -z endplate (ENP1)
    const float dz_um_n = (res_94d_n.corrected_xyz_cm[2] - tpc_pos_n[2]) * 10000.0f;
    CHECK(std::fabs(dz_um_n - 232.0f) < 1.0f);

    // ExB radial shift & B-field map correction are active
    CHECK(res_94d_p.delta_r_exb_cm > 0.0f);
    CHECK(std::fabs(res_94d_p.delta_r_bfield_cm) > 0.0f);
    CHECK(std::fabs(res_94d_p.delta_rphi_bfield_cm) > 0.0f);
  }

  // 6. ID, OD, FCA, FCB t0 and vdrift corrections
  {
    const auto id_cfg = default_wire_timing_config(WireDetector::Id, 47014);
    const auto id_res = apply_wire_t0_and_vdrift(
        WireTimingInput{
            .detector = WireDetector::Id,
            .raw_time_ns = id_cfg.t0_ns + 100.0f,
            .wire_coord_cm = 2.0f,
            .lr_sign = +1,
        },
        id_cfg);
    CHECK(id_res.valid);
    CHECK(std::fabs(id_res.corrected_time_ns - 100.0f) < 1e-4f);
    CHECK(std::fabs(id_res.drift_distance_cm - 0.40f) < 1e-4f);
    CHECK(std::fabs(id_res.coordinate_cm - 2.40f) < 1e-4f);

    const auto od_cfg = default_wire_timing_config(WireDetector::Od, 47014);
    const auto od_res = apply_wire_t0_and_vdrift(
        WireTimingInput{
            .detector = WireDetector::Od,
            .raw_time_ns = od_cfg.t0_ns + 10.0f,
            .secondary_time_ns = od_cfg.t0_secondary_ns + 14.0f,
            .wire_coord_cm = 40.0f,
            .lr_sign = +1,
            .has_dual_end_timing = true,
        },
        od_cfg);
    CHECK(od_res.valid);
    CHECK(std::fabs(od_res.z_cm - 40.0f) < 1e-3f);

    const auto fca_cfg = default_wire_timing_config(WireDetector::Fca, 47014);
    const auto fca_res = apply_wire_t0_and_vdrift(
        WireTimingInput{
            .detector = WireDetector::Fca,
            .raw_time_ns = fca_cfg.t0_ns + 40.0f,
            .wire_coord_cm = 10.0f,
            .z_cm = 160.0f,
            .lr_sign = -1,
        },
        fca_cfg);
    CHECK(fca_res.valid);
    CHECK(std::fabs(fca_res.drift_distance_cm - 0.20f) < 1e-4f);
    CHECK(std::fabs(fca_res.coordinate_cm - 9.80f) < 1e-4f);

    const auto fcb_cfg = default_wire_timing_config(WireDetector::Fcb, 47014);
    const auto fcb_bad = apply_wire_t0_and_vdrift(
        WireTimingInput{
            .detector = WireDetector::Fcb,
            .raw_time_ns = fcb_cfg.t0_ns - 5.0f,
            .wire_coord_cm = 0.0f,
        },
        fcb_cfg);
    CHECK(!fcb_bad.valid);
    CHECK(fcb_bad.drift_distance_cm == 0.0f);
  }

  // 7. HCAL P/T/isobutane correction & dead/hot tower masking
  {
    const HcalFixConfig hcal_nominal{};
    const auto e_nom = correct_hcal_tower_energy(10.0f, hcal_nominal);
    CHECK(e_nom.has_value());
    CHECK(std::fabs(*e_nom - 10.0f) < 1e-5f);

    HcalFixConfig hcal_pt;
    hcal_pt.pressure_mbar = 985.0f;  // +10 mbar -> * (1 + 0.08056)
    hcal_pt.temperature_c = 22.0f;   // +2 C -> * (1 - 0.060)
    const auto e_pt = correct_hcal_tower_energy(10.0f, hcal_pt);
    CHECK(e_pt.has_value());
    const float expected_e = 10.0f * (1.0f + 0.008056f * 10.0f) * (1.0f - 0.030f * 2.0f);
    CHECK(std::fabs(*e_pt - expected_e) < 1e-4f);

    const HcalTowerCalibration dead_tower{.status = HcalTowerStatus::Dead};
    const HcalTowerCalibration hot_tower{.status = HcalTowerStatus::Hot};
    CHECK(!correct_hcal_tower_energy(10.0f, hcal_nominal, dead_tower).has_value());
    CHECK(!correct_hcal_tower_energy(10.0f, hcal_nominal, hot_tower).has_value());
  }

  // 8. Muon chamber offset & efficiency fixer
  {
    const MuonChamberHit mu_hit{
        .subsystem = MuonSubsystem::Mub,
        .chamber_id = 42,
        .layer = 1,
        .rphi_cm = 12.0f,
        .z_cm = -80.0f,
        .sigma_rphi_cm = 0.2f,
        .sigma_z_cm = 1.0f,
        .weight = 1.0f,
    };
    const MuonChamberFix mu_fix{
        .subsystem = MuonSubsystem::Mub,
        .chamber_id = 42,
        .layer = 1,
        .rphi_offset_cm = 0.15f,
        .z_offset_cm = -0.40f,
        .anode_error_scale = 1.25f,
        .delay_error_scale = 1.10f,
        .efficiency = 0.90f,
        .is_dead = false,
    };
    const auto accepted = correct_muon_chamber(mu_hit, mu_fix, 0.50f);
    CHECK(accepted.has_value());
    CHECK(std::fabs(accepted->rphi_cm - 12.15f) < 1e-4f);
    CHECK(std::fabs(accepted->z_cm - (-80.40f)) < 1e-4f);
    CHECK(std::fabs(accepted->sigma_rphi_cm - 0.25f) < 1e-4f);
    CHECK(std::fabs(accepted->weight - 0.90f) < 1e-4f);
    CHECK(!correct_muon_chamber(mu_hit, mu_fix, 0.95f).has_value());
  }

  // 9. RICH Gladstone-Dale P/T refractive index scaling (C6F14, C5F12, C4F10)
  {
    const float n_c5_ref = adjusted_rich_refractive_index(
        RichRadiator::GasBarrelC5F12, kDefaultRichRefPressureMbar, kDefaultRichRefTempBarrelK);
    CHECK(std::fabs(n_c5_ref - kNominalIndexC5F12) < 1e-6f);

    // Higher pressure and lower temperature increase gas (n - 1) proportionally to P/T
    const float p_test = 1050.0f;
    const float t_test = 300.0f;
    const float n_c5_pt = adjusted_rich_refractive_index(
        RichRadiator::GasBarrelC5F12, p_test, t_test);
    const float expected_c5 =
        1.0f + (kNominalIndexC5F12 - 1.0f) *
                   (p_test / kDefaultRichRefPressureMbar) *
                   (kDefaultRichRefTempBarrelK / t_test);
    CHECK(std::fabs(n_c5_pt - expected_c5) < 1e-6f);

    const float n_c4_pt = adjusted_rich_refractive_index(
        RichRadiator::GasForwardC4F10, 980.0f, 305.0f);
    const float expected_c4 =
        1.0f + (kNominalIndexC4F10 - 1.0f) *
                   (980.0f / kDefaultRichRefPressureMbar) *
                   (kDefaultRichRefTempForwardK / 305.0f);
    CHECK(std::fabs(n_c4_pt - expected_c4) < 1e-6f);

    const float n_c6_warm = adjusted_rich_refractive_index(
        RichRadiator::LiquidC6F14, kDefaultRichRefPressureMbar, 323.15f);
    CHECK(n_c6_warm < kNominalIndexC6F14);

    const float theta_corr = correct_cherenkov_angle(0.060f, kNominalIndexC5F12, n_c5_pt);
    CHECK(theta_corr > 0.060f);
  }

  std::cout << "test_calibration_fixers: PASS\n";
  return 0;
}
