// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <numbers>
#include <optional>
#include <vector>

#include "delphipp/helix.hpp"
#include "delphipp/reco/pid/btag.hpp"
#include "delphipp/reco/pid/electron_id.hpp"
#include "delphipp/reco/pid/hadron_id.hpp"
#include "delphipp/reco/pid/muon_id.hpp"
#include "delphipp/reco/rich.hpp"
#include "delphipp/reco/trigger/trigger.hpp"
#include "delphipp/reco/vertex_fitter.hpp"

int main() {
  using namespace delphipp;
  using namespace delphipp::reco;
  namespace pid = delphipp::reco::pid;
  namespace trig = delphipp::reco::trigger;

  constexpr float kPi = std::numbers::pi_v<float>;

  // =========================================================================
  // 1. Electron ID (elepha32): Z -> e+ e- vs pion / early-showering pion
  // =========================================================================
  {
    const float p_ele = 45.6f;
    const float theta_ele = 1.45f;
    const float e_hpc = 45.1f;
    const auto hpc_prof = pid::expected_hpc_em_profile(e_hpc, theta_ele, 12.0f);
    const auto dedx_e = pid::bethe_bloch_dedx(p_ele, rich::Species::Electron);

    const pid::ElectronIdInput ele_in{
        .momentum_gev = p_ele,
        .charge = -1,
        .theta_rad = theta_ele,
        .ecal_energy_gev = e_hpc,
        .hcal_energy_gev = 0.25f,
        .hpc_layer_energy_gev = hpc_prof,
        .hpc_z_cm = 12.0f,
        .delta_z_cm = 0.08f,
        .delta_phi_pos_rad = 0.002f,
        .delta_phi_dir_rad = 0.002f,
        .tpc_dedx = dedx_e.dedx_mip,
        .tpc_num_wires = 160,
        .tpc_gap_cm = 0.4f,
        .rich_electron_likelihood = 0.92f,
        .rich_pion_likelihood = 0.08f,
        .is_crack = false,
        .is_conversion_veto = false,
    };
    const auto ele_res = pid::identify_electron(ele_in);
    assert(ele_res.passes_track_match);
    assert(std::fabs(ele_res.e_over_p - (e_hpc / p_ele)) < 1e-4f);
    assert(ele_res.hpc_profile_chi2 >= 0.0f && ele_res.hpc_profile_chi2 < 10.0f);
    assert(ele_res.prob_combined > 0.15f);
    assert(ele_res.tag == pid::ElectronTag::Tight);

    // Charged pion at 15 GeV/c with hadronic HCAL energy and flat HPC leakage
    const float p_pi = 15.0f;
    const auto dedx_pi = pid::bethe_bloch_dedx(p_pi, rich::Species::Pion);
    const pid::ElectronIdInput pi_in{
        .momentum_gev = p_pi,
        .charge = +1,
        .theta_rad = 1.45f,
        .ecal_energy_gev = 2.2f,
        .hcal_energy_gev = 11.5f,
        .hpc_layer_energy_gev = {0.1f, 0.2f, 0.2f, 0.3f, 0.3f, 0.3f, 0.3f, 0.3f, 0.2f},
        .hpc_z_cm = 12.0f,
        .delta_z_cm = 0.3f,
        .delta_phi_pos_rad = 0.005f,
        .delta_phi_dir_rad = 0.01f,
        .tpc_dedx = dedx_pi.dedx_mip,
        .tpc_num_wires = 160,
        .tpc_gap_cm = 0.4f,
        .rich_electron_likelihood = 0.04f,
        .rich_pion_likelihood = 0.90f,
        .is_crack = false,
        .is_conversion_veto = false,
    };
    const auto pi_res = pid::identify_electron(pi_in);
    assert(pi_res.tag == pid::ElectronTag::None);
    assert(pi_res.prob_combined < 0.01f);
  }

  // =========================================================================
  // 2. Muon ID (muflag + pxhaid): Z -> mu+ mu- vs punchthrough pion
  // =========================================================================
  {
    const std::array<float, 4> res_cm = {0.10f, -0.15f, 0.12f, -0.08f};
    const std::array<float, 4> err_cm = {0.40f, 0.40f, 0.40f, 0.40f};
    const auto chi2_sum = pid::compute_muon_hit_chi2(res_cm, err_cm, 2.4f, 4);
    assert(chi2_sum.chi2_total_pdf < 2.0f);
    assert(chi2_sum.chi2_extrap_pdf < 3.0f);

    const pid::MuonIdInput mu_barrel{
        .momentum_gev = 45.6f,
        .theta_rad = 1.40f,
        .detector = pid::MuonDetector::Mub,
        .num_layers = 4,
        .hit_pattern = 0b0011011u,  // inner (1,2) + outer (4,5)
        .chi2_total_pdf = chi2_sum.chi2_total_pdf,
        .chi2_extrap_pdf = chi2_sum.chi2_extrap_pdf,
        .hcal_layer_energy_gev = {0.40f, 0.45f, 0.35f, 0.40f},
        .hcal_total_energy_gev = 1.60f,
        .apply_hcal_mip_veto = true,
    };
    const auto mu_res = pid::identify_muon(mu_barrel);
    assert(mu_res.has_outer_hit);
    assert(mu_res.passes_hcal_mip_veto);
    assert(mu_res.hcal_veto_rv < 1.0f);
    assert(mu_res.tag == pid::MuonTag::Tight);

    // Forward muon (MUF) Standard tag
    const pid::MuonIdInput mu_fwd{
        .momentum_gev = 25.0f,
        .theta_rad = 0.35f,
        .detector = pid::MuonDetector::Muf,
        .num_layers = 3,
        .hit_pattern = 0b0111u,  // layers 1,2,3 (layer 3 is outer)
        .chi2_total_pdf = 2.8f,
        .chi2_extrap_pdf = 2.2f,
        .hcal_layer_energy_gev = {0.5f, 0.6f, 0.4f, 0.4f},
        .hcal_total_energy_gev = 1.9f,
        .apply_hcal_mip_veto = true,
    };
    const auto muf_res = pid::identify_muon(mu_fwd);
    assert(muf_res.has_outer_hit);
    assert(muf_res.passes_hcal_mip_veto);
    assert(muf_res.tag == pid::MuonTag::Standard);

    // Punchthrough pion: hits MUB chambers with small chi2 but showers in HCAL
    const pid::MuonIdInput punch_pi{
        .momentum_gev = 20.0f,
        .theta_rad = 1.50f,
        .detector = pid::MuonDetector::Mub,
        .num_layers = 3,
        .hit_pattern = 0b0001011u,
        .chi2_total_pdf = 1.6f,
        .chi2_extrap_pdf = 1.4f,
        .hcal_layer_energy_gev = {2.0f, 5.5f, 7.0f, 4.5f},
        .hcal_total_energy_gev = 19.0f,
        .apply_hcal_mip_veto = true,
    };
    const auto punch_res = pid::identify_muon(punch_pi);
    assert(!punch_res.passes_hcal_mip_veto);
    assert(punch_res.hcal_veto_rv > 1.0f);
    assert(punch_res.tag == pid::MuonTag::VeryLoose);
  }

  // =========================================================================
  // 3. Hadron ID (dedxid + pxrich): Bethe-Bloch & combined K/p/pi separation
  // =========================================================================
  {
    // Verify Bethe-Bloch 1/beta^2 ordering at p = 0.7 GeV/c: p > K > pi
    const auto bb_pi_07 = pid::bethe_bloch_dedx(0.7f, rich::Species::Pion);
    const auto bb_k_07 = pid::bethe_bloch_dedx(0.7f, rich::Species::Kaon);
    const auto bb_p_07 = pid::bethe_bloch_dedx(0.7f, rich::Species::Proton);
    assert(bb_p_07.dedx_mip > bb_k_07.dedx_mip);
    assert(bb_k_07.dedx_mip > bb_pi_07.dedx_mip);
    assert(bb_pi_07.dedx_mip > 0.9f && bb_pi_07.dedx_mip < 1.25f);

    // BBDXID pairwise separation at p = 0.7 GeV/c for a Kaon
    const float p_k_vs_pi = pid::dedx_pairwise_probability(
        0.7f, bb_k_07.dedx_mip, 160, rich::Species::Kaon, rich::Species::Pion);
    assert(p_k_vs_pi > 0.95f);

    // Combined TPC dE/dx + liquid RICH Kaon identification at p = 2.5 GeV/c
    const float p_kaon = 2.5f;
    const auto th_kaon = rich::expected_cherenkov_angle(
        p_kaon,
        rich::kSpeciesMassGeV[static_cast<std::size_t>(rich::Species::Kaon)],
        rich::kRefractiveIndexLiquid);
    assert(th_kaon.has_value());
    const auto bb_k_25 = pid::bethe_bloch_dedx(p_kaon, rich::Species::Kaon);

    const pid::HadronIdInput kaon_in{
        .momentum_gev = p_kaon,
        .tpc_dedx = bb_k_25.dedx_mip,
        .tpc_num_wires = 160,
        .tpc_gap_cm = 0.4f,
        .tpc_edge_region = false,
        .rich_liquid = rich::RichRingMeasurement{
            .momentum_gev = p_kaon,
            .refractive_index = rich::kRefractiveIndexLiquid,
            .measured_theta_rad = *th_kaon,
            .sigma_theta_rad = 0.005f,
            .num_photons = 14,
        },
        .rich_gas = std::nullopt,
        .rich_chi2_override = std::nullopt,
        .gaussianize_dedx_pull = true,
        .is_simulation = false,
    };
    const auto kaon_res = pid::identify_hadron(kaon_in);
    assert(kaon_res.dedx_quality == 4);
    assert(kaon_res.best_species == rich::Species::Kaon);
    assert(kaon_res.prob_kaon_vs_pion > 0.90f);
    assert(kaon_res.kaon_tag == 3);

    // Combined TPC dE/dx + liquid RICH Proton identification at p = 3.0 GeV/c
    const float p_prot = 3.0f;
    const auto th_prot = rich::expected_cherenkov_angle(
        p_prot,
        rich::kSpeciesMassGeV[static_cast<std::size_t>(rich::Species::Proton)],
        rich::kRefractiveIndexLiquid);
    assert(th_prot.has_value());
    const auto bb_p_30 = pid::bethe_bloch_dedx(p_prot, rich::Species::Proton);

    const pid::HadronIdInput prot_in{
        .momentum_gev = p_prot,
        .tpc_dedx = bb_p_30.dedx_mip,
        .tpc_num_wires = 160,
        .tpc_gap_cm = 0.4f,
        .tpc_edge_region = false,
        .rich_liquid = rich::RichRingMeasurement{
            .momentum_gev = p_prot,
            .refractive_index = rich::kRefractiveIndexLiquid,
            .measured_theta_rad = *th_prot,
            .sigma_theta_rad = 0.005f,
            .num_photons = 12,
        },
        .rich_gas = std::nullopt,
        .rich_chi2_override = std::nullopt,
        .gaussianize_dedx_pull = true,
        .is_simulation = false,
    };
    const auto prot_res = pid::identify_hadron(prot_in);
    assert(prot_res.best_species == rich::Species::Proton);
    assert(prot_res.prob_proton_vs_pion > 0.90f);
    assert(prot_res.prob_proton_vs_kaon > 0.85f);
    assert(prot_res.proton_tag == 3);
  }

  // =========================================================================
  // 4. 3D B-tagging (aabtagxx): displaced b-jet tracks vs prompt u-jet tracks
  // =========================================================================
  {
    const tracking::VertexFitResult pv{
        .xyz_cm = {0.0f, 0.0f, 0.0f},
        .sigma_cm = {0.0020f, 0.0020f, 0.0032f},
        .chi2 = 4.0f,
        .ndf = 8,
        .num_tracks = 6,
    };
    const Vec3 thrust{1.0f, 0.0f, 0.0f};
    const Vec3 jet1{1.0f, 0.0f, 0.0f};
    const Vec3 jet2{-1.0f, 0.0f, 0.0f};
    const pid::BTagConfig bcfg{};

    // Prompt u-jet event: 6 tracks with ~1 sigma resolution jitter
    std::vector<pid::BTagTrackInput> prompt_tracks = {
        {.track = {.eps = -0.0020f, .z0 = -0.0030f, .theta = 1.45f, .phi0 = 0.20f, .inv_r = 0.0004f},
         .sigma_rphi_cm = 0.0035f, .sigma_z_cm = 0.0070f, .n_vd_rphi = 3, .n_vd_z = 2, .jet_axis = jet1},
        {.track = {.eps = +0.0025f, .z0 = +0.0040f, .theta = 1.65f, .phi0 = -0.25f, .inv_r = -0.0005f},
         .sigma_rphi_cm = 0.0035f, .sigma_z_cm = 0.0070f, .n_vd_rphi = 3, .n_vd_z = 2, .jet_axis = jet1},
        {.track = {.eps = -0.0018f, .z0 = +0.0035f, .theta = 1.50f, .phi0 = 0.30f, .inv_r = 0.0003f},
         .sigma_rphi_cm = 0.0035f, .sigma_z_cm = 0.0070f, .n_vd_rphi = 3, .n_vd_z = 2, .jet_axis = jet1},
        {.track = {.eps = -0.0022f, .z0 = -0.0045f, .theta = 1.55f, .phi0 = kPi + 0.22f, .inv_r = 0.0004f},
         .sigma_rphi_cm = 0.0035f, .sigma_z_cm = 0.0070f, .n_vd_rphi = 3, .n_vd_z = 2, .jet_axis = jet2},
        {.track = {.eps = +0.0019f, .z0 = +0.0030f, .theta = 1.48f, .phi0 = kPi - 0.18f, .inv_r = -0.0004f},
         .sigma_rphi_cm = 0.0035f, .sigma_z_cm = 0.0070f, .n_vd_rphi = 3, .n_vd_z = 2, .jet_axis = jet2},
        {.track = {.eps = -0.0015f, .z0 = +0.0025f, .theta = 1.62f, .phi0 = kPi + 0.28f, .inv_r = 0.0005f},
         .sigma_rphi_cm = 0.0035f, .sigma_z_cm = 0.0070f, .n_vd_rphi = 3, .n_vd_z = 2, .jet_axis = jet2},
    };
    const auto prompt_res = pid::evaluate_event_btag(prompt_tracks, pv, thrust, bcfg);
    assert(prompt_res.num_selected_tracks == 6);
    assert(prompt_res.p_event_btag > 0.1f);
    assert(prompt_res.p_event_all > 0.1f);

    // Displaced b-jet event: secondary vertices at x = +0.25 cm (jet1) and x = -0.25 cm (jet2)
    auto make_displaced_track = [&](float vx_dec, float phi0, float theta,
                                    float inv_r, const Vec3& jet) {
      const float eps = -vx_dec * std::sin(phi0);
      const float z0 = -vx_dec * std::cos(phi0) / std::tan(theta);
      return pid::BTagTrackInput{
          .track = {.eps = eps, .z0 = z0, .theta = theta, .phi0 = phi0, .inv_r = inv_r},
          .sigma_rphi_cm = 0.0035f,
          .sigma_z_cm = 0.0070f,
          .n_vd_rphi = 3,
          .n_vd_z = 2,
          .jet_axis = jet,
      };
    };

    std::vector<pid::BTagTrackInput> bjet_tracks = {
        make_displaced_track(+0.25f, +0.24f, 1.38f, +0.0004f, jet1),
        make_displaced_track(+0.25f, -0.28f, 1.75f, -0.0005f, jet1),
        make_displaced_track(+0.20f, +0.18f, 1.42f, +0.0003f, jet1),
        make_displaced_track(-0.25f, kPi + 0.22f, 1.40f, +0.0004f, jet2),
        make_displaced_track(-0.25f, kPi - 0.26f, 1.72f, -0.0004f, jet2),
        {.track = {.eps = +0.0010f, .z0 = -0.0020f, .theta = 1.55f, .phi0 = kPi + 0.10f, .inv_r = 0.0004f},
         .sigma_rphi_cm = 0.0035f, .sigma_z_cm = 0.0070f, .n_vd_rphi = 3, .n_vd_z = 2, .jet_axis = jet2},
    };
    const auto bjet_res = pid::evaluate_event_btag(bjet_tracks, pv, thrust, bcfg);
    assert(bjet_res.num_selected_tracks == 6);
    assert(bjet_res.num_positive_tracks >= 5);
    assert(bjet_res.tracks[0].significance_3d > 5.0f);
    assert(bjet_res.tracks[3].significance_3d > 5.0f);
    assert(bjet_res.p_event_btag < 1.0e-3f);
    assert(bjet_res.p_hemisphere_btag[0] < 1.0e-3f);
    assert(bjet_res.p_hemisphere_btag[1] < 1.0e-3f);
  }

  // =========================================================================
  // 5. T1-T4 Trigger (trigger.car & fcabtrig.car)
  // =========================================================================
  {
    // 5a. Hadronic Z -> q qbar
    trig::TriggerDetectorInput had_in{};
    had_in.id_hits = 12;
    had_in.od_hits = 10;
    had_in.tof_hits = 6;
    had_in.hpc_clusters = 8;
    had_in.hpc_high_clusters = 3;
    had_in.hpc_energy_gev = 28.0f;
    had_in.tpc_tracks = 18;
    had_in.tpc_back_to_back = true;
    had_in.hcal_clusters = 5;
    had_in.hcal_energy_gev = 24.0f;
    had_in.charged_multiplicity = 18;
    had_in.visible_energy_gev = 82.0f;
    had_in.max_track_p_gev = 19.0f;

    const auto had_dec = trig::evaluate_trigger(had_in);
    assert(had_dec.t1.id_od);
    assert(had_dec.t1.t1_accept);
    assert(had_dec.t2.tpc_confirm);
    assert(had_dec.t2.hcal_confirm);
    assert(had_dec.t2.t2_accept);
    assert(had_dec.t3_accept);
    assert(had_dec.t4_accept);
    assert(trig::has_stream(had_dec.streams, trig::PhysicsStream::HadronicZ));
    assert(had_dec.primary_stream == trig::PhysicsStream::HadronicZ);

    // 5b. Dilepton Z -> mu+ mu-
    trig::TriggerDetectorInput dimu_in{};
    dimu_in.id_hits = 2;
    dimu_in.od_hits = 2;
    dimu_in.tof_hits = 2;
    dimu_in.tof_back_to_back = true;
    dimu_in.mub_hits = 2;
    dimu_in.mub_back_to_back = true;
    dimu_in.tpc_tracks = 2;
    dimu_in.tpc_back_to_back = true;
    dimu_in.hcal_clusters = 2;
    dimu_in.hcal_energy_gev = 3.2f;
    dimu_in.charged_multiplicity = 2;
    dimu_in.visible_energy_gev = 91.2f;
    dimu_in.max_track_p_gev = 45.6f;

    const auto dimu_dec = trig::evaluate_trigger(dimu_in);
    assert(dimu_dec.t1.mub_coinc);
    assert(dimu_dec.t1.t1_accept);
    assert(dimu_dec.t2.mub_confirm);
    assert(dimu_dec.t2.t2_accept);
    assert(dimu_dec.t3_accept);
    assert(dimu_dec.t4_accept);
    assert(trig::has_stream(dimu_dec.streams, trig::PhysicsStream::Dilepton));
    assert(!trig::has_stream(dimu_dec.streams, trig::PhysicsStream::HadronicZ));

    // 5c. Bhabha (barrel HPC + SAT luminosity)
    trig::TriggerDetectorInput bhabha_in{};
    bhabha_in.id_hits = 2;
    bhabha_in.od_hits = 2;
    bhabha_in.hpc_clusters = 2;
    bhabha_in.hpc_high_clusters = 2;
    bhabha_in.hpc_energy_gev = 88.5f;
    bhabha_in.hpc_back_to_back = true;
    bhabha_in.tpc_tracks = 2;
    bhabha_in.tpc_back_to_back = true;
    bhabha_in.hcal_energy_gev = 0.4f;
    bhabha_in.charged_multiplicity = 2;
    bhabha_in.visible_energy_gev = 90.0f;
    bhabha_in.max_track_p_gev = 45.0f;

    const auto bhabha_dec = trig::evaluate_trigger(bhabha_in);
    assert(bhabha_dec.t1.hpc_bb);
    assert(bhabha_dec.t2.hpc_confirm);
    assert(bhabha_dec.t4_accept);
    assert(trig::has_stream(bhabha_dec.streams, trig::PhysicsStream::Bhabha));

    trig::TriggerDetectorInput sat_in{};
    sat_in.sat_bhabha_coincidence = true;
    sat_in.visible_energy_gev = 89.0f;
    const auto sat_dec = trig::evaluate_trigger(sat_in);
    assert(sat_dec.t1.sat_bhabha);
    assert(sat_dec.t2.emf_sat_confirm);
    assert(sat_dec.t4_accept);
    assert(trig::has_stream(sat_dec.streams, trig::PhysicsStream::Bhabha));

    // 5d. GammaGamma (soft two-photon event)
    trig::TriggerDetectorInput gaga_in{};
    gaga_in.id_hits = 3;
    gaga_in.od_hits = 2;
    gaga_in.tof_hits = 2;
    gaga_in.hpc_clusters = 2;
    gaga_in.hpc_energy_gev = 2.4f;
    gaga_in.tpc_tracks = 4;
    gaga_in.charged_multiplicity = 4;
    gaga_in.visible_energy_gev = 6.8f;
    gaga_in.max_track_p_gev = 2.5f;

    const auto gaga_dec = trig::evaluate_trigger(gaga_in);
    assert(gaga_dec.t1.t1_accept);
    assert(gaga_dec.t2.t2_accept);
    assert(gaga_dec.t3_accept);
    assert(gaga_dec.t4_accept);
    assert(trig::has_stream(gaga_dec.streams, trig::PhysicsStream::GammaGamma));

    // 5e. FCA-FCB forward trigger (fcabtrig.car T1 & T2 sum in [93, 100])
    trig::TriggerDetectorInput fcab_in{};
    fcab_in.fcab.forward.fca_hit_groups = 2;
    fcab_in.fcab.forward.fcb_plane_hits = {4, 3, 4};
    fcab_in.fcab.forward.t1b_groups = {{30}, {32}, {34}};  // 30+32+34 = 96 in [93, 100]
    fcab_in.emf_low_fw = true;
    fcab_in.emf_low_bw = true;
    fcab_in.emf_energy_fw_gev = 42.0f;
    fcab_in.emf_energy_bw_gev = 43.0f;
    fcab_in.charged_multiplicity = 2;
    fcab_in.visible_energy_gev = 85.0f;

    const auto fcab_dec = trig::evaluate_trigger(fcab_in);
    assert(fcab_dec.fcab.t1_forward > 0);
    assert(fcab_dec.fcab.t2_forward == 1);
    assert(fcab_dec.t1.fwd_maj);
    assert(fcab_dec.t2.fcab_confirm);
    assert(fcab_dec.t4_accept);
    assert(trig::has_stream(fcab_dec.streams, trig::PhysicsStream::Bhabha));

    // 5f. Noise rejection (uncoincident single hit)
    trig::TriggerDetectorInput noise_in{};
    noise_in.id_hits = 1;
    const auto noise_dec = trig::evaluate_trigger(noise_in);
    assert(!noise_dec.t1.t1_accept);
    assert(!noise_dec.t2.t2_accept);
    assert(!noise_dec.t3_accept);
    assert(!noise_dec.t4_accept);
  }

  std::cout << "test_pid_btag_trigger: PASS\n";
  return 0;
}
