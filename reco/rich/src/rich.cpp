// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "delphipp/reco/rich.hpp"

#include <cmath>
#include <limits>
#include <numbers>

namespace delphipp::reco::rich {

std::vector<RichPhotonHit> unpack_raw_blocklet(std::span<const std::int32_t> words,
                                               std::int32_t detector_id) {
  std::vector<RichPhotonHit> hits;
  if (words.size() < 4) {
    return hits;
  }
  const std::size_t start =
      (words[0] >= 4 && static_cast<std::size_t>(words[0]) <= words.size()) ? 4u : 0u;
  constexpr float kTwoPi = 2.0f * std::numbers::pi_v<float>;
  for (std::size_t i = start; i < words.size(); ++i) {
    const std::uint32_t u = static_cast<std::uint32_t>(words[i]);
    const std::uint32_t tdc = u & 0x3FFu;
    if (tdc == 0u) {
      continue;
    }
    const std::int32_t sector = static_cast<std::int32_t>((u >> 10) & 0x1Fu) % 24 + 1;
    const float theta_c_rad = 0.02f + static_cast<float>(tdc % 640u) * 0.001f;
    const float phi_rad = (static_cast<float>(sector - 1) + 0.5f) * (kTwoPi / 24.0f);
    const float z_cm = (detector_id == kDetectorIdBrich)
                           ? (static_cast<float>((u >> 16) & 0xFFu) - 128.0f)
                           : (((u >> 15) & 1u) ? -185.0f : 185.0f);
    hits.push_back(RichPhotonHit{
        .detector_id = detector_id,
        .sector = sector,
        .phi_rad = phi_rad,
        .theta_c_rad = theta_c_rad,
        .z_cm = z_cm,
    });
  }
  return hits;
}

std::optional<float> expected_cherenkov_angle(float momentum_gev,
                                              float mass_gev,
                                              float refractive_index) noexcept {
  if (momentum_gev <= 0.0f || refractive_index <= 1.0f) {
    return std::nullopt;
  }
  const double p = static_cast<double>(momentum_gev);
  const double m = static_cast<double>(mass_gev);
  const double beta = p / std::hypot(p, m);
  const double cos_theta = 1.0 / (static_cast<double>(refractive_index) * beta);
  if (cos_theta > 1.0 || cos_theta < -1.0) {
    return std::nullopt;
  }
  return static_cast<float>(std::acos(cos_theta));
}

RichLikelihoodResult evaluate_hypotheses(const RichRingMeasurement& ring) noexcept {
  RichLikelihoodResult res{};
  float best_chi2 = std::numeric_limits<float>::infinity();
  std::size_t best_idx = static_cast<std::size_t>(Species::Pion);
  const float sig = (ring.sigma_theta_rad > 1e-6f) ? ring.sigma_theta_rad : 0.005f;

  for (std::size_t i = 0; i < 5; ++i) {
    const auto exp_th = expected_cherenkov_angle(
        ring.momentum_gev, kSpeciesMassGeV[i], ring.refractive_index);
    if (!exp_th.has_value()) {
      res.expected_theta_rad[i] = -1.0f;
      // If photons were observed but hypothesis is below threshold, penalize strongly
      res.chi2[i] = (ring.num_photons > 0) ? 100.0f : 0.0f;
    } else {
      res.expected_theta_rad[i] = *exp_th;
      const float diff = (ring.measured_theta_rad - *exp_th) / sig;
      res.chi2[i] = diff * diff;
    }
    if (res.chi2[i] < best_chi2) {
      best_chi2 = res.chi2[i];
      best_idx = i;
    }
  }
  res.best_species = static_cast<Species>(best_idx);
  return res;
}

TanagraId write_ter_bank(TanagraStore& store,
                         const RichRingMeasurement& ring,
                         const RichLikelihoodResult& res,
                         std::int32_t detector_id) {
  return store.add_ter(TerBank{
      .detector_id = detector_id,
      .submod_id = static_cast<std::int32_t>(res.best_species),
      .descriptor = 128,
      .fwords = {ring.measured_theta_rad, ring.sigma_theta_rad,
                 res.chi2[0], res.chi2[1], res.chi2[2], res.chi2[3], res.chi2[4]},
  });
}

}  // namespace delphipp::reco::rich
