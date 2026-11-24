// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "delphipp/reco/rich.hpp"

#include <cmath>
#include <limits>

namespace delphipp::reco::rich {

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
