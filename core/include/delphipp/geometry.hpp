// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// DELPHI 3D coordinate/vector transformation matching DDAPP SUBROUTINE DTRANS
// (simana_v94d/src/car/ddapp600.car, lines 19127-19180).

#pragma once

#include <array>
#include <span>

namespace delphipp {

using Vec3 = std::array<float, 3>;
using DtransMatrix = std::array<float, 12>;

// MODE = 0 (point):  vnew[i] = rmatr[i] + sum_j rmatr[3 + 3*i + j] * vold[j]
// MODE = 1 (vector): vnew[i] =            sum_j rmatr[3 + 3*i + j] * vold[j]
[[nodiscard]] inline Vec3 dtrans(const Vec3& vold,
                                 std::span<const float, 12> rmatr,
                                 int mode = 0) noexcept {
  Vec3 vnew{};
  for (std::size_t i = 0; i < 3; ++i) {
    double acc = 0.0;
    for (std::size_t j = 0; j < 3; ++j) {
      acc += static_cast<double>(rmatr[3 + 3 * i + j]) *
             static_cast<double>(vold[j]);
    }
    if (mode != 1) {
      acc += static_cast<double>(rmatr[i]);
    }
    vnew[i] = static_cast<float>(acc);
  }
  return vnew;
}

}  // namespace delphipp
