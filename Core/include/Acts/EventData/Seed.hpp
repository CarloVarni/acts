// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <array>
#include <limits>

namespace Acts {

template <typename external_spacepoint_t, std::size_t N = 3ul>
class Seed {
  static_assert(N >= 3ul, "A seed needs at least 3 space points");

 public:
  using value_type = external_spacepoint_t;
  static constexpr std::size_t DIM = N;

  template <typename... args_t>
  Seed(const args_t&... points)
    requires(sizeof...(points) == N) &&
            (std::same_as<external_spacepoint_t, args_t> && ...);

  void setVertexZ(float vertex);
  void setQuality(float seedQuality);

  const std::array<const external_spacepoint_t*, N>& sp() const;
  float z() const;
  float seedQuality() const;

 private:
  template <typename... args_t>
  void storeValues(const external_spacepoint_t& value, const args_t&... others);

  std::array<const external_spacepoint_t*, N> m_spacepoints{};
  float m_vertexZ{0.f};
  float m_seedQuality{-std::numeric_limits<float>::infinity()};
};

}  // namespace Acts

#include "Acts/EventData/Seed.ipp"
