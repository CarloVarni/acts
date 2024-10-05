// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <array>
#include <limits>

namespace Acts {

enum class OwningPolicy : bool {Owner, Viewer};
  
template <typename external_spacepoint_t, std::size_t N = 3ul, OwningPolicy owningPolicy = OwningPolicy::Viewer>
  requires(N >= 3)
class Seed {
 public:
  static constexpr std::size_t DIM = N;
  static constexpr OwningPolicy owning_policy = owningPolicy;
  using value_type = typename std::conditional<owningPolicy == OwningPolicy::Viewer,
					       const external_spacepoint_t*,
					       const external_spacepoint_t>::type;
  using element_type = external_spacepoint_t;

  Seed() = delete;

  // lvalues                                                                                                                                                                                  
  template <typename ... args_t>
  requires(sizeof...(args_t) == N)
  Seed(const args_t&... points)
    requires(owningPolicy == OwningPolicy::Viewer);
  
  // rvalues                                                                                                                                                                                  
  template <typename ... args_t>
  requires(sizeof...(args_t) == N)
  Seed(const args_t&&... points)
    requires(owningPolicy == OwningPolicy::Viewer) = delete;
  
  // any value                                                                                                                                                                                
  template <typename ... args_t>
  requires(sizeof...(args_t) == N)
  Seed(args_t&& ... points)
    requires(owningPolicy == OwningPolicy::Owner);

  void setVertexZ(float vertex);
  void setQuality(float seedQuality);

  const std::array<value_type, N>& sp() const;
  float z() const;
  float seedQuality() const;

 private:
  std::array<value_type, N> m_spacepoints{};
  float m_vertexZ{0.f};
  float m_seedQuality{-std::numeric_limits<float>::infinity()};
};

}  // namespace Acts

#include "Acts/EventData/Seed.ipp"
