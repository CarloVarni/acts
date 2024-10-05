// -*- C++ -*-
// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

namespace Acts {

  template <typename external_spacepoint_t, std::size_t N, OwningPolicy owningPolicy>
  requires(N >= 3)    
    template <typename ... args_t>
    requires(sizeof...(args_t) == N)
    Seed<external_spacepoint_t, N, owningPolicy>::Seed(const args_t&... points)
    requires(owningPolicy == OwningPolicy::Viewer)
    : m_spacepoints({&points...}) {}
  
  template <typename external_spacepoint_t, std::size_t N, OwningPolicy owningPolicy>
  requires(N >= 3)
    template <typename ... args_t>
    requires(sizeof...(args_t) == N)
    Seed<external_spacepoint_t, N, owningPolicy>::Seed(args_t&& ... points)
    requires(owningPolicy == OwningPolicy::Owner)
    : m_spacepoints({std::forward<args_t>(points)...})
  {}
  
  template <typename external_spacepoint_t, std::size_t N, OwningPolicy owningPolicy>
  requires(N >= 3)
    void Seed<external_spacepoint_t, N, owningPolicy>::setVertexZ(float vertex) {
  m_vertexZ = vertex;
}

  template <typename external_spacepoint_t, std::size_t N, OwningPolicy owningPolicy>
  requires(N >= 3)
    void Seed<external_spacepoint_t, N, owningPolicy>::setQuality(float seedQuality) {
  m_seedQuality = seedQuality;
}

  template <typename external_spacepoint_t, std::size_t N, OwningPolicy owningPolicy>
  requires(N >= 3)
    const std::array<typename Seed<external_spacepoint_t, N, owningPolicy>::value_type, N>&
    Seed<external_spacepoint_t, N, owningPolicy>::sp() const {
  return m_spacepoints;
}

  template <typename external_spacepoint_t, std::size_t N, OwningPolicy owningPolicy>
  requires(N >= 3)
    float Seed<external_spacepoint_t, N, owningPolicy>::z() const {
  return m_vertexZ;
}

  template <typename external_spacepoint_t, std::size_t N, OwningPolicy owningPolicy>
  requires(N >= 3)
    float Seed<external_spacepoint_t, N, owningPolicy>::seedQuality() const {
  return m_seedQuality;
}

}  // namespace Acts
