// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

namespace Acts {

// Implementation
template <typename container_t, bool read_only>
SpacePointProxy<container_t, read_only>::SpacePointProxy(
    typename SpacePointProxy<container_t, read_only>::ContainerType& container,
    std::size_t index)
    : m_container(container), m_index(index) {}

template <typename container_t, bool read_only>
template <bool, typename>
inline typename SpacePointProxy<container_t, read_only>::ValueType&
SpacePointProxy<container_t, read_only>::sp() {
  return container().sp(m_index);
}

template <typename container_t, bool read_only>
inline typename SpacePointProxy<container_t, read_only>::ValueType&
SpacePointProxy<container_t, read_only>::sp() const {
  return container().sp(m_index);
}

template <typename container_t, bool read_only>
inline std::size_t SpacePointProxy<container_t, read_only>::index() const {
  return m_index;
}

template <typename container_t, bool read_only>
inline float SpacePointProxy<container_t, read_only>::x() const {
  return container().x(m_index);
}

template <typename container_t, bool read_only>
inline float SpacePointProxy<container_t, read_only>::y() const {
  return container().y(m_index);
}

template <typename container_t, bool read_only>
inline float SpacePointProxy<container_t, read_only>::z() const {
  return container().z(m_index);
}

template <typename container_t, bool read_only>
inline float SpacePointProxy<container_t, read_only>::phi() const {
  return container().phi(m_index);
}

template <typename container_t, bool read_only>
inline float SpacePointProxy<container_t, read_only>::radius() const {
  return container().radius(m_index);
}

template <typename container_t, bool read_only>
inline float SpacePointProxy<container_t, read_only>::varianceR() const {
  return container().varianceR(m_index);
}

template <typename container_t, bool read_only>
inline float SpacePointProxy<container_t, read_only>::varianceZ() const {
  return container().varianceZ(m_index);
}

template <typename container_t, bool read_only>
inline const float& SpacePointProxy<container_t, read_only>::quality() const {
  return container().quality(m_index);
}

template <typename container_t, bool read_only>
inline const float& SpacePointProxy<container_t, read_only>::deltaR() const {
  return container().deltaR(m_index);
}

template <typename container_t, bool read_only>
inline void SpacePointProxy<container_t, read_only>::setQuality(
    const float& value) const {
  // this is a mutable quantity
  container().setQuality(m_index, value);
}

template <typename container_t, bool read_only>
inline void SpacePointProxy<container_t, read_only>::setDeltaR(
    const float& value) const {
  // this is a mutable quantity
  container().setDeltaR(m_index, value);
}

template <typename container_t, bool read_only>
inline const Acts::Vector3&
SpacePointProxy<container_t, read_only>::topStripVector() const {
  return container().topStripVector(m_index);
}

template <typename container_t, bool read_only>
inline const Acts::Vector3&
SpacePointProxy<container_t, read_only>::bottomStripVector() const {
  return container().bottomStripVector(m_index);
}

template <typename container_t, bool read_only>
inline const Acts::Vector3&
SpacePointProxy<container_t, read_only>::stripCenterDistance() const {
  return container().stripCenterDistance(m_index);
}

template <typename container_t, bool read_only>
inline const Acts::Vector3&
SpacePointProxy<container_t, read_only>::topStripCenterPosition() const {
  return container().topStripCenterPosition(m_index);
}
  
template <typename container_t, bool read_only>
inline typename SpacePointProxy<container_t, read_only>::ContainerType&
SpacePointProxy<container_t, read_only>::container() const {
  return *m_container;
}

template <typename container_t, bool read_only>
template <bool, typename>
inline typename SpacePointProxy<container_t, read_only>::ContainerType&
SpacePointProxy<container_t, read_only>::container() {
  return *m_container;
}

}  // namespace Acts
