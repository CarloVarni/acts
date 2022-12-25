// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

namespace Acts {

template <typename external_space_point_t>
CandidatesForMiddleSp<external_space_point_t>::CandidatesForMiddleSp() {}

template <typename external_space_point_t>
inline void CandidatesForMiddleSp<external_space_point_t>::setMaxElements(
    std::size_t n_low, std::size_t n_high) {
  m_max_size_high = n_high;
  m_max_size_low = n_low;

  // protection against default numbers
  // it may cause std::bad_alloc if we don't protect
  if (n_high == std::numeric_limits<std::size_t>::max() or
      n_low == std::numeric_limits<std::size_t>::max()) {
    return;
  }

  m_storage_high.reserve(n_high);
  m_storage_low.reserve(n_low);
}

template <typename external_space_point_t>
inline void CandidatesForMiddleSp<external_space_point_t>::pop(
    std::vector<value_type>& storage, std::size_t& current_size) {
  storage[0] = storage[current_size - 1];
  bubbledw(storage, 0, --current_size);
}

template <typename external_space_point_t>
inline bool CandidatesForMiddleSp<external_space_point_t>::exists(
    const std::size_t& n, const std::size_t& max_size) const {
  return n < max_size;
}

template <typename external_space_point_t>
inline float CandidatesForMiddleSp<external_space_point_t>::weight(
    const std::vector<value_type>& storage, std::size_t n) const {
  return storage[n].weight;
}

template <typename external_space_point_t>
inline void CandidatesForMiddleSp<external_space_point_t>::clear() {
  // do not clear max size, this is set only once
  m_n_high = 0;
  m_n_low = 0;
  // no need to clean storage
}

template <typename external_space_point_t>
void CandidatesForMiddleSp<external_space_point_t>::push(
   external_space_point_t& SpB, external_space_point_t& SpM, external_space_point_t& SpT,
    float weight, float zOrigin, bool isQuality) {
  if (isQuality) {
    return push(m_storage_high, m_n_high, m_max_size_high, SpB, SpM, SpT,
                weight, zOrigin, isQuality);
  }
  return push(m_storage_low, m_n_low, m_max_size_low, SpB, SpM, SpT,
              weight, zOrigin, isQuality);
}

template <typename external_space_point_t>
void CandidatesForMiddleSp<external_space_point_t>::push(
    std::vector<value_type>& storage, std::size_t& n, const std::size_t& n_max,
    external_space_point_t& SpB, external_space_point_t& SpM, external_space_point_t& SpT, 
    float weight, float zOrigin, bool isQuality) {
  if (n_max == 0) {
    return;
  }

  // if there is still space, add anything
  if (n < n_max) {
    addToCollection(storage, n, n_max,
                    value_type(SpB, SpM, SpT, weight, zOrigin, isQuality));
    return;
  }

  // if no space, replace one if quality is enough
  // compare to element with lower weight
  const auto& lower_weight = weight(storage, 0);
  if (weight <= lower_weight) {
    return;
  }

  // remove element with lower weight and add this one
  pop(storage, n);
  addToCollection(storage, n, n_max,
                  value_type(SpB, SpM, SpT, weight, zOrigin, isQuality));
}

template <typename external_space_point_t>
void CandidatesForMiddleSp<external_space_point_t>::addToCollection(
    std::vector<value_type>& storage, std::size_t& n, const std::size_t& n_max,
    value_type&& element) {
  // adds elements to the end of the collection
  // function called when space in storage is not full
  if (storage.size() == n_max) {
    storage[n] = std::move(element);
  } else {
    storage.push_back(std::move(element));
  }
  bubbleup(storage, n++);
}

template <typename external_space_point_t>
void CandidatesForMiddleSp<external_space_point_t>::bubbledw(
    std::vector<value_type>& storage, std::size_t n, std::size_t actual_size) {
  while (n < actual_size) {
    // left child : 2 * n + 1
    // right child: 2 * n + 2
    float current = weight(storage, n);
    std::size_t left_child = 2 * n + 1;
    std::size_t right_child = 2 * n + 2;

    // no left child, we do nothing
    if (not exists(left_child, actual_size)) {
      break;
    }

    float weight_left_child = weight(storage, left_child);

    std::size_t selected_child = left_child;
    float selected_weight = weight_left_child;

    if (exists(right_child, actual_size)) {
      float weight_right_child = weight(storage, right_child);
      if (weight_right_child <= weight_left_child) {
        selected_child = right_child;
        selected_weight = weight_right_child;
      }
    }

    if (selected_weight >= current) {
      break;
    }

    std::swap(storage[n], storage[selected_child]);
    n = selected_child;
  }  // while loop
}

template <typename external_space_point_t>
void CandidatesForMiddleSp<external_space_point_t>::bubbleup(
    std::vector<value_type>& storage, std::size_t n) {
  while (n != 0) {
    // parent: (n - 1) / 2;
    // this works because it is an integer operation
    std::size_t parent_idx = (n - 1) / 2;

    float weight_current = weight(storage, n);
    float weight_parent = weight(storage, parent_idx);

    if (weight_parent <= weight_current) {
      break;
    }

    std::swap(storage[n], storage[parent_idx]);
    n = parent_idx;
  }
}

template <typename external_space_point_t>
std::vector<typename CandidatesForMiddleSp<external_space_point_t>::value_type>
CandidatesForMiddleSp<external_space_point_t>::storage() const {
  // this will retrieve the entire storage, first high and then low quality
  // the resulting vector is not sorted!
  std::vector<value_type> output;
  output.reserve(m_n_high + m_n_low);

  for (std::size_t idx(0); idx < m_n_high; ++idx) {
    output.push_back(std::move(m_storage_high[idx]));
  }

  for (std::size_t idx(0); idx < m_n_low; ++idx) {
    output.push_back(std::move(m_storage_low[idx]));
  }

  return output;
}

template <typename external_space_point_t>
bool CandidatesForMiddleSp<external_space_point_t>::greaterSort(
    const value_type& i1, const value_type& i2) {
  if (i1.weight != i2.weight) {
    return i1.weight > i2.weight;
  }

  // This is for the case when the weights from different seeds
  // are same. This makes cpu & cuda results same

  const auto& bottom_l1 = i1.bottom;
  const auto& medium_l1 = i1.middle;
  const auto& top_l1 = i1.top;

  const auto& bottom_l2 = i2.bottom;
  const auto& top_l2 = i2.top;

  // medium is the same for all candidates
  float sum_medium =
      medium_l1->y() * medium_l1->y() + medium_l1->z() * medium_l1->z();

  float seed1_sum = sum_medium;
  float seed2_sum = sum_medium;

  seed1_sum +=
      bottom_l1->y() * bottom_l1->y() + bottom_l1->z() * bottom_l1->z();
  seed1_sum += top_l1->y() * top_l1->y() + top_l1->z() * top_l1->z();

  seed2_sum +=
      bottom_l2->y() * bottom_l2->y() + bottom_l2->z() * bottom_l2->z();
  seed2_sum += top_l2->y() * top_l2->y() + top_l2->z() * top_l2->z();

  return seed1_sum > seed2_sum;
}

}  // namespace Acts
