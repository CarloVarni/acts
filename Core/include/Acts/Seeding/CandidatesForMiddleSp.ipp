// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <iostream>

namespace Acts {

template <typename external_space_point_t>
CandidatesForMiddleSp<external_space_point_t>::CandidatesForMiddleSp()
    : m_SpB(CandidatesForMiddleSp<external_space_point_t>::default_value),
      m_SpM(CandidatesForMiddleSp<external_space_point_t>::default_value) {}

template <typename external_space_point_t>
void CandidatesForMiddleSp<external_space_point_t>::push(
    typename CandidatesForMiddleSp<external_space_point_t>::sp_type& SpT,
    float weight, float zOrigin, bool isQuality) {

  if (isQuality) {
    return push(m_storage_high, m_n_high, m_max_size_high,
    	        m_SpB, m_SpM, SpT, weight, zOrigin, isQuality);
  }
  return push(m_storage_low, m_n_low, m_max_size_low,
              m_SpB, m_SpM, SpT, weight, zOrigin, isQuality);
}

template <typename external_space_point_t>
void CandidatesForMiddleSp<external_space_point_t>::push(
    std::vector<value_type>& storage, std::size_t& n, const std::size_t& n_max,
    sp_type& SpB, sp_type& SpM, sp_type& SpT,
    float weight, float zOrigin, bool isQuality) {

  if (n_max == 0) {
    return;
  }

  // if there is still space, add anything
  if (n < n_max) {
    addToCollection(storage, n, n_max,
    		    std::make_tuple(SpB, SpM, SpT, weight, zOrigin, isQuality) );
    return;
  }

  // if no space, replace one if quality is enough
  // compare to element with lower weight
  const auto& lower_weight = top(storage);
  if (weight <= lower_weight) {
    return;
  }

  // remove element with lower weight and add this one
  pop(storage, n);
  addToCollection(storage, n, n_max,
  		  std::make_tuple(SpB, SpM, SpT, weight, zOrigin, isQuality) );
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
    storage.push_back( std::move(element) ); 
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
   } // while loop

}

template <typename external_space_point_t>
void CandidatesForMiddleSp<external_space_point_t>::bubbleup(
    std::vector<value_type>& storage, std::size_t n) {

  while( n != 0) {
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
void CandidatesForMiddleSp<external_space_point_t>::pop(
    std::vector<value_type>& storage, std::size_t& current_size) {
  storage[0] = storage[current_size - 1];
  --current_size;
  bubbledw(storage, 0, current_size);
}

template <typename external_space_point_t>
std::vector<typename CandidatesForMiddleSp<external_space_point_t>::value_type>
CandidatesForMiddleSp<external_space_point_t>::storage() const {
  // this will retrieve the entire storage, first high and then low quality
  // the resulting vector is not sorted!
  std::vector<value_type> output;
  output.reserve(m_n_high + m_n_low);

  for (std::size_t idx(0); idx < m_n_high; idx++) {
    const auto& [bottom, middle, top, weight, zOrigin, isQuality] =
        m_storage_high[idx];
    output.emplace_back(bottom, middle, top, weight, zOrigin, isQuality);
  }

  for (std::size_t idx(0); idx < m_n_low; idx++) {
    const auto& [bottom, middle, top, weight, zOrigin, isQuality] =
        m_storage_low[idx];
    output.emplace_back(bottom, middle, top, weight, zOrigin, isQuality);
  }

  // sort output according to weight and sps
  // should we collect inputs according to this criterion instead?
  std::sort(output.begin(), output.end(),
            CandidatesForMiddleSp<external_space_point_t>::greaterSort);
  return output;
}

template <typename external_space_point_t>
bool CandidatesForMiddleSp<external_space_point_t>::greaterSort(
    const value_type& i1, const value_type& i2) {
  const auto& [bottom_l1, medium_l1, top_l1, weight_l1, zOrigin_l1,
               isQuality_l1] = i1;
  const auto& [bottom_l2, medium_l2, top_l2, weight_l2, zOrigin_l2,
               isQuality_l2] = i2;

  if (weight_l1 != weight_l2) {
    return weight_l1 > weight_l2;
  }

  // This is for the case when the weights from different seeds
  // are same. This makes cpu & cuda results same

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
