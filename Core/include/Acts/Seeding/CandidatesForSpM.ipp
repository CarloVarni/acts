// This file is part of the Acts project.AA
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

namespace Acts {

  template<typename external_space_point_t>
  CandidatesForSpM<external_space_point_t>::CandidatesForSpM()
    : m_max_size_high(0),
      m_max_size_low(0),
      m_n_high(0),
      m_n_low(0),	
      m_SpB(CandidatesForSpM<external_space_point_t>::default_value),
      m_SpM(CandidatesForSpM<external_space_point_t>::default_value)
  {}

  template<typename external_space_point_t>
  void CandidatesForSpM<external_space_point_t>::push(typename CandidatesForSpM<external_space_point_t>::sp_type& SpT,
			      float weight, float zOrigin,
			      bool isQuality)
  {
    auto& storage = isQuality ? m_storage_high : m_storage_low;
    const std::size_t& current_max_size = isQuality ? m_max_size_high : m_max_size_low;
    std::size_t& current_size = isQuality ? m_n_high : m_n_low;
    
    // if there is still space, add anything
    if (current_size < current_max_size) {
      addToCollection(storage,
		      m_SpB, SpT, weight, zOrigin,
		      isQuality);
      return;
    }

    // if no space, replace one if quality is enough
    // compare to element with lower weight
    const auto& lower_weight = top(storage);
    if (weight <= lower_weight)
      return;
    
    // remove element with lower weight and add this one  
    pop(storage, current_size, current_max_size);
    insertToCollection(storage,
    		       m_SpB, SpT, weight, zOrigin,
		       isQuality);
  }

  template<typename external_space_point_t>
  void CandidatesForSpM<external_space_point_t>::addToCollection(std::vector< typename CandidatesForSpM<external_space_point_t>::value_type >& storage,
       					 typename CandidatesForSpM<external_space_point_t>::sp_type& SpB,
					 typename CandidatesForSpM<external_space_point_t>::sp_type& SpT,
					 float weight, float zOrigin,
					 bool isQuality)
  {
    // adds elements to the end of the collection
    // function called when space in storage is not full
    auto toAdd = std::make_tuple(SpB, SpT, weight, zOrigin);
    storage.push_back( toAdd );
    std::size_t added_index = isQuality ? m_n_high++ : m_n_low++;
    bubbleup(storage, added_index);
  }  

  template<typename external_space_point_t>
  void CandidatesForSpM<external_space_point_t>::insertToCollection(std::vector< typename CandidatesForSpM<external_space_point_t>::value_type >& storage,
                                            typename CandidatesForSpM<external_space_point_t>::sp_type& SpB,
					    typename CandidatesForSpM<external_space_point_t>::sp_type& SpT,
					    float weight, float zOrigin,
					    bool isQuality)
  {
    auto toAdd = std::make_tuple(SpB, SpT, weight, zOrigin);
    std::size_t added_index = isQuality ? m_n_high++ : m_n_low++;
    storage[added_index] = toAdd;
    bubbleup(storage, added_index);
  }

  template<typename external_space_point_t>
  void CandidatesForSpM<external_space_point_t>::bubbledw(std::vector< value_type >& storage,
       							  std::size_t n,
							  std::size_t max_size)
  {
    // left child : 2 * n + 1
    // right child: 2 * n + 2
    float current = weight(storage, n);
    std::size_t left_child = 2 * n + 1;
    std::size_t right_child = 2 * n + 2;
    
    // no left child, we stop
    if (not exists(left_child, max_size)) return;
    // no right child, left wins
    if (not exists(right_child, max_size)) {
      if (weight(storage, left_child) < current) {
	std::swap(storage[n], storage[left_child]);
	return bubbledw(storage, left_child, max_size);
      }
    }
    // both childs
    // left is smaller
    if (weight(storage, left_child) < weight(storage, right_child)) {
      if (weight(storage, left_child) < current) {
	std::swap(storage[n], storage[left_child]);
	return bubbledw(storage, left_child, max_size);
      }
    }
    // right is smaller
    if (weight(storage, right_child) < current) {
      std::swap(storage[n], storage[right_child]);
      return bubbledw(storage, right_child, max_size);
    }
    
  }

  template<typename external_space_point_t>	
  void CandidatesForSpM<external_space_point_t>::bubbleup(std::vector< value_type >& storage,
       							  std::size_t n)
  {
    if (n == 0) return;
    
    // parent: (n - 1) / 2;
    std::size_t parent_idx = (n - 1) / 2;

    std::size_t weight_current = weight(storage, n);
    std::size_t weight_parent = weight(storage, parent_idx);
    if (weight_parent <= weight_current)
      { return; }
    
    std::swap(storage[n], storage[parent_idx]);
    bubbleup(storage, parent_idx);
  }

  template<typename external_space_point_t>
  void CandidatesForSpM<external_space_point_t>::pop(std::vector< value_type >& storage,
       						     std::size_t& current_size,
						     const std::size_t& current_max_size)
  {
    storage[0] = storage[--current_size];
    bubbledw(storage, 0, current_max_size);
  }
  
} //namespace

