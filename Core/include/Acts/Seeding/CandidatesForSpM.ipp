// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

namespace Acts {

template<typename comparator_t>
CandidatesForSpM<comparator_t>::CandidatesForSpM(const comparator_t& cmp, std::size_t n)
  : m_max_size(n),
    m_storage(cmp)
{}

template<typename comparator_t>
void CandidatesForSpM<comparator_t>::push(std::vector< std::tuple< std::size_t, float, float, bool > >& registry,
					  std::size_t SpT, float weight, float zOrigin, bool isQuality)
{
  std::cout << "pushing to " << m_storage.size() <<"/"<<m_max_size<<std::endl;
  // if there is still space, add anything
  if (m_storage.size() < m_max_size) {
    addToCollection(registry, SpT, weight, zOrigin, isQuality);
    return;
  }
  
  // if no space, replace one if quality is enough
  // compare to element with lower weight
  std::size_t index_lower_weight = m_storage.top();
  auto& lower_element = registry[index_lower_weight];
  if (weight <= std::get<Components::WEIGHT>(lower_element))
    return;
  
  // remove element with lower weight and add this one
  addToCollectionWithIndex(registry, SpT, weight, zOrigin, isQuality,
			   index_lower_weight);
}

template<typename comparator_t>
void CandidatesForSpM<comparator_t>::addToCollection(std::vector< std::tuple< std::size_t, float, float, bool > >& registry,
						     std::size_t SpT, float weight, float zOrigin, bool isQuality)
{
  auto toAdd = std::make_tuple(SpT, weight, zOrigin, isQuality);
  registry.push_back( toAdd );
  m_storage.push( registry.size() - 1 );
}  

template<typename comparator_t>
void CandidatesForSpM<comparator_t>::addToCollectionWithIndex(std::vector< std::tuple< std::size_t, float, float, bool > >& registry,
							      std::size_t SpT, float weight, float zOrigin, bool isQuality,
							      std::size_t index)
{
  auto toAdd = std::make_tuple(SpT, weight, zOrigin, isQuality);
  registry[index] = toAdd;
  m_storage.pop();
  m_storage.push( index );
}
    
} //namespace
