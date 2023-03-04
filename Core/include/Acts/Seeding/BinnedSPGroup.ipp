// -*- C++ -*-
// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Binned SP Group Iterator

template <typename external_spacepoint_t>
Acts::BinnedSPGroupIterator<external_spacepoint_t>::BinnedSPGroupIterator(Acts::BinnedSPGroup<external_spacepoint_t>& group,
									  std::size_t index)
  : m_group(group),
    m_max_localBins( m_group->m_grid->numLocalBins() )
{
  if (index == m_group->m_grid->size()) {
    m_current_localBins = m_max_localBins;
  }
  
  // Go to the next not-empty bin
  findNotEmptyBin();
}


template <typename external_spacepoint_t>
inline 
Acts::BinnedSPGroupIterator<external_spacepoint_t>&
Acts::BinnedSPGroupIterator<external_spacepoint_t>::operator++()
{
  // Increase the position by one
  // if we were on the edge, go up one phi bin and reset z bin
  if (++m_current_localBins[INDEX::Z] == m_max_localBins[INDEX::Z]) {
    ++m_current_localBins[INDEX::PHI];
    m_current_localBins[INDEX::Z] = 0;
  }
  
  // Get the next not-empty bin in the grid
  findNotEmptyBin();
  return *this;
}

template <typename external_spacepoint_t>
inline	
bool
Acts::BinnedSPGroupIterator<external_spacepoint_t>::operator==(const Acts::BinnedSPGroupIterator<external_spacepoint_t>& other) const {
  return m_group.ptr == other.m_group.ptr and
    m_current_localBins[INDEX::PHI] == m_max_localBins[INDEX::PHI] and
    m_current_localBins[INDEX::Z] == m_max_localBins[INDEX::Z];
}

template <typename external_spacepoint_t>
inline
bool
Acts::BinnedSPGroupIterator<external_spacepoint_t>::operator!=(const Acts::BinnedSPGroupIterator<external_spacepoint_t>& other) const
{ return not (*this == other); }

template <typename external_spacepoint_t>
std::tuple< std::vector<Acts::InternalSpacePoint<external_spacepoint_t>*>,
	    std::vector<Acts::InternalSpacePoint<external_spacepoint_t>*>,
	    std::vector<Acts::InternalSpacePoint<external_spacepoint_t>*> >
Acts::BinnedSPGroupIterator<external_spacepoint_t>::operator*() {
  // Retrieve here - this is the heavy lifting
  // Less expensive then doing it in the operator++
  std::vector<Acts::InternalSpacePoint<external_spacepoint_t>*> bottoms;
  std::vector<Acts::InternalSpacePoint<external_spacepoint_t>*> middles;
  std::vector<Acts::InternalSpacePoint<external_spacepoint_t>*> tops; 

  std::size_t nBottoms = 0;
  std::size_t nTops = 0;

  // Global Index
  std::size_t global_index = m_group->m_grid->globalBinFromLocalBins({m_current_localBins[INDEX::PHI], m_group->m_bins[m_current_localBins[INDEX::Z]]});

  // Middles
  auto& collection_middles = m_group->m_grid->at(global_index);
  middles.reserve(collection_middles.size());
  for (auto& sp : collection_middles)
    middles.push_back( sp.get() );
  
  // Bottoms
  auto bottomBinIndices = m_group->m_bottomBinFinder->findBins(m_current_localBins[INDEX::PHI],
							       m_group->m_bins[m_current_localBins[INDEX::Z]],
							       m_group->m_grid.get());

  // Get n bottoms
  for (auto idx : bottomBinIndices) {
    nBottoms += m_group->m_grid->at(idx).size();
  }
  bottoms.reserve(nBottoms);
  
  for (auto idx : bottomBinIndices) {
    auto& collection_bottoms = m_group->m_grid->at(idx);
    for (auto& el : collection_bottoms) {
      bottoms.push_back( el.get() );
    }
  }
  
  // Tops
  auto topBinIndices = m_group->m_topBinFinder->findBins(m_current_localBins[INDEX::PHI],
							 m_group->m_bins[m_current_localBins[INDEX::Z]],
							 m_group->m_grid.get());
  
  for (auto idx : topBinIndices) {
    nTops += m_group->m_grid->at(idx).size();
  }
  
  for (auto idx : topBinIndices) {
    auto& collection_tops = m_group->m_grid->at(idx);
    for (auto& el : collection_tops) {
      tops.push_back( el.get() );
    }
  }

  return std::make_tuple( std::move(bottoms),
			  std::move(middles),
			  std::move(tops) );
}

template <typename external_spacepoint_t>
inline
void
Acts::BinnedSPGroupIterator<external_spacepoint_t>::findNotEmptyBin()
{
  // Iterate on the grid till we find a not-empty bin
  // We start from the current bin configuration and move forward
  for (std::size_t phiBin(m_current_localBins[INDEX::PHI]);
       phiBin < m_max_localBins[INDEX::PHI];
       ++phiBin) {
    for (std::size_t zBin(m_current_localBins[INDEX::Z]);
	 zBin < m_max_localBins[INDEX::Z];
	 ++zBin) {
      
      std::size_t zBinIndex = m_group->m_bins.size() == 0
	? zBin
	: m_group->m_bins[zBin];
      std::size_t index = m_group->m_grid->globalBinFromLocalBins({phiBin, zBinIndex});
      
      // Check if there are entries in this bin
      if (m_group->m_grid->at(index).size() == 0) {
	continue;
      }
      
      // Set the new current bins
      m_current_localBins[INDEX::PHI] = phiBin;
      m_current_localBins[INDEX::Z] = zBin;
      return;
    }
    // Reset z-index
    m_current_localBins[INDEX::Z] = 0;
  }
  
  // Could find nothing ... setting this to end()
  m_current_localBins = m_max_localBins;
}




// Binned SP Group
template <typename external_spacepoint_t>
template <typename spacepoint_iterator_t, typename callable_t>
Acts::BinnedSPGroup<external_spacepoint_t>::BinnedSPGroup(
    spacepoint_iterator_t spBegin, spacepoint_iterator_t spEnd,
    callable_t&& toGlobal,
    std::shared_ptr<const Acts::BinFinder<external_spacepoint_t>> botBinFinder,
    std::shared_ptr<const Acts::BinFinder<external_spacepoint_t>> tBinFinder,
    std::unique_ptr<SpacePointGrid<external_spacepoint_t>> grid,
    Acts::Extent& rRangeSPExtent,
    const SeedFinderConfig<external_spacepoint_t>& config,
    const SeedFinderOptions& options) {
  if (not config.isInInternalUnits) {
    throw std::runtime_error(
        "SeedFinderConfig not in ACTS internal units in BinnedSPGroup");
  }
  if (not options.isInInternalUnits) {
    throw std::runtime_error(
        "SeedFinderOptions not in ACTS internal units in BinnedSPGroup");
  }
  static_assert(
      std::is_same<
          typename std::iterator_traits<spacepoint_iterator_t>::value_type,
          const external_spacepoint_t*>::value,
      "Iterator does not contain type this class was templated with");

  // get region of interest (or full detector if configured accordingly)
  float phiMin = config.phiMin;
  float phiMax = config.phiMax;
  float zMin = config.zMin;
  float zMax = config.zMax;

  // sort by radius
  // add magnitude of beamPos to rMax to avoid excluding measurements
  // create number of bins equal to number of millimeters rMax
  // (worst case minR: configured minR + 1mm)
  // binSizeR allows to increase or reduce numRBins if needed
  size_t numRBins = static_cast<size_t>((config.rMax + options.beamPos.norm()) /
                                        config.binSizeR);
  std::vector<
      std::vector<std::unique_ptr<InternalSpacePoint<external_spacepoint_t>>>>
      rBins(numRBins);
  for (spacepoint_iterator_t it = spBegin; it != spEnd; it++) {
    if (*it == nullptr) {
      continue;
    }
    const external_spacepoint_t& sp = **it;
    const auto& [spPosition, variance] =
        toGlobal(sp, config.zAlign, config.rAlign, config.sigmaError);

    float spX = spPosition[0];
    float spY = spPosition[1];
    float spZ = spPosition[2];

    // store x,y,z values in extent
    rRangeSPExtent.extend({spX, spY, spZ});

    if (spZ > zMax || spZ < zMin) {
      continue;
    }
    float spPhi = std::atan2(spY, spX);
    if (spPhi > phiMax || spPhi < phiMin) {
      continue;
    }

    auto isp = std::make_unique<InternalSpacePoint<external_spacepoint_t>>(
        sp, spPosition, options.beamPos, variance);
    // calculate r-Bin index and protect against overflow (underflow not
    // possible)
    size_t rIndex = static_cast<size_t>(isp->radius() / config.binSizeR);
    // if index out of bounds, the SP is outside the region of interest
    if (rIndex >= numRBins) {
      continue;
    }
    rBins[rIndex].push_back(std::move(isp));
  }

  // if requested, it is possible to force sorting in R for each (z, phi) grid
  // bin
  if (config.forceRadialSorting) {
    for (auto& rbin : rBins) {
      std::sort(
          rbin.begin(), rbin.end(),
          [](std::unique_ptr<InternalSpacePoint<external_spacepoint_t>>& a,
             std::unique_ptr<InternalSpacePoint<external_spacepoint_t>>& b) {
            return a->radius() < b->radius();
          });
    }
  }

  // fill rbins into grid such that each grid bin is sorted in r
  // space points with delta r < rbin size can be out of order is sorting is not
  // requested
  for (auto& rbin : rBins) {
    for (auto& isp : rbin) {
      Acts::Vector2 spLocation(isp->phi(), isp->z());
      std::vector<std::unique_ptr<InternalSpacePoint<external_spacepoint_t>>>&
          bin = grid->atPosition(spLocation);
      bin.push_back(std::move(isp));
    }
  }
  m_grid = std::move(grid);
  m_bottomBinFinder = botBinFinder;
  m_topBinFinder = tBinFinder;

  m_bins = config.zBinsCustomLooping;
}



template <typename external_spacepoint_t>
inline
size_t
Acts::BinnedSPGroup<external_spacepoint_t>::size() const
{ return m_grid->size(); }

template <typename external_spacepoint_t>
inline
Acts::BinnedSPGroupIterator<external_spacepoint_t>
Acts::BinnedSPGroup<external_spacepoint_t>::begin()
{ return {*this, 0}; }

template <typename external_spacepoint_t>
inline
Acts::BinnedSPGroupIterator<external_spacepoint_t>
Acts::BinnedSPGroup<external_spacepoint_t>::end()
{ return {*this, m_grid->size()}; }

