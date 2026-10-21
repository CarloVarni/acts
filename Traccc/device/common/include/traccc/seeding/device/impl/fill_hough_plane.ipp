#pragma once

#include "traccc/device/global_index.hpp"
#include "detray/geometry/tracking_surface.hpp"
#include "Acts/Definitions/Units.hpp"
#include "vecmem/memory/device_atomic_ref.hpp"

namespace traccc::device {

TRACCC_HOST_DEVICE
int quant(scalar min, scalar max, unsigned nSteps, scalar val) {
  return static_cast<int>((val - min) / (max - min) * nSteps);
}

TRACCC_HOST_DEVICE
scalar transform(const scalar qopt_value,
                 const scalar r,
                 const scalar phi) {
    // Need to retrieve it from somewhere
    scalar magnetic_field_z = 2.f * Acts::UnitConstants::T;
	    
    const scalar d0 = 0;
    scalar phi_value = asin(0.50f * r * magnetic_field_z * qopt_value - d0 / r) + phi;
    ::remainderf(phi_value, 2.00f * std::numbers::pi_v<scalar>);
    return phi_value;
}

TRACCC_HOST_DEVICE
void ranges(const edm::houghplane_collection<device::N_COLS, device::N_ROWS>::device& container,
            const scalar r,
            const scalar phi,
            unsigned int* min_ranges,
            unsigned int* max_ranges) {
    for (unsigned int i(0); i<device::N_COLS; ++i) {
        scalar min_value_phi = transform(container.edges_x()[i], r, phi);
        scalar max_value_phi = transform(container.edges_x()[i + 1], r, phi);
        if (min_value_phi > max_value_phi) std::swap(min_value_phi, max_value_phi);
        
        // out of bounds
        if (min_value_phi < container.edges_y()[0] ||
            max_value_phi > container.edges_y()[device::N_ROWS]) {
            min_ranges[i] = 0;
            max_ranges[i] = 0;
            return;
        }
        
        unsigned int phi_bin_min = quant(container.edges_y()[0], 
                                         container.edges_y()[device::N_ROWS],
                                         device::N_ROWS, 
                                         min_value_phi);
        unsigned int phi_bin_max = quant(container.edges_y()[0], 
                                         container.edges_y()[device::N_ROWS],
                                         device::N_ROWS, 
                                         max_value_phi) + 1;
        
        // Clamp bins
        if (phi_bin_max > device::N_ROWS) phi_bin_max = device::N_ROWS;
        min_ranges[i] = phi_bin_min;
        max_ranges[i] = phi_bin_max;
    }
}

TRACCC_HOST_DEVICE
void fill_hough_plane(const global_index_t globalIndex,
                edm::spacepoint_collection::const_view::size_type n_spacepoints,
                const edm::spacepoint_collection::const_view& spacepoints_view,
                const vecmem::data::vector_view<unsigned long long>& spacepoint_wedge_association_view,
                const vecmem::data::vector_view<unsigned long long>& spacepoint_detector_id_view,
		const vecmem::data::vector_view<scalar>& spacepoint_r_view,
		const vecmem::data::vector_view<scalar>& spacepoint_phi_view,
                edm::houghplane_collection<device::N_COLS, device::N_ROWS>::view houghplanes_view) {
    if (globalIndex > device::N_COLS * device::N_ROWS) return;
    
    edm::houghplane_collection<device::N_COLS, device::N_ROWS>::device houghplanes(houghplanes_view);
    const edm::spacepoint_collection::const_device spacepoints(spacepoints_view);

    const vecmem::device_vector<unsigned long long> spacepoint_wedge_association(spacepoint_wedge_association_view);
    const vecmem::device_vector<unsigned long long> spacepoint_detector_id(spacepoint_detector_id_view);
    const vecmem::device_vector<scalar> spacepoint_r(spacepoint_r_view);
    const vecmem::device_vector<scalar> spacepoint_phi(spacepoint_phi_view);

    const unsigned int x_bin = globalIndex / device::N_COLS;
    const unsigned int y_bin  = globalIndex % device::N_ROWS;

    // given the bin, take the edges in both dimensions
    const scalar min_x = houghplanes.edges_x()[x_bin];
    const scalar max_x = houghplanes.edges_x()[x_bin + 1];
    const scalar min_y = houghplanes.edges_y()[y_bin];
    const scalar max_y = houghplanes.edges_y()[y_bin + 1];

    std::array<unsigned long long, 25> binContents {};
    for (unsigned int i(0); i<n_spacepoints; ++i) {
      if (spacepoint_wedge_association[i] == 0ULL) continue;

      scalar min_value_y = transform(min_x, spacepoint_r[i], spacepoint_phi[i]);
      scalar max_value_y = transform(max_x, spacepoint_r[i], spacepoint_phi[i]);

      if (min_value_y > max_y) continue;
      if (max_value_y < min_y) continue;
      
      for (unsigned int wedge_idx(0); wedge_idx<25; ++wedge_idx) {
          if ((spacepoint_wedge_association[i] & (1ULL << wedge_idx)) == 0ULL) continue;
          binContents[wedge_idx] |= spacepoint_detector_id[i];
      }
    }

    for (unsigned int i(0); i<25; ++i) {
	if (binContents[i] == 0ULL) continue;
	houghplanes[i].plane()[globalIndex] = binContents[i];
    }
}

} // namespace
