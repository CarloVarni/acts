#pragma once

#include "traccc/device/global_index.hpp"
#include <limits>


namespace traccc::device {

template <typename T>
TRACCC_HOST_DEVICE 
inline void sort_row(vecmem::device_vector<T> row) {
    const auto n = row.capacity();
    for (unsigned int i = 1; i < n; ++i) {
        const T value = row[i];
        unsigned int j = i;
        while (j > 0 && value < row[j - 1]) {
            row[j] = row[j - 1];
            --j;
        }
        row[j] = value;
    }
}

TRACCC_HOST_DEVICE
void make_seeds(const global_index_t globalIndex,
                const unsigned int n_images,
                const edm::spacepoint_collection::const_view& spacepoints_view,
                const vecmem::data::vector_view<const scalar>& spacepoint_r_view,
                const vecmem::data::vector_view<const unsigned int>& nn_outcome_view,
                const vecmem::data::jagged_vector_view<const global_index_t>& hits_view,
                vecmem::data::jagged_vector_view<scalar>& cot_theta_matrix_view,
                edm::seed_collection::view& seeds_view) {
    if (globalIndex >= n_images) return;
    
    // does not pass the NN selection
    const vecmem::device_vector<const unsigned int> nn_outcome(nn_outcome_view);
    if (nn_outcome[globalIndex] == 0) return;
 
    const vecmem::jagged_device_vector<const global_index_t> hits(hits_view);
    if (hits[globalIndex].capacity() < 3) return;

    const auto& associated_hits = hits[globalIndex];

    // create seeds
    edm::seed_collection::device seeds(seeds_view);

    assert(associated_hits.capacity() >= 3);
    seeds[globalIndex].bottom_index() = associated_hits[0];
    seeds[globalIndex].middle_index() = associated_hits[(associated_hits.capacity() - 1) / 2];
    seeds[globalIndex].top_index() = associated_hits[associated_hits.capacity() - 1];
    seeds[globalIndex].quality() = 0.00f;
}

} // namespace
