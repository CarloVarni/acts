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
void select_hits(const global_index_t globalIndex,
                 const unsigned int n_images,
                 const vecmem::data::vector_view<const unsigned int>& nn_outcome_view,
                 const vecmem::data::jagged_vector_view<const scalar>& cot_theta_matrix_view) {
    if (globalIndex >= n_images) return;
    
    // does not pass the NN selection
    const vecmem::device_vector<const unsigned int> nn_outcome(nn_outcome_view);
    if (nn_outcome[globalIndex] == 0) return;
 
}

} // namespace
