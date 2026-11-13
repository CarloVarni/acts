#pragma once

#include "traccc/device/global_index.hpp"
#include <limits>


namespace traccc::device {

TRACCC_HOST_DEVICE
inline scalar rounded_cot_theta(const scalar zA, const scalar rA,
                                const scalar zB, const scalar rB) {
    static constexpr scalar granularity = 20.f;
    const scalar cot_theta = (zB - zA) / (rB - rA);
    return std::round(cot_theta * granularity) / granularity;
}

TRACCC_HOST_DEVICE
void compute_cottheta_matrix(const global_index_t globalIndex,
                const global_index_t n_col,
                const global_index_t n_row,
                const unsigned int n_images,
                const edm::spacepoint_collection::const_view& spacepoints_view,
                const vecmem::data::vector_view<const scalar>& spacepoint_r_view,
                const vecmem::data::vector_view<const unsigned int>& nn_outcome_view,
                const vecmem::data::jagged_vector_view<const global_index_t>& hits_view,
                vecmem::data::jagged_vector_view<scalar>& cot_theta_matrix_view) {
    if (globalIndex >= n_images) return;
    
    // does not pass the NN selection
    const vecmem::device_vector<const unsigned int> nn_outcome(nn_outcome_view);
    if (nn_outcome[globalIndex] == 0) return;
 
    const vecmem::jagged_device_vector<const global_index_t> hits(hits_view);
    const auto& associated_hits = hits[globalIndex];
    unsigned int n_hits = associated_hits.capacity();
    if (n_hits < 3) return;
    if (n_col >= n_hits) return;
    if (n_row >= n_hits) return;

    static constexpr scalar default_value = std::numeric_limits<scalar>::max();

    const edm::spacepoint_collection::const_device spacepoints(spacepoints_view);
    const vecmem::device_vector<const scalar> spacepoint_r(spacepoint_r_view); 
    vecmem::jagged_device_vector<scalar> cot_theta_matrix(cot_theta_matrix_view);

    // compute cot-theta matrix
    unsigned int matrix_index = n_col * n_hits + n_row;
    if (n_col == n_row) {
        cot_theta_matrix[globalIndex][matrix_index] = default_value;
        return;    
    }

    const unsigned int i_index = associated_hits[n_col];
    const scalar ri = spacepoint_r[i_index];
    const scalar zi = spacepoints[i_index].z();

    const unsigned int j_index = associated_hits[n_row];
    const scalar rj = spacepoint_r[j_index];
    const scalar zj = spacepoints[j_index].z();
    
    const scalar cot = rounded_cot_theta(zi, ri, zj, rj);
    if (cot != cot) {
        cot_theta_matrix[globalIndex][matrix_index] = default_value;
        return;
    }
  
    if (cot == scalar{0} && (zi > 550.f || zi < -550.f)) {
        cot_theta_matrix[globalIndex][matrix_index] = default_value;
        return;
    }
 
    cot_theta_matrix[globalIndex][matrix_index] = cot;
}

} // namespace
