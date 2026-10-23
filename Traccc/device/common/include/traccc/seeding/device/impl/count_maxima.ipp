#pragma once

#include "traccc/device/global_index.hpp"
#include "vecmem/memory/device_atomic_ref.hpp"


namespace traccc::device {

	
TRACCC_HOST_DEVICE
void count_maxima(const global_index_t globalIndex,
		const edm::houghplane_collection<unsigned short, device::N_COLS, device::N_ROWS>::view houghplanes_view,
                vecmem::data::vector_view<unsigned int>& counter_view) {
	if (globalIndex >= device::N_COLS) return;
	const edm::houghplane_collection<unsigned short, device::N_COLS, device::N_ROWS>::device houghplanes(houghplanes_view);
    vecmem::device_vector<unsigned int> counter(counter_view);

    const global_index_t column_base = globalIndex * device::N_ROWS;
    for (unsigned int wedge_idx(0); wedge_idx<counter_view.size(); ++wedge_idx) {
        const auto& plane = houghplanes.plane()[wedge_idx];

        unsigned int n_maxima = 0u;
        for (unsigned int n_row = 0; n_row < device::N_ROWS; ++n_row) {
            global_index_t current_index = column_base + n_row;
            unsigned int current_value = plane[current_index];
            if (current_value == 0) continue;
            
            const unsigned int previous_row = (n_row == 0) 
                ? device::N_ROWS - 1 
                : n_row - 1;
            const unsigned int next_row = (n_row == device::N_ROWS - 1) 
                ? 0 
                : n_row + 1;

            // Neighbours within the current column.
            if (current_value < plane[column_base + previous_row] ||
                current_value < plane[column_base + next_row]) {
                continue;
            }

            // Three neighbours in the previous column.
            if (globalIndex > 0) {
                const global_index_t left_base = column_base - device::N_ROWS;

                if (current_value < plane[left_base + previous_row] ||
                    current_value < plane[left_base + n_row] ||
                    current_value < plane[left_base + next_row]) continue;
            }

            // Three neighbours in the next column.
            if (globalIndex + 1 < device::N_COLS) {
                const global_index_t right_base = column_base + device::N_ROWS;

                if (current_value < plane[right_base + previous_row] ||
                    current_value < plane[right_base + n_row] ||
                    current_value < plane[right_base + next_row]) continue;
            }

            ++n_maxima;
        } // loop on rows

        vecmem::device_atomic_ref<unsigned int> maxima_counter(counter[wedge_idx]);
        maxima_counter.fetch_add(n_maxima);
    } // loop on wedges
}

} // namespace
