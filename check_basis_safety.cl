// OpenCL kernel for checking basis safety directly over the basis list

kernel void check_basis_safety(
    global const int* basis_flat_idx,        // Flattened indices of basis elements
    global const int* next_state_table,      // Precomputed transitions (total_states * 3)
    global const int* basis_list,            // Basis element coordinates (basis_size * state_dim)
    int basis_list_size,                     // Number of basis elements
    global int* unsafe_flags,                // Output: 1 if unsafe, 0 otherwise
    int total_states,
    int state_dim
) {
    int basis_idx = get_global_id(0);
    if (basis_idx >= basis_list_size) return;

    int flat_idx = basis_flat_idx[basis_idx];
    if (flat_idx < 0 || flat_idx >= total_states) {
        unsafe_flags[basis_idx] = 1;
        return;
    }

    int next_state[3] = {0, 0, 0};  // MAX_STATE_DIM = 3
    for (int i = 0; i < state_dim; ++i) {
        next_state[i] = next_state_table[flat_idx * 3 + i];
    }

    if (next_state[0] == -1) {
        unsafe_flags[basis_idx] = 1;
        return;
    }

    int unsafe = 1;
    for (int i = 0; i < basis_list_size; ++i) {
        int dominated = 1;
        for (int j = 0; j < state_dim; ++j) {
            if (next_state[j] > basis_list[i * state_dim + j]) {
                dominated = 0;
                break;
            }
        }
        if (dominated) {
            unsafe = 0;
            break;
        }
    }

    unsafe_flags[basis_idx] = unsafe;
}
