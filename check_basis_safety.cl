// OpenCL kernel for checking basis safety with bitmap early-exit optimization
#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable

inline bool getBit(global const uchar* bitmap, int idx) {
    return (bitmap[idx / 8] & (1 << (idx % 8))) != 0;
}

inline void setBit(global uchar* bitmap, int idx) {
    int byte_idx = idx / 8;
    int bit_idx = idx % 8;
    // Use atomic_or on int-aligned address for portability (works on macOS and Linux)
    int word_idx = byte_idx / 4;
    int shift = (byte_idx % 4) * 8 + bit_idx;
    atom_or((global int*)(bitmap + word_idx * 4), (1 << shift));
}

kernel void check_basis_safety(
    global const uchar* basis_bitmap,        // Input: which states are in basis (1 bit per state)
    global const int* next_state_table,      // Input: precomputed transitions (total_states * 3)
    global const int* basis_list,            // Input: basis element coordinates (basis_size * state_dim)
    int basis_list_size,                     // Input: number of basis elements
    global uchar* unsafe_bitmap,             // Output: which basis states are unsafe
    int total_states,
    int state_dim
) {
    int state_idx = get_global_id(0);
    
    if (state_idx >= total_states) return;
    
    // Early exit: skip if not in basis
    if (!getBit(basis_bitmap, state_idx)) return;
    
    // Get next state from transition table (stored with stride 3)
    int next_state[3];  // MAX_STATE_DIM = 3
    for (int i = 0; i < state_dim; ++i) {
        next_state[i] = next_state_table[state_idx * 3 + i];
    }
    
    // Check if next_state is invalid
    if (next_state[0] == -1) {
        setBit(unsafe_bitmap, state_idx);
        return;
    }
    
    // Check if next_state is dominated by any basis element
    bool safe = false;
    for (int i = 0; i < basis_list_size; ++i) {
        bool dominated = true;
        for (int j = 0; j < state_dim; ++j) {
            if (next_state[j] > basis_list[i * state_dim + j]) {
                dominated = false;
                break;
            }
        }
        if (dominated) {
            safe = true;
            break;
        }
    }
    
    if (!safe) {
        setBit(unsafe_bitmap, state_idx);
    }
}
