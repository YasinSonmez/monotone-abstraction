
// Config structure for scalars
typedef struct {
    int basis_list_size;
    int total_states;
    int state_dim;
} SafetyConfig;

// OpenCL kernel for checking basis safety directly over the basis list
__kernel void check_basis_safety(
    __global char* transitions_bag,    // XU bag: next_state_table
    __global char* lists_bag,          // RO bag: basis_flat_idx followed by basis_list
    __global char* output_bag,         // RW bag: unsafe_flags
    __constant char* config_bag        // Config bag: SafetyConfig struct
) {
    // 1. Parse Config
    __constant SafetyConfig* config = (__constant SafetyConfig*)config_bag;
    int basis_list_size = config->basis_list_size;
    int total_states = config->total_states;
    int state_dim = config->state_dim;

    // 2. Map Bags to Pointers
    __global int* next_state_table = (__global int*)transitions_bag;
    
    // lists_bag contains basis_flat_idx first, then basis_list
    __global int* basis_flat_idx = (__global int*)lists_bag;
    // basis_list starts after basis_list_size integers
    __global int* basis_list = (__global int*)(lists_bag + (sizeof(int) * basis_list_size));
    
    __global int* unsafe_flags = (__global int*)output_bag;

    // 3. Original Logic
    const int basis_idx = get_global_id(0);
    if (basis_idx >= basis_list_size) return;

    const int flat_idx = basis_flat_idx[basis_idx];
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
