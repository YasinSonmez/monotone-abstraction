#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <cstring>
#include <vector>

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

// ============================================================================
// CONFIGURATION
// ============================================================================

constexpr double DT_2D = 0.1, DT_3D = 0.4;
constexpr double H_MIN_2D = 5.0, H_MAX_2D = 80.0, V_MIN_2D = 0.0, V_MAX_2D = 20.0;
constexpr double H_MIN_3D = 0.0, H_MAX_3D = 80.0, V_MIN_3D = 0.0, V_MAX_3D = 20.0;
constexpr double H_RES_2D = 0.2, V_RES_2D = 0.1, H_RES_3D = 0.8, V_RES_3D = 0.4;
///TODO: @MK: instead of reading 2D/3D from command line, what about letting the user define it here as STATE_DIM instead of MAX_STATE_DIM? You could keep both logics in the code thought. I am just referring to the input.
constexpr int MAX_STATE_DIM = 3;
constexpr int MAX_BASIS_ELEMENTS = 2000;

// ============================================================================
// UTILITIES
// ============================================================================

std::string loadKernelSource(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error(std::string("Failed to open: ") + filename);
    }
    std::stringstream buffer;
    buffer << file.rdbuf();
    return buffer.str();
}

void checkError(const cl_int error, const char* message) {
    if (error != CL_SUCCESS) {
        std::cerr << "Error (" << error << "): " << message << std::endl;
        throw std::runtime_error(message);
    }
}

// ============================================================================
// MONOTONE ABSTRACTION
// ============================================================================

class MonotoneAbstraction {
private:
    // State space
    double x_range_min[MAX_STATE_DIM], x_range_max[MAX_STATE_DIM], x_res[MAX_STATE_DIM];
    int x_numCells[MAX_STATE_DIM], x_priority[MAX_STATE_DIM];
    
    // Safe set basis
    int safe_set_basis[MAX_BASIS_ELEMENTS * MAX_STATE_DIM];
    int safe_set_flat_indices[MAX_BASIS_ELEMENTS];
    int safe_set_size;
    
    // Precomputed transitions
    ///TODO: @MK: maybe a std::vector, or at least a smart pointer?
    int* next_state_table;
    int total_states;
    
    // Configuration
    bool is_3d;
    double dt;
    
    // OpenCL
    cl_context context;
    cl_command_queue queue;
    cl_kernel kernel_trans, kernel_safety;
    cl_program program_trans, program_safety;
    
public:
    MonotoneAbstraction(bool use_3d) 
        : safe_set_size(0), is_3d(use_3d), next_state_table(nullptr),
          context(nullptr), queue(nullptr), 
          kernel_trans(nullptr), kernel_safety(nullptr),
          program_trans(nullptr), program_safety(nullptr) {
        memset(safe_set_basis, 0, sizeof(safe_set_basis));
        memset(safe_set_flat_indices, 0, sizeof(safe_set_flat_indices));
        dt = is_3d ? DT_3D : DT_2D;
    }
    
    ~MonotoneAbstraction() {
        if (next_state_table) delete[] next_state_table;
        if (kernel_trans) clReleaseKernel(kernel_trans);
        if (kernel_safety) clReleaseKernel(kernel_safety);
        if (program_trans) clReleaseProgram(program_trans);
        if (program_safety) clReleaseProgram(program_safety);
        if (queue) clReleaseCommandQueue(queue);
        if (context) clReleaseContext(context);
    }
    
    void setupRanges() {
        ///TODO: @MK: same as above comment regarding STATE_DIM.
        int state_dim = is_3d ? 3 : 2;
        
        if (is_3d) {
            x_range_min[0] = H_MIN_3D; x_range_max[0] = H_MAX_3D; x_res[0] = H_RES_3D;
            x_range_min[1] = V_MIN_3D; x_range_max[1] = V_MAX_3D; x_res[1] = V_RES_3D;
            x_range_min[2] = V_MIN_3D; x_range_max[2] = V_MAX_3D; x_res[2] = V_RES_3D;
            x_priority[0] = 0; x_priority[1] = 1; x_priority[2] = 0;
        } else {
            x_range_min[0] = H_MIN_2D; x_range_max[0] = H_MAX_2D; x_res[0] = H_RES_2D;
            x_range_min[1] = V_MIN_2D; x_range_max[1] = V_MAX_2D; x_res[1] = V_RES_2D;
            x_priority[0] = 0; x_priority[1] = 1;
        }
        
        total_states = 1;
        for (int i = 0; i < state_dim; ++i) {
            x_numCells[i] = std::ceil((x_range_max[i] - x_range_min[i]) / x_res[i]) + 1;
            total_states *= x_numCells[i];
        }
    }
    
    void initializeOpenCL() {
        cl_int err;
        cl_platform_id platform;
        cl_device_id device;
        
        // Get platform and device
        err = clGetPlatformIDs(1, &platform, nullptr);
        checkError(err, "Get platform failed");
        
        err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 1, &device, nullptr);
        if (err != CL_SUCCESS) {
            err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_CPU, 1, &device, nullptr);
            checkError(err, "Get device failed");
        }
        
        // Create context and queue
        context = clCreateContext(nullptr, 1, &device, nullptr, nullptr, &err);
        checkError(err, "Create context failed");
        
        queue = clCreateCommandQueue(context, device, 0, &err);
        checkError(err, "Create queue failed");
        
        // Build kernels
        buildKernel("precompute_transitions.cl", "precompute_transitions", program_trans, kernel_trans, device);
        buildKernel("check_basis_safety.cl", "check_basis_safety", program_safety, kernel_safety, device);
    }
    
    void initializeSafeSet() {
        safe_set_size = 0;
        int state_dim = is_3d ? 3 : 2;
        for (int i = 0; i < state_dim; ++i) {
            safe_set_basis[i] = x_numCells[i];
        }
        safe_set_size = 1;
        safe_set_flat_indices[0] = flattenIndex(&safe_set_basis[0]);
    }
    
    void precomputeNextStateTableGPU() {
        int state_dim = is_3d ? 3 : 2;
        next_state_table = new int[total_states * MAX_STATE_DIM];
        
        // Try cache
        const char* cache_file = is_3d ? "transition_cache_3d.bin" : "transition_cache_2d.bin";
        std::ifstream cache_in(cache_file, std::ios::binary);
        if (cache_in.good()) {
            std::cout << "Loading transition table from cache: " << cache_file << std::endl << std::flush;
            int cached_states;
            cache_in.read(reinterpret_cast<char*>(&cached_states), sizeof(int));
            if (cached_states == total_states) {
                cache_in.read(reinterpret_cast<char*>(next_state_table), 
                             sizeof(int) * total_states * MAX_STATE_DIM);
                cache_in.close();
                return;
            }
            cache_in.close();
        }
        
        // Compute on GPU
        std::cout << "Creating buffers..." << std::endl << std::flush;
        cl_int err;
        cl_mem d_total_states = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(int), &total_states, &err);
        checkError(err, "Create total_states buffer failed");
        cl_mem d_state_dim = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(int), &state_dim, &err);
        checkError(err, "Create state_dim buffer failed");
        cl_mem d_is_3d_int = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(bool), &is_3d, &err);
        checkError(err, "Create is_3d buffer failed");
        cl_mem d_dt = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(double), &dt, &err);
        checkError(err, "Create dt buffer failed");
        cl_mem d_x_range_min = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(double) * MAX_STATE_DIM, x_range_min, &err);
        checkError(err, "Create x_range_min buffer failed");
        cl_mem d_x_range_max = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(double) * MAX_STATE_DIM, x_range_max, &err);
        checkError(err, "Create x_range_max buffer failed");
        cl_mem d_x_res = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(double) * MAX_STATE_DIM, x_res, &err);
        checkError(err, "Create x_res buffer failed");
        cl_mem d_x_numCells = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(int) * MAX_STATE_DIM, x_numCells, &err);
        checkError(err, "Create x_numCells buffer failed");
        cl_mem d_x_priority = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(int) * MAX_STATE_DIM, x_priority, &err);
        checkError(err, "Create x_priority buffer failed");
        cl_mem d_next_state_table = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(int) * total_states * MAX_STATE_DIM, nullptr, &err);
        checkError(err, "Create next_state_table buffer failed");
        
        ///TODO: @MK: I'd suggest collecting all read-only vars into one struct and passing a single buffer to reduce the number of kernel arguments (may also reduce memory due to padding).
        std::cout << "Setting kernel arguments..." << std::endl;
        int is_3d_int = is_3d ? 1 : 0;
        err = clSetKernelArg(kernel_trans, 0, sizeof(cl_mem), &d_total_states);
        checkError(err, "Set kernel arg 0 failed");
        err = clSetKernelArg(kernel_trans, 1, sizeof(cl_mem), &d_state_dim);
        checkError(err, "Set kernel arg 1 failed");
        err = clSetKernelArg(kernel_trans, 2, sizeof(cl_mem), &d_is_3d_int);
        checkError(err, "Set kernel arg 2 failed");
        err = clSetKernelArg(kernel_trans, 3, sizeof(cl_mem), &d_dt);
        checkError(err, "Set kernel arg 3 failed");
        err = clSetKernelArg(kernel_trans, 4, sizeof(cl_mem), &d_x_range_min);
        checkError(err, "Set kernel arg 4 failed");
        err = clSetKernelArg(kernel_trans, 5, sizeof(cl_mem), &d_x_range_max);
        checkError(err, "Set kernel arg 5 failed");
        err = clSetKernelArg(kernel_trans, 6, sizeof(cl_mem), &d_x_res);
        checkError(err, "Set kernel arg 6 failed");
        err = clSetKernelArg(kernel_trans, 7, sizeof(cl_mem), &d_x_numCells);
        checkError(err, "Set kernel arg 7 failed");
        err = clSetKernelArg(kernel_trans, 8, sizeof(cl_mem), &d_x_priority);
        checkError(err, "Set kernel arg 8 failed");
        err = clSetKernelArg(kernel_trans, 9, sizeof(cl_mem), &d_next_state_table);
        checkError(err, "Set kernel arg 9 failed");
        
        std::cout << "Launching kernel..." << std::endl;
        size_t global = ((total_states + 255) / 256) * 256;
        clEnqueueNDRangeKernel(queue, kernel_trans, 1, nullptr, &global, nullptr, 0, nullptr, nullptr);
        clFinish(queue);
        clEnqueueReadBuffer(queue, d_next_state_table, CL_TRUE, 0,
            sizeof(int) * total_states * MAX_STATE_DIM, next_state_table, 0, nullptr, nullptr);
        
        // Cleanup
        std::cout << "Releasing buffers..." << std::endl;
        clReleaseMemObject(d_total_states);
        clReleaseMemObject(d_state_dim);
        clReleaseMemObject(d_is_3d_int);
        clReleaseMemObject(d_dt);
        clReleaseMemObject(d_x_range_min);
        clReleaseMemObject(d_x_range_max);
        clReleaseMemObject(d_x_res);
        clReleaseMemObject(d_x_numCells);
        clReleaseMemObject(d_x_priority);
        clReleaseMemObject(d_next_state_table);
        
        // Save cache
        std::ofstream cache_out(cache_file, std::ios::binary);
        if (cache_out.good()) {
            cache_out.write(reinterpret_cast<const char*>(&total_states), sizeof(int));
            cache_out.write(reinterpret_cast<const char*>(next_state_table), 
                           sizeof(int) * total_states * MAX_STATE_DIM);
        }
    }
    
    double computeSafeSetGPU(int& out_iterations, int& out_basis_elements) {
        out_iterations = 0;
        int state_dim = is_3d ? 3 : 2;
        
        // Allocate GPU buffers ONCE (no reallocation overhead)
        // Use pinned memory (CL_MEM_ALLOC_HOST_PTR) for faster transfers
        std::cout << "Allocating GPU buffers..." << std::endl << std::flush;
        ///TODO: @MK: error handling for OpenCL calls missing here and elsewhere.
        cl_int err;
        cl_mem d_total_states = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(int), &total_states, &err);
        checkError(err, "Create total_states buffer failed");
        cl_mem d_state_dim = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(int), &state_dim, &err);
        checkError(err, "Create state_dim buffer failed");
        cl_mem d_next_state_table = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(int) * total_states * MAX_STATE_DIM, next_state_table, &err);
        checkError(err, "Create next_state_table buffer failed");
        cl_mem d_safe_set_size = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(int), &safe_set_size, &err);
        checkError(err, "Create safe_set_size buffer failed");
        cl_mem d_basis_list = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR, sizeof(int) * MAX_BASIS_ELEMENTS * MAX_STATE_DIM, nullptr, &err);
        checkError(err, "Create basis_list buffer failed");
        cl_mem d_basis_flat = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR, sizeof(int) * MAX_BASIS_ELEMENTS, nullptr, &err);
        checkError(err, "Create basis_flat buffer failed");
        cl_mem d_unsafe_flags = clCreateBuffer(context, CL_MEM_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR, sizeof(int) * MAX_BASIS_ELEMENTS, nullptr, &err);
        checkError(err, "Create unsafe_flags buffer failed");
        
        // Allocate CPU buffers
        int unsafe_flags[MAX_BASIS_ELEMENTS];
        int neighbor_buffer[MAX_BASIS_ELEMENTS * MAX_STATE_DIM * MAX_STATE_DIM];
        unsigned char unsafe_mask[MAX_BASIS_ELEMENTS];
        unsigned char* seen_neighbors = new unsigned char[total_states];
        
        // Set static kernel args
        /// TODO: @MK: returned error handling missing here.
        std::cout << "Setting kernel arguments..." << std::endl << std::flush;
        err = clSetKernelArg(kernel_safety, 0, sizeof(cl_mem), &d_basis_flat);
        checkError(err, "Set kernel arg 0 failed");
        err = clSetKernelArg(kernel_safety, 1, sizeof(cl_mem), &d_next_state_table);
        checkError(err, "Set kernel arg 1 failed");
        err = clSetKernelArg(kernel_safety, 2, sizeof(cl_mem), &d_basis_list);
        checkError(err, "Set kernel arg 2 failed");
        err = clSetKernelArg(kernel_safety, 4, sizeof(cl_mem), &d_unsafe_flags);
        checkError(err, "Set kernel arg 4 failed");
        err = clSetKernelArg(kernel_safety, 5, sizeof(cl_mem), &d_total_states);
        checkError(err, "Set kernel arg 5 failed");
        err = clSetKernelArg(kernel_safety, 6, sizeof(cl_mem), &d_state_dim);
        checkError(err, "Set kernel arg 6 failed");
        
        auto compute_start = std::chrono::high_resolution_clock::now();
        
        // Fixed-point iteration
        std::cout << "Starting fixed-point iterations..." << std::endl << std::flush;
        while (true) {
            out_iterations++;
            if (safe_set_size == 0) break;
            
            // Pipelined GPU operations: upload -> kernel -> readback
            size_t basis_bytes = sizeof(int) * safe_set_size * state_dim;
            size_t flat_bytes = sizeof(int) * safe_set_size;
            
            // Async uploads
            ///TODO: @MK: are events really needed here? You could just use clFinish after each step to simplify the code. You are not measuring performance here anyway.
            cl_event upload_events[2];
            err = clEnqueueWriteBuffer(queue, d_basis_list, CL_FALSE, 0, basis_bytes, safe_set_basis, 0, nullptr, &upload_events[0]);
            checkError(err, "Enqueue write buffer d_basis_list failed");
            err = clEnqueueWriteBuffer(queue, d_basis_flat, CL_FALSE, 0, flat_bytes, safe_set_flat_indices, 0, nullptr, &upload_events[1]);
            checkError(err, "Enqueue write buffer d_basis_flat failed");
            
            // Kernel waits for uploads
            err = clSetKernelArg(kernel_safety, 3, sizeof(cl_mem), &d_safe_set_size);
            checkError(err, "Set kernel arg 3 failed");
            size_t global = ((safe_set_size + 127) / 128) * 128;
            cl_event kernel_event;
            err = clEnqueueNDRangeKernel(queue, kernel_safety, 1, nullptr, &global, nullptr, 2, upload_events, &kernel_event);
            checkError(err, "Enqueue NDRange kernel failed");
            
            // Readback waits for kernel
            cl_event read_event;
            err = clEnqueueReadBuffer(queue, d_unsafe_flags, CL_FALSE, 0, flat_bytes, unsafe_flags, 1, &kernel_event, &read_event);
            checkError(err, "Enqueue read buffer failed");
            
            // Wait for pipeline
            clWaitForEvents(1, &read_event);
            clReleaseEvent(upload_events[0]);
            clReleaseEvent(upload_events[1]);
            clReleaseEvent(kernel_event);
            clReleaseEvent(read_event);
            
            // Update safe set
            int added = updateSafeSet(unsafe_flags, unsafe_mask, seen_neighbors, neighbor_buffer, state_dim);
            if (added == 0) break;
        }
        
        auto compute_end = std::chrono::high_resolution_clock::now();
        out_basis_elements = safe_set_size;
        
        delete[] seen_neighbors;
        clReleaseMemObject(d_total_states);
        clReleaseMemObject(d_state_dim);
        clReleaseMemObject(d_next_state_table);
        clReleaseMemObject(d_basis_list);
        clReleaseMemObject(d_basis_flat);
        clReleaseMemObject(d_unsafe_flags);
        
        return std::chrono::duration<double, std::milli>(compute_end - compute_start).count();
    }
    
private:
    void buildKernel(const std::string& file, const std::string& name, cl_program& prog, cl_kernel& kern, cl_device_id dev) {
        cl_int err;
        std::string src = loadKernelSource(file);
        const char* src_str = src.c_str();
        size_t src_size = src.length();
        
        prog = clCreateProgramWithSource(context, 1, &src_str, &src_size, &err);
        checkError(err, "Create program failed");
        
        err = clBuildProgram(prog, 1, &dev, nullptr, nullptr, nullptr);
        if (err != CL_SUCCESS) {
            size_t log_size;
            clGetProgramBuildInfo(prog, dev, CL_PROGRAM_BUILD_LOG, 0, nullptr, &log_size);
            char* log = new char[log_size];
            clGetProgramBuildInfo(prog, dev, CL_PROGRAM_BUILD_LOG, log_size, log, nullptr);
            std::cerr << log << std::endl;
            delete[] log;
            throw std::runtime_error("Build failed");
        }
        
        kern = clCreateKernel(prog, name.c_str(), &err);
        checkError(err, "Create kernel failed");
    }
    
    ///TODO: @MK: We should discuss how to parallelize this as well, or parts of it at least.
    int updateSafeSet(int* unsafe_flags, unsigned char* unsafe_mask, 
                      unsigned char* seen_neighbors, int* neighbor_buffer, int state_dim) {
        // Clear masks
        std::fill_n(unsafe_mask, safe_set_size, static_cast<unsigned char>(0));
        memset(seen_neighbors, 0, total_states);
        
        // Generate unique neighbors from unsafe elements
        int neighbor_count = 0;
        for (int i = 0; i < safe_set_size; ++i) {
            if (!unsafe_flags[i]) continue;
            unsafe_mask[i] = 1;
            
            // Generate neighbors by decrementing each dimension
            for (int j = 0; j < state_dim; ++j) {
                int val = safe_set_basis[i * MAX_STATE_DIM + j];
                if (val == 1) continue;
                
                int coords[MAX_STATE_DIM];
                for (int k = 0; k < state_dim; ++k) {
                    coords[k] = safe_set_basis[i * MAX_STATE_DIM + k] - (j == k ? 1 : 0);
                }
                
                // Deduplicate using flat index
                int flat_idx = flattenIndex(coords);
                if (!seen_neighbors[flat_idx]) {
                    seen_neighbors[flat_idx] = 1;
                    int* neighbor = &neighbor_buffer[neighbor_count * MAX_STATE_DIM];
                    for (int k = 0; k < state_dim; ++k) {
                        neighbor[k] = coords[k];
                    }
                    neighbor_count++;
                }
            }
        }

        ///TODO: @MK: [TO_DISCUSS] The code from the beginning down to here (clearing mem + the first loop) can be done in a separate kernel the GPU.
        
        ///TODO: @MK: [TO_DISCUSS] This next loop seems like a Run-length Encoding (RLE)-like pattern. Can be also be parallelized on GPU.
        // Remove unsafe elements
        int write_pos = 0;
        for (int i = 0; i < safe_set_size; ++i) {
            if (!unsafe_mask[i]) {
                if (write_pos != i) {
                    for (int j = 0; j < state_dim; ++j) {
                        safe_set_basis[write_pos * MAX_STATE_DIM + j] = safe_set_basis[i * MAX_STATE_DIM + j];
                    }
                    safe_set_flat_indices[write_pos] = safe_set_flat_indices[i];
                }
                write_pos++;
            }
        }
        safe_set_size = write_pos;
        

        ///TODO: @MK: [TO_DISCUSS] This is to be kept serial due to the possible non-deterministic order of adding new neighbors. Or we need to do more research on the theory side regarding order of additions.

        // Add new neighbors
        int added = 0;
        for (int ni = 0; ni < neighbor_count; ++ni) {
            int* neighbor = &neighbor_buffer[ni * MAX_STATE_DIM];
            if (!xInSafeSet(neighbor)) {
                if (safe_set_size >= MAX_BASIS_ELEMENTS) break;
                for (int j = 0; j < state_dim; ++j) {
                    safe_set_basis[safe_set_size * MAX_STATE_DIM + j] = neighbor[j];
                }
                safe_set_flat_indices[safe_set_size] = flattenIndex(neighbor);
                safe_set_size++;
                added++;
            }
        }
        return added;
    }
    
    bool xInSafeSet(const int* x_idx) const {
        if (x_idx[0] == -1) return false;
        int state_dim = is_3d ? 3 : 2;
        for (int i = safe_set_size - 1; i >= 0; --i) {
            bool dominated = true;
            for (int j = 0; j < state_dim; ++j) {
                if (x_idx[j] > safe_set_basis[i * MAX_STATE_DIM + j]) {
                    dominated = false;
                    break;
                }
            }
            if (dominated) return true;
        }
        return false;
    }
    
    int flattenIndex(const int* idx) const {
        int result = 0, multiplier = 1;
        int state_dim = is_3d ? 3 : 2;
        for (int i = 0; i < state_dim; ++i) {
            result += (idx[i] - 1) * multiplier;
            multiplier *= x_numCells[i];
        }
        return result;
    }
};

// ============================================================================
// MAIN
// ============================================================================

int main(int argc, char* argv[]) {
    if (argc != 2 || (std::string(argv[1]) != "2d" && std::string(argv[1]) != "3d")) {
        std::cout << "Usage: " << argv[0] << " [2d|3d]" << std::endl;
        return 1;
    }
    
    bool is_3d = (std::string(argv[1]) == "3d");
    const int NUM_RUNS = 10;
    
    try {
        MonotoneAbstraction abs(is_3d);
        
        // Setup
        auto setup_start = std::chrono::high_resolution_clock::now();
        abs.setupRanges();
        abs.initializeOpenCL();
        abs.precomputeNextStateTableGPU();
        auto setup_end = std::chrono::high_resolution_clock::now();
        double setup_ms = std::chrono::duration<double, std::milli>(setup_end - setup_start).count();
        std::cout << "\nSetup: " << setup_ms << " ms\n" << std::endl;
        
        // Run multiple times
        double total_time = 0;
        int total_iter = 0;
        for (int run = 0; run < NUM_RUNS; ++run) {
            std::cout << "Run " << (run + 1) << "/" << NUM_RUNS << ": ";
            abs.initializeSafeSet();
            int iterations, basis_elements;
            double time_ms = abs.computeSafeSetGPU(iterations, basis_elements);
            total_time += time_ms;
            total_iter += iterations;
            std::cout << iterations << " iterations, " 
                      << basis_elements << " basis, " 
                      << time_ms << " ms" << std::endl;
        }
        
        std::cout << "\nAverage: " << total_iter / NUM_RUNS << " iterations, " 
                  << total_time / NUM_RUNS << " ms" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
