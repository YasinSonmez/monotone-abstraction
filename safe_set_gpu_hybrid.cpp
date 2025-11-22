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

const double DT_2D = 0.1;
const double H_MIN_2D = 5.0, H_MAX_2D = 80.0;
const double V_MIN_2D = 0.0, V_MAX_2D = 20.0;
const double H_RES_2D = 0.2, V_RES_2D = 0.1;

const double DT_3D = 0.4;
const double H_MIN_3D = 0.0, H_MAX_3D = 80.0;
const double V_MIN_3D = 0.0, V_MAX_3D = 20.0;
const double H_RES_3D = 0.8, V_RES_3D = 0.4;

const int MAX_STATE_DIM = 3;
const int MAX_BASIS_ELEMENTS = 2000;
const int MAX_GRID_SIZE_2D = 400 * 200;
const int MAX_GRID_SIZE_3D = 101 * 51 * 51;

// ============================================================================
// UTILITIES
// ============================================================================

std::string loadKernelSource(const char* filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error(std::string("Failed to open kernel file: ") + filename);
    }
    std::stringstream buffer;
    buffer << file.rdbuf();
    return buffer.str();
}

void checkError(cl_int error, const char* message) {
    if (error != CL_SUCCESS) {
        std::cerr << "OpenCL Error (" << error << "): " << message << std::endl;
        throw std::runtime_error(message);
    }
}

// ============================================================================
// CORE IMPLEMENTATION
// ============================================================================

struct ComputeStats {
    int iterations;
    int basis_elements;
    double upload_ms;
    double kernel_ms;
    double readback_ms;
    double cpu_ms;
    double total_ms;
};

class MonotoneAbstraction {
private:
    double x_range_min[MAX_STATE_DIM], x_range_max[MAX_STATE_DIM];
    double x_res[MAX_STATE_DIM];
    int x_numCells[MAX_STATE_DIM];
    int x_priority[MAX_STATE_DIM];
    
    int safe_set_basis[MAX_BASIS_ELEMENTS * MAX_STATE_DIM];
    int safe_set_size;
    int safe_set_flat_indices[MAX_BASIS_ELEMENTS];
    int flatten_multipliers[MAX_STATE_DIM];
    int total_states;
    
    int* next_state_table;
    int max_grid_size;
    
    bool is_3d;
    double dt;
    
    cl_platform_id platform;
    cl_device_id device;
    cl_context context;
    cl_command_queue queue;
    cl_program program_trans, program_safety;
    cl_kernel kernel_trans, kernel_safety;
    
public:
        MonotoneAbstraction(bool use_3d) 
                : safe_set_size(0), is_3d(use_3d), next_state_table(nullptr),
                    platform(nullptr), device(nullptr), context(nullptr), queue(nullptr), 
          program_trans(nullptr), program_safety(nullptr), kernel_trans(nullptr), kernel_safety(nullptr) {
        memset(safe_set_basis, 0, sizeof(safe_set_basis));
                memset(safe_set_flat_indices, 0, sizeof(safe_set_flat_indices));
                memset(flatten_multipliers, 0, sizeof(flatten_multipliers));
        max_grid_size = is_3d ? MAX_GRID_SIZE_3D : MAX_GRID_SIZE_2D;
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
        if (is_3d) {
            x_range_min[0] = H_MIN_3D; x_range_max[0] = H_MAX_3D;
            x_range_min[1] = V_MIN_3D; x_range_max[1] = V_MAX_3D;
            x_range_min[2] = V_MIN_3D; x_range_max[2] = V_MAX_3D;
            x_res[0] = H_RES_3D; x_res[1] = V_RES_3D; x_res[2] = V_RES_3D;
            x_priority[0] = 0; x_priority[1] = 1; x_priority[2] = 0;
        } else {
            x_range_min[0] = H_MIN_2D; x_range_max[0] = H_MAX_2D;
            x_range_min[1] = V_MIN_2D; x_range_max[1] = V_MAX_2D;
            x_res[0] = H_RES_2D; x_res[1] = V_RES_2D;
            x_priority[0] = 0; x_priority[1] = 1;
        }
        
        int state_dim = is_3d ? 3 : 2;
        for (int i = 0; i < state_dim; ++i) {
            x_numCells[i] = std::ceil((x_range_max[i] - x_range_min[i]) / x_res[i]) + 1;
        }
        
        total_states = x_numCells[0];
        for (int i = 1; i < state_dim; ++i) {
            total_states *= x_numCells[i];
        }
        
        flatten_multipliers[0] = 1;
        for (int i = 1; i < state_dim; ++i) {
            flatten_multipliers[i] = flatten_multipliers[i - 1] * x_numCells[i - 1];
        }
    }
    
    void initializeOpenCL() {
        cl_int err;
        
        err = clGetPlatformIDs(1, &platform, nullptr);
        checkError(err, "Failed to get platform");
        
        err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 1, &device, nullptr);
        if (err != CL_SUCCESS) {
            err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_CPU, 1, &device, nullptr);
            checkError(err, "Failed to get device");
        }
        
        context = clCreateContext(nullptr, 1, &device, nullptr, nullptr, &err);
        checkError(err, "Failed to create context");
        
    queue = clCreateCommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, &err);
        checkError(err, "Failed to create command queue");
        
        // Load transition kernel
        std::string kernel_trans_src = loadKernelSource("precompute_transitions.cl");
        const char* trans_str = kernel_trans_src.c_str();
        size_t trans_size = kernel_trans_src.length();
        program_trans = clCreateProgramWithSource(context, 1, &trans_str, &trans_size, &err);
        checkError(err, "Failed to create transition program");
        
        err = clBuildProgram(program_trans, 1, &device, nullptr, nullptr, nullptr);
        if (err != CL_SUCCESS) {
            size_t log_size;
            clGetProgramBuildInfo(program_trans, device, CL_PROGRAM_BUILD_LOG, 0, nullptr, &log_size);
            char* log = new char[log_size];
            clGetProgramBuildInfo(program_trans, device, CL_PROGRAM_BUILD_LOG, log_size, log, nullptr);
            std::cerr << "Transition kernel build log:\n" << log << std::endl;
            delete[] log;
            checkError(err, "Failed to build transition program");
        }
        
        kernel_trans = clCreateKernel(program_trans, "precompute_transitions", &err);
        checkError(err, "Failed to create transition kernel");
        
        // Load safety check kernel
        std::string kernel_safety_src = loadKernelSource("check_basis_safety.cl");
        const char* safety_str = kernel_safety_src.c_str();
        size_t safety_size = kernel_safety_src.length();
        program_safety = clCreateProgramWithSource(context, 1, &safety_str, &safety_size, &err);
        checkError(err, "Failed to create safety program");
        
        err = clBuildProgram(program_safety, 1, &device, nullptr, nullptr, nullptr);
        if (err != CL_SUCCESS) {
            size_t log_size;
            clGetProgramBuildInfo(program_safety, device, CL_PROGRAM_BUILD_LOG, 0, nullptr, &log_size);
            char* log = new char[log_size];
            clGetProgramBuildInfo(program_safety, device, CL_PROGRAM_BUILD_LOG, log_size, log, nullptr);
            std::cerr << "Safety kernel build log:\n" << log << std::endl;
            delete[] log;
            checkError(err, "Failed to build safety program");
        }
        
        kernel_safety = clCreateKernel(program_safety, "check_basis_safety", &err);
        checkError(err, "Failed to create safety kernel");
    }
    
    void precomputeNextStateTableGPU() {
        int state_dim = is_3d ? 3 : 2;
        next_state_table = new int[total_states * MAX_STATE_DIM];
        
        // Try to load from cache first
        const char* cache_file = is_3d ? "transition_cache_3d.bin" : "transition_cache_2d.bin";
        std::ifstream cache_in(cache_file, std::ios::binary);
        
        if (cache_in.good()) {
            int cached_total_states;
            cache_in.read(reinterpret_cast<char*>(&cached_total_states), sizeof(int));
            
            if (cached_total_states == total_states) {
                cache_in.read(reinterpret_cast<char*>(next_state_table), 
                             sizeof(int) * total_states * MAX_STATE_DIM);
                cache_in.close();
                return;
            }
            cache_in.close();
        }
        
        // Cache miss - compute on GPU
        cl_int err;
        cl_mem d_x_range_min = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
            sizeof(double) * MAX_STATE_DIM, x_range_min, &err);
        cl_mem d_x_range_max = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
            sizeof(double) * MAX_STATE_DIM, x_range_max, &err);
        cl_mem d_x_res = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
            sizeof(double) * MAX_STATE_DIM, x_res, &err);
        cl_mem d_x_numCells = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
            sizeof(int) * MAX_STATE_DIM, x_numCells, &err);
        cl_mem d_x_priority = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
            sizeof(int) * MAX_STATE_DIM, x_priority, &err);
        cl_mem d_next_state_table = clCreateBuffer(context, CL_MEM_WRITE_ONLY,
            sizeof(int) * total_states * MAX_STATE_DIM, nullptr, &err);
        
        int is_3d_int = is_3d ? 1 : 0;
        clSetKernelArg(kernel_trans, 0, sizeof(int), &total_states);
        clSetKernelArg(kernel_trans, 1, sizeof(int), &state_dim);
        clSetKernelArg(kernel_trans, 2, sizeof(int), &is_3d_int);
        clSetKernelArg(kernel_trans, 3, sizeof(double), &dt);
        clSetKernelArg(kernel_trans, 4, sizeof(cl_mem), &d_x_range_min);
        clSetKernelArg(kernel_trans, 5, sizeof(cl_mem), &d_x_range_max);
        clSetKernelArg(kernel_trans, 6, sizeof(cl_mem), &d_x_res);
        clSetKernelArg(kernel_trans, 7, sizeof(cl_mem), &d_x_numCells);
        clSetKernelArg(kernel_trans, 8, sizeof(cl_mem), &d_x_priority);
        clSetKernelArg(kernel_trans, 9, sizeof(cl_mem), &d_next_state_table);
        
        size_t global_work_size = total_states;
        size_t local_work_size = 256;
        global_work_size = ((global_work_size + local_work_size - 1) / local_work_size) * local_work_size;
        
        cl_event trans_event;
        clEnqueueNDRangeKernel(queue, kernel_trans, 1, nullptr, &global_work_size, &local_work_size, 0, nullptr, &trans_event);
        clWaitForEvents(1, &trans_event);
        clEnqueueReadBuffer(queue, d_next_state_table, CL_TRUE, 0,
            sizeof(int) * total_states * MAX_STATE_DIM, next_state_table, 0, nullptr, nullptr);
    cl_ulong trans_start = 0, trans_end = 0;
    clGetEventProfilingInfo(trans_event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &trans_start, nullptr);
    clGetEventProfilingInfo(trans_event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &trans_end, nullptr);
    double trans_ms = (trans_end - trans_start) * 1e-6;
        clReleaseEvent(trans_event);
    std::cout << "Transition kernel: " << trans_ms << " ms" << std::endl;
        
        clReleaseMemObject(d_x_range_min);
        clReleaseMemObject(d_x_range_max);
        clReleaseMemObject(d_x_res);
        clReleaseMemObject(d_x_numCells);
        clReleaseMemObject(d_x_priority);
        clReleaseMemObject(d_next_state_table);
        
        // Save to cache
        std::ofstream cache_out(cache_file, std::ios::binary);
        if (cache_out.good()) {
            cache_out.write(reinterpret_cast<const char*>(&total_states), sizeof(int));
            cache_out.write(reinterpret_cast<const char*>(next_state_table), 
                           sizeof(int) * total_states * MAX_STATE_DIM);
            cache_out.close();
        }
    }
    
    void initializeSafeSet() {
        safe_set_size = 0;
        int state_dim = is_3d ? 3 : 2;
        for (int i = 0; i < state_dim; ++i) {
            safe_set_basis[safe_set_size * MAX_STATE_DIM + i] = x_numCells[i];
        }
        safe_set_size = 1;
        safe_set_flat_indices[0] = flattenIndex(&safe_set_basis[0]);
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
    
    ComputeStats computeSafeSetGPU() {
        int iter = 0;
        int state_dim = is_3d ? 3 : 2;

        cl_int err;
        cl_mem d_next_state_table = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
            sizeof(int) * total_states * MAX_STATE_DIM, next_state_table, &err);
        cl_mem d_basis_list = clCreateBuffer(context, CL_MEM_READ_ONLY,
            sizeof(int) * MAX_BASIS_ELEMENTS * MAX_STATE_DIM, nullptr, &err);
        cl_mem d_basis_flat = clCreateBuffer(context, CL_MEM_READ_ONLY,
            sizeof(int) * MAX_BASIS_ELEMENTS, nullptr, &err);
        cl_mem d_unsafe_flags = clCreateBuffer(context, CL_MEM_READ_WRITE,
            sizeof(int) * MAX_BASIS_ELEMENTS, nullptr, &err);

    int unsafe_flags[MAX_BASIS_ELEMENTS];
        int unsafe_indices[MAX_BASIS_ELEMENTS];
        int neighbor_buffer[MAX_BASIS_ELEMENTS * MAX_STATE_DIM * MAX_STATE_DIM];
    unsigned char unsafe_mask[MAX_BASIS_ELEMENTS];
        unsigned char* seen_neighbors = new unsigned char[total_states];

        size_t local_work_size = 128;

        clSetKernelArg(kernel_safety, 0, sizeof(cl_mem), &d_basis_flat);
        clSetKernelArg(kernel_safety, 1, sizeof(cl_mem), &d_next_state_table);
        clSetKernelArg(kernel_safety, 2, sizeof(cl_mem), &d_basis_list);
        clSetKernelArg(kernel_safety, 4, sizeof(cl_mem), &d_unsafe_flags);
        clSetKernelArg(kernel_safety, 5, sizeof(int), &total_states);
        clSetKernelArg(kernel_safety, 6, sizeof(int), &state_dim);

        double total_upload_ms = 0.0;
        double total_kernel_ms = 0.0;
        double total_read_ms = 0.0;
        double total_cpu_ms = 0.0;
        
        auto compute_start = std::chrono::high_resolution_clock::now();

        while (true) {
            iter++;
            if (safe_set_size == 0) break;

            size_t basis_bytes = sizeof(int) * safe_set_size * state_dim;
            size_t flat_bytes = sizeof(int) * safe_set_size;

            auto upload_start = std::chrono::high_resolution_clock::now();
            cl_event write_events[2];
            clEnqueueWriteBuffer(queue, d_basis_list, CL_FALSE, 0, basis_bytes, safe_set_basis, 0, nullptr, &write_events[0]);
            clEnqueueWriteBuffer(queue, d_basis_flat, CL_FALSE, 0, flat_bytes, safe_set_flat_indices, 0, nullptr, &write_events[1]);
            clWaitForEvents(2, write_events);
            auto upload_end = std::chrono::high_resolution_clock::now();
            total_upload_ms += std::chrono::duration<double, std::milli>(upload_end - upload_start).count();
            for (int we = 0; we < 2; ++we) {
                clReleaseEvent(write_events[we]);
            }

            clSetKernelArg(kernel_safety, 3, sizeof(int), &safe_set_size);

            size_t global_work_size = ((safe_set_size + local_work_size - 1) / local_work_size) * local_work_size;
            cl_event safety_event;
            clEnqueueNDRangeKernel(queue, kernel_safety, 1, nullptr, &global_work_size, &local_work_size, 0, nullptr, &safety_event);
            clWaitForEvents(1, &safety_event);
            cl_ulong safety_start = 0, safety_end = 0;
            clGetEventProfilingInfo(safety_event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &safety_start, nullptr);
            clGetEventProfilingInfo(safety_event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &safety_end, nullptr);
            total_kernel_ms += (safety_end - safety_start) * 1e-6;
            clReleaseEvent(safety_event);

            auto read_start = std::chrono::high_resolution_clock::now();
            clEnqueueReadBuffer(queue, d_unsafe_flags, CL_TRUE, 0, flat_bytes, unsafe_flags, 0, nullptr, nullptr);
            auto read_end = std::chrono::high_resolution_clock::now();
            total_read_ms += std::chrono::duration<double, std::milli>(read_end - read_start).count();

            auto cpu_start = std::chrono::high_resolution_clock::now();

            std::fill_n(unsafe_mask, safe_set_size, static_cast<unsigned char>(0));
            memset(seen_neighbors, 0, total_states);

            int unsafe_count = 0;
            int neighbor_count = 0;
            for (int i = 0; i < safe_set_size; ++i) {
                if (!unsafe_flags[i]) continue;

                unsafe_indices[unsafe_count++] = i;
                unsafe_mask[i] = 1;

                int idx = i;
                for (int j = 0; j < state_dim; ++j) {
                    int val = safe_set_basis[idx * MAX_STATE_DIM + j];
                    if (val == 1) continue;

                    int neighbor_coords[MAX_STATE_DIM];
                    for (int k = 0; k < state_dim; ++k) {
                        neighbor_coords[k] = safe_set_basis[idx * MAX_STATE_DIM + k] - (j == k ? 1 : 0);
                    }
                    
                    int flat_idx = flattenIndex(neighbor_coords);
                    if (!seen_neighbors[flat_idx]) {
                        seen_neighbors[flat_idx] = 1;
                        int* neighbor = &neighbor_buffer[neighbor_count * MAX_STATE_DIM];
                        for (int k = 0; k < state_dim; ++k) {
                            neighbor[k] = neighbor_coords[k];
                        }
                        neighbor_count++;
                    }
                }
            }

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

            int added_count = 0;
            for (int ni = 0; ni < neighbor_count; ++ni) {
                int* neighbor = &neighbor_buffer[ni * MAX_STATE_DIM];
                if (!xInSafeSet(neighbor)) {
                    if (safe_set_size >= MAX_BASIS_ELEMENTS) break;
                    for (int j = 0; j < state_dim; ++j) {
                        safe_set_basis[safe_set_size * MAX_STATE_DIM + j] = neighbor[j];
                    }
                    safe_set_flat_indices[safe_set_size] = flattenIndex(neighbor);
                    safe_set_size++;
                    added_count++;
                }
            }

            auto cpu_end = std::chrono::high_resolution_clock::now();
            total_cpu_ms += std::chrono::duration<double, std::milli>(cpu_end - cpu_start).count();

            if (added_count == 0) {
                break;
            }
        }
        
        auto compute_end = std::chrono::high_resolution_clock::now();
        double total_compute_ms = std::chrono::duration<double, std::milli>(compute_end - compute_start).count();

        delete[] seen_neighbors;
        
        clReleaseMemObject(d_next_state_table);
        clReleaseMemObject(d_basis_list);
        clReleaseMemObject(d_basis_flat);
        clReleaseMemObject(d_unsafe_flags);

        return ComputeStats{iter, safe_set_size, total_upload_ms, total_kernel_ms, 
                           total_read_ms, total_cpu_ms, total_compute_ms};
    }
    
    void printSafeSet(int max_elements = 10) const {
        std::cout << "\nSafe Set Basis Elements:" << std::endl;
        int state_dim = is_3d ? 3 : 2;
        for (int i = 0; i < std::min(max_elements, safe_set_size); ++i) {
            double state[MAX_STATE_DIM];
            int basis[MAX_STATE_DIM];
            for (int j = 0; j < state_dim; ++j) {
                basis[j] = safe_set_basis[i * MAX_STATE_DIM + j];
            }
            getPriorityStateAtIdx(basis, state);
            
            std::cout << "Basis " << i << ": [";
            for (int j = 0; j < state_dim; ++j) {
                std::cout << state[j];
                if (j < state_dim - 1) std::cout << ", ";
            }
            std::cout << "]" << std::endl;
        }
    }
    
private:
    void getPriorityStateAtIdx(const int* idx, double* val) const {
        int state_dim = is_3d ? 3 : 2;
        for (int i = 0; i < state_dim; ++i) {
            if (x_priority[i] == 1) {
                val[i] = x_range_min[i] + (idx[i] - 1) * x_res[i];
            } else {
                val[i] = x_range_max[i] - (idx[i] - 1) * x_res[i];
            }
        }
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
    
    void unflattenIndex(int flat_idx, int* idx) const {
        int state_dim = is_3d ? 3 : 2;
        for (int i = 0; i < state_dim; ++i) {
            idx[i] = (flat_idx % x_numCells[i]) + 1;
            flat_idx /= x_numCells[i];
        }
    }
};

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cout << "Usage: " << argv[0] << " [2d|3d]" << std::endl;
        return 1;
    }
    
    std::string mode = argv[1];
    bool is_3d = (mode == "3d");
    
    if (mode != "2d" && mode != "3d") {
        std::cout << "Error: Mode must be '2d' or '3d'" << std::endl;
        return 1;
    }
    
    const int NUM_RUNS = 10;
    
    try {
        MonotoneAbstraction abs(is_3d);
        
        // One-time setup
        auto setup_start = std::chrono::high_resolution_clock::now();
        abs.setupRanges();
        abs.initializeOpenCL();
        abs.precomputeNextStateTableGPU();
        auto setup_end = std::chrono::high_resolution_clock::now();
        double setup_ms = std::chrono::duration<double, std::milli>(setup_end - setup_start).count();
        
        std::cout << "\n==================================" << std::endl;
        std::cout << "Setup phase: " << setup_ms << " ms" << std::endl;
        std::cout << "==================================" << std::endl;
        
        // Run computation 10 times
        std::vector<ComputeStats> all_stats;
        
        for (int run = 0; run < NUM_RUNS; ++run) {
            std::cout << "\n==================================" << std::endl;
            std::cout << "Run " << (run + 1) << "/" << NUM_RUNS << std::endl;
            std::cout << "==================================" << std::endl;
            
            auto run_start = std::chrono::high_resolution_clock::now();
            abs.initializeSafeSet();
            auto stats = abs.computeSafeSetGPU();
            auto run_end = std::chrono::high_resolution_clock::now();
            double run_total_ms = std::chrono::duration<double, std::milli>(run_end - run_start).count();
            
            all_stats.push_back(stats);
            
            std::cout << "  Iterations: " << stats.iterations << std::endl;
            std::cout << "  Basis elements: " << stats.basis_elements << std::endl;
            std::cout << "  Upload time: " << stats.upload_ms << " ms" << std::endl;
            std::cout << "  Kernel time: " << stats.kernel_ms << " ms" << std::endl;
            std::cout << "  Readback time: " << stats.readback_ms << " ms" << std::endl;
            std::cout << "  CPU update time: " << stats.cpu_ms << " ms" << std::endl;
            std::cout << "  Compute total: " << stats.total_ms << " ms" << std::endl;
            std::cout << "  Run total (with init): " << run_total_ms << " ms" << std::endl;
        }
        
        // Compute and display averages
        std::cout << "\n==================================" << std::endl;
        std::cout << "Summary Statistics (avg over " << NUM_RUNS << " runs)" << std::endl;
        std::cout << "==================================" << std::endl;
        
        double avg_upload = 0.0, avg_kernel = 0.0, avg_readback = 0.0;
        double avg_cpu = 0.0, avg_total = 0.0;
        int avg_iterations = 0;
        
        for (const auto& stats : all_stats) {
            avg_upload += stats.upload_ms;
            avg_kernel += stats.kernel_ms;
            avg_readback += stats.readback_ms;
            avg_cpu += stats.cpu_ms;
            avg_total += stats.total_ms;
            avg_iterations += stats.iterations;
        }
        
        avg_upload /= NUM_RUNS;
        avg_kernel /= NUM_RUNS;
        avg_readback /= NUM_RUNS;
        avg_cpu /= NUM_RUNS;
        avg_total /= NUM_RUNS;
        avg_iterations /= NUM_RUNS;
        
        std::cout << "  Avg iterations: " << avg_iterations << std::endl;
        std::cout << "  Avg upload time: " << avg_upload << " ms" << std::endl;
        std::cout << "  Avg kernel time: " << avg_kernel << " ms" << std::endl;
        std::cout << "  Avg readback time: " << avg_readback << " ms" << std::endl;
        std::cout << "  Avg CPU update time: " << avg_cpu << " ms" << std::endl;
        std::cout << "  Avg compute total: " << avg_total << " ms" << std::endl;
        std::cout << "  Setup time: " << setup_ms << " ms" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
