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

inline void setBit(unsigned char* bitmap, int idx) {
    bitmap[idx / 8] |= (1 << (idx % 8));
}

inline bool getBit(const unsigned char* bitmap, int idx) {
    return (bitmap[idx / 8] & (1 << (idx % 8))) != 0;
}

inline void clearBit(unsigned char* bitmap, int idx) {
    bitmap[idx / 8] &= ~(1 << (idx % 8));
}

// ============================================================================
// CORE IMPLEMENTATION
// ============================================================================

class MonotoneAbstraction {
private:
    double x_range_min[MAX_STATE_DIM], x_range_max[MAX_STATE_DIM];
    double x_res[MAX_STATE_DIM];
    int x_numCells[MAX_STATE_DIM];
    int x_priority[MAX_STATE_DIM];
    
    int safe_set_basis[MAX_BASIS_ELEMENTS * MAX_STATE_DIM];
    int safe_set_size;
    
    unsigned char* basis_bitmap;
    int bitmap_size;
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
        : safe_set_size(0), is_3d(use_3d), basis_bitmap(nullptr), next_state_table(nullptr),
          platform(nullptr), device(nullptr), context(nullptr), queue(nullptr), 
          program_trans(nullptr), program_safety(nullptr), kernel_trans(nullptr), kernel_safety(nullptr) {
        memset(safe_set_basis, 0, sizeof(safe_set_basis));
        max_grid_size = is_3d ? MAX_GRID_SIZE_3D : MAX_GRID_SIZE_2D;
        dt = is_3d ? DT_3D : DT_2D;
    }
    
    ~MonotoneAbstraction() {
        if (next_state_table) delete[] next_state_table;
        if (basis_bitmap) delete[] basis_bitmap;
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
        
        bitmap_size = (total_states + 7) / 8;
        basis_bitmap = new unsigned char[bitmap_size];
        memset(basis_bitmap, 0, bitmap_size);
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
        
        queue = clCreateCommandQueue(context, device, 0, &err);
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
        
        clEnqueueNDRangeKernel(queue, kernel_trans, 1, nullptr, &global_work_size, &local_work_size, 0, nullptr, nullptr);
        clEnqueueReadBuffer(queue, d_next_state_table, CL_TRUE, 0,
            sizeof(int) * total_states * MAX_STATE_DIM, next_state_table, 0, nullptr, nullptr);
        
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
        buildBitmapFromBasis();
    }
    
    void buildBitmapFromBasis() {
        memset(basis_bitmap, 0, bitmap_size);
        int state_dim = is_3d ? 3 : 2;
        
        // Only mark the actual basis elements (not all dominated states)
        // This makes early-exit extremely effective (~1000 threads vs 262K threads)
        for (int i = 0; i < safe_set_size; ++i) {
            int flat_idx = 0, multiplier = 1;
            for (int j = 0; j < state_dim; ++j) {
                flat_idx += (safe_set_basis[i * MAX_STATE_DIM + j] - 1) * multiplier;
                multiplier *= x_numCells[j];
            }
            setBit(basis_bitmap, flat_idx);
        }
    }
    
    bool xInSafeSet(const int* x_idx) const {
        if (x_idx[0] == -1) return false;
        
        int state_dim = is_3d ? 3 : 2;
        for (int i = 0; i < safe_set_size; ++i) {
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
    
    void computeSafeSetGPU() {
        int iter = 0;
        int state_dim = is_3d ? 3 : 2;
        
        // Create GPU buffers once (reuse across iterations)
        cl_int err;
        cl_mem d_basis_bitmap = clCreateBuffer(context, CL_MEM_READ_ONLY, bitmap_size, nullptr, &err);
        cl_mem d_next_state_table = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
            sizeof(int) * total_states * MAX_STATE_DIM, next_state_table, &err);
        cl_mem d_basis_list = clCreateBuffer(context, CL_MEM_READ_ONLY, 
            sizeof(int) * MAX_BASIS_ELEMENTS * MAX_STATE_DIM, nullptr, &err);
        cl_mem d_unsafe_bitmap = clCreateBuffer(context, CL_MEM_READ_WRITE, bitmap_size, nullptr, &err);
        
        unsigned char* unsafe_bitmap = new unsigned char[bitmap_size];
        
        // Pre-allocate fixed arrays (avoid vector allocations)
        int unsafe_indices[MAX_BASIS_ELEMENTS];
        int neighbor_buffer[MAX_BASIS_ELEMENTS * MAX_STATE_DIM * MAX_STATE_DIM];
        
        size_t local_work_size = 256;
        size_t global_work_size = ((total_states + local_work_size - 1) / local_work_size) * local_work_size;
        
        // Set constant kernel args once
        clSetKernelArg(kernel_safety, 1, sizeof(cl_mem), &d_next_state_table);
        clSetKernelArg(kernel_safety, 5, sizeof(int), &total_states);
        clSetKernelArg(kernel_safety, 6, sizeof(int), &state_dim);
        
        while (true) {
            iter++;
            if (safe_set_size == 0) break;
            
            // Phase 1: GPU marks unsafe states
            memset(unsafe_bitmap, 0, bitmap_size);  // Clear bitmap - kernel only sets bits, never clears
            
            clEnqueueWriteBuffer(queue, d_basis_bitmap, CL_FALSE, 0, bitmap_size, basis_bitmap, 0, nullptr, nullptr);
            clEnqueueWriteBuffer(queue, d_basis_list, CL_FALSE, 0, 
                sizeof(int) * safe_set_size * MAX_STATE_DIM, safe_set_basis, 0, nullptr, nullptr);
            clEnqueueWriteBuffer(queue, d_unsafe_bitmap, CL_FALSE, 0, bitmap_size, unsafe_bitmap, 0, nullptr, nullptr);
            
            clSetKernelArg(kernel_safety, 0, sizeof(cl_mem), &d_basis_bitmap);
            clSetKernelArg(kernel_safety, 2, sizeof(cl_mem), &d_basis_list);
            clSetKernelArg(kernel_safety, 3, sizeof(int), &safe_set_size);
            clSetKernelArg(kernel_safety, 4, sizeof(cl_mem), &d_unsafe_bitmap);
            
            clEnqueueNDRangeKernel(queue, kernel_safety, 1, nullptr, &global_work_size, &local_work_size, 0, nullptr, nullptr);
            clEnqueueReadBuffer(queue, d_unsafe_bitmap, CL_TRUE, 0, bitmap_size, unsafe_bitmap, 0, nullptr, nullptr);
            
            // Phase 2: CPU identifies unsafe basis elements (use fixed array)
            int unsafe_count = 0;
            for (int i = 0; i < safe_set_size; ++i) {
                int flat_idx = 0, multiplier = 1;
                for (int j = 0; j < state_dim; ++j) {
                    flat_idx += (safe_set_basis[i * MAX_STATE_DIM + j] - 1) * multiplier;
                    multiplier *= x_numCells[j];
                }
                if (getBit(unsafe_bitmap, flat_idx)) {
                    unsafe_indices[unsafe_count++] = i;
                }
            }
            
            // Phase 3: Collect neighbors from unsafe bases (use fixed array)
            int neighbor_count = 0;
            for (int ui = 0; ui < unsafe_count; ++ui) {
                int idx = unsafe_indices[ui];
                for (int j = 0; j < state_dim; ++j) {
                    int val = safe_set_basis[idx * MAX_STATE_DIM + j];
                    if (val == 1) continue;
                    
                    int* neighbor = &neighbor_buffer[neighbor_count * MAX_STATE_DIM];
                    for (int k = 0; k < state_dim; ++k) {
                        neighbor[k] = safe_set_basis[idx * MAX_STATE_DIM + k] - (j == k ? 1 : 0);
                    }
                    neighbor_count++;
                }
            }
            
            // Remove unsafe bases (mark-and-compact)
            int write_pos = 0;
            for (int i = 0; i < safe_set_size; ++i) {
                bool is_unsafe = false;
                for (int ui = 0; ui < unsafe_count; ++ui) {
                    if (unsafe_indices[ui] == i) {
                        is_unsafe = true;
                        break;
                    }
                }
                if (!is_unsafe) {
                    if (write_pos != i) {
                        for (int j = 0; j < state_dim; ++j) {
                            safe_set_basis[write_pos * MAX_STATE_DIM + j] = safe_set_basis[i * MAX_STATE_DIM + j];
                        }
                    }
                    write_pos++;
                }
            }
            safe_set_size = write_pos;
            
            // Add unique neighbors
            int added_count = 0;
            for (int ni = 0; ni < neighbor_count; ++ni) {
                int* neighbor = &neighbor_buffer[ni * MAX_STATE_DIM];
                if (!xInSafeSet(neighbor)) {
                    if (safe_set_size >= MAX_BASIS_ELEMENTS) break;
                    for (int j = 0; j < state_dim; ++j) {
                        safe_set_basis[safe_set_size * MAX_STATE_DIM + j] = neighbor[j];
                    }
                    safe_set_size++;
                    added_count++;
                }
            }
            
            if (added_count == 0) break;
            buildBitmapFromBasis();
        }
        
        delete[] unsafe_bitmap;
        clReleaseMemObject(d_basis_bitmap);
        clReleaseMemObject(d_next_state_table);
        clReleaseMemObject(d_basis_list);
        clReleaseMemObject(d_unsafe_bitmap);
        
        std::cout << "Converged in " << iter << " iterations with " << safe_set_size << " basis elements" << std::endl;
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
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    try {
        MonotoneAbstraction abs(is_3d);
        abs.setupRanges();
        abs.initializeOpenCL();
        abs.precomputeNextStateTableGPU();
        abs.initializeSafeSet();
        abs.computeSafeSetGPU();
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        
        std::cout << "Total: " << duration.count() << " ms" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
