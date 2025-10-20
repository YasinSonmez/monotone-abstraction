#include <iostream>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <cstring>
#include <thread>
#include <mutex>
#include <atomic>
#include <unordered_map>

/**
 * Cached Safe Set Computation for Vehicle Following
 * Implements full dynamics lookup table for maximum performance
 */

// ============================================================================
// CONFIGURATION
// ============================================================================

// 2D Configuration (matching MATLAB)
const double DT_2D = 0.1;
const double H_MIN_2D = 5.0, H_MAX_2D = 80.0;
const double V_MIN_2D = 0.0, V_MAX_2D = 20.0;
const double H_RES_2D = 0.2, V_RES_2D = 0.1;
const double U_MIN_2D = -3.0, U_MAX_2D = 3.0;
const double U_RES_2D = 0.1;

// 3D Configuration  
const double DT_3D = 0.4;
const double H_MIN_3D = 0.0, H_MAX_3D = 80.0;
const double V_MIN_3D = 0.0, V_MAX_3D = 20.0;
const double H_RES_3D = 0.8, V_RES_3D = 0.4;
const double T_RES_3D = 150.0;

// 3D Vehicle parameters
const double T_MAX = 1200.0, T_BRAKE_MIN = -1800.0, T_BRAKE_MAX = -2400.0;
const double M_MIN = 2000.0, M_MAX = 2500.0;
const double R_W_MIN = 0.30, R_W_MAX = 0.35;
const double ALPHA_MIN = 300.0, ALPHA_MAX = 350.0;
const double BETA_MIN = 0.10, BETA_MAX = 0.25;
const double GAMMA_MIN = 0.30, GAMMA_MAX = 0.65;

// Maximum sizes for pre-allocation
const int MAX_STATE_DIM = 3;
const int MAX_INPUT_DIM = 1;
const int MAX_DISTURBANCE_DIM = 1;
const int MAX_BASIS_ELEMENTS = 2000;

// Optimization: Reduce RK4 steps for faster computation
const int RK4_STEPS_OPTIMIZED = 200;  // Reduced from 1000 to 200 (5x faster)

// ============================================================================
// HASH FUNCTION FOR CACHE KEYS
// ============================================================================

struct TransitionKey {
    int x_idx[MAX_STATE_DIM];
    int u_idx[MAX_INPUT_DIM];
    int w_idx[MAX_DISTURBANCE_DIM];
    int state_dim;
    
    bool operator==(const TransitionKey& other) const {
        for (int i = 0; i < state_dim; ++i) {
            if (x_idx[i] != other.x_idx[i]) return false;
        }
        for (int i = 0; i < MAX_INPUT_DIM; ++i) {
            if (u_idx[i] != other.u_idx[i]) return false;
        }
        for (int i = 0; i < MAX_DISTURBANCE_DIM; ++i) {
            if (w_idx[i] != other.w_idx[i]) return false;
        }
        return true;
    }
};

struct TransitionKeyHash {
    std::size_t operator()(const TransitionKey& key) const {
        std::size_t hash = 0;
        for (int i = 0; i < key.state_dim; ++i) {
            hash ^= std::hash<int>()(key.x_idx[i]) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        }
        for (int i = 0; i < MAX_INPUT_DIM; ++i) {
            hash ^= std::hash<int>()(key.u_idx[i]) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        }
        for (int i = 0; i < MAX_DISTURBANCE_DIM; ++i) {
            hash ^= std::hash<int>()(key.w_idx[i]) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        }
        return hash;
    }
};

// ============================================================================
// CORE IMPLEMENTATION
// ============================================================================

class MonotoneDynamics {
private:
    double dt;
    bool is_3d;
    
public:
    MonotoneDynamics(int /*state_dim*/, int /*input_dim*/, int /*disturbance_dim*/, double timestep, bool use_3d)
        : dt(timestep), is_3d(use_3d) {}
    
    void nextState(const double* x, const double* u, const double* w, double* x_plus) const {
        if (is_3d) {
            solve3DDynamics(x, u, w, x_plus);
        } else {
            // 2D dynamics: h(t+dt) = h(t) + v_L*dt - v(t)*dt - 0.5*a*dt^2, v(t+dt) = v(t) + a*dt
            x_plus[0] = x[0] + w[0]*dt - x[1]*dt - 0.5*u[0]*dt*dt;
            x_plus[1] = x[1] + u[0]*dt;
        }
    }
    
private:
    void solve3DDynamics(const double* x, const double* u, const double* w, double* x_plus) const {
        // Choose worst-case parameters for ego vehicle
        double R_w_ego = (u[0] > 0) ? R_W_MIN : R_W_MAX;
        double M_ego = ((u[0] / R_w_ego - ALPHA_MIN) > 0) ? M_MIN : M_MAX;
        double a = (1.0 / M_ego) * (u[0] / R_w_ego - ALPHA_MIN);
        double b = -(1.0 / M_MAX) * BETA_MIN;
        double c = -(1.0 / M_MAX) * GAMMA_MIN;
        
        // Choose worst-case parameters for lead vehicle
        double R_w_L = (w[0] > 0) ? R_W_MAX : R_W_MIN;
        double M_L = ((w[0] / R_w_L - ALPHA_MAX) > 0) ? M_MAX : M_MIN;
        double a_L = (1.0 / M_L) * (w[0] / R_w_L - ALPHA_MAX);
        double b_L = -(1.0 / M_MIN) * BETA_MAX;
        double c_L = -(1.0 / M_MIN) * GAMMA_MAX;
        
        // Solve vehicle dynamics using RK4 with reduced steps
        double pos_plus, vel_plus, posL_plus, velL_plus;
        solveODE(x[1], a, b, c, false, pos_plus, vel_plus);
        solveODE(x[2], a_L, b_L, c_L, true, posL_plus, velL_plus);
        
        x_plus[0] = x[0] + posL_plus - pos_plus;
        x_plus[1] = std::max(V_MIN_3D, std::min(vel_plus, V_MAX_3D));
        x_plus[2] = std::max(V_MIN_3D, std::min(velL_plus, V_MAX_3D));
    }
    
    void solveODE(double v0, double a, double b, double c, bool is_lead, double& pos, double& vel) const {
        const int num_steps = RK4_STEPS_OPTIMIZED;  // Reduced from 1000 to 200
        const double step_size = dt / num_steps;
        pos = 0; 
        vel = v0;
        
        for (int i = 0; i < num_steps; ++i) {
            double k1_pos, k1_vel, k2_pos, k2_vel, k3_pos, k3_vel, k4_pos, k4_vel;
            
            f(pos, vel, a, b, c, is_lead, k1_pos, k1_vel);
            f(pos + step_size/2 * k1_pos, vel + step_size/2 * k1_vel, a, b, c, is_lead, k2_pos, k2_vel);
            f(pos + step_size/2 * k2_pos, vel + step_size/2 * k2_vel, a, b, c, is_lead, k3_pos, k3_vel);
            f(pos + step_size * k3_pos, vel + step_size * k3_vel, a, b, c, is_lead, k4_pos, k4_vel);
            
            pos += step_size/6 * (k1_pos + 2*k2_pos + 2*k3_pos + k4_pos);
            vel += step_size/6 * (k1_vel + 2*k2_vel + 2*k3_vel + k4_vel);
        }
    }
    
    void f(double /*pos*/, double vel, double a, double b, double c, bool is_lead, double& dpos, double& dvel) const {
        double dvdt = a + b * vel + c * vel * vel;
        if (vel <= 0.0 && dvdt < 0.0) dvdt = 0.0;
        if (is_lead && vel >= V_MAX_3D && dvdt > 0.0) dvdt = 0.0;
        dpos = vel;
        dvel = dvdt;
    }
};

class MonotoneAbstraction {
private:
    MonotoneDynamics dyn;
    
    // Fixed-size arrays instead of STL containers
    double x_range_min[MAX_STATE_DIM], x_range_max[MAX_STATE_DIM];
    double u_range_min[MAX_INPUT_DIM], u_range_max[MAX_INPUT_DIM];
    double w_range_min[MAX_DISTURBANCE_DIM], w_range_max[MAX_DISTURBANCE_DIM];
    double x_res[MAX_STATE_DIM], u_res[MAX_INPUT_DIM], w_res[MAX_DISTURBANCE_DIM];
    int x_numCells[MAX_STATE_DIM], u_numCells[MAX_INPUT_DIM], w_numCells[MAX_DISTURBANCE_DIM];
    int x_priority[MAX_STATE_DIM], u_priority[MAX_INPUT_DIM], w_priority[MAX_DISTURBANCE_DIM];
    
    // Safe set basis as flattened 2D array: [basis_index][dimension]
    int safe_set_basis[MAX_BASIS_ELEMENTS * MAX_STATE_DIM];
    int safe_set_size;  // Current number of basis elements
    
    // NEW: Full dynamics lookup table using hash map
    std::unordered_map<TransitionKey, int*, TransitionKeyHash> dynamics_cache;
    mutable std::mutex cache_mutex;
    
    // Statistics
    mutable std::atomic<long long> cache_hits{0};
    mutable std::atomic<long long> cache_misses{0};
    
    bool is_3d;
    bool use_parallel;
    
public:
    MonotoneAbstraction(const MonotoneDynamics& dynamics, bool use_3d, bool parallel = false) 
        : dyn(dynamics), safe_set_size(0), is_3d(use_3d), use_parallel(parallel) {
        // Initialize arrays to zero
        memset(safe_set_basis, 0, sizeof(safe_set_basis));
    }
    
    ~MonotoneAbstraction() {
        // Clean up dynamically allocated cache entries
        for (auto& entry : dynamics_cache) {
            delete[] entry.second;
        }
    }
    
    void setupRanges() {
        if (is_3d) {
            x_range_min[0] = H_MIN_3D; x_range_max[0] = H_MAX_3D;
            x_range_min[1] = V_MIN_3D; x_range_max[1] = V_MAX_3D;
            x_range_min[2] = V_MIN_3D; x_range_max[2] = V_MAX_3D;
            u_range_min[0] = T_BRAKE_MIN; u_range_max[0] = T_MAX;
            w_range_min[0] = T_BRAKE_MAX; w_range_max[0] = T_MAX;
            x_res[0] = H_RES_3D; x_res[1] = V_RES_3D; x_res[2] = V_RES_3D;
            u_res[0] = T_RES_3D;
            w_res[0] = T_RES_3D;
            x_priority[0] = 0; x_priority[1] = 1; x_priority[2] = 0;
            u_priority[0] = 1;
            w_priority[0] = 0;
        } else {
            x_range_min[0] = H_MIN_2D; x_range_max[0] = H_MAX_2D;
            x_range_min[1] = V_MIN_2D; x_range_max[1] = V_MAX_2D;
            u_range_min[0] = U_MIN_2D; u_range_max[0] = U_MAX_2D;
            w_range_min[0] = V_MIN_2D; w_range_max[0] = V_MAX_2D;
            x_res[0] = H_RES_2D; x_res[1] = V_RES_2D;
            u_res[0] = U_RES_2D;
            w_res[0] = V_RES_2D;
            x_priority[0] = 0; x_priority[1] = 1;
            u_priority[0] = 1;
            w_priority[0] = 0;
        }
        
        // Compute number of cells
        int state_dim = is_3d ? 3 : 2;
        for (int i = 0; i < state_dim; ++i) {
            x_numCells[i] = std::ceil((x_range_max[i] - x_range_min[i]) / x_res[i]) + 1;
        }
        for (int i = 0; i < MAX_INPUT_DIM; ++i) {
            u_numCells[i] = std::ceil((u_range_max[i] - u_range_min[i]) / u_res[i]) + 1;
        }
        for (int i = 0; i < MAX_DISTURBANCE_DIM; ++i) {
            w_numCells[i] = std::ceil((w_range_max[i] - w_range_min[i]) / w_res[i]) + 1;
        }
    }
    
    void initializeSafeSet() {
        safe_set_size = 0;
        if (is_3d) {
            // 3D: Start empty for precomputation
        } else {
            // 2D: Initialize with entire state space
            int state_dim = 2;
            for (int i = 0; i < state_dim; ++i) {
                safe_set_basis[safe_set_size * MAX_STATE_DIM + i] = x_numCells[i];
            }
            safe_set_size = 1;
        }
    }
    
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
    
    bool getStateIdx(const double* val, int* idx) const {
        int state_dim = is_3d ? 3 : 2;
        const double threshold = 1e-5;
        
        for (int i = 0; i < state_dim; ++i) {
            if (val[i] < x_range_min[i] - threshold || val[i] > x_range_max[i] + threshold) {
                idx[0] = -1;  // Invalid state
                return false;
            }
            
            if (x_priority[i] == 1) {
                idx[i] = std::ceil((val[i] - x_range_min[i]) / x_res[i] - threshold) + 1;
            } else {
                idx[i] = std::ceil((x_range_max[i] - val[i]) / x_res[i] - threshold) + 1;
            }
        }
        return true;
    }
    
    bool xInSafeSet(const int* x_idx) const {
        if (x_idx[0] == -1) return false;
        
        int state_dim = is_3d ? 3 : 2;
        // Check most likely candidates first (reverse order)
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
    
    // NEW: Improved transition state with proper caching
    void getTransitionState(const int* x_idx, const int* u_idx, const int* w_idx, int* x_plus_idx) {
        int state_dim = is_3d ? 3 : 2;
        
        // Create cache key
        TransitionKey key;
        key.state_dim = state_dim;
        for (int i = 0; i < state_dim; ++i) {
            key.x_idx[i] = x_idx[i];
        }
        for (int i = 0; i < MAX_INPUT_DIM; ++i) {
            key.u_idx[i] = u_idx[i];
        }
        for (int i = 0; i < MAX_DISTURBANCE_DIM; ++i) {
            key.w_idx[i] = w_idx[i];
        }
        
        // Check cache first (thread-safe)
        {
            std::lock_guard<std::mutex> lock(cache_mutex);
            auto it = dynamics_cache.find(key);
            if (it != dynamics_cache.end()) {
                cache_hits++;
                for (int j = 0; j < state_dim; ++j) {
                    x_plus_idx[j] = it->second[j];
                }
                return;
            }
        }
        
        cache_misses++;
        
        // Compute transition
        double x_val[MAX_STATE_DIM], u_val[MAX_INPUT_DIM], w_val[MAX_DISTURBANCE_DIM];
        double x_plus_val[MAX_STATE_DIM];
        
        getPriorityStateAtIdx(x_idx, x_val);
        getPriorityInputAtIdx(u_idx, u_val);
        getPriorityDisturbanceAtIdx(w_idx, w_val);
        
        dyn.nextState(x_val, u_val, w_val, x_plus_val);
        getStateIdx(x_plus_val, x_plus_idx);
        
        // Cache the result (thread-safe)
        {
            std::lock_guard<std::mutex> lock(cache_mutex);
            int* cached_result = new int[MAX_STATE_DIM];
            for (int j = 0; j < state_dim; ++j) {
                cached_result[j] = x_plus_idx[j];
            }
            dynamics_cache[key] = cached_result;
        }
    }
    
    void removeBasisIdx(const int* basis_idx) {
        int state_dim = is_3d ? 3 : 2;
        
        // Find and remove the basis element
        for (int i = 0; i < safe_set_size; ++i) {
            bool match = true;
            for (int j = 0; j < state_dim; ++j) {
                if (safe_set_basis[i * MAX_STATE_DIM + j] != basis_idx[j]) {
                    match = false;
                    break;
                }
            }
            if (match) {
                // Use memmove for faster array shifting
                if (i < safe_set_size - 1) {
                    memmove(&safe_set_basis[i * MAX_STATE_DIM], 
                           &safe_set_basis[(i + 1) * MAX_STATE_DIM],
                           (safe_set_size - i - 1) * MAX_STATE_DIM * sizeof(int));
                }
                safe_set_size--;
                
                // Add new basis elements
                for (int j = 0; j < state_dim; ++j) {
                    if (basis_idx[j] == 1) continue;
                    int candidate[MAX_STATE_DIM];
                    for (int k = 0; k < state_dim; ++k) {
                        candidate[k] = basis_idx[k];
                    }
                    candidate[j] -= 1;
                    
                    if (!xInSafeSet(candidate)) {
                        addBasisIdx(candidate);
                    }
                }
                break;
            }
        }
    }
    
    void addBasisIdx(const int* x_idx) {
        if (safe_set_size >= MAX_BASIS_ELEMENTS) return;
        if (xInSafeSet(x_idx)) return;
        
        int state_dim = is_3d ? 3 : 2;
        for (int i = 0; i < state_dim; ++i) {
            safe_set_basis[safe_set_size * MAX_STATE_DIM + i] = x_idx[i];
        }
        safe_set_size++;
    }
    
    void precomputeBoundary(bool use_precomp) {
        if (!is_3d) {
            std::cout << "2D precomputation completed. Basis size: " << safe_set_size << std::endl;
            return;
        }
        
        if (!use_precomp) {
            // 3D without precomputation: Initialize with entire state space
            safe_set_basis[0] = x_numCells[0];
            safe_set_basis[1] = x_numCells[1];
            safe_set_basis[2] = x_numCells[2];
            safe_set_size = 1;
            std::cout << "3D without precomputation completed. Basis size: " << safe_set_size << std::endl;
            return;
        }
        
        int point_count = 0, added_count = 0;
        for (double h_point = H_MIN_3D; h_point <= H_MAX_3D; h_point += H_RES_3D) {
            for (double v_lead_point = V_MIN_3D; v_lead_point <= V_MAX_3D; v_lead_point += V_RES_3D) {
                double v_follow_point = outerApproxBoundary(h_point, v_lead_point);
                v_follow_point = std::min(v_follow_point, V_MAX_3D);
                
                double x_bound[3] = {h_point, v_follow_point, v_lead_point};
                int x_bound_idx[3];
                getStateIdx(x_bound, x_bound_idx);
                
                if (std::fmod(v_follow_point, V_RES_3D) > 1e-5 && x_bound_idx[1] > 1) {
                    x_bound_idx[1] -= 1;
                }
                
                int old_size = safe_set_size;
                addBasisIdx(x_bound_idx);
                if (safe_set_size > old_size) added_count++;
                point_count++;
            }
        }
        
        std::cout << "3D precomputation: Processed " << point_count << " points, added " << added_count << " new basis elements" << std::endl;
        std::cout << "3D precomputation completed. Basis size: " << safe_set_size << std::endl;
    }
    
    void computeSafeSet() {
        std::cout << "Starting " << (use_parallel ? "parallel" : "sequential") << " invariant set computation..." << std::endl;
        
        int iter = 0;
        while (true) {
            int unsafe_basis_indices[MAX_BASIS_ELEMENTS * MAX_STATE_DIM];
            int unsafe_count = 0;
            iter++;
            
            if (safe_set_size == 0) break;
            
            if (use_parallel) {
                // Parallel processing of basis elements
                std::vector<bool> unsafe_flags(safe_set_size, false);
                std::mutex unsafe_mutex;
                
                auto process_basis = [&](int start, int end) {
                    for (int i = start; i < end; ++i) {
                        int x_idx[MAX_STATE_DIM];
                        for (int j = 0; j < (is_3d ? 3 : 2); ++j) {
                            x_idx[j] = safe_set_basis[i * MAX_STATE_DIM + j];
                        }
                        
                        double x_val[MAX_STATE_DIM];
                        getPriorityStateAtIdx(x_idx, x_val);
                        
                        bool is_unsafe = false;
                        
                        if (is_3d) {
                            // 3D: Use worst-case input and disturbance
                            int u_idx[1] = {1};  // minimum input
                            int w_idx[1] = {w_numCells[0]};  // maximum disturbance
                            
                            int x_plus_idx[MAX_STATE_DIM];
                            getTransitionState(x_idx, u_idx, w_idx, x_plus_idx);
                            if (!xInSafeSet(x_plus_idx)) {
                                is_unsafe = true;
                            }
                        } else {
                            // 2D: Use MATLAB algorithm
                            double velocity = x_val[1];
                            
                            // Use maximum braking force: u_val = max(a_min, -10*velocity)
                            double u_val = std::max(U_MIN_2D, -10.0 * velocity);
                            double u_val_vec[1] = {u_val};
                            int u_idx[1];
                            getInputIdx(u_val_vec, u_idx);
                            
                            // Lead vehicle not moving: w_val = v_min
                            double w_val = V_MIN_2D;
                            double w_val_vec[1] = {w_val};
                            int w_idx[1];
                            getDisturbanceIdx(w_val_vec, w_idx);
                            
                            int x_plus_idx[MAX_STATE_DIM];
                            getTransitionState(x_idx, u_idx, w_idx, x_plus_idx);
                            
                            // Analytical safety check
                            double headway = x_val[0];
                            double stop_dist = headway - (velocity * velocity) / (2.0 * std::abs(u_val));
                            bool should_be_safe = (stop_dist >= 5.0);
                            
                            if (!xInSafeSet(x_plus_idx)) {
                                is_unsafe = true;
                                if (should_be_safe) {
                                    std::lock_guard<std::mutex> lock(unsafe_mutex);
                                    std::cout << "Warning: Point should be safe! stop_dist=" << stop_dist << std::endl;
                                }
                            }
                        }
                        
                        unsafe_flags[i] = is_unsafe;
                    }
                };
                
                // Create threads
                unsigned int num_threads = std::thread::hardware_concurrency();
                if (num_threads == 0) num_threads = 4;
                
                std::vector<std::thread> threads;
                int chunk_size = safe_set_size / num_threads;
                
                for (unsigned int t = 0; t < num_threads; ++t) {
                    int start = t * chunk_size;
                    int end = (t == num_threads - 1) ? safe_set_size : (t + 1) * chunk_size;
                    threads.emplace_back(process_basis, start, end);
                }
                
                // Wait for all threads to complete
                for (auto& thread : threads) {
                    thread.join();
                }
                
                // Collect unsafe basis indices
                for (int i = 0; i < safe_set_size; ++i) {
                    if (unsafe_flags[i]) {
                        for (int j = 0; j < (is_3d ? 3 : 2); ++j) {
                            unsafe_basis_indices[unsafe_count * MAX_STATE_DIM + j] = safe_set_basis[i * MAX_STATE_DIM + j];
                        }
                        unsafe_count++;
                    }
                }
            } else {
                // Sequential processing (original algorithm)
                for (int i = 0; i < safe_set_size; ++i) {
                    int x_idx[MAX_STATE_DIM];
                    for (int j = 0; j < (is_3d ? 3 : 2); ++j) {
                        x_idx[j] = safe_set_basis[i * MAX_STATE_DIM + j];
                    }
                    
                    double x_val[MAX_STATE_DIM];
                    getPriorityStateAtIdx(x_idx, x_val);
                    
                    bool is_unsafe = false;
                    
                    if (is_3d) {
                        // 3D: Use worst-case input and disturbance
                        int u_idx[1] = {1};  // minimum input
                        int w_idx[1] = {w_numCells[0]};  // maximum disturbance
                        
                        int x_plus_idx[MAX_STATE_DIM];
                        getTransitionState(x_idx, u_idx, w_idx, x_plus_idx);
                        if (!xInSafeSet(x_plus_idx)) {
                            is_unsafe = true;
                        }
                    } else {
                        // 2D: Use MATLAB algorithm
                        double velocity = x_val[1];
                        
                        // Use maximum braking force: u_val = max(a_min, -10*velocity)
                        double u_val = std::max(U_MIN_2D, -10.0 * velocity);
                        double u_val_vec[1] = {u_val};
                        int u_idx[1];
                        getInputIdx(u_val_vec, u_idx);
                        
                        // Lead vehicle not moving: w_val = v_min
                        double w_val = V_MIN_2D;
                        double w_val_vec[1] = {w_val};
                        int w_idx[1];
                        getDisturbanceIdx(w_val_vec, w_idx);
                        
                        int x_plus_idx[MAX_STATE_DIM];
                        getTransitionState(x_idx, u_idx, w_idx, x_plus_idx);
                        
                        // Analytical safety check
                        double headway = x_val[0];
                        double stop_dist = headway - (velocity * velocity) / (2.0 * std::abs(u_val));
                        bool should_be_safe = (stop_dist >= 5.0);
                        
                        if (!xInSafeSet(x_plus_idx)) {
                            is_unsafe = true;
                            if (should_be_safe) {
                                std::cout << "Warning: Point should be safe! stop_dist=" << stop_dist << std::endl;
                            }
                        }
                    }
                    
                    if (is_unsafe) {
                        for (int j = 0; j < (is_3d ? 3 : 2); ++j) {
                            unsafe_basis_indices[unsafe_count * MAX_STATE_DIM + j] = x_idx[j];
                        }
                        unsafe_count++;
                    }
                }
            }
            
            if (unsafe_count == 0) {
                std::cout << "Converged! All basis elements are safe." << std::endl;
                break;
            } else {
                for (int i = 0; i < unsafe_count; ++i) {
                    int basis_idx[MAX_STATE_DIM];
                    for (int j = 0; j < (is_3d ? 3 : 2); ++j) {
                        basis_idx[j] = unsafe_basis_indices[i * MAX_STATE_DIM + j];
                    }
                    removeBasisIdx(basis_idx);
                }
            }
        }
        
        std::cout << (use_parallel ? "Parallel" : "Sequential") << " safe set computation completed in " << iter << " iterations." << std::endl;
        std::cout << "Final safe set has " << safe_set_size << " basis elements." << std::endl;
        
        // Print cache statistics
        long long total_lookups = cache_hits.load() + cache_misses.load();
        if (total_lookups > 0) {
            double hit_rate = 100.0 * cache_hits.load() / total_lookups;
            std::cout << "\nDynamics Cache Statistics:" << std::endl;
            std::cout << "  Cache hits: " << cache_hits.load() << std::endl;
            std::cout << "  Cache misses: " << cache_misses.load() << std::endl;
            std::cout << "  Hit rate: " << hit_rate << "%" << std::endl;
            std::cout << "  Cache size: " << dynamics_cache.size() << " entries" << std::endl;
        }
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
    double outerApproxBoundary(double h, double v_L) const {
        if (h < H_MIN_3D) return 0.0;
        if (v_L < 0) return 0.0;
        
        double a_lead_brake = (1.0 / M_MIN) * (std::abs(T_BRAKE_MAX) / R_W_MIN + ALPHA_MAX);
        double a_ego_brake = (1.0 / M_MAX) * (std::abs(T_BRAKE_MIN) / R_W_MAX + ALPHA_MIN + BETA_MIN * V_MAX_3D + GAMMA_MIN * V_MAX_3D * V_MAX_3D);
        
        double v_squared = 2 * a_lead_brake * (h + v_L * v_L / (2 * a_ego_brake) - H_MIN_3D);
        return (v_squared < 0) ? 0.0 : std::sqrt(v_squared);
    }
    
    void getPriorityInputAtIdx(const int* idx, double* val) const {
        for (int i = 0; i < MAX_INPUT_DIM; ++i) {
            if (u_priority[i] == 1) {
                val[i] = u_range_min[i] + (idx[i] - 1) * u_res[i];
            } else {
                val[i] = u_range_max[i] - (idx[i] - 1) * u_res[i];
            }
        }
    }
    
    void getPriorityDisturbanceAtIdx(const int* idx, double* val) const {
        for (int i = 0; i < MAX_DISTURBANCE_DIM; ++i) {
            if (w_priority[i] == 1) {
                val[i] = w_range_min[i] + (idx[i] - 1) * w_res[i];
            } else {
                val[i] = w_range_max[i] - (idx[i] - 1) * w_res[i];
            }
        }
    }
    
    bool getInputIdx(const double* val, int* idx) const {
        const double threshold = 1e-5;
        
        for (int i = 0; i < MAX_INPUT_DIM; ++i) {
            if (val[i] < u_range_min[i] - threshold || val[i] > u_range_max[i] + threshold) {
                idx[0] = -1;  // Invalid input
                return false;
            }
            
            if (u_priority[i] == 1) {
                idx[i] = std::ceil((val[i] - u_range_min[i]) / u_res[i] - threshold) + 1;
            } else {
                idx[i] = std::ceil((u_range_max[i] - val[i]) / u_res[i] - threshold) + 1;
            }
        }
        return true;
    }
    
    bool getDisturbanceIdx(const double* val, int* idx) const {
        const double threshold = 1e-5;
        
        for (int i = 0; i < MAX_DISTURBANCE_DIM; ++i) {
            if (val[i] < w_range_min[i] - threshold || val[i] > w_range_max[i] + threshold) {
                idx[0] = -1;  // Invalid disturbance
                return false;
            }
            
            if (w_priority[i] == 1) {
                idx[i] = std::ceil((val[i] - w_range_min[i]) / w_res[i] - threshold) + 1;
            } else {
                idx[i] = std::ceil((w_range_max[i] - val[i]) / w_res[i] - threshold) + 1;
            }
        }
        return true;
    }
};

int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cout << "Usage: " << argv[0] << " [2d|3d] [precomp|noprecomp] [seq|parallel]" << std::endl;
        std::cout << "  precomp    - Enable precomputation (3D only)" << std::endl;
        std::cout << "  noprecomp  - Disable precomputation" << std::endl;
        std::cout << "  seq        - Sequential computation" << std::endl;
        std::cout << "  parallel   - Parallel computation" << std::endl;
        return 1;
    }
    
    std::string mode = argv[1];
    std::string precomp = argv[2];
    std::string exec_mode = argv[3];
    bool is_3d = (mode == "3d");
    bool use_precomp = (precomp == "precomp");
    bool use_parallel = (exec_mode == "parallel");
    
    if (mode != "2d" && mode != "3d") {
        std::cout << "Error: Mode must be '2d' or '3d'" << std::endl;
        return 1;
    }
    
    if (precomp != "precomp" && precomp != "noprecomp") {
        std::cout << "Error: Precomputation must be 'precomp' or 'noprecomp'" << std::endl;
        return 1;
    }
    
    if (exec_mode != "seq" && exec_mode != "parallel") {
        std::cout << "Error: Execution mode must be 'seq' or 'parallel'" << std::endl;
        return 1;
    }
    
    std::cout << "=== " << (is_3d ? "3D" : "2D") << " Safe Set Computation (CACHED) ===" << std::endl;
    std::cout << "Computing invariant safe sets using " << (use_parallel ? "parallel" : "sequential") << " monotone abstractions..." << std::endl;
    std::cout << "Optimizations: Full dynamics caching with hash table lookup" << std::endl;
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    try {
        // Create dynamics and abstraction
        MonotoneDynamics dyn(is_3d ? 3 : 2, 1, 1, is_3d ? DT_3D : DT_2D, is_3d);
        MonotoneAbstraction abs(dyn, is_3d, use_parallel);
        abs.setupRanges();
        abs.initializeSafeSet();
        abs.precomputeBoundary(use_precomp);
        abs.computeSafeSet();
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        
        std::cout << "\n" << (use_parallel ? "Parallel" : "Sequential") << " computation completed in " << duration.count() << " ms" << std::endl;
        abs.printSafeSet();
        std::cout << "\n=== Safe Set Computation Complete ===" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}

