#include <iostream>
#include <vector>
#include <unordered_set>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <thread>
#include <mutex>
#include <atomic>

/**
 * Safe Set Computation for Vehicle Following
 * Supports both sequential and parallel execution
 * 
 * Usage: ./safe_set [2d|3d] [precomp|noprecomp] [seq|parallel]
 * Example: ./safe_set 3d precomp parallel
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
    
    std::vector<double> nextState(const std::vector<double>& x, 
                                  const std::vector<double>& u, 
                                  const std::vector<double>& w) const {
        if (is_3d) {
            return solve3DDynamics(x, u, w);
        } else {
            // 2D dynamics: h(t+dt) = h(t) + v_L*dt - v(t)*dt - 0.5*a*dt^2, v(t+dt) = v(t) + a*dt
            std::vector<double> x_plus(2);
            x_plus[0] = x[0] + w[0]*dt - x[1]*dt - 0.5*u[0]*dt*dt;
            x_plus[1] = x[1] + u[0]*dt;
            return x_plus;
        }
    }
    
private:
    std::vector<double> solve3DDynamics(const std::vector<double>& x, 
                                       const std::vector<double>& u, 
                                       const std::vector<double>& w) const {
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
        
        // Solve vehicle dynamics using RK4
        auto [pos_plus, vel_plus] = solveODE(x[1], a, b, c, false);
        auto [posL_plus, velL_plus] = solveODE(x[2], a_L, b_L, c_L, true);
        
        std::vector<double> x_plus(3);
        x_plus[0] = x[0] + posL_plus - pos_plus;
        x_plus[1] = std::max(V_MIN_3D, std::min(vel_plus, V_MAX_3D));
        x_plus[2] = std::max(V_MIN_3D, std::min(velL_plus, V_MAX_3D));
        return x_plus;
    }
    
    std::pair<double, double> solveODE(double v0, double a, double b, double c, bool is_lead) const {
        const int num_steps = 1000;
        const double step_size = dt / num_steps;
        double pos = 0, vel = v0;
        
        for (int i = 0; i < num_steps; ++i) {
            auto k1 = f(pos, vel, a, b, c, is_lead);
            auto k2 = f(pos + step_size/2 * k1.first, vel + step_size/2 * k1.second, a, b, c, is_lead);
            auto k3 = f(pos + step_size/2 * k2.first, vel + step_size/2 * k2.second, a, b, c, is_lead);
            auto k4 = f(pos + step_size * k3.first, vel + step_size * k3.second, a, b, c, is_lead);
            
            pos += step_size/6 * (k1.first + 2*k2.first + 2*k3.first + k4.first);
            vel += step_size/6 * (k1.second + 2*k2.second + 2*k3.second + k4.second);
        }
        return {pos, vel};
    }
    
    std::pair<double, double> f(double /*pos*/, double vel, double a, double b, double c, bool is_lead) const {
        double dvdt = a + b * vel + c * vel * vel;
        if (vel <= 0.0 && dvdt < 0.0) dvdt = 0.0;
        if (is_lead && vel >= V_MAX_3D && dvdt > 0.0) dvdt = 0.0;
        return {vel, dvdt};
    }
};

class MonotoneAbstraction {
private:
    MonotoneDynamics dyn;
    std::vector<std::pair<double, double>> x_range, u_range, w_range;
    std::vector<double> x_res, u_res, w_res;
    std::vector<int> x_numCells, u_numCells, w_numCells;
    std::vector<int> x_priority, u_priority, w_priority;
    std::vector<std::vector<int>> safe_set_basis;
    std::unordered_map<int, std::vector<int>> transition_cache;
    bool is_3d;
    bool use_parallel;
    mutable std::mutex cache_mutex;
    
public:
    MonotoneAbstraction(const MonotoneDynamics& dynamics, bool use_3d, bool parallel = false) 
        : dyn(dynamics), is_3d(use_3d), use_parallel(parallel) {
        if (use_parallel) {
            unsigned int num_threads = std::thread::hardware_concurrency();
            if (num_threads == 0) num_threads = 8; // fallback
            std::cout << "Parallel computation initialized with " << num_threads << " threads" << std::endl;
        }
    }
    
    void setupRanges() {
        if (is_3d) {
            x_range = {{H_MIN_3D, H_MAX_3D}, {V_MIN_3D, V_MAX_3D}, {V_MIN_3D, V_MAX_3D}};
            u_range = {{T_BRAKE_MIN, T_MAX}};
            w_range = {{T_BRAKE_MAX, T_MAX}};
            x_res = {H_RES_3D, V_RES_3D, V_RES_3D};
            u_res = {T_RES_3D};
            w_res = {T_RES_3D};
            x_priority = {0, 1, 0};
            u_priority = {1};
            w_priority = {0};
        } else {
            x_range = {{H_MIN_2D, H_MAX_2D}, {V_MIN_2D, V_MAX_2D}};
            u_range = {{U_MIN_2D, U_MAX_2D}};
            w_range = {{V_MIN_2D, V_MAX_2D}};
            x_res = {H_RES_2D, V_RES_2D};
            u_res = {U_RES_2D};
            w_res = {V_RES_2D};
            x_priority = {0, 1};
            u_priority = {1};
            w_priority = {0};
        }
        
        // Compute number of cells
        x_numCells.resize(x_res.size());
        u_numCells.resize(u_res.size());
        w_numCells.resize(w_res.size());
        for (size_t i = 0; i < x_res.size(); ++i) {
            x_numCells[i] = std::ceil((x_range[i].second - x_range[i].first) / x_res[i]) + 1;
        }
        for (size_t i = 0; i < u_res.size(); ++i) {
            u_numCells[i] = std::ceil((u_range[i].second - u_range[i].first) / u_res[i]) + 1;
        }
        for (size_t i = 0; i < w_res.size(); ++i) {
            w_numCells[i] = std::ceil((w_range[i].second - w_range[i].first) / w_res[i]) + 1;
        }
    }
    
    void initializeSafeSet() {
        safe_set_basis.clear();
        if (is_3d) {
            // 3D: Start empty for precomputation
        } else {
            // 2D: Initialize with entire state space
            safe_set_basis.push_back(x_numCells);
        }
    }
    
    std::vector<double> getPriorityStateAtIdx(const std::vector<int>& idx) const {
        std::vector<double> val(idx.size());
        for (size_t i = 0; i < idx.size(); ++i) {
            if (x_priority[i] == 1) {
                val[i] = x_range[i].first + (idx[i] - 1) * x_res[i];
            } else {
                val[i] = x_range[i].second - (idx[i] - 1) * x_res[i];
            }
        }
        return val;
    }
    
    std::vector<double> getPriorityInputAtIdx(const std::vector<int>& idx) const {
        std::vector<double> val(idx.size());
        for (size_t i = 0; i < idx.size(); ++i) {
            if (u_priority[i] == 1) {
                val[i] = u_range[i].first + (idx[i] - 1) * u_res[i];
            } else {
                val[i] = u_range[i].second - (idx[i] - 1) * u_res[i];
            }
        }
        return val;
    }
    
    std::vector<double> getPriorityDisturbanceAtIdx(const std::vector<int>& idx) const {
        std::vector<double> val(idx.size());
        for (size_t i = 0; i < idx.size(); ++i) {
            if (w_priority[i] == 1) {
                val[i] = w_range[i].first + (idx[i] - 1) * w_res[i];
            } else {
                val[i] = w_range[i].second - (idx[i] - 1) * w_res[i];
            }
        }
        return val;
    }
    
    std::vector<int> getStateIdx(const std::vector<double>& val) const {
        std::vector<int> idx(val.size());
        const double threshold = 1e-5;
        
        for (size_t i = 0; i < val.size(); ++i) {
            if (val[i] < x_range[i].first - threshold || val[i] > x_range[i].second + threshold) {
                return {-1};
            }
            
            if (x_priority[i] == 1) {
                idx[i] = std::ceil((val[i] - x_range[i].first) / x_res[i] - threshold) + 1;
            } else {
                idx[i] = std::ceil((x_range[i].second - val[i]) / x_res[i] - threshold) + 1;
            }
        }
        return idx;
    }
    
    std::vector<int> getInputIdx(const std::vector<double>& val) const {
        std::vector<int> idx(val.size());
        const double threshold = 1e-5;
        
        for (size_t i = 0; i < val.size(); ++i) {
            if (val[i] < u_range[i].first - threshold || val[i] > u_range[i].second + threshold) {
                return {-1};
            }
            
            if (u_priority[i] == 1) {
                idx[i] = std::ceil((val[i] - u_range[i].first) / u_res[i] - threshold) + 1;
            } else {
                idx[i] = std::ceil((u_range[i].second - val[i]) / u_res[i] - threshold) + 1;
            }
        }
        return idx;
    }
    
    std::vector<int> getDisturbanceIdx(const std::vector<double>& val) const {
        std::vector<int> idx(val.size());
        const double threshold = 1e-5;
        
        for (size_t i = 0; i < val.size(); ++i) {
            if (val[i] < w_range[i].first - threshold || val[i] > w_range[i].second + threshold) {
                return {-1};
            }
            
            if (w_priority[i] == 1) {
                idx[i] = std::ceil((val[i] - w_range[i].first) / w_res[i] - threshold) + 1;
            } else {
                idx[i] = std::ceil((w_range[i].second - val[i]) / w_res[i] - threshold) + 1;
            }
        }
        return idx;
    }
    
    bool xInSafeSet(const std::vector<int>& x_idx) const {
        if (x_idx.size() == 1 && x_idx[0] == -1) return false;
        
        for (const auto& basis : safe_set_basis) {
            bool dominated = true;
            for (size_t i = 0; i < x_idx.size(); ++i) {
                if (x_idx[i] > basis[i]) {
                    dominated = false;
                    break;
                }
            }
            if (dominated) return true;
        }
        return false;
    }
    
    std::vector<int> getTransitionState(const std::vector<int>& x_idx, 
                                       const std::vector<int>& u_idx, 
                                       const std::vector<int>& w_idx) {
        int flat_idx = flattenIndex(x_idx);
        
        if (use_parallel) {
            std::lock_guard<std::mutex> lock(cache_mutex);
            if (transition_cache.find(flat_idx) != transition_cache.end()) {
                return transition_cache[flat_idx];
            }
        } else {
            if (transition_cache.find(flat_idx) != transition_cache.end()) {
                return transition_cache[flat_idx];
            }
        }
        
        auto x_val = getPriorityStateAtIdx(x_idx);
        auto u_val = getPriorityInputAtIdx(u_idx);
        auto w_val = getPriorityDisturbanceAtIdx(w_idx);
        
        auto x_plus_val = dyn.nextState(x_val, u_val, w_val);
        auto x_plus_idx = getStateIdx(x_plus_val);
        
        if (use_parallel) {
            std::lock_guard<std::mutex> lock(cache_mutex);
            transition_cache[flat_idx] = x_plus_idx;
        } else {
            transition_cache[flat_idx] = x_plus_idx;
        }
        return x_plus_idx;
    }
    
    void removeBasisIdx(const std::vector<int>& basis_idx) {
        auto it = std::find(safe_set_basis.begin(), safe_set_basis.end(), basis_idx);
        if (it != safe_set_basis.end()) {
            safe_set_basis.erase(it);
            
            for (size_t i = 0; i < basis_idx.size(); ++i) {
                if (basis_idx[i] == 1) continue;
                std::vector<int> candidate = basis_idx;
                candidate[i] -= 1;
                if (!xInSafeSet(candidate)) {
                    safe_set_basis.push_back(candidate);
                }
            }
        }
    }
    
    void addBasisIdx(const std::vector<int>& x_idx) {
        if (!xInSafeSet(x_idx)) {
            safe_set_basis.push_back(x_idx);
        }
    }
    
    void precomputeBoundary(bool use_precomp) {
        if (!is_3d) {
            std::cout << "2D precomputation completed. Basis size: " << safe_set_basis.size() << std::endl;
            return;
        }
        
        if (!use_precomp) {
            // 3D without precomputation: Initialize with entire state space
            safe_set_basis.push_back(x_numCells);
            std::cout << "3D without precomputation completed. Basis size: " << safe_set_basis.size() << std::endl;
            return;
        }
        
        int point_count = 0, added_count = 0;
        for (double h_point = H_MIN_3D; h_point <= H_MAX_3D; h_point += H_RES_3D) {
            for (double v_lead_point = V_MIN_3D; v_lead_point <= V_MAX_3D; v_lead_point += V_RES_3D) {
                double v_follow_point = outerApproxBoundary(h_point, v_lead_point);
                v_follow_point = std::min(v_follow_point, V_MAX_3D);
                
                std::vector<double> x_bound = {h_point, v_follow_point, v_lead_point};
                auto x_bound_idx = getStateIdx(x_bound);
                
                if (std::fmod(v_follow_point, V_RES_3D) > 1e-5 && x_bound_idx[1] > 1) {
                    x_bound_idx[1] -= 1;
                }
                
                size_t old_size = safe_set_basis.size();
                addBasisIdx(x_bound_idx);
                if (safe_set_basis.size() > old_size) added_count++;
                point_count++;
            }
        }
        
        std::cout << "3D precomputation: Processed " << point_count << " points, added " << added_count << " new basis elements" << std::endl;
        std::cout << "3D precomputation completed. Basis size: " << safe_set_basis.size() << std::endl;
    }
    
    void computeSafeSet() {
        std::cout << "Starting " << (use_parallel ? "parallel" : "sequential") << " invariant set computation..." << std::endl;
        
        int iter = 0;
        while (true) {
            std::vector<std::vector<int>> unsafe_basis_indices;
            iter++;
            
            std::cout << "Iteration: " << iter << std::endl;
            
            size_t basis_size = getBasisSize();
            if (basis_size == 0) break;
            
            if (use_parallel) {
                // Parallel processing of basis elements
                std::vector<bool> unsafe_flags(basis_size, false);
                std::mutex unsafe_mutex;
                
                auto process_basis = [&](size_t start, size_t end) {
                    for (size_t i = start; i < end; ++i) {
                        auto x_idx = getBasisIdx(i);
                        auto x_val = getPriorityStateAtIdx(x_idx);
                        
                        bool is_unsafe = false;
                        
                        if (is_3d) {
                            // 3D: Use worst-case input and disturbance (original efficient approach)
                            std::vector<int> u_idx = {1};  // minimum input
                            std::vector<int> w_idx = {static_cast<int>(w_numCells[0])};  // maximum disturbance
                            
                            auto x_plus_idx = getTransitionState(x_idx, u_idx, w_idx);
                            if (!xInSafeSet(x_plus_idx)) {
                                is_unsafe = true;
                            }
                        } else {
                            // 2D: Use MATLAB algorithm - specific input/disturbance selection
                            double velocity = x_val[1];
                            
                            // Use maximum braking force: u_val = max(a_min, -10*velocity)
                            double u_val = std::max(U_MIN_2D, -10.0 * velocity);
                            auto u_val_vec = getPriorityInputAtIdx({1}); // Get input at index 1
                            u_val_vec[0] = u_val;
                            auto u_idx = getInputIdx(u_val_vec);
                            
                            // Lead vehicle not moving: w_val = v_min
                            double w_val = V_MIN_2D;
                            auto w_val_vec = getPriorityDisturbanceAtIdx({1}); // Get disturbance at index 1
                            w_val_vec[0] = w_val;
                            auto w_idx = getDisturbanceIdx(w_val_vec);
                            
                            auto x_plus_idx = getTransitionState(x_idx, u_idx, w_idx);
                            
                            // Analytical safety check: stop_dist = headway - velocity^2/(2*abs(u_val))
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
                size_t chunk_size = basis_size / num_threads;
                
                for (unsigned int t = 0; t < num_threads; ++t) {
                    size_t start = t * chunk_size;
                    size_t end = (t == num_threads - 1) ? basis_size : (t + 1) * chunk_size;
                    threads.emplace_back(process_basis, start, end);
                }
                
                // Wait for all threads to complete
                for (auto& thread : threads) {
                    thread.join();
                }
                
                // Collect unsafe basis indices
                for (size_t i = 0; i < basis_size; ++i) {
                    if (unsafe_flags[i]) {
                        unsafe_basis_indices.push_back(getBasisIdx(i));
                    }
                }
            } else {
                // Sequential processing (original algorithm)
                for (size_t i = 0; i < getBasisSize(); ++i) {
                    auto x_idx = getBasisIdx(i);
                    auto x_val = getPriorityStateAtIdx(x_idx);
                    
                    if (is_3d) {
                        // 3D: Use worst-case input and disturbance (original efficient approach)
                        std::vector<int> u_idx = {1};  // minimum input
                        std::vector<int> w_idx = {static_cast<int>(w_numCells[0])};  // maximum disturbance
                        
                        auto x_plus_idx = getTransitionState(x_idx, u_idx, w_idx);
                        if (!xInSafeSet(x_plus_idx)) {
                            unsafe_basis_indices.push_back(x_idx);
                        }
                    } else {
                        // 2D: Use MATLAB algorithm - specific input/disturbance selection
                        double headway = x_val[0];
                        double velocity = x_val[1];
                        
                        // Use maximum braking force: u_val = max(a_min, -10*velocity)
                        double u_val = std::max(U_MIN_2D, -10.0 * velocity);
                        auto u_val_vec = getPriorityInputAtIdx({1}); // Get input at index 1
                        u_val_vec[0] = u_val;
                        auto u_idx = getInputIdx(u_val_vec);
                        
                        // Lead vehicle not moving: w_val = v_min
                        double w_val = V_MIN_2D;
                        auto w_val_vec = getPriorityDisturbanceAtIdx({1}); // Get disturbance at index 1
                        w_val_vec[0] = w_val;
                        auto w_idx = getDisturbanceIdx(w_val_vec);
                        
                        auto x_plus_idx = getTransitionState(x_idx, u_idx, w_idx);
                        
                        // Analytical safety check: stop_dist = headway - velocity^2/(2*abs(u_val))
                        double stop_dist = headway - (velocity * velocity) / (2.0 * std::abs(u_val));
                        bool should_be_safe = (stop_dist >= 5.0);
                        
                        if (!xInSafeSet(x_plus_idx)) {
                            unsafe_basis_indices.push_back(x_idx);
                            if (should_be_safe) {
                                std::cout << "Warning: Point should be safe! stop_dist=" << stop_dist << std::endl;
                            }
                        }
                    }
                }
            }
            
            if (unsafe_basis_indices.empty()) {
                std::cout << "Converged! All basis elements are safe." << std::endl;
                break;
            } else {
                std::cout << "Unsafe basis elements: " << unsafe_basis_indices.size() << std::endl;
                for (const auto& unsafe_idx : unsafe_basis_indices) {
                    removeBasisIdx(unsafe_idx);
                }
            }
        }
        
        std::cout << (use_parallel ? "Parallel" : "Sequential") << " safe set computation completed in " << iter << " iterations." << std::endl;
        std::cout << "Final safe set has " << getBasisSize() << " basis elements." << std::endl;
    }
    
    size_t getBasisSize() const { return safe_set_basis.size(); }
    std::vector<int> getBasisIdx(int basis_idx) const { return safe_set_basis[basis_idx]; }
    
    void printSafeSet(int max_elements = 10) const {
        std::cout << "\nSafe Set Basis Elements:" << std::endl;
        for (int i = 0; i < std::min(max_elements, static_cast<int>(getBasisSize())); ++i) {
            auto basis = getBasisIdx(i);
            auto state = getPriorityStateAtIdx(basis);
            std::cout << "Basis " << i << ": [";
            for (size_t j = 0; j < state.size(); ++j) {
                std::cout << state[j];
                if (j < state.size() - 1) std::cout << ", ";
            }
            std::cout << "]" << std::endl;
        }
    }
    
private:
    int flattenIndex(const std::vector<int>& idx) const {
        int result = 0, multiplier = 1;
        for (size_t i = 0; i < idx.size(); ++i) {
            result += (idx[i] - 1) * multiplier;
            multiplier *= x_numCells[i];
        }
        return result;
    }
    
    double outerApproxBoundary(double h, double v_L) const {
        if (h < H_MIN_3D) return 0.0;
        if (v_L < 0) return 0.0;
        
        double a_lead_brake = (1.0 / M_MIN) * (std::abs(T_BRAKE_MAX) / R_W_MIN + ALPHA_MAX);
        double a_ego_brake = (1.0 / M_MAX) * (std::abs(T_BRAKE_MIN) / R_W_MAX + ALPHA_MIN + BETA_MIN * V_MAX_3D + GAMMA_MIN * V_MAX_3D * V_MAX_3D);
        
        double v_squared = 2 * a_lead_brake * (h + v_L * v_L / (2 * a_ego_brake) - H_MIN_3D);
        return (v_squared < 0) ? 0.0 : std::sqrt(v_squared);
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
    
    std::cout << "=== " << (is_3d ? "3D" : "2D") << " Safe Set Computation for Vehicle Following ===" << std::endl;
    std::cout << "Computing invariant safe sets using " << (use_parallel ? "parallel" : "sequential") << " monotone abstractions..." << std::endl;
    
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
