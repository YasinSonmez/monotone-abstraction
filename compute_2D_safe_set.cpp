#include <iostream>
#include <vector>
#include <unordered_set>
#include <cmath>
#include <algorithm>
#include <chrono>

/**
 * C++ Implementation of 2D Safe Set Computation for Vehicle Following
 * 
 * This is a clean, fast, object-oriented implementation of the MATLAB
 * compute_2D_safe_set.m example, computing invariant safe sets for
 * adaptive cruise control using monotone abstractions.
 */

// Forward declarations
class MonotoneDynamics;
class MonotoneAbstraction;

/**
 * MonotoneDynamics: Represents a monotone dynamical system
 * 
 * State: [h, v] where h = headway, v = ego velocity
 * Input: a (ego acceleration)
 * Disturbance: v_L (lead velocity)
 * 
 * Dynamics: h(t+dt) = h(t) + v_L*dt - v(t)*dt - 0.5*a*dt^2
 *           v(t+dt) = v(t) + a*dt
 */
class MonotoneDynamics {
private:
    int n_x, n_u, n_w;
    double dt;
    
public:
    MonotoneDynamics(int state_dim, int input_dim, int disturbance_dim, double timestep = 0.1)
        : n_x(state_dim), n_u(input_dim), n_w(disturbance_dim), dt(timestep) {}
    
    // Vehicle following dynamics with exact discretization
    std::vector<double> nextState(const std::vector<double>& x, 
                                  const std::vector<double>& u, 
                                  const std::vector<double>& w) const {
        if (x.size() != static_cast<size_t>(n_x) || u.size() != static_cast<size_t>(n_u) || w.size() != static_cast<size_t>(n_w)) {
            throw std::invalid_argument("Invalid dimensions in nextState");
        }
        
        std::vector<double> x_plus(n_x);
        // h(t+dt) = h(t) + v_L*dt - v(t)*dt - 0.5*a*dt^2
        x_plus[0] = x[0] + w[0]*dt - x[1]*dt - 0.5*u[0]*dt*dt;
        // v(t+dt) = v(t) + a*dt
        x_plus[1] = x[1] + u[0]*dt;
        
        return x_plus;
    }
    
    int getStateDim() const { return n_x; }
    int getInputDim() const { return n_u; }
    int getDisturbanceDim() const { return n_w; }
};

/**
 * MonotoneAbstraction: Grid-based abstraction for monotone systems
 * 
 * Implements the core algorithm for computing invariant safe sets
 * using boundary basis representation and priority directions.
 */
class MonotoneAbstraction {
private:
    MonotoneDynamics dyn;
    
    // Ranges and resolutions
    std::vector<std::pair<double, double>> x_range, u_range, w_range;
    std::vector<double> x_res, u_res, w_res;
    std::vector<int> x_numCells, u_numCells, w_numCells;
    
    // Priority directions (1 = larger values more unsafe, 0 = smaller values more unsafe)
    std::vector<int> x_priority, u_priority, w_priority;
    
    // Safe set representation using boundary basis
    std::vector<std::vector<int>> safe_set_basis;
    
    // Cached transitions for performance
    std::unordered_map<int, std::vector<int>> transition_cache;
    
    // Helper functions
    int flattenIndex(const std::vector<int>& idx) const {
        int result = 0;
        int multiplier = 1;
        for (size_t i = 0; i < idx.size(); ++i) {
            result += (idx[i] - 1) * multiplier;
            multiplier *= x_numCells[i];
        }
        return result;
    }
    
    std::vector<int> unflattenIndex(int flat_idx) const {
        std::vector<int> idx(x_numCells.size());
        for (size_t i = 0; i < x_numCells.size(); ++i) {
            idx[i] = (flat_idx % x_numCells[i]) + 1;
            flat_idx /= x_numCells[i];
        }
        return idx;
    }
    
    std::vector<double> getPriorityValue(const std::vector<int>& idx, 
                                         const std::vector<std::pair<double, double>>& range,
                                         const std::vector<double>& res,
                                         const std::vector<int>& priority) const {
        std::vector<double> val(idx.size());
        for (size_t i = 0; i < idx.size(); ++i) {
            if (priority[i] == 1) {
                // Larger values have priority
                val[i] = range[i].first + (idx[i] - 1) * res[i];
            } else {
                // Smaller values have priority
                val[i] = range[i].second - (idx[i] - 1) * res[i];
            }
        }
        return val;
    }
    
    std::vector<int> getIndex(const std::vector<double>& val,
                             const std::vector<std::pair<double, double>>& range,
                             const std::vector<double>& res,
                             const std::vector<int>& priority) const {
        std::vector<int> idx(val.size());
        const double threshold = 1e-5;
        
        for (size_t i = 0; i < val.size(); ++i) {
            if (val[i] < range[i].first - threshold || val[i] > range[i].second + threshold) {
                return {-1}; // Out of bounds
            }
            
            if (priority[i] == 1) {
                idx[i] = std::ceil((val[i] - range[i].first) / res[i] - threshold) + 1;
            } else {
                idx[i] = std::ceil((range[i].second - val[i]) / res[i] - threshold) + 1;
            }
        }
        return idx;
    }
    
public:
    MonotoneAbstraction(const MonotoneDynamics& dynamics) : dyn(dynamics) {}
    
    void setXRange(const std::vector<std::pair<double, double>>& range) {
        x_range = range;
    }
    
    void setURange(const std::vector<std::pair<double, double>>& range) {
        u_range = range;
    }
    
    void setWRange(const std::vector<std::pair<double, double>>& range) {
        w_range = range;
    }
    
    void setXRes(const std::vector<double>& res) {
        x_res = res;
        x_numCells.resize(res.size());
        for (size_t i = 0; i < res.size(); ++i) {
            x_numCells[i] = std::ceil((x_range[i].second - x_range[i].first) / res[i]) + 1;
        }
    }
    
    void setURes(const std::vector<double>& res) {
        u_res = res;
        u_numCells.resize(res.size());
        for (size_t i = 0; i < res.size(); ++i) {
            u_numCells[i] = std::ceil((u_range[i].second - u_range[i].first) / res[i]) + 1;
        }
    }
    
    void setWRes(const std::vector<double>& res) {
        w_res = res;
        w_numCells.resize(res.size());
        for (size_t i = 0; i < res.size(); ++i) {
            w_numCells[i] = std::ceil((w_range[i].second - w_range[i].first) / res[i]) + 1;
        }
    }
    
    void setXPriority(const std::vector<int>& priority) {
        x_priority = priority;
    }
    
    void setUPriority(const std::vector<int>& priority) {
        u_priority = priority;
    }
    
    void setWPriority(const std::vector<int>& priority) {
        w_priority = priority;
    }
    
    void initializeSafeSet() {
        safe_set_basis.clear();
        safe_set_basis.push_back(x_numCells); // Start with entire state space
    }
    
    std::vector<double> getPriorityStateAtIdx(const std::vector<int>& idx) const {
        return getPriorityValue(idx, x_range, x_res, x_priority);
    }
    
    std::vector<double> getPriorityInputAtIdx(const std::vector<int>& idx) const {
        return getPriorityValue(idx, u_range, u_res, u_priority);
    }
    
    std::vector<double> getPriorityDisturbanceAtIdx(const std::vector<int>& idx) const {
        return getPriorityValue(idx, w_range, w_res, w_priority);
    }
    
    std::vector<int> getStateIdx(const std::vector<double>& val) const {
        return getIndex(val, x_range, x_res, x_priority);
    }
    
    std::vector<int> getInputIdx(const std::vector<double>& val) const {
        return getIndex(val, u_range, u_res, u_priority);
    }
    
    std::vector<int> getDisturbanceIdx(const std::vector<double>& val) const {
        return getIndex(val, w_range, w_res, w_priority);
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
        
        // Check cache first
        if (transition_cache.find(flat_idx) != transition_cache.end()) {
            return transition_cache[flat_idx];
        }
        
        // Compute transition
        auto x_val = getPriorityStateAtIdx(x_idx);
        auto u_val = getPriorityInputAtIdx(u_idx);
        auto w_val = getPriorityDisturbanceAtIdx(w_idx);
        
        auto x_plus_val = dyn.nextState(x_val, u_val, w_val);
        auto x_plus_idx = getStateIdx(x_plus_val);
        
        // Cache the result
        transition_cache[flat_idx] = x_plus_idx;
        return x_plus_idx;
    }
    
    void removeBasisIdx(const std::vector<int>& basis_idx) {
        // Find and remove the basis element
        auto it = std::find(safe_set_basis.begin(), safe_set_basis.end(), basis_idx);
        if (it != safe_set_basis.end()) {
            safe_set_basis.erase(it);
            
            // Add new elements to basis if needed
            for (size_t i = 0; i < basis_idx.size(); ++i) {
                if (basis_idx[i] == 1) continue; // Out of bounds
                
                std::vector<int> candidate = basis_idx;
                candidate[i] -= 1;
                
                if (!xInSafeSet(candidate)) {
                    safe_set_basis.push_back(candidate);
                }
            }
        }
    }
    
    size_t getBasisSize() const {
        return safe_set_basis.size();
    }
    
    std::vector<int> getBasisIdx(int basis_idx) const {
        return safe_set_basis[basis_idx];
    }
    
    // Main invariant set algorithm
    void computeSafeSet() {
        std::cout << "Starting invariant set computation..." << std::endl;
        
        int iter = 0;
        while (true) {
            std::vector<std::vector<int>> unsafe_basis_indices;
            iter++;
            std::cout << "Iteration " << iter << std::endl;
            
            for (size_t i = 0; i < getBasisSize(); ++i) {
                // Test points on the boundary of the safe set
                auto x_idx = getBasisIdx(i);
                auto x_val = getPriorityStateAtIdx(x_idx);
                
                // Use max possible braking force
                double u_val = std::max(-3.0, -10.0 * x_val[1]);
                auto u_idx = getInputIdx({u_val});
                
                // Assume front vehicle not moving (worst case)
                auto w_idx = getDisturbanceIdx({0.0});
                
                auto x_plus_idx = getTransitionState(x_idx, u_idx, w_idx);
                
                // Check if the point should be safe analytically
                double headway = x_val[0];
                double velocity = x_val[1];
                double stop_dist = headway - velocity * velocity / (2.0 * std::abs(u_val));
                
                // Check if the point is safe
                if (!xInSafeSet(x_plus_idx)) {
                    unsafe_basis_indices.push_back(x_idx);
                    if (stop_dist >= 5.0) {
                        // Should be safe!
                    }
                }
            }
            
            if (unsafe_basis_indices.empty()) {
                std::cout << "Converged! All basis elements are safe." << std::endl;
                break;
            } else {
                for (const auto& unsafe_idx : unsafe_basis_indices) {
                    removeBasisIdx(unsafe_idx);
                }
            }
        }
        
        std::cout << "Safe set computation completed in " << iter << " iterations." << std::endl;
        std::cout << "Final safe set has " << getBasisSize() << " basis elements." << std::endl;
    }
    
    // Print safe set for visualization
    void printSafeSet() const {
        std::cout << "\nSafe Set Basis Elements:" << std::endl;
        for (size_t i = 0; i < getBasisSize(); ++i) {
            auto basis = getBasisIdx(i);
            auto state = getPriorityStateAtIdx(basis);
            std::cout << "Basis " << i << ": [" << state[0] << ", " << state[1] << "]" << std::endl;
        }
    }
};

/**
 * Main function demonstrating 2D safe set computation
 */
int main() {
    std::cout << "=== 2D Safe Set Computation for Vehicle Following ===" << std::endl;
    std::cout << "Computing invariant safe sets using monotone abstractions..." << std::endl;
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Create vehicle following dynamics
    MonotoneDynamics veh_follow_dyn(2, 1, 1, 0.1);
    
    // Create abstraction
    MonotoneAbstraction veh_follow_abs(veh_follow_dyn);
    
    // Motion constraints
    const double h_min = 5.0;   // m
    const double h_max = 80.0;  // m
    const double v_min = 0.0;   // m/s
    const double v_max = 20.0; // m/s
    const double a_min = -3.0;  // m/s^2
    const double a_max = 3.0;  // m/s^2
    
    // Set ranges
    veh_follow_abs.setXRange({{h_min, h_max}, {v_min, v_max}});
    veh_follow_abs.setURange({{a_min, a_max}});
    veh_follow_abs.setWRange({{v_min, v_max}});
    
    // Set resolutions
    const double h_res = 0.2;
    const double v_res = 0.1;
    const double a_res = 0.1;
    
    veh_follow_abs.setXRes({h_res, v_res});
    veh_follow_abs.setURes({a_res});
    veh_follow_abs.setWRes({v_res});
    
    // Set priorities
    veh_follow_abs.setXPriority({0, 1}); // smaller headway more unsafe, larger velocity more unsafe
    veh_follow_abs.setUPriority({1});    // larger acceleration more unsafe
    veh_follow_abs.setWPriority({0});    // smaller velocity more unsafe
    
    // Initialize and compute safe set
    veh_follow_abs.initializeSafeSet();
    veh_follow_abs.computeSafeSet();
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    
    std::cout << "\nComputation completed in " << duration.count() << " ms" << std::endl;
    
    // Print results
    veh_follow_abs.printSafeSet();
    
    std::cout << "\n=== Safe Set Computation Complete ===" << std::endl;
    
    return 0;
}
