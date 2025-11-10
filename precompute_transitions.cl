/**
 * OpenCL Kernel for Pre-computing Worst-Case State Transitions
 * Each work item computes the next state for one grid cell in parallel
 */

// Vehicle dynamics parameters (3D)
#define T_MAX 1200.0
#define T_BRAKE_MIN -1800.0
#define T_BRAKE_MAX -2400.0
#define M_MIN 2000.0
#define M_MAX 2500.0
#define R_W_MIN 0.30
#define R_W_MAX 0.35
#define ALPHA_MIN 300.0
#define ALPHA_MAX 350.0
#define BETA_MIN 0.10
#define BETA_MAX 0.25
#define GAMMA_MIN 0.30
#define GAMMA_MAX 0.65
#define V_MIN_3D 0.0
#define V_MAX_3D 20.0

// 2D parameters
#define U_MIN_2D -3.0
#define V_MIN_2D 0.0

// ODE solver parameters
#define NUM_STEPS 1000

/**
 * ODE right-hand side: computes derivatives for position and velocity
 */
void f(double vel, double a, double b, double c, bool is_lead, double* dpos, double* dvel) {
    double dvdt = a + b * vel + c * vel * vel;
    if (vel <= 0.0 && dvdt < 0.0) dvdt = 0.0;
    if (is_lead && vel >= V_MAX_3D && dvdt > 0.0) dvdt = 0.0;
    *dpos = vel;
    *dvel = dvdt;
}

/**
 * RK4 ODE solver for vehicle dynamics
 */
void solveODE(double v0, double a, double b, double c, bool is_lead, double dt, 
              double* pos, double* vel) {
    const double step_size = dt / NUM_STEPS;
    *pos = 0.0;
    *vel = v0;
    
    for (int i = 0; i < NUM_STEPS; ++i) {
        double k1_pos, k1_vel, k2_pos, k2_vel, k3_pos, k3_vel, k4_pos, k4_vel;
        
        f(*vel, a, b, c, is_lead, &k1_pos, &k1_vel);
        f(*vel + step_size/2 * k1_vel, a, b, c, is_lead, &k2_pos, &k2_vel);
        f(*vel + step_size/2 * k2_vel, a, b, c, is_lead, &k3_pos, &k3_vel);
        f(*vel + step_size * k3_vel, a, b, c, is_lead, &k4_pos, &k4_vel);
        
        *pos += step_size/6 * (k1_pos + 2*k2_pos + 2*k3_pos + k4_pos);
        *vel += step_size/6 * (k1_vel + 2*k2_vel + 2*k3_vel + k4_vel);
    }
}

/**
 * Solve 3D vehicle dynamics with worst-case parameters
 */
void solve3DDynamicsWorstCase(double x0, double x1, double x2, double u_worst, double w_worst, 
                              double dt, double* x_plus) {
    // Choose worst-case parameters for ego vehicle (braking)
    double R_w_ego = (u_worst > 0) ? R_W_MIN : R_W_MAX;
    double M_ego = ((u_worst / R_w_ego - ALPHA_MIN) > 0) ? M_MIN : M_MAX;
    double a = (1.0 / M_ego) * (u_worst / R_w_ego - ALPHA_MIN);
    double b = -(1.0 / M_MAX) * BETA_MIN;
    double c = -(1.0 / M_MAX) * GAMMA_MIN;
    
    // Choose worst-case parameters for lead vehicle
    double R_w_L = (w_worst > 0) ? R_W_MAX : R_W_MIN;
    double M_L = ((w_worst / R_w_L - ALPHA_MAX) > 0) ? M_MAX : M_MIN;
    double a_L = (1.0 / M_L) * (w_worst / R_w_L - ALPHA_MAX);
    double b_L = -(1.0 / M_MIN) * BETA_MAX;
    double c_L = -(1.0 / M_MIN) * GAMMA_MAX;
    
    // Solve vehicle dynamics using RK4
    double pos_plus, vel_plus, posL_plus, velL_plus;
    solveODE(x1, a, b, c, false, dt, &pos_plus, &vel_plus);
    solveODE(x2, a_L, b_L, c_L, true, dt, &posL_plus, &velL_plus);
    
    x_plus[0] = x0 + posL_plus - pos_plus;
    x_plus[1] = fmax((double)V_MIN_3D, fmin((double)vel_plus, (double)V_MAX_3D));
    x_plus[2] = fmax((double)V_MIN_3D, fmin((double)velL_plus, (double)V_MAX_3D));
}

/**
 * Convert flattened index to multi-dimensional indices
 */
void unflattenIndex(int flat_idx, int state_dim, __global const int* x_numCells, int* idx) {
    if (state_dim == 2) {
        idx[1] = (flat_idx % x_numCells[1]) + 1;
        idx[0] = (flat_idx / x_numCells[1]) + 1;
    } else {  // 3D
        idx[2] = (flat_idx % x_numCells[2]) + 1;
        int temp = flat_idx / x_numCells[2];
        idx[1] = (temp % x_numCells[1]) + 1;
        idx[0] = (temp / x_numCells[1]) + 1;
    }
}

/**
 * Convert indices to physical state values using priority directions
 */
void getPriorityStateAtIdx(int idx0, int idx1, int idx2, int state_dim,
                           __global const double* x_range_min,
                           __global const double* x_range_max,
                           __global const double* x_res,
                           __global const int* x_priority,
                           double* val) {
    int idx[3];
    idx[0] = idx0; idx[1] = idx1; idx[2] = idx2;
    for (int i = 0; i < state_dim; ++i) {
        if (x_priority[i] == 1) {
            val[i] = x_range_min[i] + (idx[i] - 1) * x_res[i];
        } else {
            val[i] = x_range_max[i] - (idx[i] - 1) * x_res[i];
        }
    }
}

/**
 * Convert physical state values to indices
 */
void getStateIdx(double val0, double val1, double val2, int state_dim,
                 __global const double* x_range_min,
                 __global const double* x_range_max,
                 __global const double* x_res,
                 __global const int* x_priority,
                 __global const int* x_numCells,
                 int* idx) {
    double val[3];
    val[0] = val0; val[1] = val1; val[2] = val2;
    for (int i = 0; i < state_dim; ++i) {
        double val_clamped = fmax((double)x_range_min[i], fmin((double)val[i], (double)x_range_max[i]));
        
        if (x_priority[i] == 1) {
            idx[i] = (int)floor((val_clamped - x_range_min[i]) / x_res[i]) + 1;
        } else {
            idx[i] = (int)floor((x_range_max[i] - val_clamped) / x_res[i]) + 1;
        }
        // Clamp to valid range
        if (idx[i] < 1) idx[i] = 1;
        if (idx[i] > x_numCells[i]) idx[i] = x_numCells[i];
    }
}

/**
 * Main kernel: Compute worst-case next state for each grid cell
 * 
 * @param flat_idx: Global work item ID (one per state)
 * @param total_states: Total number of states in grid
 * @param state_dim: Dimensionality (2 or 3)
 * @param is_3d: Whether using 3D dynamics
 * @param dt: Time step
 * @param x_range_min: Minimum values for each dimension
 * @param x_range_max: Maximum values for each dimension
 * @param x_res: Resolution for each dimension
 * @param x_numCells: Number of cells in each dimension
 * @param x_priority: Priority direction for each dimension (0=decreasing, 1=increasing)
 * @param next_state_table: Output buffer (total_states * 3 integers)
 */
__kernel void precompute_transitions(
    const int total_states,
    const int state_dim,
    const int is_3d,
    const double dt,
    __global const double* x_range_min,
    __global const double* x_range_max,
    __global const double* x_res,
    __global const int* x_numCells,
    __global const int* x_priority,
    __global int* next_state_table
) {
    int flat_idx = get_global_id(0);
    
    // Bounds check
    if (flat_idx >= total_states) return;
    
    // Unflatten index to get multi-dimensional indices
    int x_idx[3];
    unflattenIndex(flat_idx, state_dim, x_numCells, x_idx);
    
    // Convert indices to physical state values
    double x_val[3];
    getPriorityStateAtIdx(x_idx[0], x_idx[1], x_idx[2], state_dim, 
                          x_range_min, x_range_max, x_res, x_priority, x_val);
    
    // Compute worst-case next state
    double x_plus_val[3];
    if (is_3d) {
        // 3D worst-case: minimum input, maximum disturbance
        double u_worst = T_BRAKE_MIN;
        double w_worst = T_BRAKE_MAX;
        solve3DDynamicsWorstCase(x_val[0], x_val[1], x_val[2], u_worst, w_worst, dt, x_plus_val);
    } else {
        // 2D worst-case
        double velocity = x_val[1];
        double u_worst = fmax((double)U_MIN_2D, -10.0 * velocity);
        double w_worst = V_MIN_2D;
        
        x_plus_val[0] = x_val[0] + w_worst*dt - x_val[1]*dt - 0.5*u_worst*dt*dt;
        x_plus_val[1] = x_val[1] + u_worst*dt;
        x_plus_val[2] = 0.0;  // Unused in 2D
    }
    
    // Convert next state values to indices
    int x_plus_idx[3];
    getStateIdx(x_plus_val[0], x_plus_val[1], x_plus_val[2], state_dim, 
                x_range_min, x_range_max, x_res, x_priority, x_numCells, x_plus_idx);
    
    // Store result in output buffer
    for (int j = 0; j < state_dim; ++j) {
        next_state_table[flat_idx * 3 + j] = x_plus_idx[j];
    }
}
