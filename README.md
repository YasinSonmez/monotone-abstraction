# Monotone Abstractions for Vehicle Safety Sets

This repository contains implementations for computing, verifying, and visualizing safety (invariant) sets and maximal safety controllers for monotone dynamical systems arising in two traffic scenarios:

- Adaptive cruise control (ACC) vehicle-following
- Unprotected left turn with oncoming traffic (two variants: ego-first, oncoming-first)

It provides:
- **C++ GPU-accelerated implementation** (OpenCL) for fast safe set computation
- MATLAB library for defining monotone dynamics and creating grid-based monotone abstractions (`library/`)
- Scenario-specific constructors and constants
- Scripts to compute 2D/3D safe sets (from scratch), optionally using precomputation to speed up convergence
- Plotting utilities and example precomputed data


## Quick Start (Recommended)

### Run C++ GPU-Accelerated Version

**Requirements:**
- C++ compiler with C++17 support (g++, clang++)
- OpenCL support (available on macOS, Linux, Windows)
- macOS: Built-in OpenCL framework
- Linux: Install `ocl-icd-opencl-dev` or `opencl-headers`

**Run the clean, optimized version in 3D mode:**

```bash
./build_and_run.sh safe_set_gpu_hybrid_clean 3d
```

**Other options:**

```bash
# Run in 2D mode
./build_and_run.sh safe_set_gpu_hybrid_clean 2d

# Run the original version with detailed output
./build_and_run.sh safe_set_gpu_hybrid 3d
```

**What it does:**
- Automatically compiles the C++ source with optimizations
- Runs 10 iterations of the safe set computation
- Reports detailed timing breakdown (GPU upload, kernel execution, readback, CPU update)
- Uses GPU acceleration for maximum performance
- Caches transition computations for faster repeated runs

**Performance:** The GPU implementation is significantly faster than MATLAB (minutes vs hours for 3D problems).

## MATLAB Requirements

- MATLAB R2020a or newer recommended (tested with R2020a/R2021+)
- Optimization Toolbox not required; uses base MATLAB plus `ode45`
- Enough RAM/CPU time for grid algorithms (3D sets can run for hours)

Optional (for figures): the plotting helpers under `plot_scripts/` and bundled `freezeColors` utilities are included.


## Repository layout

### C++ GPU Implementation (Root Directory)

- **`safe_set_gpu_hybrid_clean.cpp`**: Clean, professional 500-line GPU-accelerated implementation (recommended)
- **`safe_set_gpu_hybrid.cpp`**: Original GPU implementation with detailed output
- **`build_and_run.sh`**: Universal build and run script for C++ files
- **`precompute_transitions.cl`**: OpenCL kernel for GPU transition computation
- **`check_basis_safety.cl`**: OpenCL kernel for GPU safety checking
- **`transition_cache_*.bin`**: Binary cache files for precomputed transitions for 3D

### MATLAB Library

- `library/`
  - `@monotone_dyn/monotone_dyn.m`: lightweight wrapper representing a monotone dynamical system; dispatches to a user-provided dynamics function and enforces dimensional checks.
  - `@monotone_abstraction/monotone_abstraction.m`: grid abstraction and invariant-set algorithms, including boundary-basis representation, transitions, caching, and membership checks.
  - Additional scenario-specific classes under `library/` are provided for completeness but are not direct entry points.

- `acc-example/`
  - `2D_safe_set/compute_2D_safe_set.m`: minimal 2D toy version (headway and ego speed, disturbance is lead speed) demonstrating invariant-set computation.
  - `3D_safe_set/`: “hard-collision” ACC scenario with 3D state `[h; v; v_L]` and disturbance `a_L` (lead acceleration)
    - `veh_follow_const.m`: modeling and grid constants (ranges, resolution, sampling, torque bounds)
    - `veh_follow_expr.m`: dynamics using a parametric drag model, integrated via `ode45`
    - `get_3D_monotone_abs.m`: constructs the `monotone_abstraction` object and sets ranges, resolutions, and priority directions
    - `compute_3D_safe_set.m`: main script to compute the invariant (safe) set boundary basis; includes optional precomputation using `helper_functions/outer_approx_boundary.m`
  - `3D_safe_set_soft_collisions/`: ACC variant allowing “soft” collisions bounded by maximum allowable relative velocity
    - `veh_follow_const_v2.m`, `veh_follow_expr_v2.m`, `get_3D_monotone_abs_v2.m`, `compute_3D_safe_set_v2.m`
  - `plot_scripts/`: plotting helpers for ACC safe sets
  - `helper_functions/`: utilities such as `outer_approx_boundary.m` and `unsafe_collision_check.m`
  - `data/`: example precomputed `.mat` and `.fig` results

- `turn-example/`
  - `3D_safe_set_ego_first/`: unprotected left turn where the ego vehicle goes first
    - `veh_turn_const_v3.m`, `veh_turn_expr_v3.m`, `get_3D_monotone_abs_v3.m`
    - `compute_3D_safe_set_v3.m`: computes the ego-first safe set
    - `compute_controller_ego_first.m`: constructs a maximal safety controller from stored safe sets
  - `3D_safe_set_oncoming_first/`: oncoming vehicle has priority
    - `veh_turn_const_v4.m`, `veh_turn_expr_v4.m`, `get_3D_monotone_abs_v4.m`
    - `compute_3D_safe_set_v4.m`: computes the oncoming-first safe set
    - `compute_controller_oncoming_first.m`: constructs a maximal safety controller from stored safe sets
  - `plot_scripts/`: plotting helpers for both turn scenarios
  - `data/` and `data_old/`: example precomputed results


## Concepts in brief

- Monotone abstraction: The continuous dynamics are discretized on a grid with per-dimension priority directions specifying which end of each axis is “more unsafe.” This allows representing down-closed safe sets by their boundary basis.
- Safe set (robust control invariant set): The set of states from which, under worst-case disturbances, there exists a control that keeps the next state inside the set. Here we compute the maximal such set under specified priority discretization and constraints.
- Basis iteration: The algorithm initializes with the whole state space as a single basis element, then iteratively removes boundary elements whose worst-case successor exits the set, adding their immediate predecessors until convergence.


## How to run

### C++ GPU Implementation (Recommended for Performance)

The GPU-accelerated C++ implementation provides the fastest computation times. It includes:
- **`safe_set_gpu_hybrid_clean.cpp`**: Clean, well-documented 500-line implementation (recommended)
- **`safe_set_gpu_hybrid.cpp`**: Original version with additional verbose output

#### Running with the build script:

```bash
# Recommended: Run the clean version in 3D mode
./build_and_run.sh safe_set_gpu_hybrid_clean 3d

# Run in 2D mode
./build_and_run.sh safe_set_gpu_hybrid_clean 2d

# Run the original version
./build_and_run.sh safe_set_gpu_hybrid 3d
```

The script automatically:
1. Cleans previous builds
2. Compiles with `-O3` optimization and OpenCL support
3. Runs the program with specified arguments
4. Reports performance metrics per iteration

#### Manual compilation:

```bash
# macOS
g++ -std=c++17 -O3 -o safe_set_gpu_hybrid_clean safe_set_gpu_hybrid_clean.cpp -framework OpenCL

# Linux
g++ -std=c++17 -O3 -o safe_set_gpu_hybrid_clean safe_set_gpu_hybrid_clean.cpp -lOpenCL

# Run
./safe_set_gpu_hybrid_clean 3d
```

#### Output:

The program reports:
- Setup time (OpenCL initialization, transition table precomputation)
- Per-run statistics (10 runs by default):
  - Iterations to convergence
  - Number of basis elements
  - Time breakdown: Upload, Kernel execution, Readback, CPU update
- Average statistics across all runs

#### Features:

- **Transition caching**: First run precomputes transitions, subsequent runs load from cache
- **GPU acceleration**: OpenCL kernels for safety checks
- **Neighbor deduplication**: Efficient on-the-fly deduplication using flat indices
- **Clean code**: Professional, modular, well-commented implementation

---

### MATLAB Implementation

General notes:
- Use MATLAB's Current Folder to the repository root or provide absolute paths from below.
- Scripts that compute 3D safe sets can be long-running (from tens of minutes to several hours depending on grid resolution and machine). Consider starting with the 2D example or using precomputation sections (uncommented) to accelerate convergence.
- **For faster computation, use the C++ GPU implementation above.**

### 0) Open MATLAB in the repo

- Set Current Folder to `/Users/cusgadmin/Projects/monotone-abstraction`
- Add subfolders to the path (once per session):

```matlab
addpath(genpath('/Users/cusgadmin/Projects/monotone-abstraction'));
```


### 1) Quick demo: 2D ACC safe set

Run:

```matlab
cd('/Users/cusgadmin/Projects/monotone-abstraction/acc-example/2D_safe_set');
compute_2D_safe_set;
```

This prints iterations and converges quickly, demonstrating the mechanics on a small grid.


### 2) ACC: 3D safe set (hard collisions)

Entry points:
- Build abstraction and compute safe set: `acc-example/3D_safe_set/compute_3D_safe_set.m`
- Constructors/constants: `get_3D_monotone_abs.m`, `veh_follow_const.m`, `veh_follow_expr.m`

Run from MATLAB:

```matlab
cd('/Users/cusgadmin/Projects/monotone-abstraction/acc-example/3D_safe_set');
[con, veh_follow_abs] = get_3D_monotone_abs();
veh_follow_abs.initialize_safe_set();
% Optional precompute boundary candidates (uncomment in script to enable)
compute_3D_safe_set;  % may take hours on default resolution
```

After convergence, `veh_follow_abs.safe_set_basis` holds the boundary. To plot, use or adapt `acc-example/plot_scripts/plot_3D_safe_set.m`. That script assumes a variable `veh_follow_abs` in the workspace and uses `veh_follow_const_v2()`; either change to `veh_follow_const()` or switch to the v2 scenario (soft collisions) accordingly.

Minimal plotting example for the hard-collision set present in workspace as `veh_follow_abs`:

```matlab
cd('/Users/cusgadmin/Projects/monotone-abstraction/acc-example/plot_scripts');
plot_3D_safe_set;  % ensure it references the correct constants function
```

Tip: To accelerate, uncomment the precomputation loop in `compute_3D_safe_set.m` which seeds the basis using the analytic over-approximation in `helper_functions/outer_approx_boundary.m`.


### 3) ACC: 3D safe set (soft collisions allowed)

Entry points:
- Build abstraction and compute set: `acc-example/3D_safe_set_soft_collisions/compute_3D_safe_set_v2.m`
- Constructors/constants: `get_3D_monotone_abs_v2.m`, `veh_follow_const_v2.m`, `veh_follow_expr_v2.m`

Run:

```matlab
cd('/Users/cusgadmin/Projects/monotone-abstraction/acc-example/3D_safe_set_soft_collisions');
[con, veh_follow_abs] = get_3D_monotone_abs_v2();
veh_follow_abs.initialize_safe_set();
compute_3D_safe_set_v2;  % long-running; uses optional precompute block if you uncomment it
```

Plotting:

```matlab
cd('/Users/cusgadmin/Projects/monotone-abstraction/acc-example/plot_scripts');
plot_3D_safe_set;  % this script was authored for v2 constants; adjust if needed
```

You can also open example figures in `acc-example/data/` to compare with expected results.


### 4) Unprotected left turn: ego-first scenario

Entry points:
- Build abstraction and compute set: `turn-example/3D_safe_set_ego_first/compute_3D_safe_set_v3.m`
- Constructors/constants: `get_3D_monotone_abs_v3.m`, `veh_turn_const_v3.m`, `veh_turn_expr_v3.m`

Run:

```matlab
cd('/Users/cusgadmin/Projects/monotone-abstraction/turn-example/3D_safe_set_ego_first');
[con, veh_turn_abs] = get_3D_monotone_abs_v3();
veh_turn_abs.initialize_safe_set();
compute_3D_safe_set_v3;  % ~75 iterations reported in script comments
save('veh_turn_abs_ego_first_latest.mat','veh_turn_abs');
```

Plotting:

```matlab
cd('/Users/cusgadmin/Projects/monotone-abstraction/turn-example/plot_scripts');
veh_turn_abs = load('/Users/cusgadmin/Projects/monotone-abstraction/turn-example/data/veh_turn_abs_ego_first_5_7_2021.mat','veh_turn_abs').veh_turn_abs;  % or use the one from your workspace
plot_3D_safe_set_ego_first;
```

Maximal safety controller:

```matlab
cd('/Users/cusgadmin/Projects/monotone-abstraction/turn-example/3D_safe_set_ego_first');
% expects a saved abstraction with controllable set stored per script header
% adjust filename to the artifact you saved above
save('veh_turn_abs_ego_first_5_7_2021.mat','veh_turn_abs');
compute_controller_ego_first;  % produces ego_first_controller array in workspace
```


### 5) Unprotected left turn: oncoming-first scenario

Entry points:
- Build abstraction and compute set: `turn-example/3D_safe_set_oncoming_first/compute_3D_safe_set_v4.m`
- Constructors/constants: `get_3D_monotone_abs_v4.m`, `veh_turn_const_v4.m`, `veh_turn_expr_v4.m`

Run:

```matlab
cd('/Users/cusgadmin/Projects/monotone-abstraction/turn-example/3D_safe_set_oncoming_first');
[con, veh_turn_abs] = get_3D_monotone_abs_v4();
veh_turn_abs.initialize_safe_set();
compute_3D_safe_set_v4;  % ~77 iterations reported in script comments
save('veh_turn_abs_oncoming_first_latest.mat','veh_turn_abs');
```

Plotting:

```matlab
cd('/Users/cusgadmin/Projects/monotone-abstraction/turn-example/plot_scripts');
veh_turn_abs = load('/Users/cusgadmin/Projects/monotone-abstraction/turn-example/data/veh_turn_abs_oncoming_first_5_7_2021.mat','veh_turn_abs').veh_turn_abs;  % or use the one from your run
plot_3D_safe_set_oncoming_first;
```

Maximal safety controller:

```matlab
cd('/Users/cusgadmin/Projects/monotone-abstraction/turn-example/3D_safe_set_oncoming_first');
save('veh_turn_abs_oncoming_first_5_7_2021.mat','veh_turn_abs');
compute_controller_oncoming_first;  % produces oncoming_first_controller array
```


## Performance tips

### MATLAB Performance:
- Grid resolution strongly impacts runtime and memory. Reduce `*_res` values (coarser grids) in constants files to speed up, at the cost of fidelity.
- Enable the precomputation blocks in `compute_*` scripts to seed the safe-set boundary using analytic outer approximations. This can dramatically cut iterations.
- The transition cache in `monotone_abstraction` speeds up repeated successor queries; avoid changing ranges/resolutions mid-run.


## Reproducing paper-quality figures

- ACC variants: compare your computed sets against `.fig`/`.mat` in `acc-example/data/` and `acc-example/data_old/`.
- Turn variants: use `turn-example/data/` artifacts as references and the provided plotting scripts to match viewing angles and axes.


## Troubleshooting

### C++ Implementation:
- **Missing OpenCL**: If compilation fails with OpenCL errors:
  - macOS: OpenCL is built-in, ensure Xcode Command Line Tools are installed
  - Linux: Install `ocl-icd-opencl-dev` (Debian/Ubuntu) or `ocl-icd` (Arch)
  - Windows: Install GPU vendor's SDK (NVIDIA CUDA, AMD APP SDK, or Intel OpenCL)
- **Build script permission denied**: Run `chmod +x build_and_run.sh`
- **Kernel files not found**: Ensure `precompute_transitions.cl` and `check_basis_safety.cl` are in the same directory
- **Cache files**: Delete `transition_cache_*.bin` files if you modify grid parameters

### MATLAB Implementation:
- Dimension errors typically arise if a dynamics function signature or dimensions do not match the `monotone_dyn` settings. Ensure your `get_*_monotone_abs*.m` constructor and dynamics function agree on `(n_x, n_u, n_w)`.
- If plotting scripts reference `veh_follow_const_v2()` but you ran the v1 ACC scenario, switch the call to `veh_follow_const()` or run the v2 scenario.
- Long runs show iteration counters and progress; MATLAB may clear the command window (`clc`) inside loops by design in these scripts.


## License

If not otherwise specified, this repository is provided under an academic/research-friendly license. See individual folders for any third-party licenses (e.g., `freezeColors`).

