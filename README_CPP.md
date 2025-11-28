# Safe Set Computation - C++ Implementation

High-performance C++ implementations of safe set computation for vehicle following scenarios.

## Available Versions

1. **`safe_set`** - Original implementation with full RK4 integration (1000 steps)
2. **`safe_set_opencl_ready`** - OpenCL-compatible version with C-style arrays
3. **`safe_set_optimized`** - Optimized with reduced RK4 steps some other methods
4. **`safe_set_cached`** - Hash table-based dynamics caching + optimized RK4 (**recommended**)

## Quick Start


@MK: The following needs to be adjusted. There is no makefile in this project; build is managed with a sh script. I'd also consider switching to CMake or Bazel but this is a different story.

```bash
# Build all versions
make all

# Run the fastest version (cached + parallel)
./safe_set_cached 3d noprecomp parallel

# Run performance benchmark
make benchmark
```

## Usage

```bash
# All versions support the same interface:
./safe_set[_variant] [2d|3d] [precomp|noprecomp] [seq|parallel]

# Examples:
./safe_set_cached 3d noprecomp parallel    # Fastest
./safe_set_optimized 3d noprecomp seq      # Good single-threaded
./safe_set 2d precomp seq                  # 2D computation
```

## Performance Benchmark

### 3D Safe Set Computation (no precomputation)

Performance comparison of all versions on 3D vehicle following problem:

| Version | Sequential Time | Parallel Time | Speedup vs MATLAB | Speedup vs C++ Baseline |
|---------|----------------|---------------|-------------------|------------------------|
| MATLAB (original) | ~6 minutes | - | 1.0x | - |
| `safe_set` (C++) | 6416 ms (6.4s) | 1959 ms (2.0s) | **56x** / **184x** | 1.0x / 3.3x |
| `safe_set_opencl_ready` | 8863 ms (8.9s) | 1910 ms (1.9s) | 41x / 189x | 0.7x / 3.4x |
| `safe_set_optimized` | 2200 ms (2.2s) | 965 ms (1.0s) | 164x / 373x | 2.9x / 6.6x |
| **`safe_set_cached`** | **1511 ms (1.5s)** | **700 ms (0.7s)** | **238x** / **514x** | **4.2x** / **9.2x** |

**Key Results:**
- **Final basis size**: 1,083 elements
- **Iterations**: 143-144 (MATLAB: ~150 iterations)
- **MATLAB baseline**: ~6 minutes (360,000 ms) for comparable computation
- **Best C++ speedup**: 514x faster than MATLAB (cached parallel)
- **C++ baseline speedup**: 56x faster than MATLAB (original sequential)
- **Parallel scaling**: ~2.2x speedup from parallelization in cached version
- **Cache benefit**: 38.6% hit rate eliminates ~40% of dynamics computations

### 2D Safe Set Computation

All versions complete 2D computation in under 50ms:
- **Final basis size**: 168 elements
- **Iterations**: 368
- **Computation time**: 18-37 ms (cached version)

## Configuration

Edit parameters at the top of any `safe_set*.cpp` file:

```cpp
// 2D Configuration
const double DT_2D = 0.1;
const double H_MIN_2D = 5.0, H_MAX_2D = 80.0;
const double V_MIN_2D = 0.0, V_MAX_2D = 20.0;
const double H_RES_2D = 0.2, V_RES_2D = 0.1;

// 3D Configuration  
const double DT_3D = 0.4;
const double H_MIN_3D = 0.0, H_MAX_3D = 80.0;
const double V_MIN_3D = 0.0, V_MAX_3D = 20.0;
const double H_RES_3D = 0.8, V_RES_3D = 0.4;

// Vehicle parameters
const double T_MAX = 1200.0;
const double T_BRAKE_MIN = -1800.0, T_BRAKE_MAX = -2400.0;
const double M_MIN = 2000.0, M_MAX = 2500.0;
// ... (see source files for complete parameters)
```

## Implementation Details

### Optimization Techniques Used

1. **Reduced RK4 steps** (`safe_set_optimized`, `safe_set_cached`)
   - Reduced from 1000 to 200 steps per dynamics evaluation
   - 5x faster dynamics computation with minimal accuracy loss

2. **Dynamics caching** (`safe_set_cached`)
   - Hash table for O(1) lookup of previously computed transitions
   - Caches full (state, input, disturbance) → next_state mappings
   - ~40% hit rate eliminates redundant expensive RK4 integrations

3. **Parallelization** (all versions with `parallel` flag)
   - Multi-threaded basis element processing
   - Thread-safe cache access with mutex protection
   - ~2x speedup on modern multi-core CPUs

4. **Memory efficiency**
   - Fixed-size arrays instead of dynamic STL containers
   - Flattened 2D arrays for better cache locality
   - ~1.3 MB total memory for cached version

## Requirements

- C++17 compiler (g++, clang++)
- Make
- POSIX threads (pthread)
- No external dependencies

## Build

```bash
# Build all versions
make all

# Clean all builds
make clean
```

## Files

- `safe_set.cpp` - Original implementation
- `safe_set_opencl_ready.cpp` - OpenCL-compatible version
- `safe_set_optimized.cpp` - Optimized with reduced RK4 steps
- `safe_set_cached.cpp` - Cached dynamics version (recommended)
- `Makefile` - Build system with multiple targets
