# C++ Implementation of 2D Safe Set Computation

This is a clean, fast, object-oriented C++ implementation of the MATLAB `compute_2D_safe_set.m` example, computing invariant safe sets for adaptive cruise control using monotone abstractions.

## Features

- **Clean OOD Design**: Separate classes for `MonotoneDynamics` and `MonotoneAbstraction`
- **Fast Performance**: Optimized with caching and efficient data structures
- **Single File**: Complete implementation in one file for easy compilation
- **Memory Efficient**: Uses `std::vector` and `std::unordered_map` for optimal performance
- **No Dependencies**: Pure C++17 standard library

## Algorithm

The implementation computes invariant safe sets for vehicle following dynamics:

- **State**: `[h, v]` where `h` = headway, `v` = ego velocity  
- **Input**: `a` (ego acceleration)
- **Disturbance**: `v_L` (lead velocity)

**Dynamics**:
```
h(t+dt) = h(t) + v_L*dt - v(t)*dt - 0.5*a*dt^2
v(t+dt) = v(t) + a*dt
```

The algorithm uses:
- Grid-based monotone abstractions
- Boundary basis representation for safe sets
- Priority directions for efficient computation
- Cached transitions for performance

## Compilation

```bash
# Using make
make

# Or directly with g++
g++ -std=c++17 -O3 -Wall -Wextra compute_2D_safe_set.cpp -o compute_2D_safe_set
```

## Usage

```bash
# Run the computation
./compute_2D_safe_set

# Or using make
make run
```

## Performance

- **Computation Time**: ~43ms on modern hardware
- **Iterations**: 368 iterations to convergence
- **Safe Set Size**: 168 basis elements
- **Memory**: Efficient caching with `std::unordered_map`

## Output

The program outputs:
- Iteration progress
- Final safe set basis elements with [headway, velocity] coordinates
- Computation time
- Number of basis elements

## Example Output

```
=== 2D Safe Set Computation for Vehicle Following ===
Computing invariant safe sets using monotone abstractions...
Starting invariant set computation...
...
Converged! All basis elements are safe.
Safe set computation completed in 368 iterations.
Final safe set has 168 basis elements.

Computation completed in 43 ms
```

## Key Classes

### `MonotoneDynamics`
- Represents the vehicle following dynamics
- Handles state transitions with exact discretization
- Validates input dimensions

### `MonotoneAbstraction`  
- Grid-based abstraction for monotone systems
- Implements boundary basis representation
- Caches transitions for performance
- Computes invariant safe sets using iterative algorithm

## Parameters

- **State Range**: h ∈ [5, 80] m, v ∈ [0, 20] m/s
- **Input Range**: a ∈ [-3, 3] m/s²  
- **Resolution**: h_res = 0.2m, v_res = 0.1 m/s, a_res = 0.1 m/s²
- **Timestep**: dt = 0.1s
- **Priorities**: Smaller headway and larger velocity more unsafe
