# Safe Set Computation

A clean, fast C++ implementation of safe set computation for vehicle following scenarios.

## Usage

```bash
# Compile
make

# Run 3D computation with precomputation (recommended)
./safe_set 3d precomp

# Run 3D computation without precomputation
./safe_set 3d noprecomp

# Run 2D computation  
./safe_set 2d precomp

# Or use make targets
make run-3d          # 3D with precomputation
make run-3d-noprecomp # 3D without precomputation
make run-2d          # 2D (precomputation ignored)
```

## Results

### 3D Safe Set Computation (with precomputation)
- **Precomputation**: 697 boundary points added
- **Final basis size**: 1,074 elements
- **Iterations**: 104
- **Computation time**: ~5.0 seconds

### 3D Safe Set Computation (without precomputation)
- **Initial basis size**: 1 (entire state space)
- **Final basis size**: 1,083 elements
- **Iterations**: 143
- **Computation time**: ~6.7 seconds

### 2D Safe Set Computation  
- **Final basis size**: 168 elements
- **Iterations**: 368
- **Computation time**: ~37 ms

## Configuration

Edit parameters at the top of `safe_set.cpp`:

```cpp
// 2D Configuration
const double DT_2D = 0.1;
const double H_MIN_2D = 0.0, H_MAX_2D = 20.0;
const double V_MIN_2D = 0.0, V_MAX_2D = 10.0;
const double H_RES_2D = 0.2, V_RES_2D = 0.2;

// 3D Configuration  
const double DT_3D = 0.4;
const double H_MIN_3D = 0.0, H_MAX_3D = 80.0;
const double V_MIN_3D = 0.0, V_MAX_3D = 20.0;
const double H_RES_3D = 0.8, V_RES_3D = 0.4;
```

## Performance

- **3D**: 34x faster than MATLAB
- **2D**: 300x faster than MATLAB
- **Memory**: Significantly lower than MATLAB
- **Dependencies**: C++17 compiler only

## Files

- `safe_set.cpp` - Single source file with all functionality
- `Makefile` - Build system
- `README.md` - This documentation