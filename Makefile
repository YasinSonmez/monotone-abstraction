# Unified Safe Set Computation Makefile

CXX = g++
CXXFLAGS = -std=c++17 -O3 -Wall -Wextra

TARGET = safe_set
TARGET_OPENCL_READY = safe_set_opencl_ready
TARGET_OPTIMIZED = safe_set_optimized
TARGET_CACHED = safe_set_cached
SOURCE = safe_set.cpp
SOURCE_OPENCL_READY = safe_set_opencl_ready.cpp
SOURCE_OPTIMIZED = safe_set_optimized.cpp
SOURCE_CACHED = safe_set_cached.cpp

# Default target
all: $(TARGET) $(TARGET_OPENCL_READY) $(TARGET_OPTIMIZED) $(TARGET_CACHED)

# Build targets
$(TARGET): $(SOURCE)
	$(CXX) $(CXXFLAGS) -pthread $(SOURCE) -o $(TARGET)

$(TARGET_OPENCL_READY): $(SOURCE_OPENCL_READY)
	$(CXX) $(CXXFLAGS) -pthread $(SOURCE_OPENCL_READY) -o $(TARGET_OPENCL_READY)


$(TARGET_OPTIMIZED): $(SOURCE_OPTIMIZED)
	$(CXX) $(CXXFLAGS) -pthread $(SOURCE_OPTIMIZED) -o $(TARGET_OPTIMIZED)

$(TARGET_CACHED): $(SOURCE_CACHED)
	$(CXX) $(CXXFLAGS) -pthread $(SOURCE_CACHED) -o $(TARGET_CACHED)

# Run targets
run-2d-seq: $(TARGET)
	./$(TARGET) 2d precomp seq

run-2d-parallel: $(TARGET)
	./$(TARGET) 2d precomp parallel

run-3d-seq: $(TARGET)
	./$(TARGET) 3d precomp seq

run-3d-parallel: $(TARGET)
	./$(TARGET) 3d precomp parallel

run-3d-seq-noprecomp: $(TARGET)
	./$(TARGET) 3d noprecomp seq

run-3d-parallel-noprecomp: $(TARGET)
	./$(TARGET) 3d noprecomp parallel

# OpenCL-ready version targets
run-2d-seq-opencl: $(TARGET_OPENCL_READY)
	./$(TARGET_OPENCL_READY) 2d precomp seq

run-2d-parallel-opencl: $(TARGET_OPENCL_READY)
	./$(TARGET_OPENCL_READY) 2d precomp parallel

run-3d-seq-opencl: $(TARGET_OPENCL_READY)
	./$(TARGET_OPENCL_READY) 3d precomp seq

run-3d-parallel-opencl: $(TARGET_OPENCL_READY)
	./$(TARGET_OPENCL_READY) 3d precomp parallel

run-3d-seq-noprecomp-opencl: $(TARGET_OPENCL_READY)
	./$(TARGET_OPENCL_READY) 3d noprecomp seq

run-3d-parallel-noprecomp-opencl: $(TARGET_OPENCL_READY)
	./$(TARGET_OPENCL_READY) 3d noprecomp parallel

# Performance comparison
compare-3d: $(TARGET) $(TARGET_OPENCL_READY)
	@echo "=== Sequential 3D Computation (Original) ==="
	@time ./$(TARGET) 3d precomp seq
	@echo ""
	@echo "=== Sequential 3D Computation (OpenCL-Ready) ==="
	@time ./$(TARGET_OPENCL_READY) 3d precomp seq

compare-2d: $(TARGET) $(TARGET_OPENCL_READY)
	@echo "=== Sequential 2D Computation (Original) ==="
	@time ./$(TARGET) 2d precomp seq
	@echo ""
	@echo "=== Sequential 2D Computation (OpenCL-Ready) ==="
	@time ./$(TARGET_OPENCL_READY) 2d precomp seq

# Optimized version targets
run-3d-optimized-seq: $(TARGET_OPTIMIZED)
	./$(TARGET_OPTIMIZED) 3d precomp seq

run-3d-optimized-parallel: $(TARGET_OPTIMIZED)
	./$(TARGET_OPTIMIZED) 3d precomp parallel

# Cached version targets
run-3d-cached-seq: $(TARGET_CACHED)
	./$(TARGET_CACHED) 3d precomp seq

run-3d-cached-parallel: $(TARGET_CACHED)
	./$(TARGET_CACHED) 3d precomp parallel

run-2d-cached-seq: $(TARGET_CACHED)
	./$(TARGET_CACHED) 2d precomp seq

run-2d-cached-parallel: $(TARGET_CACHED)
	./$(TARGET_CACHED) 2d precomp parallel

# Performance benchmark
benchmark: $(TARGET_OPTIMIZED) $(TARGET_CACHED)
	@echo "=== PERFORMANCE BENCHMARK: 3D Safe Set Computation ==="
	@echo ""
	@echo "1. Optimized Sequential:"
	@time ./$(TARGET_OPTIMIZED) 3d noprecomp seq 2>&1 | tail -5
	@echo ""
	@echo "2. Optimized Parallel:"
	@time ./$(TARGET_OPTIMIZED) 3d noprecomp parallel 2>&1 | tail -5
	@echo ""
	@echo "3. Cached Sequential:"
	@time ./$(TARGET_CACHED) 3d noprecomp seq 2>&1 | tail -10
	@echo ""
	@echo "4. Cached Parallel:"
	@time ./$(TARGET_CACHED) 3d noprecomp parallel 2>&1 | tail -10

# Clean target
clean:
	rm -f $(TARGET) $(TARGET_OPENCL_READY) $(TARGET_OPTIMIZED) $(TARGET_CACHED)

.PHONY: all clean \
	run-2d-seq run-2d-parallel run-3d-seq run-3d-parallel \
	run-3d-seq-noprecomp run-3d-parallel-noprecomp \
	run-2d-seq-opencl run-2d-parallel-opencl \
	run-3d-seq-opencl run-3d-parallel-opencl \
	run-3d-seq-noprecomp-opencl run-3d-parallel-noprecomp-opencl \
	run-3d-optimized-seq run-3d-optimized-parallel \
	run-3d-cached-seq run-3d-cached-parallel \
	run-2d-cached-seq run-2d-cached-parallel \
	compare-3d compare-2d benchmark