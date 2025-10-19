# Unified Safe Set Computation Makefile

CXX = g++
CXXFLAGS = -std=c++17 -O3 -Wall -Wextra

TARGET = safe_set
SOURCE = safe_set.cpp

# Default target
all: $(TARGET)

# Build target
$(TARGET): $(SOURCE)
	$(CXX) $(CXXFLAGS) -pthread $(SOURCE) -o $(TARGET)

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

# Performance comparison
compare-3d:
	@echo "=== Sequential 3D Computation ==="
	@time ./$(TARGET) 3d precomp seq
	@echo ""
	@echo "=== Parallel 3D Computation ==="
	@time ./$(TARGET) 3d precomp parallel

compare-2d:
	@echo "=== Sequential 2D Computation ==="
	@time ./$(TARGET) 2d precomp seq
	@echo ""
	@echo "=== Parallel 2D Computation ==="
	@time ./$(TARGET) 2d precomp parallel

# Clean target
clean:
	rm -f $(TARGET)

.PHONY: all run-2d-seq run-2d-parallel run-3d-seq run-3d-parallel run-3d-seq-noprecomp run-3d-parallel-noprecomp compare-3d compare-2d clean