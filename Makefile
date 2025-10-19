# Safe Set Computation Makefile

CXX = g++
CXXFLAGS = -std=c++17 -O3 -Wall -Wextra

TARGET = safe_set
SOURCE = safe_set.cpp

# Default target
all: $(TARGET)

# Build target
$(TARGET): $(SOURCE)
	$(CXX) $(CXXFLAGS) $(SOURCE) -o $(TARGET)

# Run targets
run-2d: $(TARGET)
	./$(TARGET) 2d precomp

run-3d: $(TARGET)
	./$(TARGET) 3d precomp

run-3d-noprecomp: $(TARGET)
	./$(TARGET) 3d noprecomp

# Clean target
clean:
	rm -f $(TARGET)

.PHONY: all run-2d run-3d run-3d-noprecomp clean