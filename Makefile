# Makefile for 2D Safe Set Computation

CXX = g++
CXXFLAGS = -std=c++17 -O3 -Wall -Wextra
TARGET = compute_2D_safe_set
SOURCE = compute_2D_safe_set.cpp

$(TARGET): $(SOURCE)
	$(CXX) $(CXXFLAGS) $(SOURCE) -o $(TARGET)

clean:
	rm -f $(TARGET)

run: $(TARGET)
	./$(TARGET)

.PHONY: clean run
