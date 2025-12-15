#!/bin/bash

# Generic script to build and run C++ files with OpenCL
# Usage: ./build_and_run.sh <cpp_file_name> <args...>
# Example: ./build_and_run.sh safe_set_gpu_hybrid 3d

set -e  # Exit on error

if [ $# -lt 2 ]; then
    echo "Usage: $0 <cpp_file_name> <args...>"
    echo "Example: $0 safe_set_gpu_hybrid 3d"
    exit 1
fi

CPP_FILE=$1
shift  # Remove first argument, rest are program arguments
PROGRAM_ARGS="$@"

# Add .cpp extension if not present
if [[ ! "$CPP_FILE" == *.cpp ]]; then
    CPP_SOURCE="${CPP_FILE}.cpp"
    EXECUTABLE="$CPP_FILE"
else
    CPP_SOURCE="$CPP_FILE"
    EXECUTABLE="${CPP_FILE%.cpp}"
fi

# Check if source file exists
if [ ! -f "$CPP_SOURCE" ]; then
    echo "Error: Source file '$CPP_SOURCE' not found"
    exit 1
fi

echo "==================================="
echo "Cleaning previous builds..."
echo "==================================="
rm -f "$EXECUTABLE"

echo ""
echo "==================================="
echo "Compiling $CPP_SOURCE..."
echo "==================================="
clang++ -std=c++17 -O3 -o "$EXECUTABLE" "$CPP_SOURCE" -framework OpenCL

echo ""
echo "==================================="
echo "Running: ./$EXECUTABLE $PROGRAM_ARGS"
echo "==================================="
./"$EXECUTABLE" $PROGRAM_ARGS

echo ""
echo "==================================="
echo "Done!"
echo "==================================="

