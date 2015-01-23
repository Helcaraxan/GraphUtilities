cmake_minimum_required(VERSION 2.6)
project(GraphUtilities)

# Define some compile options passed onto the code by #define's
option(ENABLE_STATISTICS
	"Gather statistics from each query about the algorithm's effectiveness" ON)
option(ENABLE_BENCHMARKS
	"Perform benchmark measures on indexing and queries" ON)
option(ENABLE_RETRO_LABELS
	"During labeling perform a 3rd and 4th reversed pass for stricter labeling" OFF)

# Generate the header file according to the previously defined options
if(NOT MAX_THREADS)
	set(MAX_THREADS "2")
endif(NOT MAX_THREADS)

configure_file(
	"${PROJECT_SOURCE_DIR}/include/graph.hpp.in"
	"${PROJECT_BINARY_DIR}/include/graph.hpp"
	)

# Set the compile options depending on the build-type
set(CMAKE_CXX_FLAGS_DEBUG "-std=c++11 -pedantic -Wall -pthread -g3")
set(CMAKE_CXX_FLAGS_RELEASE "-std=c++11 -pedantic -Wall -pthread -O3")

# Include the code for building
add_subdirectory(src)

enable_testing()
add_subdirectory(test)

