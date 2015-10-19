cmake_minimum_required(VERSION 2.6)
project(GraphUtilities)

# Set specific ENV variables for idchire
site_name(HOSTNAME)
if (HOSTNAME STREQUAL "uv2000")
	message(STATUS "Compiling @ ${HOSTNAME} - Setting ENV variables")
	set(ENV{PATH} "/opt/gcc/gcc-4.9.2/bin:$PATH")
	set(ENV{CPP_INCLUDE_DIR "/opt/papi/papi-5.2.0/include")
	set(CMAKE_C_COMPILER "gcc-4.9.2")
	set(CMAKE_CXX_COMPILER "g++-4.9.2")
endif()

# Set the default build type when necessary
if (NOT CMAKE_BUILD_TYPE)
	message(STATUS "Setting build-type to 'Release' as none was specified")
	set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)


# Define some compile options passed onto the code by #define's
option(ENABLE_STATISTICS
	"Gather statistics from each query about the algorithm's effectiveness" ON)
option(ENABLE_BENCHMARKS
	"Perform benchmark measures on indexing and queries" OFF)
option(ENABLE_TLS
	"Use thread local storage for DFS identification" OFF)
option(ENABLE_COARSEN_STATISTICS
	"Gather statistics from graph coarsening operations" OFF)

if (NOT AUTO_THRESHOLD)
	set(AUTO_THRESHOLD "10")
endif(NOT AUTO_THRESHOLD)

# Set the compile options depending on the build-type
set(CMAKE_CXX_FLAGS_DEBUG "-std=c++11 -pedantic -Wall -pthread -g3 -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-std=c++11 -pedantic -Wall -pthread -g -O3")


# Create a target for removing the files installed by the 'install' target
add_custom_target(uninstall
	COMMAND ${CMAKE_COMMAND} -E remove_directory "${CMAKE_INSTALL_PREFIX}/include/graph-utilities"
	COMMAND ${CMAKE_COMMAND} -E remove_directory "${CMAKE_INSTALL_PREFIX}/lib/libGraphUtilities.so"
	COMMAND ${CMAKE_COMMAND} -E remove_directory "${CMAKE_INSTALL_PREFIX}/bin/QueryDriver"
	)


# Include the code for building the library
add_subdirectory(tools)
add_subdirectory(src)

# Include the code for building the driver and testing
enable_testing()
add_subdirectory(test)
