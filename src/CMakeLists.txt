include_directories("${PROJECT_BINARY_DIR}/include")

link_directories("${PROJECT_BINARY_DIR}/lib")


# Try to locate the required PAPI library
find_library(PAPI_LIB libpapi.so)

# Adapt to the specific case of idchire
if (HOSTNAME STREQUAL "uv2000")
	set(PAPI_LIB "/opt/papi/5.2.0/lib/libpapi.so")
endif()

if (PAPI_LIB STREQUAL "PAPI_LIB-NOTFOUND")
	message(FATAL_ERROR "Could not find the PAPI library. Please install it before continuing.")
endif()


# Add the library
add_library(
	GraphUtilities SHARED
	metis-extra.cpp
	queries.cpp
	graph.cpp
	)

target_link_libraries(GraphUtilities
	${PAPI_LIB}
	mtmetis
	tbb
	m
	z
	)


# Add the install parameters
install(TARGETS GraphUtilities
	LIBRARY DESTINATION "${PROJECT_SOURCE_DIR}/lib"
	)

install(
	FILES
		"${PROJECT_BINARY_DIR}/include/graph.hpp"
		"${PROJECT_BINARY_DIR}/include/queries.hpp"
	DESTINATION "${PROJECT_SOURCE_DIR}/include"
	)
