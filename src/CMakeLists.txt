include_directories("${PROJECT_BINARY_DIR}/include")

link_directories("${PROJECT_BINARY_DIR}/lib")


# Generate the header files according to the previously defined options
configure_file(
	"${PROJECT_SOURCE_DIR}/src/include/graph.hpp.in"
	"${PROJECT_BINARY_DIR}/include/graph-utilities/graph.hpp"
	)

configure_file(
	"${PROJECT_SOURCE_DIR}/src/include/queries.hpp.in"
	"${PROJECT_BINARY_DIR}/include/graph-utilities/queries.hpp"
	)

configure_file(
	"${PROJECT_SOURCE_DIR}/src/include/metis-extra.hpp.in"
	"${PROJECT_BINARY_DIR}/include/graph-utilities/metis-extra.hpp"
	)

configure_file(
	"${PROJECT_SOURCE_DIR}/src/include/support.hpp.in"
	"${PROJECT_BINARY_DIR}/include/graph-utilities/support.hpp"
	)


# Add the library
add_library(
	GraphUtilities SHARED
	metis-extra.cpp
	queries.cpp
	graph.cpp
	)

target_link_libraries(GraphUtilities
	mtmetis
	tbb
	m
	z
	)


# Add the install parameters
install(TARGETS GraphUtilities
	LIBRARY DESTINATION "lib"
	)

install(
	FILES
		"${PROJECT_BINARY_DIR}/include/graph-utilities/graph.hpp"
		"${PROJECT_BINARY_DIR}/include/graph-utilities/queries.hpp"
		"${PROJECT_BINARY_DIR}/include/graph-utilities/support.hpp"
		"${PROJECT_BINARY_DIR}/include/graph-utilities/metis-extra.hpp"
	DESTINATION "include/graph-utilities"
	)
