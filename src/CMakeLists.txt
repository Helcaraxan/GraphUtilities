include_directories("${PROJECT_BINARY_DIR}/include")

link_directories("${PROJECT_BINARY_DIR}/lib")


# Generate the header files according to the previously defined options
file(GLOB SOURCE_HEADER_LIST ${PROJECT_SOURCE_DIR}/src/include/*.hpp.in)
foreach(SOURCE_HEADER ${SOURCE_HEADER_LIST})
	get_filename_component(HEADER ${SOURCE_HEADER} NAME_WE)

	configure_file(
		"${PROJECT_SOURCE_DIR}/src/include/${HEADER}.hpp.in"
		"${PROJECT_BINARY_DIR}/include/graph-utilities/${HEADER}.hpp"
		)
endforeach(SOURCE_HEADER)


# Add the library
file(GLOB SOURCE_LIST ${PROJECT_SOURCE_DIR}/src/*.cpp)
add_library(
	GraphUtilities SHARED
	${SOURCE_LIST}
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
		${PROJECT_BINARY_DIR}/include/graph-utilities/*.hpp
	DESTINATION "include/graph-utilities"
	)
