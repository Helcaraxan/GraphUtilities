include_directories("${PROJECT_BINARY_DIR}/include")

link_directories("${PROJECT_BINARY_DIR}/lib")


# Generate the public header files according to the previously defined options
file(GLOB PUBLIC_HEADER_LIST ${PROJECT_SOURCE_DIR}/src/include/*.hpp.in)
foreach(PUBLIC_HEADER ${PUBLIC_HEADER_LIST})
	get_filename_component(HEADER ${PUBLIC_HEADER} NAME_WE)

	configure_file(
		"${PROJECT_SOURCE_DIR}/src/include/${HEADER}.hpp.in"
		"${PROJECT_BINARY_DIR}/include/graph-utilities/${HEADER}.hpp"
		)
endforeach(PUBLIC_HEADER)

# Generate the private header files according to the previously defined options
file(GLOB PRIVATE_HEADER_LIST ${PROJECT_SOURCE_DIR}/src/include/implementation/*.hpp.in)
foreach(PRIVATE_HEADER ${PRIVATE_HEADER_LIST})
	get_filename_component(HEADER ${PRIVATE_HEADER} NAME_WE)

	configure_file(
		"${PROJECT_SOURCE_DIR}/src/include/implementation/${HEADER}.hpp.in"
		"${PROJECT_BINARY_DIR}/include/graph-utilities/implementation/${HEADER}.hpp"
		)
endforeach(PRIVATE_HEADER)


# Add the library
file(GLOB SOURCE_LIST ${PROJECT_SOURCE_DIR}/src/*.cpp)
add_library(
	GraphUtilities SHARED
	${SOURCE_LIST}
	)

target_link_libraries(GraphUtilities
	${PROJECT_SOURCE_DIR}/lib/libpatoh.a
	tbb
	m
	z
	)


# Add the install parameters
install(TARGETS GraphUtilities
	LIBRARY DESTINATION "lib"
	)

file(GLOB PUBLIC_HEADERS ${PROJECT_BINARY_DIR}/include/graph-utilities/*.hpp)
install(
	FILES ${PUBLIC_HEADERS}
	DESTINATION "include/graph-utilities"
	)
