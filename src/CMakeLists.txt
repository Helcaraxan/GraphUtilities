include_directories("${PROJECT_BINARY_DIR}/include")

add_library(
		GraphUtilities SHARED
		graph.cpp
		)

install(TARGETS GraphUtilities
		LIBRARY DESTINATION "lib"
	   )
