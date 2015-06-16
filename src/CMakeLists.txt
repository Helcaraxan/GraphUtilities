include_directories("${PROJECT_BINARY_DIR}/include")


add_library(
		GraphUtilities SHARED
		queries.cpp
		graph.cpp
		)

install(TARGETS GraphUtilities
		LIBRARY DESTINATION "lib"
	    )
