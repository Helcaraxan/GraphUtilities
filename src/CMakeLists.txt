include_directories("${PROJECT_BINARY_DIR}/include")

add_library(
		GraphUtilities SHARED
		metis-extra.cpp
		queries.cpp
		graph.cpp
		)

target_link_libraries(
	GraphUtilities
	mtmetis
	m
	z
	)

install(TARGETS GraphUtilities
		LIBRARY DESTINATION "lib"
	    )
