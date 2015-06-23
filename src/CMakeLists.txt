include_directories(
	"${PROJECT_BINARY_DIR}/include"
	"${PROJECT_SOURCE_DIR}/include"
	)

link_directories("${PROJECT_SOURCE_DIR}/lib")

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
