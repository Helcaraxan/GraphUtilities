include_directories("../include")

add_library(
		GraphUtilities SHARED
		graph.cpp
		)

install(TARGETS GraphUtilities
		LIBRARY DESTINATION "${PROJECT_SOURCE_DIR}/lib"
	   )
