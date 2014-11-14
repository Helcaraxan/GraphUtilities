include_directories("../include")

add_library(
		GraphUtilities SHARED
		graph.cpp
		random.cpp
		)

install(TARGETS GraphUtilities
		LIBRARY DESTINATION "lib"
	   )
