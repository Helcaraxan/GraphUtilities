include_directories("${PROJECT_BINARY_DIR}/include")
link_directories("${PROJECT_BINARY_DIR}/src")


# Set the libraries necessary to build the driver
set(USED_LIBRARIES GraphUtilities /home/dvanamst/opt/lib/libpapi.so stdc++ tbb)


# Define the driver executable
add_executable(QueryDriver driver.cpp)

target_link_libraries(QueryDriver ${USED_LIBRARIES})


# Define the test directories
set(GRAPH_DIR ${PROJECT_SOURCE_DIR}/test/graphs)
set(QUERY_DIR ${PROJECT_SOURCE_DIR}/test/queries)


# Macro to add all tests on a specified graph
macro(graph_test filepath)
	get_filename_component(graph ${filepath} NAME_WE)
	get_filename_component(graphfile ${filepath} NAME)

	# Random query DFS test
	add_test(${graph}-random-DFS
		QueryDriver
		--input=${GRAPH_DIR}/${graphfile} --test=${QUERY_DIR}/${graph}.test --search=DFS --unique-edges --batch)

	# Random query BBFS test
	add_test(${graph}-random-BBFS
		QueryDriver
		--input=${GRAPH_DIR}/${graphfile} --test=${QUERY_DIR}/${graph}.test --search=BBFS --unique-edges --batch)

	# Positive query test with automatic DFS / BBFS choice
	add_test(${graph}-positive-automatic
		QueryDriver
		--input=${GRAPH_DIR}/${graphfile} --test=${QUERY_DIR}/${graph}_positive.test --unique-edges --batch)

	# Deep positive query test with automatic DFS / BBFS choice
	add_test(${graph}-deep-positive-automatic
		QueryDriver
		--input=${GRAPH_DIR}/${graphfile} --test=${QUERY_DIR}/${graph}_positive.test  --unique-edges --batch)

	# Set properties for random query tests
	set_tests_properties(${graph}-random-DFS ${graph}-random-BBFS
		PROPERTIES
		FAIL_REGULAR_EXPRESSIONS "*Erroneous result for*"
		TIMEOUT 300)

	# Set properties for positive query tests
	set_tests_properties(${graph}-positive-automatic ${graph}-deep-positive-automatic
		PROPERTIES
		FAIL_REGULAR_EXPRESSIONS "*Erroneous result for*"
		PASS_REGULAR_EXPRESSIONS "*- Number of negative answers   : 0*"
		TIMEOUT 300)
endmacro(graph_test)


# Test 1. Usage message
add_test(QueryDriverUsage QueryDriver --help)

set_tests_properties(QueryDriverUsage
		PROPERTIES
		PASS_REGULAR_EXPRESSION "Usage:*")

# Test 2. Graph tests
file(GLOB GRAPHS ${GRAPH_DIR}/*)

foreach(testgraph ${GRAPHS})
	graph_test(${testgraph})
endforeach(testgraph)

