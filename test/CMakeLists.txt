include_directories("${PROJECT_BINARY_DIR}/include")

link_directories(
	"${PROJECT_SOURCE_DIR}/lib"
	"${PROJECT_BINARY_DIR}/src"
	)


## Reachability Query driver

# Define the driver executable
add_executable(ReachabilityDriver reachabilityEval.cpp)

target_link_libraries(ReachabilityDriver
	GraphUtilities
	stdc++
	tbb
	)

# Set the install parameters
install(TARGETS ReachabilityDriver
	RUNTIME DESTINATION "bin"
	)


# Define the test directories
set(GRAPH_DIR ${PROJECT_SOURCE_DIR}/test/graphs)
set(QUERY_DIR ${PROJECT_SOURCE_DIR}/test/queries)


# Macro to add all tests on a specified graph
macro(graph_test filepath)
	get_filename_component(graph ${filepath} NAME_WE)
	get_filename_component(graphfile ${filepath} NAME)

	# Random query DFS test
	add_test(${graph}-random-DFS
		ReachabilityDriver
		--input=${GRAPH_DIR}/${graphfile} --test=${QUERY_DIR}/${graph}.test --search=DFS --unique-edges --batch)

	# Random query BBFS test
	add_test(${graph}-random-BBFS
		ReachabilityDriver
		--input=${GRAPH_DIR}/${graphfile} --test=${QUERY_DIR}/${graph}.test --search=BBFS --unique-edges --batch)

	# Positive query test with automatic DFS / BBFS choice
	add_test(${graph}-positive-automatic
		ReachabilityDriver
		--input=${GRAPH_DIR}/${graphfile} --test=${QUERY_DIR}/${graph}_positive.test --unique-edges --batch)

	# Deep positive query test with automatic DFS / BBFS choice
	add_test(${graph}-deep-positive-automatic
		ReachabilityDriver
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
add_test(ReachabilityDriverUsage ReachabilityDriver --help)

set_tests_properties(ReachabilityDriverUsage
		PROPERTIES
		PASS_REGULAR_EXPRESSION "Usage:*")

# Test 2. Graph tests
file(GLOB GRAPHS ${GRAPH_DIR}/*.gra)

foreach(testgraph ${GRAPHS})
	graph_test(${testgraph})
endforeach(testgraph)


## Partition evaluation driver

# Define the driver executable
add_executable(PartitionDriver partitionEval.cpp)

target_link_libraries(PartitionDriver
	GraphUtilities
	stdc++
	tbb
	)

# Set the install parameters
install(TARGETS PartitionDriver
	RUNTIME DESTINATION "bin"
	)
