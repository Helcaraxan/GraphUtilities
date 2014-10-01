cmake_minimum_required(VERSION 2.6)
project(GraphUtilities)

set(CMAKE_CXX_FLAGS "-std=c++11 -pedantic -Wall -g3")
link_directories("${CMAKE_HOME_DIR}/lib")

add_subdirectory(src)


