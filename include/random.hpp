/*
 * random.hpp
 */

#ifndef RANDOM_HPP
#define RANDOM_HPP

#include <cstdlib>

#include "graph.hpp"


bool alternativeConnectivity(Graph * graph, Vertex * source, Vertex * target);
Graph * randomDAG(int vertices, int edges, int sources, int sinks);

#endif//RANDOM_HPP
