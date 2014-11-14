/*
 * random.cpp
 */

#include <set>
#include <list>
#include <vector>
#include <chrono>
#include <random>
#include <cstdlib>

#include "random.hpp"

using namespace std;


Graph * randomDAG(int vertices, int edges, int sources, int sinks) {
	int i, j;
	int vertex;
	int edgeCount;
	int counter = 0;
	double edgeMean = ((double) edges) / ((double) vertices);
	Vertex * sourceVertex, * targetVertex;

	default_random_engine generator(chrono::system_clock::now().time_since_epoch().count());
	poisson_distribution<int> arityDistribution(edgeMean - 1.0);
	uniform_int_distribution<int> intDistribution;
	uniform_int_distribution<int> sourceDistribution(0, vertices - sinks - 1);
	
	Graph * graph = new Graph();

	if ((edges < (vertices - 1)) || (edges > (vertices * (vertices - 1))))
		exit(EXIT_FAILURE);

	// Create all the vertices first
	for (i = 0; i < vertices; i++)
		graph->addVertex();

	// Do a first pass to create edges
	for (i = 0; i < vertices - sinks; i++) {
		int min = i < sources ? sources : i + 1;
		sourceVertex = graph->getVertexFromId(i);

		edgeCount = arityDistribution(generator) + 1;
		if (edgeCount > (vertices - i - 1))
			edgeCount = vertices - i - 1;

		intDistribution = uniform_int_distribution<int>(min, vertices - 1);
		for (j = 0; j < edgeCount; j++) {
			vertex = intDistribution(generator);
			targetVertex = graph->getVertexFromId(vertex);

			graph->addEdge(sourceVertex, targetVertex);
			counter++;
		}
	}

	// Verify that no sources have appeared
	for (i = sources; i < vertices; i++) {
		if (graph->getVertexFromId(i)->getNumberOfPredecessors() == 0) {
			intDistribution = uniform_int_distribution<int> (0, i - 1);
			vertex = intDistribution(generator);

			sourceVertex = graph->getVertexFromId(vertex);
			targetVertex = graph->getVertexFromId(i);
			graph->addEdge(sourceVertex, targetVertex);
			counter++;
		}
	}

	// Adjust the number of edges to the initial request
  while (counter < edges) {
    vertex = sourceDistribution(generator);
    sourceVertex = graph->getVertexFromId(vertex);

    intDistribution = uniform_int_distribution<int>(vertex + 1, vertices - 1);
    vertex = intDistribution(generator);
    targetVertex = graph->getVertexFromId(vertex);
  }

  while (counter > edges) {
    vertex = sourceDistribution(generator);
    sourceVertex = graph->getVertexFromId(vertex);

    if ((edgeCount = sourceVertex->getNumberOfSuccessors()) > 1) {
      for (auto it = sourceVertex->successors_begin(); it != sourceVertex->successors_end(); ++it) {
        if (((*it)->getNumberOfPredecessors() > 1) &&
            (graph->indirectPathExists(sourceVertex, *it))) {
          graph->removeEdge(sourceVertex, *it);
          counter--;
          break;
        }
      }
    }
  }
	
  return graph;
}

