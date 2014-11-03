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


bool alternativeConnectivity(Graph * graph, Vertex * source, Vertex * target) {
  Vertex * currVertex;
  list<Vertex *> toVisit;
  set<Vertex *> scheduled;

  // Fill up the initial toVisit list
  for (auto it = source->successors.begin(); it != source->successors.end(); ++it) {
    if (*it != target) {
      toVisit.push_back(*it);
      scheduled.insert(*it);
    }
  }

  while (!toVisit.empty()) {
    currVertex = toVisit.front();
    toVisit.pop_front();

    for (auto it = currVertex->successors.begin(); it != currVertex->successors.end(); ++it) {
      if (*it == target)
        return true;

      if (((*it)->id < target->id) && (scheduled.count(*it) == 0)) {
        toVisit.push_front(*it);
        scheduled.insert(*it);
      }
    }
  }

  return false;
}


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
		sourceVertex = graph->vertices[i];

		edgeCount = arityDistribution(generator) + 1;
		if (edgeCount > (vertices - i - 1))
			edgeCount = vertices - i - 1;

		intDistribution = uniform_int_distribution<int>(min, vertices - 1);
		for (j = 0; j < edgeCount; j++) {
			vertex = intDistribution(generator);
			targetVertex = graph->vertices[vertex];

			graph->addEdge(sourceVertex, targetVertex);
			counter++;
		}
	}

	// Verify that no sources have appeared
	for (i = sources; i < vertices; i++) {
		if (graph->vertices[i]->predecessors.size() == 0) {
			intDistribution = uniform_int_distribution<int> (0, i - 1);
			vertex = intDistribution(generator);

			sourceVertex = graph->vertices[vertex];
			targetVertex = graph->vertices[i];
			graph->addEdge(sourceVertex, targetVertex);
			counter++;
		}
	}

	// Adjust the number of edges to the initial request
  while (counter < edges) {
    vertex = sourceDistribution(generator);
    sourceVertex = graph->vertices[vertex];

    intDistribution = uniform_int_distribution<int>(vertex + 1, vertices - 1);
    vertex = intDistribution(generator);
    targetVertex = graph->vertices[vertex];
  }

  while (counter > edges) {
    vertex = sourceDistribution(generator);
    sourceVertex = graph->vertices[vertex];

    if ((edgeCount = sourceVertex->successors.size()) > 1) {
      for (i = 0; i < edgeCount; i++) {
        targetVertex = sourceVertex->successors[i];
        if ((targetVertex->predecessors.size() > 1) &&
            (alternativeConnectivity(graph, sourceVertex, targetVertex))) {
          graph->removeEdge(sourceVertex, targetVertex);
          counter--;
          break;
        }
      }
    }
  }
	
  return graph;
}

