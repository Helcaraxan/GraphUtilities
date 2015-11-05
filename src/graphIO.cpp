#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "graph-utilities/defs.hpp"
#include "graph-utilities/graph.hpp"
#include "graph-utilities/vertex.hpp"
#include "graph-utilities/support.hpp"

using namespace std;


// Parser functions

Graph *
Graph::createFromDotFile(const char * fileName, bool noDoubleEdges) {
  char dump[128];
  int source, target, maxId;
  fstream input(fileName, fstream::in);
  Graph * graph = NULL;

  graph = new Graph();

  if (!input.good()) {
    cerr << "ERROR: Could not open the Dot input file.\n";
    exit(EXIT_FAILURE);
  }

  input.getline(dump, 127);
  if (!strstr(dump, "digraph")) {
    cerr << "ERROR: The supplied file is not a graph in Dot format.\n";
    exit(EXIT_FAILURE);
  }

  while (input.good()) {
    input.getline(dump, 127);

    if (strchr(dump, '}')) {
      break;
    }

    if (strchr(dump, '>')) {
      sscanf(dump, "%d -> %d", &source, &target);
      maxId = source < target ? target : source;
      while (graph->getVertexCount() <= (unsigned) maxId)
        graph->addVertexUnsafe(graph->threadCount);

      if (noDoubleEdges)
        graph->addEdgeUnsafe(graph->getVertexFromId(source),
            graph->getVertexFromId(target));
      else
        graph->addEdge(graph->getVertexFromId(source),
            graph->getVertexFromId(target));
    } else {
      sscanf(dump, "%d", &source);
      while (graph->getVertexCount() <= (unsigned) source)
        graph->addVertex();
    }
  }
  input.close();
  graph->indexed = false;

  // Make sure the graph is a DAG
  if (!graph->condenseGraph(true))
    graph->condensed = true;
  else
    graph->condensed = false;

  return graph;
}


Graph *
Graph::createFromGraFile(const char * fileName, bool noDoubleEdges) {
  char dump[128];
  int source, target, lineNumber;
  string barTitle = "Parsing graph ";
  fstream input(fileName, fstream::in);
  Graph *  graph = NULL;

  graph = new Graph();

  if (!input.good()) {
    cerr << "ERROR: Could not open the Gra input file.\n";
    exit(EXIT_FAILURE);
  }

  // Throw first lines away while they are text
  do {
    input.getline(dump, 127);
  } while (!isdigit(dump[0]));


  // Get the number of vertices / lines to read
  lineNumber = atoi(dump);
  configureProgressBar(&barTitle, lineNumber);

  // Create all vertices
  for (int i = 0; i < lineNumber; i++)
    graph->addVertexUnsafe(graph->threadCount);

  // Parse the adjacency list
  for (int i = 0; i < lineNumber; i++) {
    // Get the source at the start of the line
    input.get(dump, 127, ' ');
    source = atoi(dump);

    while (true) {
      if (input.get() != ' ') {
        cerr << "ERROR: Could not correctly read the graph definition.\n";
        exit(EXIT_FAILURE);
      }

      if (input.peek() == '#') {
        input.get(dump, 127);
        break;
      }

      input.get(dump, 127, ' ');
      target = atoi(dump);

      if (noDoubleEdges)
        graph->addEdgeUnsafe(graph->getVertexFromId(source),
            graph->getVertexFromId(target));
      else
        graph->addEdge(graph->getVertexFromId(source),
            graph->getVertexFromId(target));
    }

    resultProgressBar(i + 1);
  }
  cout << "\n";

  input.close();
  graph->indexed = false;

  // Make sure the graph is a DAG
  if (!graph->condenseGraph(true))
    graph->condensed = true;
  else
    graph->condensed = false;

  return graph;
}


// Serialization

void
Graph::printToFile(fstream &printStream, bool withWeights) {
  printStream << "graph_for_greach\n" << getVertexCount();
  for (auto nodeIt = vertices.begin(), end = vertices.end();
      nodeIt != end; ++nodeIt) {
    if (*nodeIt == NULL)
      continue;

    printStream << "\n" << (*nodeIt)->id << ":";

    for (auto succIt = (*nodeIt)->succ_begin(), succEnd = (*nodeIt)->succ_end();
        succIt != succEnd; ++succIt)
      printStream << " " << (*succIt)->id;

    printStream << " #";

    if (!withWeights)
      continue;

    printStream << "\n" << (*nodeIt)->weight << ":";

    for (auto weightIt = (*nodeIt)->successorWeights.begin(),
        weightEnd = (*nodeIt)->successorWeights.end();
        weightIt != weightEnd; ++weightIt)
      printStream << " " << (*weightIt);

    printStream << " #";
  }
}


