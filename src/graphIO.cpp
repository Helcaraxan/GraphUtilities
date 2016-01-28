#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "graph-utilities/implementation/support.hpp"
#include "graph-utilities/implementation/graphImpl.hpp"

using namespace std;


// Parser functions

Graph *
parseDotFile(const char * fileName, bool noDoubleEdges) {
  char dump[128];
  int source, target, maxId;
  fstream input(fileName, fstream::in);
  GraphImpl * graph = new GraphImpl();

  if (!input.good()) {
    cerr << "ERROR: Could not open the Dot input file.\n";
    exit(EXIT_FAILURE);
  }

  input.getline(dump, 127);
  if (!strstr(dump, "digraph")) {
    cerr << "ERROR: The supplied file is not a graph in Dot format.\n";
    exit(EXIT_FAILURE);
  }

  cout << "Parsing graph...";
  while (input.good()) {
    input.getline(dump, 127);

    if (strchr(dump, '}')) {
      break;
    }

    if (strchr(dump, '>')) {
      sscanf(dump, "%d -> %d", &source, &target);
      maxId = source < target ? target : source;
      while (graph->getVertexCount() <= (unsigned) maxId)
        graph->addVertexUnsafe(graph->getThreadCount());

      if (noDoubleEdges)
        graph->addEdgeUnsafe(graph->vertices[source], graph->vertices[target]);
      else
        graph->addEdge(graph->vertices[source], graph->vertices[target]);
    } else {
      sscanf(dump, "%d", &source);
      while (graph->getVertexCount() <= (unsigned) source)
        graph->addVertex();
    }
  }
  input.close();
  graph->setIndexed(false);

  cout << " DONE" << endl;

  // Make sure the graph is a DAG
  if (graph->isDAG())
    graph->setCondensed(true);
  else
    graph->setCondensed(false);

  return graph;
}


Graph *
parseGraFile(const char * fileName, bool noDoubleEdges) {
  char dump[128];
  int source, target, lineNumber;
  string barTitle = "Parsing graph ";
  fstream input(fileName, fstream::in);
  GraphImpl * graph = new GraphImpl();

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
    graph->addVertexUnsafe(graph->getThreadCount());

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
        graph->addEdgeUnsafe(graph->vertices[source], graph->vertices[target]);
      else
        graph->addEdge(graph->vertices[source], graph->vertices[target]);
    }

    resultProgressBar(i + 1);
  }
  cout << "\n";

  input.close();
  graph->setIndexed(false);

  // Make sure the graph is a DAG
  if (!graph->isDAG())
    graph->setCondensed(true);
  else
    graph->setCondensed(false);

  return graph;
}


// Serialization

void
GraphImpl::printToFile(const char * fileName, bool withWeights) {
  fstream printStream;

  printStream.open(fileName, ios_base::out);

  printStream << "graph_for_greach\n" << getVertexCount();
  for (auto nodeIt = vertices.begin(), end = vertices.end();
      nodeIt != end; ++nodeIt) {
    if (*nodeIt == nullptr)
      continue;

    printStream << "\n" << (*nodeIt)->getId() << ":";

    for (auto succIt = (*nodeIt)->succBegin(), succEnd = (*nodeIt)->succEnd();
        succIt != succEnd; ++succIt)
      printStream << " " << (*succIt)->getId();

    printStream << " #";

    if (!withWeights)
      continue;

    printStream << "\n" << (*nodeIt)->getWeight() << ":";

    for (int i = 0, e = (*nodeIt)->getSuccessorCount(); i < e; i++)
      printStream << " " << (*nodeIt)->getSuccessorWeight(i);

    printStream << " #";
  }

  printStream.close();
}
