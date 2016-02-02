#include <cstring>
#include <fstream>
#include <iomanip>
#include <sstream>
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
  int source = -1, target = -1, vertexCount = -1, lineNumber = 0;
  string rawLine;
  string barTitle = "Parsing graph ";
  fstream input(fileName, fstream::in);
  stringstream line;
  GraphImpl * graph = new GraphImpl();

  if (!input.good()) {
    cerr << "ERROR: Could not open the Gra input file.\n";
    exit(EXIT_FAILURE);
  }

  // Throw first lines away while they are text
  do {
    getline(input, rawLine);
    lineNumber++;
  } while (!isdigit(rawLine[0]));


  // Get the number of vertices (lines) to read
  line = stringstream(rawLine);
  line >> vertexCount;
  configureProgressBar(&barTitle, vertexCount);

  // Create all vertices
  for (int i = 0; i < vertexCount; i++)
    graph->addVertexUnsafe(graph->getThreadCount());

  // Parse the adjacency list
  for (int i = 0; i < vertexCount; i++) {
    char peak;

    // Get new line from input and exit if EOF is reached
    getline(input, rawLine);
    if (input.eof())
      break;

    // Get source node ID from the start of the line
    lineNumber++;
    line = stringstream(rawLine);
    line >> source;

    if (line.fail()) {
      cerr << "ERROR: Incorrect source in input file at line " << lineNumber
        << "\n.";
      exit(EXIT_FAILURE);
    }

    // Iterate over the successors of the source
    while (true) {
      // Fast-forward to the next successor or end-of-list
      while (true) {
        line.get(peak);

        if (isdigit(peak) || peak == '#') {
          line.putback(peak);
          break;
        }
      }

      // End-of-list reached. Ignore the rest of the line (comment)
      if (peak == '#')
        break;

      line >> target;
      if (line.fail()) {
        cerr << "ERROR: Incorrect target in input file at line " << lineNumber
          << "\n.";
        exit(EXIT_FAILURE);
      }

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
