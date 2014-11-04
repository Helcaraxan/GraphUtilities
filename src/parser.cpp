#include <cstdio>
#include <cstring>
#include <cstdlib>

#include "parser.hpp"


using namespace std;

Graph *
graphFromDotFile(char * fileName) {
  char dump[128];
  int source, target;
  FILE * input = fopen(fileName, "r");
  Graph * graph = new Graph();

  if (input == NULL) {
    fprintf(stderr, "Error while opening input Dot file.\n");
    exit(EXIT_FAILURE);
  }

  fgets(dump, 127, input);
  if (!strstr(dump, "digraph")) {
    fprintf(stderr, "Error - the supplied file is not a graph in Dot format.\n");
    exit(EXIT_FAILURE);
  }

  while (!feof(input)) {
    fgets(dump, 127, input);

    if (strchr(dump, '}')) {
      fclose(input);
      return graph;
    }

    if (strchr(dump, '>')) {
      sscanf(dump, "%d -> %d", &source, &target);
      while (graph->vertices.size() <= (unsigned) target)
        graph->addVertex();

      graph->addEdge(graph->vertices[source], graph->vertices[target]);
    } else {
      sscanf(dump, "%d", &source);
      while (graph->vertices.size() <= (unsigned) source)
        graph->addVertex();
    }
  }

  return graph;
}
