#include "graph-utilities/defs.hpp"
#include "graph-utilities/graph.hpp"
#include "graph-utilities/vertex.hpp"

using namespace std;


// Partitioning queries

int
Graph::getPartitionCount() {
  return 1;
}


Graph::PartitionSet *
getPartitionSet(int count) {
  return NULL;
}


// Internal partitioning query functions

void
#ifdef ENABLE_TLS
Graph::processPartitionQuery(Query * query) {
#else // ENABLE_TLS
Graph::processPartitionQuery(int threadId, Query * query) {
#endif // ENABLE_TLS
  // TODO
}


