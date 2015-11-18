#include "graph-utilities/implementation/graphImpl.hpp"

using namespace std;


// Partitioning queries

int
GraphImpl::getPartitionCount() const {
  return 1;
}


const Vertex::PartitionArray&
GraphImpl::getPartitionSet(void) const {
  return sccSets;
}


// Internal partitioning query functions

void
#if defined(ENABLE_TLS)
GraphImpl::processPartitionQuery(PartitionQueryImpl * query) {
#else // ENABLE_TLS
GraphImpl::processPartitionQuery(PartitionQueryImpl * query, int threadId) {
#endif // ENABLE_TLS
  // TODO
}
