#include <iostream>

#include "graph-utilities/implementation/graphImpl.hpp"

using namespace std;


// File local functions
static PartitionNode *
computeHierarchicalPartitioning(const GraphImpl * graph,
    vector<int>& scheduling, int start = -1, int end = -1) {
  PartitionNode * result = NULL;

  if (start == -1) {
    start = 0;
    end = scheduling.size() - 1;
  }

  if (start == end) {
    result = new PartitionNode(scheduling[start]);
    result->maxLive = 1;
  } else {
    int IOcost = 0;
    int intermediatePos = (start + end) / 2 + (start + end) % 2;
    set<int> targetSet(scheduling.begin() + intermediatePos,
        scheduling.begin() + end);

    PartitionNode * first = computeHierarchicalPartitioning(graph, scheduling,
        start, intermediatePos - 1);
    PartitionNode * second = computeHierarchicalPartitioning(graph, scheduling,
        intermediatePos, end);
    
    result = new PartitionNode(first, second);

    // Compute the IO size between the two children
    for (int i = start; i < intermediatePos; i++) {
      Vertex * curr = graph->getVertex(scheduling[i]);

      for (int j = 0, e = curr->getSuccessorCount(); j < e; j++) {
        if (targetSet.count(curr->getSuccessor(j)->getId())) {
          IOcost++;
          break;
        }
      }
    }

    // Compute the maxLive of this partition level as :
    // maxLive = max(maxLive(first), maxLive(second), IO(first, second))
    result->maxLive = first->maxLive;
    if (second->maxLive > result->maxLive)
      result->maxLive = second->maxLive;

    if (IOcost > result->maxLive)
      result->maxLive = IOcost;
  }

  return result;
}


static int
assignPartitions(PartitionNode * node, int memorySize,
    map<int, int>& partitions, int partitionId = 0) {
  if (node->id != -1) {
    partitions[node->id] = partitionId;
  } else {
    partitionId =
      assignPartitions(node->first, memorySize, partitions, partitionId);

    if (node->maxLive > memorySize)
      partitionId++;

    partitionId =
      assignPartitions(node->second, memorySize, partitions, partitionId);

    if (node->maxLive > memorySize)
      partitionId++;
  }

  return partitionId;
}


// Partitioning queries

int
GraphImpl::getPartitionCount() const {
  return 1;
}


const Vertex::PartitionArray&
GraphImpl::getPartitionSet(void) const {
  return sccSets;
}


int
GraphImpl::getSchedulingCost(vector<int>& scheduling, int memorySize) const {
  int cost = 0;
  map<int, int> partitions;
  PartitionNode * partitionRoot = NULL;

  // First check the validity of the scheduling
  if (!checkSchedulingValidity(scheduling))
    return 0;

  // Construct the default hierarchical partitioning of the scheduling.
  partitionRoot = computeHierarchicalPartitioning(this, scheduling);

  // Obtain a partition assignement for the specified memory-size
  assignPartitions(partitionRoot, memorySize, partitions);

  // Compute the IO cost for each partition
  for (auto it = partitions.begin(), end = partitions.end(); it != end; ++it) {
    set<int> targetSet;
    Vertex * curr = getVertex(it->first);

    for (int i = 0, e = curr->getSuccessorCount(); i < e; i++)
      targetSet.insert(curr->getSuccessor(i)->getId());

    cost += targetSet.size();
  }

  return cost;
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


bool
GraphImpl::checkSchedulingValidity(vector<int>& scheduling) const {
  set<int> idSet;

  if (scheduling.size() != getVertexCount()) {
    cerr << "WARNING: The supplied scheduling in 'getSchedulingCost' is not "
      << "valid. It does not contain the same number of vertices as the "
      << "graph.\n";
    return false;
  }

  for (int i = 0, e = scheduling.size(); i < e; i++) {
    Vertex * curr = getVertex(scheduling[i]);

    if (!curr) {
      cerr << "WARNING: The supplied scheduling in 'getSchedulingCost' contains"
        << " an ID that does not correspond to an existing vertex.\n";
      return false;
    }

    if (idSet.count(scheduling[i])) {
      cerr << "WARNING: The supplied scheduling in 'getSchedulingCost' contains"
        << " two times the same ID (" << scheduling[i] << ").\n";
      return false;
    }

    for (auto j = 0, e2 = curr->getSuccessorCount(); j < e2; j++) {
      if (idSet.count(curr->getSuccessor(j)->getId())) {
        cerr << "WARNING: The supplied scheduling in 'getSchedulingCost' is "
          << "illegal. ID (" << curr->getSuccessor(j)->getId() << ") should be "
          << "scheduling after ID (" << scheduling[i] << ").\n";
        return false;
      }
    }

    idSet.insert(scheduling[i]);
  }

  return true;
}
