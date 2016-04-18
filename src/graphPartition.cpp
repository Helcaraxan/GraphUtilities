#include <atomic>
#include <chrono>
#include <thread>
#include <fstream>
#include <iostream>
#include <algorithm>

#include <tbb/concurrent_queue.h>

#include "graph-utilities/implementation/graphImpl.hpp"
#include "graph-utilities/implementation/patoh.hpp"
#include "graph-utilities/implementation/support.hpp"

using namespace std;


// Global variables for graph partitioning

atomic<int> partCount{0};
bool taskPoison = false;
tbb::concurrent_queue<PartitionTask *> taskQueue;
PartitionImpl * part = nullptr;
PartitionMethod method = UndefinedPartitionMethod;


// Internal main partitioning functions

void
#if defined(ENABLE_TLS)
GraphImpl::processPartitionQuery(PartitionQueryImpl * query) {
#else // ENABLE_TLS
GraphImpl::processPartitionQuery(PartitionQueryImpl * query, int threadId) {
#endif // ENABLE_TLS
  method = query->getMethod();

  switch (query->getMethod()) {
    case Greedy:
      partitionGreedy(query);
      break;

    case Convexify:
      partitionConvexify(query);
      break;

    case MaxDistance:
      partitionMaxDistance(query);
      break;

    case PaToH:
      partitionPaToH(query);
      break;

    case UndefinedPartitionMethod:
      query->setError(true);
      break;
  }

  resultQueue.push(query);
}


HGraph *
GraphImpl::getHGraph(Vertex::IdSet& idSet) const {
  int idx = 0;
  HGraph * hyperGraph = new HGraph();

  hyperGraph->cellCount = idSet.size();
  hyperGraph->backwardMap.resize(idSet.size());

  for (auto it = idSet.begin(), end = idSet.end(); it != end; ++it) {
    Vertex * curr = getVertex(*it);

    hyperGraph->forwardMap[*it] = idx;
    hyperGraph->backwardMap[idx++] = *it;

    for (int i = 0, e = curr->getSuccessorCount(); i < e; i++) {
      if (idSet.count(curr->getSuccessor(i)->getId())) {
        hyperGraph->netCount++;
        break;
      }
    }
  }

  hyperGraph->xpins.push_back(0);

  for (auto it = idSet.begin(), end = idSet.end(); it != end; ++it) {
    Vertex * curr = getVertex(*it);

    hyperGraph->pins.push_back(hyperGraph->forwardMap[*it]);
    for (int i = 0, e = curr->getSuccessorCount(); i < e; i++) {
      int succ = curr->getSuccessor(i)->getId();

      if (idSet.count(succ))
        hyperGraph->pins.push_back(hyperGraph->forwardMap[succ]);
    }

    hyperGraph->xpins.push_back(hyperGraph->pins.size());
  }

  return hyperGraph;
}


void
partitionWorker(const GraphImpl * graph) {
  PartitionTask * task = nullptr;

  while (true) {
    while (!taskQueue.try_pop(task)) {
      if (taskPoison)
        return;
    }

    switch (method) {
      case Greedy:
        cerr << "ERROR: A partition-worker received a Greedy query.\n";
        exit(EXIT_FAILURE);
        break;

      case Convexify:
        graph->convBisect(task->node);

        resultProgressBar(partCount);
        if (partCount == (int) graph->getVertexCount())
          return;

        break;

      case MaxDistance:
        graph->maxDistBisect(task->node, task->sourceSet, task->sinkSet);
        break;

      case PaToH:
        graph->patohBisect(task->node);

        resultProgressBar(partCount);
        if (partCount == (int) graph->getVertexCount())
          return;

        break;

      case UndefinedPartitionMethod:
        return;
    }

    delete task;
  }
}
