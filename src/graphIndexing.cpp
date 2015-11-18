#include <iostream>

#include "graph-utilities/implementation/graphImpl.hpp"

using namespace std;


// Local perf measurement function

#if defined(__amd64)
static __inline__ uint64_t getClock(void) {
  uint64_t a, d;
  __asm__ volatile ("rdtsc" : "=a" (a), "=d" (d));
  return (d << 32) | a;
}
#else // __amd64
#error "The GraphUtilities library only supports x86-64 architectures"
#endif // __amd64


// Graph traversal methods

Vertex *
GraphImpl::getNextDFS(bool postOrder, bool reverse) {
  static bool reverseMemory = false;
  static VertexImpl * currentVertex = NULL;
  static VertexImpl::iterator originIt = sources.begin();
  VertexImpl * lastVertex = NULL;
  pair<bool, VertexImpl *> next(false, NULL);

  if (reverse != reverseMemory) {
    discoverExtremities();
    reverseMemory = reverse;
    originIt = reverse ? sinks.begin() : sources.begin();
  }

  while (true) {
    // When at the end reinitialize the static variables and signal the end
    if ((originIt ==
          (reverse ? sinks.end() : sources.end())) && (currentVertex == NULL)) {
      discoverExtremities();
      originIt = reverse ? sinks.begin() : sources.begin();
      return NULL;
    }

    lastVertex = currentVertex;
    if (!currentVertex) {
      currentVertex = *originIt;
      originIt++;

      return currentVertex;
    }

    do {
      next = currentVertex->getNextDFS(lastVertex, postOrder, reverse);
      lastVertex = currentVertex;
      currentVertex = next.second;
    } while (!next.first);

    if (currentVertex)
      break;
  }

  return currentVertex;
}


// Indexing

void
GraphImpl::setIndexed(bool value) {
  indexed = value;
}


void
GraphImpl::labelVertices(bool reverse) {
  int traversalMethod = 0;
  int currLabel = 0;
  VertexImpl * nextVertex = NULL;
  VertexImpl::Array * order = reverse ? &predecessorQueue : &successorQueue;

  // Compose the method used for post-order labeling
  if (reverse)
    traversalMethod |= 0x01;

  traversalMethod |= (indexMethod & 0x02);

  while (true) {
    if (!(nextVertex = dynamic_cast<VertexImpl *>(getNextDFS(true, reverse))))
      break;

    order->push_back(nextVertex);

    if (reverse)
      nextVertex->clearSuccessors();
    else
      nextVertex->clearPredecessors();
  }

  // Refill the predecessors or successors of all the vertices
  for (auto it = order->begin(), end = order->end(); it != end; ++it) {
    if (reverse) {
      for (int i = 0, e = (*it)->getPredecessorCount(); i < e; i++) {
        Vertex * predVertex = (*it)->getPredecessor(i);
        int predWeight = (*it)->getPredecessorWeight(i);

        VertexImpl * predVertexImpl = dynamic_cast<VertexImpl *>(predVertex);
        predVertexImpl->addSuccessorUnsafe(*it, predWeight);
      }
    } else {
      for (int i = 0, e = (*it)->getSuccessorCount(); i < e; i++) {
        Vertex * succVertex = (*it)->getSuccessor(i);
        int succWeight = (*it)->getSuccessorWeight(i);

        VertexImpl * succVertexImpl = dynamic_cast<VertexImpl *>(succVertex);
        succVertexImpl->addPredecessorUnsafe(*it, succWeight);
      }
    }
  }

  // Label the vertices in reverse post-order
  for (auto it = order->rbegin(), end = order->rend(); it != end; ++it) {
    if (reverse)
      (*it)->setOrders(-1, currLabel++);
    else
      (*it)->setOrders(currLabel++, -1);
  }
}


void
GraphImpl::indexGraph() {
  unique_lock<mutex> indexLock(indexMutex);

  startGlobalOperation();

  if (indexed)
    return;

  // Make sure we are indexing a DAG
  condenseToDAG();

  if (indexMethod == UndefinedIndexMethod) {
    cerr << "ERROR: Unknown indexing method. Aborting.\n";
    exit(EXIT_FAILURE);
  }

#if defined(ENABLE_BENCHMARKS)
  uint64_t startClock = getClock();
#endif // ENABLE_BENCHMARKS

  // First traversal
  labelVertices(false);

  // Second traversal
  labelVertices(true);

  // When requested reorder the successors and predecessors in all vertices
  if (indexMethod & 0x04) {
    VertexImpl * curr;

    while (!successorQueue.empty()) {
      curr = successorQueue.back();
      successorQueue.pop_back();

      for (auto it = curr->predBegin(), end = curr->predEnd();
          it != end; ++it)
        (*it)->addSuccessorUnsafe(curr, curr->getPredecessorWeight(*it));

      curr->clearSuccessors();
    }

    while (!predecessorQueue.empty()) {
      curr = predecessorQueue.back();
      predecessorQueue.pop_back();

      for (auto it = curr->succBegin(), end = curr->succEnd();
          it != end; ++it)
        (*it)->addPredecessorUnsafe(curr, curr->getSuccessorWeight(*it));

      curr->clearPredecessors();
    }
  } else {
    successorQueue.clear();
    predecessorQueue.clear();
  }

#if defined(ENABLE_BENCHMARKS)
  cyclesSpentIndexing = getClock() - startClock;
#endif // ENABLE_BENCHMARKS

  indexed = true;

  stopGlobalOperation();
}


// Maintenance

void
GraphImpl::discoverExtremities() {
  sources.clear();
  sinks.clear();

  for (auto it = vertices.begin(), end = vertices.end(); it != end; ++it) {
    if (*it == NULL)
      continue;

    if ((*it)->getPredecessorCount() == 0)
      sources.push_back(*it);

    if ((*it)->getSuccessorCount() == 0)
      sinks.push_back(*it);
  }
}
