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
GraphImpl::getNextDFS(TraversalOrder order, TraversalDirection direction) {
  // Static variables for statefull processing
  static Graph::TraversalOrder orderMemory = PostOrder;
  static Graph::TraversalDirection directionMemory = Forward;
  static VertexImpl * currentVertex = NULL;
  static VertexImpl::iterator originIt = sources.begin();

  // Local temporary variables
  VertexImpl * lastVertex = NULL;
  pair<bool, VertexImpl *> next(false, NULL);

  // When necessary reset the state
  if ((order != NAOrder) || (direction != NADirection)) {
    if (order != NAOrder)
      orderMemory = order;

    if (direction != NADirection)
      directionMemory = direction;

    for (auto it = vertices.begin(), end = vertices.end(); it != end; ++it)
      if (*it)
        (*it)->resetDFS();

    discoverExtremities();

    if (directionMemory == Forward)
      originIt = sources.begin();
    else
      originIt = sinks.begin();
  }

  while (true) {
    // When at the end reinitialize the static variables and signal the end
    if (directionMemory == Forward) {
      if ((originIt == sources.end()) && (currentVertex == NULL)) {
        discoverExtremities();

        originIt = sources.begin();
        return NULL;
      }
    } else if (directionMemory == Backward) {
      if ((originIt == sinks.end()) && (currentVertex == NULL)) {
        discoverExtremities();

        originIt = sinks.begin();
        return NULL;
      }
    }

    lastVertex = currentVertex;
    if (!currentVertex) {
      currentVertex = *originIt;
      originIt++;

      return currentVertex;
    }

    do {
      next =
        currentVertex->getNextDFS(lastVertex, orderMemory, directionMemory);
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
GraphImpl::labelVertices(TraversalDirection direction) {
  int traversalMethod = 0;
  int currLabel = 0;
  VertexImpl * nextVertex = NULL;
  VertexImpl::Array * order = NULL;
  
  switch(direction) {
    case Forward: order = &successorQueue; break;
    case Backward: order = &predecessorQueue; break;
    default:
      exit(EXIT_FAILURE);
      break;
  }

  // Compose the method used for post-order labeling
  if (direction == Backward)
    traversalMethod |= 0x01;

  traversalMethod |= (indexMethod & 0x02);

  nextVertex = dynamic_cast<VertexImpl *>(getNextDFS(PostOrder, direction));
  while (true) {
    if (!nextVertex)
      break;

    order->push_back(nextVertex);

    if (direction == Forward)
      nextVertex->clearPredecessors();
    else
      nextVertex->clearSuccessors();

    nextVertex = dynamic_cast<VertexImpl *>(getNextDFS());
  }

  // Refill the predecessors or successors of all the vertices
  for (auto it = order->begin(), end = order->end(); it != end; ++it) {
    if (direction == Forward) {
      for (int i = 0, e = (*it)->getSuccessorCount(); i < e; i++) {
        Vertex * succVertex = (*it)->getSuccessor(i);
        int succWeight = (*it)->getSuccessorWeight(i);

        VertexImpl * succVertexImpl = dynamic_cast<VertexImpl *>(succVertex);
        succVertexImpl->addPredecessorUnsafe(*it, succWeight);
      }
    } else {
      for (int i = 0, e = (*it)->getPredecessorCount(); i < e; i++) {
        Vertex * predVertex = (*it)->getPredecessor(i);
        int predWeight = (*it)->getPredecessorWeight(i);

        VertexImpl * predVertexImpl = dynamic_cast<VertexImpl *>(predVertex);
        predVertexImpl->addSuccessorUnsafe(*it, predWeight);
      }
    }
  }

  // Label the vertices in inverse post-order
  for (auto it = order->rbegin(), end = order->rend(); it != end; ++it) {
    if (direction == Forward)
      (*it)->setOrders(currLabel++, -1);
    else
      (*it)->setOrders(-1, currLabel++);
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
  labelVertices(Forward);

  // Second traversal
  labelVertices(Backward);

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
