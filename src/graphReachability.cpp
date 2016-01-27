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


// Reachability query configuration

SearchMethod
GraphImpl::getMethod() const {
  return searchMethod.load(memory_order_acquire);
}


void
GraphImpl::setMethod(SearchMethod method) {
  searchMethod.store(method, memory_order_release);
}


// Internal reachability query functions

void
#if defined(ENABLE_TLS)
GraphImpl::processReachabilityQuery(ReachabilityQueryImpl * query) {
  VertexImpl::checkDFSId(vertices.size());

  switch (query->getMethod()) {
    case DFS:
      areConnectedDFS(query);
      break;

    case BBFS:
      areConnectedBBFS(query);
      break;

    case NoLabels:
      areConnectedNoLabels(query);
      break;
#else // ENABLE_TLS
GraphImpl::processReachabilityQuery(ReachabilityQueryImpl * query,
    int threadId) {
  switch (query->getMethod()) {
    case DFS:
      areConnectedDFS(query, threadId);
      break;

    case BBFS:
      areConnectedBBFS(query, threadId);
      break;

    case NoLabels:
      areConnectedNoLabels(query, threadId);
      break;
#endif // ENABLE_TLS

    case UndefinedSearchMethod:
      query->setError(true);
      break;
  }

  if (query->getInternal() == nullptr) {
#if defined(ENABLE_STATISTICS)
    registerQueryStatistics(query);
#endif // ENABLE_STATISTICS

    resultQueue.push(query);
  }
}


/* Use the previously done indexation to answer to the query */
void
#if defined(ENABLE_TLS)
GraphImpl::areConnectedDFS(ReachabilityQueryImpl * query) {
#else // ENABLE_TLS
GraphImpl::areConnectedDFS(ReachabilityQueryImpl * query, int threadId) {
#endif // ENABLE_TLS
  uint64_t cycleCount = 0;
  uint64_t timestamp = 0;
  VertexImpl * curr;
  VertexImpl::Array searchArray;
  pair<int, int> sourceOrders;
  pair<int, int> targetOrders;

#if defined(ENABLE_BENCHMARKS)
  unique_lock<mutex> counterLock(counterMutex, defer_lock);
  cycleCount = getClock();
#else // ENABLE_BENCHMARKS
  if (query->getInternal() != nullptr)
    cycleCount = getClock();
#endif // ENABLE_BENCHMARKS

  query->setAnswer(false);

  // Are U and V the same vertex?
  if (query->getSource() == query->getTarget()) {
#if defined(ENABLE_STATISTICS)
    query->addSearchedNodes();
    query->getPath().push_back(query->getSourceI());
#endif // ENABLE_STATISTICS
    query->setAnswer(true);
    goto end;
  }

  // Can V be a descendant of U in the standard graph or U a descendant of V in
  // the reverse graph?
  sourceOrders = query->getSourceI()->getOrders();
  targetOrders = query->getTargetI()->getOrders();

  if ((sourceOrders.first > targetOrders.first) ||
      (sourceOrders.second < targetOrders.second))
    goto end;

  // Do a DFS on the subgraph specified by both orders to get the final answer
  timestamp = chrono::high_resolution_clock::now().time_since_epoch().count();

  searchArray.push_back(query->getSourceI());

  while (!searchArray.empty()) {
    if (query->getCancel()) {
      query->setError(true);
      goto cancel;
    }

    curr = searchArray.back();

#if defined(ENABLE_STATISTICS)
    if (!query->getPath().empty() && (curr == query->getPath().back())) {
      query->getPath().pop_back();
      searchArray.pop_back();
      continue;
    }

    query->addSearchedNodes();
    query->getPath().push_back(curr);
#else // ENABLE_STATISTICS
    searchArray.pop_back();
#endif // ENABLE_STATISTICS

    targetOrders = query->getTargetI()->getOrders();
    for (auto it = curr->succBegin(), end = curr->succEnd();
        it != end; ++it) {
      sourceOrders = (*it)->getOrders();

      if ((indexMethod & 0x04) && (sourceOrders.first > targetOrders.first))
        break;

      if (*it == query->getTarget()) {
#if defined(ENABLE_STATISTICS)
        query->addSearchedNodes();
        query->getPath().push_back(query->getTargetI());
#endif // ENABLE_STATISTICS

        query->setAnswer(true);
        goto end;
      }

#if defined(ENABLE_TLS)
      if ((*it)->getDFSId() != timestamp) {
        (*it)->setDFSId(timestamp);
#else // ENABLE_TLS
      if ((*it)->getDFSId(threadId) != timestamp) {
        (*it)->setDFSId(threadId, timestamp);
#endif // ENABLE_TLS
        if ((sourceOrders.first < targetOrders.first) &&
            (sourceOrders.second > targetOrders.second))
          searchArray.push_back(*it);
      }
    }
  }

end:
#if defined(ENABLE_BENCHMARKS)
  cycleCount = getClock() - cycleCount;
  counterLock.lock();
  cyclesSpentQuerying += cycleCount;
  counterLock.unlock();

  if (query->getInternal() != nullptr) {
#else // ENABLE_BENCHMARKS
  if (query->getInternal() != nullptr) {
    cycleCount = getClock() - cycleCount;
#endif // ENABLE_BENCHMARKS

    // Set the query time in the case of an internal query
    unique_lock<mutex> lock(query->getInternal()->guard);

    pair<ReachabilityQueryImpl *, uint64_t> result(query, cycleCount);
    query->getInternal()->results.push_back(result);

    if (query->getInternal()->results.size() == 2)
      query->getInternal()->signal.notify_all();
  }

#if defined(ENABLE_STATISTICS) || defined(ENABLE_BENCHMARKS)
  queryCount++;
#endif //  ENABLE_STATISTICS || ENABLE_BENCHMARKS

  return;

cancel:
#if defined(ENABLE_BENCHMARKS)
  cycleCount = getClock() - cycleCount;

  if (query->getInternal() != nullptr) {
#else // ENABLE_BENCHMARKS
  if (query->getInternal() != nullptr) {
    cycleCount = getClock() - cycleCount;
#endif // ENABLE_BENCHMARKS

    // Set the query time in the case of an internal query
    unique_lock<mutex> lock(query->getInternal()->guard);

    pair<ReachabilityQueryImpl *, uint64_t> result(query, cycleCount);
    query->getInternal()->results.push_back(result);

    if (query->getInternal()->results.size() == 2)
      query->getInternal()->signal.notify_all();
  }

  query->setError(true);

  return;
}


void
#if defined(ENABLE_TLS)
GraphImpl::areConnectedBBFS(ReachabilityQueryImpl * query) {
#else // ENABLE_TLS
GraphImpl::areConnectedBBFS(ReachabilityQueryImpl * query, int threadId) {
#endif
  uint64_t cycleCount = 0;
  uint64_t forwardId = 0, backwardId = 0;
  VertexImpl * curr;
  VertexImpl::Queue searchQueueForward;
  VertexImpl::Queue searchQueueBackward;
  pair<int, int> sourceOrders;
  pair<int, int> targetOrders;

#if defined(ENABLE_BENCHMARKS)
  unique_lock<mutex> counterLock(counterMutex, defer_lock);
  cycleCount = getClock();
#else // ENABLE_BENCHMARKS
  if (query->getInternal() == nullptr)
    cycleCount = getClock();
#endif // ENABLE_BENCHMARKS

  query->setAnswer(false);

  // Are U and V the same vertex?
  if (query->getSource() == query->getTarget()) {
#if defined(ENABLE_STATISTICS)
    query->addSearchedNodes();
#endif // ENABLE_STATISTICS
    query->setAnswer(true);
    goto end;
  }

  // Can V be a descendant of U in the standard graph or U a descendant of V in
  // the reverse graph?
  sourceOrders = query->getSourceI()->getOrders();
  targetOrders = query->getTargetI()->getOrders();

  if ((sourceOrders.first > targetOrders.first) ||
      (sourceOrders.second < targetOrders.second))
    goto end;

  // Do a BBFS on the subgraph specified by both orders to get the final answer
  forwardId = chrono::high_resolution_clock::now().time_since_epoch().count();
  backwardId = forwardId + 1;

#if defined(ENABLE_TLS)
  query->getSourceI()->setDFSId(forwardId);
  query->getTargetI()->setDFSId(backwardId);
#else // ENABLE_TLS
  query->getSourceI()->setDFSId(threadId, forwardId);
  query->getTargetI()->setDFSId(threadId, backwardId);
#endif // ENABLE_TLS

  searchQueueForward.push(query->getSourceI());
  searchQueueBackward.push(query->getTargetI());

  while (!searchQueueForward.empty() || !searchQueueBackward.empty()) {
    if (query->getCancel()) {
      query->setError(true);
      goto cancel;
    }

    if (!searchQueueForward.empty()) {
      curr = searchQueueForward.front();
      searchQueueForward.pop();

#if defined(ENABLE_STATISTICS)
      query->addSearchedNodes();
#endif // ENABLE_STATISTICS

      targetOrders = query->getTargetI()->getOrders();
      for (auto it = curr->succBegin(), end = curr->succEnd();
          it != end; ++it) {
        sourceOrders = (*it)->getOrders();

        if ((indexMethod & 0x04) &&
            (sourceOrders.first > targetOrders.first))
          break;

#if defined(ENABLE_TLS)
        if ((*it)->getDFSId() == backwardId) {
#else // ENABLE_TLS
        if ((*it)->getDFSId(threadId) == backwardId) {
#endif // ENABLE_TLS
          query->setAnswer(true);
          goto end;
        }

#if defined(ENABLE_TLS)
        if ((*it)->getDFSId() != forwardId) {
          (*it)->setDFSId(forwardId);
#else // ENABLE_TLS
        if ((*it)->getDFSId(threadId) != forwardId) {
          (*it)->setDFSId(threadId, forwardId);
#endif
          if ((sourceOrders.first < targetOrders.first) &&
              (targetOrders.second > targetOrders.second))
            searchQueueForward.push(*it);
        }
      }
    }

    if (!searchQueueBackward.empty()) {
      curr = searchQueueBackward.front();
      searchQueueBackward.pop();

#if defined(ENABLE_STATISTICS)
      query->addSearchedNodes();
#endif // ENABLE_STATISTICS

      sourceOrders = query->getSourceI()->getOrders();
      for (auto it = curr->predBegin(), end = curr->predEnd();
          it != end; ++it) {
        targetOrders = (*it)->getOrders();

        if ((indexMethod & 0x04) &&
            (sourceOrders.first > targetOrders.first ))
          break;

#if defined(ENABLE_TLS)
        if ((*it)->getDFSId() == forwardId) {
#else // ENABLE_TLS
        if ((*it)->getDFSId(threadId) == forwardId) {
#endif // ENABLE_TLS
          query->setAnswer(true);
          goto end;
        }

#if defined(ENABLE_TLS)
        if ((*it)->getDFSId() != backwardId) {
          (*it)->setDFSId(backwardId);
#else // ENABLE_TLS
        if ((*it)->getDFSId(threadId) != backwardId) {
          (*it)->setDFSId(threadId, backwardId);
#endif // ENABLE_TLS
          if ((sourceOrders.first < targetOrders.first) &&
              (sourceOrders.second > targetOrders.second))
            searchQueueBackward.push(*it);
        }
      }
    }
  }

end:
#if defined(ENABLE_BENCHMARKS)
  cycleCount = getClock() - cycleCount;
  counterLock.lock();
  cyclesSpentQuerying += cycleCount;
  counterLock.unlock();

  // Set the query time
  if (query->getInternal() != nullptr) {
#else // ENABLE_BENCHMARKS
  if (query->getInternal() != nullptr) {
    cycleCount = getClock() - cycleCount;
#endif // ENABLE_BENCHMARKS

    // Set the query time in the case of an internal query
    unique_lock<mutex> lock(query->getInternal()->guard);

    pair<ReachabilityQueryImpl *, uint64_t> result(query, cycleCount);
    query->getInternal()->results.push_back(result);

    if (query->getInternal()->results.size() == 2)
      query->getInternal()->signal.notify_all();
  }

#if defined(ENABLE_STATISTICS) || defined(ENABLE_BENCHMARKS)
  queryCount++;
#endif //  ENABLE_STATISTICS || ENABLE_BENCHMARKS

  return;

cancel:
#if defined(ENABLE_BENCHMARKS)
  cycleCount = getClock() - cycleCount;
  if (query->getInternal() != nullptr) {
#else // ENABLE_BENCHMARKS
  if (query->getInternal() != nullptr) {
    cycleCount = getClock() - cycleCount;
#endif // ENABLE_BENCHMARKS

    // Set the query time in the case of an internal query
    unique_lock<mutex> lock(query->getInternal()->guard);

    pair<ReachabilityQueryImpl *, uint64_t> result(query, cycleCount);
    query->getInternal()->results.push_back(result);

    if (query->getInternal()->results.size() == 2)
      query->getInternal()->signal.notify_all();
  }

  query->setError(true);

  return;
}


void
#if defined(ENABLE_TLS)
GraphImpl::areConnectedNoLabels(ReachabilityQueryImpl * query) {
#else // ENABLE_TLS
GraphImpl::areConnectedNoLabels(ReachabilityQueryImpl * query, int threadId) {
#endif // ENABLE_TLS
  uint64_t timestamp;
  VertexImpl * curr;
  VertexImpl::Stack toVisit;

  timestamp = chrono::high_resolution_clock::now().time_since_epoch().count();

  query->setAnswer(false);

  toVisit.push(query->getSourceI());
#if defined(ENABLE_TLS)
  query->getSourceI()->setDFSId(timestamp);
#else // ENABLE_TLS
  query->getSourceI()->setDFSId(threadId, timestamp);
#endif // ENABLE_TLS

  while (!toVisit.empty()) {
    if (query->getCancel()) {
      query->setError(true);
      return;
    }

    curr = toVisit.top();
    toVisit.pop();

    for (auto it = curr->succBegin(), end = curr->succEnd();
        it != end; ++it) {
      if (*it == query->getSourceI()) {
        query->setAnswer(true);
        return;
      }

#if defined(ENABLE_TLS)
      if ((*it)->getDFSId() != timestamp) {
        (*it)->setDFSId(timestamp);
#else // ENABLE_TLS
      if ((*it)->getDFSId(threadId) != timestamp) {
        (*it)->setDFSId(threadId, timestamp);
#endif // ENABLE_TLS
        toVisit.push(*it);
      }
    }
  }
}
