#include "graph-utilities/defs.hpp"
#include "graph-utilities/graph.hpp"
#include "graph-utilities/vertex.hpp"

using namespace std;


// Local perf measurement function

#ifdef __amd64
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
Graph::getMethod() {
  return searchMethod.load(memory_order_acquire);
}


void
Graph::setMethod(SearchMethod method) {
  searchMethod.store(method, memory_order_release);
}


// Internal reachability query functions

void
#ifdef ENABLE_TLS
Graph::processReachabilityQuery(Query * query) {
  Vertex::checkDFSId(vertices.size());

  switch (query->query.reachability->getMethod()) {
    case DFS:
      areConnectedDFS(query->query.reachability);
      break;

    case BBFS:
      areConnectedBBFS(query->query.reachability);
      break;

    case NoLabels:
      areConnectedNoLabels(query->query.reachability);
      break;
#else // ENABLE_TLS
Graph::processReachabilityQuery(int threadId, Query * query) {
  switch (query->query.reachability->getMethod()) {
    case DFS:
      areConnectedDFS(query->query.reachability, threadId);
      break;

    case BBFS:
      areConnectedBBFS(query->query.reachability, threadId);
      break;

    case NoLabels:
      areConnectedNoLabels(query->query.reachability, threadId);
      break;
#endif // ENABLE_TLS

    case UndefinedSearchMethod:
      query->query.reachability->setError(true);
      break;
  }

  if (query->query.reachability->getInternal() == NULL) {
#ifdef ENABLE_STATISTICS
    registerQueryStatistics(query->query.reachability);
#endif // ENABLE_STATISTICS

    resultQueue.push(query);
  }
}


/* Use the previously done indexation to answer to the query */
void
#ifdef ENABLE_TLS
Graph::areConnectedDFS(ReachabilityQuery * query) {
#else // ENABLE_TLS
Graph::areConnectedDFS(ReachabilityQuery * query, int threadId) {
#endif // ENABLE_TLS
  uint64_t cycleCount = 0;
  uint64_t timestamp = 0;
  Vertex * curr;
  vector<Vertex *> searchStack;

#ifdef ENABLE_BENCHMARKS
  unique_lock<mutex> benchmarkLock(benchmarkMutex, defer_lock);
  cycleCount = getClock();
#else // ENABLE_BENCHMARKS
  if (query->getInternal() != NULL)
    cycleCount = getClock();
#endif // ENABLE_BENCHMARKS

  query->setAnswer(false);

  // Are U and V the same vertex?
  if (query->getSource() == query->getTarget()) {
#ifdef ENABLE_STATISTICS
    query->searchedNodes++;
    query->path.push_back(query->getSource());
#endif // ENABLE_STATISTICS
    query->setAnswer(true);
    goto end;
  }

  // Can V be a descendant of U in the standard graph or U a descendant of V in
  // the reverse graph?
  if ((query->getSource()->orderLabel > query->getTarget()->orderLabel) ||
      (query->getTarget()->revOrderLabel > query->getSource()->revOrderLabel))
    goto end;

  // Do a DFS on the subgraph specified by both orders to get the final answer
  timestamp = chrono::high_resolution_clock::now().time_since_epoch() / chrono::nanoseconds(1);

  searchStack.push_back(query->getSource());

  while (!searchStack.empty()) {
    if (query->getCancel()) {
      query->setError(true);
      goto cancel;
    }

    curr = searchStack.back();

#ifdef ENABLE_STATISTICS
    if (!query->path.empty() && (curr == query->path.back())) {
      query->path.pop_back();
      searchStack.pop_back();
      continue;
    }

    query->searchedNodes++;
    query->path.push_back(curr);
#else // ENABLE_STATISTICS
    searchStack.pop_back();
#endif // ENABLE_STATISTICS

    for (auto it = curr->successors.begin(), end = curr->successors.end(); it != end; ++it) {
      if ((indexMethod & 0x04) && ((*it)->orderLabel > query->getTarget()->orderLabel))
        break;

      if (*it == query->getTarget()) {
#ifdef ENABLE_STATISTICS
        query->searchedNodes++;
        query->path.push_back(query->getTarget());
#endif // ENABLE_STATISTICS

        query->setAnswer(true);
        goto end;
      }

#ifdef ENABLE_TLS
      if ((*it)->getDFSId() != timestamp) {
        (*it)->setDFSId(timestamp);
#else // ENABLE_TLS
      if ((*it)->getDFSId(threadId) != timestamp) {
        (*it)->setDFSId(threadId, timestamp);
#endif // ENABLE_TLS
        if (((*it)->orderLabel < query->getTarget()->orderLabel) &&
            ((*it)->revOrderLabel > query->getTarget()->revOrderLabel))
          searchStack.push_back(*it);
      }
    }
  }

end:
#ifdef ENABLE_BENCHMARKS
  cycleCount = getClock() - cycleCount;
  benchmarkLock.lock();
  cyclesSpentQuerying += cycleCount;
  queryNumber++;
  benchmarkLock.unlock();

  if (query->getInternal() != NULL) {
#else // ENABLE_BENCHMARKS
  if (query->getInternal() != NULL) {
    cycleCount = getClock() - cycleCount;
#endif // ENABLE_BENCHMARKS

    // Set the query time in the case of an internal query
    unique_lock<mutex> lock(query->getInternal()->guard);

    pair<ReachabilityQuery *, uint64_t> result(query, cycleCount);
    query->getInternal()->results.push_back(result);

    if (query->getInternal()->results.size() == 2)
      query->getInternal()->signal.notify_all();
  }

  return;

cancel:
#ifdef ENABLE_BENCHMARKS
  cycleCount = getClock() - cycleCount;

  if (query->getInternal() != NULL) {
#else // ENABLE_BENCHMARKS
  if (query->getInternal() != NULL) {
    cycleCount = getClock() - cycleCount;
#endif // ENABLE_BENCHMARKS

    // Set the query time in the case of an internal query
    unique_lock<mutex> lock(query->getInternal()->guard);

    pair<ReachabilityQuery *, uint64_t> result(query, cycleCount);
    query->getInternal()->results.push_back(result);

    if (query->getInternal()->results.size() == 2)
      query->getInternal()->signal.notify_all();
  }

  query->setError(true);

  return;
}


void
#ifdef ENABLE_TLS
Graph::areConnectedBBFS(ReachabilityQuery * query) {
#else // ENABLE_TLS
Graph::areConnectedBBFS(ReachabilityQuery * query, int threadId) {
#endif
  uint64_t cycleCount = 0;
  uint64_t forwardId = 0, backwardId = 0;
  Vertex * curr;
  queue<Vertex *> searchQueueForward;
  queue<Vertex *> searchQueueBackward;

#ifdef ENABLE_BENCHMARKS
  unique_lock<mutex> benchmarkLock(benchmarkMutex, defer_lock);
  cycleCount = getClock();
#else // ENABLE_BENCHMARKS
  if (query->getInternal() == NULL)
    cycleCount = getClock();
#endif // ENABLE_BENCHMARKS

  query->setAnswer(false);

  // Are U and V the same vertex?
  if (query->getSource() == query->getTarget()) {
#ifdef ENABLE_STATISTICS
    query->searchedNodes++;
#endif // ENABLE_STATISTICS
    query->setAnswer(true);
    goto end;
  }

  // Can V be a descendant of U in the standard graph or U a descendant of V in
  // the reverse graph?
  if ((query->getSource()->orderLabel > query->getTarget()->orderLabel) ||
      (query->getTarget()->revOrderLabel > query->getSource()->revOrderLabel))
    goto end;

  // Do a BBFS on the subgraph specified by both orders to get the final answer
  forwardId = chrono::high_resolution_clock::now().time_since_epoch().count();
  backwardId = forwardId + 1;

#ifdef ENABLE_TLS
  query->getSource()->setDFSId(forwardId);
  query->getTarget()->setDFSId(backwardId);
#else // ENABLE_TLS
  query->getSource()->setDFSId(threadId, forwardId);
  query->getTarget()->setDFSId(threadId, backwardId);
#endif // ENABLE_TLS

  searchQueueForward.push(query->getSource());
  searchQueueBackward.push(query->getTarget());

  while (!searchQueueForward.empty() || !searchQueueBackward.empty()) {
    if (query->getCancel()) {
      query->setError(true);
      goto cancel;
    }

    if (!searchQueueForward.empty()) {
      curr = searchQueueForward.front();
      searchQueueForward.pop();

#ifdef ENABLE_STATISTICS
      query->searchedNodes++;
#endif // ENABLE_STATISTICS

      for (auto it = curr->succ_begin(), end = curr->succ_end();
          it != end; ++it) {
        if ((indexMethod & 0x04) &&
            ((*it)->orderLabel > query->getTarget()->orderLabel))
          break;

#ifdef ENABLE_TLS
        if ((*it)->getDFSId() == backwardId) {
#else // ENABLE_TLS
        if ((*it)->getDFSId(threadId) == backwardId) {
#endif // ENABLE_TLS
          query->setAnswer(true);
          goto end;
        }

#ifdef ENABLE_TLS
        if ((*it)->getDFSId() != forwardId) {
          (*it)->setDFSId(forwardId);
#else // ENABLE_TLS
        if ((*it)->getDFSId(threadId) != forwardId) {
          (*it)->setDFSId(threadId, forwardId);
#endif
          if (((*it)->orderLabel < query->getTarget()->orderLabel) &&
              ((*it)->revOrderLabel > query->getTarget()->revOrderLabel))
            searchQueueForward.push(*it);
        }
      }
    }

    if (!searchQueueBackward.empty()) {
      curr = searchQueueBackward.front();
      searchQueueBackward.pop();

#ifdef ENABLE_STATISTICS
      query->searchedNodes++;
#endif // ENABLE_STATISTICS

      for (auto it = curr->pred_begin(), end = curr->pred_end();
          it != end; ++it) {
        if ((indexMethod & 0x04) &&
            ((*it)->orderLabel < query->getSource()->orderLabel))
          break;

#ifdef ENABLE_TLS
        if ((*it)->getDFSId() == forwardId) {
#else // ENABLE_TLS
        if ((*it)->getDFSId(threadId) == forwardId) {
#endif // ENABLE_TLS
          query->setAnswer(true);
          goto end;
        }

#ifdef ENABLE_TLS
        if ((*it)->getDFSId() != backwardId) {
          (*it)->setDFSId(backwardId);
#else // ENABLE_TLS
        if ((*it)->getDFSId(threadId) != backwardId) {
          (*it)->setDFSId(threadId, backwardId);
#endif // ENABLE_TLS
          if (((*it)->orderLabel > query->getSource()->orderLabel) &&
              ((*it)->revOrderLabel < query->getSource()->revOrderLabel))
            searchQueueBackward.push(*it);
        }
      }
    }
  }

end:
#ifdef ENABLE_BENCHMARKS
  cycleCount = getClock() - cycleCount;
  benchmarkLock.lock();
  cyclesSpentQuerying += cycleCount;
  queryNumber++;
  benchmarkLock.unlock();

  // Set the query time
  if (query->getInternal() != NULL) {
#else // ENABLE_BENCHMARKS
  if (query->getInternal() != NULL) {
    cycleCount = getClock() - cycleCount;
#endif // ENABLE_BENCHMARKS

    // Set the query time in the case of an internal query
    unique_lock<mutex> lock(query->getInternal()->guard);

    pair<ReachabilityQuery *, uint64_t> result(query, cycleCount);
    query->getInternal()->results.push_back(result);

    if (query->getInternal()->results.size() == 2)
      query->getInternal()->signal.notify_all();
  }

  return;

cancel:
#ifdef ENABLE_BENCHMARKS
  cycleCount = getClock() - cycleCount;
  if (query->getInternal() != NULL) {
#else // ENABLE_BENCHMARKS
  if (query->getInternal() != NULL) {
    cycleCount = getClock() - cycleCount;
#endif // ENABLE_BENCHMARKS

    // Set the query time in the case of an internal query
    unique_lock<mutex> lock(query->getInternal()->guard);

    pair<ReachabilityQuery *, uint64_t> result(query, cycleCount);
    query->getInternal()->results.push_back(result);

    if (query->getInternal()->results.size() == 2)
      query->getInternal()->signal.notify_all();
  }

  query->setError(true);

  return;
}


void
#ifdef ENABLE_TLS
Graph::areConnectedNoLabels(ReachabilityQuery * query) {
#else // ENABLE_TLS
Graph::areConnectedNoLabels(ReachabilityQuery * query, int threadId) {
#endif // ENABLE_TLS
  uint64_t timestamp;
  stack<Vertex *> toVisit;
  Vertex * curr;

  timestamp = chrono::high_resolution_clock::now().time_since_epoch().count();

  query->setAnswer(false);

  toVisit.push(query->getSource());
#ifdef ENABLE_TLS
  query->getSource()->setDFSId(timestamp);
#else // ENABLE_TLS
  query->getSource()->setDFSId(threadId, timestamp);
#endif // ENABLE_TLS

  while (!toVisit.empty()) {
    if (query->getCancel()) {
      query->setError(true);
      return;
    }

    curr = toVisit.top();
    toVisit.pop();

    for (auto it = curr->succ_begin(), end = curr->succ_end();
        it != end; ++it) {
      if (*it == query->getSource()) {
        query->setAnswer(true);
        return;
      }

#ifdef ENABLE_TLS
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
