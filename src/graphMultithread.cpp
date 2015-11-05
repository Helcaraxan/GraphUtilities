#include "graph-utilities/defs.hpp"
#include "graph-utilities/graph.hpp"
#include "graph-utilities/vertex.hpp"

using namespace std;


// Multi-thread management

void
Graph::startWorker() {
  unique_lock<mutex> queryWaitLock(queryWaitMutex);

  if (noQueries)
    queryWaitCondition.wait(queryWaitLock);

  activeThreads++;
}


void
Graph::stopWorker() {
  unique_lock<mutex> queryWaitLock(queryWaitMutex, defer_lock);
  unique_lock<mutex> globalWaitLock(globalWaitMutex, defer_lock);

  queryWaitLock.lock();
  activeThreads--;
  queryWaitLock.unlock();

  if (noQueries && (activeThreads == 0)) {
    globalWaitLock.lock();
    globalWaitCondition.notify_one();
    globalWaitLock.unlock();
  }
}


void
Graph::startGlobalOperation() {
  unique_lock<mutex> queryWaitLock(queryWaitMutex, defer_lock);
  unique_lock<mutex> globalWaitLock(globalWaitMutex, defer_lock);

  queryWaitLock.lock();
  noQueries = true;
  queryWaitLock.unlock();

  globalWaitLock.lock();
  if (activeThreads > 0)
    globalWaitCondition.wait(globalWaitLock);

  globalWaitLock.unlock();

#ifdef ENABLE_TLS
  Vertex::checkDFSId(vertices.size());
#endif // ENABLE_TLS
}


void
Graph::stopGlobalOperation() {
  unique_lock<mutex> queryWaitLock(queryWaitMutex);

#ifdef ENABLE_TLS
  Vertex::checkDFSId(vertices.size());
#endif // ENABLE_TLS

  noQueries = false;

  queryWaitCondition.notify_all();
}


void
Graph::enableGraph() {
  startGlobalOperation();

  if (!enabled) {
    for (int i = 0; i < threadCount; i++)
      queryThreads.push_back(new thread(queryWorker, this, i));

    enabled = true;
  }

  stopGlobalOperation();
}


void
Graph::disableGraph() {
  if (threadCount > 0) {
    endOfQueries();

    while (!queryThreads.empty()) {
      queryThreads.back()->join();

      delete queryThreads.back();
      queryThreads.pop_back();
    }

    threadCount = 0;
  }
}


// Worker thread main function

void
queryWorker(Graph * graph, int threadId) {
  Query * query = NULL;

  while (true) {
    while (!graph->jobQueue.try_pop(query)) {}
      // Do nothing

    if (graph->threadShutdown)
      return;

    graph->startWorker();

    while (true) {
      switch (query->type) {
#ifdef ENABLE_TLS
        case Reachability:
          graph->processReachabilityQuery(query);
          break;

        case Partition:
          graph->processPartitionQuery(query);
          break;
#else // ENABLE_TLS
        case Reachability:
          graph->processReachabilityQuery(threadId, query);
          break;

        case Partition:
          graph->processPartitionQuery(threadId, query);
          break;
#endif // ENABLE_TLS
      }

      if (!graph->jobQueue.try_pop(query))
        break;
    }

    graph->stopWorker();
  }
}
