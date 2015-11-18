#include "graph-utilities/implementation/graphImpl.hpp"
#include "graph-utilities/implementation/queriesImpl.hpp"

using namespace std;


// Multi-thread management

int
GraphImpl::getThreadCount() const {
  return threadCount;
}


void
GraphImpl::startWorker() {
  unique_lock<mutex> queryWaitLock(queryWaitMutex);

  if (noQueries)
    queryWaitCondition.wait(queryWaitLock);

  activeThreads++;
}


void
GraphImpl::stopWorker() {
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
GraphImpl::startGlobalOperation() {
  unique_lock<mutex> queryWaitLock(queryWaitMutex, defer_lock);
  unique_lock<mutex> globalWaitLock(globalWaitMutex, defer_lock);

  queryWaitLock.lock();
  noQueries = true;
  queryWaitLock.unlock();

  globalWaitLock.lock();
  if (activeThreads > 0)
    globalWaitCondition.wait(globalWaitLock);

  globalWaitLock.unlock();

#if defined(ENABLE_TLS)
  VertexImpl::checkDFSId(vertices.size());
#endif // ENABLE_TLS
}


void
GraphImpl::stopGlobalOperation() {
  unique_lock<mutex> queryWaitLock(queryWaitMutex);

#if defined(ENABLE_TLS)
  VertexImpl::checkDFSId(vertices.size());
#endif // ENABLE_TLS

  noQueries = false;

  queryWaitCondition.notify_all();
}


// Worker thread main function

void
queryWorker(GraphImpl * graph, int threadId) {
  Query * query = NULL;

  while (true) {
    while (!graph->jobQueue.try_pop(query)) {}
      // Do nothing

    if (graph->threadShutdown)
      return;

    graph->startWorker();

    while (true) {
      PartitionQueryImpl * pQuery = dynamic_cast<PartitionQueryImpl *>(query);
      ReachabilityQueryImpl * rQuery =
        dynamic_cast<ReachabilityQueryImpl *>(query);

#if defined(ENABLE_TLS)
      if (pQuery)
        graph->processPartitionQuery(pQuery);
      else if (rQuery)
        graph->processReachabilityQuery(rQuery);
#else // ENABLE_TLS
      if (pQuery)
        graph->processPartitionQuery(pQuery, threadId);
      else if (rQuery)
        graph->processReachabilityQuery(rQuery, threadId);
#endif // ENABLE_TLS

      if (!graph->jobQueue.try_pop(query))
        break;
    }

    graph->stopWorker();
  }
}
