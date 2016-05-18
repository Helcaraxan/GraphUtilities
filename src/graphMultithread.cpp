// GraphUtilities library - Directed graph manipulation and querying
// Copyright (C) 2016 - Duco van Amstel
//
// For more license information please see the 'LICENSE' file at the top-level
// of the source code tree.

#include <chrono>
#include <iostream>

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


uint64_t
GraphImpl::getTimestamp() {
  unique_lock<mutex> timestampLock(timestampMutex);

  return chrono::high_resolution_clock::now().time_since_epoch().count();
}

// Worker thread main function

void
queryWorker(GraphImpl * graph, int threadId) {
  Query * query = nullptr;
  CoarsenQueryImpl * cQuery = nullptr;
  PartitionQueryImpl * pQuery = nullptr;
  ReachabilityQueryImpl * rQuery = nullptr;

  while (true) {
    while (!graph->jobQueue.try_pop(query)) {}
      // Do nothing

    if (!query)
      return;

    graph->startWorker();

    while (true) {
      if (!query) {
        graph->stopWorker();
        return;
      }

#if defined(ENABLE_TLS)
      if ((cQuery = dynamic_cast<CoarsenQueryImpl *>(query)))
        graph->processCoarsenQuery(cQuery);
      else if ((pQuery = dynamic_cast<PartitionQueryImpl *>(query)))
        graph->processPartitionQuery(pQuery);
      else if ((rQuery = dynamic_cast<ReachabilityQueryImpl *>(query)))
        graph->processReachabilityQuery(rQuery);
#else // ENABLE_TLS
      if ((cQuery = dynamic_cast<CoarsenQueryImpl *>(query)))
        graph->processCoarsenQuery(cQuery, threadId);
      else if ((pQuery = dynamic_cast<PartitionQueryImpl *>(query)))
        graph->processPartitionQuery(pQuery, threadId);
      else if ((rQuery = dynamic_cast<ReachabilityQueryImpl *>(query)))
        graph->processReachabilityQuery(rQuery, threadId);
#endif // ENABLE_TLS
      else
        cerr << "WARNING: Received a Query of unknown type for processing. It "
          << "was discarded.\n";

      if (!graph->jobQueue.try_pop(query))
        break;
    }

    graph->stopWorker();
  }
}
