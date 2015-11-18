#include <iostream>

#include "graph-utilities/implementation/graphImpl.hpp"

using namespace std;


// Generic query management

void
GraphImpl::pushQuery(Query * query) {
  static int DFSwin = 0;
  static int BBFSwin = 0;
  ReachabilityQueryImpl * rQuery = dynamic_cast<ReachabilityQueryImpl *>(query);

  unique_lock<mutex> methodLock(methodMutex, defer_lock);

  // Verify that the graph has been indexed
  if (!indexed) {
    indexGraph();
    methodLock.lock();
    setMethod(UndefinedSearchMethod);
    DFSwin = 0;
    BBFSwin = 0;
    methodLock.unlock();
  }

  if (queryThreads.size() == 0)
    enableQueries();

  // If a method is specified use it
  if (!rQuery || (rQuery->getMethod() != UndefinedSearchMethod)) {
    jobQueue.push(query);
    return;
  }

  // Else select the method automatically. If there is no
  // preferred method yet then find one.
  if (getMethod() == UndefinedSearchMethod) {
    ReachabilityQueryImpl * DFSQuery =
      new ReachabilityQueryImpl(rQuery->getSourceI(), rQuery->getTargetI(),
          DFS);
    ReachabilityQueryImpl * BBFSQuery =
      new ReachabilityQueryImpl(rQuery->getSourceI(), rQuery->getTargetI(),
          BBFS);
    ReachabilityQueryImpl::InternalQueue * internal =
      new ReachabilityQueryImpl::InternalQueue();

    DFSQuery->setInternal(internal);
    BBFSQuery->setInternal(internal);

    unique_lock<mutex> lock(internal->guard);
    jobQueue.push(DFSQuery);
    jobQueue.push(BBFSQuery);

    // Wait for the results
    internal->signal.wait(lock);

    if (internal->results.size() < 2) {
      cerr << "ERROR: Did not receive results for an internal query" << endl;
      exit(EXIT_FAILURE);
    }

    // Check for errors in the queries
    if (DFSQuery->getError() || BBFSQuery->getError()) {
      cerr << "ERROR: Could not correctly process BBFS or DFS search." << endl;
      exit(EXIT_FAILURE);
    }

    // Check for coherency between queries
    if (DFSQuery->getAnswer() != BBFSQuery->getAnswer()) {
      cerr << "ERROR: Incoherent results between BBFS and DFS search." << endl;
      exit(EXIT_FAILURE);
    }

    // Record the winner of the query contest
    int winningIndex = 0;
    if (internal->results[0].second > internal->results[1].second)
      winningIndex = 1;

    if (internal->results[winningIndex].first == DFSQuery) {
      DFSwin++;
      rQuery->setAnswer(DFSQuery->getAnswer());
      rQuery->setMethod(DFS);

#if defined(ENABLE_STATISTICS)
      registerQueryStatistics(DFSQuery);
#endif // ENABLE_STATISTICS
    } else {
      BBFSwin++;
      rQuery->setAnswer(BBFSQuery->getAnswer());
      rQuery->setMethod(BBFS);

#if defined(ENABLE_STATISTICS)
      registerQueryStatistics(BBFSQuery);
#endif // ENABLE_STATISTICS
    }

    // Push the real query to the results
    resultQueue.push(query);

    // If necessary register the global query contest winner as method
    methodLock.lock();
    if (getMethod() == UndefinedSearchMethod) {
      if (DFSwin >= AUTO_THRESHOLD)
        setMethod(DFS);
      else if (BBFSwin >= AUTO_THRESHOLD)
        setMethod(BBFS);
    }
    methodLock.unlock();

    // Clean-up the internal queries
    delete internal;
    delete DFSQuery;
    delete BBFSQuery;
  } else {
    rQuery->setMethod(getMethod());

    jobQueue.push(query);
  }

  return;
}


Query *
GraphImpl::pullResult() {
  Query * result = NULL;

  while (!resultQueue.try_pop(result)) {}
    // Do nothing

  return result;
}


void
GraphImpl::enableQueries() {
  startGlobalOperation();

  threadShutdown = false;

  for (int i = 0; i < threadCount; i++)
    queryThreads.push_back(new thread(queryWorker, this, i));

  stopGlobalOperation();
}


void
GraphImpl::disableQueries() {
  startGlobalOperation();

  threadShutdown = true;

  for (int i = 0; i < threadCount; i++)
    jobQueue.push(NULL);

  while (threadCount > 0) {
    queryThreads.back()->join();

    delete queryThreads.back();
    queryThreads.pop_back();

    threadCount--;
  }

  stopGlobalOperation();
}
