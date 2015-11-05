#include <iostream>

#include "graph-utilities/defs.hpp"
#include "graph-utilities/graph.hpp"
#include "graph-utilities/vertex.hpp"

using namespace std;


// Generic query management

void
Graph::pushQuery(Query * query) {
  static int DFSwin = 0;
  static int BBFSwin = 0;

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

  if (!enabled)
    enableGraph();

  // If a method is specified use it
  if ((query->type != Reachability) ||
      (query->query.reachability->getMethod() != UndefinedSearchMethod)) {
    jobQueue.push(query);
    return;
  }

  // Else select the method automatically. If there is no
  // preferred method yet then find one.
  if (getMethod() == UndefinedSearchMethod) {
    ReachabilityQuery * reachQuery = query->query.reachability;
    ReachabilityQuery * DFSReachQuery =
      new ReachabilityQuery(reachQuery->getSource(),
          reachQuery->getTarget(), DFS);
    ReachabilityQuery * BBFSReachQuery =
      new ReachabilityQuery(reachQuery->getSource(),
          reachQuery->getTarget(), BBFS);
    ReachabilityQuery::InternalQueue * internal =
      new ReachabilityQuery::InternalQueue();

    DFSReachQuery->setInternal(internal);
    BBFSReachQuery->setInternal(internal);

    Query * DFSQuery = new Query(DFSReachQuery);
    Query * BBFSQuery = new Query(BBFSReachQuery);

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
    if (DFSReachQuery->getError() || BBFSReachQuery->getError()) {
      cerr << "ERROR: Could not correctly process BBFS or DFS search." << endl;
      exit(EXIT_FAILURE);
    }

    // Check for coherency between queries
    if (DFSReachQuery->getAnswer() != BBFSReachQuery->getAnswer()) {
      cerr << "ERROR: Incoherent results between BBFS and DFS search." << endl;
      exit(EXIT_FAILURE);
    }

    // Record the winner of the query contest
    int winningIndex = 0;
    if (internal->results[0].second > internal->results[1].second)
      winningIndex = 1;

    if (internal->results[winningIndex].first == DFSReachQuery) {
      DFSwin++;
      query->query.reachability->setAnswer(DFSReachQuery->getAnswer());
      query->query.reachability->setMethod(DFS);

#ifdef ENABLE_STATISTICS
      registerQueryStatistics(DFSReachQuery);
#endif // ENABLE_STATISTICS
    } else {
      BBFSwin++;
      query->query.reachability->setAnswer(BBFSReachQuery->getAnswer());
      query->query.reachability->setMethod(BBFS);

#ifdef ENABLE_STATISTICS
      registerQueryStatistics(BBFSReachQuery);
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
    delete DFSQuery;
    delete BBFSQuery;
    delete internal;
  } else {
    query->query.reachability->setMethod(getMethod());

    jobQueue.push(query);
  }

  return;
}


Query *
Graph::pullResult() {
  Query * result = NULL;

  while (!resultQueue.try_pop(result)) {}
    // Do nothing

  return result;
}


void
Graph::endOfQueries() {
  threadShutdown = true;

  for (int i = 0; i < threadCount; i++)
    jobQueue.push(NULL);
}


