// GraphUtilities library - Directed graph manipulation and querying
// Copyright (C) 2016 - Duco van Amstel
//
// For more license information please see the 'LICENSE' file at the top-level
// of the source code tree.

#include <fstream>
#include <iostream>

#include "graph-utilities/implementation/support.hpp"
#include "graph-utilities/implementation/graphImpl.hpp"

using namespace std;


// Coarsening queries

void
#if defined(ENABLE_TLS)
GraphImpl::processCoarsenQuery(CoarsenQueryImpl * query) {
#else // ENABLe_TLS
GraphImpl::processCoarsenQuery(CoarsenQueryImpl * query, int threadId) {
#endif // ENABLE_TLS
  indexGraph();

  startGlobalOperation();

  switch (query->getMethod()) {
    case GreedyCoarsen:
      query->setAnswer(coarsenGreedy(query));
      break;

    case UndefinedCoarsenMethod:
      cerr << "ERROR: Called coarsening method with undefined method" << endl;
      exit(EXIT_FAILURE);
  }

  stopGlobalOperation();
}


// Internal coarsening functions

static void
processVertex(VertexImpl * curr, VertexImpl::Set &readySet,
    VertexImpl::Set &processedSet) {
  VertexImpl::Set::iterator it = readySet.find(curr);

  if (it == readySet.end()) {
    cerr << "ERROR: Tried to process a vertex that is not in the readySet.\n";
    exit(EXIT_FAILURE);
  }

  readySet.erase(it);
  processedSet.insert(curr);
}


bool
GraphImpl::addToReadySet(VertexImpl * curr, VertexImpl::Set &readySet,
    VertexImpl::Set &processedSet) {
  if (readySet.count(curr))
    return true;

  for (auto it = curr->predBegin(), end = curr->predEnd(); it != end; ++it)
    if (!processedSet.count(*it))
      return false;

  readySet.insert(curr);
  return true;
}


GraphImpl *
GraphImpl::coarsenGreedy(CoarsenQueryImpl * query) {
  int progressCount = 0;
  string barTitle = "Greedy coarsening ";
  VertexImpl * curr = nullptr;
  VertexImpl * newVertex = nullptr;
  GraphImpl * coarseGraph = new GraphImpl();
  map<VertexImpl *, int> localMapping;
  VertexImpl::Set readySet;
  VertexImpl::Set vertexGroup;
  VertexImpl::Set processedSet;
  VertexImpl::Queue workQueue;
  VertexImpl::Queue localWorkQueue;
  fstream mappingStream;

#if defined(ENABLE_COARSEN_STATISTICS)
  double totalWorkListSize = 0;
  double workListCount = 0;
  double totalSuccessorCount = 0;
  double illegalSuccessorCount = 0;
#endif // ENABLE_COARSE_STATISTICS

  configureProgressBar(&barTitle, getVertexCount());

  discoverExtremities();

  for (auto it = sources.begin(), end = sources.end(); it != end; ++it) {
    addToReadySet(*it, readySet, processedSet);
    workQueue.push(*it);
  }

  while (!workQueue.empty()) {
    curr = workQueue.front();
    workQueue.pop();

    if (processedSet.count(curr))
      continue;

#if defined(ENABLE_COARSEN_STATISTICS)
    workListCount += 1;
#endif // ENABLE_COARSEN_STATISTICS

    vertexGroup.clear();
    localWorkQueue.push(curr);

    while (!localWorkQueue.empty()) {
      curr = localWorkQueue.front();
      localWorkQueue.pop();

#if defined(ENABLE_COARSEN_STATISTICS)
      totalWorkListSize += 1;
      totalSuccessorCount += curr->getSuccessorCount();
#endif // ENABLE_COARSEN_STATISTICS

      processVertex(curr, readySet, processedSet);
      vertexGroup.insert(curr);

      resultProgressBar(++progressCount);

      for (auto it = curr->succBegin(), end = curr->succEnd();
          it != end; ++it) {
        if (addToReadySet(*it, readySet, processedSet))
          localWorkQueue.push(*it);
#if defined(ENABLE_COARSEN_STATISTICS)
        else
          illegalSuccessorCount += 1;
#endif // ENABLE_COARSEN_STATISTICS
      }

      if (vertexGroup.size() >= (unsigned int) query->getFactor())
        break;
    }

    while (!localWorkQueue.empty()) {
      curr = localWorkQueue.front();
      localWorkQueue.pop();
      workQueue.push(curr);
    }

    int newWeight = 0;
    newVertex = coarseGraph->addVertexUnsafe(0);
    for (auto it = vertexGroup.begin(), end = vertexGroup.end();
        it != end; ++it) {
      query->getMap()[(*it)->getId()] = newVertex->getId();
      newWeight += (*it)->getWeight();

      for (auto predIt = (*it)->predBegin(), predEnd = (*it)->predEnd();
          predIt != predEnd; ++predIt) {
        VertexImpl::Set predecessorSet;

        if (!vertexGroup.count(*predIt)) {
          VertexImpl * localPred =
            coarseGraph->vertices[query->getMap()[(*predIt)->getId()]];

          if (predecessorSet.find(localPred) != predecessorSet.end())
            continue;
          else
            predecessorSet.insert(localPred);

          auto localIt = localMapping.find(localPred);
          if (localIt == localMapping.end())
            localMapping[localPred] = (*it)->getPredecessorWeight(*predIt);
          else
            localIt->second += (*it)->getPredecessorWeight(*predIt);
        }
      }

      for (auto mapIt = localMapping.begin(), mapEnd = localMapping.end();
          mapIt != mapEnd; ++mapIt)
        coarseGraph->addEdgeUnsafe(mapIt->first, newVertex, mapIt->second);

      localMapping.clear();
    }

    newVertex->setWeight(newWeight);
  }

  cout << endl;

#if defined(ENABLE_COARSEN_STATISTICS)
  cout << "Average coarsen ratio: ";
  cout << ((double) getVertexCount()) / (double) coarseGraph->getVertexCount();
  cout << "\nWorklist count: " << workListCount << endl;
  cout << "\nAverage worklist size: " << totalWorkListSize / workListCount;
  cout << "\nAverage successor count: ";
  cout << totalSuccessorCount / totalWorkListSize;
  cout << "\nAverage illegal successor count: ";
  cout << illegalSuccessorCount / totalWorkListSize << endl;
#endif // ENABLE_COARSEN_STATISTICS

  return coarseGraph;
}
