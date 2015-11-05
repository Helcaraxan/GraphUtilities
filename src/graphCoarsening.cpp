#include <fstream>
#include <iostream>

#include "graph-utilities/defs.hpp"
#include "graph-utilities/graph.hpp"
#include "graph-utilities/vertex.hpp"
#include "graph-utilities/support.hpp"

using namespace std;


// Coarsening queries

Graph *
Graph::coarsen(CoarsenMethod method, int factor, int secondaryFactor,
    map<int, int> &vertexMap) {
  Graph * coarsenedGraph = NULL;

  indexGraph();

  startGlobalOperation();

  switch (method) {
    case Greedy:
      coarsenedGraph = coarsenGreedy(factor, vertexMap);
      break;

    case UndefinedCoarsenMethod:
      cerr << "ERROR: Called coarsening method with undefined method" << endl;
      exit(EXIT_FAILURE);
  }

  stopGlobalOperation();

  return coarsenedGraph;
}


// Internal coarsening functions

static void
processVertex(Vertex * curr,
    set<Vertex *> &readySet, set<Vertex *> &processedSet) {
  set<Vertex *>::iterator it = readySet.find(curr);

  if (it == readySet.end()) {
    cerr << "ERROR: Tried to process a vertex that is not in the readySet.\n";
    exit(EXIT_FAILURE);
  }

  readySet.erase(it);
  processedSet.insert(curr);
}


bool
Graph::addToReadySet(Vertex * curr,
    set<Vertex *> &readySet, set<Vertex *> &processedSet) {
  if (readySet.count(curr))
    return true;

  for (auto it = curr->predecessors.begin(), end = curr->predecessors.end();
      it != end; ++it) {
    if (!processedSet.count(*it))
      return false;
  }

  readySet.insert(curr);
  return true;
}


Graph *
Graph::coarsenGreedy(int factor, map<int, int> &vertexMap) {
  int progressCount = 0;
  string barTitle = "Greedy coarsening ";
  Vertex * curr = NULL;
  Vertex * newVertex = NULL;
  Graph * coarseGraph = new Graph(false);
  map<Vertex *, int> localMapping;
  set<Vertex *> readySet;
  set<Vertex *> vertexGroup;
  set<Vertex *> processedSet;
  queue<Vertex *> workQueue;
  queue<Vertex *> localWorkQueue;
  fstream mappingStream;

#ifdef ENABLE_COARSEN_STATISTICS
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

#ifdef ENABLE_COARSEN_STATISTICS
    workListCount += 1;
#endif // ENABLE_COARSEN_STATISTICS

    vertexGroup.clear();
    localWorkQueue.push(curr);

    while (!localWorkQueue.empty()) {
      curr = localWorkQueue.front();
      localWorkQueue.pop();

#ifdef ENABLE_COARSEN_STATISTICS
      totalWorkListSize += 1;
      totalSuccessorCount += curr->getSuccessorCount();
#endif // ENABLE_COARSEN_STATISTICS

      processVertex(curr, readySet, processedSet);
      vertexGroup.insert(curr);

      resultProgressBar(++progressCount);

      for (auto it = curr->succ_begin(), end = curr->succ_end();
          it != end; ++it) {
        if (addToReadySet(*it, readySet, processedSet))
          localWorkQueue.push(*it);
#ifdef ENABLE_COARSEN_STATISTICS
        else
          illegalSuccessorCount += 1;
#endif // ENABLE_COARSEN_STATISTICS
      }

      if (vertexGroup.size() >= (unsigned int) factor)
        break;
    }

    while (!localWorkQueue.empty()) {
      curr = localWorkQueue.front();
      localWorkQueue.pop();
      workQueue.push(curr);
    }

    newVertex = coarseGraph->addVertex();
    for (auto it = vertexGroup.begin(), end = vertexGroup.end();
        it != end; ++it) {
      vertexMap[(*it)->id] = newVertex->id;
      newVertex->weight += (*it)->weight;

      for (auto predIt = (*it)->pred_begin(), predEnd = (*it)->pred_end();
          predIt != predEnd; ++predIt) {
        set<Vertex *> predecessorSet;

        if (!vertexGroup.count(*predIt)) {
          Vertex * localPred = coarseGraph->vertices[vertexMap[(*predIt)->id]];

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
  }

  cout << endl;

#ifdef ENABLE_COARSEN_STATISTICS
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
