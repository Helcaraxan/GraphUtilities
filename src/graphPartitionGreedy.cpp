// GraphUtilities library - Directed graph manipulation and querying
// Copyright (C) 2016 - Duco van Amstel
//
// For more license information please see the 'LICENSE' file at the top-level
// of the source code tree.

#include <set>
#include <list>
#include <vector>
#include <iostream>

#include "graph-utilities/implementation/graphImpl.hpp"
#include "graph-utilities/implementation/support.hpp"

using namespace std;


// Global state for a Greedy partitioning

static PartitionImpl * part = nullptr;
static set<VertexImpl *> * liveSet = nullptr;
static set<VertexImpl *> readyVertices;
static double readyBalance = 1.0;
static list<VertexImpl *> readySuccessors;
static list<VertexImpl *> readyNeighbours;
static double successorsTaken = 0;
static double neighboursTaken = 0;
static vector<bool> processedVertices;


static void
updateReadyVertices(VertexImpl * lastVertex) {
  for (auto it = lastVertex->succBegin(), end = lastVertex->succEnd();
      it != end; ++it) {
    bool ready = true;

    for (auto it2 = (*it)->predBegin(), end2 = (*it)->predEnd();
        it2 != end2; ++it2) {
      if (!processedVertices[(*it2)->getId()]) {
        ready = false;
        break;
      }
    }

    if (ready)
      readyVertices.insert(*it);
  }
}


static bool
updateLiveSet(VertexImpl * nextVertex, int memorySize) {
  set<VertexImpl *> * newLiveSet = new set<VertexImpl *>(*liveSet);

  for (auto it = nextVertex->predBegin(), end = nextVertex->predEnd();
      it != end; ++it) {
    bool addVertex = false;

    for (auto it2 = (*it)->succBegin(), end2 = (*it)->succEnd();
        it2 != end2; ++it2) {
      if (!processedVertices[(*it2)->getId()]) {
        addVertex = true;
        break;
      }
    }

    if (addVertex)
      newLiveSet->insert(*it);
    else
      newLiveSet->erase(*it);
  }

  if (newLiveSet->size() > (unsigned int) memorySize) {
    delete newLiveSet;
    return false;
  }

  delete liveSet;
  liveSet = newLiveSet;
  return true;
}


static VertexImpl *
selectNextVertex(VertexImpl * lastVertex) {
  VertexImpl * next = nullptr;
  set<VertexImpl *> newNeighbours;
  set<VertexImpl *> newSuccessors;

  for (auto it = lastVertex->succBegin(), end = lastVertex->succEnd();
      it != end; ++it) {
    if (readyVertices.count(*it))
      newSuccessors.insert(*it);

    for (auto it2 = (*it)->predBegin(), end2 = (*it)->predEnd();
      it2 != end2; ++it2) {
      if (readyVertices.count(*it2))
        newNeighbours.insert(*it2);
    }
  }

  for (auto it = readySuccessors.begin(), end = readySuccessors.end();
      it != end; ++it) {
    if (*it == lastVertex) {
      it = readySuccessors.erase(it);
      it--;
    } else {
      newSuccessors.erase(*it);
    }
  }

  for (auto it = readyNeighbours.begin(), end = readyNeighbours.end();
      it != end; ++it) {
    if (*it == lastVertex) {
      it = readyNeighbours.erase(it);
      it--;
    } else {
      newNeighbours.erase(*it);
    }
  }

  for (auto it = newSuccessors.begin(), end = newSuccessors.end();
      it != end; ++it)
    readySuccessors.push_back(*it);

  for (auto it = newNeighbours.begin(), end = newNeighbours.end();
      it != end; ++it)
    readyNeighbours.push_back(*it);

  if ((neighboursTaken < successorsTaken * readyBalance) &&
      !readyNeighbours.empty()) {
    next = readyNeighbours.front();
    readyNeighbours.pop_front();
    neighboursTaken += 1;
  } else if (!readySuccessors.empty()) {
    next = readySuccessors.front();
    readySuccessors.pop_front();
    successorsTaken += 1;
  } else {
    next = *readyVertices.begin();
  }

  return next;
}


void
GraphImpl::partitionGreedy(PartitionQueryImpl * query) const {
  int partCount = 0;
  IaP<PartitionNodeImpl> ptr = 0;
  string barTitle = "Partitioning ";

  // Prepare the global state and the progress bar
  processedVertices.resize(vertices.size(), false);

  partCount = 0;
  configureProgressBar(&barTitle, getVertexCount());
  resultProgressBar(partCount);

  // Construct the target Partition instance
  part = new PartitionImpl(vertices.size());
  ptr = part->root;
  for (int i = 0, e = vertices.size(); i < e; i++) {
    if (vertices[i]) {
      part->representants[i] = ptr;
      part->root->id = i;
      ptr = i;
    } else {
      part->representants[i] = -1;
    }
  }

  // Construct the initial readyQueue
  for (auto it = sources.begin(), end = sources.end(); it != end; ++it)
    readyVertices.insert(*it);

  // Construct the partition scheme
  while (!readyVertices.empty()) {
    // Start a new partition
    PartitionNodeImpl * node = new PartitionNodeImpl(part->root);
    VertexImpl * nextVertex = *readyVertices.begin();

    ptr = node;

    // Reinitialize the global state
    readySuccessors.clear();
    readyNeighbours.clear();
    successorsTaken = 0;
    neighboursTaken = 0;

    if (liveSet)
      liveSet->clear();
    else
      liveSet = new set<VertexImpl *>();

    // Iteratively add vertices to the partition
    while (!readyVertices.empty()
        && updateLiveSet(nextVertex, query->getMemorySize())) {
      part->representants[nextVertex->getId()] = ptr;
      node->id = nextVertex->getId();
      ptr = nextVertex->getId();

      processedVertices[nextVertex->getId()] = true;

      readyVertices.erase(nextVertex);

      updateReadyVertices(nextVertex);

      nextVertex = selectNextVertex(nextVertex);

      resultProgressBar(partCount++);
    }
  }

  resultProgressBar(partCount++);
  cout << "\n";

  part->method = Greedy;
  query->setPartition(part);
  part = nullptr;
}
