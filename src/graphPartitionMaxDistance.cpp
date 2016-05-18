// GraphUtilities library - Directed graph manipulation and querying
// Copyright (C) 2016 - Duco van Amstel
//
// For more license information please see the 'LICENSE' file at the top-level
// of the source code tree.

#include <iostream>

#include <tbb/concurrent_queue.h>

#include "graph-utilities/implementation/graphImpl.hpp"
#include "graph-utilities/implementation/support.hpp"

using namespace std;


// Global variables (declared in main file - graphPartition.cpp)

extern bool taskPoison;
extern atomic<int> partCount;
extern PartitionImpl * part;
extern tbb::concurrent_queue<PartitionTask *> taskQueue;


// File local variables

vector<int> unions;
vector<int> sizes;
vector<int> sourceMax;
vector<int> sinkMax;


// File local functions

static int
findUnion(int id) {
  while (id != unions[id])
    id = unions[id];

  return id;
}


static void
mergeUnion(int i, int j) {
  i = findUnion(i);
  j = findUnion(j);

  if (i == j)
    return;
  
  if (sizes[i] < sizes[j]) {
    unions[i] = j;
    sizes[j] += sizes[i];
  } else {
    unions[j] = i;
    sizes[i] += sizes[j];
  }
}


// Internal Partition (MaxDistance method) implementation functions

void
GraphImpl::maxDistBisect(PartitionNodeImpl * parent, Vertex::IdSet * sourceSet,
    Vertex::IdSet * sinkSet) const {
  /* Create the target set of vertex IDs */
  Vertex::IdSet idSet;
  IaP<PartitionNodeImpl> id = parent->id;

  while (id.isInteger()) {
    idSet.insert(id);
    id = part->representants[id];
  }

  /* Check that the current target set for bisection is actually a connected
   * sub-graph through union-find */
  for (auto it = idSet.begin(), end = idSet.end(); it != end; ++it) {
    unions[*it] = *it;
    sizes[*it] = 1;
  }

  for (auto it = idSet.begin(), end = idSet.end(); it != end; ++it) {
    Vertex * curr = getVertex(*it);

    for (int i = 0, e = curr->getSuccessorCount(); i < e; i++) {
      int succ = curr->getSuccessor(i)->getId();

      if (idSet.count(succ))
        mergeUnion(*it, succ);
    }
  }

  int target = findUnion(*idSet.begin());
  if (sizes[target] != (int) idSet.size()) {
    /* If the sub-graph is not connected then perform a trivial multi-section
     * into the connected components */
    map<int, PartitionNodeImpl *> subParts;
    map<int, Vertex::IdSet *> subSources;
    map<int, Vertex::IdSet *> subSinks;

    for (auto it = idSet.begin(), end = idSet.end(); it != end; ++it) {
      target = findUnion(*it);

      auto mapIt = subParts.find(target);
      if (mapIt == subParts.end()) {
        subParts[target] = new PartitionNodeImpl(parent);
        subSources[target] = new Vertex::IdSet();
        subSinks[target] = new Vertex::IdSet();
        part->nodeCount++;

        subParts[target]->id = *it;
        part->representants[*it] = subParts[target];

        subParts[target]->size++;
      } else {
        part->representants[*it] = mapIt->second->id;
        mapIt->second->id = *it;

        mapIt->second->size++;
      }

      if (sourceSet->count(*it))
        subSources[target]->insert(*it);

      if (sinkSet->count(*it))
        subSinks[target]->insert(*it);
    }

    delete sourceSet;
    delete sinkSet;

    parent->id = -1;

    for (auto mapIt = subParts.begin(), mapEnd = subParts.end();
        mapIt != mapEnd; ++mapIt) {
      if (part->representants[mapIt->second->id].isInteger())
        taskQueue.push(new PartitionTask(mapIt->second,
            subSources[mapIt->first], subSinks[mapIt->first]));
      else
        partCount++;
    }
  } else {
    /* If the sub-graph is connected then perform a bisection based on the
     * sourceMax and sinkMax values of each node */
    stack<Vertex *> todo;
    Vertex::IdSet * subSources = new Vertex::IdSet();
    Vertex::IdSet * subSinks = new Vertex::IdSet();

    // First reset the sourceMax and sinkMax values for the target set
    for (auto it = idSet.begin(), end = idSet.end(); it != end; ++it) {
      sourceMax[*it] = 0;
      sinkMax[*it] = 0;
    }

    // Perform a DFS to recompute the sourceMax values
    for (auto it = sourceSet->begin(), end = sourceSet->end(); it != end; ++it)
      todo.push(getVertex(*it));

    while (!todo.empty()) {
      Vertex * curr = todo.top();
      todo.pop();

      for (int i = 0, e = curr->getSuccessorCount(); i < e; i++) {
        Vertex * succ = curr->getSuccessor(i);

        if (idSet.count(succ->getId())) {
          if (sourceMax[succ->getId()] < sourceMax[curr->getId()] + 1) {
            sourceMax[succ->getId()] = sourceMax[curr->getId()] + 1;
            todo.push(succ);
          }
        }
      }
    }

    // Perform a DFS to recompute the sinkMax values and select potential new
    // sources and sinks
    for (auto it = sinkSet->begin(), end = sinkSet->end(); it != end; ++it)
      todo.push(getVertex(*it));

    while (!todo.empty()) {
      Vertex * curr = todo.top();
      todo.pop();

      for (int i = 0, e = curr->getPredecessorCount(); i < e; i++) {
        Vertex * pred = curr->getPredecessor(i);

        if (idSet.count(pred->getId())) {
          if (sinkMax[pred->getId()] < sinkMax[curr->getId()] + 1) {
            sinkMax[pred->getId()] = sinkMax[pred->getId()] + 1;
            todo.push(pred);
          }

          if ((sourceMax[curr->getId()] >= sinkMax[curr->getId()]) &&
              (sourceMax[pred->getId()] < sinkMax[pred->getId()])) {
            subSources->insert(curr->getId());
            subSinks->insert(pred->getId());
          }
        }
      }
    }

    // Eliminate false positives from the potential new sources and sinks
    auto setIt = subSources->begin();
    auto setEnd = subSources->end();
    while (setIt != setEnd) {
      if (sourceMax[*setIt] < sinkMax[*setIt]) {
        setIt = subSources->erase(setIt);
      } else {
        bool erased = false;
        Vertex * curr = getVertex(*setIt);

        for (int i = 0, e = curr->getPredecessorCount(); i < e; i++) {
          Vertex * pred = curr->getPredecessor(i);

          if (idSet.count(pred->getId())) {
            if (sourceMax[pred->getId()] >= sinkMax[pred->getId()]) {
              setIt = subSources->erase(setIt);
              erased = true;
              break;
            }
          }
        }

        if (!erased)
          setIt++;
      }
    }

    setIt = subSinks->begin();
    setEnd = subSinks->end();
    while (setIt != setEnd) {
      if (sourceMax[*setIt] >= sinkMax[*setIt]) {
        setIt = subSinks->erase(setIt);
      } else {
        bool erased = false;
        Vertex * curr = getVertex(*setIt);

        for (int i = 0, e = curr->getSuccessorCount(); i < e; i++) {
          Vertex * succ = curr->getSuccessor(i);

          if (idSet.count(succ->getId())) {
            if (sourceMax[succ->getId()] < sinkMax[succ->getId()]) {
              setIt = subSinks->erase(setIt);
              erased = true;
              break;
            }
          }
        }

        if (!erased)
          setIt++;
      }
    }
    
    // Divide the target set into two partitions
    PartitionNodeImpl * childA = new PartitionNodeImpl(parent);
    PartitionNodeImpl * childB = new PartitionNodeImpl(parent);
    IaP<PartitionNodeImpl> ptrA = childA;
    IaP<PartitionNodeImpl> ptrB = childB;
    part->nodeCount++;
    part->nodeCount++;

    for (auto it = idSet.begin(), end = idSet.end(); it != end; ++it) {
      if (sourceMax[*it] < sinkMax[*it]) {
        part->representants[*it] = ptrA;
        childA->id = *it;
        ptrA = *it;

        childA->size++;
      } else {
        part->representants[*it] = ptrB;
        childB->id = *it;
        ptrB = *it;

        childB->size++;
      }
    }

    parent->id = -1;

    // Recurse on the two new partitions if they are not atomic
    if (part->representants[childA->id].isInteger())
      taskQueue.push(new PartitionTask(childA, sourceSet, subSinks));
    else
      partCount++;

    if (part->representants[childB->id].isInteger())
      taskQueue.push(new PartitionTask(childB, subSources, sinkSet));
    else
      partCount++;
  }
}


void
GraphImpl::partitionMaxDistance(PartitionQueryImpl * query) const {
  vector<thread *> workers;
  Vertex::IdSet * sourceSet = new Vertex::IdSet();
  Vertex::IdSet * sinkSet = new Vertex::IdSet();
  IaP<PartitionNodeImpl> ptr = 0;
  string barTitle = "Partitioning ";

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

  // Construct the root source and sink sets
  for (auto it = sources.begin(), end = sources.end(); it != end; ++it)
    sourceSet->insert((*it)->getId());

  for (auto it = sinks.begin(), end = sinks.end(); it != end; ++it)
    sinkSet->insert((*it)->getId());

  // Adjust the size of union-find and sourceMax/sinkMax vectors
  sizes.resize(vertices.size());
  unions.resize(vertices.size());
  sourceMax.resize(vertices.size());
  sinkMax.resize(vertices.size());

  // Create the required number of worker threads and perform partitioning
  taskPoison = false;
  taskQueue.push(new PartitionTask(part->root, sourceSet, sinkSet));
  for (int i = 0, e = query->getThreadCount(); i < e; i++)
    workers.push_back(new thread(partitionWorker, this));

  // Monitor the partition count actively
  do {
    resultProgressBar(partCount);
    this_thread::sleep_for(chrono::milliseconds(100));
  } while (partCount < (int) getVertexCount());

  resultProgressBar(partCount);
  cout << "\n";

  // Terminate workers when partitioning is finished
  taskPoison = true;
  while (!workers.empty()) {
    workers.back()->join();

    delete workers.back();
    workers.pop_back();
  }

  // Associate the new partition with the query that requested it
  part->method = MaxDistance;
  query->setPartition(part);
  part = nullptr;
}
