// GraphUtilities library - Directed graph manipulation and querying
// Copyright (C) 2016 - Duco van Amstel
//
// For more license information please see the 'LICENSE' file at the top-level
// of the source code tree.

#include <iostream>

#include "graph-utilities/implementation/graphImpl.hpp"
#include "graph-utilities/implementation/patoh.hpp"
#include "graph-utilities/implementation/support.hpp"

using namespace std;


// Global variables (declared in main file - graphPartion.cpp)

extern atomic<int> partCount;
extern PartitionImpl * part;
extern tbb::concurrent_queue<PartitionTask *> taskQueue;


// Internal Partition (PaToH method) implementation functions

void
GraphImpl::patohBisect(PartitionNodeImpl * parent) const {
  HGraph * hyperGraph = nullptr;
  Vertex::IdSet idSet;
  vector<int> patohPart;

  // Create the target idSet and the associated hyper-graph
  IaP<PartitionNodeImpl> ptr = parent->id;
  while (ptr.isInteger()) {
    idSet.insert(ptr);
    ptr = part->representants[ptr];
  }

  hyperGraph = getHGraph(idSet);
  patohPart.resize(idSet.size(), -1);

  // Perform the partioning with the prefixed vertex assignements
  int cut;
  int partitionWeights[2];
  vector<int> weights;

  if (hyperGraph->cellCount > hyperGraph->netCount)
    weights.resize(hyperGraph->cellCount, 1);
  else
    weights.resize(hyperGraph->netCount, 1);

  PPaToH_Parameters pargs = new PaToH_Parameters;
  PaToH_Initialize_Parameters(pargs, PATOH_CONPART, PATOH_SUGPARAM_DEFAULT);
  pargs->_k = 2;

  PaToH_Alloc(pargs, hyperGraph->cellCount, hyperGraph->netCount, 1,
      weights.data(), weights.data(), hyperGraph->xpins.data(),
      hyperGraph->pins.data());
  
  PaToH_Part(pargs, hyperGraph->cellCount, hyperGraph->netCount, 1, 1,
      weights.data(), weights.data(), hyperGraph->xpins.data(),
      hyperGraph->pins.data(), nullptr, patohPart.data(), partitionWeights,
      &cut);

  PaToH_Free();

  // Add the two new partitions to the partioning hierarchy
  int memberA = -1;
  int memberB = -1;
  PartitionNodeImpl * childA = new PartitionNodeImpl(parent);
  PartitionNodeImpl * childB = new PartitionNodeImpl(parent);
  part->nodeCount++;
  part->nodeCount++;

  ptr = parent->id;
  do {
    int idx = hyperGraph->forwardMap[ptr];
    IaP<PartitionNodeImpl> newPtr = part->representants[ptr];

    if (patohPart[idx] == 0) {
      if (memberA == -1)
        part->representants[ptr] = childA;
      else
        part->representants[ptr] = memberA;

      memberA = ptr;
      childA->size++;
    } else if (patohPart[idx] == 1) { 
      if (memberB == -1)
        part->representants[ptr] = childB;
      else
        part->representants[ptr] = memberB;

      memberB = ptr;
      childB->size++;
    } else {
      cerr << "ERROR: Found node not in partition after convex bisection.\n";
      exit(EXIT_FAILURE);
    }

    ptr = newPtr;
  } while (ptr.isInteger());

  childA->id = memberA;
  childB->id = memberB;
  parent->id = -1;

  if (part->representants[memberA].isInteger())
    taskQueue.push(new PartitionTask(childA, nullptr, nullptr));
  else
    partCount++;

  if (part->representants[memberB].isInteger())
    taskQueue.push(new PartitionTask(childB, nullptr, nullptr));
  else
    partCount++;
}


void
GraphImpl::partitionPaToH(PartitionQueryImpl * query) const {
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

  // Perform the partitioning
  taskQueue.push(new PartitionTask(part->root, nullptr, nullptr));
  partitionWorker(this);
  cout << "\n";

  // Associate the new partition with the query that requested it
  part->method = PaToH;
  query->setPartition(part);
  part = nullptr;
}
