#include <iostream>

#include "graph-utilities/implementation/graphImpl.hpp"
#include "graph-utilities/implementation/patoh.hpp"
#include "graph-utilities/implementation/support.hpp"

using namespace std;


// Global variables (declared in main file - graphPartion.cpp)

extern atomic<int> partCount;
extern PartitionImpl * part;
extern tbb::concurrent_queue<PartitionTask *> taskQueue;


// Internal Partition (Convexify method) implementation functions

void
GraphImpl::convBisect(PartitionNodeImpl * parent) const {
  int cut;
  int safeCount = 0;
  int partitionWeights[2];
  bool forceProgress = false;
  HGraph * hyperGraph = nullptr;
  Vertex::IdSet idSet;
  vector<int> weights;
  vector<int> convPart;
  vector<int> patohPart;
  vector<int> predCount;
  vector<int> succCount;

  // Create the target idSet and the associated hyper-graph
  IaP<PartitionNodeImpl> ptr = parent->id;
  while (ptr.isInteger()) {
    idSet.insert(ptr);
    ptr = part->representants[ptr];
  }

  hyperGraph = getHGraph(idSet);
  convPart.resize(idSet.size(), 2);
  patohPart.resize(idSet.size(), -1);
  predCount.resize(idSet.size(), 0);
  succCount.resize(idSet.size(), 0);

  if (hyperGraph->cellCount > hyperGraph->netCount)
    weights.resize(hyperGraph->cellCount, 1);
  else
    weights.resize(hyperGraph->netCount, 1);

  // Initialize the successor and predecessor counts
  for (auto it = idSet.begin(), end = idSet.end(); it != end; ++it) {
    int id = hyperGraph->forwardMap[*it];
    Vertex * curr = getVertex(*it);

    for (int i = 0, e = curr->getPredecessorCount(); i < e; i++) {
      if (idSet.count(curr->getPredecessor(i)->getId()))
        predCount[id]++;
    }

    for (int i = 0, e = curr->getSuccessorCount(); i < e; i++) {
      if (idSet.count(curr->getSuccessor(i)->getId()))
        succCount[id]++;
    }
  }

  // We iterate partitioning untill there are no unsafe vertices left
  while (!idSet.empty()) {
    // We process the current sources and sinks of the unsafe vertices
    int partId = 0;
    queue<int> checkQueue;

    for (auto it = idSet.begin(), end = idSet.end(); it != end; ++it) {
      int id = hyperGraph->forwardMap[*it];

      if (!succCount[id] && !predCount[id]) {
        if (forceProgress) {
          patohPart[id] = partId;
          convPart[id] = partId;
          partId = 1 - partId;
        }
        checkQueue.push(id);
      } else if (!predCount[id]) {
        if (forceProgress) {
          patohPart[id] = 0;
          convPart[id] = 0;
        }
        checkQueue.push(id);
      } else if (!succCount[id]) {
        if (forceProgress) {
          patohPart[id] = 1;
          convPart[id] = 1;
        }
        checkQueue.push(id);
      } else {
        patohPart[id] = -1;
      }

      if (!forceProgress)
        patohPart[id] = -1;
    }

    // Perform the partioning with the prefixed vertex assignements
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

    // Check for new safe vertices and update the partition vectors accordingly
    int newSafeCount = 0;
    while (!checkQueue.empty()) {
      int id = checkQueue.front();
      checkQueue.pop();

      Vertex * curr = getVertex(hyperGraph->backwardMap[id]);

      if ((patohPart[id] == 0) && (predCount[id] == 0)) {
        idSet.erase(hyperGraph->backwardMap[id]);
        convPart[id] = 0;
        newSafeCount++;

        for (int i = 0, e = curr->getSuccessorCount(); i < e; i++) {
          int succ = curr->getSuccessor(i)->getId();
          int target = hyperGraph->forwardMap[succ];

          if (idSet.count(succ)) {
            predCount[target]--;

            if (predCount[target] == 0)
              checkQueue.push(target);
          }
        }
      } else if ((patohPart[id] == 1) && (succCount[id] == 0)) {
        idSet.erase(hyperGraph->backwardMap[id]);
        convPart[id] = 1;
        newSafeCount++;

        for (int i = 0, e = curr->getPredecessorCount(); i < e; i++) {
          int pred = curr->getPredecessor(i)->getId();
          int target = hyperGraph->forwardMap[pred];

          if (idSet.count(pred)) {
            succCount[target]--;

            if (succCount[target] == 0)
              checkQueue.push(target);
          }
        }
      }
    }

    if (safeCount == newSafeCount) {
      forceProgress = true;
    } else {
      forceProgress = false;
      safeCount = newSafeCount;
    }
  }

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

    if (convPart[idx] == 0) {
      if (memberA == -1)
        part->representants[ptr] = childA;
      else
        part->representants[ptr] = memberA;

      memberA = ptr;
      childA->size++;
    } else if (convPart[idx] == 1) { 
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
GraphImpl::partitionConvexify(PartitionQueryImpl * query) const {
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
  part->method = Convexify;
  query->setPartition(part);
  part = nullptr;
}
