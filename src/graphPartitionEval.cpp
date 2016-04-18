#include <fstream>
#include <iostream>

#include "graph-utilities/implementation/graphImpl.hpp"
#include "graph-utilities/implementation/patoh.hpp"
#include "graph-utilities/implementation/support.hpp"

using namespace std;


// Partition cost evaluation Interface functions

double
GraphImpl::getPartitionCost(const Partition * partition, int memorySize,
    IOComplexityType type, const char * tileFile) const {
  double cost = 0;
  const PartitionImpl * internalPartition =
    dynamic_cast<const PartitionImpl *>(partition);

  // If not done yet compute the PartitionNode costs
  if (internalPartition->root->maxLive == -1)
    computePartitionNodeCosts(const_cast<PartitionImpl *>(internalPartition));

  // Open the tile dump file if provided
  fstream tileStream;
  if (tileFile)
    tileStream.open(tileFile, ios_base::out);

  // Traverse the partition tree for the IO complexity evaluation
  queue<PartitionNodeImpl *> workQueue;
  PartitionNodeImpl * curr =
    dynamic_cast<const PartitionImpl *>(partition)->root;

  workQueue.push(curr);
  while (!workQueue.empty()) {
    curr = workQueue.front();
    workQueue.pop();

    if (curr->maxLive <= memorySize) {
      if (type == AvgLoadStore)
        cost += curr->exportCost;

      if (tileFile) {
        Vertex::IdSet idSet;

        partition->represents(curr, idSet);
        
        tileStream << "{\n";
        for (auto it = idSet.begin(), end = idSet.end(); it != end; ++it)
          tileStream << "\t" << *it << "\n";

        tileStream << "}" << endl;
      }

      continue;
    }

    cost += curr->ioCost;

    for (int i = 0, e = curr->children.size(); i < e; i++)
      workQueue.push(curr->children[i]);
  }

  switch (type) {
    case TotalLoads:
      return cost;

    case AvgLoadStore:
      cost /= getVertexCount();
      if (cost < (double) 10 / (double) getVertexCount())
        return 0.0;
      else
        return cost;

    default:
      break;
  }

  // Close the tile dump file
  if (tileFile)
    tileStream.close();

  return -1.0;
}


int
GraphImpl::getCutCost(int partitionCount) const {
  int cost = 0;
  HGraph * hyperGraph = nullptr;
  vector<int> patohPart;
  vector<Vertex::IdSet *> idSets;

  // Prepare the ID sets
  for (int i = 0; i < partitionCount; i++)
    idSets.push_back(new Vertex::IdSet());

  // Create the graphs ID set and the associated hyper-graph
  for (int i = 0, e = vertices.size(); i < e; i++) {
    if (vertices[i])
      idSets.front()->insert(i);
  }

  hyperGraph = getHGraph(*idSets.front());
  patohPart.resize(idSets.front()->size(), -1);
  idSets.front()->clear();

  // Perform the partioning with the PaToH library
  int cut = 0;
  vector<int> weights;
  vector<int> partitionWeights(partitionCount);

  if (hyperGraph->cellCount > hyperGraph->netCount)
    weights.resize(hyperGraph->cellCount, 1);
  else
    weights.resize(hyperGraph->netCount, 1);

  PPaToH_Parameters pargs = new PaToH_Parameters;
  PaToH_Initialize_Parameters(pargs, PATOH_CONPART, PATOH_SUGPARAM_DEFAULT);
  pargs->_k = partitionCount;

  PaToH_Alloc(pargs, hyperGraph->cellCount, hyperGraph->netCount, 1,
      weights.data(), weights.data(), hyperGraph->xpins.data(),
      hyperGraph->pins.data());

  PaToH_Part(pargs, hyperGraph->cellCount, hyperGraph->netCount, 1, 1,
      weights.data(), weights.data(), hyperGraph->xpins.data(),
      hyperGraph->pins.data(), nullptr, patohPart.data(),
      partitionWeights.data(), &cut);

  PaToH_Free();

  for (int i = 0, e = patohPart.size(); i < e; i++)
    idSets[patohPart[i]]->insert(hyperGraph->backwardMap[i]);

  // Clean-up memory
  delete hyperGraph;

  // Compute the cut cost
  for (int i = 0; i < partitionCount; i++) {
    for (auto it = idSets[i]->begin(), end = idSets[i]->end();
        it != end; ++it) {
      Vertex * curr = getVertex(*it);
      vector<int> succs(partitionCount, 0);

      for (int j = 0, e = curr->getSuccessorCount(); j < e; j++) {
        int target = curr->getSuccessor(j)->getId();

        for (int k = 0; k < partitionCount; k++) {
          if (idSets[k]->count(target)) {
            succs[k] = 1;
            break;
          }
        }
      }

      for (int j = 1; j < partitionCount; j++) {
        if (j != i)
          cost += succs[j];
      }
    }
  }

  for (auto it = partitions.begin(), end = partitions.end(); it != end; ++it)
    (*it)->clear();

  return cost;
}


int
GraphImpl::getCutCost(const Partition * partition) const {
  int cost = 0;
  int partitionCount = 0;
  vector<Vertex::IdSet *> partitions;

  // Convert the partition scheme
  const PartitionNode * root = partition->getRoot();
  for (int i = 0, e = root->getChildCount(); i < e; i++) {
    const PartitionNode * node = root->getChild(i);

    if (node->getChildCount() > 0) {
      cerr << "ERROR: Cannot compute cut-cost of a non k-way partitioning.\n";
      goto end;
    }

    partitions.push_back(new Vertex::IdSet());
    partition->represents(node, *partitions.back());
  }

  // Compute the cut cost
  partitionCount = partitions.size();
  for (int i = 0; i < partitionCount; i++) {
    for (auto it = partitions[i]->begin(), end = partitions[i]->end();
        it != end; ++it) {
      Vertex * curr = getVertex(*it);
      vector<int> succs(partitionCount, 0);

      for (int j = 0, e = curr->getSuccessorCount(); j < e; j++) {
        int target = curr->getSuccessor(j)->getId();

        for (int k = 0; k < partitionCount; k++) {
          if (partitions[k]->count(target)) {
            succs[k] = 1;
            break;
          }
        }
      }

      for (int j = 1; j < partitionCount; j++) {
        if (j != i)
          cost += succs[j];
      }
    }
  }

end:
  for (auto it = partitions.begin(), end = partitions.end(); it != end; ++it)
    (*it)->clear();

  return cost;
}


int
GraphImpl::getCutCost(const Partition * partition, int partitionCount) const {
  int cost = 0;
  int startIdx = 0, endIdx = 0;
  vector<int> schedule;
  vector<Vertex::IdSet *> partitions;

  // Iteratively create the partitions
  partition->extractSchedule(schedule);
  for (int i = 0; i < partitionCount; i++) {
    partitions.push_back(new Vertex::IdSet());

    endIdx = startIdx + schedule.size() / partitionCount;
    if (endIdx > (int) schedule.size())
      endIdx = schedule.size();

    for (int j = startIdx; j < endIdx; j++)
      partitions.back()->insert(schedule[j]);

    startIdx = endIdx;
  }

  // Compute the cut cost
  for (int i = 0; i < partitionCount; i++) {
    for (auto it = partitions[i]->begin(), end = partitions[i]->end();
        it != end; ++it) {
      Vertex * curr = getVertex(*it);
      vector<int> succs(partitionCount, 0);

      for (int j = 0, e = curr->getSuccessorCount(); j < e; j++) {
        int target = curr->getSuccessor(j)->getId();

        for (int k = 0; k < partitionCount; k++) {
          if (partitions[k]->count(target)) {
            succs[k] = 1;
            break;
          }
        }
      }

      for (int j = 1; j < partitionCount; j++) {
        if (j != i)
          cost += succs[j];
      }
    }
  }

  for (auto it = partitions.begin(), end = partitions.end(); it != end; ++it)
    (*it)->clear();

  return cost;
}


// Internal partition cost evaluation functions

void
GraphImpl::computePartitionNodeCosts(PartitionImpl * partition) const {
  vector<int> schedule;
  stack<PartitionNodeImpl *> workQueue;
  Vertex * curr, * pred, * succ;
  PartitionNodeImpl * sourceNode, * targetNode, * subTree;

  partition->extractSchedule(schedule);

  // Compute the IO and export costs
  for (auto it = schedule.begin(), end = schedule.end(); it != end; ++it) {
    set<PartitionNodeImpl *> targetNodes;
    vector<PartitionNodeImpl *> parents;

    curr = getVertex(*it);
    targetNode = partition->representants[*it];

    for (int i = 0, e = curr->getPredecessorCount(); i < e; i++) {
      pred = curr->getPredecessor(i);
      sourceNode = partition->representants[pred->getId()];

      subTree = partition->getSubTree(sourceNode, targetNode);
      subTree->ioCost++;
    }

    for (int i = 0, e = curr->getSuccessorCount(); i < e; i++) {
      succ = curr->getSuccessor(i);
      sourceNode = partition->representants[succ->getId()];

      subTree = partition->getSubTree(sourceNode, targetNode);
      targetNodes.insert(subTree);
    }

    while (targetNode) {
      parents.push_back(targetNode);
      targetNode = targetNode->parent;
    }

    for (auto it2 = parents.rbegin(), end2 = parents.rend();
        it2 != end2; ++it2) {
      if (targetNodes.count(*it2)) {
        (*it2)->exportCost++;
        break;
      }
    }
  }

  workQueue.push(partition->root);
  while (!workQueue.empty()) {
    sourceNode = workQueue.top();

    if (sourceNode->id != -1) {
      workQueue.pop();

      sourceNode->maxLive = 1;

      if (sourceNode->parent)
        sourceNode->parent->evaluatedChildren++;
    } else if (sourceNode->evaluatedChildren == 0) {
      for (int i = 0, e = sourceNode->children.size(); i < e; i++)
        workQueue.push(sourceNode->children[i]);
    } else {
      workQueue.pop();

      if (sourceNode->evaluatedChildren != (int) sourceNode->children.size()) {
        cerr << "ERROR: Did not evaluate all children for a PartitionNode.\n";
        exit(EXIT_FAILURE);
      }

      sourceNode->maxLive = sourceNode->ioCost;
      for (int i = 0, e = sourceNode->children.size(); i < e; i++) {
        if (sourceNode->children[i]->maxLive > sourceNode->maxLive)
          sourceNode->maxLive = sourceNode->children[i]->maxLive;
      }

      if (sourceNode->parent)
        sourceNode->parent->evaluatedChildren++;
    }
  }
}
