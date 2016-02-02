#include <iostream>

#include "graph-utilities/implementation/graphImpl.hpp"
#include "graph-utilities/implementation/support.hpp"

using namespace std;


// Global variables (declared in main file - graphPartion.cpp)

extern atomic<int> partCount;
extern PartitionImpl * part;


// Internal Partition (from scheduling) implementation functions

const Partition *
GraphImpl::computeConvexPartition(vector<int>& schedule) const {
  string barTitle = "Partitioning "; 
  queue<PartitionNodeImpl *> workQueue;
  PartitionImpl * partition = new PartitionImpl(vertices.size());

  // Check the provided schedule's validity
  if (!checkSchedule(schedule)) {
    cerr << "Given schedule is invalid. Can not compute convex partition."
      << endl;
    exit(EXIT_FAILURE);
  }

  partCount = 0;
  configureProgressBar(&barTitle, getVertexCount());
  resultProgressBar(partCount);
  
  // Set up the partition structure
  partition->root->id = schedule.back();
  partition->representants[schedule.front()] = partition->root;
  for (int i = 1, e = schedule.size(); i < e; i++)
    partition->representants[schedule[i]] = schedule[i - 1];

  // Iteratively split partitions in two
  workQueue.push(partition->root);
  while (!workQueue.empty()) {
    int memberCount = 0;
    PartitionNodeImpl * child = nullptr;
    PartitionNodeImpl * curr = workQueue.front();
    workQueue.pop();

    // Count number of vertices in the current partition
    IaP<PartitionNodeImpl> id = curr->id;
    do {
      memberCount++;
      id = partition->representants[id];
    } while (id.isInteger());
    
    // Create the first child partition
    child = new PartitionNodeImpl(curr);
    partition->nodeCount++;

    child->size = memberCount / 2;
    child->id = curr->id;
    id = curr->id;
    for (int i = 1, e = memberCount / 2; i < e; i++)
      id = partition->representants[id];

    IaP<PartitionNodeImpl> nextStart = partition->representants[id];

    partition->representants[id] = child;

    if (child->size > 1)
      workQueue.push(child);
    else
      resultProgressBar(++partCount);
    
    // Create the second child partition
    child = new PartitionNodeImpl(curr);
    partition->nodeCount++;

    child->size = memberCount / 2 + memberCount % 2;
    child->id = nextStart;
    id = nextStart;
    while (id.isInteger())
      id = partition->representants[id];

    if (child->size > 1)
      workQueue.push(child);
    else
      resultProgressBar(++partCount);

    curr->id = -1;
  }

  cout << "\n";

  return partition;
}


bool
GraphImpl::checkSchedule(vector<int>& schedule) const {
  set<int> idSet;

  if (schedule.size() != getVertexCount()) {
    cerr << "WARNING: The supplied schedule is not valid. It does not contain "
      << "the same number of vertices as the graph.\n";
    return false;
  }

  for (int i = 0, e = schedule.size(); i < e; i++) {
    Vertex * curr = getVertex(schedule[i]);

    if (!curr) {
      cerr << "WARNING: The supplied schedule contains an ID that does not "
        << "correspond to an existing vertex.\n";
      return false;
    }

    if (idSet.count(schedule[i])) {
      cerr << "WARNING: The supplied schedule contains two times the same ID ("
        << schedule[i] << ").\n";
      return false;
    }

    for (auto j = 0, e2 = curr->getSuccessorCount(); j < e2; j++) {
      if (idSet.count(curr->getSuccessor(j)->getId())) {
        cerr << "WARNING: The supplied schedule is illegal. ID ("
          << curr->getSuccessor(j)->getId() << ") should be scheduled after ID "
          << "(" << schedule[i] << ").\n";
        return false;
      }
    }

    idSet.insert(schedule[i]);
  }

  return true;
}
