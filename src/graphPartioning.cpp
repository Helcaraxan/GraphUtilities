#include <iostream>
#include <algorithm>

#include "graph-utilities/implementation/graphImpl.hpp"
#include "graph-utilities/implementation/patoh.hpp"
#include "graph-utilities/implementation/support.hpp"

using namespace std;


// File local variables
static int partCount = 0;
static vector<int> unions;
static vector<int> sizes;
static vector<int> sourceMax;
static vector<int> sinkMax;
static queue<PartitionNodeImpl *> workQueue;
static PartitionImpl * part = nullptr;


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


// Partitioning queries

const Partition *
GraphImpl::computeConvexPartition(vector<int>& schedule) const {
  queue<PartitionNodeImpl *> workQueue;
  PartitionImpl * partition = new PartitionImpl(schedule.size());

  // Check the provided schedule's validity
  if (!checkSchedule(schedule)) {
    cerr << "Given schedule is invalid. Can not compute convex partition."
      << endl;
    exit(EXIT_FAILURE);
  }
  
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

    child->id = curr->id;
    id = curr->id;
    for (int i = 1, e = memberCount / 2; i < e; i++)
      id = partition->representants[id];

    IaP<PartitionNodeImpl> nextStart = partition->representants[id];

    partition->representants[id] = child;
    curr->children.push_back(child);

    if ((memberCount / 2) > 1)
      workQueue.push(child);
    
    // Create the second child partition
    child = new PartitionNodeImpl(curr);

    child->id = nextStart;
    id = nextStart;
    while (id.isInteger())
      id = partition->representants[id];

    partition->representants[id] = child;
    curr->children.push_back(child);

    if ((memberCount / 2 + memberCount % 2) > 1)
      workQueue.push(child);

    curr->id = -1;
  }

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


int
GraphImpl::getPartitionCost(const Partition * partition, int memorySize) const {
  int cost = 0;
  queue<PartitionNodeImpl *> workQueue;
  PartitionNodeImpl * curr =
    dynamic_cast<const PartitionImpl *>(partition)->root;

  // Compute the max-live and IO-costs for all nodes
  computePartitionNodeCosts(curr);

  workQueue.push(curr);

  while (!workQueue.empty()) {
    curr = workQueue.front();
    workQueue.pop();

    if (curr->maxLive <= memorySize)
      continue;

    cost += curr->ioCost;

    for (int i = 0, e = curr->children.size(); i < e; i++)
      workQueue.push(curr->children[i]);
  }

  return cost;
}


// Internal partitioning query functions

void
GraphImpl::computePartitionNodeCosts(PartitionNodeImpl * node,
    Vertex::IdSet * group) const {
  vector<Vertex::IdSet *> childSets;

  // In the case of a leaf node just instantiate the right values
  if (node->id != -1) {
    node->maxLive = 1;
    node->ioCost = 0;

    if (group)
      group->insert(node->id);

    return;
  }

  // Recursively compute the values for the child nodes
  for (int i = 0, e = node->children.size(); i < e; i++) {
    childSets.push_back(new Vertex::IdSet());
    computePartitionNodeCosts(node->children[i], childSets.back());
  }

  // Evaluate communications between the child nodes
  for (int i = 0, e = node->children.size(); i < e; i++) {
    if (node->maxLive < node->children[i]->maxLive)
      node->maxLive = node->children[i]->maxLive;
    
    // Account for communications to each other child node seperatly
    vector<int> ioCosts(e - i - 1, 0);
    for (auto it = childSets[i]->begin(), end = childSets[i]->end();
        it != end; ++it) {
      Vertex * curr = getVertex(*it);

      vector<int> localCosts(e - i - 1, 0);
      for (int j = 0, e2 = curr->getSuccessorCount(); j < e2; j++) {
        int target = curr->getSuccessor(j)->getId();

        for (int k = i + 1; k < e; k++) {
          if (childSets[k]->count(target)) {
            localCosts[k - i - 1] = 1;
            break;
          }
        }
      }

      for (int j = 0; j < e - i - 1; j++)
        ioCosts[j] += localCosts[j];
    }

    // Update the node's values
    for (int j = 0; j < e - i - 1; j++) {
      node->ioCost += ioCosts[j];

      if (node->maxLive < ioCosts[j])
        node->maxLive = ioCosts[j];
    }
  }

  // Except for the root node call gather the IDs of the vertices member of the
  // current node's partition
  if (group) {
    for (int i = 0, e = childSets.size(); i < e; i++) {
      group->insert(childSets[i]->begin(), childSets[i]->end());
      delete childSets[i];
    }
  }
}


void
#if defined(ENABLE_TLS)
GraphImpl::processPartitionQuery(PartitionQueryImpl * query) {
#else // ENABLE_TLS
GraphImpl::processPartitionQuery(PartitionQueryImpl * query, int threadId) {
#endif // ENABLE_TLS
  switch (query->getMethod()) {
    case Convexify:
      partitionConvexify(query);
      break;

    case MaxDistance:
      partitionMaxDistance(query);
      break;

    case UndefinedPartitionMethod:
      query->setError(true);
      break;
  }

  resultQueue.push(query);
}


HGraph *
GraphImpl::getHGraph(Vertex::IdSet& idSet) const {
  int idx = 0;
  HGraph * hyperGraph = new HGraph();

  hyperGraph->cellCount = idSet.size();
  hyperGraph->backwardMap.resize(idSet.size());

  for (auto it = idSet.begin(), end = idSet.end(); it != end; ++it) {
    Vertex * curr = getVertex(*it);

    hyperGraph->forwardMap[*it] = idx;
    hyperGraph->backwardMap[idx++] = *it;

    for (int i = 0, e = curr->getSuccessorCount(); i < e; i++) {
      if (idSet.count(curr->getSuccessor(i)->getId())) {
        hyperGraph->netCount++;
        break;
      }
    }
  }

  hyperGraph->xpins.push_back(0);

  for (auto it = idSet.begin(), end = idSet.end(); it != end; ++it) {
    Vertex * curr = getVertex(*it);

    hyperGraph->pins.push_back(hyperGraph->forwardMap[*it]);
    for (int i = 0, e = curr->getSuccessorCount(); i < e; i++) {
      int succ = curr->getSuccessor(i)->getId();

      if (idSet.count(succ))
        hyperGraph->pins.push_back(hyperGraph->forwardMap[succ]);
    }

    hyperGraph->xpins.push_back(hyperGraph->pins.size());
  }

  return hyperGraph;
}


void
GraphImpl::convBisect(PartitionNodeImpl * parent) {
  HGraph * hyperGraph = nullptr;
  Vertex::IdSet idSet;
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


  // We iterate partitioning untill there are no unsafe vertices left
  while (!idSet.empty()) {
    // We assign the current sources and sinks to the two target partitions
    int partId = 0;
    queue<int> checkQueue;

    for (auto it = idSet.begin(), end = idSet.end(); it != end; ++it) {
      int source = hyperGraph->forwardMap[*it];
      Vertex * curr = getVertex(*it);

      predCount[source] = 0;
      for (int i = 0, e = curr->getPredecessorCount(); i < e; i++) {
        if (idSet.count(curr->getPredecessor(i)->getId()))
          predCount[source]++;
      }

      succCount[source] = 0;
      for (int i = 0, e = curr->getSuccessorCount(); i < e; i++) {
        if (idSet.count(curr->getSuccessor(i)->getId()))
          succCount[source]++;
      }

      if (!succCount[source] && !predCount[source]) {
        patohPart[source] = partId;
        convPart[source] = partId;
        partId = 1 - partId;
        checkQueue.push(source);
      } else if (!predCount[source]) {
        patohPart[source] = 0;
        convPart[source] = 0;
        checkQueue.push(source);
      } else if (!succCount[source]) {
        patohPart[source] = 1;
        convPart[source] = 1;
        checkQueue.push(source);
      } else {
        patohPart[source] = -1;
      }
    }

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

    // Check for new safe vertices and update the patohPart vector accordingly
    while (!checkQueue.empty()) {
      int id = checkQueue.front();
      checkQueue.pop();

      // This vertex is safe so remove it from the partition set
      idSet.erase(hyperGraph->backwardMap[id]);

      Vertex * curr = getVertex(hyperGraph->backwardMap[id]);

      if (convPart[id] == 0) {
        for (int i = 0, e = curr->getSuccessorCount(); i < e; i++) {
          int succ = curr->getSuccessor(i)->getId();
          int target = hyperGraph->forwardMap[succ];

          if (idSet.count(succ)) {
            predCount[target]--;

            if ((convPart[target] == 2) && (predCount[target] == 0)) {
              patohPart[target] = 0;
              convPart[target] = 0;
              checkQueue.push(target);
            }
          }
        }
      } else if (convPart[id] == 1) {
        for (int i = 0, e = curr->getPredecessorCount(); i < e; i++) {
          int pred = curr->getPredecessor(i)->getId();
          int target = hyperGraph->forwardMap[pred];

          if (idSet.count(pred)) {
            succCount[target]--;

            if ((convPart[target] == 2) && (succCount[target] == 0)) {
              patohPart[target] = 1;
              convPart[target] = 1;
              checkQueue.push(target);
            }
          }
        }
      } else {
        exit(EXIT_FAILURE);
      }

    }
  }

  // Add the two new partitions to the partioning hierarchy
  int memberA = -1;
  int memberB = -1;
  PartitionNodeImpl * childA = new PartitionNodeImpl(parent);
  PartitionNodeImpl * childB = new PartitionNodeImpl(parent);

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
    } else if (convPart[idx] == 1) { 
      if (memberB == -1)
        part->representants[ptr] = childB;
      else
        part->representants[ptr] = memberB;

      memberB = ptr;
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
    workQueue.push(childA);
  else
    resultProgressBar(++partCount);

  if (part->representants[memberB].isInteger())
    workQueue.push(childB);
  else
    resultProgressBar(++partCount);
}


void
GraphImpl::maxDistBisect(PartitionNodeImpl * parent, Vertex::IdSet& sourceSet,
    Vertex::IdSet& sinkSet) {
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

        subParts[target]->id = *it;
        part->representants[*it] = subParts[target];
      } else {
        part->representants[*it] = mapIt->second->id;
        mapIt->second->id = *it;
      }

      if (sourceSet.count(*it))
        subSources[target]->insert(*it);

      if (sinkSet.count(*it))
        subSinks[target]->insert(*it);
    }

    for (auto mapIt = subParts.begin(), mapEnd = subParts.end();
        mapIt != mapEnd; ++mapIt) {
      if (part->representants[mapIt->second->id].isInteger())
        maxDistBisect(mapIt->second, *subSources[mapIt->first],
            *subSinks[mapIt->first]);
      else
        resultProgressBar(++partCount);

      delete subSources[mapIt->first];
      delete subSinks[mapIt->first];
    }

    parent->id = -1;
  } else {
    /* If the sub-graph is connected then perform a bisection based on the
     * sourceMax and sinkMax values of each node */
    stack<Vertex *> todo;
    Vertex::IdSet subSources;
    Vertex::IdSet subSinks;

    // First reset the sourceMax and sinkMax values for the target set
    for (auto it = idSet.begin(), end = idSet.end(); it != end; ++it) {
      sourceMax[*it] = 0;
      sinkMax[*it] = 0;
    }

    // Perform a DFS to recompute the sourceMax values
    for (auto it = sourceSet.begin(), end = sourceSet.end(); it != end; ++it)
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
    for (auto it = sinkSet.begin(), end = sinkSet.end(); it != end; ++it)
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
            subSources.insert(curr->getId());
            subSinks.insert(pred->getId());
          }
        }
      }
    }

    // Eliminate false positives from the potential new sources and sinks
    auto setIt = subSources.begin(), setEnd = subSources.end();
    while (setIt != setEnd) {
      if (sourceMax[*setIt] < sinkMax[*setIt]) {
        setIt = subSources.erase(setIt);
      } else {
        bool erased = false;
        Vertex * curr = getVertex(*setIt);

        for (int i = 0, e = curr->getPredecessorCount(); i < e; i++) {
          Vertex * pred = curr->getPredecessor(i);

          if (idSet.count(pred->getId())) {
            if (sourceMax[pred->getId()] >= sinkMax[pred->getId()]) {
              setIt = subSources.erase(setIt);
              erased = true;
              break;
            }
          }
        }

        if (!erased)
          setIt++;
      }
    }

    setIt = subSinks.begin();
    setEnd = subSinks.end();
    while (setIt != setEnd) {
      if (sourceMax[*setIt] >= sinkMax[*setIt]) {
        setIt = subSinks.erase(setIt);
      } else {
        bool erased = false;
        Vertex * curr = getVertex(*setIt);

        for (int i = 0, e = curr->getSuccessorCount(); i < e; i++) {
          Vertex * succ = curr->getSuccessor(i);

          if (idSet.count(succ->getId())) {
            if (sourceMax[succ->getId()] < sinkMax[succ->getId()]) {
              setIt = subSinks.erase(setIt);
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
    PartitionNodeImpl * subPartA = new PartitionNodeImpl(parent);
    PartitionNodeImpl * subPartB = new PartitionNodeImpl(parent);
    IaP<PartitionNodeImpl> ptrA = subPartA;
    IaP<PartitionNodeImpl> ptrB = subPartB;

    for (auto it = idSet.begin(), end = idSet.end(); it != end; ++it) {
      if (sourceMax[*it] < sinkMax[*it]) {
        part->representants[*it] = ptrA;
        subPartA->id = *it;
        ptrA = *it;
      } else {
        part->representants[*it] = ptrB;
        subPartB->id = *it;
        ptrB = *it;
      }
    }

    // Recurse on the two new partitions if they are not atomic
    if (part->representants[subPartA->id].isInteger())
      maxDistBisect(subPartA, sourceSet, subSinks);
    else
      resultProgressBar(++partCount);

    sourceSet.clear();
    subSinks.clear();

    if (part->representants[subPartB->id].isInteger())
      maxDistBisect(subPartB, subSources, sinkSet);
    else
      resultProgressBar(++partCount);

    parent->id = -1;
  }
}


void
GraphImpl::partitionConvexify(PartitionQueryImpl * query) {
  IaP<PartitionNodeImpl> ptr = 0;
  string barTitle = "Partitioning ";

  partCount = 0;
  configureProgressBar(&barTitle, vertices.size());
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

  workQueue.push(part->root);
  while (!workQueue.empty()) {
    PartitionNodeImpl * parent = workQueue.front();
    workQueue.pop();

    convBisect(parent);
  }

  cout << "\n";

  query->setPartition(part);
  part = nullptr;
}


void
GraphImpl::partitionMaxDistance(PartitionQueryImpl * query) {
  Vertex::IdSet sourceSet;
  Vertex::IdSet sinkSet;
  IaP<PartitionNodeImpl> ptr = 0;
  string barTitle = "Partitioning ";

  partCount = 0;
  configureProgressBar(&barTitle, vertices.size());
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

  // Construct the root souce and sink sets
  for (auto it = sources.begin(), end = sources.end(); it != end; ++it)
    sourceSet.insert((*it)->getId());

  for (auto it = sinks.begin(), end = sinks.end(); it != end; ++it)
    sinkSet.insert((*it)->getId());

  // Adjust the size of union-find and sourceMax/sinkMax vectors
  sizes.resize(vertices.size());
  unions.resize(vertices.size());
  sourceMax.resize(vertices.size());
  sinkMax.resize(vertices.size());

  // Recursively bisect the graph
  maxDistBisect(part->root, sourceSet, sinkSet);

  cout << "\n";

  // Associate the new partition with the query that requested it
  query->setPartition(part);
  part = nullptr;
}
