#include <iostream>

#include "graph-utilities/implementation/graphImpl.hpp"
#include "graph-utilities/implementation/patoh.hpp"

using namespace std;


// File local functions

static vector<int> *
patohBisect(HGraph * hyperGraph) {
  int cut;
  int * ones = NULL;
  int * tpart = NULL;
  int pweight[2] = { 0, 0};
  float tweight[2] = { .5, .5};

  PPaToH_Parameters pargs = new PaToH_Parameters;
  PaToH_Initialize_Parameters(pargs, PATOH_CONPART, PATOH_SUGPARAM_DEFAULT);

  ones = new int[2 * hyperGraph->nhedges];
  tpart = new int[hyperGraph->nvtxs];

  fill(ones, ones + hyperGraph->nvtxs, 1);

  PaToH_Alloc(pargs, hyperGraph->nvtxs, hyperGraph->nhedges, 1, ones, ones,
      hyperGraph->eptr, hyperGraph->eind);

  pargs->_k = 2;
  PaToH_Part(pargs, hyperGraph->nvtxs, hyperGraph->nhedges, 1, 0, ones, ones,
      hyperGraph->eptr, hyperGraph->eind, tweight, tpart, pweight, &cut);

  delete[] ones;
  delete pargs;
  PaToH_Free();

  vector<int> * part = new vector<int>(tpart, tpart + hyperGraph->nvtxs);
  delete tpart;
  delete hyperGraph;

  return part;
}


static PartitionNode *
computeHierarchicalPartitioning(const GraphImpl * graph,
    vector<int>& scheduling, int start = -1, int end = -1) {
  PartitionNode * result = NULL;

  if (start == -1) {
    start = 0;
    end = scheduling.size() - 1;
  }

  if (start == end) {
    result = new PartitionNode(scheduling[start]);
    result->maxLive = 1;
  } else {
    int IOcost = 0;
    int intermediatePos = (start + end) / 2 + (start + end) % 2;
    set<int> targetSet(scheduling.begin() + intermediatePos,
        scheduling.begin() + end);

    PartitionNode * first = computeHierarchicalPartitioning(graph, scheduling,
        start, intermediatePos - 1);
    PartitionNode * second = computeHierarchicalPartitioning(graph, scheduling,
        intermediatePos, end);
    
    result = new PartitionNode(first, second);

    // Compute the IO size between the two children
    for (int i = start; i < intermediatePos; i++) {
      Vertex * curr = graph->getVertex(scheduling[i]);

      for (int j = 0, e = curr->getSuccessorCount(); j < e; j++) {
        if (targetSet.count(curr->getSuccessor(j)->getId())) {
          IOcost++;
          break;
        }
      }
    }

    // Compute the maxLive of this partition level as :
    // maxLive = max(maxLive(first), maxLive(second), IO(first, second))
    result->maxLive = first->maxLive;
    if (second->maxLive > result->maxLive)
      result->maxLive = second->maxLive;

    if (IOcost > result->maxLive)
      result->maxLive = IOcost;
  }

  return result;
}


static int
assignPartitions(PartitionNode * node, int memorySize,
    map<int, int>& partitions, int partitionId = 0) {
  if (node->id != -1) {
    partitions[node->id] = partitionId;
  } else {
    partitionId =
      assignPartitions(node->first, memorySize, partitions, partitionId);

    if (node->maxLive > memorySize)
      partitionId++;

    partitionId =
      assignPartitions(node->second, memorySize, partitions, partitionId);

    if (node->maxLive > memorySize)
      partitionId++;
  }

  return partitionId;
}


static int
findUnion(int id, vector<int>& unions) {
  while (id != unions[id])
    id = unions[id];

  return id;
}


static void
mergeUnion(int i, int j, vector<int>& unions, vector<int>& sizes) {
  i = findUnion(i, unions);
  j = findUnion(j, unions);

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

int
GraphImpl::getPartitionLevels() const {
  return partitions.size();
}


const vector<Vertex::PartitionArray *>&
GraphImpl::getPartitions() const {
  return partitions;
}


int
GraphImpl::getSchedulingCost(vector<int>& scheduling, int memorySize) const {
  int cost = 0;
  map<int, int> partitions;
  PartitionNode * partitionRoot = NULL;

  // First check the validity of the scheduling
  if (!checkSchedulingValidity(scheduling))
    return 0;

  // Construct the default hierarchical partitioning of the scheduling.
  partitionRoot = computeHierarchicalPartitioning(this, scheduling);

  // Obtain a partition assignement for the specified memory-size
  assignPartitions(partitionRoot, memorySize, partitions);

  // Compute the IO cost for each partition
  for (auto it = partitions.begin(), end = partitions.end(); it != end; ++it) {
    set<int> targetSet;
    Vertex * curr = getVertex(it->first);

    for (int i = 0, e = curr->getSuccessorCount(); i < e; i++)
      targetSet.insert(curr->getSuccessor(i)->getId());

    cost += targetSet.size();
  }

  return cost;
}


// Internal partitioning query functions

void
#if defined(ENABLE_TLS)
GraphImpl::processPartitionQuery(PartitionQueryImpl * query) {
#else // ENABLE_TLS
GraphImpl::processPartitionQuery(PartitionQueryImpl * query, int threadId) {
#endif // ENABLE_TLS
  
}


HGraph *
GraphImpl::getHGraph(Vertex::IdSet& idSet) const {
  int vCount = idSet.size();
  int eCount = 0;
  int ind;
  int idx;
  HGraph * hyperGraph = new HGraph();

  for (auto it = idSet.begin(), end = idSet.end(); it != end; ++it) {
    Vertex * curr = getVertex(*it);

    for (int i = 0, e = curr->getSuccessorCount(); i < e; i++)
      eCount += idSet.count(curr->getSuccessor(i)->getId());
  }

  hyperGraph->nvtxs = vCount;
  hyperGraph->nhedges = eCount;

  hyperGraph->eptr = new int[vCount + 2];
  hyperGraph->eind = new int[vCount + eCount + 2];

  hyperGraph->eptr[0] = 0;

  ind = 0;
  idx = 0;
  for (auto it = idSet.begin(), end = idSet.end(); it != end; ++it) {
    Vertex * curr = getVertex(*it);

    hyperGraph->eind[ind++] = idx;

    for (int i = 0, e = curr->getSuccessorCount(); i < e; i++)
      ind += idSet.count(curr->getSuccessor(i)->getId());

    hyperGraph->eptr[idx + 1] = ind;
    hyperGraph->forwardMap[*it] = idx;
    hyperGraph->backwardMap.push_back(*it);

    idx++;
  }

  ind = 0;
  idx = 0;
  for (int i = 0; i < vCount; i++) {
    int id = hyperGraph->eind[ind++];
    Vertex * curr = getVertex(hyperGraph->backwardMap[id]);

    for (int i = 0, e = curr->getSuccessorCount(); i < e; i++) {
      int succId = curr->getSuccessor(i)->getId();
      if (idSet.count(succId))
        hyperGraph->eind[ind++] = hyperGraph->forwardMap[succId];
    }
  }

  return hyperGraph;
}


void
GraphImpl::convBisect(Vertex::IdSet& idSet, int level) {
  int nbCuts[2] = { 0, 0};
  Vertex * curr = NULL;
  queue<int> ready;
  vector<int> nbPred;
  vector<int> nbSucc;
  vector<int> convPart;
  HGraph * hyperGraph = NULL;

  // Bail out if the source set is of size 1
  if (idSet.size() <= 1) {
    return;
  } else {
    while ((int) partitions.size() <= level)
      partitions.push_back(new Vertex::PartitionArray());
  }

  nbPred.resize(idSet.size(), 0);
  nbSucc.resize(idSet.size(), 0);
  convPart.resize(idSet.size(), 2);
  hyperGraph = getHGraph(idSet);

  // We determine the best order for the two parts
  vector<int> * part = patohBisect(hyperGraph);

  for (auto it = hyperGraph->forwardMap.begin(),
      end = hyperGraph->forwardMap.end(); it != end; ++it) {
    int pId = (*part)[it->second];
    Vertex * curr = getVertex(it->first);

    for(int i = 0, e = curr->getSuccessorCount(); i < e; i++) {
      int tId = (*part)[hyperGraph->forwardMap[curr->getSuccessor(i)->getId()]];
      if (pId != tId)
        nbCuts[pId]++;
    }
  }

  // If there are more cuts from 1 to 0, we switch parts
  if (nbCuts[1] > nbCuts[0]) {
    for (int i = 0, e = idSet.size(); i < e; i++)
      (*part)[i] = 1 - (*part)[i];
  }
  
  for (auto it = hyperGraph->forwardMap.begin(),
      end = hyperGraph->forwardMap.end(); it != end; ++it) {
    curr = getVertex(it->first);

    if (((*part)[it->second] == 0) && (curr->getPredecessorCount() == 0)) {
      ready.push(it->second);
      convPart[it->second] = 0;
    }

    if (((*part)[it->second] == 1) && (curr->getSuccessorCount() == 0)) {
      ready.push(it->second);
      convPart[it->second] = 1;
    }

    nbPred[it->second] = curr->getPredecessorCount();
    nbSucc[it->second] = curr->getSuccessorCount();
  }

  int nbOk[2] = {0, 0};
  while (!ready.empty()) {
    int id = ready.front();
    curr = getVertex(hyperGraph->backwardMap[id]);
    ready.pop();

    nbOk[(*part)[id]]++;
    if ((*part)[id] == 0) {
      for (int i = 0, e = curr->getSuccessorCount(); i < e; i++) {
        int target = hyperGraph->forwardMap[curr->getSuccessor(i)->getId()];

        if ((*part)[id] == (*part)[target]) {
          nbPred[target]--;
          if ((convPart[target] == 2) && (nbPred[target] == 0)) {
            ready.push(target);
            convPart[target] = (*part)[target];
          }
        }
      }
    } else {
      for (int i = 0, e = curr->getPredecessorCount(); i < e; i++) {
        int target = hyperGraph->forwardMap[curr->getPredecessor(i)->getId()];

        if ((*part)[id] == (*part)[target]) {
          nbSucc[target]--;
          if ((convPart[target] == 2) && (nbSucc[target] == 0)) {
            ready.push(target);
            convPart[target] = (*part)[target];
          }
        }
      }
    }
  }

  if ((nbOk[0] + nbOk[1]) < (int) idSet.size()) {
    set<int> subIdSet;
    set<int> targetSubIdSet;
    vector<int> idMap;

    for (int i = 0, e = idSet.size(); i < e; i++) {
      if (convPart[i] == 2)
        subIdSet.insert(hyperGraph->backwardMap[i]);
    }

    for (int i = 0, e = subIdSet.size(); i < e; i++)
      targetSubIdSet.insert(i);

    GraphImpl * subGraph = getSubGraph(subIdSet, idMap);
    subGraph->convBisect(targetSubIdSet, 1);

    bool invert = false;
    const vector<Vertex::PartitionArray *>& subPartitions =
      subGraph->getPartitions();
    
    const Vertex::IdSet * subSetA = subPartitions.front()->front();
    const Vertex::IdSet * subSetB = subPartitions.front()->back();

    for (auto it = subSetA->begin(), end = subSetA->end(); it != end; ++it) {
      Vertex * curr = subGraph->getVertex(*it);

      for (int i = 0, e = curr->getPredecessorCount(); i < e; i++) {
        if (!subSetA->count(curr->getPredecessor(i)->getId())) {
          invert = true;
          break;
        }
      }

      if (invert) {
        subSetA = subPartitions.front()->back();
        subSetB = subPartitions.front()->front();

        break;
      }
    }

    for (auto it = subSetA->begin(), end = subSetA->end(); it != end; ++it)
      convPart[idMap[*it]] = 0;

    for (auto it = subSetB->begin(), end = subSetB->end(); it != end; ++it)
      convPart[idMap[*it]] = 1;

    delete subGraph;
  }

  // Add the two new partitions to the partioning hierarchy
  Vertex::IdSet * setA = new Vertex::IdSet();
  Vertex::IdSet * setB = new Vertex::IdSet();
  for (int i = 0, e = idSet.size(); i < e; i++) {
    if (convPart[i] == 0) {
      setA->insert(hyperGraph->backwardMap[i]);
    } else if (convPart[i] == 1) { 
      setB->insert(hyperGraph->backwardMap[i]);
    } else {
      cerr << "ERROR: Found a node not in partition after convex-bisection.\n";
      exit(EXIT_FAILURE);
    }
  }

  partitions[level]->push_back(setA);
  partitions[level]->push_back(setB);
}


void
GraphImpl::maxDistBisect(Vertex::IdSet& idSet, vector<int>& sourceMax,
    vector<int>& sinkMax, int level) {
  // Bail out if the source set is of size 1
  if (idSet.size() == 1) {
    return;
  } else {
    while ((int) partitions.size() <= level)
      partitions.push_back(new Vertex::PartitionArray());
  }

  /* First check for connectiviy in the graph through union-find */
  vector<int> unions(*idSet.end() + 1, -1);
  vector<int> sizes(*idSet.end() + 1, 1);

  // Initialize
  int id = 0;
  for (auto it = idSet.begin(), end = idSet.end(); it != end; ++it)
    unions[*it] = id++;

  // Perform merging
  for (auto it = idSet.begin(), end = idSet.end(); it != end; ++it) {
    Vertex * curr = getVertex(*it);

    for (int i = 0, e = curr->getSuccessorCount(); i < e; i++) {
      int succ = curr->getSuccessor(i)->getId();

      if (idSet.count(succ))
        mergeUnion(*it, succ, unions, sizes);
    }
  }

  // If the graph is not connected perform trivial multisection
  int target = findUnion(*idSet.begin(), unions);
  if (sizes[target] < (int) idSet.size()) {
    map<int, Vertex::IdSet *> subSets;

    for (auto it = idSet.begin(), end = idSet.end(); it != end; ++it) {
      target = findUnion(*it, unions);

      if (!subSets.count(target)) {
        Vertex::IdSet * newSet = new Vertex::IdSet();

        subSets[target] = newSet;
        partitions[level]->push_back(newSet);
      }

      subSets[target]->insert(*it);
    }

    return;
  }

  /* Perform max-distance based bisection for a connected graph */
  Vertex::IdSet * setA = new Vertex::IdSet();
  Vertex::IdSet * setB = new Vertex::IdSet();

  for (auto it = idSet.begin(), end = idSet.end(); it != end; ++it) {
    if (sourceMax[*it] < sinkMax[*it]) {
      setA->insert(*it);
      sinkMax[*it] = sinkMax[*it] / 2;
    } else {
      setB->insert(*it);
      sourceMax[*it] = sourceMax[*it] / 2 + sourceMax[*it] % 2;
    }
  }

  partitions[level]->push_back(setA);
  partitions[level]->push_back(setB);
}


void
GraphImpl::partitionConvex(void) {
  // Erase existing partitions
  partitions.clear();

  Vertex::IdSet * topSet = new Vertex::IdSet();

  for (auto it = vertices.begin(), end = vertices.end(); it != end; ++it) {
    if (*it)
      topSet->insert((*it)->getId());
  }

  partitions.push_back(new Vertex::PartitionArray);
  partitions[0]->push_back(topSet);

  int level = 0;
  while (level < (int) partitions.size()) {
    for (auto it = partitions[level]->begin(), end = partitions[level]->end();
        it != end; ++it)
      convBisect(**it, level + 1);

    level++;
  }
}


void
GraphImpl::partitionMaxDistance(void) {
  vector<int> sourceMax;
  vector<int> sinkMax;

  // Erase existing partitions
  partitions.clear();

  Vertex::IdSet * topSet = new Vertex::IdSet();

  for (auto it = vertices.begin(), end = vertices.end(); it != end; ++it) {
    if (*it)
      topSet->insert((*it)->getId());
  }

  partitions.push_back(new Vertex::PartitionArray);
  partitions[0]->push_back(topSet);

  sourceMax.resize(vertices.size() + 1, -1);
  sinkMax.resize(vertices.size() + 1, -1);

  /* Perform the initial computations of the max-distances */
  stack<Vertex *> todo;

  // Construct an artificial root node for all sources and sinks
  VertexImpl root(vertices.size());
  for (auto it = sources.begin(), end = sources.end(); it != end; ++it)
    root.addSuccessor(*it);

  for (auto it = sinks.begin(), end = sinks.end(); it != end; ++it)
    root.addPredecessor(*it);

  // Perform a forward DFS for sourceMax
  todo.push(&root);
  while (!todo.empty()) {
    Vertex * curr = todo.top();
    todo.pop();

    for (int i = 0, e = curr->getSuccessorCount(); i < e; i++) {
      Vertex * succ = curr->getSuccessor(i);
      
      if (sourceMax[succ->getId()] < sourceMax[curr->getId()] + 1) {
        sourceMax[succ->getId()] = sourceMax[curr->getId()] + 1;
        todo.push(succ);
      }
    }
  }

  // Perform a backward DFS for sinkMax
  todo.push(&root);
  while (!todo.empty()) {
    Vertex * curr = todo.top();
    todo.pop();

    for (int i = 0, e = curr->getPredecessorCount(); i < e; i++) {
      Vertex * pred = curr->getPredecessor(i);
      
      if (sinkMax[pred->getId()] < sinkMax[pred->getId()] + 1) {
        sinkMax[pred->getId()] = sinkMax[pred->getId()] + 1;
        todo.push(pred);
      }
    }
  }

  /* Perform the real hierarchical partitioning */
  int level = 0;
  while (level < (int) partitions.size()) {
    for (auto it = partitions[level]->begin(), end = partitions[level]->end();
        it != end; ++it)
      maxDistBisect(**it, sourceMax, sinkMax, level + 1);

    level++;
  }
}

bool
GraphImpl::checkSchedulingValidity(vector<int>& scheduling) const {
  set<int> idSet;

  if (scheduling.size() != getVertexCount()) {
    cerr << "WARNING: The supplied scheduling in 'getSchedulingCost' is not "
      << "valid. It does not contain the same number of vertices as the "
      << "graph.\n";
    return false;
  }

  for (int i = 0, e = scheduling.size(); i < e; i++) {
    Vertex * curr = getVertex(scheduling[i]);

    if (!curr) {
      cerr << "WARNING: The supplied scheduling in 'getSchedulingCost' contains"
        << " an ID that does not correspond to an existing vertex.\n";
      return false;
    }

    if (idSet.count(scheduling[i])) {
      cerr << "WARNING: The supplied scheduling in 'getSchedulingCost' contains"
        << " two times the same ID (" << scheduling[i] << ").\n";
      return false;
    }

    for (auto j = 0, e2 = curr->getSuccessorCount(); j < e2; j++) {
      if (idSet.count(curr->getSuccessor(j)->getId())) {
        cerr << "WARNING: The supplied scheduling in 'getSchedulingCost' is "
          << "illegal. ID (" << curr->getSuccessor(j)->getId() << ") should be "
          << "scheduling after ID (" << scheduling[i] << ").\n";
        return false;
      }
    }

    idSet.insert(scheduling[i]);
  }

  return true;
}
