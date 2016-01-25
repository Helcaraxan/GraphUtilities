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
static PartitionImpl * part = NULL;
static stack<Vertex::IdSet *> subSetA;
static stack<Vertex::IdSet *> subSetB;


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
  delete[] tpart;

  return part;
}


void
computeHierarchicalPartition(const GraphImpl * graph, vector<int>& scheduling,
    int start = -1, int end = -1) {
  if (part == NULL) {
    start = 0;
    end = scheduling.size() - 1;

    part = new PartitionImpl(scheduling.size());
    part->root->id = 0;
    for (int i = 0, e = scheduling.size(); i < e; i++)
      part->representants[i] = i + 1;

    part->representants[end] = part->root;
  }

  IaP<PartitionNodeImpl> ptr = part->representants[start];
  while (ptr.isInteger())
    ptr = part->representants[ptr];

  PartitionNodeImpl * parent = ptr;
  if (start == end) {
    PartitionNodeImpl * node = new PartitionNodeImpl(parent);

    node->id = start;
    node->maxLive = 0;
    part->representants[start] = node;
    parent->id = -1;
  } else {
    int IOcost = 0;
    int intermediatePos = (start + end) / 2 + (start + end) % 2;
    set<int> targetSet(scheduling.begin() + intermediatePos,
        scheduling.begin() + end);
    PartitionNodeImpl * first = new PartitionNodeImpl(parent);
    PartitionNodeImpl * second = new PartitionNodeImpl(parent);

    parent->id = -1;

    part->representants[start] = first;
    part->representants[intermediatePos] = second;

    computeHierarchicalPartition(graph, scheduling, start, intermediatePos - 1);
    computeHierarchicalPartition(graph, scheduling, intermediatePos, end);

    // Compute the IO size between the two children
    for (int i = start; i < intermediatePos; i++) {
      Vertex * curr = graph->getVertex(scheduling[i]);

      for (int j = 0, e = curr->getSuccessorCount(); j < e; j++) {
        if (targetSet.count(curr->getSuccessor(j)->getId()))
          IOcost += curr->getSuccessorWeight(j);
      }
    }

    // Compute the maxLive of this partition level as :
    // maxLive = max(maxLive(first), maxLive(second), IO(first, second))
    parent->maxLive = first->maxLive;
    if (second->maxLive > parent->maxLive)
      parent->maxLive = second->maxLive;

    if (IOcost > parent->maxLive)
      parent->maxLive = IOcost;
  }
}


static int
assignTiles(PartitionNodeImpl * node, int memorySize, vector<int>& scheduling,
    map<int, int>& tileMap, int tileId = 0) {
  if (node->id != -1) {
    tileMap[scheduling[node->id]] = tileId;
  } else {
    for (auto it = node->children.begin(), end = node->children.end();
        it != end; ++it) {
      tileId = assignTiles(*it, memorySize, scheduling, tileMap, tileId);

      if ((*it)->maxLive > memorySize)
        tileId++;
    }
  }

  return tileId;
}


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

int
GraphImpl::getSchedulingCost(vector<int>& scheduling, int memorySize) const {
  int cost = 0;
  map<int, int> tileMap;

  // First check the validity of the scheduling
  if (!checkSchedulingValidity(scheduling))
    return 0;

  // Construct the default hierarchical partition of the scheduling.
  computeHierarchicalPartition(this, scheduling);

  // Obtain a tile assignement for the specified memory-size
  assignTiles(part->root, memorySize, scheduling, tileMap);

  // Compute the IO cost for each tile
  for (auto it = tileMap.begin(), end = tileMap.end(); it != end; ++it) {
    set<int> targetSet;
    Vertex * curr = getVertex(it->first);

    for (int i = 0, e = curr->getSuccessorCount(); i < e; i++) {
      if (tileMap[curr->getSuccessor(i)->getId()] != it->second)
        cost += curr->getSuccessorWeight(i);
    }
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
  int vCount = idSet.size();
  int eCount = 0;
  int ind = 0;
  int idx = 0;
  HGraph * hyperGraph = new HGraph();

  hyperGraph->backwardMap.resize(idSet.size());
  for (auto it = idSet.begin(), end = idSet.end(); it != end; ++it) {
    Vertex * curr = getVertex(*it);

    hyperGraph->forwardMap[*it] = idx;
    hyperGraph->backwardMap[idx++] = *it;

    for (int i = 0, e = curr->getSuccessorCount(); i < e; i++)
      eCount += idSet.count(curr->getSuccessor(i)->getId());
  }

  hyperGraph->nvtxs = vCount;
  hyperGraph->nhedges = vCount;

  hyperGraph->eptr = new int[vCount + 2];
  hyperGraph->eind = new int[vCount + eCount + 2];

  hyperGraph->eptr[0] = 0;

  idx = 0;
  for (auto it = idSet.begin(), end = idSet.end(); it != end; ++it) {
    Vertex * curr = getVertex(*it);

    hyperGraph->eind[ind++] = hyperGraph->forwardMap[*it];
    for (int i = 0, e = curr->getSuccessorCount(); i < e; i++) {
      int succ = curr->getSuccessor(i)->getId();

      if (idSet.count(succ))
        hyperGraph->eind[ind++] = hyperGraph->forwardMap[succ];
    }

    hyperGraph->eptr[++idx] = ind;
  }

  return hyperGraph;
}


void
GraphImpl::convBisect(PartitionNodeImpl * parent) {
  int nbCuts[2] = { 0, 0};
  Vertex * curr = NULL;
  HGraph * hyperGraph = NULL;
  Vertex::IdSet idSet;
  queue<int> ready;
  vector<int> nbPred;
  vector<int> nbSucc;
  vector<int> convPart;

  // Create the target idSet and the associated hyper-graph
  IaP<PartitionNodeImpl> ptr = parent->id;
  while (ptr.isInteger()) {
    idSet.insert(ptr);
    ptr = part->representants[ptr];
  }

  hyperGraph = getHGraph(idSet);
  nbPred.resize(idSet.size(), 0);
  nbSucc.resize(idSet.size(), 0);
  convPart.resize(idSet.size(), 2);

  // We determine the best order for the two parts
  vector<int> * pPart = patohBisect(hyperGraph);

  for (auto it = hyperGraph->forwardMap.begin(),
      end = hyperGraph->forwardMap.end(); it != end; ++it) {
    int pId = (*pPart)[it->second];
    Vertex * curr = getVertex(it->first);

    for(int i = 0, e = curr->getSuccessorCount(); i < e; i++) {
      int tId = (*pPart)[hyperGraph->forwardMap[curr->getSuccessor(i)->getId()]];
      if (pId != tId)
        nbCuts[pId]++;
    }
  }

  // If there are more cuts from 1 to 0, we switch parts
  if (nbCuts[1] > nbCuts[0]) {
    for (int i = 0, e = idSet.size(); i < e; i++)
      (*pPart)[i] = 1 - (*pPart)[i];
  }
  
  for (auto it = hyperGraph->forwardMap.begin(),
      end = hyperGraph->forwardMap.end(); it != end; ++it) {
    curr = getVertex(it->first);

    if (((*pPart)[it->second] == 0) && (curr->getPredecessorCount() == 0)) {
      ready.push(it->second);
      convPart[it->second] = 0;
    }

    if (((*pPart)[it->second] == 1) && (curr->getSuccessorCount() == 0)) {
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

    nbOk[(*pPart)[id]]++;
    if ((*pPart)[id] == 0) {
      for (int i = 0, e = curr->getSuccessorCount(); i < e; i++) {
        int target = hyperGraph->forwardMap[curr->getSuccessor(i)->getId()];

        if ((*pPart)[id] == (*pPart)[target]) {
          nbPred[target]--;
          if ((convPart[target] == 2) && (nbPred[target] == 0)) {
            ready.push(target);
            convPart[target] = (*pPart)[target];
          }
        }
      }
    } else {
      for (int i = 0, e = curr->getPredecessorCount(); i < e; i++) {
        int target = hyperGraph->forwardMap[curr->getPredecessor(i)->getId()];

        if ((*pPart)[id] == (*pPart)[target]) {
          nbSucc[target]--;
          if ((convPart[target] == 2) && (nbSucc[target] == 0)) {
            ready.push(target);
            convPart[target] = (*pPart)[target];
          }
        }
      }
    }
  }

  if ((nbOk[0] + nbOk[1]) < (int) idSet.size()) {
    set<int> subIdSet;
    set<int> targetSubIdSet;
    vector<int> idMap;

    IaP<PartitionNodeImpl> ptr = parent;
    for (int i = 0, e = idSet.size(); i < e; i++) {
      if (convPart[i] == 2) {
        int id = hyperGraph->backwardMap[i];
        part->representants[id] = ptr;
        parent->id = id;
        ptr = id;
      }
    }

    subSetA.push(new Vertex::IdSet());
    subSetB.push(new Vertex::IdSet());

    convBisect(parent);

    bool invert = false;
    for (auto it = subSetA.top()->begin(), end = subSetA.top()->end();
        it != end; ++it) {
      Vertex * curr = getVertex(*it);

      for (int i = 0, e = curr->getPredecessorCount(); i < e; i++) {
        if (!subSetA.top()->count(curr->getPredecessor(i)->getId())) {
          invert = true;
          break;
        }
      }

      if (invert)
        break;
    }

    for (auto it = subSetA.top()->begin(), end = subSetA.top()->end();
        it != end; ++it)
      convPart[idMap[*it]] = invert ? 1 : 0;

    for (auto it = subSetB.top()->begin(), end = subSetB.top()->end();
        it != end; ++it)
      convPart[idMap[*it]] = invert ? 0 : 1;

    delete subSetA.top();
    delete subSetB.top();

    subSetA.pop();
    subSetB.pop();
  }

  // Add the two new partitions to the partioning hierarchy
  if (!subSetA.empty()) {
    for (int i = 0, e = idSet.size(); i < e; i++) {
      if (convPart[i] == 0)
        subSetA.top()->insert(hyperGraph->backwardMap[i]);
      else
        subSetB.top()->insert(hyperGraph->backwardMap[i]);
    }
  } else {
    int memberA = -1;
    int memberB = -1;
    PartitionNodeImpl * childA = new PartitionNodeImpl(parent);
    PartitionNodeImpl * childB = new PartitionNodeImpl(parent);

    for (int i = 0, e = idSet.size(); i < e; i++) {
      int idx = hyperGraph->backwardMap[i];

      if (convPart[i] == 0) {
        if (memberA == -1)
          part->representants[idx] = childA;
        else
          part->representants[idx] = memberA;

        memberA = idx;
      } else if (convPart[i] == 1) { 
        if (memberB == -1)
          part->representants[idx] = childB;
        else
          part->representants[idx] = memberB;

        memberB = idx;
      } else {
        cerr << "ERROR: Found node not in partition after convex-bisection.\n";
        exit(EXIT_FAILURE);
      }
    }

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

  query->setPartition(part);
  part = NULL;
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
  part = NULL;
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
