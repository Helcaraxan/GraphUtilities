#include <map>
#include <iostream>

#include "graph-utilities/defs.hpp"
#include "graph-utilities/graph.hpp"
#include "graph-utilities/vertex.hpp"

using namespace std;


// Local variables

static map<int, int> remappedVertices;


// Modificators

Vertex *
Graph::addVertex(void) {
  Vertex * newVertex = NULL;

  startGlobalOperation();

  newVertex = addVertexUnsafe(threadCount);

  stopGlobalOperation();

  return newVertex;
}


void
Graph::removeVertex(Vertex * v) {
  int remapTarget = -1;
  auto remapping = remappedVertices.find(v->id);

  if (remapping != remappedVertices.end())
    remapTarget = remapping->second;

  for (int i = 0, e = v->getPredecessorCount(); i < e; i++) {
    v->getPredecessor(i)->removeSuccessor(v);
    edgeCount--;
  }

  for (int i = 0, e = v->getSuccessorCount(); i < e; i++) {
    v->getSuccessor(i)->removePredecessor(v);
    edgeCount--;
  }

  if (remapTarget != -1) {
    for (auto it = remappedVertices.begin(), end = remappedVertices.end();
        it != end; ++it) {
      if (it->second == v->id)
        it->second = remapTarget;
    }
  }

  vertexCount--;

  vertices[v->id] = NULL;
  delete v;
}


Vertex *
Graph::mergeVertices(set<Vertex *> &s) {
  Vertex * target = addVertexUnsafe(threadCount);
  UserDataInterface * mergedData = NULL;
  startGlobalOperation();

  int targetWeight = 0;
  map<Vertex *, int> predWeightMap;
  map<Vertex *, int> succWeightMap;

  for (auto mergeIt = s.begin(), mergeEnd = s.end();
      mergeIt != mergeEnd; ++mergeIt) {
    auto weightIt = (*mergeIt)->predecessorWeights.begin();
    for (auto predIt = (*mergeIt)->pred_begin(),
        predEnd = (*mergeIt)->pred_end(); predIt != predEnd; ++predIt) {
      if (s.find(*predIt) == s.end()) {
        auto mapIt = predWeightMap.find(*predIt);
        if (mapIt != predWeightMap.end())
          mapIt->second += *weightIt;
        else
          predWeightMap[*predIt] = *weightIt;
      }

      weightIt++;
    }

    weightIt = (*mergeIt)->successorWeights.begin();
    for (auto succIt = (*mergeIt)->succ_begin(),
        succEnd = (*mergeIt)->succ_end(); succIt != succEnd; ++succIt) {
      if (s.find(*succIt) == s.end()) {
        auto mapIt = succWeightMap.find(*succIt);
        if (mapIt != succWeightMap.end())
          mapIt->second += *weightIt;
        else
          succWeightMap[*succIt] = *weightIt;
      }

      weightIt++;
    }

    targetWeight += (*mergeIt)->weight;

    auto remapIt = remappedVertices.find((*mergeIt)->id);
    if (remapIt != remappedVertices.end())
      remapIt->second = target->id;
    else
      remappedVertices[(*mergeIt)->id] = target->id;

    if ((*mergeIt)->userData != NULL) {
      if (mergedData == NULL)
        mergedData = (*mergeIt)->userData->clone();
      else
        mergedData->merge((*mergeIt)->userData);
    }

    removeVertex(*mergeIt);
  }

  target->weight = targetWeight;

  for (auto mapIt = predWeightMap.begin(), mapEnd = predWeightMap.end();
      mapIt != mapEnd; ++mapIt)
    addEdgeUnsafe(mapIt->first, target, mapIt->second);

  for (auto mapIt = succWeightMap.begin(), mapEnd = succWeightMap.end();
      mapIt != mapEnd; ++mapIt)
    addEdgeUnsafe(target, mapIt->first, mapIt->second);

  if (mergedData)
    target->setUserData(mergedData);

  indexed = false;
  stopGlobalOperation();
  
  return target;
}


bool
Graph::addEdge(Vertex * source, Vertex * target, int weight) {
  startGlobalOperation();

  if (!source->addSuccessor(target, weight))
    return false;

  target->addPredecessor(source, weight);
  edgeCount++;
  indexed = false;
  condensed = false;

  stopGlobalOperation();
  return true;
}


bool
Graph::removeEdge(Vertex * source, Vertex * target) {
  startGlobalOperation();

  if (!source->removeSuccessor(target)) {
    stopGlobalOperation();
    return false;
  }

  source->removePredecessor(source);
  edgeCount--;

  stopGlobalOperation();

  return true;
}


void
Graph::setIndexMethod(IndexMethod newMethod) {
  indexMethod = newMethod;
}


// Access

unsigned int
Graph::getEdgeCount() const {
  return edgeCount;
}


unsigned int
Graph::getVertexCount() const {
  return vertexCount;
}


Vertex *
Graph::getVertexFromId(int id) const {
  auto remapping = remappedVertices.find(id);

  if (remapping != remappedVertices.end())
    return vertices.at(remapping->second);
  else
    return vertices.at(id);
}


IndexMethod
Graph::getIndexMethod() const {
  return indexMethod;
}


// Checks

void
Graph::verifyVertexIds() const {
  int i = 0;
  int count = getVertexCount();

  for (auto it = vertices.begin(), end = vertices.end(); it != end; ++it) {
    if (*it == NULL) {
      cerr << "Vertex slot n°" << i << " was NULL" << endl;
    } else if ((*it)->id > count) {
      cerr << "Vertex slot n°" << i << " has an abnormal ID (" << (*it)->id;
      cerr << ")" << endl;
    }

    i++;
  }
}


// Maintenance

Vertex *
Graph::addVertexUnsafe(int threadCount) {
  int id = vertices.size();
#ifdef ENABLE_TLS
  Vertex * newVertex = new Vertex(id);
#else // ENABLE_TLS
  Vertex * newVertex = new Vertex(threadCount, id);
#endif // ENABLE_TLS

  vertices.push_back(newVertex);
  indexed = false;
  vertexCount++;

  return newVertex;
}


bool
Graph::addEdgeUnsafe(Vertex * source, Vertex * target, int weight) {
  source->addSuccessorUnsafe(target, weight);
  target->addPredecessorUnsafe(source, weight);
  edgeCount++;
  return true;
}


