#include <map>
#include <iostream>

#include "graph-utilities/implementation/graphImpl.hpp"

using namespace std;


// Local variables

static map<int, int> remappedVertices;


// Modificators

Vertex *
GraphImpl::addVertex(int weight) {
  Vertex * newVertex = NULL;

  startGlobalOperation();

  newVertex = addVertexUnsafe(weight);

  stopGlobalOperation();

  return newVertex;
}


void
GraphImpl::removeVertex(Vertex * vertex) {
  int remapTarget = -1;
  auto remapping = remappedVertices.find(vertex->getId());

  if (remapping != remappedVertices.end())
    remapTarget = remapping->second;

  for (int i = 0, e = vertex->getPredecessorCount(); i < e; i++) {
    vertex->getPredecessor(i)->removeSuccessor(vertex);
    edgeCount--;
  }

  for (int i = 0, e = vertex->getSuccessorCount(); i < e; i++) {
    vertex->getSuccessor(i)->removePredecessor(vertex);
    edgeCount--;
  }

  if (remapTarget != -1) {
    for (auto it = remappedVertices.begin(), end = remappedVertices.end();
        it != end; ++it) {
      if (it->second == vertex->getId())
        it->second = remapTarget;
    }
  }

  vertexCount--;

  vertices[vertex->getId()] = NULL;
  delete vertex;
}


Vertex *
GraphImpl::mergeVertices(const Vertex::Set &vertexSet) {
  VertexImpl * target = addVertexUnsafe(threadCount);
  VertexImpl::Set internalSet;
  UserDataInterface * mergedData = NULL;
  startGlobalOperation();

  int targetWeight = 0;
  map<VertexImpl *, int> predWeightMap;
  map<VertexImpl *, int> succWeightMap;

  for (auto setIt = vertexSet.begin(), setEnd = vertexSet.end();
      setIt != setEnd; ++setIt)
    internalSet.insert(dynamic_cast<VertexImpl *>(*setIt));

  for (auto setIt = internalSet.begin(), setEnd = internalSet.end();
      setIt != setEnd; ++setIt) {
    for (int i = 0, e = (*setIt)->getPredecessorCount(); i < e; i++) {
      VertexImpl * pred = (*setIt)->getPredecessorI(i);
      if (vertexSet.find(pred) == vertexSet.end()) {
        auto mapIt = predWeightMap.find(pred);
        if (mapIt != predWeightMap.end())
          mapIt->second += (*setIt)->getPredecessorWeight(i);
        else
          predWeightMap[pred] = (*setIt)->getPredecessorWeight(i);
      }
    }

    for (int i = 0, e = (*setIt)->getSuccessorCount(); i < e; i++) {
      VertexImpl * succ = (*setIt)->getSuccessorI(i);
      if (vertexSet.find(succ) == vertexSet.end()) {
        auto mapIt = succWeightMap.find(succ);
        if (mapIt != succWeightMap.end())
          mapIt->second += (*setIt)->getSuccessorWeight(i);
        else
          succWeightMap[succ] = (*setIt)->getSuccessorWeight(i);
      }
    }

    targetWeight += (*setIt)->getWeight();

    auto remapIt = remappedVertices.find((*setIt)->getId());
    if (remapIt != remappedVertices.end())
      remapIt->second = target->getId();
    else
      remappedVertices[(*setIt)->getId()] = target->getId();

    if ((*setIt)->getUserData() != NULL) {
      if (mergedData == NULL)
        mergedData = (*setIt)->getUserData()->clone();
      else
        mergedData->merge((*setIt)->getUserData());
    }

    removeVertex(*setIt);
  }

  target->setWeight(targetWeight);

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
GraphImpl::addEdge(Vertex * source, Vertex * target, int weight) {
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
GraphImpl::removeEdge(Vertex * source, Vertex * target) {
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
GraphImpl::setIndexMethod(IndexMethod newMethod) {
  indexMethod = newMethod;
}


// Access

Vertex *
GraphImpl::getVertex(int id) const {
  auto remapping = remappedVertices.find(id);

  if (remapping != remappedVertices.end())
    return vertices.at(remapping->second);
  else
    return vertices.at(id);
}


unsigned int
GraphImpl::getEdgeCount() const {
  return edgeCount;
}


unsigned int
GraphImpl::getVertexCount() const {
  return vertexCount;
}


IndexMethod
GraphImpl::getIndexMethod() const {
  return indexMethod;
}


VertexImpl::Array&
GraphImpl::getSourceVertices() {
  return sources;
}


VertexImpl::Array&
GraphImpl::getSinkVertices() {
  return sinks;
}


// Checks

void
GraphImpl::verifyVertexIds() const {
  int i = 0;
  int count = getVertexCount();

  for (auto it = vertices.begin(), end = vertices.end(); it != end; ++it) {
    if (*it == NULL) {
      cerr << "Vertex slot n°" << i << " was NULL\n";
    } else if ((*it)->getId() > count) {
      cerr << "Vertex slot n°" << i << " has an abnormal ID ("
        << (*it)->getId() << ")\n";
    }

    i++;
  }
}


// Maintenance

VertexImpl *
GraphImpl::addVertexUnsafe(int weight) {
  int id = vertices.size();
  VertexImpl * newVertex = new VertexImpl(id, weight);

  vertices.push_back(newVertex);
  indexed = false;
  vertexCount++;

  return newVertex;
}


bool
GraphImpl::addEdgeUnsafe(VertexImpl * source, VertexImpl * target, int weight) {
  source->addSuccessorUnsafe(target, weight);
  target->addPredecessorUnsafe(source, weight);
  edgeCount++;
  return true;
}


Graph *
createEmptyGraph(void) {
  return new GraphImpl();
}
