#include "graph-utilities/implementation/queriesImpl.hpp"

using namespace std;


/* CoarsenQueryImpl Interface implementation */
bool
CoarsenQueryImpl::getError() const {
  return error;
}


int
CoarsenQueryImpl::getFactor() const {
  return factor;
}


int
CoarsenQueryImpl::getSecondaryFactor() const {
  return secondaryFactor;
}


Graph *
CoarsenQueryImpl::getAnswer() {
  return coarsenedGraph;
}


CoarsenMethod
CoarsenQueryImpl::getMethod() const {
  return method;
}


map<int, int>&
CoarsenQueryImpl::getMap() {
  return vertexMap;
}


// Modifications
void
CoarsenQueryImpl::setError(bool value) {
  error = value;
}


void
CoarsenQueryImpl::setFactor(int newFactor) {
  factor = newFactor;
}


void
CoarsenQueryImpl::setSecondaryFactor(int newFactor) {
  secondaryFactor = newFactor;
}


void
CoarsenQueryImpl::setMethod(CoarsenMethod newMethod) {
  method = newMethod;
}


CoarsenQuery *
createCQuery() {
  return new CoarsenQueryImpl();
}


/* CoarsenQueryImpl implementations */

// Modifications
void
CoarsenQueryImpl::setAnswer(GraphImpl * newGraph) {
  coarsenedGraph = newGraph;
}


/* PartitionQuery Interface implementations */
// Access

bool
PartitionQueryImpl::getError() const {
  return error;
}


PartitionMethod
PartitionQueryImpl::getMethod() const {
  return method;
}


int
PartitionQueryImpl::getThreadCount() const {
  return threadCount;
}


int
PartitionQueryImpl::getMemorySize() const {
  return memorySize;
}


const Partition *
PartitionQueryImpl::getPartition() const {
  return partition;
}


// Modifications

void
PartitionQueryImpl::setMethod(PartitionMethod newMethod) {
  method = newMethod;
}


void
PartitionQueryImpl::setThreadCount(int newCount) {
  threadCount = newCount;
}


void
PartitionQueryImpl::setMemorySize(int newSize) {
  memorySize = newSize;
}


/* PartitionQueryImpl implementations */
// Modifications

void
PartitionQueryImpl::setError(bool value) {
  error = value;
}


void
PartitionQueryImpl::setPartition(PartitionImpl * part) {
  partition = part;
}


PartitionQuery *
createPQuery() {
  return new PartitionQueryImpl();
}


/* ReachabilityQuery Interface implementations */
// Access

Vertex *
ReachabilityQueryImpl::getSource() const {
  return source;
}


Vertex *
ReachabilityQueryImpl::getTarget() const {
  return target;
}


bool
ReachabilityQueryImpl::getAnswer() const {
  return answer;
}


bool
ReachabilityQueryImpl::getError() const {
  return error;
}


SearchMethod
ReachabilityQueryImpl::getMethod() const {
  return method;
}


// Modifications

void
ReachabilityQueryImpl::setSource(Vertex * sourceVertex) {
  VertexImpl * castSource = dynamic_cast<VertexImpl *>(sourceVertex);

  if (castSource)
    source = castSource;
}


void
ReachabilityQueryImpl::setTarget(Vertex * targetVertex) {
  VertexImpl * castTarget = dynamic_cast<VertexImpl *>(targetVertex);

  if (castTarget)
    target = castTarget;
}


void
ReachabilityQueryImpl::setMethod(SearchMethod value) {
  method = value;
}


ReachabilityQuery *
createRQuery() {
  return new ReachabilityQueryImpl();
}


/* ReachabilityQueryImpl implementations */
// Access

VertexImpl *
ReachabilityQueryImpl::getSourceI() const {
  return source;
}


VertexImpl *
ReachabilityQueryImpl::getTargetI() const {
  return target;
}


bool
ReachabilityQueryImpl::getCancel() const {
  return cancel;
}


ReachabilityQueryImpl::InternalQueue *
ReachabilityQueryImpl::getInternal() const {
  return internal;
}


#if defined(ENABLE_STATISTICS)
VertexImpl::Array&
ReachabilityQueryImpl::getPath() {
  return path;
}


int
ReachabilityQueryImpl::getSearchedNodes() {
  return searchedNodes;
}
#endif // ENABLE_STATISTICS


// Modifications

void
ReachabilityQueryImpl::setError(bool value) {
  error = value;
}


void
ReachabilityQueryImpl::setAnswer(bool value) {
  answer = value;
}


void
ReachabilityQueryImpl::setCancel(bool value) {
  cancel = value;
}


void
ReachabilityQueryImpl::setInternal(
    ReachabilityQueryImpl::InternalQueue * value) {
  internal = value;
}


#if defined(ENABLE_STATISTICS)
void
ReachabilityQueryImpl::nodeSearched() {
  searchedNodes++;
}
#endif // ENABLE_STATISTICS
