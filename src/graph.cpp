#include <map>
#include <queue>
#include <chrono>
#include <climits>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <algorithm>

#include <unistd.h>

#include "graph-utilities/graph.hpp"

using namespace std;


// Local global variables

#ifdef ENABLE_BENCHMARKS
static PAPI_dmem_info_t memoryInfo;
static long long baseMemoryUsage;
#endif // ENABLE_BENCHMARKS

static map<int, int> remappedVertices;

// Semaphore class (passive wait variant)

void
semaphore::post() {
  unique_lock<mutex> lck(mtx);
  atomic_fetch_add_explicit(&count, 1, memory_order_acq_rel);
  cv.notify_one();
}


void
semaphore::wait() {
  unique_lock<mutex> lck(mtx);
  while (count.load(memory_order_acquire) == 0)
    cv.wait(lck);

  atomic_fetch_add_explicit(&count, -1, memory_order_acq_rel);
}


// Spin-semaphore class (passive wait variant)

void
spin_semaphore::initialize(int memberCount) {
  unique_lock<mutex> lck(mtx);

  if (cvs) {
    cerr << "ERROR: Trying to initialize an already initialized spin-semaphore." << endl;
    exit(EXIT_FAILURE);
  }

  cvs = new condition_variable[memberCount];
  members = memberCount;
  maxIdx = 0;
  currIdx = 0;
  overflow = 0;
}


void
spin_semaphore::reset() {
  unique_lock<mutex> lck(mtx, defer_lock);

  if (cvs) {
    while (currIdx != maxIdx) {
      lck.lock();
      cvs[currIdx].notify_all();
      lck.unlock();
      currIdx = (currIdx + 1) % members;
    }
  }

  lck.lock();
  maxIdx = 0;
  currIdx = 0;
  members = 0;
  overflow = 0;

  if (cvs) {
    delete[] cvs;
    cvs = NULL;
  }

  lck.unlock();
}


void
spin_semaphore::post() {
  unique_lock<mutex> lck(mtx);

  if (currIdx == maxIdx) {
    overflow++;
  } else {
    cvs[currIdx].notify_all();
    currIdx = (currIdx + 1) % members;
  }
}


void
spin_semaphore::wait() {
  unique_lock<mutex> lck(mtx);

  if (overflow > 0) {
    overflow--;
  } else {
    int myIdx = maxIdx;

    maxIdx = (maxIdx + 1) % members;
    cvs[myIdx].wait(lck);
  }
}


// Modificators

bool
Vertex::addPredecessor(Vertex * pred, int weight) {
  if (pred == this)
    return false;

  for (auto it = predecessors.begin(), end = predecessors.end(); it != end; ++it) {
    if (*it == pred)
      return false;
  }

  predecessors.push_back(pred);
  predecessorWeights.push_back(weight);
  predecessorCount++;
  return true;
}


bool
Vertex::addSuccessor(Vertex * succ, int weight) {
  if (succ == this)
    return false;

  for (auto it = successors_begin(), end = successors_end(); it != end; ++it) {
    if (*it == succ)
      return false;
  }

  successors.push_back(succ);
  successorWeights.push_back(weight);
  successorCount++;
  return true;
}


bool
Vertex::removePredecessor(Vertex * pred) {
  auto weightIt = predecessorWeights.begin();
  for (auto predIt = predecessors.begin(), end = predecessors.end(); predIt != end; ++predIt) {
    if (*predIt == pred) {
      predecessors.erase(predIt);
      predecessorWeights.erase(weightIt);
      predecessorCount--;
      return true;
    }

    ++weightIt;
  }

  return false;
}


bool
Vertex::removeSuccessor(Vertex * succ) {
  auto weightIt = successorWeights.begin();
  for (auto succIt = successors_begin(), end = successors_end(); succIt != end; ++succIt) {
    if (*succIt == succ) {
      successors.erase(succIt);
      successorWeights.erase(weightIt);
      successorCount--;
      return true;
    }

    ++weightIt;
  }

  return false;
}


void
Vertex::clearPredecessors() {
  predecessors.clear();
  predecessorWeights.clear();
}


void
Vertex::clearSuccessors() {
  successors.clear();
  successorWeights.clear();
}


// Access

Vertex *
Vertex::getPredecessor(int idx) const {
  return predecessors.at(idx);
}


Vertex *
Vertex::getSuccessor(int idx) const {
  return successors.at(idx);
}


int
Vertex::getPredecessorCount() const {
  return predecessorCount;
}


int
Vertex::getSuccessorCount() const {
  return successorCount;
}


int
Vertex::getPredecessorWeight(Vertex * pred) const {
  auto weightIt = predecessorWeights.begin();
  for (auto predIt = predecessors.begin(), end = predecessors.end(); predIt != end; ++predIt) {
    if (*predIt == pred)
      return *weightIt;

    ++weightIt;
  }

  return -1;
}


int
Vertex::getSuccessorWeight(Vertex * succ) const {
  auto weightIt = successorWeights.begin();
  for (auto succIt = successors.begin(), end = successors.end(); succIt != end; ++succIt) {
    if (*succIt == succ)
      return *weightIt;

    ++weightIt;
  }

  return -1;
}


int
Vertex::getPredecessorWeight(int idx) const {
  if (getPredecessorCount() <= idx)
    return -1;

  return predecessorWeights[idx];
}


int
Vertex::getSuccessorWeight(int idx) const {
  if (getSuccessorCount() <= idx)
    return -1;

  return successorWeights[idx];
}


// Iterators

Vertex::iterator
Vertex::predecessors_begin() {
  return predecessors.begin();
}


Vertex::iterator
Vertex::predecessors_end() {
  return predecessors.end();
}


Vertex::iterator
Vertex::successors_begin() {
  return successors.begin();
}


Vertex::iterator
Vertex::successors_end() {
  return successors.end();
}


// Indexing

#ifdef ENABLE_TLS
void
Vertex::setDFSId(uint64_t newId) {
  DFSId->at(id) = newId;
}
#else // ENABLE_TLS
void
Vertex::setDFSId(int idx, uint64_t newId) {
  DFSId[idx] = newId;
}
#endif // ENABLE_TLS


#ifdef ENABLE_TLS
uint64_t
Vertex::getDFSId() {
  return DFSId->at(id);
}
#else // ENABLE_TLS
uint64_t
Vertex::getDFSId(int idx) {
  return DFSId[idx];
}
#endif // ENABLE_TLS


/* Performs a "visit" of the Vertex when coming from the pred Vertex.
 * The boolean reverse indicates if the graph traversal is on the normal graph
 * or on the inverse graph.
 * The main effect is the reordering of the successor/predecessor list according
 * to the order in which the parent vertices access this Vertex instance during
 * graph traversal.
 */
void
Vertex::visit(Vertex * pred, int method) {
  // Reinitialize in case we are at the start of a new labeling traversal
  if ((inVisits + outVisits) == (predecessorCount + successorCount)) {
    inVisits = 0;
    outVisits = 0;
  }

  // Register the first visiting node
  if (inVisits == 0)
    firstVisit = pred;

  // If this is the first vertex to be visited (source vertex) then return
  if (!pred)
    return;

  // Clear the vector that should be reordered when requested
  if (((inVisits + outVisits) == 0) && (method & 0x02)) {
    if (method & 0x01)
      clearSuccessors();
    else
      clearPredecessors();
  }

  // When necessary push the visiting parent on the reordered vector
  if (method & 0x02) {
    if (method & 0x01) {
      successors.push_back(pred);
      successorWeights.push_back(pred->getPredecessorWeight(this));
    } else {
      predecessors.push_back(pred);
      predecessorWeights.push_back(pred->getSuccessorWeight(this));
    }
  }

  inVisits++;
}


/* Compute the next Vertex to be visited in a post-order traversal of the graph
 * from the current Vertex instance. This allows for an iterative traversal and
 * not a recursive one in order to prevent stack-overflows on large graphs
 */
Vertex *
Vertex::createPostOrder(vector<Vertex *> * postOrder, int method) {
  int visitIndex, visitNumber;
  vector<Vertex *> * nextVertices;

  // Determine what is up and what is down depending on the reverse boolean
  nextVertices = (method & 0x01) ? &predecessors : &successors;
  visitNumber = (method & 0x01) ? predecessorCount : successorCount;

  // Visit all the child vertices
  while (outVisits < visitNumber) {
    visitIndex = (method & 0x04) ? visitNumber - outVisits - 1 : outVisits;
    (*nextVertices)[visitIndex]->visit(this, method);

    outVisits++;

    if ((*nextVertices)[visitIndex]->inVisits == 1)
      return (*nextVertices)[visitIndex];
  }

  // Push the vertex on the stack only when all children have been visited
  postOrder->push_back(this);

  // Return the predecessor from which the first visit to this vertex was made
  return firstVisit;
}


// Modificators

void
Vertex::addPredecessorUnsafe(Vertex * pred, int weight) {
  predecessors.push_back(pred);
  predecessorWeights.push_back(weight);
}


void
Vertex::addSuccessorUnsafe(Vertex * succ, int weight) {
  successors.push_back(succ);
  successorWeights.push_back(weight);
}


// Parser

/* This function parses a subset of Dot files. Only edges and vertices can be
 * defined without any attributes
 */
Graph *
Graph::createFromDotFile(const char * fileName, bool noDoubleEdges) {
  char dump[128];
  int source, target, maxId;
  fstream input(fileName, fstream::in);
  Graph * graph = NULL;

#ifdef ENABLE_BENCHMARKS
  PAPI_get_dmem_info(&memoryInfo);
  baseMemoryUsage = memoryInfo.resident;
#endif // ENABLE_BENCHMARKS

  graph = new Graph();

  if (!input.good()) {
    cerr << "ERROR: Could not open the Dot input file.\n";
    exit(EXIT_FAILURE);
  }

  input.getline(dump, 127);
  if (!strstr(dump, "digraph")) {
    cerr << "ERROR: The supplied file is not a graph in Dot format.\n";
    exit(EXIT_FAILURE);
  }

  while (input.good()) {
    input.getline(dump, 127);

    if (strchr(dump, '}')) {
      break;
    }

    if (strchr(dump, '>')) {
      sscanf(dump, "%d -> %d", &source, &target);
      maxId = source < target ? target : source;
      while (graph->getVertexCount() <= (unsigned) maxId)
        graph->addVertexUnsafe(graph->threadCount);

      if (noDoubleEdges)
        graph->addEdgeUnsafe(graph->getVertexFromId(source), graph->getVertexFromId(target));
      else
        graph->addEdge(graph->getVertexFromId(source), graph->getVertexFromId(target));
    } else {
      sscanf(dump, "%d", &source);
      while (graph->getVertexCount() <= (unsigned) source)
        graph->addVertex();
    }
  }
  input.close();
  graph->indexed = false;

  // Make sure the graph is a DAG
  if (!graph->hasCycles())
    graph->condensed = true;
  else
    graph->condensed = false;

#ifdef ENABLE_BENCHMARKS
  PAPI_get_dmem_info(&memoryInfo);
  graph->graphMemoryUsage = memoryInfo.resident - baseMemoryUsage;
#endif // ENABLE_BENCHMARKS

  return graph;
}


Graph *
Graph::createFromGraFile(const char * fileName, bool noDoubleEdges) {
  char dump[128];
  int source, target, lineNumber;
  string barTitle = "Parsing graph ";
  fstream input(fileName, fstream::in);
  Graph *  graph = NULL;

#ifdef ENABLE_BENCHMARKS
  PAPI_get_dmem_info(&memoryInfo);
  baseMemoryUsage = memoryInfo.resident;
#endif // ENABLE_BENCHMARKS

  graph = new Graph();

  if (!input.good()) {
    cerr << "ERROR: Could not open the Gra input file.\n";
    exit(EXIT_FAILURE);
  }

  // Throw first lines away while they are text
  do {
    input.getline(dump, 127);
  } while (!isdigit(dump[0]));


  // Get the number of vertices / lines to read
  lineNumber = atoi(dump);
  configureProgressBar(&barTitle, lineNumber);

  // Create all vertices
  for (int i = 0; i < lineNumber; i++)
    graph->addVertexUnsafe(graph->threadCount);

  // Parse the adjacency list
  for (int i = 0; i < lineNumber; i++) {
    // Get the source at the start of the line
    input.get(dump, 127, ' ');
    source = atoi(dump);

    while (true) {
      if (input.get() != ' ') {
        cerr << "ERROR: Could not correctly read the graph definition.\n";
        exit(EXIT_FAILURE);
      }

      if (input.peek() == '#') {
        input.get(dump, 127);
        break;
      }

      input.get(dump, 127, ' ');
      target = atoi(dump);

      if (noDoubleEdges)
        graph->addEdgeUnsafe(graph->getVertexFromId(source), graph->getVertexFromId(target));
      else
        graph->addEdge(graph->getVertexFromId(source), graph->getVertexFromId(target));
    }

    resultProgressBar(i + 1);
  }
  cout << "\n";

  input.close();
  graph->indexed = false;

  // Make sure the graph is a DAG
  if (!graph->hasCycles())
    graph->condensed = true;
  else
    graph->condensed = false;

#ifdef ENABLE_BENCHMARKS
  PAPI_get_dmem_info(&memoryInfo);
  graph->graphMemoryUsage = memoryInfo.resident - baseMemoryUsage;
#endif // ENABLE_BENCHMARKS

  return graph;
}


// Serialization
void
Graph::printToFile(fstream &printStream, bool withWeights) {
  printStream << "graph_for_greach\n" << getVertexCount();
  for (auto nodeIt = vertices.begin(), end = vertices.end(); nodeIt != end; ++nodeIt) {
    if (*nodeIt == NULL)
      continue;

    printStream << "\n" << (*nodeIt)->id << ":";

    for (auto succIt = (*nodeIt)->successors_begin(), succEnd = (*nodeIt)->successors_end();
        succIt != succEnd; ++succIt)
      printStream << " " << (*succIt)->id;

    printStream << " #";

    if (!withWeights)
      continue;

    printStream << "\n" << (*nodeIt)->weight << ":";

    for (auto weightIt = (*nodeIt)->successorWeights.begin(),
        weightEnd = (*nodeIt)->successorWeights.end(); weightIt != weightEnd; ++weightIt)
      printStream << " " << (*weightIt);

    printStream << " #";
  }
}


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
    for (auto it = remappedVertices.begin(), end = remappedVertices.end(); it != end; ++it) {
      if (it->second == v->id)
        it->second = remapTarget;
    }
  }

  vertexCount--;

  vertices[v->id] = NULL;
  delete v;
}


void
Graph::mergeVertices(set<Vertex *> &s) {
  Vertex * target = addVertexUnsafe(threadCount);
  startGlobalOperation();

  int targetWeight = 0;
  map<Vertex *, int> predWeightMap;
  map<Vertex *, int> succWeightMap;

  for (auto mergeIt = s.begin(), mergeEnd = s.end(); mergeIt != mergeEnd; ++mergeIt) {
    auto weightIt = (*mergeIt)->predecessorWeights.begin();
    for (auto predIt = (*mergeIt)->predecessors_begin(), predEnd = (*mergeIt)->predecessors_end();
        predIt != predEnd; ++predIt) {
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
    for (auto succIt = (*mergeIt)->successors_begin(), succEnd = (*mergeIt)->successors_end();
        succIt != succEnd; ++succIt) {
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

    removeVertex(*mergeIt);
  }

  target->weight = targetWeight;

  for (auto mapIt = predWeightMap.begin(), mapEnd = predWeightMap.end(); mapIt != mapEnd; ++mapIt)
    addEdgeUnsafe(mapIt->first, target, mapIt->second);

  for (auto mapIt = succWeightMap.begin(), mapEnd = succWeightMap.end(); mapIt != mapEnd; ++mapIt)
    addEdgeUnsafe(target, mapIt->first, mapIt->second);

  indexed = false;
  stopGlobalOperation();
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


void
Graph::condenseGraph(bool dumpSCC) {
  int label = 0;
  stack<Vertex *> unclassified;
  stack<Vertex *> undetermined;
  vector<int> labels(vertices.size(), -1);

  // Erase any existing SCCs
  for (auto it = sccSets.begin(), end = sccSets.end(); it != end; ++it)
    delete *it;

  sccSets.clear();

  // Loop over all vertices and apply the Path-based SCC algorithm
  for (auto it = vertices.begin(), end = vertices.end(); it != end; ++it) {
    if ((*it == NULL) || (labels[(*it)->id] != -1))
      continue;

    label = sccVertexVisit(*it, label, labels, unclassified, undetermined);
  }

  // Verify that all valid vertices have been visited
  if (label != (int) getVertexCount())
    cerr << "Did not visit all vertices during graph condensation!" << endl;

  // When requested dump the SCCs that were found
  if (dumpSCC) {
    cerr << "DEBUG : SCC dump\n----" << endl;
    for (auto it = sccSets.begin(), end = sccSets.end(); it != end; ++it) {
      for (auto setIt = (*it)->begin(), setEnd = (*it)->end(); setIt != setEnd; ++setIt)
        cerr << (*setIt)->id << " ";

      cerr << "\n----" << endl;
    }
  }

  // Merge the SCCs that have more than one vertex
  for (auto it = sccSets.begin(), end = sccSets.end(); it != end; ++it)
    if ((*it)->size() > 1)
      mergeVertices(**it);
}



// Access

graph_t *
Graph::getMetisGraph() {
  int adjIdx = 0;
  graph_t * metisGraph = new graph_t;

  metisGraph->nvtxs = getVertexCount();
  metisGraph->xadj = new mtmetis_adj_t[metisGraph->nvtxs];
  metisGraph->adjncy = new mtmetis_vtx_t[getEdgeCount()];

  metisGraph->xadj[0] = adjIdx;

  for (auto it = vertices.begin(), end = vertices.end(); it != end; ++it) {
    if (*it == NULL)
      continue;

    for (auto it2 = (*it)->predecessors_begin(), end2 = (*it)->predecessors_end();
        it2 != end2; ++it2)
      metisGraph->adjncy[adjIdx++] = (*it)->id;

    for (auto it2 = (*it)->successors_begin(), end2 = (*it)->successors_end();
        it2 != end2; ++it2)
      metisGraph->adjncy[adjIdx++] = (*it)->id;

    metisGraph->xadj[(*it)->id + 1] = adjIdx;
  }

  return metisGraph;
}


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

list<Vertex *> *
Graph::hasCycles() {
  list<Vertex *> * retValue = new list<Vertex *>();
  vector<Vertex *> * cyclePath = NULL;
  globalTimestamp =
    chrono::high_resolution_clock::now().time_since_epoch() / chrono::nanoseconds(1);

  startGlobalOperation();

  discoverExtremities();

  // Iterate over the sources and check
  for (auto it = sources.begin(); it != sources.end(); ++it) {
    if ((cyclePath = hasCyclesFromSource(*it)))
      break;
  }

  if (cyclePath) {
    Vertex * start = cyclePath->back();
    cyclePath->pop_back();

    do {
      retValue->push_front(cyclePath->back());
      cyclePath->pop_back();
    } while (retValue->front() != start);

    delete cyclePath;

    stopGlobalOperation();
    return retValue;
  }

  delete retValue;

  stopGlobalOperation();

  return NULL;
}


// Reachability queries

void
Graph::pushQuery(Query * query) {
  static int DFSwin = 0;
  static int BBFSwin = 0;

  unique_lock<mutex> methodLock(methodMutex, defer_lock);
  unique_lock<mutex> internalResultLock(internalResultMutex, defer_lock);

  // Verify that the graph has been indexed
  if (!indexed) {
    indexGraph();
    methodLock.lock();
    setMethod(UndefinedSearchMethod);
    DFSwin = 0;
    BBFSwin = 0;
    methodLock.unlock();
  }

  // If a method is specified use it
  if ((query->type != Reachability) ||
      (query->query.reachability->getMethod() != UndefinedSearchMethod)) {
    jobQueue.push(query);
    return;
  }

  // Else select the method automatically. If there is no
  // preferred method yet then find one.
  if (getMethod() == UndefinedSearchMethod) {
    int results = 0;
    ReachabilityQuery * reachQuery = query->query.reachability;
    ReachabilityQuery * DFSReachQuery =
      new ReachabilityQuery(reachQuery->getSource(), reachQuery->getTarget(), DFS);
    ReachabilityQuery * BBFSReachQuery =
      new ReachabilityQuery(reachQuery->getSource(), reachQuery->getTarget(), BBFS);

    DFSReachQuery->setInternal(1);
    BBFSReachQuery->setInternal(1);

    Query * DFSQuery = new Query(Reachability, DFSReachQuery);
    Query * BBFSQuery = new Query(Reachability, BBFSReachQuery);

    jobQueue.push(DFSQuery);
    jobQueue.push(BBFSQuery);

    // Wait for the results
    while (results < 2) {
      internalResultSpinSemaphore.wait();
      internalResultLock.lock();

      auto it = internalResultQueue.begin();
      for (auto end = internalResultQueue.end(); it != end; ++it) {
        if (((*it) == DFSQuery) || ((*it) == BBFSQuery)) {
          results++;
          internalResultQueue.erase(it);
          break;
        }
      }

      if (it == internalResultQueue.end())
        internalResultSpinSemaphore.post();

      internalResultLock.unlock();
    }

    // Check for errors in the queries
    if (DFSReachQuery->getError() || BBFSReachQuery->getError()) {
      cerr << "ERROR: Could not correctly process BBFS or DFS search." << endl;
      exit(EXIT_FAILURE);
    }

    // Check for coherency between queries
    if (DFSReachQuery->getAnswer() != BBFSReachQuery->getAnswer()) {
      cerr << "ERROR: Incoherent results between BBFS and DFS search." << endl;
      exit(EXIT_FAILURE);
    }

    // Find the fastest search method and register the win
    methodLock.lock();
    if (DFSReachQuery->getInternal() < BBFSReachQuery->getInternal()) {
      if ((++DFSwin >= AUTO_THRESHOLD) && (getMethod() == UndefinedSearchMethod))
        setMethod(DFS);

      query->query.reachability->setAnswer(DFSReachQuery->getAnswer());
      query->query.reachability->setMethod(DFS);

#ifdef ENABLE_STATISTICS
      registerQueryStatistics(DFSReachQuery);
#endif // ENABLE_STATISTICS
    } else {
      if ((++BBFSwin >= AUTO_THRESHOLD) && (getMethod() == UndefinedSearchMethod))
        setMethod(BBFS);

      query->query.reachability->setAnswer(BBFSReachQuery->getAnswer());
      query->query.reachability->setMethod(BBFS);

#ifdef ENABLE_STATISTICS
      registerQueryStatistics(BBFSReachQuery);
#endif // ENABLE_STATISTICS
    }
    methodLock.unlock();

    // Clean-up the internal queries
    delete DFSQuery;
    delete BBFSQuery;

    resultQueue.push(query);
  } else {
    query->query.reachability->setMethod(getMethod());

    jobQueue.push(query);
  }

  return;
}


Query *
Graph::pullResult() {
  Query * result = NULL;

  while (!resultQueue.try_pop(result)) {}
    // Do nothing

  return result;
}


void
Graph::endOfQueries() {
  threadShutdown = true;

  for (int i = 0; i < threadCount; i++)
    jobQueue.push(NULL);
}


// Partitioning queries
int
Graph::getPartitionCount() {
  return 1;
}


Graph::PartitionSet *
getPartitionSet(int count) {
  return NULL;
}


// Graph coarsening

Graph *
Graph::coarsen(CoarsenMethod method, int factor, int secondaryFactor, map<int, int> &vertexMap) {
  Graph * coarsenedGraph = NULL;

  indexGraph();

  startGlobalOperation();

  switch (method) {
    case Greedy:
      coarsenedGraph = coarsenGreedy(factor, vertexMap);
      break;

    case EdgeRedux:
      coarsenedGraph = coarsenEdgeRedux(factor, vertexMap);
      break;

    case ApproxIteration:
      coarsenedGraph = coarsenApproxIteration(factor, secondaryFactor, vertexMap);
      break;

    case UndefinedCoarsenMethod:
      cerr << "ERROR: Called coarsening method with undefined method" << endl;
      exit(EXIT_FAILURE);
  }

  stopGlobalOperation();

  return coarsenedGraph;
}


// Benchmark statistics

bool
Graph::statisticsAreEnabled() {
#ifdef ENABLE_STATISTICS
  return true;
#else // ENABLE_STATISTICS
  return false;
#endif // ENABLE_STATISTICS
}


uintmax_t
Graph::getQueryCount() {
#ifdef ENABLE_STATISTICS
  return queryCount;
#else // ENABLE_STATISTICS
  return 0;
#endif // ENABLE_STATISTICS
}


uintmax_t
Graph::getPositiveQueryCount() {
#ifdef ENABLE_STATISTICS
  return positiveQueryCount;
#else // ENABLE_STATISTICS
  return 0;
#endif // ENABLE_STATISTICS
}


uintmax_t
Graph::getNegativeQueryCount() {
#ifdef ENABLE_STATISTICS
  return negativeQueryCount;
#else // ENABLE_STATISTICS
  return 0;
#endif // ENABLE_STATISTICS
}


uintmax_t
Graph::getShortNegativeQueryCount() {
#ifdef ENABLE_STATISTICS
  return shortNegativeQueryCount;
#else // ENABLE_STATISTICS
  return 0;
#endif // ENABLE_STATISTICS
}


double
Graph::getPositiveQueryOverhead() {
#ifdef ENABLE_STATISTICS
  return positiveQueryOverhead;
#else // ENABLE_STATISTICS
  return 0.0L;
#endif // ENABLE_STATISTICS
}


double
Graph::getNegativeQueryOverhead() {
#ifdef ENABLE_STATISTICS
  return negativeQueryOverhead;
#else // ENABLE_STATISTICS
  return 0.0L;
#endif // ENABLE_STATISTICS
}


void
Graph::printStatistics(ostream &os) {
#ifdef ENABLE_STATISTICS
  double shortFraction = ((double) shortNegativeQueryCount / ((double) negativeQueryCount));
  os << "\n---\nStatistics:\n";
  os << "General statistics\n";
  os << "Number of vertices: " << getVertexCount() << "\n";
  os << "Number of edges: " << edgeCount << "\n";
  os << "Number of performed queries: " << queryCount << "\n";

  os << "\nPositive query statistics:\n";
  os << "- Number of positive answers: " << positiveQueryCount << "\n";
  if (positiveQueryCount && (getMethod() == DFS)) {
    os.precision(3);
    os << "- Average DFS length / path-length ratio: " << positiveQueryOverhead << "\n";
  }

  os << "\nNegative query statistics:\n";
  os << "- Number of negative answers   : " << negativeQueryCount << "\n";
  if (negativeQueryCount) {
    os << "- Number of immediate negatives: " << shortNegativeQueryCount;
    os.precision(4);
    os << " (" << shortFraction * 100 << "%).\n";
    if (getMethod() == DFS) {
      os.precision(3);
      os << "- Average DFS length / graph-size ratio: " << negativeQueryOverhead << "\n";
    }
  }
  os << "---\n\n";
#else // ENABLE_STATISTICS
  os << "WARNING: Statistics gathering has not been enabled at compile time.\n";
  os << "WARNING: To enable statistics invoke cmake with '-DENABLE_STATISTICS=1'.\n\n";
#endif // ENABLE_STATISTICS
}


// Benchmarking

bool
Graph::benchmarksAreEnabled(void) {
#ifdef ENABLE_BENCHMARKS
  return true;
#else // ENABLE_BENCHMARKS
  return false;
#endif // ENABLE_BENCHMARKS
}


long long
Graph::getQueryNumber(void) {
#ifdef ENABLE_BENCHMARKS
  return queryNumber;
#else // ENABLE_BENCHMARKS
  return 0;
#endif // ENABLE_BENCHMARKS
}


long long
Graph::getCyclesSpentIndexing(void) {
#ifdef ENABLE_BENCHMARKS
  return cyclesSpentIndexing;
#else // ENABLE_BENCHMARKS
  return 0;
#endif // ENABLE_BENCHMARKS
}


long long
Graph::getCyclesSpentQuerying(void) {
#ifdef ENABLE_BENCHMARKS
  return cyclesSpentQuerying;
#else // ENABLE_BENCHMARKS
  return 0;
#endif // ENABLE_BENCHMARKS
}


void
Graph::printBenchmarks(ostream &os) {
#ifdef ENABLE_BENCHMARKS
  int indexMemoryUsage = (8 * getVertexCount()) / 1024;
  double indexOverhead = ((double) indexMemoryUsage) / ((double) graphMemoryUsage - indexMemoryUsage);
  os << "\n---\nBenchmarking\n";
  os << "Speed performances:\n";
  os << "Cycles spent indexing: " << cyclesSpentIndexing << "\n";
  os << "Cycles spent querying: " << cyclesSpentQuerying << "\n";
  os << "-> average of " << cyclesSpentQuerying / queryNumber << " cycles per query\n\n";
  os << "Memory usage: " << graphMemoryUsage << " kilo-bytes.\n";
  os << "Index memory requirement: " << indexMemoryUsage << " kilo-bytes\n";
  os.precision(2);
  os << "Index memory overhead: " << indexOverhead * 100 << "%\n---\n\n";
#else // ENABLE_BENCHMARKS
  os << "WARNING: Benchmarking has not been enabled at compile time.\n";
  os << "WARNING: To enable benchmarks invoke cmake with '-DENABLE_BENCHMARKS=1'.\n\n";
#endif // ENABLE_BENCHMARKS
}


// Indexing

void
Graph::labelVertices(bool retro, bool reverse) {
  int traversalMethod = 0;
  int currLabel = 0;
  Vertex * nextVertex;
  vector<Vertex *> * startVertices;;
  vector<Vertex *> * postOrder = reverse ? &predecessorQueue : &successorQueue;

#ifndef ENABLE_RETRO_LABELS
  // In case retro labels have not been compiled we should abort
  if (retro)
    return;
#endif // ENABLE_RETRO_LABELS

  startVertices = reverse ? &sinks : &sources;

  // Compose the method used for post-order labeling
  if (reverse)
    traversalMethod |= 0x01;

  if (retro)
    traversalMethod |= 0x04;
  else
    traversalMethod |= (indexMethod & 0x02);

  // Loop while there are unlabeled sources / sinks
  for (auto it = startVertices->begin(), end = startVertices->end(); it != end; ++it) {
    // Initialize
    nextVertex = *it;
    (*it)->visit(NULL, traversalMethod);

    // Run the DFS and use an iterative method to prevent memory overflow in large graphs
    while (nextVertex)
      nextVertex = nextVertex->createPostOrder(postOrder, traversalMethod);
  }

  // Label the vertices in reverse post-order
  for (auto it = postOrder->rbegin(), end = postOrder->rend(); it != end; ++it) {
    if (retro) {
#ifdef ENABLE_RETRO_LABELS
      if (reverse)
        (*it)->retroReverseOrderLabel = currLabel++;
      else
        (*it)->retroOrderLabel = currLabel++;;
#endif // ENABLE_RETRO_LABELS
    } else {
      if (reverse)
        (*it)->reverseOrderLabel = currLabel++;
      else
        (*it)->orderLabel = currLabel++;
    }
  }
}


void
Graph::indexGraph() {
  unique_lock<mutex> indexLock(indexMutex);

  startGlobalOperation();

  if (indexed)
    return;

  // Make sure we are indexing a DAG
  if (!condensed)
    condenseGraph();

  if (indexMethod == UndefinedIndexMethod) {
    cerr << "ERROR: Unknown indexing method. Aborting.\n";
    exit(EXIT_FAILURE);
  }

#ifdef ENABLE_BENCHMARKS
  int event = PAPI_TOT_CYC;
  PAPI_start_counters(&event, 1);
#endif // ENABLE_BENCHMARKS

  // First traversal
  labelVertices(false, false);

  // Second traversal
  labelVertices(false, true);

  // When requested reorder the successors and predecessors in all vertices
  if (indexMethod & 0x04) {
    Vertex * curr;

    while (!successorQueue.empty()) {
      curr = successorQueue.back();
      successorQueue.pop_back();

      for (auto it = curr->predecessors_begin(), end = curr->predecessors_end(); it != end; ++it)
        (*it)->addSuccessorUnsafe(curr, curr->getPredecessorWeight(*it));

      curr->clearSuccessors();
    }

    while (!predecessorQueue.empty()) {
      curr = predecessorQueue.back();
      predecessorQueue.pop_back();

      for (auto it = curr->successors_begin(), end = curr->successors_end(); it != end; ++it)
        (*it)->addPredecessorUnsafe(curr, curr->getSuccessorWeight(*it));

      curr->clearPredecessors();
    }
  } else {
    successorQueue.clear();
    predecessorQueue.clear();
  }
#ifdef ENABLE_RETRO_LABELS
  // Third traversal
  labelVertices(true, false);

  // Fourth traversal
  labelVertices(true, true);

  successorQueue.clear();
  predecessorQueue.clear();
#endif // ENABLE_RETRO_LABELS

#ifdef ENABLE_BENCHMARKS
  long long indexCycles;
  PAPI_stop_counters(&indexCycles, 1);
  cyclesSpentIndexing = indexCycles;
#endif // ENABLE_BENCHMARKS

  indexed = true;

  stopGlobalOperation();
}


// Maintenance

void
Graph::discoverExtremities() {
  sources.clear();
  sinks.clear();

  for (auto it = vertices.begin(), end = vertices.end(); it != end; ++it) {
    if (*it == NULL)
      continue;

    if ((*it)->predecessors.empty())
      sources.push_back(*it);

    if ((*it)->successors.empty())
      sinks.push_back(*it);
  }
}


// Multi-threading
void
Graph::startQuery() {
  unique_lock<mutex> queryWaitLock(queryWaitMutex);

  if (noQueries)
    queryWaitCondition.wait(queryWaitLock);

  activeThreads++;
}


void
Graph::stopQuery() {
  unique_lock<mutex> queryWaitLock(queryWaitMutex, defer_lock);
  unique_lock<mutex> globalWaitLock(globalWaitMutex, defer_lock);

  queryWaitLock.lock();
  activeThreads--;
  queryWaitLock.unlock();

  if (noQueries && (activeThreads == 0)) {
    globalWaitLock.lock();
    globalWaitCondition.notify_one();
    globalWaitLock.unlock();
  }
}


void
Graph::startGlobalOperation() {
  unique_lock<mutex> queryWaitLock(queryWaitMutex, defer_lock);
  unique_lock<mutex> globalWaitLock(globalWaitMutex, defer_lock);

  queryWaitLock.lock();
  noQueries = true;
  queryWaitLock.unlock();

  globalWaitLock.lock();
  if (activeThreads > 0)
    globalWaitCondition.wait(globalWaitLock);

  globalWaitLock.unlock();

#ifdef ENABLE_TLS
  if (!DFSId)
    DFSId = new vector<uint64_t>();

  if (DFSId->size() < vertices.size())
    DFSId->resize(vertices.size(), 0);
#endif // ENABLE_TLS
}


void
Graph::stopGlobalOperation() {
  unique_lock<mutex> queryWaitLock(queryWaitMutex);

#ifdef ENABLE_TLS
  if (!DFSId)
    DFSId = new vector<uint64_t>();

  if (DFSId->size() < vertices.size())
    DFSId->resize(vertices.size(), 0);
#endif // ENABLE_TLS

  noQueries = false;

  queryWaitCondition.notify_all();
}


void
Graph::startWorkers() {
  if (threadsActive)
    return;

  for (int i = 0; i < threadCount; i++)
    queryThreads.push_back(new thread(queryWorker, this, i));

  threadsActive = true;
}


void
Graph::enableGraph(int numberOfThreads) {
  if ((enabled) || (numberOfThreads == 0))
    return;

  threadCount = numberOfThreads;
  PAPI_library_init(PAPI_VER_CURRENT);
  internalResultSpinSemaphore.initialize(numberOfThreads);

  startWorkers();

  enabled = true;
}


void
Graph::disableGraph() {
  if (threadCount > 0) {
    endOfQueries();

    while (!queryThreads.empty()) {
      queryThreads.back()->join();

      delete queryThreads.back();
      queryThreads.pop_back();
    }

    internalResultSpinSemaphore.reset();

    PAPI_shutdown();

    threadCount = 0;
  }
}


// Maintenance

vector<Vertex *> *
Graph::hasCyclesFromSource(Vertex * source) {
  set<Vertex *> DFSSet;
  vector<Vertex *> * DFSPath = new vector<Vertex *>();
  stack<Vertex *> toVisit;
  Vertex * curr = NULL;

  toVisit.push(source);
#ifdef ENABLE_TLS
  source->setDFSId(globalTimestamp);
#else // ENABLE_TLS
  source->setDFSId(0, globalTimestamp);
#endif // ENABLE_TLS

  while (!toVisit.empty()) {
    curr = toVisit.top();
    if (!DFSPath->empty() && (curr == DFSPath->back())) {
      DFSPath->pop_back();
      DFSSet.erase(curr);
      toVisit.pop();
      continue;
    }

    DFSPath->push_back(curr);
    DFSSet.insert(curr);

    for (auto it = curr->successors.begin(), end = curr->successors.end(); it != end; ++it) {
      if (DFSSet.count(*it)) {
        DFSPath->push_back(*it);
        return DFSPath;
      }

#ifdef ENABLE_TLS
      if ((*it)->getDFSId() != globalTimestamp) {
        (*it)->setDFSId(globalTimestamp);
#else // ENABLE_TLS
      if ((*it)->getDFSId(0) != globalTimestamp) {
        (*it)->setDFSId(0, globalTimestamp);
#endif // ENABLE_TLS
        toVisit.push(*it);
      }
    }
  }

  delete DFSPath;

  return NULL;
}


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


// Internal condense functions

int
Graph::sccVertexVisit(Vertex * target, int label, vector<int> &labels,
    stack<Vertex *> &unclassified, stack<Vertex *> &undetermined) {
  if (labels[target->id] != -1)
    cerr << "Revisited node " << target->id << " during graph condensation!" << endl;

  labels[target->id] = label++;
  unclassified.push(target);
  undetermined.push(target);

  for (auto it = target->successors_begin(), end = target->successors_end(); it != end; ++it) {
    if (labels[(*it)->id] == -1) {
      label = sccVertexVisit(*it, label, labels, unclassified, undetermined);
    } else {
      while (labels[undetermined.top()->id] > labels[(*it)->id])
        undetermined.pop();
    }
  }

  if (target == undetermined.top()) {
    sccSets.push_back(new set<Vertex *>);

    while (unclassified.top() != target) {
      sccSets.back()->insert(unclassified.top());
      unclassified.pop();
    }

    sccSets.back()->insert(unclassified.top());
    unclassified.pop();

    undetermined.pop();
  }

  return label;
}


// Internal reachability query functions

void
#ifdef ENABLE_TLS
Graph::processReachabilityQuery(Query * query) {
  unique_lock<mutex> internalResultLock(internalResultMutex, defer_lock);

#ifdef ENABLE_TLS
  if (!DFSId)
    DFSId = new vector<uint64_t>();

  if (DFSId->size() < vertices.size())
    DFSId->resize(vertices.size(), 0);
#endif // ENABLE_TLS

  switch (query->query.reachability->getMethod()) {
    case DFS:
      areConnectedDFS(query->query.reachability);
      break;

    case BBFS:
      areConnectedBBFS(query->query.reachability);
      break;

    case Approx:
      areConnectedApprox(query->query.reachability);
      break;

    case NoLabels:
      areConnectedNoLabels(query->query.reachability);
      break;
#else // ENABLE_TLS
Graph::processReachabilityQuery(int threadId, Query * query) {
  unique_lock<mutex> internalResultLock(internalResultMutex, defer_lock);

  switch (query->query.reachability->getMethod()) {
    case DFS:
      areConnectedDFS(query->query.reachability, threadId);
      break;

    case BBFS:
      areConnectedBBFS(query->query.reachability, threadId);
      break;

    case Approx:
      areConnectedApprox(query->query.reachability, threadId);
      break;

    case NoLabels:
      areConnectedNoLabels(query->query.reachability, threadId);
      break;
#endif // ENABLE_TLS

    case UndefinedSearchMethod:
      query->query.reachability->setError(true);
      break;
  }

  if (query->query.reachability->getInternal() != 0) {
    internalResultLock.lock();
    internalResultQueue.push_back(query);
    internalResultSpinSemaphore.post();
    internalResultLock.unlock();
  } else {
#ifdef ENABLE_STATISTICS
    registerQueryStatistics(query->query.reachability);
#endif // ENABLE_STATISTICS

    resultQueue.push(query);
  }
}


/* Use the previously done indexation to answer to the query */
void
#ifdef ENABLE_TLS
Graph::areConnectedDFS(ReachabilityQuery * query) {
#else // ENABLE_TLS
Graph::areConnectedDFS(ReachabilityQuery * query, int threadId) {
#endif // ENABLE_TLS
  int event = PAPI_TOT_CYC;
  long long counterValue;
  uint64_t timestamp;
  Vertex * curr;
  vector<Vertex *> searchStack;

#ifdef ENABLE_BENCHMARKS
  unique_lock<mutex> benchmarkLock(benchmarkMutex, defer_lock);
  PAPI_start_counters(&event, 1);
#else // ENABLE_BENCHMARKS
  if (query->getInternal() == 1)
    PAPI_start_counters(&event, 1);
#endif // ENABLE_BENCHMARKS

  query->setAnswer(false);

  // Are U and V the same vertex?
  if (query->getSource() == query->getTarget()) {
#ifdef ENABLE_STATISTICS
    query->searchedNodes++;
    query->path.push_back(query->getSource());
#endif // ENABLE_STATISTICS
    query->setAnswer(true);
    goto end;
  }

  // Can V be a descendant of U in the standard graph or U a descendant of V in
  // the reverse graph?
  if ((query->getSource()->orderLabel > query->getTarget()->orderLabel) ||
      (query->getTarget()->reverseOrderLabel > query->getSource()->reverseOrderLabel))
    goto end;

#ifdef ENABLE_RETRO_LABELS
  if ((query->getSource()->retroOrderLabel > query->getTarget()->retroOrderLabel) ||
      (query->getTarget()->retroReverseOrderLabel > query->getSource()->retroReverseOrderLabel))
    goto end;
#endif // ENABLE_RETRO_LABELS

  // Do a DFS on the subgraph specified by both orders to get the final answer
  timestamp = chrono::high_resolution_clock::now().time_since_epoch() / chrono::nanoseconds(1);

  searchStack.push_back(query->getSource());

  while (!searchStack.empty()) {
    if (query->getCancel()) {
      query->setError(true);
      goto cancel;
    }

    curr = searchStack.back();

#ifdef ENABLE_STATISTICS
    if (!query->path.empty() && (curr == query->path.back())) {
      query->path.pop_back();
      searchStack.pop_back();
      continue;
    }

    query->searchedNodes++;
    query->path.push_back(curr);
#else // ENABLE_STATISTICS
    searchStack.pop_back();
#endif // ENABLE_STATISTICS

    for (auto it = curr->successors.begin(), end = curr->successors.end(); it != end; ++it) {
      if ((indexMethod & 0x04) && ((*it)->orderLabel > query->getTarget()->orderLabel))
        break;

      if (*it == query->getTarget()) {
#ifdef ENABLE_STATISTICS
        query->searchedNodes++;
        query->path.push_back(query->getTarget());
#endif // ENABLE_STATISTICS

        query->setAnswer(true);
        goto end;
      }

#ifdef ENABLE_TLS
      if ((*it)->getDFSId() != timestamp) {
        (*it)->setDFSId(timestamp);
#else // ENABLE_TLS
      if ((*it)->getDFSId(threadId) != timestamp) {
        (*it)->setDFSId(threadId, timestamp);
#endif // ENABLE_TLS

#ifdef ENABLE_RETRO_LABELS
        if (((*it)->orderLabel < query->getTarget()->orderLabel) &&
            ((*it)->reverseOrderLabel > query->getTarget()->reverseOrderLabel) &&
            ((*it)->retroOrderLabel < query->getTarget()->retroOrderLabel) &&
            ((*it)->retroReverseOrderLabel > query->getTarget()->retroReverseOrderLabel))
#else // ENABLE_RETRO_LABELS
        if (((*it)->orderLabel < query->getTarget()->orderLabel) &&
            ((*it)->reverseOrderLabel > query->getTarget()->reverseOrderLabel))
#endif // ENABLE_RETRO_LABELS
          searchStack.push_back(*it);
      }
    }
  }

end:
#ifdef ENABLE_BENCHMARKS
  PAPI_stop_counters(&counterValue, 1);
  benchmarkLock.lock();
  cyclesSpentQuerying += counterValue;
  query->umber++;
  benchmarkLock.unlock();

  // Set the query->time
  if (query->getInternal() == 1)
    query->setInternal(counterValue);
#else // ENABLE_BENCHMARKS
  // In the case of an internal query->register the query->time
  if (query->getInternal() == 1) {
    PAPI_stop_counters(&counterValue, 1);
    query->setInternal(counterValue);
  }
#endif // ENABLE_BENCHMARKS

  return;

cancel:
#ifdef ENABLE_BENCHMARKS
  PAPI_stop_counters(&counterValue, 1);
#endif // ENABLE_BENCHMARKS

  query->setError(true);
  return;
}


void
#ifdef ENABLE_TLS
Graph::areConnectedBBFS(ReachabilityQuery * query) {
#else // ENABLE_TLS
Graph::areConnectedBBFS(ReachabilityQuery * query, int threadId) {
#endif
  int event = PAPI_TOT_CYC;
  long long counterValue;
  uint64_t forwardId, backwardId;
  Vertex * curr;
  queue<Vertex *> searchQueueForward;
  queue<Vertex *> searchQueueBackward;

#ifdef ENABLE_BENCHMARKS
  unique_lock<mutex> benchmarkLock(benchmarkMutex, defer_lock);
  PAPI_start_counters(&event, 1);
#else // ENABLE_BENCHMARKS
  if (query->getInternal() == 1)
    PAPI_start_counters(&event, 1);
#endif // ENABLE_BENCHMARKS

  query->setAnswer(false);

  // Are U and V the same vertex?
  if (query->getSource() == query->getTarget()) {
#ifdef ENABLE_STATISTICS
    query->searchedNodes++;
#endif // ENABLE_STATISTICS
    query->setAnswer(true);
    goto end;
  }

  // Can V be a descendant of U in the standard graph or U a descendant of V in
  // the reverse graph?
  if ((query->getSource()->orderLabel > query->getTarget()->orderLabel) ||
      (query->getTarget()->reverseOrderLabel > query->getSource()->reverseOrderLabel))
    goto end;

#ifdef ENABLE_RETRO_LABELS
  if ((query->getSource()->retroOrderLabel > query->getTarget()->retroOrderLabel) ||
      (query->getTarget()->retroReverseOrderLabel > query->getSource()->retroReverseOrderLabel))
    goto end;
#endif // ENABLE_RETRO_LABELS

  // Do a BBFS on the subgraph specified by both orders to get the final answer
  forwardId = chrono::high_resolution_clock::now().time_since_epoch() / chrono::nanoseconds(1);
  backwardId = forwardId + 1;

#ifdef ENABLE_TLS
  query->getSource()->setDFSId(forwardId);
  query->getTarget()->setDFSId(backwardId);
#else // ENABLE_TLS
  query->getSource()->setDFSId(threadId, forwardId);
  query->getTarget()->setDFSId(threadId, backwardId);
#endif // ENABLE_TLS

  searchQueueForward.push(query->getSource());
  searchQueueBackward.push(query->getTarget());

  while (!searchQueueForward.empty() || !searchQueueBackward.empty()) {
    if (query->getCancel()) {
      query->setError(true);
      goto cancel;
    }

    if (!searchQueueForward.empty()) {
      curr = searchQueueForward.front();
      searchQueueForward.pop();

#ifdef ENABLE_STATISTICS
      query->searchedNodes++;
#endif // ENABLE_STATISTICS

      for (auto it = curr->successors.begin(), end = curr->successors.end(); it != end; ++it) {
        if ((indexMethod & 0x04) && ((*it)->orderLabel > query->getTarget()->orderLabel))
          break;

#ifdef ENABLE_TLS
        if ((*it)->getDFSId() == backwardId) {
#else // ENABLE_TLS
        if ((*it)->getDFSId(threadId) == backwardId) {
#endif // ENABLE_TLS
          query->setAnswer(true);
          goto end;
        }

#ifdef ENABLE_TLS
        if ((*it)->getDFSId() != forwardId) {
          (*it)->setDFSId(forwardId);
#else // ENABLE_TLS
        if ((*it)->getDFSId(threadId) != forwardId) {
          (*it)->setDFSId(threadId, forwardId);
#endif

#ifdef ENABLE_RETRO_LABELS
          if (((*it)->orderLabel < query->getTarget()->orderLabel) &&
              ((*it)->reverseOrderLabel > query->getTarget()->reverseOrderLabel) &&
              ((*it)->retroOrderLabel < query->getTarget()->retroOrderLabel) &&
              ((*it)->retroReverseOrderLabel > query->getTarget()->retroReverseOrderLabel))
#else // ENABLE_RETRO_LABELS
          if (((*it)->orderLabel < query->getTarget()->orderLabel) &&
              ((*it)->reverseOrderLabel > query->getTarget()->reverseOrderLabel))
#endif // ENABLE_RETRO_LABELS
            searchQueueForward.push(*it);
        }
      }
    }

    if (!searchQueueBackward.empty()) {
      curr = searchQueueBackward.front();
      searchQueueBackward.pop();

#ifdef ENABLE_STATISTICS
      query->searchedNodes++;
#endif // ENABLE_STATISTICS

      for (auto it = curr->predecessors.begin(), end = curr->predecessors.end(); it != end; ++it) {
        if ((indexMethod & 0x04) && ((*it)->orderLabel < query->getSource()->orderLabel))
          break;

#ifdef ENABLE_TLS
        if ((*it)->getDFSId() == forwardId) {
#else // ENABLE_TLS
        if ((*it)->getDFSId(threadId) == forwardId) {
#endif // ENABLE_TLS
          query->setAnswer(true);
          goto end;
        }

#ifdef ENABLE_TLS
        if ((*it)->getDFSId() != backwardId) {
          (*it)->setDFSId(backwardId);
#else // ENABLE_TLS
        if ((*it)->getDFSId(threadId) != backwardId) {
          (*it)->setDFSId(threadId, backwardId);
#endif // ENABLE_TLS

#ifdef ENABLE_RETRO_LABELS
          if (((*it)->orderLabel > query->getSource()->orderLabel) &&
              ((*it)->reverseOrderLabel < query->getSource()->reverseOrderLabel) &&
              ((*it)->retroOrderLabel > query->getSource()->retroOrderLabel) &&
              ((*it)->retroReverseOrderLabel < query->getSource()->retroReverseOrderLabel))
#else // ENABLE_RETRO_LABELS
          if (((*it)->orderLabel > query->getSource()->orderLabel) &&
              ((*it)->reverseOrderLabel < query->getSource()->reverseOrderLabel))
#endif // ENABLE_RETRO_LABELS
            searchQueueBackward.push(*it);
        }
      }
    }
  }

end:
#ifdef ENABLE_BENCHMARKS
  PAPI_stop_counters(&counterValue, 1);
  benchmarkLock.lock();
  cyclesSpentQuerying += counterValue;
  query->umber++;
  benchmarkLock.unlock();

  // Set the query time
  if (query->getInternal() == 1)
    query->setInternal(counterValue);
#else // ENABLE_BENCHMARKS
  // In the case of an internal query register the query time
  if (query->getInternal() == 1) {
    PAPI_stop_counters(&counterValue, 1);
    query->setInternal(counterValue);
  }
#endif // ENABLE_BENCHMARKS

  return;

cancel:
#ifdef ENABLE_BENCHMARKS
  PAPI_stop_counters(&counterValue, 1);
#endif // ENABLE_BENCHMARKS

  query->setError(true);
  return;
}


void
#ifdef ENABLE_TLS
Graph::areConnectedApprox(ReachabilityQuery * query) {
#else // ENABLE_TLS
Graph::areConnectedApprox(ReachabilityQuery * query, int threadId) {
#endif // ENABLE_TLS
  if (query->getSource() == query->getTarget()) {
#ifdef ENABLE_STATISTICS
    query->searchedNodes++;
#endif // ENABLE_STATISTICS
    query->setAnswer(true);
    return;
  }

  if ((query->getSource()->orderLabel > query->getTarget()->orderLabel) ||
      (query->getTarget()->reverseOrderLabel > query->getSource()->reverseOrderLabel)) {
    query->setAnswer(false);
    return;
  }

#ifdef ENABLE_RETRO_LABELS
  if ((query->getSource()->retroOrderLabel > query->getTarget()->retroOrderLabel) ||
      (query->getTarget()->retroReverseOrderLabel > query->getSource()->retroReverseOrderLabel)) {
    query->setAnswer(false);
    return;
  }
#endif // ENABLE_RETRO_LABELS

  query->setAnswer(true);
}


void
#ifdef ENABLE_TLS
Graph::areConnectedNoLabels(ReachabilityQuery * query) {
#else // ENABLE_TLS
Graph::areConnectedNoLabels(ReachabilityQuery * query, int threadId) {
#endif // ENABLE_TLS
  uint64_t timestamp;
  stack<Vertex *> toVisit;
  Vertex * curr;

  timestamp = chrono::high_resolution_clock::now().time_since_epoch() / chrono::nanoseconds(1);

  query->setAnswer(false);

  toVisit.push(query->getSource());
#ifdef ENABLE_TLS
  query->getSource()->setDFSId(timestamp);
#else // ENABLE_TLS
  query->getSource()->setDFSId(threadId, timestamp);
#endif // ENABLE_TLS

  while (!toVisit.empty()) {
    if (query->getCancel()) {
      query->setError(true);
      return;
    }

    curr = toVisit.top();
    toVisit.pop();

    for (auto it = curr->successors.begin(), end = curr->successors.end(); it != end; ++it) {
      if (*it == query->getSource()) {
        query->setAnswer(true);
        return;
      }

#ifdef ENABLE_TLS
      if ((*it)->getDFSId() != timestamp) {
        (*it)->setDFSId(timestamp);
#else // ENABLE_TLS
      if ((*it)->getDFSId(threadId) != timestamp) {
        (*it)->setDFSId(threadId, timestamp);
#endif // ENABLE_TLS
        toVisit.push(*it);
      }
    }
  }
}


SearchMethod
Graph::getMethod() {
  return searchMethod.load(memory_order_acquire);
}


void
Graph::setMethod(SearchMethod method) {
  searchMethod.store(method, memory_order_release);
}


// Internal partitioning query functions

void
#ifdef ENABLE_TLS
Graph::processPartitionQuery(Query * query) {
#else // ENABLE_TLS
Graph::processPartitionQuery(int threadId, Query * query) {
#endif // ENABLE_TLS
  // TODO
}


// Internal coarsening functions

static void
processVertex(Vertex * curr, set<Vertex *> &readySet, set<Vertex *> &processedSet) {
  set<Vertex *>::iterator it = readySet.find(curr);

  if (it == readySet.end()) {
    cerr << "ERROR: Tried to process a vertex that is not in the readySet." << endl;
    exit(EXIT_FAILURE);
  }

  readySet.erase(it);
  processedSet.insert(curr);
}


bool
Graph::addToReadySet(Vertex * curr, set<Vertex *> &readySet, set<Vertex *> &processedSet) {
  if (readySet.count(curr))
    return true;

  for (auto it = curr->predecessors.begin(), end = curr->predecessors.end(); it != end; ++it) {
    if (!processedSet.count(*it))
      return false;
  }

  readySet.insert(curr);
  return true;
}


Graph *
Graph::coarsenGreedy(int factor, map<int, int> &vertexMap) {
  int progressCount = 0;
  string barTitle = "Greedy coarsening ";
  Vertex * curr = NULL;
  Vertex * newVertex = NULL;
  Graph * coarseGraph = new Graph(false);
  map<Vertex *, int> localMapping;
  set<Vertex *> readySet;
  set<Vertex *> vertexGroup;
  set<Vertex *> processedSet;
  queue<Vertex *> workQueue;
  queue<Vertex *> localWorkQueue;
  fstream mappingStream;

#ifdef ENABLE_COARSEN_STATISTICS
  double totalWorkListSize = 0;
  double workListCount = 0;
  double totalSuccessorCount = 0;
  double illegalSuccessorCount = 0;
#endif // ENABLE_COARSE_STATISTICS

  configureProgressBar(&barTitle, getVertexCount());

  for (auto it = sources.begin(), end = sources.end(); it != end; ++it) {
    addToReadySet(*it, readySet, processedSet);
    workQueue.push(*it);
  }

  while (!workQueue.empty()) {
    curr = workQueue.front();
    workQueue.pop();

    if (processedSet.count(curr))
      continue;

#ifdef ENABLE_COARSEN_STATISTICS
    workListCount += 1;
#endif // ENABLE_COARSEN_STATISTICS

    vertexGroup.clear();
    localWorkQueue.push(curr);

    while (!localWorkQueue.empty()) {
      curr = localWorkQueue.front();
      localWorkQueue.pop();

#ifdef ENABLE_COARSEN_STATISTICS
      totalWorkListSize += 1;
      totalSuccessorCount += curr->getSuccessorCount();
#endif // ENABLE_COARSEN_STATISTICS

      processVertex(curr, readySet, processedSet);
      vertexGroup.insert(curr);

      resultProgressBar(++progressCount);

      for (auto it = curr->successors_begin(), end = curr->successors_end(); it != end; ++it) {
        if (addToReadySet(*it, readySet, processedSet))
          localWorkQueue.push(*it);
#ifdef ENABLE_COARSEN_STATISTICS
        else
          illegalSuccessorCount += 1;
#endif // ENABLE_COARSEN_STATISTICS
      }

      if (vertexGroup.size() >= (unsigned int) factor)
        break;
    }

    while (!localWorkQueue.empty()) {
      curr = localWorkQueue.front();
      localWorkQueue.pop();
      workQueue.push(curr);
    }

    newVertex = coarseGraph->addVertex();
    for (auto it = vertexGroup.begin(), end = vertexGroup.end(); it != end; ++it) {
      vertexMap[(*it)->id] = newVertex->id;
      newVertex->weight += (*it)->weight;

      for (auto predIt = (*it)->predecessors_begin(), predEnd = (*it)->predecessors_end();
          predIt != predEnd; ++predIt) {
        set<Vertex *> predecessorSet;

        if (!vertexGroup.count(*predIt)) {
          Vertex * localPred = coarseGraph->vertices[vertexMap[(*predIt)->id]];

          if (predecessorSet.find(localPred) != predecessorSet.end())
            continue;
          else
            predecessorSet.insert(localPred);

          auto localIt = localMapping.find(localPred);
          if (localIt == localMapping.end())
            localMapping[localPred] = (*it)->getPredecessorWeight(*predIt);
          else
            localIt->second += (*it)->getPredecessorWeight(*predIt);
        }
      }

      for (auto mapIt = localMapping.begin(), mapEnd = localMapping.end();
          mapIt != mapEnd; ++mapIt)
        coarseGraph->addEdgeUnsafe(mapIt->first, newVertex, mapIt->second);

      localMapping.clear();
    }
  }

  cout << endl;

#ifdef ENABLE_COARSEN_STATISTICS
  cout << "Average coarsen ratio: ";
  cout << ((double) getVertexCount()) / ((double) coarseGraph->getVertexCount()) << endl;
  cout << "Worklist count: " << workListCount << endl;
  cout << "Average worklist size: " << totalWorkListSize / workListCount << endl;
  cout << "Average successor count: " << totalSuccessorCount / totalWorkListSize << endl;
  cout << "Average illegal successor count: " << illegalSuccessorCount / totalWorkListSize << endl;
#endif // ENABLE_COARSEN_STATISTICS

  return coarseGraph;
}


static bool
reduxOrder(const Vertex * lhs, const Vertex * rhs, int f) {
  static int factor = 1;
  int lhsWeight = 0;
  int rhsWeight = 0;

  if (f != 0) {
    factor = f - 1;
    return true;
  }

  for (int i = 0, e = lhs->getSuccessorCount(); i < e; i++)
    lhsWeight += lhs->getSuccessorWeight(i);

  for (int i = 0, e = rhs->getSuccessorCount(); i < e; i++)
    rhsWeight += rhs->getSuccessorWeight(i);

  if (lhsWeight >= rhsWeight) {
    if (rhsWeight > factor)
      return false;
    else if (lhsWeight > factor)
      return true;
    else
      return false;
  } else {
    return !reduxOrder(rhs, lhs, 0);
  }
}


class edgeReduxOrder {
public:
  bool operator() (const Vertex * lhs, const Vertex * rhs) const {
    return reduxOrder(lhs, rhs, 0);
  }
};


class edgeReduxOrderReverse {
public:
  bool operator() (const Vertex * lhs, const Vertex * rhs) const {
    return !reduxOrder(lhs, rhs, 0);
  }
};


Graph *
Graph::coarsenEdgeRedux(int factor, map<int, int> &vertexMap) {
  int progressCount = 0;
  string barTitle = "Edge Redux coarsening ";
  Vertex * curr = NULL, * search = NULL, * source = NULL;
  Vertex * newVertex = NULL;
  Graph * coarseGraph = new Graph(false);
  ReachabilityQuery * query = new ReachabilityQuery(NULL, NULL, DFS);
  map<int, set<int> > retroMapping;
  set<Vertex *> vertexGroup;
  set<Vertex *> processedSet;
  priority_queue<Vertex *, vector<Vertex *>, edgeReduxOrder> workQueue;
  fstream mappingStream;

#ifdef ENABLE_COARSEN_STATISTICS
  double totalWorkListSize = 0;
  double workListCount = 0;
  double totalSuccessorCount = 0;
  double illegalSuccessorCount = 0;
  double convexSearchCount = 0;
  double convexSearchAbortsIllegal = 0;
  double convexSearchAbortsMaxSize = 0;
#endif // ENABLE_COARSEN_STATISTICS

  reduxOrder(NULL, NULL, factor);

  configureProgressBar(&barTitle, getVertexCount());

  // Push nodes in reverse ID order. This allows for nodes of same order to be
  // listed in topological order instead of a detrimental reverse topological
  // order.
  for (auto it = vertices.rbegin(), end = vertices.rend(); it != end; ++it) {
    if (*it)
      workQueue.push(*it);
  }

  while (!workQueue.empty()) {
    priority_queue<Vertex *, vector<Vertex *>, edgeReduxOrderReverse> localWorkQueue;

    curr = workQueue.top();
    workQueue.pop();

    resultProgressBar(++progressCount);

    if (processedSet.count(curr))
      continue;

    vertexGroup.clear();

    localWorkQueue.push(curr);
    source = curr;

#ifdef ENABLE_COARSEN_STATISTICS
    workListCount += 1;
#endif // ENABLE_COARSE_STATISTICS

    while (!localWorkQueue.empty()) {
      curr = localWorkQueue.top();
      localWorkQueue.pop();

#ifdef ENABLE_COARSEN_STATISTICS
      totalWorkListSize += 1;
      totalSuccessorCount += curr->getSuccessorCount();
#endif // ENABLE_COARSE_STATISTICS

      if (vertexGroup.count(curr))
        continue;

      if (curr != source) {
        bool backTrace = false;
        set<Vertex *> shadowVertexGroup;
        set<Vertex *> newlyProcessedSet;
        queue<Vertex *> searchQueue;

        shadowVertexGroup = vertexGroup;

        for (auto predIt = curr->predecessors_begin(), predEnd = curr->predecessors_end();
            predIt != predEnd; ++predIt) {
          if (!vertexGroup.count(*predIt))
            searchQueue.push(*predIt);
        }

        while (!searchQueue.empty()) {
          search = searchQueue.front();
          searchQueue.pop();

          if (vertexGroup.count(search))
            continue;

#ifdef ENABLE_COARSEN_STATISTICS
          convexSearchCount += 1;
#endif // ENABLE_COARSE_STATISTICS

          query->source = source;
          query->target = search;
          query->answer = false;
#ifdef ENABLE_STATISTICS
          query->path.clear();
          query->searchedNodes = 0;
#endif // ENABLE_STATISTICS

#ifdef ENABLE_TLS
          areConnectedDFS(query);
#else // ENABLE_TLS
          areConnectedDFS(query, 0);
#endif // ENABLE_TLS

          if (query->getAnswer()) {
#ifdef ENABLE_STATISTICS
            for (auto pathIt = query->path.begin(), pathEnd = query->path.end();
                pathIt != pathEnd; ++pathIt) {
              if (processedSet.count(*pathIt) && !vertexGroup.count(*pathIt)) {
#ifdef ENABLE_COARSEN_STATISTICS
                convexSearchAbortsIllegal += 1;
#endif // ENABLE_COARSE_STATISTICS
                backTrace = true;
                break;
              }

              vertexGroup.insert(*pathIt);
              processedSet.insert(*pathIt);
              newlyProcessedSet.insert(*pathIt);
              for (auto predIt = (*pathIt)->predecessors_begin(), predEnd = (*pathIt)->predecessors_end();
                  predIt != predEnd; ++predIt) {
                if (!vertexGroup.count(*predIt))
                  searchQueue.push(*predIt);
              }

              if (vertexGroup.size() > 2 * (unsigned int) factor) {
#ifdef ENABLE_COARSEN_STATISTICS
                convexSearchAbortsMaxSize += 1;
#endif // ENABLE_COARSE_STATISTICS
                backTrace = true;
                break;
              }
            }
#else // ENABLE_STATISTICS
            backTrace = true;
#endif // ENABLE_STATISTICS
          }

          if (backTrace) {
            vertexGroup = shadowVertexGroup;

            for (auto it = newlyProcessedSet.begin(), end = newlyProcessedSet.end(); it != end; ++it)
              processedSet.erase(*it);

            break;
          }
        }

        if (backTrace)
          continue;
      }

      vertexGroup.insert(curr);
      processedSet.insert(curr);

      if (vertexGroup.size() >= (unsigned int) factor)
        break;

      for (auto succIt = curr->successors_begin(), succEnd = curr->successors_end();
          succIt != succEnd; ++succIt) {
        if (!processedSet.count(*succIt))
          localWorkQueue.push(*succIt);
#ifdef ENABLE_COARSEN_STATISTICS
        else
          illegalSuccessorCount += 1;
#endif // ENABLE_COARSE_STATISTICS
      }
    }

    newVertex = coarseGraph->addVertex();
    for (auto it = vertexGroup.begin(), end = vertexGroup.end(); it != end; ++it) {
      vertexMap[(*it)->id] = newVertex->id;
      retroMapping[newVertex->id].insert((*it)->id);
      newVertex->weight += (*it)->weight;
    }
  }

  for (auto it = retroMapping.begin(), end = retroMapping.end(); it != end; ++it) {
    map<Vertex *, int> localMapping;
    Vertex * source = coarseGraph->getVertexFromId(it->first);

    for (auto setIt = it->second.begin(), setEnd = it->second.end();
        setIt != setEnd; ++setIt) {
      Vertex * origin = getVertexFromId(*setIt);
      set<Vertex *> predecessorSet;

      for (auto predIt = origin->predecessors_begin(), predEnd = origin->predecessors_end();
          predIt != predEnd; ++predIt) {
        if (!it->second.count((*predIt)->id)) {
          Vertex * localPred = coarseGraph->getVertexFromId(vertexMap[(*predIt)->id]);

          if (predecessorSet.count(localPred))
            continue;
          else
            predecessorSet.insert(localPred);

          auto localIt = localMapping.find(localPred);
          if (localIt == localMapping.end())
            localMapping[localPred] = origin->getPredecessorWeight(*predIt);
          else
            localIt->second += origin->getPredecessorWeight(*predIt);
        }
      }
    }

    for (auto mapIt = localMapping.begin(), mapEnd = localMapping.end(); mapIt != mapEnd; ++mapIt)
      coarseGraph->addEdgeUnsafe(mapIt->first, source, mapIt->second);
  }

  cout << endl;

#ifdef ENABLE_COARSEN_STATISTICS
  double includedNodes = 0L;
  for (auto mapIt = retroMapping.begin(), mapEnd = retroMapping.end(); mapIt != mapEnd; ++mapIt)
    includedNodes += (double) mapIt->second.size();

  includedNodes /= (double) retroMapping.size();

  cout << "Average coarsen ratio: " << includedNodes << endl;
  cout << "Worklist count: " << workListCount << endl;
  cout << "Average worklist size: " << totalWorkListSize / workListCount << endl;
  cout << "Average successor count: " << totalSuccessorCount / totalWorkListSize << endl;
  cout << "Average illegal successor count: " << illegalSuccessorCount / totalWorkListSize << endl;
  cout << "Convex search count: " << convexSearchCount << endl;
  cout << "Illegal path aborts: " << convexSearchAbortsIllegal << endl;
  cout << "MaxSize path aborts: " << convexSearchAbortsMaxSize << endl;
#endif // ENABLE_COARSE_STATISTICS

  delete query;

  return coarseGraph;
}


void
Graph::gatherEdgeInformation(set<Vertex *> &vertexGroup, Graph * coarseGraph,
    map<int, int> &IDMap, map<Vertex *, int> &weightMap, bool backward) {
  set<Vertex *> neighbourSet;

  for (auto setIt = vertexGroup.begin(), setEnd = vertexGroup.end(); setIt != setEnd; ++setIt) {
    vector<Vertex *> &neighbours = backward ? (*setIt)->predecessors : (*setIt)->successors;
    neighbourSet.clear();

    for (auto it = neighbours.begin(), end = neighbours.end(); it != end; ++it) {
      if (vertexGroup.count(*it))
        continue;

      auto localIt = IDMap.find((*it)->id);
      if (localIt != IDMap.end()) {
        Vertex * target = coarseGraph->getVertexFromId(localIt->second);

        if (neighbourSet.count(target))
          continue;
        else
          neighbourSet.insert(target);

        int weight =
          backward ? (*setIt)->getPredecessorWeight(*it) : (*setIt)->getSuccessorWeight(*it);
        auto weightIt = weightMap.find(target);
        if (weightIt == weightMap.end())
          weightMap[target] = weight;
        else
          weightIt->second += weight;
      }
    }
  }
}

Graph *
Graph::coarsenApproxIteration(int factor, int secondaryFactor, map<int, int> &vertexMap) {
  int iterations = 0;
  string barTitle = "Iterative coarsening - level ";
  map<int, set<int> > retroMapping;
  Graph * sourceGraph = NULL;
  Graph * coarseGraph = this;
  ReachabilityQuery * query = new ReachabilityQuery(NULL, NULL, Approx);
  fstream mappingStream;

#ifdef ENABLE_COARSEN_STATISTICS
  double totalWorkListSize = 0;
  double workListCount = 0;
  double convexSearchCount = 0;
  double convexSearchAborts= 0;
#endif // ENABLE_COARSEN_STATISTICS

  // Compute the number of necessary iterations
  while (factor > 1) {
    iterations++;
    factor /= secondaryFactor;
  }

  // Pre-fill the vertex mappings
  for (unsigned int i = 0, e = getVertexCount(); i < e; i++)
    retroMapping[i].insert(i);

  for (int i = 1; i <= iterations; i++) {
    int progressCount = 0;
    uint64_t iterationId =
        chrono::high_resolution_clock::now().time_since_epoch() / chrono::nanoseconds(1);
    map<int, Vertex *> sourceMap;
    map<int, set<int> > newRetroMapping;
    vector<Vertex *> workQueue;

    // Forward graphs one level
    sourceGraph = coarseGraph;
    coarseGraph = new Graph(false);
    vertexMap.clear();

    // Index the source graph
    sourceGraph->indexGraph();

    // Prepare the progress bar
    string localTitle = barTitle + to_string(i) + "/" + to_string(iterations) + " ";
    configureProgressBar(&localTitle, sourceGraph->getVertexCount());

    // Randomize the workQueue
    for (auto it = sourceGraph->vertices.begin(), end = sourceGraph->vertices.end();
        it != end; ++it) {
      if (*it)
        workQueue.push_back(*it);
    }
    random_shuffle(workQueue.begin(), workQueue.end());

    for (auto it = workQueue.begin(), end = workQueue.end(); it != end; ++it) {
      Vertex * newVertex = NULL;
      set<Vertex *> vertexGroup;
      set<int> predecessorSources;
      set<int> successorSources;
      map<Vertex *, int> localMapping;
      queue<Vertex *> localWorkQueue;

      resultProgressBar(++progressCount);

#ifdef ENABLE_TLS
      if ((*it)->getDFSId() == iterationId)
        continue;
      else
        (*it)->setDFSId(iterationId);
#else // ENABLE_TLS
      if ((*it)->getDFSId(0) == iterationId)
        continue;
      else
        (*it)->setDFSId(0, iterationId);
#endif // ENABLE_TLS

      vertexGroup.insert(*it);

      for (auto succIt = (*it)->successors_begin(), succEnd = (*it)->successors_end();
          succIt != succEnd; ++succIt) {
        if (!vertexMap.count((*succIt)->id))
          localWorkQueue.push(*succIt);
      }

#ifdef ENABLE_COARSEN_STATISTICS
      workListCount += 1;
#endif // ENABLE_COARSEN_STATISTICS

      while (!localWorkQueue.empty()) {
        bool skipVertex = false;
        Vertex * next = localWorkQueue.front();
        localWorkQueue.pop();

#ifdef ENABLE_TLS
        if (next->getDFSId() == iterationId)
#else // ENABLE_TLS
        if (next->getDFSId(0) == iterationId)
#endif // ENABLE_TLS
          continue;

#ifdef ENABLE_COARSEN_STATISTICS
        totalWorkListSize += 1;
#endif // ENABLE_COARSEN_STATISTICS

        for (auto predIt = next->predecessors_begin(), predEnd = next->predecessors_end();
            predIt != predEnd; ++predIt) {
          if (vertexGroup.count(*predIt))
            continue;

#ifdef ENABLE_COARSEN_STATISTICS
          convexSearchCount += 1;
#endif // ENABLE_COARSEN_STATISTICS

          query->source = *it;
          query->target = *predIt;
#ifdef ENABLE_STATISTICS
          query->searchedNodes = 0;
          query->path.clear();
#endif // ENABLE_STATISTICS

#ifdef ENABLE_TLS
          sourceGraph->areConnectedApprox(query);
#else // ENABLE_TLS
          sourceGraph->areConnectedApprox(query, 0);
#endif // ENABLE_TLS

          if ((skipVertex = query->getAnswer())) {
#ifdef ENABLE_COARSEN_STATISTICS
            convexSearchAborts += 1;
#endif // ENABLE_COARSEN_STATISTICS
            break;
          }
        }

        if (skipVertex)
          continue;

        vertexGroup.insert(next);

#ifdef ENABLE_TLS
        next->setDFSId(iterationId);
#else // ENABLE_TLS
        next->setDFSId(0, iterationId);
#endif // ENABLE_TLS

        for (auto succIt = next->successors_begin(), succEnd = next->successors_end();
            succIt != succEnd; ++succIt) {
          if (!vertexMap.count((*succIt)->id))
            localWorkQueue.push(*succIt);
        }

        if (vertexGroup.size() >= (unsigned int) secondaryFactor)
          break;
      }

      newVertex = coarseGraph->addVertex();
      sourceMap[newVertex->id] = *it;
      newRetroMapping.insert(pair<int, set<int> >(newVertex->id, set<int>()));
      for (auto setIt = vertexGroup.begin(), setEnd = vertexGroup.end(); setIt != setEnd; ++setIt) {
        set<int> &retroSet = retroMapping[(*setIt)->id];
        vertexMap[(*setIt)->id] = newVertex->id;
        newRetroMapping[newVertex->id].insert(retroSet.begin(), retroSet.end());
      }

      // Create predecessors edges where necessary
      localMapping.clear();
      gatherEdgeInformation(vertexGroup, coarseGraph, vertexMap, localMapping, true);

      for (auto mapIt = localMapping.begin(), mapEnd = localMapping.end();
          mapIt != mapEnd; ++mapIt)
        coarseGraph->addEdgeUnsafe(mapIt->first, newVertex, mapIt->second);

      // Create successors edges where necessary
      localMapping.clear();
      gatherEdgeInformation(vertexGroup, coarseGraph, vertexMap, localMapping, false);

      for (auto mapIt = localMapping.begin(), mapEnd = localMapping.end();
          mapIt != mapEnd; ++mapIt)
        coarseGraph->addEdgeUnsafe(newVertex, mapIt->first, mapIt->second);
    }

    retroMapping = newRetroMapping;
    coarseGraph->indexMethod = sourceGraph->indexMethod;

    // When necessary clean-up the intermediary graph
    if (sourceGraph != this)
      delete sourceGraph;
  }

  // Construct the forward mapping of vertex IDs
  vertexMap.clear();
  for (auto mapIt = retroMapping.begin(), mapEnd = retroMapping.end(); mapIt != mapEnd; ++mapIt) {
    for (auto setIt = mapIt->second.begin(), setEnd = mapIt->second.end(); setIt != setEnd; ++setIt)
      vertexMap[*setIt] = mapIt->first;
  }

  cout << endl;

#ifdef ENABLE_COARSEN_STATISTICS
  double includedNodes = 0L;
  for (auto mapIt = retroMapping.begin(), mapEnd = retroMapping.end(); mapIt != mapEnd; ++mapIt)
    includedNodes += (double) mapIt->second.size();

  includedNodes /= (double) retroMapping.size();

  cout << "Average coarsen ratio: " << includedNodes << endl;
  cout << "Worklist count: " << workListCount << endl;
  cout << "Average worklist size: " << totalWorkListSize / workListCount << endl;
  cout << "Convex search count: " << convexSearchCount << endl;
  cout << "Convex search aborts: " << convexSearchAborts << endl;
#endif // ENABLE_COARSE_STATISTICS

  return coarseGraph;
}


#ifdef ENABLE_STATISTICS
// Internal statistics maintenance

void
Graph::registerQueryStatistics(ReachabilityQuery * query) {
  double coefficient = 1.0;
  double overhead = 1.0;
  unique_lock<mutex> statisticsLock(statisticsMutex);

  queryCount++;

  if (query->getAnswer()) {
    positiveQueryCount++;

    coefficient = 1.0 / ((double) positiveQueryCount);
    if (query->getMethod() == DFS)
      overhead = ((double) query->searchedNodes) / ((double) query->path.size());

    positiveQueryOverhead *= (1.0 - coefficient);
    positiveQueryOverhead += coefficient * overhead;
  } else {
    negativeQueryCount++;

    if (query->searchedNodes > 0) {
      coefficient = 1.0 / ((double) negativeQueryCount);
      overhead = ((double) query->searchedNodes) / ((double) getVertexCount());

      negativeQueryOverhead *= (1.0 - coefficient);
      negativeQueryOverhead += coefficient * overhead;
    } else {
      shortNegativeQueryCount++;
    }
  }
}

#endif // ENABLE_STATISTICS


// Local help and work functions

unsigned long int getThreadID(void) {
  hash<thread::id> h;
  return h(this_thread::get_id());
}

mutex printMutex;

void
queryWorker(Graph * graph, int threadId) {
  Query * query = NULL;

  /*
  unique_lock<mutex> printLock(printMutex);
  cerr << "Initializing PAPI for thread " << threadId << " with system-id ";
  cerr << getThreadID() << endl;
  printLock.unlock();

  int errCode = PAPI_register_thread();
  if (errCode != PAPI_OK) {
    cerr << "ERROR: Could not initialize worker thread for PAPI measurements in thread ";
    cerr << threadId << ". Error was : ";
    switch (errCode) {
      case PAPI_ENOMEM:
        cerr << "PAPI_ENOMEM" << endl;
        break;
      case PAPI_ESYS:
        cerr << "PAPI_ESYS" << endl;
        break;
      case PAPI_ESBSTR:
        cerr << "PAPI_ESBSTR" << endl;
        break;
      default:
        cerr << "Unknown (" << errCode << ")" << endl;
    }
    exit(EXIT_FAILURE);
  }
  */

  while (true) {
    while (!graph->jobQueue.try_pop(query)) {}
      // Do nothing

    if (graph->threadShutdown)
      return;

    graph->startQuery();

    switch (query->type) {
#ifdef ENABLE_TLS
      case Reachability:
        graph->processReachabilityQuery(query);
        break;

      case Partition:
        graph->processPartitionQuery(query);
        break;
#else // ENABLE_TLS
      case Reachability:
        graph->processReachabilityQuery(threadId, query);
        break;

      case Partition:
        graph->processPartitionQuery(threadId, query);
        break;
#endif // ENABLE_TLS
    }

    graph->stopQuery();
  }

  //PAPI_unregister_thread();
}
