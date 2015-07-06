#include <map>
#include <set>
#include <queue>
#include <stack>
#include <chrono>
#include <climits>
#include <cstring>
#include <fstream>
#include <iostream>

#include <unistd.h>

#include "graph.hpp"

using namespace std;


// Local global variables

#ifdef ENABLE_BENCHMARKS
static PAPI_dmem_info_t memoryInfo;
static long long baseMemoryUsage;
#endif // ENABLE_BENCHMARKS


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
Vertex::addPredecessor(Vertex * pred) {
  if (pred == this)
    return false;

  for (auto it = predecessors.begin(), end = predecessors.end(); it != end; ++it) {
    if (*it == pred)
      return false;
  }
  
  predecessors.push_back(pred);
  predecessorCount++;
  return true;
}


bool
Vertex::addSuccessor(Vertex * succ) {
  if (succ == this)
    return false;

  for (auto it = successors.begin(), end = successors.end(); it != end; ++it) {
    if (*it == succ)
      return false;
  }

  successors.push_back(succ);
  successorCount++;
  return true;
}


bool
Vertex::removePredecessor(Vertex * pred) {
  for (auto it = predecessors.begin(), end = predecessors.end(); it != end; ++it) {
    if (*it == pred) {
      predecessors.erase(it);
      predecessorCount--;
      return true;
    }
  }

  return false;
}


bool
Vertex::removeSuccessor(Vertex * succ) {
  for (auto it = successors.begin(), end = successors.end(); it != end; ++it) {
    if (*it == succ) {
      successors.erase(it);
      successorCount--;
      return true;
    }
  }

  return false;
}


int
Vertex::getNumberOfPredecessors() {
  return predecessorCount;
}


int
Vertex::getNumberOfSuccessors() {
  return successorCount;
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

void
Vertex::setDFSId(int idx, uint64_t id) {
  DFSId[idx].store(id, memory_order_release);
}


uint64_t
Vertex::getDFSId(int idx) {
  return DFSId[idx].load(memory_order_acquire);
}


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
  if (((inVisits + outVisits) == 0) && (method & 0x02))
    (method & 0x01) ? successors.clear() : predecessors.clear();

  // When necessary push the visiting parent on the reordered vector
  if (method & 0x02)
    (method & 0x01) ? successors.push_back(pred) : predecessors.push_back(pred);

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
      while (graph->vertices.size() <= (unsigned) maxId)
        graph->addVertexUnsafe();

      if (noDoubleEdges)
        graph->addEdgeUnsafe(graph->vertices[source], graph->vertices[target]);
      else
        graph->addEdge(graph->vertices[source], graph->vertices[target]);
    } else {
      sscanf(dump, "%d", &source);
      while (graph->vertices.size() <= (unsigned) source)
        graph->addVertex();
    }
  }

  input.close();
  graph->indexed = false;
  graph->condensed = false;

  // Make sure the graph is a DAG
  graph->condenseGraph();

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

  // Get the first line out of the way
  input.getline(dump, 127);

  // Get the number of vertices / lines to read
  input >> lineNumber;
  configureProgressBar(&barTitle, lineNumber);

  // Create all vertices
  for (int i = 0; i < lineNumber; i++)
    graph->addVertexUnsafe();

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
        graph->addEdgeUnsafe(graph->vertices[source], graph->vertices[target]);
      else
        graph->addEdge(graph->vertices[source], graph->vertices[target]);
    }

    resultProgressBar(i + 1);
  }
  cout << "\n";

  input.close();
  graph->indexed = false;
  graph->condensed = false;

  // Make sure the graph is a DAG
  graph->condenseGraph();

#ifdef ENABLE_BENCHMARKS
  PAPI_get_dmem_info(&memoryInfo);
  graph->graphMemoryUsage = memoryInfo.resident - baseMemoryUsage;
#endif // ENABLE_BENCHMARKS

  return graph;
}


// Modificators

Vertex *
Graph::addVertex(void) {
  startGlobalOperation();

  int id = vertices.size();;
  Vertex * newVertex = new Vertex(id);

  vertices.push_back(newVertex);
  indexed = false;

  stopGlobalOperation();

  return newVertex;
}


void
Graph::removeVertex(Vertex * v) {
  vertices.back()->id = v->id;
  vertices[v->id] = vertices.back();
  vertices.pop_back();
  delete v;
}


void
Graph::mergeVertices(Vertex * s, Vertex * t) {
  startGlobalOperation();

  for (auto it = s->predecessors.begin(), end = s->predecessors.end(); it != end; ++it)
    t->addPredecessor(*it);

  for (auto it = s->successors.begin(), end = s->successors.end(); it != end; ++it)
    t->addSuccessor(*it);

  removeVertex(s);
  indexed = false;

  stopGlobalOperation();
}


bool
Graph::addEdge(Vertex * source, Vertex * target) {
  startGlobalOperation();

  if (!source->addSuccessor(target))
    return false;

  target->addPredecessor(source);
  edgeCount++;
  indexed = false;
  condensed = false;

  stopGlobalOperation();
  return true;
}


bool
Graph::removeEdge(Vertex * source, Vertex * target) {
  startGlobalOperation();

  if (!source->removeSuccessor(target))
    return false;

  source->removePredecessor(source);
  edgeCount--;

  return true;

  startGlobalOperation();
}


void
Graph::setIndexMethod(IndexMethod newMethod) {
  indexMethod = newMethod;
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
Graph::getEdgeCount() {
  return edgeCount;
}


unsigned int
Graph::getVertexCount() {
  return vertices.size();
}


Vertex *
Graph::getVertexFromId(int id) {
  return vertices[id];
}


IndexMethod
Graph::getIndexMethod() {
  return indexMethod;
}


// Queries

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

  for (int i = 0; i < MAX_THREADS; i++)
    jobQueue.push(NULL);
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
  os << "Number of vertices: " << vertices.size() << "\n";
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
  int indexMemoryUsage = (8 * vertices.size()) / 1024;
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

  if (indexed)
    return;

  // Make sure we are indexing a DAG
  condenseGraph();

  // Make sure we discovered the sources
  discoverExtremities();

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

      for (auto it = curr->predecessors.begin(), end = curr->predecessors.end(); it != end; ++it)
        (*it)->successors.push_back(curr);

      curr->successors.clear();
    }

    while (!predecessorQueue.empty()) {
      curr = predecessorQueue.back();
      predecessorQueue.pop_back();

      for (auto it = curr->successors.begin(), end = curr->successors.end(); it != end; ++it)
        (*it)->predecessors.push_back(curr);

      curr->predecessors.clear();
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
}


// Maintenance

void
Graph::discoverExtremities() {
  sources.clear();
  sinks.clear();

  for (auto it = vertices.begin(), end = vertices.end(); it != end; ++it) {
    if ((*it)->predecessors.empty())
      sources.push_back(*it);

    if ((*it)->successors.empty())
      sinks.push_back(*it);
  }
}


void
Graph::condenseGraph() {
  if (condensed)
    return;

  // Prepare for the upcoming DFSs
  globalTimestamp =
    chrono::high_resolution_clock::now().time_since_epoch() / chrono::nanoseconds(1);

  // Iterate over the vertices and condense
  for (auto it = vertices.begin(); it != vertices.end(); ++it) {
    if ((*it)->getDFSId(0) != globalTimestamp)
      condenseFromSource(*it);
  }

  condensed = true;
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
}


void
Graph::stopGlobalOperation() {
  unique_lock<mutex> queryWaitLock(queryWaitMutex);

  noQueries = false;

  queryWaitCondition.notify_all();
}


// Maintenance

void
Graph::condenseFromSource(Vertex * source) {
  set<Vertex *> DFSSet;
  vector<Vertex *> DFSPath;
  stack<Vertex *> toVisit;
  map<Vertex *, Vertex *> mergeMap;
  Vertex * curr = NULL;

  toVisit.push(source);
  source->setDFSId(0, globalTimestamp);

  while (!toVisit.empty()) {
    curr = toVisit.top();
    if (!DFSPath.empty() && (curr == DFSPath.back())) {
      DFSPath.pop_back();
      DFSSet.erase(curr);
      toVisit.pop();
      continue;
    }

    // Insert the current vertex in to the DFS path
    DFSPath.push_back(curr);
    DFSSet.insert(curr);

    // Iterate over the successors of the current vertex
    for (auto it = curr->successors.begin(), end = curr->successors.end(); it != end; ++it) {
      if (DFSSet.count(*it)) {
        // Mark the cycle for merging into the first vertex of the cycle
        auto it2 = DFSPath.rbegin();
        while (*it2 != *it) {
          if (mergeMap.count(*it2) == 0)
            mergeMap[*it2] = *it;

          it2++;
        }
      } else {
        // See if we need to visit this successor
        if ((*it)->getDFSId(0) != globalTimestamp) {
          (*it)->setDFSId(0, globalTimestamp);
          toVisit.push(*it);
        }
      }
    }
  }

  // Perform the necessary merging
  for (auto it = mergeMap.begin(), end = mergeMap.end(); it != end; ++it) {
    // Find the final merge target
    Vertex * mergeTarget = it->second;
    while (mergeMap.count(mergeTarget) == 1)
      mergeTarget = mergeMap[mergeTarget];

    // Merge the two vertices
    mergeVertices(it->first, mergeTarget);
  }
}


Vertex *
Graph::addVertexUnsafe() {
  int id = vertices.size();;
  Vertex * newVertex = new Vertex(id);

  vertices.push_back(newVertex);
  indexed = false;

  return newVertex;
}


bool
Graph::addEdgeUnsafe(Vertex * source, Vertex * target) {
  startGlobalOperation();
  source->successors.push_back(target);
  source->successorCount++;
  target->predecessors.push_back(source);
  target->predecessorCount++;
  edgeCount++;
  stopGlobalOperation();
  return true;
}


// Internal reachability query functions

void
Graph::processReachabilityQuery(int threadId, Query * query) {
  vector<Vertex *> searchStack;
  unique_lock<mutex> internalResultLock(internalResultMutex, defer_lock);

  switch (query->query.reachability->getMethod()) {
    case DFS:
      searchStack.clear();
      areConnectedDFS(query->query.reachability, threadId, searchStack);
      break;

    case BBFS:
      areConnectedBBFS(query->query.reachability, threadId);
      break;

    case NoLabels:
      areConnectedNoLabels(query->query.reachability, threadId);
      break;

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
Graph::areConnectedDFS(ReachabilityQuery * query, int threadId, vector<Vertex *> &searchStack) {
  int event = PAPI_TOT_CYC;
  long long counterValue;
  uint64_t timestamp;
  Vertex * curr;

#ifdef ENABLE_BENCHMARKS
  unique_lock<mutex> benchmarkLock(benchmarkMutex, defer_lock);
  PAPI_start_counters(&event, 1);
#else // ENABLE_BENCHMARKS
  if (query->getInternal() == 1)
    PAPI_start_counters(&event, 1);
#endif // ENABLE_BENCHMARKS

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
    if (query->getCancel())
      goto cancel;

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

      if ((*it)->getDFSId(threadId) != timestamp) {
        (*it)->setDFSId(threadId, timestamp);

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
Graph::areConnectedBBFS(ReachabilityQuery * query, int threadId) {
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
      (query->getTarget()->retroReverseOrderLabel > query->getSource()->retroReverseOrderLabel)) {
    goto end;
#endif // ENABLE_RETRO_LABELS

  // Do a BBFS on the subgraph specified by both orders to get the final answer
  forwardId = chrono::high_resolution_clock::now().time_since_epoch() / chrono::nanoseconds(1);
  backwardId = forwardId + 1;

  query->getSource()->setDFSId(threadId, forwardId);
  query->getTarget()->setDFSId(threadId, backwardId);

  searchQueueForward.push(query->getSource());
  searchQueueBackward.push(query->getTarget());

  while (!searchQueueForward.empty() || !searchQueueBackward.empty()) {
    if (query->getCancel())
      goto cancel;

    if (!searchQueueForward.empty()) {
      curr = searchQueueForward.front();
      searchQueueForward.pop();

#ifdef ENABLE_STATISTICS
      query->searchedNodes++;
#endif // ENABLE_STATISTICS

      for (auto it = curr->successors.begin(), end = curr->successors.end(); it != end; ++it) {
        if ((indexMethod & 0x04) && ((*it)->orderLabel > query->getTarget()->orderLabel))
          break;

        if ((*it)->getDFSId(threadId) == backwardId) {
          query->setAnswer(true);
          goto end;
        }

        if ((*it)->getDFSId(threadId) != forwardId) {
          (*it)->setDFSId(threadId, forwardId);

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

        if ((*it)->getDFSId(threadId) == forwardId) {
          query->setAnswer(true);
          goto end;
        }

        if ((*it)->getDFSId(threadId) != backwardId) {
          (*it)->setDFSId(threadId, backwardId);

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
Graph::areConnectedNoLabels(ReachabilityQuery * query, int threadId) {
  uint64_t timestamp;
  stack<Vertex *> toVisit;
  Vertex * curr;

  timestamp = chrono::high_resolution_clock::now().time_since_epoch() / chrono::nanoseconds(1);

  toVisit.push(query->getSource());
  query->getSource()->setDFSId(threadId, timestamp);

  while (!toVisit.empty()) {
    if (query->getCancel())
      return;

    curr = toVisit.top();
    toVisit.pop();

    for (auto it = curr->successors.begin(), end = curr->successors.end(); it != end; ++it) {
      if (*it == query->getSource()) {
        query->setAnswer(true);
        return;
      }

      if ((*it)->getDFSId(threadId) != timestamp) {
        (*it)->setDFSId(threadId, timestamp);
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
Graph::processPartitionQuery(int threadId, Query * query) {
  // TODO
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
      overhead = ((double) query->searchedNodes) / ((double) vertices.size());

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


void
queryWorker(Graph * graph, int threadId) {
  Query * query = NULL;

  if (PAPI_thread_init(getThreadID) != PAPI_OK) {
    cerr << "ERROR: Could not initialize worker thread for PAPI measurements." << endl;
    exit(EXIT_FAILURE);
  }

  while (true) {
    while (!graph->jobQueue.try_pop(query)) {}
      // Do nothing

    if (graph->threadShutdown)
      return;

    graph->startQuery();

    switch (query->type) {
      case Reachability:
        graph->processReachabilityQuery(threadId, query);
        break;

      case Partition:
        graph->processPartitionQuery(threadId, query);
        break;
    }

    graph->stopQuery();
  }
}
