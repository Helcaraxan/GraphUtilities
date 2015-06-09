#include <map>
#include <set>
#include <queue>
#include <stack>
#include <chrono>
#include <climits>
#include <cstring>
#include <fstream>
#include <iostream>

#include <papi.h>
#include <unistd.h>

#include "graph.hpp"

using namespace std;


// Local global variables

#ifdef ENABLE_BENCHMARKS
static PAPI_dmem_info_t memoryInfo;
static long long baseMemoryUsage;
#endif // ENABLE_BENCHMARKS


// Semaphore class (active wait variant)

void
semaphore::post() {
  unique_lock<mutex> lck(mtx);
  count++;
  cv.notify_one();
}


void
semaphore::wait() {
  unique_lock<mutex> lck(mtx);
  while (count == 0)
    cv.wait(lck);

  count--;
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


// Query declarations

Vertex *
Graph::Query::getSource() {
  return source;
}


Vertex *
Graph::Query::getTarget() {
  return target;
}


bool
Graph::Query::getAnswer() {
  return answer;
}


bool
Graph::Query::isError() {
  return error;
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

  cout << "Start parsing the graph definition from a Dot file.\n";
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

  cout << "Finished parsing the graph from a Dot file.\n";
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

  cout << "Start parsing the graph definition from a Gra file.\n\n";

  // Get the first line out of the way
  input.getline(dump, 127);

  // Get the number of vertices / lines to read
  input >> lineNumber;

  // Create all vertices
  for (int i = 0; i < lineNumber; i++)
    graph->addVertexUnsafe();

  // Parse the adjacency list
  for (int i = 0; i < lineNumber; i++) {
    if (((i % 100000) == 0) && (i != 0))
      cerr << "Done " << i << " nodes.\n";

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
  }

  cout << "Finished parsing the graph from a Gra file.\n\n";
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
Graph::setIndexMethod(Graph::IndexMethod newMethod) {
  indexMethod = newMethod;
}


// Access

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


Graph::IndexMethod
Graph::getIndexMethod() {
  return indexMethod;
}


// Queries

void
Graph::pushQuery(Query * query) {
  static int DFSwin = 0;
  static int BBFSwin = 0;
  static Graph::SearchMethod preferredMethod = Undefined;

  unique_lock<mutex> jobLock(jobMutex, defer_lock);
  unique_lock<mutex> methodLock(methodMutex, defer_lock);
  unique_lock<mutex> resultLock(resultMutex, defer_lock);
  unique_lock<mutex> internalResultLock(internalResultMutex, defer_lock);

  // Verify that the graph has been indexed
  if (!indexed) {
    indexGraph();
    methodLock.lock();
    preferredMethod = Undefined;
    DFSwin = 0;
    BBFSwin = 0;
    methodLock.unlock();
  }

  // If a method is specified use it
  if (query->method != Undefined) {
    jobLock.lock();
    jobQueue.push_back(query);
    jobSemaphore.post();
    jobLock.unlock();
    return;
  }

  // Else select the method automatically. If there is no
  // preferred method yet then find one.
  if (preferredMethod == Undefined) {
    Query * result = NULL;
    Query DFSQuery(query->source, query->target, DFS);
    Query BBFSQuery(query->source, query->target, BBFS);

    DFSQuery.internal = true;
    BBFSQuery.internal = true;

    jobLock.lock();
    jobQueue.push_back(&DFSQuery);
    jobQueue.push_back(&BBFSQuery);
    jobSemaphore.post();
    jobSemaphore.post();
    jobLock.unlock();

    // Wait for the first result
    while ((result != &DFSQuery) && (result != &BBFSQuery)) {
      internalResultSemaphore.wait();
      internalResultLock.lock();

      auto it = internalResultQueue.begin();
      for (auto end = internalResultQueue.end(); it != end; ++it) {
        if (((*it) == &DFSQuery) || ((*it) == &BBFSQuery)) {
          result = *it;
          internalResultQueue.erase(it);
          break;
        }
      }

      if (it == internalResultQueue.end())
        internalResultSemaphore.post();

      internalResultLock.unlock();
    }

    // Cancel the jobs to stop the worker threads
    DFSQuery.cancel = true;
    BBFSQuery.cancel = true;

    // Only register winning search method on positive queries
    if (result->method == DFS) {
      methodLock.lock();
      DFSwin++;
      methodLock.unlock();

      query->answer = DFSQuery.answer;
      query->method = DFS;
    } else if (result->method == BBFS) {
      methodLock.lock();
      BBFSwin++;
      methodLock.unlock();

      query->answer = BBFSQuery.answer;
      query->method = BBFS;
    }

#ifdef ENABLE_STATISTICS
    registerQueryStatistics(result);
#endif // ENABLE_STATISTICS

    resultLock.lock();
    resultQueue.push_back(query);
    resultSemaphore.post();
    resultLock.unlock();

    if (preferredMethod == Undefined) {
      methodLock.lock();

      if (DFSwin >= 10) {
        preferredMethod = DFS;
        cerr << "\nMulti-thread calibration chose DFS as search-method.\n";
      } else if (BBFSwin >= 10) {
        preferredMethod = BBFS;
        cerr << "\nMulti-thread calibration chose BBFS as search-method.\n";
      }

      methodLock.unlock();
    }

    // Clean-up
    if (result == &DFSQuery)
      result = &BBFSQuery;
    else if (result == &BBFSQuery)
      result = &DFSQuery;

    if (!result->error) {
      while (result) {
        internalResultSemaphore.wait();
        internalResultLock.lock();

        auto it = internalResultQueue.begin();
        for (auto end = internalResultQueue.end(); it != end; ++it) {
          if ((*it) == result) {
            internalResultQueue.erase(it);
            result = NULL;
            break;
          }
        }

        if (it == internalResultQueue.end())
          internalResultSemaphore.post();

        internalResultLock.unlock();
      }
    }
  } else {
    query->method = preferredMethod;

    jobLock.lock();
    jobQueue.push_back(query);
    jobSemaphore.post();
    jobLock.unlock();
  }

  return;
}


Graph::Query *
Graph::pullResult() {
  Query * result = NULL;
  unique_lock<mutex> resultLock(resultMutex, defer_lock);

  resultSemaphore.wait();
  resultLock.lock();
  result = resultQueue.front();
  resultQueue.pop_front();
  resultLock.unlock();

  return result;
}


void
Graph::endOfQueries() {
  threadShutdown = true;

  for (int i = 0; i < MAX_THREADS; i++)
    jobSemaphore.post();
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
  if (positiveQueryCount) {
    os.precision(3);
    os << "- Average DFS length / path-length ratio: " << positiveQueryOverhead << "\n"; 
  }

  os << "\nNegative query statistics:\n";
  os << "- Number of negative answers   : " << negativeQueryCount << "\n";
  if (negativeQueryCount) {
    os << "- Number of immediate negatives: " << shortNegativeQueryCount;
    os.precision(4);
    os << " (" << shortFraction * 100 << "%).\n";
    os.precision(3);
    os << "- Average DFS length / graph-size ratio: " << negativeQueryOverhead << "\n";
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

  if (indexMethod == Graph::UndefinedMethod) {
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

  cout << "Condensing the graph.\n";
  cout << "Initial graph has " << vertices.size();
  cout << " vertices and " << edgeCount << " edges.\n\n";

  // Prepare for the upcoming DFSs
  DFSId[0]++;

  // Iterate over the vertices and condense
  for (auto it = vertices.begin(); it != vertices.end(); ++it) {
    if ((*it)->DFSId[0] != DFSId[0])
      condenseFromSource(*it);
  }

  condensed = true;

  cout << "Condensed graph has " << vertices.size();
  cout << " vertices and " << edgeCount << " edges.\n\n";
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


// Internal query functions

/* Use the previously done indexation to answer to the query */
void
Graph::areConnectedDFS(Query * query, int threadId) {
  Vertex * curr;
  stack<Vertex *> searchStack;

#ifdef ENABLE_BENCHMARKS
  int event = PAPI_TOT_CYC;
  long long counterValue;
  PAPI_start_counters(&event, 1);
  unique_lock<mutex> benchmarkLock(benchmarkMutex, defer_lock);
#endif // ENABLE_BENCHMARKS

  // Are U and V the same vertex?
  if (query->source == query->target) {
#ifdef ENABLE_STATISTICS
    query->searchedNodes++;
    query->path.push_back(query->source);
#endif // ENABLE_STATISTICS
    query->answer = true;
    goto end;
  }

  // Can V be a descendant of U in the standard graph or U a descendant of V in
  // the reverse graph?
  if ((query->source->orderLabel > query->target->orderLabel) ||
      (query->target->reverseOrderLabel > query->source->reverseOrderLabel))
    goto end;

#ifdef ENABLE_RETRO_LABELS
  if ((query->source->retroOrderLabel > query->target->retroOrderLabel) ||
      (query->target->retroReverseOrderLabel > query->source->retroReverseOrderLabel))
    goto end;
#endif // ENABLE_RETRO_LABELS

  // Do a DFS on the subgraph specified by both orders to get the final answer
  DFSId[threadId] =
    chrono::high_resolution_clock::now().time_since_epoch() / chrono::nanoseconds(1);

  searchStack.push(query->source);

  while (!searchStack.empty()) {
    if (query->cancel)
      goto cancel;

    curr = searchStack.top();

#ifdef ENABLE_STATISTICS
    if (!query->path.empty() && (curr == query->path.back())) {
      query->path.pop_back();
      searchStack.pop();
      continue;
    }

    query->searchedNodes++;
    query->path.push_back(curr);
#else // ENABLE_STATISTICS
    searchStack.pop();
#endif // ENABLE_STATISTICS

    for (auto it = curr->successors.begin(), end = curr->successors.end(); it != end; ++it) {
      if ((indexMethod & 0x04) && ((*it)->orderLabel > query->target->orderLabel))
        break;

      if (*it == query->target) {
#ifdef ENABLE_STATISTICS
        query->searchedNodes++;
        query->path.push_back(query->target);
#endif // ENABLE_STATISTICS

        query->answer = true;
        goto end;
      }

      if ((*it)->DFSId[threadId] != DFSId[threadId]) {
        (*it)->DFSId[threadId] = DFSId[threadId];

#ifdef ENABLE_RETRO_LABELS
        if (((*it)->orderLabel < query->target->orderLabel) &&
            ((*it)->reverseOrderLabel > query->target->reverseOrderLabel) &&
            ((*it)->retroOrderLabel < query->target->retroOrderLabel) &&
            ((*it)->retroReverseOrderLabel > query->target->retroReverseOrderLabel))
#else // ENABLE_RETRO_LABELS
        if (((*it)->orderLabel < query->target->orderLabel) &&
            ((*it)->reverseOrderLabel > query->target->reverseOrderLabel))
#endif // ENABLE_RETRO_LABELS
          searchStack.push(*it);
      }
    }
  }

end:
#ifdef ENABLE_BENCHMARKS
  PAPI_stop_counters(&counterValue, 1);
  benchmarkLock.lock();
  cyclesSpentQuerying += counterValue;
  queryNumber++;
  benchmarkLock.unlock();
#endif // ENABLE_BENCHMARKS

  return;

cancel:
#ifdef ENABLE_BENCHMARKS
  PAPI_stop_counters(&counterValue, 1);
#endif // ENABLE_BENCHMARKS
  
  query->error = true;
  return;
}


void
Graph::areConnectedBBFS(Query * query, int threadId) {
  uint64_t forwardId, backwardId;
  Vertex * curr;
  queue<Vertex *> searchQueueForward;
  queue<Vertex *> searchQueueBackward;

#ifdef ENABLE_BENCHMARKS
  int event = PAPI_TOT_CYC;
  long long counterValue;
  PAPI_start_counters(&event, 1);
  unique_lock<mutex> benchmarkLock(benchmarkMutex, defer_lock);
#endif // ENABLE_BENCHMARKS

  // Are U and V the same vertex?
  if (query->source == query->target) {
#ifdef ENABLE_STATISTICS
    query->searchedNodes++;
#endif // ENABLE_STATISTICS
    query->answer = true;
    goto end;
  }

  // Can V be a descendant of U in the standard graph or U a descendant of V in
  // the reverse graph?
  if ((query->source->orderLabel > query->target->orderLabel) ||
      (query->target->reverseOrderLabel > query->source->reverseOrderLabel)) {
#ifdef ENABLE_STATISTICS
    shortNegativeQueryCount++;
#endif //ENABLE_STATISTICS
    goto end;
  }

#ifdef ENABLE_RETRO_LABELS
  if ((query->source->retroOrderLabel > query->target->retroOrderLabel) ||
      (query->target->retroReverseOrderLabel > query->source->retroReverseOrderLabel)) {
#ifdef ENABLE_STATISTICS
    shortNegativeQueryCount++;
#endif // ENABLE_STATISTICS
    goto end;
  }
#endif // ENABLE_RETRO_LABELS

  // Do a BBFS on the subgraph specified by both orders to get the final answer
  forwardId = chrono::high_resolution_clock::now().time_since_epoch() / chrono::nanoseconds(1);
  backwardId = forwardId + 1;

  query->source->DFSId[threadId] = forwardId;
  query->target->DFSId[threadId] = backwardId;

  searchQueueForward.push(query->source);
  searchQueueBackward.push(query->target);

  while (!searchQueueForward.empty() || !searchQueueBackward.empty()) {
    if (query->cancel)
      goto cancel;

    if (!searchQueueForward.empty()) {
      curr = searchQueueForward.front();
      searchQueueForward.pop();

#ifdef ENABLE_STATISTICS
      query->searchedNodes++;
#endif // ENABLE_STATISTICS

      for (auto it = curr->successors.begin(), end = curr->successors.end(); it != end; ++it) {
        if ((indexMethod & 0x04) && ((*it)->orderLabel > query->target->orderLabel))
          break;

        if ((*it)->DFSId[threadId] == backwardId) {
          query->answer = true;
          goto end;
        }

        if ((*it)->DFSId[threadId] != forwardId) {
          (*it)->DFSId[threadId] = forwardId;

#ifdef ENABLE_RETRO_LABELS
          if (((*it)->orderLabel < query->target->orderLabel) &&
              ((*it)->reverseOrderLabel > query->target->reverseOrderLabel) &&
              ((*it)->retroOrderLabel < query->target->retroOrderLabel) &&
              ((*it)->retroReverseOrderLabel > query->target->retroReverseOrderLabel))
#else // ENABLE_RETRO_LABELS
          if (((*it)->orderLabel < query->target->orderLabel) &&
              ((*it)->reverseOrderLabel > query->target->reverseOrderLabel))
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
        if ((indexMethod & 0x04) && ((*it)->orderLabel < query->source->orderLabel))
          break;

        if ((*it)->DFSId[threadId] == forwardId) {
          query->answer = true;
          goto end;
        }

        if ((*it)->DFSId[threadId] != backwardId) {
          (*it)->DFSId[threadId] = backwardId;

#ifdef ENABLE_RETRO_LABELS
          if (((*it)->orderLabel > query->source->orderLabel) &&
              ((*it)->reverseOrderLabel < query->source->reverseOrderLabel) &&
              ((*it)->retroOrderLabel > query->source->retroOrderLabel) &&
              ((*it)->retroReverseOrderLabel < query->source->retroReverseOrderLabel))
#else // ENABLE_RETRO_LABELS
          if (((*it)->orderLabel > query->source->orderLabel) &&
              ((*it)->reverseOrderLabel < query->source->reverseOrderLabel))
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
  queryNumber++;
  benchmarkLock.unlock();
#endif // ENABLE_BENCHMARKS

  return;

cancel:
#ifdef ENABLE_BENCHMARKS
  PAPI_stop_counters(&counterValue, 1);
#endif // ENABLE_BENCHMARKS
  
  query->error = true;
  return;
}


void
Graph::areConnectedNoLabels(Query * query, int threadId) {
  stack<Vertex *> toVisit;
  Vertex * curr;

  DFSId[threadId] =
    chrono::high_resolution_clock::now().time_since_epoch() / chrono::nanoseconds(1);

  toVisit.push(query->source);
  query->source->DFSId[threadId] = DFSId[threadId];

  while (!toVisit.empty()) {
    if (query->cancel)
      return;

    curr = toVisit.top();
    toVisit.pop();

    for (auto it = curr->successors.begin(), end = curr->successors.end(); it != end; ++it) {
      if (*it == query->source) {
        query->answer = true;
        return;
      }

      if ((*it)->DFSId[threadId] != DFSId[threadId]) {
        (*it)->DFSId[threadId] = DFSId[threadId];
        toVisit.push(*it);
      }
    }
  }
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
  source->DFSId[0] = DFSId[0];

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
        if ((*it)->DFSId[0] != DFSId[0]) {
          (*it)->DFSId[0] = DFSId[0];
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


#ifdef ENABLE_STATISTICS
// Internal statistics maintenance

void
Graph::registerQueryStatistics(Query * query) {
  double coefficient = 1.0;
  double overhead = 1.0;
  unique_lock<mutex> statisticsLock(statisticsMutex);

  queryCount++;

  if (query->answer) {
    positiveQueryCount++;

    coefficient = 1.0 / ((double) positiveQueryCount);
    if (query->method == DFS)
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


// Local work function

void
queryWorker(Graph * graph, int threadId) {
  Graph::Query * query = NULL;
  unique_lock<mutex> jobLock(graph->jobMutex, defer_lock);
  unique_lock<mutex> resultLock(graph->resultMutex, defer_lock);
  unique_lock<mutex> internalResultLock(graph->internalResultMutex, defer_lock);

  while (true) {
    graph->jobSemaphore.wait();

    if (graph->threadShutdown)
      return;

    jobLock.lock();
    query = graph->jobQueue.front();
    graph->jobQueue.pop_front();
    jobLock.unlock();

    graph->startQuery();

    switch (query->method) {
      case Graph::DFS:
        graph->areConnectedDFS(query, threadId);
        break;

      case Graph::BBFS:
        graph->areConnectedBBFS(query, threadId);
        break;

      case Graph::NoLabels:
        graph->areConnectedNoLabels(query, threadId);
        break;

      case Graph::Undefined:
        query->error = true;
        break;
    }

    graph->stopQuery();

    if (query->internal) {
      internalResultLock.lock();
      graph->internalResultQueue.push_back(query);
      graph->internalResultSemaphore.post();
      internalResultLock.unlock();
    } else {
#ifdef ENABLE_STATISTICS
      graph->registerQueryStatistics(query);
#endif // ENABLE_STATISTICS

      resultLock.lock();
      graph->resultQueue.push_back(query);
      graph->resultSemaphore.post();
      resultLock.unlock();
    }
  }
}
