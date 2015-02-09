#include <map>
#include <set>
#include <queue>
#include <stack>
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

thread_local int threadId;

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
    cerr << "Error while opening input Dot file.\n";
    exit(EXIT_FAILURE);
  }

  input.getline(dump, 127);
  if (!strstr(dump, "digraph")) {
    cerr << "Error - the supplied file is not a graph in Dot format.\n";
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
    cerr << "Error while opening input Dot file.\n";
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
        cerr << "Error while reading the graph definition.\n";
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
  // Verify that the graph has been indexed
  if (!indexed) {
    indexGraph();
    pthread_mutex_lock(&methodMutex);
    preferredMethod = Undefined;
    DFSwin = 0;
    BBFSwin = 0;
    pthread_mutex_unlock(&methodMutex);
  }

  // If a method is specified use it
  if (query->method != Undefined) {
    pthread_mutex_lock(&jobMutex);
    jobQueue.push_back(query);
    sem_post(&jobSemaphore);
    pthread_mutex_unlock(&jobMutex);
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

    pthread_mutex_lock(&jobMutex);
    jobQueue.push_back(&DFSQuery);
    jobQueue.push_back(&BBFSQuery);
    sem_post(&jobSemaphore);
    sem_post(&jobSemaphore);
    pthread_mutex_unlock(&jobMutex);

    // Wait for the first result
    while ((result != &DFSQuery) && (result != &BBFSQuery)) {
      sem_wait(&internalResultSemaphore);
      pthread_mutex_lock(&internalResultMutex);

      auto it = internalResultQueue.begin();
      for (auto end = internalResultQueue.end(); it != end; ++it) {
        if (((*it) == &DFSQuery) || ((*it) == &BBFSQuery)) {
          result = *it;
          internalResultQueue.erase(it);
          break;
        }
      }

      if (it == internalResultQueue.end())
        sem_post(&internalResultSemaphore);

      pthread_mutex_unlock(&internalResultMutex);
    }

    // Cancel the jobs to stop the worker threads
    DFSQuery.cancel = true;
    BBFSQuery.cancel = true;

    // Only register winning search method on positive queries
    if (result == &DFSQuery) {
      pthread_mutex_lock(&methodMutex);
      DFSwin++;
      pthread_mutex_unlock(&methodMutex);

      query->answer = DFSQuery.answer;
      query->method = DFS;
    } else if (result == &BBFSQuery) {
      pthread_mutex_lock(&methodMutex);
      BBFSwin++;
      pthread_mutex_unlock(&methodMutex);

      query->answer = BBFSQuery.answer;
      query->method = BBFS;
    }

    pthread_mutex_lock(&resultMutex);
    resultQueue.push_back(query);
    sem_post(&resultSemaphore);
    pthread_mutex_unlock(&resultMutex);

    if (DFSwin >= 10) {
      preferredMethod = DFS;
      cerr << "\nMulti-thread calibration chose DFS as search-method.\n";
    } else if (BBFSwin >= 10) {
      preferredMethod = BBFS;
      cerr << "\nMulti-thread calibration chose BBFS as search-method.\n";
    }

    // Clean-up
    if (result == &DFSQuery)
      result = &BBFSQuery;
    else if (result == &BBFSQuery)
      result = &DFSQuery;

    if (!result->error) {
      while (result) {
        sem_wait(&internalResultSemaphore);
        pthread_mutex_lock(&internalResultMutex);

        auto it = internalResultQueue.begin();
        for (auto end = internalResultQueue.end(); it != end; ++it) {
          if ((*it) == result) {
            internalResultQueue.erase(it);
            result = NULL;
            break;
          }
        }

        if (it == internalResultQueue.end())
          sem_post(&internalResultSemaphore);

        pthread_mutex_unlock(&internalResultMutex);
      }
    }
  } else {
    query->method = preferredMethod;

    pthread_mutex_lock(&jobMutex);
    jobQueue.push_back(query);
    sem_post(&jobSemaphore);
    pthread_mutex_unlock(&jobMutex);
  }

  return;
}


Graph::Query *
Graph::pullResult() {
  Query * result = NULL;

  sem_wait(&resultSemaphore);
  pthread_mutex_lock(&resultMutex);
  result = resultQueue.front();
  resultQueue.pop_front();
  pthread_mutex_unlock(&resultMutex);

  return result;
}


void
Graph::endOfQueries() {
  threadShutdown = true;

  for (int i = 0; i < MAX_THREADS; i++)
    sem_post(&jobSemaphore);
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
  os << "WARNING: To enable statistics uncomment the #define in the header file 'graph.h'.\n\n";
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
  os << "WARNING: To enable benchmarks uncomment the #define in the header file 'graph.h'.\n\n";
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
  pthread_mutex_lock(&indexMutex);

  if (indexed)
    return;

  // Make sure we are indexing a DAG
  condenseGraph();

  // Make sure we discovered the sources
  discoverExtremities();

  if (indexMethod == Graph::UndefinedMethod) {
    cerr << "Unknown indexing method. Aborting.\n";
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

  pthread_mutex_unlock(&indexMutex);
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
  pthread_mutex_lock(&queryWaitMutex);

  if (noQueries)
    pthread_cond_wait(&queryWaitCondition, &queryWaitMutex);

  activeThreads++;
  pthread_mutex_unlock(&queryWaitMutex);
}


void
Graph::stopQuery() {
  pthread_mutex_lock(&queryWaitMutex);
  activeThreads--;
  pthread_mutex_unlock(&queryWaitMutex);

  if (noQueries && (activeThreads == 0)) {
    pthread_mutex_lock(&globalWaitMutex);
    pthread_cond_signal(&globalWaitCondition);
    pthread_mutex_unlock(&globalWaitMutex);
  }
}


void
Graph::startGlobalOperation() {
  pthread_mutex_lock(&queryWaitMutex);
  noQueries = true;
  pthread_mutex_unlock(&queryWaitMutex);

  pthread_mutex_lock(&globalWaitMutex);
  if (activeThreads > 0)
    pthread_cond_wait(&globalWaitCondition, &globalWaitMutex);

  pthread_mutex_unlock(&globalWaitMutex);
}


void
Graph::stopGlobalOperation() {
  pthread_mutex_lock(&queryWaitMutex);
  noQueries = false;

  pthread_cond_broadcast(&queryWaitCondition);
  pthread_mutex_unlock(&queryWaitMutex);
}


// Internal query functions

/* Use the previously done indexation to answer to the query */
void *
Graph::areConnectedDFS(void * arg) {
  Query * query = (Query *) arg;
  Vertex * curr;
  stack<Vertex *> searchStack;

#ifdef ENABLE_STATISTICS
  vector<Vertex *> path;
  uintmax_t searchedNodes = 0;
#endif // ENABLE_STATISTICS

#ifdef ENABLE_BENCHMARKS
  int event = PAPI_TOT_CYC;
  long long counterValue;
  PAPI_start_counters(&event, 1);
#endif // ENABLE_BENCHMARKS

  // Are U and V the same vertex?
  if (query->source == query->target) {
#ifdef ENABLE_STATISTICS
    searchedNodes++;
    path.push_back(query->source);
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
  searchStack.push(query->source);
  DFSId[threadId]++;

  while (!searchStack.empty()) {
    if (query->cancel)
      goto cancel;

    curr = searchStack.top();

#ifdef ENABLE_STATISTICS
    if (!path.empty() && (curr == path.back())) {
      path.pop_back();
      searchStack.pop();
      continue;
    }

    searchedNodes++;
    path.push_back(curr);
#else // ENABLE_STATISTICS
    searchStack.pop();
#endif // ENABLE_STATISTICS

    for (auto it = curr->successors.begin(), end = curr->successors.end(); it != end; ++it) {
      if ((indexMethod & 0x04) && ((*it)->orderLabel > query->target->orderLabel))
        break;

      if (*it == query->target) {
#ifdef ENABLE_STATISTICS
        searchedNodes++;
        path.push_back(query->target);
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
  pthread_mutex_lock(&benchmarkMutex);
  cyclesSpentQuerying += counterValue;
  queryNumber++;
  pthread_mutex_unlock(&benchmarkMutex);
#endif // ENABLE_BENCHMARKS

#ifdef ENABLE_STATISTICS
  registerQueryStatistics(query->answer, path.size(), searchedNodes);
#endif // ENABLE_STATISTICS

  checkIdOverflow();

  return arg;

cancel:
#ifdef ENABLE_BENCHMARKS
  PAPI_stop_counters(&counterValue, 1);
#endif // ENABLE_BENCHMARKS
  
  checkIdOverflow();
  query->error = true;

  return arg;
}


void *
Graph::areConnectedBBFS(void * arg) {
  int forwardId, backwardId;
  Vertex * curr;
  Query * query = (Query *) arg;
  queue<Vertex *> searchQueueForward;
  queue<Vertex *> searchQueueBackward;

#ifdef ENABLE_BENCHMARKS
  int event = PAPI_TOT_CYC;
  long long counterValue;
  PAPI_start_counters(&event, 1);
#endif // ENABLE_BENCHMARKS

  // Are U and V the same vertex?
  if (query->source == query->target) {
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

  // Do a DFS on the subgraph specified by both orders to get the final answer
  searchQueueForward.push(query->source);
  searchQueueBackward.push(query->target);
  forwardId = ++DFSId[threadId];
  backwardId = ++DFSId[threadId];
  query->source->DFSId[threadId] = forwardId;
  query->target->DFSId[threadId] = backwardId;

  while (!searchQueueForward.empty() || !searchQueueBackward.empty()) {
    if (query->cancel)
      goto cancel;

    if (!searchQueueForward.empty()) {
      curr = searchQueueForward.front();
      searchQueueForward.pop();

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
  pthread_mutex_lock(&benchmarkMutex);
  cyclesSpentQuerying += counterValue;
  queryNumber++;
  pthread_mutex_unlock(&benchmarkMutex);
#endif // ENABLE_BENCHMARKS

#ifdef ENABLE_STATISTICS
  registerQueryStatistics(query->answer, 0, 0);
#endif // ENABLE_STATISTICS

  checkIdOverflow();

  return arg;

cancel:
#ifdef ENABLE_BENCHMARKS
  PAPI_stop_counters(&counterValue, 1);
#endif // ENABLE_BENCHMARKS
  
  checkIdOverflow();
  query->error = true;

  return arg;
}


void *
Graph::areConnectedNoLabels(void * arg) {
  stack<Vertex *> toVisit;
  Query * query = (Query *) arg;
  Vertex * curr;

  DFSId[threadId]++;

  toVisit.push(query->source);
  query->source->DFSId[threadId] = DFSId[threadId];

  while (!toVisit.empty()) {
    if (query->cancel)
      return arg;

    curr = toVisit.top();
    toVisit.pop();

    for (auto it = curr->successors.begin(), end = curr->successors.end(); it != end; ++it) {
      if (*it == query->source) {
        query->answer = true;
        return arg;
      }
        ;

      if ((*it)->DFSId[threadId] != DFSId[threadId]) {
        (*it)->DFSId[threadId] = DFSId[threadId];
        toVisit.push(*it);
      }
    }
  }

  checkIdOverflow();

  return arg;
}


void
Graph::checkIdOverflow() {
  if (threadId >= MAX_THREADS) {
    cerr << "Thread ID higher than maximum thread count." << endl;
    return;
  }

  if (DFSId[threadId] > (UCHAR_MAX - 2)) {
    DFSId[threadId] = 0;
    for (auto it = vertices.begin(), end = vertices.end(); it != end; ++it)
      (*it)->DFSId[threadId] = 0;
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
Graph::registerQueryStatistics(bool result, unsigned int pathSize, uintmax_t searchedNodes) {
  double coefficient, overhead;

  pthread_mutex_lock(&statisticsMutex);

  queryCount++;

  if (result) {
    positiveQueryCount++;

    coefficient = 1.0 / ((double) positiveQueryCount);
    overhead = ((double) searchedNodes) / ((double) pathSize);

    positiveQueryOverhead *= (1.0 - coefficient);
    positiveQueryOverhead += coefficient * overhead;
  } else {
    negativeQueryCount++;

    if (searchedNodes) {
      coefficient = 1.0 / ((double) negativeQueryCount);
      overhead = ((double) searchedNodes) / ((double) vertices.size());

      negativeQueryOverhead *= (1.0 - coefficient);
      negativeQueryOverhead += coefficient * overhead;
    } else {
      shortNegativeQueryCount++;
    }
  }

  pthread_mutex_unlock(&statisticsMutex);
}

#endif // ENABLE_STATISTICS


// Local work function

void *
queryWorker(void * args) {
  struct QueryWorkerArgs * queryArgs = (struct QueryWorkerArgs *) args;
  Graph::Query * query = NULL;
  Graph * graph = queryArgs->graphArg;
  threadId = queryArgs->threadIdArg;

  while (true) {
    sem_wait(&graph->jobSemaphore);

    if (graph->threadShutdown)
      return NULL;

    pthread_mutex_lock(&graph->jobMutex);
    query = graph->jobQueue.front();
    graph->jobQueue.pop_front();
    pthread_mutex_unlock(&graph->jobMutex);

    graph->startQuery();

    switch (query->method) {
      case Graph::DFS:
        graph->areConnectedDFS((void *) query);
        break;

      case Graph::BBFS:
        graph->areConnectedBBFS((void *) query);
        break;

      case Graph::NoLabels:
        graph->areConnectedNoLabels((void *) query);
        break;

      case Graph::Undefined:
        query->error = true;
        break;
    }

    graph->stopQuery();

    if (query->internal) {
      pthread_mutex_lock(&graph->internalResultMutex);
      graph->internalResultQueue.push_back(query);
      sem_post(&graph->internalResultSemaphore);
      pthread_mutex_unlock(&graph->internalResultMutex);
    } else {
      pthread_mutex_lock(&graph->resultMutex);
      graph->resultQueue.push_back(query);
      sem_post(&graph->resultSemaphore);
      pthread_mutex_unlock(&graph->resultMutex);
    }

  }

  return NULL;
}
