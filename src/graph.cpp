#include <set>
#include <list>
#include <cstring>
#include <fstream>
#include <iostream>

#include "graph.hpp"

#ifdef ENABLE_PAPI_BENCHMARKS
#include <papi.h>
#endif // ENABLE_PAPI_BENCHMARKS

using namespace std;


// Constructors & destructors

Vertex::Vertex(int i) :
	id(i),
	orderLabel(-1),
	reverseOrderLabel(-1),
	inVisits(0),
	outVisits(0),
	queryID(0)
{}


Vertex::~Vertex(void)
{}


Graph::Graph(void) :
	indexed(false),
  edgeCount(0),
	queryID(0)
{
#ifdef ENABLE_STATISTICS
  statisticsEnabled = true;
  queryCount = 0;
  positiveQueryCount = 0;
  negativeQueryCount = 0;
  shortNegativeQueryCount = 0;
  positiveQueryOverhead = 0.0L;
  negativeQueryOverhead = 0.0L;
#else // ENABLE_STATISTICS
  statisticsEnabled = false;
#endif // ENABLE_STATISTICS

#ifdef ENABLE_PAPI_BENCHMARKS
  queryNumber = 0;
  benchmarkEvents[0] = PAPI_TOT_CYC;

  if (PAPI_num_counters() < 1)
    papiBenchmarksEnabled = false;
  else
    papiBenchmarksEnabled = true;
#else // ENABLE_PAPI_BENCHMARKS
  papiBenchmarksEnabled = false;
#endif // ENABLE_PAPI_BENCHMARKS
}


Graph::~Graph(void) {
	for (auto it = vertices.begin(); it != vertices.end(); ++it)
		delete *it;

	vertices.clear();
}


// Modificators

bool
Vertex::addPredecessor(Vertex * pred) {
	for (auto it = predecessors.begin(); it != predecessors.end(); ++it) {
		if (*it == pred)
			return false;
	}
	
	predecessors.push_back(pred);
	return true;
}


bool
Vertex::addSuccessor(Vertex * succ) {
	for (auto it = successors.begin(); it != successors.end(); ++it) {
		if (*it == succ)
			return false;
	}

	successors.push_back(succ);
	return true;
}


bool
Vertex::removePredecessor(Vertex * pred) {
  for (auto it = predecessors.begin(); it != predecessors.end(); ++it) {
    if (*it == pred) {
      predecessors.erase(it);
      return true;
    }
  }

  return false;
}


bool
Vertex::removeSuccessor(Vertex * succ) {
  for (auto it = successors.begin(); it != successors.end(); ++it) {
    if (*it == succ) {
      successors.erase(it);
      return true;
    }
  }

  return false;
}


int
Vertex::getNumberOfPredecessors() {
  return predecessors.size();
}


int
Vertex::getNumberOfSuccessors() {
  return successors.size();
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
 * to the order in which the parent nodes access this Vertex instance during
 * graph traversal.
 */
void
Vertex::visit(Vertex * pred, bool reverse) {
	vector<Vertex *> * orderedVertices;

	// Reinitialize in case of reverse post-order traversal
	if (reverse && ((inVisits + outVisits) == (int) (successors.size() + predecessors.size()))) {
		inVisits = 0;
		outVisits = 0;
	}

  // Determine which Vertex vector should be reordered
	orderedVertices = reverse ? &successors : &predecessors;

  // If this is the first node to be visited there is no reordering
	if (!pred)
    return;
    
  // If the ordering is OK do nothing except registering the current visit
  if ((*orderedVertices)[inVisits] == pred) {
    inVisits++;
		return;
  }

  // When necessary reorder the Vertex vector as to put the parent provoking the
  // n-th visit of this Vertex on the n-th spot
  for (auto it = orderedVertices->begin(); it != orderedVertices->end(); ++it) {
    if (*it == pred) {
			*it = (*orderedVertices)[inVisits];
			(*orderedVertices)[inVisits++] = pred;
			break;
		}
	}
}


/* Compute the next Vertex to be visited in a post-order traversal of the graph
 * from the current Vertex instance. This allow for an iterative traversal and
 * not a recursive one in order to prevent stack-overflows on large graphs
 */
Vertex *
Vertex::createPostOrder(stack<Vertex *> * postOrder, bool reverse) {
	vector<Vertex *> * upVertices;
	vector<Vertex *> * downVertices;

  /* Determine what is up and what is down depending on the reverse boolean */
	upVertices = reverse ? &successors : &predecessors;
	downVertices = reverse ? &predecessors : &successors;

  /* Visit all the child nodes */
	while (outVisits < (int) downVertices->size()) {
		(*downVertices)[outVisits]->visit(this, reverse);

		if ((*downVertices)[outVisits]->inVisits == 1)
			return (*downVertices)[outVisits++];
		else
			outVisits++;
	}

  /* Push the node on the post-order stack only when all children have been
   * visited */
	postOrder->push(this);

  /* Return the predecessor from which the first visit to this node was made as
   * next Vertex */
  if (upVertices->size() == 0)
		return NULL;
	else
		return (*upVertices)[0];
}


// Parser

/* This function parses a subset of Dot files. Only edges and nodes can be
 * defined without any attributes
 */
Graph *
Graph::createFromDotFile(const char * fileName) {
  char dump[128];
  int source, target, maxId;
  fstream input(fileName, fstream::in);
  Graph * graph = new Graph();
  
  if (input == NULL) {
    fprintf(stderr, "Error while opening input Dot file.\n");
    exit(EXIT_FAILURE);
  }

  input.getline(dump, 127);
  if (!strstr(dump, "digraph")) {
    fprintf(stderr, "Error - the supplied file is not a graph in Dot format.\n");
    exit(EXIT_FAILURE);
  }

  while (input.good()) {
    input.getline(dump, 127);

    if (strchr(dump, '}')) {
      input.close();
      return graph;
    }

    if (strchr(dump, '>')) {
      sscanf(dump, "%d -> %d", &source, &target);
      maxId = source < target ? target : source;
      while (graph->vertices.size() <= (unsigned) maxId)
        graph->addVertex();

      graph->addEdge(graph->vertices[source], graph->vertices[target]);
    } else {
      sscanf(dump, "%d", &source);
      while (graph->vertices.size() <= (unsigned) source)
        graph->addVertex();
    }
  }

  return graph;
}


// Modificators

Vertex *
Graph::addVertex(void) {
	int id = vertices.size();;
	Vertex * newVertex = new Vertex(id);

	vertices.push_back(newVertex);

	return newVertex;
}


bool
Graph::addEdge(Vertex * source, Vertex * target) {
	if (!source->addSuccessor(target))
		return false;

	target->addPredecessor(source);
  edgeCount++;
	return true;
}


bool
Graph::removeEdge(Vertex * source, Vertex * target) {
  if (!source->removeSuccessor(target))
    return false;

  source->removePredecessor(source);
  edgeCount--;
  return true;
}


// Access

uintmax_t
Graph::getEdgeCount() {
  return edgeCount;
}


uintmax_t
Graph::getVertexCount() {
  return vertices.size();
}


Vertex *
Graph::getVertexFromId(int id) {
  return vertices[id];
}


// Queries

/* Use the previously done indexation to answer to the query */
vector<Vertex *> * Graph::areConnected(Vertex * u, Vertex * v, vector<Vertex *> * path) {
	Vertex * curr;
	stack<Vertex *> searchStack;

#ifdef ENABLE_STATISTICS
  uintmax_t searchedNodes = 0;
#endif // ENABLE_STATISTICS

  // Verify that the graph has been indexed
  if (!indexed)
    indexGraph();

  // Clear the specified path container
  path->clear();

#ifdef ENABLE_PAPI_BENCHMARKS
  PAPI_start_counters(benchmarkEvents, 1);
#endif // ENABLE_PAPI_BENCHMARKS

	// Are U and V the same vertex?
	if (u == v) {
    path->push_back(u);
    goto positiveEnd;
  }

	// Can V be a descendant of U in the standard graph or U a descendant of V in
  // the reverse graph?
	if ((u->orderLabel > v->orderLabel) ||
      (v->reverseOrderLabel > u->reverseOrderLabel))
    goto negativeEnd;

	// Do a DFS on the subgraph specified by both orders to get the final answer
	searchStack.push(u);
	queryID++;

	while (!searchStack.empty()) {
		curr = searchStack.top();
    if (path->size() && (curr == path->back())) {
      path->pop_back();
      searchStack.pop();
      continue;
    }

#ifdef ENABLE_STATISTICS
    searchedNodes++;
#endif // ENABLE_STATISTICS

    path->push_back(curr);
		for (auto it = curr->successors.begin(); it !=curr->successors.end(); ++it) {
      if (*it == v) {
        path->push_back(v);
        goto positiveEnd;
      }

			if ((*it)->queryID != queryID) {
				(*it)->queryID = queryID;

				if (((*it)->orderLabel < v->orderLabel) && ((*it)->reverseOrderLabel > v->reverseOrderLabel))
					searchStack.push(*it);
			}
		}
	}

#ifdef ENABLE_PAPI_BENCHMARKS
  long long counterValue;
  PAPI_stop_counters(&counterValue, 1);
  cyclesSpentQuerying += counterValue;
  queryNumber++;
#endif // ENABLE_PAPI_BENCHMARKS

negativeEnd:
#ifdef ENABLE_STATISTICS
  registerQueryStatistics(NULL, searchedNodes);
#endif // ENABLE_STATISTICS
	return NULL;

positiveEnd:
#ifdef ENABLE_STATISTICS
  registerQueryStatistics(path, searchedNodes + 1);
#endif // ENABLE_STATISTICS
  return path;
}


bool
Graph::indirectPathExists(Vertex * u, Vertex * v) {
  Vertex * currVertex;
  list<Vertex *> toVisit;
  set<Vertex *> scheduled;

  // Fill up the initial toVisit list
  for (auto it = u->successors.begin(); it != u->successors.end(); ++it) {
    if (*it != v) {
      toVisit.push_back(*it);
      scheduled.insert(*it);
    }
  }

  while (!toVisit.empty()) {
    currVertex = toVisit.front();
    toVisit.pop_front();

    for (auto it = currVertex->successors.begin(); it != currVertex->successors.end(); ++it) {
      if (*it == v)
        return true;

      if (scheduled.count(*it) == 0) {
        toVisit.push_front(*it);
        scheduled.insert(*it);
      }
    }
  }

  return false;
}


// Benchmark statistics

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
  return papiBenchmarksEnabled;
}


long long
Graph::getQueryNumber(void) {
#ifdef ENABLE_PAPI_BENCHMARKS
  return queryNumber;
#else // ENABLE_PAPI_BENCHMARKS
  return 0;
#endif // ENABLE_PAPI_BENCHMARKS
}


long long
Graph::getCyclesSpentIndexing(void) {
#ifdef ENABLE_PAPI_BENCHMARKS
  return cyclesSpentIndexing;
#else // ENABLE_PAPI_BENCHMARKS
  return 0;
#endif // ENABLE_PAPI_BENCHMARKS
}


long long
Graph::getCyclesSpentQuerying(void) {
#ifdef ENABLE_PAPI_BENCHMARKS
  return cyclesSpentQuerying;
#else // ENABLE_PAPI_BENCHMARKS
  return 0;
#endif // ENABLE_PAPI_BENCHMARKS
}


void
Graph::printBenchmarks(ostream &os) {
#ifdef ENABLE_PAPI_BENCHMARKS
  os << "\n---\nBenchmarking\n";
  os << "Cycles spent indexing: " << cyclesSpentIndexing << "\n";
  os << "Cycles spent querying: " << cyclesSpentQuerying << "\n";
  os << "-> average of " << cyclesSpentQuerying / queryNumber << " cycles per query.\n---\n\n";
#else // ENABLE_PAPI_BENCHMARKS
  os << "WARNING: Benchmarking has not been enabled at compile time.\n";
  os << "WARNING: To enable benchmarks uncomment the #define in the header file 'graph.h'.\n\n";
#endif // ENABLE_PAPI_BENCHMARKS
}


// Indexing

void
Graph::labelVertices(bool reverse) {
	int currLabel = 0;
	stack<Vertex *> postOrder;
	Vertex * nextVertex;
	vector<Vertex *> * startVertices;;

	startVertices = reverse ? &sinks : &sources;

	// Loop while there are unlabeled sources / sinks
	for (auto it = startVertices->begin(); it != startVertices->end(); ++it) {
		// Initialize
		nextVertex = *it;
		(*it)->visit(NULL, reverse);

		// Run the DFS
		while (nextVertex) {
			// In case of the first forward labeling create the sinks list
			if (!reverse && (nextVertex->successors.size() == 0))
				sinks.push_back(nextVertex);

			// Use an iterative method to prevent memory overflow in large graphs
			nextVertex = nextVertex->createPostOrder(&postOrder, reverse);
		}
	}

  // Label the nodes in reverse post-order
  while (!postOrder.empty()) {
    if (reverse)
      postOrder.top()->reverseOrderLabel = currLabel++;
    else
      postOrder.top()->orderLabel = currLabel++;;

    postOrder.pop();
  }
}


void
Graph::indexGraph() {
#ifdef ENABLE_PAPI_BENCHMARKS
  PAPI_start_counters(benchmarkEvents, 1);
#endif // ENABLE_PAPI_BENCHMARKS

	// Create the vector with the source nodes
	for (auto it = vertices.begin(); it != vertices.end(); ++it) {
		if ((*it)->predecessors.size() == 0)
			sources.push_back(*it);
	}

	// Perform the forward graph ordering
	labelVertices(false);

	// Perform the reverse graph ordering
	labelVertices(true);

  indexed = true;

#ifdef ENABLE_PAPI_BENCHMARKS
  PAPI_stop_counters(&cyclesSpentIndexing, 1);
#endif // ENABLE_PAPI_BENCHMARKS
}


#ifdef ENABLE_STATISTICS
void
Graph::registerQueryStatistics(vector<Vertex *> * path, uintmax_t searchedNodes) {
  double coefficient, overhead;

  queryCount++;

  if (path) {
    positiveQueryCount++;

    coefficient = 1.0 / ((double) positiveQueryCount);
    overhead = ((double) searchedNodes) / ((double) path->size());

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
}
#endif // ENABLE_STATISTICS
