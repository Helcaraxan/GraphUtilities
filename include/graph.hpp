/*
 * graph.hpp
 */

/* Forward declarations */
class Edge;
class Vertex;
class Graph;


#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <stack>
#include <vector>

//#define ENABLE_RETRO_LABELS
#define ENABLE_STATISTICS
#define ENABLE_PAPI_BENCHMARKS

using namespace std;

class Vertex;
class Graph;


/* Class declarations */
class Vertex {
friend class Graph;

// Local types
public:
  typedef vector<Vertex *>::iterator iterator;

// Data members
public:
	int id; /* A unique Vertex ID */

private:
	int orderLabel; /* Ordering in a reverse post-order traversal of the Graph */
	int reverseOrderLabel; /* Ordering in a reverse post-order traversal of the inverse Graph */
#ifdef ENABLE_RETRO_LABELS
  int retroOrderLabel;
  int retroReverseOrderLabel;
#endif // ENABLE_RETRO_LABELS

	int inVisits; /* Variable used in reverse post-order */
	int outVisits; /* Variable used in reverse post-order */
	int DFSId; /* Variable to distinguish between seperate DFSs */

  Vertex * firstVisit;

  int predecessorCount;
  int successorCount;
	vector<Vertex *> predecessors;
	vector<Vertex *> successors;

// Function members
public:
	// Constructors & Destructor
	Vertex(int i);
	~Vertex(void);

	// Modificators
	bool addPredecessor(Vertex * pred);
	bool addSuccessor(Vertex * succ);
  bool removePredecessor(Vertex * pred);
  bool removeSuccessor(Vertex * succ);

  // Access
  int getNumberOfPredecessors(void);
  int getNumberOfSuccessors(void);

  // Iterators
  iterator predecessors_begin(void);
  iterator predecessors_end(void);

  iterator successors_begin(void);
  iterator successors_end(void);

private:
	// Indexing
	void visit(Vertex * pred, int method); /* Function used in post-order traversal */
	Vertex * createPostOrder(stack<Vertex *> * postOrder, int method); /* Idem */
};


class Graph {
// Local types
public:
  /* The IndexMethod values obey the following rules:
   * 0x01 (first bit) Indicates if predecessors are rescheduled on traversal
   * 0x02 (second bit) Indicates if successors are rescheduled on traversal
   * 0x04 (third bit) Indicates if successors are ordered with their labels
   */
  typedef enum {
    ShortIndexing = 0x00,
    SuccessorOrder = 0x02,
    Standard = 0x03,
    DFSShortCut = 0x05,
    UndefinedMethod = 0x10
  } IndexMethod;

// Data members
protected:
	bool indexed; /* Indicates if the ordering of the graph has been done */
  bool condensed;
  IndexMethod indexMethod;
  unsigned int edgeCount;

	vector<Vertex *> vertices;
	vector<Vertex *> sources;
	vector<Vertex *> sinks;

  // Global ID to distinguish between seperate DFSs
	int DFSId;

private:
  bool statisticsEnabled;
  bool papiBenchmarksEnabled;

#ifdef ENABLE_STATISTICS
  // Counters
  uintmax_t queryCount;
  uintmax_t positiveQueryCount;
  uintmax_t negativeQueryCount;
  uintmax_t shortNegativeQueryCount;

  // DFS overhead statistics
  double positiveQueryOverhead;
  double negativeQueryOverhead;
#endif // ENABLE_STATISTICS

#ifdef ENABLE_PAPI_BENCHMARKS
  // Counters
  long long queryNumber;
  long long cyclesSpentIndexing;
  long long cyclesSpentQuerying;
  long long graphMemoryUsage;

  // PAPI administration
  int benchmarkEvents[1];
#endif // ENABLE_PAPI_BENCHMARKS

// Function members
public:
	// Constructors & Destructor
	Graph(void);
	~Graph(void);

  // Parser functions
  static Graph * createFromDotFile(const char * fileName, bool noDoubleEdges = false);
  static Graph * createFromGraFile(const char * fileName, bool noDoubleEdges = false);

	// Modificators
	Vertex * addVertex(void);
  void removeVertex(Vertex * v);
  void mergeVertices(Vertex * s, Vertex * t);
	bool addEdge(Vertex * s, Vertex * t);
  bool removeEdge(Vertex * s, Vertex * t);
  void setIndexMethod(IndexMethod newMethod);

  // Access
  unsigned int getEdgeCount(void);
  unsigned int getVertexCount(void);
  Vertex * getVertexFromId(int id);
  IndexMethod getIndexMethod(void);

	// Queries
  /* The query function returning NULL or the path connecting the vertices */
	vector<Vertex *> * areConnected(Vertex * u, Vertex * v, vector<Vertex *> * path);

  /* Find if two vertices are connected through a multi-hop path */
  bool indirectPathExists(Vertex * u, Vertex * v);

  /* Do a query without considering labels : standard DFS */
  bool areConnectedDFS(Vertex * u, Vertex * v);
  bool areConnectedBBFS(Vertex * u, Vertex * V);

  // Statistics
  bool statisticsAreEnabled(void);
  uintmax_t getQueryCount(void);
  uintmax_t getPositiveQueryCount(void);
  uintmax_t getNegativeQueryCount(void);
  uintmax_t getShortNegativeQueryCount(void);
  double getPositiveQueryOverhead(void);
  double getNegativeQueryOverhead(void);
  void printStatistics(ostream &os);

  // Benchmarking
  bool benchmarksAreEnabled(void);
  long long getQueryNumber(void);
  long long getCyclesSpentIndexing(void);
  long long getCyclesSpentQuerying(void);
  void printBenchmarks(ostream &os);

protected:
	// Indexing
	void labelVertices(bool retro, bool reverse);
	void indexGraph();

  // Maintenance
  void discoverExtremities(void);
  void condenseGraph(void);

private:
  // Maintenance
  void condenseFromSource(Vertex * source);

#ifdef ENABLE_STATISTICS
  // Maintenance
  bool addEdgeUnsafe(Vertex * source, Vertex * target);

  // Internal statistic maintenance
  void registerQueryStatistics(vector<Vertex *> * path, uintmax_t searchedNodes);
#endif // ENABLE_STATISTICS
};

#endif // GRAPH_HPP
