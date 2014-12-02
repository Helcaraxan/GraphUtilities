/*
 * graph.hpp
 */

/* Forward declarations */
class Edge;
class Vertex;
class Graph;


#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <vector>

//#define ENABLE_RETRO_LABELS
//#define ENABLE_STATISTICS
//#define ENABLE_PAPI_BENCHMARKS

using namespace std;

class Vertex;
class Graph;
class DFSlauncher;
class BBFSlauncher;


// Class declarations
class Vertex {
friend class Graph;

// Local types
public:
  typedef vector<Vertex *>::iterator iterator;

// Data members
public:
	int id;

private:
	int orderLabel;
	int reverseOrderLabel;
#ifdef ENABLE_RETRO_LABELS
  int retroOrderLabel;
  int retroReverseOrderLabel;
#endif // ENABLE_RETRO_LABELS

	int inVisits;
	int outVisits;
	int DFSId;

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
	void visit(Vertex * pred, int method);
	Vertex * createPostOrder(vector<Vertex *> * postOrder, int method);
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
    LabelOrder = 0x07,
    UndefinedMethod = 0x10
  } IndexMethod;

  typedef enum {
    DFS,
    BBFS,
    NoLabels,
    Undefined
  } SearchMethod;

private:
  struct Query_t {
    Vertex * source;
    Vertex * target;
    bool timed;
    bool answer;
    long long searchTime;

    Query_t(Vertex * u, Vertex * v) :
      source(u),
      target(v),
      timed(false),
      answer(false),
      searchTime(-1)
    {}
  };

  typedef struct Query_t Query;

// Data members
protected:
	bool indexed;
  bool condensed;
  IndexMethod indexMethod;
  unsigned int edgeCount;

	vector<Vertex *> vertices;
	vector<Vertex *> sources;
	vector<Vertex *> sinks;

  // Global ID to distinguish between seperate DFSs
	int DFSId;

private:
  // Flags
  bool statisticsEnabled;
  bool papiBenchmarksEnabled;

  // Reordering vectors
  vector<Vertex *> successorQueue;
  vector<Vertex *> predecessorQueue;

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
  bool areConnected(Vertex * u, Vertex * v);
  bool areConnected(Vertex * u, Vertex * v, SearchMethod method);
  bool indirectPathExists(Vertex * u, Vertex * v);

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
	void indexGraph(void);

  // Maintenance
  void discoverExtremities(void);
  void condenseGraph(void);

private:
  // Maintenance
  void condenseFromSource(Vertex * source);

  // Internal query functions
  void * areConnectedDFS(void * arg);
  void * areConnectedBBFS(void * arg);
  void * areConnectedNoLabels(void * arg);

  // Declare launcher classes as friends
  friend class DFSlauncher;
  friend class BBFSlauncher;

  // Maintenance
  bool addEdgeUnsafe(Vertex * source, Vertex * target);

#ifdef ENABLE_STATISTICS
  // Internal statistic maintenance
  void registerQueryStatistics(bool result, unsigned int pathSize, uintmax_t searchedNodes);
#endif // ENABLE_STATISTICS
};


class DFSlauncher {
public:
  void * runArg;
  Graph * myGraph;

  DFSlauncher(Graph * graph, void * arg);

  void * runDFS(void);
};


class BBFSlauncher {
public:
  void * runArg;
  Graph * myGraph;

  BBFSlauncher(Graph * graph, void * arg);

  void * runBBFS(void);
};


void * launchDFS(void * helper);
void * launchBBFS(void * helper);


#endif // GRAPH_HPP
