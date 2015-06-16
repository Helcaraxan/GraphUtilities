#include <atomic>
#include <vector>

using namespace std;


// Forward declarations

class Vertex;
class Graph;


// Enumeration types

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
  UndefinedIndexMethod = 0x10
} IndexMethod;


typedef enum {
  DFS,
  BBFS,
  NoLabels,
  UndefinedSearchMethod
} SearchMethod;


typedef enum {
  UndefinedPartitionMethod
} PartitionMethod;


typedef enum {
  Reachability,
  Partition
} QueryType;



// The Reachability variant of the Query class

class ReachabilityQuery {
friend class Graph;

// Data members
private:
  atomic<Vertex *> source;
  atomic<Vertex *> target;
  atomic<bool> answer{false};
  atomic<bool> error{false};
  atomic<bool> cancel{false};
  atomic<uint64_t> internal{0};
  atomic<SearchMethod> method{UndefinedSearchMethod};

#ifdef ENABLE_STATISTICS
  vector<Vertex *> path;
  atomic<uintmax_t> searchedNodes{0};
#endif // ENABLE_STATISTICS


// Function members
public:
  ReachabilityQuery(Vertex * u, Vertex * v, SearchMethod searchMethod = UndefinedSearchMethod) {
    source.store(u, memory_order_release);
    target.store(v, memory_order_release);
    method.store(searchMethod, memory_order_release);
  }

  Vertex * getSource(void);
  Vertex * getTarget(void);
  bool getAnswer(void);
  bool getError(void);
  SearchMethod getMethod(void);


private:
  void setAnswer(bool value);
  void setError(bool value);
  void setCancel(bool value);
  void setInternal(uint64_t value);
  void setMethod(SearchMethod value);
  bool getCancel(void);
  uint64_t getInternal(void);
};



// The Partition variant of the Query class

class PartitionQuery {
friend class Graph;

// Data members
private:


// Function members
public:
  PartitionQuery() {
  }
};



// The Query class

class Query {
// Local types
public:
  union QueryUnion {
    ReachabilityQuery &reachability;
    PartitionQuery &partition;

    QueryUnion(ReachabilityQuery &q) :
      reachability(q)
    {}

    QueryUnion(PartitionQuery &q) :
      partition(q)
    {}
  };


// Data members
public:
  const QueryType type;
  const QueryUnion query;


// Function members
public:
  Query(QueryType t, ReachabilityQuery &q) :
    type(t),
    query(q)
  {}

  Query(QueryType t, PartitionQuery &q) :
    type(t),
    query(q)
  {}
};
