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

using namespace std;

class Vertex;
class Graph;


/* Class declarations */
class Vertex {
friend class Graph;

public:
	int id; /* A unique Vertex ID */

	// Constructors & Destructor
	Vertex(int i);
	~Vertex(void);

  // Local types
  typedef vector<Vertex *>::iterator iterator;

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


protected:
	int orderLabel; /* Ordering in a reverse post-order traversal of the Graph */
	int reverseOrderLabel; /* Ordering in a reverse post-order traversal of the inverse Graph */

	int inVisits; /* Variable used in reverse post-order */
	int outVisits; /* Variable used in reverse post-order */
	int queryID; /* Variable to distinguish searches for seperate reachability queries */

	vector<Vertex *> predecessors;
	vector<Vertex *> successors;

	// Indexing
	void visit(Vertex * pred, bool reverse); /* Function used in reverse post-order traversal */
	Vertex * createPostOrder(stack<Vertex *> * postOrder, bool reverse); /* Idem */
};


class Graph {
public:
	// Constructors & Destructor
	Graph(void);
	~Graph(void);

  // Parser function
  static Graph * createFromDotFile(const char * fileName);

	// Modificators
	Vertex * addVertex(void);
	bool addEdge(Vertex * s, Vertex * t);
  bool removeEdge(Vertex * s, Vertex * t);

  // Access
  int getVertexCount(void);
  Vertex * getVertexFromId(int id);

	// Queries
	bool areConnected(Vertex * u, Vertex * v); /* The query function */
  bool indirectPathExists(Vertex * u, Vertex * v); /* Find if two vertices are connected through a multi-hop path */

protected:
	bool indexed; /* Indicates if the ordering of the graph has been done */

	vector<Vertex *> vertices;
	vector<Vertex *> sources;
	vector<Vertex *> sinks;

	int queryID; /* Global ID to distinguish between seperate reachability queries */

	// Indexing
	void labelVertices(bool reverse);
	void indexGraph(void);
};

#endif // GRAPH_HPP
