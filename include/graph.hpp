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

	// Modificators
	bool addPredecessor(Vertex * pred);
	bool addSuccessor(Vertex * succ);
  bool removePredecessor(Vertex * pred);
  bool removeSuccessor(Vertex * succ);

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
	bool indexed; /* Indicates if the ordering of the graph has been done */

	// Constructors & Destructor
	Graph(void);
	~Graph(void);

  // Parser function
  void fillFromDotFile(const char * fileName);

	// Modificators
	Vertex * addVertex(void);
	bool addEdge(Vertex * s, Vertex * t);
  bool removeEdge(Vertex * s, Vertex * t);

  // Access
  Vertex * getVertexFromId(int id);

	// Queries
	bool areConnected(Vertex * u, Vertex * v); /* The query function */

protected:
	vector<Vertex *> vertices;
	vector<Vertex *> sources;
	vector<Vertex *> sinks;

	int queryID; /* Global ID to distinguish between seperate reachability queries */

	// Indexing
	void labelVertices(bool reverse);
	void indexGraph(void);
};

#endif // GRAPH_HPP
