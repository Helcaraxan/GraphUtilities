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


/* Class declarations */
class Vertex {
public:
	int id;
	int orderLabel;
	int reverseOrderLabel;

	int inVisits;
	int outVisits;

	int queryID;

	vector<Vertex *> predecessors;
	vector<Vertex *> successors;

	// Constructors & Destructor
	Vertex(int i);
	~Vertex(void);

	// Modificators
	bool addPredecessor(Vertex * pred);
	bool addSuccessor(Vertex * succ);

	// Indexing
	void visit(Vertex * pred, bool reverse);
	Vertex * createPostOrder(stack<Vertex *> * postOrder, bool reverse);
};


class Graph {
public:
	bool indexed;
	int queryID;

	vector<Vertex *> vertices;
	vector<Vertex *> sources;
	vector<Vertex *> sinks;

	// Constructors & Destructor
	Graph(void);
	~Graph(void);

	// Modificators
	Vertex * addVertex(void);
	bool addEdge(Vertex * s, Vertex * t);

	// Indexing
	void labelVertices(bool reverse);
	void indexGraph(void);

	// Queries
	bool areConnected(Vertex * u, Vertex * v);
};

#endif // GRAPH_HPP
