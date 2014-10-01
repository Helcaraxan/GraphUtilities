#include <stack>
#include <vector>
#include <iostream>

#include "graph.hpp"

using namespace std;


/* Constructors & destructors */
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
	queryID(0)
{}


Graph::~Graph(void)
{
	for (auto it = vertices.begin(); it != vertices.end(); ++it)
		delete *it;

	vertices.clear();
}


// Modificators
bool Vertex::addPredecessor(Vertex * pred)
{
	for (auto it = predecessors.begin(); it != predecessors.end(); ++it) {
		if (*it == pred)
			return false;
	}
	
	predecessors.push_back(pred);
	return true;
}


bool Vertex::addSuccessor(Vertex * succ)
{
	for (auto it = successors.begin(); it != successors.end(); ++it) {
		if (*it == succ)
			return false;
	}

	successors.push_back(succ);
	return true;
}


Vertex * Graph::addVertex(void)
{
	int id = vertices.size();;
	Vertex * newVertex = new Vertex(id);

	vertices[id] = newVertex;

	return newVertex;
}


// Indexing
void Vertex::visit(Vertex * pred, bool reverse)
{
	vector<Vertex *> * orderedVertices;

	// Reinitialize in case of reverse post-order traversal
	if (reverse && ((inVisits + outVisits) == (int) (successors.size() + predecessors.size()))) {
		inVisits = 0;
		outVisits = 0;
	}

	orderedVertices = reverse ? &successors : &predecessors;

	if (!pred || ((*orderedVertices)[inVisits] == pred))
		return;

	for (int i = inVisits; i < (int) orderedVertices->size(); i++) {
		if ((*orderedVertices)[i] == pred) {
			(*orderedVertices)[i] = (*orderedVertices)[inVisits];
			(*orderedVertices)[inVisits++] = pred;
			break;
		}
	}
}


Vertex * Vertex::createPostOrder(stack<Vertex *> * postOrder, bool reverse)
{
	vector<Vertex *> * upVertices;
	vector<Vertex *> * downVertices;

	upVertices = reverse ? &successors : &predecessors;
	downVertices = reverse ? &predecessors : &successors;

	while (outVisits < (int) downVertices->size()) {
		(*downVertices)[outVisits]->visit(this, reverse);

		if ((*downVertices)[outVisits]->inVisits == 1)
			return (*downVertices)[outVisits++];
		else
			outVisits++;
	}

	postOrder->push(this);

	if (upVertices->size() == 0)
		return NULL;
	else
		return (*upVertices)[0];
}


void Graph::labelVertices(bool reverse)
{
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
		while (!postOrder.empty())
			postOrder.pop();

		// Run the DFS
		while (nextVertex) {
			// In case of regular labeling create the sinks repository
			if (!reverse && (nextVertex->successors.size() == 0))
				sinks.push_back(nextVertex);

			// Use an iterative method to prevent memory overflow in large graphs
			nextVertex = nextVertex->createPostOrder(&postOrder, reverse);
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
}


void Graph::indexGraph()
{
	// Create the vector with the source nodes
	for (auto it = vertices.begin(); it != vertices.end(); ++it) {
		if ((*it)->predecessors.size() == 0)
			sources.push_back(*it);
	}

	// Perform the forward graph ordering
	labelVertices(false);

	// Perform the reverse graph ordering
	labelVertices(true);
}


// Queries
bool Graph::areConnected(Vertex * u, Vertex * v) {
	Vertex * curr;
	stack<Vertex *> searchStack;

	// Are U and V the same vertex?
	if (u == v)
		return true;

	// Can V be a descendant of U in the standard graph?
	if (u->orderLabel > v->orderLabel)
		return false;

	// Can U be a descendant of V in the reverse graph?
	if (v->reverseOrderLabel > u->reverseOrderLabel)
		return false;

	// Do a DFS on the subgraph specified by both orders to get the final answer
	searchStack.push(u);
	queryID++;

	while (!searchStack.empty()) {
		curr = searchStack.top();
		searchStack.pop();
		if (curr == v)
			return true;

		for (auto it = curr->successors.begin(); it !=curr->successors.end(); ++it) {
			if ((*it)->queryID != queryID) {
				(*it)->queryID = queryID;

				if (((*it)->orderLabel < v->orderLabel) && ((*it)->reverseOrderLabel < u->reverseOrderLabel))
					searchStack.push(*it);
			}
		}
	}

	return false;
}


