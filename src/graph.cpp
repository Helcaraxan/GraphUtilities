#include <stack>
#include <vector>
#include <cstring>
#include <fstream>
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


// Indexing

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
void
Graph::fillFromDotFile(const char * fileName) {
  char dump[128];
  int source, target;
  fstream input(fileName, fstream::in);
  
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
      return;
    }

    if (strchr(dump, '>')) {
      sscanf(dump, "%d -> %d", &source, &target);
      while (vertices.size() <= (unsigned) target)
        addVertex();

      addEdge(vertices[source], vertices[target]);
    } else {
      sscanf(dump, "%d", &source);
      while (vertices.size() <= (unsigned) source)
        addVertex();
    }
  }
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
	return true;
}


bool
Graph::removeEdge(Vertex * source, Vertex * target) {
  if (!source->removeSuccessor(target))
    return false;

  source->removePredecessor(source);
  return true;
}


// Access
Vertex *
Graph::getVertexFromId(int id) {
  return vertices[id];
}


// Queries

/* Use the previously done indexation to answer to the query */
bool Graph::areConnected(Vertex * u, Vertex * v) {
	Vertex * curr;
	stack<Vertex *> searchStack;

  // Verify that the graph has been indexed
  if (!indexed)
    indexGraph();

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

		for (auto it = curr->successors.begin(); it !=curr->successors.end(); ++it) {
      if (*it == v)
        return true;

			if ((*it)->queryID != queryID) {
				(*it)->queryID = queryID;

				if (((*it)->orderLabel < v->orderLabel) && ((*it)->reverseOrderLabel > v->reverseOrderLabel))
					searchStack.push(*it);
			}
		}
	}

	return false;
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
}

