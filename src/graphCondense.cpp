#include <fstream>
#include <iostream>

#include "graph-utilities/defs.hpp"
#include "graph-utilities/graph.hpp"
#include "graph-utilities/vertex.hpp"

using namespace std;


// Graph condensation


bool
Graph::condenseGraph(bool dummy, const char * dumpFile) {
  int label = 0;
  bool result = false;
  stack<Vertex *> unclassified;
  stack<Vertex *> undetermined;
  vector<pair<int, int> > labels(vertices.size(), pair<int, int>(-1, -1));
  vector<int> sccLabels(vertices.size(), -1);
  fstream outStream;

  // Erase any existing SCCs
  for (auto it = sccSets.begin(), end = sccSets.end(); it != end; ++it)
    delete *it;

  sccSets.clear();
  discoverExtremities();

  // Loop over all vertices and apply the Path-based SCC algorithm
  for (auto it = sources.begin(), end = sources.end(); it != end; ++it) {
    if ((*it == NULL) || (labels[(*it)->id].first != -1))
      continue;

    label = sccVertexVisit(*it, label, labels, unclassified, undetermined);

    if (!unclassified.empty())
      cerr << "ERROR : unclassified stack was not empty." << endl;

    if (!undetermined.empty())
      cerr << "ERROR : undetermined stack was not empty." << endl;
  }

  // Verify that all valid vertices have been visited
  if (label != (int) getVertexCount())
    cerr << "Did not visit all vertices during graph condensation!" << endl;

  // When a filename is provided prepare the stream to dump the correspondance
  // between original and merged nodes
  if (dumpFile)
    outStream.open(dumpFile, ios_base::out);

  // Merge the SCCs that have more than one vertex and dump if necessary
  for (auto it = sccSets.begin(), end = sccSets.end(); it != end; ++it) {
    if ((*it)->size() > 1) {
      if (!dummy) {
        if (dumpFile) {
          outStream << vertexCount << ":";
          for (auto setIt = (*it)->begin(), setEnd = (*it)->end();
              setIt != setEnd; ++setIt)
            outStream << " " << (*setIt)->id; 

          outStream << "\n";
        }

        mergeVertices(**it);
      }

      result = true;
    }
  }

  if (dumpFile)
    outStream.close();

  return result;
}


// Access

graph_t *
Graph::getMetisGraph() {
  int adjIdx = 0;
  graph_t * metisGraph = new graph_t;

  metisGraph->nvtxs = getVertexCount();
  metisGraph->xadj = new mtmetis_adj_t[metisGraph->nvtxs];
  metisGraph->adjncy = new mtmetis_vtx_t[getEdgeCount()];

  metisGraph->xadj[0] = adjIdx;

  for (auto it = vertices.begin(), end = vertices.end(); it != end; ++it) {
    if (*it == NULL)
      continue;

    for (auto it2 = (*it)->pred_begin(), end2 = (*it)->pred_end();
        it2 != end2; ++it2)
      metisGraph->adjncy[adjIdx++] = (*it)->id;

    for (auto it2 = (*it)->succ_begin(), end2 = (*it)->succ_end();
        it2 != end2; ++it2)
      metisGraph->adjncy[adjIdx++] = (*it)->id;

    metisGraph->xadj[(*it)->id + 1] = adjIdx;
  }

  return metisGraph;
}


// Internal condensation function

int
Graph::sccVertexVisit(Vertex * target, int label,
    vector<pair<int, int> > &labels,
    stack<Vertex *> &unclassified, stack<Vertex *> &undetermined) {
  if (labels[target->id].first != -1)
    cerr << "Revisited node " << target->id << " during graph condensation!\n";

  labels[target->id].first = label++;
  unclassified.push(target);
  undetermined.push(target);

  for (auto it = target->succ_begin(), end = target->succ_end();
      it != end; ++it) {
    if (labels[(*it)->id].first == -1) {
      label = sccVertexVisit(*it, label, labels, unclassified, undetermined);
    } else if (labels[(*it)->id].second == -1) {
      while (labels[undetermined.top()->id].first >= labels[(*it)->id].first)
        undetermined.pop();
    }
  }

  if (target == undetermined.top()) {
    int sccId = sccSets.size();
    sccSets.push_back(new set<Vertex *>);

    while (unclassified.top() != target) {
      sccSets.back()->insert(unclassified.top());
      labels[unclassified.top()->id].second = sccId;
      unclassified.pop();
    }

    sccSets.back()->insert(unclassified.top());
    labels[unclassified.top()->id].second = sccId;
    unclassified.pop();

    undetermined.pop();
  }

  return label;
}
