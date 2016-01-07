#include <fstream>
#include <iostream>

#include "graph-utilities/implementation/graphImpl.hpp"

using namespace std;


// Internal functions for graph condensation

static int
sccVertexVisit(GraphImpl * graph, VertexImpl * target, int label,
    vector<pair<int, int> > &labels, VertexImpl::Stack &unclassified,
    VertexImpl::Stack &undetermined) {
  Vertex::PartitionArray& sccSets = graph->getSCCsI();

  if (labels[target->getId()].first != -1)
    cerr << "Revisited node " << target->getId()
      << " during graph condensation!\n";

  labels[target->getId()].first = label++;
  unclassified.push(target);
  undetermined.push(target);

  for (auto it = target->succBegin(), end = target->succEnd();
      it != end; ++it) {
    if (labels[(*it)->getId()].first == -1)
      label =
        sccVertexVisit(graph, *it, label, labels, unclassified, undetermined);
    else if (labels[(*it)->getId()].second == -1)
      while (labels[undetermined.top()->getId()].first
          >= labels[(*it)->getId()].first)
        undetermined.pop();
  }

  if (target == undetermined.top()) {
    int sccId = sccSets.size();
    sccSets.push_back(new Vertex::IdSet);

    while (unclassified.top() != target) {
      sccSets.back()->insert(unclassified.top()->getId());
      labels[unclassified.top()->getId()].second = sccId;
      unclassified.pop();
    }

    sccSets.back()->insert(unclassified.top()->getId());
    labels[unclassified.top()->getId()].second = sccId;
    unclassified.pop();

    undetermined.pop();
  }

  return label;
}


static bool
checkDAG(GraphImpl * graph, bool condense, const char * correspondanceFile) {
  int label = 0;
  bool result = true;
  VertexImpl::Stack unclassified;
  VertexImpl::Stack undetermined;
  vector<pair<int, int> >
    labels(graph->getVertexCount(), pair<int, int>(-1, -1));
  vector<int> sccLabels(graph->getVertexCount(), -1);
  fstream outStream;
  Vertex::PartitionArray& sccSets = graph->getSCCsI();

  // Erase any existing SCCs
  for (auto it = sccSets.begin(), end = sccSets.end(); it != end; ++it)
    delete *it;

  sccSets.clear();
  graph->discoverExtremities();
  VertexImpl::Array& sources = graph->getSourceVertices();

  // Loop over all vertices and apply the Path-based SCC algorithm
  for (auto it = sources.begin(), end = sources.end(); it != end; ++it) {
    if (labels[(*it)->getId()].first != -1)
      continue;

    label =
      sccVertexVisit(graph, *it, label, labels, unclassified, undetermined);

    if (!unclassified.empty())
      cerr << "ERROR : unclassified stack was not empty." << endl;

    if (!undetermined.empty())
      cerr << "ERROR : undetermined stack was not empty." << endl;
  }

  // Verify that all valid vertices have been visited
  if (label != (int) graph->getVertexCount())
    cerr << "Did not visit all vertices during graph condensation!" << endl;

  // When a filename is provided prepare the stream to dump the correspondance
  // between original and merged nodes
  if (condense && (correspondanceFile != NULL))
    outStream.open(correspondanceFile, ios_base::out);

  // Merge the SCCs that have more than one vertex and dump if necessary
  for (auto it = sccSets.begin(), end = sccSets.end(); it != end; ++it) {
    if ((*it)->size() > 1) {
      if (condense) {
        Vertex::Set mergeSet;

        if (correspondanceFile != NULL)
          outStream << graph->getVertexCount() << ":";

        for (auto setIt = (*it)->begin(), setEnd = (*it)->end();
            setIt != setEnd; ++setIt) {
          mergeSet.insert(graph->getVertex(*setIt));

          if (correspondanceFile != NULL)
            outStream << " " << *setIt;
        }

        if (correspondanceFile != NULL)
          outStream << "\n";

        result = false;
        
        graph->mergeVertices(mergeSet);
      } else {
        return false;
      }
    }
  }

  if (condense && (correspondanceFile != NULL))
    outStream.close();

  return result;
}


// Graph condensation

void
GraphImpl::setCondensed(bool value) {
  condensed = value;
}


bool
GraphImpl::isCondensed() const {
  return condensed;
}


const Vertex::PartitionArray&
GraphImpl::getSCCs() {
  return sccSets;
}


Vertex::PartitionArray&
GraphImpl::getSCCsI() {
  return sccSets;
}


bool
GraphImpl::isDAG() {
  return checkDAG(this, false, NULL);
}


void
GraphImpl::condenseToDAG(const char * correspondanceFile) {
  if (!condensed)
    checkDAG(this, true, correspondanceFile);
  
  condensed = true;
}
