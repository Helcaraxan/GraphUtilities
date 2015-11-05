#include "graph-utilities/vertex.hpp"

using namespace std;


// Global DFSId variable in case of TLS
#ifdef ENABLE_TLS
thread_local vector<int> DFSId;
#endif // ENABLE_TLS


/* Public member functions */

// Modificators

bool
Vertex::addPredecessor(Vertex * pred, int weight) {
  if (pred == this)
    return false;

  for (auto it = pred_begin(), end = pred_end(); it != end; ++it) {
    if (*it == pred)
      return false;
  }

  addPredecessorUnsafe(pred, weight);
  return true;
}


bool
Vertex::addSuccessor(Vertex * succ, int weight) {
  if (succ == this)
    return false;

  for (auto it = succ_begin(), end = succ_end(); it != end; ++it) {
    if (*it == succ)
      return false;
  }

  addSuccessorUnsafe(succ, weight);

  return true;
}


bool
Vertex::removePredecessor(Vertex * pred) {
  auto weightIt = predecessorWeights.begin();
  for (auto predIt = pred_begin(), end = pred_end(); predIt != end; ++predIt) {
    if (*predIt == pred) {
      predecessors.erase(predIt);
      predecessorWeights.erase(weightIt);
      predecessorCount--;
      return true;
    }

    ++weightIt;
  }

  return false;
}


bool
Vertex::removeSuccessor(Vertex * succ) {
  auto weightIt = successorWeights.begin();
  for (auto succIt = succ_begin(), end = succ_end(); succIt != end; ++succIt) {
    if (*succIt == succ) {
      successors.erase(succIt);
      successorWeights.erase(weightIt);
      successorCount--;
      return true;
    }

    ++weightIt;
  }

  return false;
}


void
Vertex::clearPredecessors() {
  predecessors.clear();
  predecessorWeights.clear();
  predecessorCount = 0;
}


void
Vertex::clearSuccessors() {
  successors.clear();
  successorWeights.clear();
  successorCount = 0;
}


// Access

Vertex *
Vertex::getPredecessor(int idx) const {
  return predecessors.at(idx);
}


Vertex *
Vertex::getSuccessor(int idx) const {
  return successors.at(idx);
}


int
Vertex::getPredecessorCount() const {
  return predecessorCount;
}


int
Vertex::getSuccessorCount() const {
  return successorCount;
}


int
Vertex::getPredecessorWeight(Vertex * pred) const {
  auto weightIt = predecessorWeights.begin();
  for (auto predIt = pred_begin(), end = pred_end(); predIt != end; ++predIt) {
    if (*predIt == pred)
      return *weightIt;

    ++weightIt;
  }

  return -1;
}


int
Vertex::getSuccessorWeight(Vertex * succ) const {
  auto weightIt = successorWeights.begin();
  for (auto succIt = succ_begin(), end = succ_end(); succIt != end; ++succIt) {
    if (*succIt == succ)
      return *weightIt;

    ++weightIt;
  }

  return -1;
}


int
Vertex::getPredecessorWeight(int idx) const {
  if (getPredecessorCount() <= idx)
    return -1;

  return predecessorWeights[idx];
}


int
Vertex::getSuccessorWeight(int idx) const {
  if (getSuccessorCount() <= idx)
    return -1;

  return successorWeights[idx];
}


// Iterators

Vertex::const_iterator
Vertex::pred_begin() const {
  return predecessors.begin();
}


Vertex::const_iterator
Vertex::pred_end() const {
  return predecessors.end();
}


Vertex::const_iterator
Vertex::succ_begin() const {
  return successors.begin();
}


Vertex::const_iterator
Vertex::succ_end() const {
  return successors.end();
}


Vertex::iterator
Vertex::pred_begin() {
  return predecessors.begin();
}


Vertex::iterator
Vertex::pred_end() {
  return predecessors.end();
}


Vertex::iterator
Vertex::succ_begin() {
  return successors.begin();
}


Vertex::iterator
Vertex::succ_end() {
  return successors.end();
}


// User data management (local)

void
Vertex::setUserData(UserDataInterface * newData) {
  if (userData != NULL)
    delete userData;

  userData = newData;
}


UserDataInterface *
Vertex::getUserData() const {
  return userData;
}


/* Private member functions */

// Indexing

#ifdef ENABLE_TLS
void
Vertex::setDFSId(uint64_t newId) {
  DFSId.at(id) = newId;
}


uint64_t
Vertex::getDFSId() const {
  return DFSId.at(id);
}


void
Vertex::checkDFSId(size_t size) {
  if (DFSId.size() < size)
    DFSId.resize(size);
}

#else // ENABLE_TLS

void
Vertex::setDFSId(int idx, uint64_t newId) {
  DFSId[idx] = newId;
}


uint64_t
Vertex::getDFSId(int idx) const {
  return DFSId[idx];
}
#endif // ENABLE_TLS


/* Compute the next Vertex to be visited in a DFS traversal of the graph
 * from the current Vertex instance. This allows for an iterative traversal and
 * not a recursive one in order to prevent stack-overflows on large graphs. The
 * result pair gives the potential next vertice in the DFS traversal and the
 * boolean telling if its the real next vertice or if it's an intermediary one.
 */
pair<bool, Vertex *>
Vertex::getNextDFS(Vertex * origin, bool postOrder, bool reverse) {
  bool isTarget = false;
  Vertex * target = NULL;

  if (inVisits == 0) {
    // This is part of the first visit to this Vertex

    if (firstVisit == NULL) {
      // Register the origin as the source of the DFS traversal
      if (origin == this) {
        firstVisit = (Vertex *) 0x1;
      } else {
        firstVisit = origin;
      }

      target = this;
      if (!postOrder)
        isTarget = true;
      else
        isTarget = false;
    } else if (outVisits < (reverse ? predecessorCount : successorCount)) {
      // Visit successors if some have not yet been visited
      isTarget = false;
      target = reverse ? predecessors[outVisits++] : successors[outVisits++];
    } else {
      // All successors have been visited so either...

      if ((!postOrder)
          || (outVisits == (reverse ? predecessorCount : successorCount))) {
        // ... go back to the first visiting Vertex for this traversal
        if (firstVisit == (Vertex *) 0x1) {
          inVisits = 0;
          outVisits = 0;
          firstVisit = NULL;

          isTarget = true;
          target = NULL;
        } else {
          inVisits = 1;

          isTarget = false;
          target = firstVisit;
        }
      } else {
        // ... or return this Vertex as target ...
        outVisits++;

        isTarget = true;
        target = this;
      }
    }
  } else {
    // Register a visit to this Vertex but bounce back to the origin
    inVisits++;

    isTarget = false;
    target = origin;
  }

  // If this was the last visit to a non source then reinitialize the Vertex's
  // counters
  if ((firstVisit != (Vertex *) 0x1) &&
        (inVisits == (reverse ? successorCount : predecessorCount))) {
    inVisits = 0;
    outVisits = 0;
    firstVisit = NULL;
  }

  return pair<bool, Vertex *>(isTarget, target);
}


// Modificators

void
Vertex::addPredecessorUnsafe(Vertex * pred, int weight) {
  predecessors.push_back(pred);
  predecessorWeights.push_back(weight);
  predecessorCount++;
}


void
Vertex::addSuccessorUnsafe(Vertex * succ, int weight) {
  successors.push_back(succ);
  successorWeights.push_back(weight);
  successorCount++;
}
