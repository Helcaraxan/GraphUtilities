#include "graph-utilities/implementation/vertexImpl.hpp"

using namespace std;


// Global DFSId variable in case of TLS
#ifdef ENABLE_TLS
thread_local vector<int> DFSId;
#endif // ENABLE_TLS


// UserData implementation
template<class UserDataImpl>
UserDataInterface *
UserData<UserDataImpl>::clone() const {
  return new UserDataImpl(static_cast<const UserDataImpl&>(*this));
}


// User data management

void
VertexImpl::setUserData(UserDataInterface * newData) {
  if (userData != NULL)
    delete userData;

  userData = newData;
}


UserDataInterface *
VertexImpl::getUserData() const {
  return userData;
}


// Indexing

#ifdef ENABLE_TLS
void
VertexImpl::setDFSId(uint64_t newId) {
  DFSId.at(id) = newId;
}


uint64_t
VertexImpl::getDFSId() const {
  return DFSId.at(id);
}


void
VertexImpl::checkDFSId(size_t size) {
  if (DFSId.size() < size)
    DFSId.resize(size);
}

#else // ENABLE_TLS

void
VertexImpl::setDFSId(int index, uint64_t newId) {
  DFSId[index] = newId;
}


uint64_t
VertexImpl::getDFSId(int index) const {
  return DFSId[index];
}
#endif // ENABLE_TLS


pair<int, int>
VertexImpl::getOrders() const {
  return pair<int, int>(orderLabel, revOrderLabel);
}


// Graph traversal

void
VertexImpl::resetDFS() {
  inVisits = 0;
  outVisits = 0;
  firstVisit = NULL;
}

/* Compute the next Vertex to be visited in a DFS traversal of the graph
 * from the current Vertex instance. This allows for an iterative traversal and
 * not a recursive one in order to prevent stack-overflows on large graphs. The
 * result pair gives the potential next vertice in the DFS traversal and the
 * boolean telling if its the real next vertice or if it's an intermediary one.
 */
pair<bool, VertexImpl *>
VertexImpl::getNextDFS(VertexImpl * origin, Graph::TraversalOrder order,
    Graph::TraversalDirection direction) {
  int visitLimit = 0;
  bool isTarget = false;
  VertexImpl * target = NULL;

  if (inVisits == 0) {
    if (direction == Graph::TraversalDirection::Forward)
      visitLimit = successorCount;
    else if (direction == Graph::TraversalDirection::Backward)
      visitLimit = predecessorCount;

    // This is part of the first visit to this Vertex
    if (firstVisit == NULL) {
      // Register the origin as the source of the DFS traversal
      if (origin == this) {
        firstVisit = (VertexImpl *) 0x1;
      } else {
        firstVisit = origin;
      }

      target = this;
      if (order == Graph::TraversalOrder::PreOrder)
        isTarget = true;
      else if (order == Graph::TraversalOrder::PostOrder)
        isTarget = false;
    } else if (outVisits < visitLimit) {
      // Visit successors if some have not yet been visited
      isTarget = false;

      if (direction == Graph::TraversalDirection::Forward)
        target = successors[outVisits++];
      else if (direction == Graph::TraversalDirection::Backward)
        target = predecessors[outVisits++];
    } else {
      // All successors have been visited so either...
      if ((order == Graph::TraversalOrder::PreOrder)
          || (outVisits == visitLimit)) {
        // ... go back to the first visiting Vertex for this traversal
        if (firstVisit == (Vertex *) 0x1) {
          resetDFS();

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
  if (direction == Graph::TraversalDirection::Forward)
    visitLimit = predecessorCount;
  else if (direction == Graph::TraversalDirection::Backward)
    visitLimit = successorCount;

  if ((firstVisit != (Vertex *) 0x1) && (inVisits == visitLimit))
    resetDFS();

  return pair<bool, VertexImpl *>(isTarget, target);
}


// Modificators

void
VertexImpl::setWeight(int newWeight) {
  weight = newWeight;
}


bool
VertexImpl::addPredecessor(Vertex * predecessor, int weight) {
  if (predecessor == this)
    return false;

  for (auto it = predBegin(), end = predEnd(); it != end; ++it) {
    if ((Vertex *) *it == predecessor)
      return false;
  }

  addPredecessorUnsafe(dynamic_cast<VertexImpl *>(predecessor), weight);
  return true;
}


bool
VertexImpl::addSuccessor(Vertex * successor, int weight) {
  if (successor == this)
    return false;

  for (auto it = succBegin(), end = succEnd(); it != end; ++it) {
    if ((Vertex *) *it == successor)
      return false;
  }

  addSuccessorUnsafe(dynamic_cast<VertexImpl *>(successor), weight);

  return true;
}


bool
VertexImpl::removePredecessor(Vertex * predecessor) {
  auto weightIt = predecessorWeights.begin();
  for (auto predIt = predBegin(), end = predEnd(); predIt != end; ++predIt) {
    if ((Vertex *) *predIt == predecessor) {
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
VertexImpl::removeSuccessor(Vertex * successor) {
  auto weightIt = successorWeights.begin();
  for (auto succIt = succBegin(), end = succEnd(); succIt != end; ++succIt) {
    if ((Vertex *) *succIt == successor) {
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
VertexImpl::clearPredecessors() {
  predecessors.clear();
  predecessorWeights.clear();
  predecessorCount = 0;
}


void
VertexImpl::clearSuccessors() {
  successors.clear();
  successorWeights.clear();
  successorCount = 0;
}


void
VertexImpl::setPredecessorWeight(const Vertex * predecessor, int weight) {
  auto weightIt = predecessorWeights.begin();
  for (auto it = predecessors.begin(), end = predecessors.end();
      it != end; ++it) {
    if (*it == predecessor) {
      *weightIt = weight;
      break;
    }

    weightIt++;
  }
}


void
VertexImpl::setSuccessorWeight(const Vertex * successor, int weight) {
  auto weightIt = successorWeights.begin();
  for (auto it = successors.begin(), end = successors.end();
      it != end; ++it) {
    if (*it == successor) {
      *weightIt = weight;
      break;
    }

    weightIt++;
  }
}


void
VertexImpl::addPredecessorUnsafe(VertexImpl * predecessor, int weight) {
  predecessors.push_back(predecessor);
  predecessorWeights.push_back(weight);
  predecessorCount++;
}


void
VertexImpl::addSuccessorUnsafe(VertexImpl * successor, int weight) {
  successors.push_back(successor);
  successorWeights.push_back(weight);
  successorCount++;
}


void
VertexImpl::setOrders(int forward, int backward) {
  if (forward != -1)
    orderLabel = forward;
  
  if (backward != -1)
    revOrderLabel = backward;
}


// Access

int
VertexImpl::getId() const {
  return id;
}


int
VertexImpl::getWeight() const {
  return weight;
}


Vertex *
VertexImpl::getPredecessor(int index) const {
  return predecessors.at(index);
}


Vertex *
VertexImpl::getSuccessor(int index) const {
  return successors.at(index);
}


int
VertexImpl::getPredecessorCount() const {
  return predecessorCount;
}


int
VertexImpl::getSuccessorCount() const {
  return successorCount;
}


int
VertexImpl::getPredecessorWeight(const Vertex * predecessor) const {
  auto weightIt = predecessorWeights.begin();
  for (auto predIt = predBegin(), end = predEnd(); predIt != end; ++predIt) {
    if ((Vertex *) *predIt == predecessor)
      return *weightIt;

    ++weightIt;
  }

  return -1;
}


int
VertexImpl::getSuccessorWeight(const Vertex * successor) const {
  auto weightIt = successorWeights.begin();
  for (auto succIt = succBegin(), end = succEnd(); succIt != end; ++succIt) {
    if ((Vertex *) *succIt == successor)
      return *weightIt;

    ++weightIt;
  }

  return -1;
}


int
VertexImpl::getPredecessorWeight(int index) const {
  if (getPredecessorCount() <= index)
    return -1;

  return predecessorWeights[index];
}


int
VertexImpl::getSuccessorWeight(int index) const {
  if (getSuccessorCount() <= index)
    return -1;

  return successorWeights[index];
}


VertexImpl *
VertexImpl::getPredecessorI(int index) const {
  return predecessors[index];
}


VertexImpl *
VertexImpl::getSuccessorI(int index) const {
  return successors[index];
}


bool
VertexImpl::hasPredecessor(const VertexImpl * predecessor) const {
  for (auto it = predecessors.begin(), end = predecessors.end();
      it != end; ++it) {
    if (*it == predecessor)
      return true;
  }

  return false;
}


bool
VertexImpl::hasSuccessor(const VertexImpl * successor) const {
  for (auto it = successors.begin(), end = successors.end();
      it != end; ++it) {
    if (*it == successor)
      return true;
  }

  return false;
}


// Iterators

VertexImpl::const_iterator
VertexImpl::predBegin() const {
  return predecessors.begin();
}


VertexImpl::const_iterator
VertexImpl::predEnd() const {
  return predecessors.end();
}


VertexImpl::const_iterator
VertexImpl::succBegin() const {
  return successors.begin();
}


VertexImpl::const_iterator
VertexImpl::succEnd() const {
  return successors.end();
}


VertexImpl::iterator
VertexImpl::predBegin() {
  return predecessors.begin();
}


VertexImpl::iterator
VertexImpl::predEnd() {
  return predecessors.end();
}


VertexImpl::iterator
VertexImpl::succBegin() {
  return successors.begin();
}


VertexImpl::iterator
VertexImpl::succEnd() {
  return successors.end();
}


/* PartionNode Interface implementation */
// Access methods

const PartitionNode *
PartitionNodeImpl::getParent() const {
  return parent;
}


int
PartitionNodeImpl::getChildCount() const {
  return children.size();
}


const PartitionNode *
PartitionNodeImpl::getChild(int idx) const {
  return children[idx];
}


/* Partion Interface implementation */
// Access methods

void
PartitionImpl::represents(const PartitionNode * node,
    Vertex::IdSet& idSet) const {
  const PartitionNodeImpl * iNode =
    dynamic_cast<const PartitionNodeImpl *>(node);

  if (!iNode)
    return;

  idSet.insert(iNode->id);
  IaP<PartitionNodeImpl> ptr = representants[iNode->id];

  while (ptr.isInteger()) {
    idSet.insert(ptr);
    ptr = representants[ptr];
  }
}


const PartitionNode *
PartitionImpl::getRoot() const {
  return root;
}


const PartitionNode *
PartitionImpl::getLeaf(int idx) const {
  IaP<PartitionNodeImpl> ptr = representants[idx];

  while (ptr.isInteger())
    ptr = representants[ptr];

  return ptr;
}


const PartitionNode *
PartitionImpl::getSubTree(int idA, int idB) const {
  set<const PartitionNode *> candidates;

  const PartitionNode * curr = getLeaf(idA);
  while (curr) {
    candidates.insert(curr);
    curr = curr->getParent();
  }

  curr = getLeaf(idB);
  while (curr) {
    if (candidates.count(curr))
      break;
    else
      curr = curr->getParent();
  }

  return curr;
}


void
PartitionImpl::extractSchedule(vector<int>& schedule) const {
  unsigned int depth = 1;
  vector<int> visitIdx;
  PartitionNodeImpl * curr = root;

  schedule.clear();

  while (curr) {
    if (curr->id != -1) {
      schedule.push_back(curr->id);
      curr = curr->parent;
      depth--;
      continue;
    }

    if (depth > visitIdx.size())
      visitIdx.push_back(0);

    if (visitIdx.back() >= curr->getChildCount()) {
      curr = curr->parent;
      visitIdx.pop_back();
      depth--;
    } else {
      curr = curr->children[visitIdx[depth - 1]];
      visitIdx[depth - 1]++;
      depth++;
    }
  }
}
