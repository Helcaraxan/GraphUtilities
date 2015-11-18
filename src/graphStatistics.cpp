#include "graph-utilities/implementation/graphImpl.hpp"

using namespace std;


// Statistics

bool
GraphImpl::getStatisticsEnabled() {
#if defined(ENABLE_STATISTICS)
  return true;
#else // ENABLE_STATISTICS
  return false;
#endif // ENABLE_STATISTICS
}


#if defined(ENABLE_STATISTICS)
uintmax_t
GraphImpl::getQueryCount() {
  return queryCount;
}


uintmax_t
GraphImpl::getPositiveQueryCount() {
  return positiveQueryCount;
}


uintmax_t
GraphImpl::getNegativeQueryCount() {
  return negativeQueryCount;
}


uintmax_t
GraphImpl::getShortNegativeQueryCount() {
  return shortNegativeQueryCount;
}


double
GraphImpl::getPositiveQueryOverhead() {
  return positiveQueryOverhead;
}


double
GraphImpl::getNegativeQueryOverhead() {
  return negativeQueryOverhead;
}
#endif // ENABLE_STATISTICS


void
GraphImpl::printStatistics(ostream &os) {
#if defined(ENABLE_STATISTICS)
  double shortFraction =
    ((double) shortNegativeQueryCount / ((double) negativeQueryCount));
  os << "\n---\nStatistics:\n";
  os << "General statistics\n";
  os << "Number of vertices: " << getVertexCount() << "\n";
  os << "Number of edges: " << edgeCount << "\n";
  os << "Number of performed queries: " << queryCount << "\n";

  os << "\nPositive query statistics:\n";
  os << "- Number of positive answers: " << positiveQueryCount << "\n";
  if (positiveQueryCount && (getMethod() == DFS)) {
    os.precision(3);
    os << "- Average DFS length / path-length ratio: " << positiveQueryOverhead;
    os << "\n";
  }

  os << "\nNegative query statistics:\n";
  os << "- Number of negative answers   : " << negativeQueryCount << "\n";
  if (negativeQueryCount) {
    os << "- Number of immediate negatives: " << shortNegativeQueryCount;
    os.precision(4);
    os << " (" << shortFraction * 100 << "%).\n";
    if (getMethod() == DFS) {
      os.precision(3);
      os << "- Average DFS length / graph-size ratio: ";
      os << negativeQueryOverhead << "\n";
    }
  }
  os << "---\n\n";
#else // ENABLE_STATISTICS
  os << "WARNING: Statistics gathering has not been enabled at compile time.\n";
  os << "WARNING: To use statistics run cmake with '-DENABLE_STATISTICS=1'.\n";
#endif // ENABLE_STATISTICS
}


// Internal statistics maintenance

#if defined(ENABLE_STATISTICS)
void
GraphImpl::registerQueryStatistics(ReachabilityQuery * query) {
  double coefficient = 1.0;
  double overhead = 1.0;
  unique_lock<mutex> counterLock(counterMutex);
  ReachabilityQueryImpl * cast = dynamic_cast<ReachabilityQueryImpl *>(query);

  if (!cast)
    return;

  if (cast->getAnswer()) {
    positiveQueryCount++;

    coefficient = 1.0 / ((double) positiveQueryCount);
    if (cast->getMethod() == DFS)
      overhead =
        ((double) cast->getSearchedNodes()) / (double) cast->getPath().size();

    positiveQueryOverhead *= (1.0 - coefficient);
    positiveQueryOverhead += coefficient * overhead;
  } else {
    negativeQueryCount++;

    if (cast->getSearchedNodes() > 0) {
      coefficient = 1.0 / ((double) negativeQueryCount);
      overhead =
        ((double) cast->getSearchedNodes()) / (double) getVertexCount();

      negativeQueryOverhead *= (1.0 - coefficient);
      negativeQueryOverhead += coefficient * overhead;
    } else {
      shortNegativeQueryCount++;
    }
  }
}
#endif // ENABLE_STATISTICS
