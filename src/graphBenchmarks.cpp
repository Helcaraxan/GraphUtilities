// GraphUtilities library - Directed graph manipulation and querying
// Copyright (C) 2016 - Duco van Amstel
//
// For more license information please see the 'LICENSE' file at the top-level
// of the source code tree.

#include "graph-utilities/implementation/graphImpl.hpp"

using namespace std;


// Benchmarks

bool
GraphImpl::getBenchmarksEnabled(void) {
#if defined(ENABLE_BENCHMARKS)
  return true;
#else // ENABLE_BENCHMARKS
  return false;
#endif // ENABLE_BENCHMARKS
}


#if defined(ENABLE_BENCHMARKS)
long long
GraphImpl::getCyclesSpentIndexing(void) {
  return cyclesSpentIndexing;
}


long long
GraphImpl::getCyclesSpentQuerying(void) {
  return cyclesSpentQuerying;
}
#endif // ENABLE_BENCHMARKS


void
GraphImpl::printBenchmarks(ostream &os) {
#if defined(ENABLE_BENCHMARKS)
  os << "\n---\nBenchmarking\n";
  os << "Speed performances:\n";
  os << "Cycles spent indexing: " << cyclesSpentIndexing << "\n";
  os << "Cycles spent querying: " << cyclesSpentQuerying << "\n";
  os << "-> average of " << cyclesSpentQuerying / queryCount;
  os << " cycles per query\n\n";
#else // ENABLE_BENCHMARKS
  os << "WARNING: Benchmarking has not been enabled at compile time.\n";
  os << "WARNING: To use benchmarks run cmake with '-DENABLE_BENCHMARKS=1'.\n";
#endif // ENABLE_BENCHMARKS
}
