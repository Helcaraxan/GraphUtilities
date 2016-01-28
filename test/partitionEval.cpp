#include <climits>
#include <cstdlib>
#include <iostream>
#include <set>

#include <graph-utilities/graph.hpp>


int
main(int argc, char * argv[]) {
  int minLevel, maxLevel;
  double acc = 0;
  vector<int> schedule;
  Graph * graph = nullptr;
  PartitionQuery * pQuery = nullptr;
  IOComplexityType type = UndefinedIOType;


  if (string(argv[2]) == "GRA")
    graph = parseGraFile(argv[1], true);
  else
    graph = parseDotFile(argv[1], true);

  if (string(argv[3]) == "ALS")
    type = AvgLoadStore;
  else
    type = TotalLoads;

  graph->setIndexMethod(Standard);

  if (!graph->isDAG()) {
    cerr << "Error: is not a DAG\n";
    exit(EXIT_FAILURE);
  }

  graph->enableQueries();


  cout << "Performing Convexify partitioning\n";
  pQuery = createPQuery();
  pQuery->setMethod(Convexify);

  graph->pushQuery(pQuery);
  graph->pullQuery(true);

  cout << "Partitioning done\n";
  const Partition * part = pQuery->getPartition();

  if (!part) {
    cerr << "ERROR: No partition was found.\n";
    exit(EXIT_FAILURE);
  }

  minLevel = INT_MAX;
  maxLevel = 0;
  acc = 0;

  cout << "Computing partition statistics.\n";
  for (int i = 0, e = graph->getVertexCount(); i < e; i++) {
    int lCount = 0;
    const PartitionNode * node = part->getLeaf(i);
    
    while (node) {
      lCount++;
      node = node->getParent();
    }

    acc += lCount;

    if (lCount > maxLevel)
      maxLevel = lCount;

    if (lCount < minLevel)
      minLevel = lCount;
  }
  cout << "Average depth: " << acc / (double) graph->getVertexCount() << "\n";
  cout << "Minimum depth: " << minLevel << "\n";
  cout << "Maximum depth: " << maxLevel << "\n";

  cout << "Verifying scheduling... ";
  cout.flush();

  part->extractSchedule(schedule);
  if (graph->checkSchedule(schedule))
    cout << "VALID" << endl;
  else
    cout << "INVALID" << endl;


  cout << "Evaluating schedule costs:\n";
  cout << "Target memory (32)  : ";
  cout.flush();

  cout << graph->getPartitionCost(part, 32, type) << endl;

  cout << "Target memory (1k)  : ";
  cout.flush();

  cout << graph->getPartitionCost(part, 1024, type) << endl;

  cout << "Target memory (64k) : ";
  cout.flush();

  cout << graph->getPartitionCost(part, 65536, type) << endl;

  cout << "Target memory (512k): ";
  cout.flush();

  cout << graph->getPartitionCost(part, 524288, type) << endl;

  cout << "Target memory (2M)  : ";
  cout.flush();

  cout << graph->getPartitionCost(part, 2097152, type) << endl;

  delete part;
  delete pQuery;


  cout << "Performing Max-Distance partitioning\n";
  pQuery = createPQuery();
  pQuery->setMethod(MaxDistance);

  graph->pushQuery(pQuery);
  graph->pullQuery(true);

  cout << "Partitioning done\n";
  part = pQuery->getPartition();

  if (!part) {
    cerr << "ERROR: No partition was found.\n";
    exit(EXIT_FAILURE);
  }

  minLevel = INT_MAX;
  maxLevel = 0;
  acc = 0;

  cout << "Computing partition statistics.\n";
  for (int i = 0, e = graph->getVertexCount(); i < e; i++) {
    int lCount = 0;
    const PartitionNode * node = part->getLeaf(i);
    
    while (node) {
      lCount++;
      node = node->getParent();
    }

    acc += lCount;

    if (lCount > maxLevel)
      maxLevel = lCount;

    if (lCount < minLevel)
      minLevel = lCount;
  }
  cout << "Average depth: " << acc / (double) graph->getVertexCount() << "\n";
  cout << "Minimum depth: " << minLevel << "\n";
  cout << "Maximum depth: " << maxLevel << "\n";

  cout << "Verifying scheduling... ";
  cout.flush();

  part->extractSchedule(schedule);
  if (graph->checkSchedule(schedule))
    cout << "VALID" << endl;
  else
    cout << "INVALID" << endl;


  cout << "Evaluating schedule costs:\n";
  cout << "Target memory (32)  : ";
  cout.flush();

  cout << graph->getPartitionCost(part, 32, type) << endl;

  cout << "Target memory (1k)  : ";
  cout.flush();

  cout << graph->getPartitionCost(part, 1024, type) << endl;

  cout << "Target memory (64k) : ";
  cout.flush();

  cout << graph->getPartitionCost(part, 65536, type) << endl;

  cout << "Target memory (512k): ";
  cout.flush();

  cout << graph->getPartitionCost(part, 524288, type) << endl;

  cout << "Target memory (2M)  : ";
  cout.flush();

  cout << graph->getPartitionCost(part, 2097152, type) << endl;

  delete part;
  delete pQuery;


  exit(EXIT_SUCCESS);
}
