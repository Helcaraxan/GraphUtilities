#include <set>
#include <string>
#include <climits>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>

#include <getopt.h>

#include <graph-utilities/graph.hpp>

using namespace std;


// Command-line option specifications

static const struct option longopts[] = {
  {"graph",      required_argument, 0, 'g'},
  {"schedule",   required_argument, 0, 's'},
  {"tile",       required_argument, 0, 'T'},
  {"evaluation", required_argument, 0, 'e'},
  {"method",     required_argument, 0, 'm'},
  {"memory",     required_argument, 0, 'M'},
  {"threads",    required_argument, 0, 't'},
  {"help",       no_argument,       0, 'h'},
  {0,0,0,0}
};


// Help message printing

void
printHelpMessage() {
  cout << "Usage: -g|--graph=<graph> -m|--method=<method> [options]\n";
  cout << "File options:\n";
  cout << "\t-s | --schedule=<file>\t\tFile to which the schedule produced by the partition will be dumped.\n";
  cout << "\t-T | --tile=<file>\t\tFile to which to dump the tiles created for IO complexity evaluation.\n";
  cout << "Partition options:\n";
  cout << "\t-m | --method=<method>\t\tMethod from <Convexify|MaxDistance>\n";
  cout << "\t-e | --evaluation=<type>\tMethod from <TotalLoads|AvgLoadStore>\n";
  cout << "\t-M | --memory=<size>\t\tSize of the memory for which IO complexity should be computed.\n";
  cout << "\t-t | --threads=<count>\t\tNumber of worker threads to use for the partitioning (only for MaxDistance).\n";
  cout << "\nMiscellaneous options:\n";
  cout << "\t-h | --help\t\tDisplay this help message\n";
}


// Driver function

int
main(int argc, char * argv[]) {
  int c;
  int threadCount = 1;
  int memorySize = 256;
  string graphFile = "", schedFile = "", tileFile = "";
  Graph * graph = nullptr;
  PartitionQuery * pQuery = nullptr;
  PartitionMethod method = UndefinedPartitionMethod;
  IOComplexityType type = UndefinedIOType;


  // Parse command-line options
  while ((c = getopt_long(argc, argv, "g:s:T:e:m:M:t:h", longopts, nullptr)) != -1) {
    switch (c) {
      case 'g':
        graphFile = optarg;
        break;

      case 's':
        schedFile = optarg;
        break;

      case 'T':
        tileFile = optarg;
        break;

      case 'e':
        if (!strcmp(optarg, "TotalLoads")) {
          type = TotalLoads;
        } else if (!strcmp(optarg, "AvgLoadStore")) {
          type = AvgLoadStore;
        } else {
          cerr << "ERROR: Unknown argument to --evaluation." << endl;
          printHelpMessage();
          exit(EXIT_FAILURE);
        }
        break;

      case 'm':
        if (!strcmp(optarg, "Convexify")) {
          method = Convexify;
        } else if (!strcmp(optarg, "MaxDistance")) {
          method = MaxDistance;
        } else {
          cerr << "ERROR: Unknown argument to --method." << endl;
          printHelpMessage();
          exit(EXIT_FAILURE);
        }
        break;

      case 'M':
        if (!isdigit(optarg[0])) {
          cerr << "ERROR: The -M | --memory argument is not a number\n";
          printHelpMessage();
          exit(EXIT_FAILURE);
        } else {
          memorySize = atoi(optarg);
        }
        break;

      case 't':
        if (!isdigit(optarg[0])) {
          cerr << "ERROR: The -t | --threads argument is not a number\n";
          printHelpMessage();
          exit(EXIT_FAILURE);
        } else {
          threadCount = atoi(optarg);
        }
        break;

      case 'h':
        printHelpMessage();
        exit(EXIT_SUCCESS);
        break;

      case 0:
        break;

      default:
        cerr << "ERROR: GetOpt class error while parsing command-line "
          << " arguments." << endl;
        exit(EXIT_FAILURE);
    }
  }

  if (graphFile.size() == 0) {
    cerr << "ERROR: No graph file was specified.\n\n";
    printHelpMessage();
    exit(EXIT_FAILURE);
  }

  if (!graphFile.substr(graphFile.size() - 4, 4).compare(".dot")) {
    graph = parseDotFile(graphFile.c_str(), true);
  } else if (!graphFile.substr(graphFile.size() - 4, 4).compare(".gra")) {
    graph = parseGraFile(graphFile.c_str(), true);
  } else {
    cerr << "ERROR: Unknown graph-file extension. Only accepts .dot and .gra "
      << "files.\n";
    exit(EXIT_FAILURE);
  }

  if (!graph) {
    cerr << "ERROR: No graph was generated from the graph file.\n\n";
    printHelpMessage();
    exit(EXIT_FAILURE);
  }


  // Set the indexing method in case it is needed and start worker threads
  graph->setIndexMethod(Standard);
  graph->enableQueries();


  // Create a PartitionQuery instance and set the partition method & threadcount
  pQuery = createPQuery();
  pQuery->setMethod(method);
  pQuery->setThreadCount(threadCount);
  

  // Submit the query and retrieve it once completed 
  graph->pushQuery(pQuery);
  graph->pullQuery(true);


  // Retrieve the obtained Partition
  const Partition * part = pQuery->getPartition();

  if (!part) {
    cerr << "ERROR: No partition was found.\n";
    exit(EXIT_FAILURE);
  }


  // Compute statistics over the Partition
  int minLevel = INT_MAX;
  int maxLevel = 0;
  double acc = 0;

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

  cout << "\nPartition statistics:\n";
  cout << "Average depth: " << acc / (double) graph->getVertexCount() << "\n";
  cout << "Minimum depth: " << minLevel << "\n";
  cout << "Maximum depth: " << maxLevel << "\n";
  cout.flush();


  // Retrieve the scheduling that is implied by the Partition instance if
  // necessary
  vector<int> schedule;
  if ((type != UndefinedIOType) || (schedFile.size() > 0)) {
    part->extractSchedule(schedule);

    // Check the schedule for validity
    cout << "\nVerifying scheduling... ";
    if (graph->checkSchedule(schedule))
      cout << "VALID" << endl;
    else
      cout << "INVALID" << endl;
  }

  // Perform IO complexity evaluation if required
  if (type != UndefinedIOType) {
    double cost = -1;

    cout << "\nEvaluating schedule costs:\n";
    if (tileFile.size() > 0)
      cost = graph->getPartitionCost(part, memorySize, type, tileFile.c_str());
    else
      cost = graph->getPartitionCost(part, memorySize, type);

    cout << "Target memory size: " << memorySize << " 32-bit words\n";
    cout << "IO complexity: " << cost << endl;
  }

  // Dump the schedule when requested
  if (schedFile.size() > 0) {
    fstream schedStream(schedFile.c_str(), ios_base::out);

    if (!schedStream.good()) {
      cerr << "\nERROR: Could not open the specified schedule dump file.\n";
      exit(EXIT_FAILURE);
    }

    cout << "\nDumping schedule...";
    cout.flush();

    for (auto it = schedule.begin(), end = schedule.end(); it != end; ++it)
      schedStream << *it << "\n";

    schedStream.close();

    cout << " DONE" << endl;
  }

  delete part;
  delete pQuery;
  delete graph;

  exit(EXIT_SUCCESS);
}
