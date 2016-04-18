#include <set>
#include <string>
#include <climits>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <algorithm>

#include <getopt.h>

#include <graph-utilities/graph.hpp>

using namespace std;


// Externally defined variables for progress bars
extern int batchFlag;


// Command-line option specifications

static const struct option longopts[] = {
  {"graph",           required_argument,  0, 'g'},
  {"schedule",        required_argument,  0, 's'},
  {"tile",            required_argument,  0, 'T'},
  {"hint",            required_argument,  0, 'H'},
  {"evaluation",      required_argument,  0, 'e'},
  {"method",          required_argument,  0, 'm'},
  {"memory",          required_argument,  0, 'M'},
  {"threads",         required_argument,  0, 't'},
  {"partition-count", required_argument,  0, 'p'},
  {"original",        no_argument,        0, 'o'},
  {"batch",           no_argument, &batchFlag, 1},
  {"help",            no_argument,        0, 'h'},
  {0,0,0,0}
};


// Help message printing

void
printHelpMessage() {
  cout << "Usage: -g|--graph=<graph> -m|--method=<method> [options]\n";
  cout << "File options:\n";
  cout << "\t-s | --schedule=<file>\t\tFile to which the schedule produced by the partition will be dumped.\n";
  cout << "\t-T | --tile=<file>\t\tFile to which to dump the tiles created for IO complexity evaluation.\n";
  cout << "\t-H | --hint=<file>\t\tFile with a proposed scheduling of nodes for IO evaluation (works only with -e | --evaluation).\n";
  cout << "Partition options:\n";
  cout << "\t-p | --partition-count=<int>\tNumber of partitions to create with the PaToH library.\n";
  cout << "\t-m | --method=<method>\t\tMethod from <Convexify|MaxDistance|Greedy>\n";
  cout << "\t-e | --evaluation=<type>\tMethod from <TotalLoads|AvgLoadStore>\n";
  cout << "\t-M | --memory=<size>\t\tSize of the memory for which IO complexity should be computed.\n";
  cout << "\t-t | --threads=<count>\t\tNumber of worker threads to use for the partitioning (only for MaxDistance).\n";
  cout << "\nMiscellaneous options:\n";
  cout << "\t-o | --original\t\tEvaluate the original scheduling's IO complexity (works only with -e | --evaluation).\n";
  cout << "\t-b | --batch\t\tDo not print progress-bar\n";
  cout << "\t-h | --help\t\tDisplay this help message\n";
}


// Driver function

int
main(int argc, char * argv[]) {
  int c;
  int threadCount = 1;
  int partitionCount = 0;
  bool evaluateOriginal = false;
  set<int> memorySizes;
  string graphFile = "", schedFile = "", tileFile = "", hintFile = "";
  Graph * graph = nullptr;
  PartitionQuery * pQuery = nullptr;
  PartitionMethod method = UndefinedPartitionMethod;
  IOComplexityType type = UndefinedIOType;


  // Parse command-line options
  while ((c = getopt_long(argc, argv, "g:s:T:H:e:m:M:t:p:obh", longopts, nullptr)) != -1) {
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

      case 'H':
        hintFile = optarg;
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
        } else if (!strcmp(optarg, "PaToH")) {
          method = PaToH;
        } else if (!strcmp(optarg, "Greedy")) {
          method = Greedy;
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
          memorySizes.insert(atoi(optarg));
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

      case 'p':
        if (!isdigit(optarg[0])) {
          cerr << "ERROR: The -p | --partition-count argument is not a number\n";
          printHelpMessage();
          exit(EXIT_FAILURE);
        } else {
          partitionCount = atoi(optarg);
        }
        break;

      case 'o':
        evaluateOriginal = true;
        break;

      case 'b':
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

  if ((method == Greedy) && (memorySizes.size() > 1)) {
    cerr << "ERROR: Greedy method is not compatible with multiple memory "
      << "sizes.\n\n";
    exit(EXIT_FAILURE);
  }

  if (memorySizes.empty())
    memorySizes.insert(256);

  // Set the indexing method in case it is needed and start worker threads
  graph->setIndexMethod(Standard);
  graph->enableQueries();


  // Create a PartitionQuery instance and set the partition method & threadcount
  pQuery = createPQuery();
  pQuery->setMethod(method);
  pQuery->setThreadCount(threadCount);
  if (method == Greedy)
    pQuery->setMemorySize(*memorySizes.begin());
  

  // Submit the query and retrieve it once completed 
  graph->pushQuery(pQuery);
  graph->pullQuery(true);


  // Retrieve the obtained Partition
  const Partition * part = pQuery->getPartition();

  if (!part) {
    cerr << "ERROR: No partition was found.\n";
    exit(EXIT_FAILURE);
  }


  if (method != Greedy) {
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
    cout << "Maximum depth: " << maxLevel << "\n" << endl;
    cout.flush();
  }


  // Retrieve the scheduling implied by the Partition instance if necessary
  vector<int> schedule;
  if (((type != UndefinedIOType) || (schedFile.size() > 0))
      && (method != Greedy)) {
    part->extractSchedule(schedule);

    if (method != PaToH) {
      // Check the schedule for validity except for PaToH partitioning
      cout << "Verifying scheduling implied by partition... ";
      if (graph->checkSchedule(schedule)) {
        cout << "VALID\n";
      } else {
        cout << "INVALID" << endl;
        exit(EXIT_FAILURE);
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
    } else if (schedFile.size() > 0) {
      cerr << "WARNING: Can not dump schedule for PaToH partitioning. It would "
        << "not be valid.\n";
      exit(EXIT_FAILURE);
    }
  }

  // Perform IO complexity evaluation if required
  if (type != UndefinedIOType) {
    const Partition * hintPart = nullptr;
    const Partition * origPart = nullptr;

    if ((tileFile.size() > 0) && memorySizes.size() > 1) {
      cerr << "ERROR: Can not dump tiles to file with multiple specified memory "
        << " sizes.\n";
      exit(EXIT_FAILURE);
    }

    if (hintFile.size() > 0) {
      fstream hintStream(hintFile.c_str(), ios_base::in);

      for (int i = 0, e = schedule.size(); i < e; i++) {
        if (hintStream.eof()) {
          cerr << "ERROR: Hinted schedule does not contain all vertex IDs.\n";
          exit(EXIT_FAILURE);
        }

        hintStream >> schedule[i];
      }

      cout << "Verifying scheduling implied by hinted partition... ";
      if (graph->checkSchedule(schedule)) {
        cout << "VALID\n";
      } else {
        cout << "INVALID" << endl;
        exit(EXIT_FAILURE);
      }

      hintPart = graph->computeConvexPartition(schedule);
    }

    if (evaluateOriginal) {
      iota(schedule.begin(), schedule.end(), 0);

      cout << "Verifying the implicit original schedule... ";
      if (graph->checkSchedule(schedule)) {
        cout << "VALID\n";
      } else {
        cout << "INVALID" << endl;
        exit(EXIT_FAILURE);
      }

      origPart = graph->computeConvexPartition(schedule);
    }

    cout << "Evaluating partition costs:\n";
    for (const int &size : memorySizes) {
      double cost = -1;

      cout << "- Target memory size: " << size << " 32-bit words" << endl;

      if (method != Greedy) {
        if (tileFile.size() > 0)
          cost = graph->getPartitionCost(part, size, type, tileFile.c_str());
        else
          cost = graph->getPartitionCost(part, size, type);
      } else {
        cost = graph->getCutCost(part);
      }

      cout << "  Partition IO complexity: " << cost << endl;

      if (hintFile.size() > 0) {
        cost = graph->getPartitionCost(hintPart, size, type);
        cout << "  Hinted schedule IO complexity: " << cost << "\n";
      }

      if (evaluateOriginal) {
        cost = graph->getPartitionCost(origPart, size, type);
        cout << "  Original schedule IO complexity: " << cost << "\n";
      }

      cout << endl;
    }
  }

  // When required perform a comparison on k-way partitioning
  if (partitionCount) {
    int cutCost = 0;
    cout << "Using PaToH for a " << partitionCount << "-way partitioning." << endl;

    cutCost = graph->getCutCost(partitionCount);
    cout << "PaToH cut-cost: " << cutCost << "\n\n";

    cout << "Perform a " << partitionCount << "-way partitioning on the "
      << "computed partition." << endl;

    cutCost = graph->getCutCost(part, partitionCount);
    cout << "Hierarchical cut-cost: " << cutCost << endl;
  }

  // Final clean-up
  delete part;
  delete pQuery;
  delete graph;

  exit(EXIT_SUCCESS);
}
