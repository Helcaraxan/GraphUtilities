#include <map>
#include <chrono>
#include <random>
#include <string>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <iostream>

#include <getopt.h>
#include <tbb/concurrent_hash_map.h>

#include "graph.hpp"

using namespace std;


static int dryFlag = 0;
static int uniqueFlag = 0;
static int verifyFlag = 0;
static int queryNumber = 100000;
static bool poison = false;
static semaphore openQueries;
static tbb::concurrent_hash_map<Query *, bool> queryMap;


static const struct option longopts[] = {
  {"input",     required_argument, 0, 'i'},
  {"output",    required_argument, 0, 'o'},
  {"test",      required_argument, 0, 't'},
  {"method",    required_argument, 0, 'm'},
  {"search",    required_argument, 0, 's'},
  {"count",     required_argument, 0, 'c'},
  {"graph",     optional_argument, 0, 'g'},
  {"queries",   optional_argument, 0, 'q'},
  {"help",      no_argument,       0, 'h'},
  {"verify",  no_argument, &verifyFlag, 1},
  {"unique-edges", no_argument, &uniqueFlag, 1},
  {"dry",       no_argument,  &dryFlag, 1},
  {"batch",     no_argument,  &batchFlag, 1},
  {0,0,0,0}
};


void
printHelpMessage() {
  cout << "Usage: -i|--input=<input-file> [options]\n";
  cout << "File options:\n";
  cout << "\t-o | --output=<output-file>\t\tFile to which query statistics and benchmarks are dumped (default: stdout)\n";
  cout << "\t-t | --test=<test-file>\t\t\tFile from which queries are read (default: random queries)\n";
  cout << "\t-g | --graph=<graph-dump-file>\t\tFile to which to dump the condensed graph with unique edges (default: no dump)\n";
  cout << "\t-q | --queries=<query-dump-file>\tFile to which randomly generated queries should be dumpes (default: no dump)\n";
  cout << "\nQuery options:\n";
  cout << "\t-m | --method=<index-method>\tMethod from <ShortIndexing|SuccessorOrder|Standard|LabelOrder> (default: Standard)\n";
  cout << "\t-s | --search=<search-method>\tMethod from <DFS|BBFS> (default: automatically select the fastest for the graph)\n";
  cout << "\t-c | --count=<number>\t\tNumber of queries to generate (only used when --test is not used)\n";
  cout << "\nBoolean options:\n";
  cout << "\t-h | --help\t\tDisplay this help message\n";
  cout << "\t-v | --verify\t\tVerify the query results by a DFS query that ignores labeling\n";
  cout << "\t-u | --unique-edges\tDon't check for double-edges on the input graph (speeds-up parsing of large graphs)\n";
  cout << "\t-d | --dry\t\tDo not perform queries stop after graph condensation (and eventual dumping)\n";
  cout << "\t-b | --batch\t\tDo not print progress-bar\n";
}


void
queryGenerator(Graph * graph, fstream &testFile, SearchMethod method) {
  int i, a, b, res;
  string garbage;
  Query * newQuery = NULL;
  ReachabilityQuery * subQuery = NULL;
  tbb::concurrent_hash_map<Query *, bool>::accessor queryAccess;

  // Fill the queries to process...
  if (testFile.is_open()) {
    // ... from the specified file

    // Count the number of queries
    queryNumber = 0;
    while (getline(testFile, garbage))
      queryNumber++;

    testFile.clear();
    testFile.seekg(0, testFile.beg);
    configureProgressBar(NULL, queryNumber);

    // Parse the actual queries
    while (testFile.good()) {
      testFile >> a >> b >> res;

      if (testFile.eof())
        break;

      subQuery =
        new ReachabilityQuery(graph->getVertexFromId(a), graph->getVertexFromId(b), method);
      newQuery = new Query(Reachability, subQuery);

      queryMap.insert(queryAccess, pair<Query *, bool>(newQuery, (res == 0 ? false : true)));
      queryAccess.release();
      graph->pushQuery(newQuery);
      openQueries.post();
    }

    testFile.close();
  } else {
    // ... at random
    default_random_engine generator(chrono::system_clock::now().time_since_epoch().count());
    uniform_int_distribution<int> queryDistribution(0, graph->getVertexCount() - 1);

    progressBarFinish = queryNumber;
    for (i = 0; i < queryNumber; i++) {
      a = queryDistribution(generator);
      do {
        b = queryDistribution(generator);
      } while (a == b);

      if (a < b)
        subQuery =
          new ReachabilityQuery(graph->getVertexFromId(a), graph->getVertexFromId(b), method);
      else
        subQuery =
          new ReachabilityQuery(graph->getVertexFromId(b), graph->getVertexFromId(a), method);

      newQuery = new Query(Reachability, subQuery);

      queryMap.insert(queryAccess, pair<Query *, bool>(newQuery, false));
      queryAccess.release();
      graph->pushQuery(newQuery);
      openQueries.post();
    }
  }

  poison = true;
  openQueries.post();
}


void
resultAnalysis(Graph * graph, fstream &queryFile) {
  int resultCounter = 0;
  string barTitle = "Parsing results ";
  Query * nextQuery = NULL;
  tbb::concurrent_hash_map<Query *, bool>::const_accessor queryAccess;

  configureProgressBar(&barTitle, 0);

  while (true) {
    openQueries.wait();

    if (poison && (resultCounter == queryNumber))
      break;

    nextQuery = graph->pullResult();

    if (nextQuery->type != Reachability) {
      cerr << "ERROR: Encountered an unexpected partition query." << endl;
      continue;
    }

    queryMap.find(queryAccess, nextQuery);

    if (nextQuery->query.reachability->getError()) {
      cerr << "ERROR: Could not process query " << nextQuery->query.reachability->getSource()->id;
      cerr << " -> " << nextQuery->query.reachability->getTarget()->id << "\n";
    } else if (verifyFlag) {
      if (nextQuery->query.reachability->getAnswer() != queryAccess->second) {
        cerr << "ERROR: Wrong answer for query " << nextQuery->query.reachability->getSource()->id;
        cerr << " -> " << nextQuery->query.reachability->getTarget()->id << "\n";
      }
    } else if (queryFile.is_open()) {
      queryFile << nextQuery->query.reachability->getSource()->id << " ";
      queryFile << nextQuery->query.reachability->getTarget()->id << " ";

      if (queryAccess->second)
        queryFile << "1\n";
      else
        queryFile << "0\n";
    }
    queryAccess.release();
    resultProgressBar(++resultCounter);
  }

  cout << "\n";

  if (queryFile.is_open())
    queryFile.close();
}


int
main(int argc, char * argv[]) {
	Graph * graph = NULL;
  IndexMethod indexMethod = Standard;
  SearchMethod searchMethod = UndefinedSearchMethod;
  fstream outputFile, testFile, dumpFile, queryFile;
  ostream * output = &cout;
  int i, c;
  char fileName[512] = {'\0'};

  // Parse command-line options
  while ((c = getopt_long(argc, argv, "i:o:t:m:s:c:g::q::hvudb", longopts, NULL)) != -1) {
    switch (c) {
      case 'i':
        strncpy(fileName, optarg, 511);
        break;

      case 'o':
        outputFile.open(optarg, fstream::out);
        if (!outputFile.good()) {
          cerr << "ERROR: Could not open the output file.\n";
          exit(EXIT_FAILURE);
        }

        output = &outputFile;
        break;

      case 't':
        testFile.open(optarg, fstream::in);
        if (!testFile.good()) {
          cerr << "ERROR: Could not open the test file.\n";
          exit(EXIT_FAILURE);
        }
        break;

      case 'm':
        if (!strcmp(optarg, "ShortIndexing")) {
          indexMethod = ShortIndexing;
        } else if (!strcmp(optarg, "SuccessorOrder")) {
          indexMethod = SuccessorOrder;
        } else if (!strcmp(optarg, "Standard")) {
          indexMethod = Standard;
        } else if (!strcmp(optarg, "LabelOrder")) {
          indexMethod = LabelOrder;
        } else {
          cerr << "ERROR: Unknown label method specified.\n";
          printHelpMessage();
          exit(EXIT_FAILURE);
        }
        break;

      case 's':
        if (!strcmp(optarg, "DFS")) {
          searchMethod = DFS;
        } else if (!strcmp(optarg, "BBFS")) {
          searchMethod = BBFS;
        } else {
          cerr << "ERROR: Unknown search method specified.\n";
          printHelpMessage();
          exit(EXIT_FAILURE);
        }
        break;

      case 'c':
        if (!isdigit(optarg[0])) {
          cerr << "ERROR: The -c | --count argument is not a number\n";
          printHelpMessage();
          exit(EXIT_FAILURE);
        } else {
          queryNumber = atoi(optarg);
        }
        break;

      case 'g':
        dumpFile.open(optarg, fstream::out);
        if (!dumpFile.good()) {
          cerr << "ERROR: Could not open the graph dump file.\n";
          exit(EXIT_FAILURE);
        }
        break;

      case 'q':
        queryFile.open(optarg, fstream::out);
        if (!queryFile.good()) {
          cerr << "ERROR: Could not open the query dump file.\n";
          exit(EXIT_FAILURE);
        }
        break;

      case 'h':
        printHelpMessage();
        exit(EXIT_SUCCESS);
        break;

      case 0:
        break;

      default:
        cerr << "ERROR: GetOpt class error while parsing command-line arguments.\n";
        exit(EXIT_FAILURE);
    }
  }

  if (fileName[0] == '\0') {
    cerr << "ERROR: No input file was specified.\n\n";
    printHelpMessage();
    exit(EXIT_FAILURE);
  }

  if (strstr(fileName, ".dot")) {
    graph = Graph::createFromDotFile(fileName, uniqueFlag);
  } else if (strstr(fileName, ".gra")) {
    graph = Graph::createFromGraFile(fileName, uniqueFlag);
  } else {
    cerr << "ERROR: Unknown graph-file extension. Only accepts .dot and .gra files.\n";
    exit(EXIT_FAILURE);
  }

  if (!graph) {
    cerr << "ERROR: No graph was generated from the input file.\n\n";
    printHelpMessage();
    exit(EXIT_FAILURE);
  }
  
  // Set the indexing method as specified / default
  graph->setIndexMethod(indexMethod);

  // Dump the parsed graph for verification purposes
  if (dumpFile.is_open()) {
    Vertex * curr;

    dumpFile << "graph_for_greach\n" << graph->getVertexCount() << "\n";
    for (i = 0; i < (int) graph->getVertexCount(); i++) {
      dumpFile << i << ":";
      curr = graph->getVertexFromId(i);
      for (auto it = curr->successors_begin(); it != curr->successors_end(); ++it)
        dumpFile << " " << (*it)->id;

      dumpFile << " #\n";
    }
    dumpFile.close();
  }

  // Exit here in case of a 'dry' run
  if (dryFlag)
    exit(EXIT_SUCCESS);

  // Process the queries with two threads
  thread pushThread(queryGenerator, graph, ref(testFile), searchMethod);
  thread pullThread(resultAnalysis, graph, ref(queryFile));

  // Wait for the queries to be issued and treated
  pushThread.join();
  pullThread.join();
  graph->endOfQueries();
	
  // Dump the statistics of the queries...
  if (graph->statisticsAreEnabled())
    graph->printStatistics(*output);

  // ...and the benchmarks
  if (graph->benchmarksAreEnabled())
    graph->printBenchmarks(*output);

  // Clear up the queries
  for (auto it = queryMap.begin(), end = queryMap.end(); it != end; ++it)
    delete it->first;

  queryMap.clear();

  delete graph;

	//exit(EXIT_SUCCESS);
  return 0;
}
