#include <map>
#include <atomic>
#include <chrono>
#include <random>
#include <string>
#include <thread>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <iostream>

#include <getopt.h>
#include <tbb/concurrent_hash_map.h>

#include "graph-utilities/graph.hpp"
#include "graph-utilities/implementation/support.hpp"
#include "graph-utilities/implementation/semaphore.hpp"

using namespace std;


// Externally defined variables for progress bars

extern int batchFlag;
extern int progressBarFinish;


// Local variables

static int dryFlag = 0;
static int uniqueFlag = 0;
static int verifyFlag = 0;
static int queryNumber = 100000;
static tbb::concurrent_hash_map<Query *, bool> queryMap;


// Command-line option specifications

static const struct option longopts[] = {
  {"input",     required_argument, 0, 'i'},
  {"output",    required_argument, 0, 'o'},
  {"test",      required_argument, 0, 't'},
  {"graph",     required_argument, 0, 'g'},
  {"queries",   required_argument, 0, 'q'},
  {"index-method",   required_argument, 0, 'I'},
  {"search-method",  required_argument, 0, 'S'},
  {"count",     required_argument, 0, 'c'},
  {"coarsen-factor", required_argument, 0, 'f'},
  {"coarsen-method", required_argument, 0, 'C'},
  {"coarsen-file",   required_argument, 0, 'F'},
  {"help",      no_argument,       0, 'h'},
  {"verify",  no_argument,      &verifyFlag, 1},
  {"unique-edges", no_argument, &uniqueFlag, 1},
  {"dry",       no_argument,    &dryFlag,    1},
  {"batch",     no_argument,    &batchFlag,  1},
  {0,0,0,0}
};


// Driver functions

void
printHelpMessage() {
  cout << "Usage: -i|--input=<input-file> [options]\n";
  cout << "File options:\n";
  cout << "\t-o | --output=<output-file>\t\tFile to which query statistics and benchmarks are dumped\n";
  cerr << "\t\t\t\t\t\t(default: stdout)\n";
  cout << "\t-t | --test=<test-file>\t\t\tFile from which queries are read (default: random queries)\n";
  cout << "\t-g | --graph=<graph-dump-file>\t\tFile to which to dump the condensed graph with unique edges\n";
  cout << "\t\t\t\t\t\t(default: no dump)\n";
  cout << "\t-q | --queries=<query-dump-file>\tFile to which randomly generated queries should be dumped\n";
  cout << "\t\t\t\t\t\t(default: no dump)\n";
  cout << "\nQuery options:\n";
  cout << "\t-I | --index-method=<index-method>\tMethod from <ShortIndexing|SuccessorOrder|Standard|LabelOrder>\n";
  cout << "\t\t\t\t\t\t(default: Standard)\n";
  cout << "\t-S | --search-method=<search-method>\tMethod from <DFS|BBFS> (default: select fastest for the graph)\n";
  cout << "\t-c | --count=<number>\t\t\tNumber of queries to generate (only used when --test is not used)\n";
  cout << "\nCoarsening options:\n";
  cout << "\t-f | --coarsen-factor=<number>\t\tSpecify the coarsening factor that is to be applied\n";
  cout << "\t-C | --coarsen-method=<coarsen-method>\tMethod from <Greedy> (default: Greedy)\n";
  cout << "\t-F | --coarsen-file=<coarsen-dump-file>\tFile-prefix for dumping the coarsened graph and mapping\n";
  cout << "\nBoolean options:\n";
  cout << "\t-h | --help\t\tDisplay this help message\n";
  cout << "\t-v | --verify\t\tVerify the query results by a DFS query that ignores labeling\n";
  cout << "\t-u | --unique-edges\tNo check for double-edges on input graph (speeds-up parsing of large graphs)\n";
  cout << "\t-d | --dry\t\tDo not perform queries stop after graph condensation (and eventual dumping)\n";
  cout << "\t-b | --batch\t\tDo not print progress-bar\n";
}


void
queryGenerator(Graph * graph, const char * testFile, SearchMethod method) {
  int i, a, b, res;
  string garbage;
  fstream testStream;
  ReachabilityQuery * query = NULL;
  tbb::concurrent_hash_map<Query *, bool>::accessor queryAccess;

  if (queryNumber == 0)
    return;

  if (*testFile != '\0')
    testStream.open(testFile, ios_base::in);

  // Fill the queries to process...
  if (testStream.is_open()) {
    // ... from the specified file

    // Count the number of queries
    queryNumber = 0;
    while (getline(testStream, garbage))
      queryNumber++;

    testStream.clear();
    testStream.seekg(0, testStream.beg);
    configureProgressBar(NULL, queryNumber);

    // Parse the actual queries
    while (testStream.good()) {
      testStream >> a >> b >> res;

      if (testStream.eof())
        break;

      query = createRQuery();
      query->setSource(graph->getVertex(a));
      query->setTarget(graph->getVertex(b));
      query->setMethod(method);

      queryMap.insert(queryAccess,
          pair<Query *, bool>(query, (res == 0 ? false : true)));
      queryAccess.release();
      graph->pushQuery(query);
    }

    testStream.close();
  } else {
    // ... at random
    default_random_engine
      generator(chrono::system_clock::now().time_since_epoch().count());
    uniform_int_distribution<int> queryDistribution(0, graph->getVertexCount() - 1);

    progressBarFinish = queryNumber;
    for (i = 0; i < queryNumber; i++) {
      a = queryDistribution(generator);
      do {
        b = queryDistribution(generator);
      } while (a == b);

      query = createRQuery();
      query->setMethod(method);
      if (a < b) {
        query->setSource(graph->getVertex(a));
        query->setTarget(graph->getVertex(b));
      } else {
        query->setSource(graph->getVertex(b));
        query->setTarget(graph->getVertex(a));
      }

      queryMap.insert(queryAccess, pair<Query *, bool>(query, false));
      queryAccess.release();
      graph->pushQuery(query);
    }
  }
}


void
resultAnalysis(Graph * graph, const char * queryFile) {
  int resultCounter = 0;
  string barTitle = "Parsing results ";
  fstream queryStream;
  Query * query = NULL;
  ReachabilityQuery * rQuery = NULL;
  tbb::concurrent_hash_map<Query *, bool>::const_accessor queryAccess;

  if (queryNumber == 0)
    return;

  configureProgressBar(&barTitle, 0);

  if (*queryFile != '\0')
    queryStream.open(queryFile, ios_base::out);

  while (true) {
    query = graph->pullQuery(true);

    if (!query)
      break;

    rQuery = dynamic_cast<ReachabilityQuery *>(query);
    if (!rQuery) {
      cerr << "ERROR: Encountered an unexpected partition query." << endl;
      continue;
    }

    queryMap.find(queryAccess, query);

    if (rQuery->getError()) {
      cerr << "ERROR: Could not process query " << rQuery->getSource()->getId();
      cerr << " -> " << rQuery->getTarget()->getId() << "\n";
    } else if (verifyFlag) {
      if (rQuery->getAnswer() != queryAccess->second) {
        cerr << "ERROR: Wrong answer on query " << rQuery->getSource()->getId();
        cerr << " -> " << rQuery->getTarget()->getId() << "\n";
      }
    } else if (queryStream.is_open()) {
      queryStream << rQuery->getSource()->getId() << " ";
      queryStream << rQuery->getTarget()->getId() << " ";

      if (queryAccess->second)
        queryStream << "1\n";
      else
        queryStream << "0\n";
    }

    queryAccess.release();
    resultProgressBar(++resultCounter);
  }

  cout << "\n";

  if (queryStream.is_open())
    queryStream.close();
}


int
main(int argc, char * argv[]) {
  string inputFile, outputFile, coarseFile, testFile, dumpFile, queryFile;
	Graph * graph = NULL;
  CoarsenMethod coarsenMethod = Greedy;
  IndexMethod indexMethod = Standard;
  SearchMethod searchMethod = UndefinedSearchMethod;
  int c;
  int coarsenFactor = 0;
  int coarsenSecondaryFactor = 2;

  // Parse command-line options
  while ((c = getopt_long(argc, argv, "i:o:t:g:q:I:S:c:f:C:F:hvudb", longopts,
          NULL)) != -1) {
    switch (c) {
      case 'i':
        inputFile = optarg;
        break;

      case 'o':
        outputFile = optarg;
        break;

      case 't':
        testFile = optarg;
        break;

      case 'g':
        dumpFile = optarg;
        break;

      case 'q':
        queryFile = optarg;
        break;

      case 'I':
        if (!strcmp(optarg, "ShortIndexing")) {
          indexMethod = ShortIndexing;
        } else if (!strcmp(optarg, "SuccessorOrder")) {
          indexMethod = SuccessorOrder;
        } else if (!strcmp(optarg, "Standard")) {
          indexMethod = Standard;
        } else if (!strcmp(optarg, "LabelOrder")) {
          indexMethod = LabelOrder;
        } else {
          cerr << "ERROR: Unknown argument to --index-method.\n";
          printHelpMessage();
          exit(EXIT_FAILURE);
        }
        break;

      case 'S':
        if (!strcmp(optarg, "DFS")) {
          searchMethod = DFS;
        } else if (!strcmp(optarg, "BBFS")) {
          searchMethod = BBFS;
        } else {
          cerr << "ERROR: Unknown argument to --search-method.\n";
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

      case 'f':
        if (!isdigit(optarg[0])) {
          cerr << "ERROR: The -C | --coarsen-factor argument is not a number\n";
          printHelpMessage();
          exit(EXIT_FAILURE);
        } else {
          coarsenFactor = atoi(optarg);

          if (coarsenFactor < 2)
            coarsenFactor = 2;
        }
        break;

      case 'C':
        if (!strcmp(optarg, "Greedy")) {
          coarsenMethod = Greedy;
        } else {
          cerr << "ERROR: Unknown argument to --coarsen-method.\n";
          printHelpMessage();
          exit(EXIT_FAILURE);
        }
        break;

      case 'F':
        coarseFile = optarg;
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

  if (inputFile.size() == 0) {
    cerr << "ERROR: No input file was specified.\n\n";
    printHelpMessage();
    exit(EXIT_FAILURE);
  }

  if (!inputFile.substr(inputFile.size() - 4, 4).compare(".dot")) {
    graph = parseDotFile(inputFile.c_str(), uniqueFlag);
  } else if (!inputFile.substr(inputFile.size() - 4, 4).compare(".gra")) {
    graph = parseGraFile(inputFile.c_str(), uniqueFlag);
  } else {
    cerr << "ERROR: Unknown graph-file extension. Only accepts .dot and .gra "
      << "files.\n";
    exit(EXIT_FAILURE);
  }

  if (!graph) {
    cerr << "ERROR: No graph was generated from the input file.\n\n";
    printHelpMessage();
    exit(EXIT_FAILURE);
  }

  // Set the indexing method as specified / default
  graph->setIndexMethod(indexMethod);

  // Start worker threads on the graph
  graph->enableQueries();

  // When necessary create a coarsened graph and dump it
  if (coarsenFactor > 1) {
    map<int, int> vertexMap;
    Query * query = NULL;
    CoarsenQuery * cQuery = createCQuery();
    cQuery->setFactor(coarsenFactor);
    cQuery->setSecondaryFactor(coarsenSecondaryFactor);
    cQuery->setMethod(coarsenMethod);

    graph->pushQuery(cQuery);

    query = graph->pullQuery(true);

    if (!query) {
      cerr << "ERROR: Coarsen query dissapeared into the void!\n";
      exit(EXIT_FAILURE);
    }

    if (coarseFile.size() > 0) {
      fstream coarseMapStream;
      string mapFile = coarseFile + ".map";

      cQuery->getAnswer()->printToFile(coarseFile.c_str(), true);
      coarseMapStream.open(mapFile.c_str(), ios_base::out);

      coarseMapStream << "General to coarsened graph index mapping";

      map<int, int>& vertexMap = cQuery->getMap();
      for (auto it = vertexMap.begin(), end = vertexMap.end(); it != end; ++it)
        coarseMapStream << "\n" << it->first << " : " << it->second;

      coarseMapStream.close();
    }

    delete cQuery->getAnswer();
    delete cQuery;
    delete graph;

    exit(EXIT_SUCCESS);
  }

  // Dump the parsed graph for verification purposes
  if (dumpFile.size() > 0) {
    fstream dumpStream;
    Vertex * curr;

    dumpStream.open(dumpFile.c_str(), ios_base::out);

    dumpStream << "graph_for_greach\n" << graph->getVertexCount() << "\n";
    for (int i = 0, e = graph->getVertexCount(); i < e; i++) {
      if (!(curr = graph->getVertex(i)))
        continue;

      dumpStream << i << ":";

      for (int i2 = 0, e2 = curr->getSuccessorCount(); i2 < e2; i2++)
        dumpStream << " " << curr->getSuccessor(i2)->getId();

      dumpStream << " #\n";
    }

    dumpStream.close();
  }

  // Exit here in case of a 'dry' run
  if (dryFlag) {
    delete graph;
    exit(EXIT_SUCCESS);
  }

  // Process the queries with two threads
  thread pushThread(queryGenerator, graph, testFile.c_str(), searchMethod);
  thread pullThread(resultAnalysis, graph, queryFile.c_str());

  // Wait for the queries to be issued and treated
  pushThread.join();
  graph->disableQueries();
  pullThread.join();

  // Prepare output stream for statistics and benchmarks
  ostream * output = &cout;
  if (outputFile.size() > 0)
    output = new fstream(outputFile.c_str(), ios_base::out);

  // Dump the statistics of the queries...
  if (graph->getStatisticsEnabled())
    graph->printStatistics(*output);

  // ...and the benchmarks
  if (graph->getBenchmarksEnabled())
    graph->printBenchmarks(*output);

  // Clear up the queries
  for (auto it = queryMap.begin(), end = queryMap.end(); it != end; ++it)
    delete it->first;

  queryMap.clear();

  delete graph;

	//exit(EXIT_SUCCESS);
  return 0;
}
