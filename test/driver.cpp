#include <map>
#include <chrono>
#include <random>
#include <thread>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <iostream>

#include <getopt.h>
#include "tbb/concurrent_hash_map.h"

#include "graph.hpp"

using namespace std;

#define QUERY_NB 100000

int dryFlag = 0;
int uniqueFlag = 0;
int verifyFlag = 0;
bool poison = false;
semaphore openQueries;
tbb::concurrent_hash_map<Graph::Query *, bool> queryMap;
vector<Graph::Query *> results;


const struct option longopts[] = {
  {"input",     required_argument, 0, 'i'},
  {"output",    required_argument, 0, 'o'},
  {"test",      required_argument, 0, 't'},
  {"method",    required_argument, 0, 'm'},
  {"search",    required_argument, 0, 's'},
  {"graph",     optional_argument, 0, 'g'},
  {"queries",   optional_argument, 0, 'q'},
  {"help",      no_argument,       0, 'h'},
  {"verify",  no_argument, &verifyFlag, 1},
  {"unique-edges", no_argument, &uniqueFlag, 1},
  {"dry",       no_argument,  &dryFlag, 1},
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
  cout << "\nBoolean options:\n";
  cout << "\t-h | --help\t\tDisplay this help message\n";
  cout << "\t-v | --verify\t\tVerify the query results by a DFS query that ignores labeling\n";
  cout << "\t-u | --unique-edges\tDon't check for double-edges on the input graph (speeds-up parsing of large grpahs)\n";
  cout << "\t-d | --dry\t\tDo not perform queries stop after graph condensation (and eventual dumping)\n";
}


void
queryGenerator(Graph * graph, fstream * testFile, Graph::SearchMethod method) {
  int i, a, b, res;
  Graph::Query * newQuery = NULL;
  tbb::concurrent_hash_map<Graph::Query *, bool>::accessor queryAccess;

  // Fill the queries to process...
  if (testFile->is_open()) {
    cout << "Parsing queries from a test file.\n";
    // ... from the specified file
    while (testFile->good()) {
      (*testFile) >> a >> b >> res;

      newQuery = new Graph::Query(graph->getVertexFromId(a), graph->getVertexFromId(b), method);
      queryMap.insert(queryAccess, pair<Graph::Query *, bool>(newQuery, (res == 0 ? false : true)));
      queryAccess.release();
      graph->pushQuery(newQuery);
      openQueries.post();
    }

    cout << "Finished parsing queries from file.\n\n";

    testFile->close();
  } else {
    // ... at random
    default_random_engine generator(chrono::system_clock::now().time_since_epoch().count());
    uniform_int_distribution<int> queryDistribution(0, graph->getVertexCount() - 1);

    for (i = 0; i < QUERY_NB; i++) {
      a = queryDistribution(generator);
      do {
        b = queryDistribution(generator);
      } while (a == b);

      if (a < b)
        newQuery = new Graph::Query(graph->getVertexFromId(a), graph->getVertexFromId(b), method);
      else
        newQuery = new Graph::Query(graph->getVertexFromId(b), graph->getVertexFromId(a), method);

      queryMap.insert(queryAccess, pair<Graph::Query *, bool>(newQuery, false));
      queryAccess.release();
      graph->pushQuery(newQuery);
      openQueries.post();
    }
  }

  poison = true;
  openQueries.post();
  return;
}


void
resultAnalysis(Graph * graph, fstream * queryFile) {
  int resultCounter = 0;
  Graph::Query * nextQuery = NULL;
  tbb::concurrent_hash_map<Graph::Query *, bool>::const_accessor queryAccess;


  cout << "Analyzed results: " << resultCounter;
  while (true) {
    openQueries.wait();

    if (poison && (queryMap.size() == results.size()))
      break;

    nextQuery = graph->pullResult();
    queryMap.find(queryAccess, nextQuery);

    if (nextQuery->isError()) {
      cerr << "Error when processing query " << nextQuery->getSource()->id << " -> " << nextQuery->getTarget()->id << "\n";
    } else if (verifyFlag) {
      if (nextQuery->getAnswer() != queryAccess->second)
        cerr << "Wrong answer for query " << nextQuery->getSource()->id << " -> " << nextQuery->getTarget()->id << "\n";
    } else if (queryFile->is_open()) {
      (*queryFile) << nextQuery->getSource()->id << " " << nextQuery->getTarget()->id << " ";
      if (queryAccess->second)
        (*queryFile) << "1\n";
      else
        (*queryFile) << "0\n";
    }
    queryAccess.release();
    results.push_back(nextQuery);
    resultCounter++;
    cout << "\rAnalyzed results: " << resultCounter;
  }

  cout << "\n";

  if (queryFile->is_open())
    queryFile->close();
}


int
main(int argc, char * argv[]) {
	Graph * graph = NULL;
  Graph::IndexMethod indexMethod = Graph::Standard;
  Graph::SearchMethod searchMethod = Graph::Undefined;
  fstream outputFile, testFile, dumpFile, queryFile;
  ostream * output = &cout;
  int i, c;
  char fileName[512] = {'\0'};

  // Parse command-line options
  while ((c = getopt_long(argc, argv, "i:o:t:m:s:g::q::hvud", longopts, NULL)) != -1) {
    switch (c) {
      case 'i':
        strncpy(fileName, optarg, 511);
        break;

      case 'o':
        outputFile.open(optarg, fstream::out);
        if (!outputFile.good()) {
          cerr << "Error: could not open the output file.\n";
          exit(EXIT_FAILURE);
        }

        output = &outputFile;
        break;

      case 't':
        testFile.open(optarg, fstream::in);
        if (!testFile.good()) {
          cerr << "Error: could not open the test file.\n";
          exit(EXIT_FAILURE);
        }
        break;

      case 'm':
        if (!strcmp(optarg, "ShortIndexing")) {
          indexMethod = Graph::ShortIndexing;
        } else if (!strcmp(optarg, "SuccessorOrder")) {
          indexMethod = Graph::SuccessorOrder;
        } else if (!strcmp(optarg, "Standard")) {
          indexMethod = Graph::Standard;
        } else if (!strcmp(optarg, "LabelOrder")) {
          indexMethod = Graph::LabelOrder;
        } else {
          cerr << "Error: unknown label method specified.\n";
          printHelpMessage();
          exit(EXIT_FAILURE);
        }
        break;

      case 's':
        if (!strcmp(optarg, "DFS")) {
          searchMethod = Graph::DFS;
        } else if (!strcmp(optarg, "BBFS")) {
          searchMethod = Graph::BBFS;
        } else {
          cerr << "Error: unknown search method specified.\n";
          printHelpMessage();
          exit(EXIT_FAILURE);
        }
        break;

      case 'g':
        dumpFile.open(optarg, fstream::out);
        if (!dumpFile.good()) {
          cerr << "Error: could not open the graph dump file.\n";
          exit(EXIT_FAILURE);
        }
        break;

      case 'q':
        queryFile.open(optarg, fstream::out);
        if (!queryFile.good()) {
          cerr << "Error: could not open the query dump file.\n";
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
        cerr << "Error in GetOpt class while parsing command-line arguments.\n";
        exit(EXIT_FAILURE);
    }
  }

  if (fileName[0] == '\0') {
    cerr << "Error: No input file was specified.\n\n";
    printHelpMessage();
    exit(EXIT_FAILURE);
  }

  if (strstr(fileName, ".dot")) {
    graph = Graph::createFromDotFile(fileName, uniqueFlag);
  } else if (strstr(fileName, ".gra")) {
    graph = Graph::createFromGraFile(fileName, uniqueFlag);
  } else {
    cerr << "Error: Unknown graph-file extension. Only accepts .dot and .gra files.\n";
    exit(EXIT_FAILURE);
  }

  if (!graph) {
    cerr << "Error: No graph was generated from the input file.\n\n";
    printHelpMessage();
    exit(EXIT_FAILURE);
  }
  
  // Set the indexing method as specified / default
  graph->setIndexMethod(indexMethod);

  // Dump the parsed graph for verification purposes
  if (dumpFile.is_open()) {
    Vertex * curr;

    cout << "Dump file for consistency check.\n";
    dumpFile << "graph_for_greach\n" << graph->getVertexCount() << "\n";
    for (i = 0; i < (int) graph->getVertexCount(); i++) {
      dumpFile << i << ":";
      curr = graph->getVertexFromId(i);
      for (auto it = curr->successors_begin(); it != curr->successors_end(); ++it)
        dumpFile << " " << (*it)->id;

      dumpFile << " #\n";
    }

    cout << "Finished dumping file.\n\n";
    dumpFile.close();
  }

  // Exit here in case of a 'dry' run
  if (dryFlag)
    exit(EXIT_SUCCESS);

  // Process the queries with two threads
  thread pushThread(queryGenerator, graph, &testFile, searchMethod);
  thread pullThread(resultAnalysis, graph, &queryFile);

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

  // Clear up the query map
  for (auto it = queryMap.begin(), end = queryMap.end(); it != end; ++it)
    delete it->first;

  queryMap.clear();
  delete graph;

	exit(EXIT_SUCCESS);
}
