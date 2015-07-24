#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <string>
#include <queue>
#include <map>
#include <climits>

#include "graph.hpp"
#include "patoh.h"
//#include "metis_part.hpp"

using namespace std;

vector<vector<int>> parseGraFile(FILE* f)
{
   int n;
   char dump[2];
   fscanf(f, "%d\n", &n);
   vector<vector<int>> succ(n);
   for(int i=0; i<n; i++) {
      int j;
      fscanf(f, "%*d: ");
      while(fscanf(f, "%d ", &j) == 1) {
         succ[i].push_back(j);
      }
      if(fscanf(f, "%1[#]", dump) != 1) {
         printf("Wrong file format\n");
         exit(1);
      }
   }
   return succ;
}

void output_dot(const char* filename, Graph g, vector<int> part) 
{
   const char* colors[] = {"aliceblue", "antiquewhite", "aquamarine", "azure", "beige", "bisque", "black", "blanchedalmond", "blue", "blueviolet", "brown", "burlywood", "cadetblue", "chartreuse", "chocolate", "coral", "cornflowerblue", "cornsilk", "crimson", "cyan", "darkgoldenrod", "darkgreen", "darkkhaki", "darkolivegreen", "darkorange", "darkorchid", "darksalmon", "darkseagreen", "darkslateblue", "darkslategray", "darkslategrey", "darkturquoise", "darkviolet", "deeppink", "deepskyblue", "dimgray", "dimgrey", "dodgerblue", "firebrick", "floralwhite", "forestgreen", "gainsboro", "ghostwhite", "gold", "goldenrod", "green", "hotpink", "indianred", "indigo", "invis", "ivory", "khaki", "lavender", "lavenderblush", "lawngreen", "lemonchiffon", "lightblue", "lightcoral", "lightcyan", "lightgoldenrod", "lightgoldenrodyellow", "lightgray", "lightgrey", "lightpink", "lightsalmon", "lightseagreen", "lightskyblue", "lightslateblue", "lightslategray", "lightslategrey", "lightsteelblue", "lightyellow", "limegreen", "linen", "magenta", "maroon", "mediumaquamarine", "mediumblue", "mediumorchid", "mediumpurple", "mediumseagreen", "mediumslateblue", "mediumspringgreen", "mediumturquoise", "mediumvioletred", "midnightblue", "mintcream", "mistyrose", "moccasin", "navajowhite", "navy", "navyblue", "none", "oldlace", "olivedrab", "orange", "orangered", "orchid", "papayawhip", "peachpuff", "peru", "pink", "plum", "powderblue", "purple", "red", "rosybrown", "royalblue", "saddlebrown", "salmon", "sandybrown", "seagreen", "seashell", "sienna", "skyblue", "slateblue", "slategray", "slategrey", "snow", "springgreen", "steelblue", "tan", "thistle", "tomato", "transparent", "turquoise", "violet", "violetred", "wheat", "white", "whitesmoke", "yellow", "yellowgreen"};
   //randomize colors
   for(int i=0; i<2000; i++)
   {
      int a = rand() % 132, b = rand() % 132;
      swap(colors[a], colors[b]);
   }

   FILE* fo = fopen(filename, "w");
   if(!fo) exit(1);
   fprintf(fo, "digraph {\n");
   for(int i=0, n = part.size(); i<n; i++)
      fprintf(fo, "\t%d [fillcolor=%s, style=filled]\n", i, colors[part[i]]);
   for(int i=0; i<g.getNbNodes(); i++)
      for(int s : g.getSuccessors(i))
         fprintf(fo, "\t%d -> %d\n [color=%s]", i, s, part[i] == part[s] ? "black" : "red");
   fprintf(fo, "}");
}
/*
vector<int> mtmetis_bisect(Graph graph)
{
   //printf("conversion to MeTis format\n");
   graph_t* metisGraph = graph.getMetisGraph();
   //printf("MeTis bisection\n");

   double* opts = mtmetis_init_options();
   opts[MTMETIS_OPTION_NPARTS] = 2;
   opts[MTMETIS_OPTION_VERBOSITY] = 0;
   opts[MTMETIS_OPTION_NTHREADS] = 4;
   mtmetis_pid_t* part = new mtmetis_pid_t[metisGraph->nvtxs];
   int err = mtmetis_partition_explicit(metisGraph->nvtxs, metisGraph->xadj, metisGraph->adjncy, NULL, NULL, opts, part, NULL);
   if(err != MTMETIS_SUCCESS) {
      printf("mt-metis error during partitioning: %s\n", err == MTMETIS_ERROR_INVALIDINPUT ? "Invalid input" : "Not enough memory");
      exit(1);
   }
   vector<int> vpart(part, part+graph.getNbNodes());
   return vpart;
}
*/

pair<vector<int>, int> patoh_bisect(Graph graph)
{
   PPaToH_Parameters pargs = new PaToH_Parameters;
   PaToH_Initialize_Parameters(pargs, PATOH_CONPART, PATOH_SUGPARAM_DEFAULT);

   hgraph_t* hg = graph.getHyperGraph();
   /*
   printf("%d %d\n", hg->nvtxs, hg->nhedges);
   for(int i = 0; i < hg->nhedges; i++) {
      for(int j = hg->eptr[i]; j < hg->eptr[i+1]; j++)
         printf("%d ", hg->eind[j]);
      printf("\n");
   }
   */

   int* ones = new int[2*hg->nhedges];
   fill(ones, ones+hg->nvtxs, 1);
   PaToH_Alloc(pargs, hg->nvtxs, hg->nhedges, 1, ones, ones, hg->eptr, hg->eind);

   int cut;
   int * tpart = new int[hg->nvtxs];
   float tweight[] = {.5, .5};
   int pweight[] = {0, 0};
   pargs->_k = 2;
   PaToH_Part(pargs, hg->nvtxs, hg->nhedges, 1, 0, ones, ones, hg->eptr, hg->eind, tweight, tpart, pweight, &cut);
   delete [] ones;
   delete pargs;
   PaToH_Free();
   vector<int> part(tpart, tpart+hg->nvtxs);
   return make_pair(part, cut);
}
