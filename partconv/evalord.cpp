#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <string>
#include <queue>
#include <map>
#include <climits>

#include "graph.hpp"

using namespace std;

int main(int argc, char* argv[])
{
   srand(24);
   //parsing command line args
   if(argc < 2) {
      printf("No graph file given\n");
      return 1;
   }

   FILE* fgraph = fopen(argv[1], "r");
   if(!fgraph) {
      printf("Could not open file\n");
      return 1;
   }

   vector<vector<int>> v = parseGraFile(fgraph);
   Graph graph(v);
   fclose(fgraph);

   vector<int> order(graph.getNbNodes());
   if(argc < 3) {
      for(int i = 0; i < graph.getNbNodes(); i++)
         order[i] = i;
   }
   else {
      FILE* forder = fopen(argv[2], "r");
      if(!forder)
         return 1;
      for(int i = 0; i < graph.getNbNodes(); i++)
         fscanf(forder, "%d\n", &order[i]);
      fclose(forder);
   }

   vector<int> use(graph.getNbNodes(), -1);
   vector<int> reuse;
   int maxReuse = 0;
   for(int i : order) {
      for(int j : graph.getPredecessors(i)) {
         if(use[j] == -1)
            reuse.push_back(-1);
         else
            reuse.push_back(reuse.size() - use[j]);
         use[j] = reuse.size();
         maxReuse = max(maxReuse, reuse.back());
      }
   }
   vector<int> nbReuse(maxReuse+1, 0);
   for(int r : reuse)
      if(r != -1)
         nbReuse[r]++;
   for(int r = maxReuse; r>=0; r--)
      nbReuse[r] += nbReuse[r+1];
   for(int i = 0; i < maxReuse; i += 50)
      printf("%d %d\n", i, nbReuse[i]);
}
