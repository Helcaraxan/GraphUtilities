#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <string>
#include <queue>
#include <set>
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

   vector<int> trace;
   for(int i : order) {
      for(int j : graph.getPredecessors(i))
         trace.push_back(j);
      trace.push_back(i);
   }
   vector<int> reuse(trace.size(), 0);
   vector<bool> seen(graph.getNbNodes(), false);
   vector<int> stack;
   int maxReuse = 0;
   for(int i = 0, s = trace.size(); i < s; i++) {
      if(!seen[trace[i]])
         reuse[i] = -1;
      else if(!stack.empty())
         for(auto it = stack.end()-1, e = stack.begin() ; it != e; it--) {
            if(*it == trace[i]) {
               stack.erase(it);
               break;
            }
            reuse[i]++;
         }
      stack.push_back(trace[i]);
      seen[trace[i]] = true;
      maxReuse = max(maxReuse, reuse[i]);
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
