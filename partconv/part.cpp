#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <string>
#include <queue>
#include <map>
#include <set>
#include <climits>
#include <iostream>

//#include "metis_part.hpp"
#include "graph.hpp"

using namespace std;

const int NB_TRY_BISECT = 3;

struct TreeExt {
  TreeExt* succ[2];
  int id, size, nbCuts;

  TreeExt(int id = -1, int s = 0, int n = 0, TreeExt* s1 = NULL, TreeExt* s2 = NULL) : id(id), size(s), nbCuts(n) {
    succ[0] = s1;
    succ[1] = s2;
  }

  ~TreeExt() {
    for(int i = 0; i < 2; i++)
      delete succ[i];
  }
};

void print_treee_p(TreeExt* t, string p0, string p1)
{
  if(t->id != -1)
    cout << p0 << t->id << ' ' << t->nbCuts << '\n' << p1 << '\n';
  else {
    if(t->succ[0]) {
      string p00 = p0 + to_string(t->nbCuts) + "--";
      string p01 = p1 + (t->succ[1] ? "|  " : "   "); 
      print_treee_p(t->succ[0], p00, p01);
    }
    if(t->succ[1]) {
      string p10 = p1 + "|--";
      string p11 = p1 + "   ";
      print_treee_p(t->succ[1], p10, p11);
    }
  }
}

void print_treee(TreeExt* t)
{
  print_treee_p(t, "", "");
}

void partStats(Graph g, vector<int> part) {
   int nbParts = 0;
   for(int j : part)
      nbParts = max(j+1, nbParts);
   int nbCut = 0;
   int nbCut0 = 0;
   for(int i = 0; i < g.getNbNodes(); i++)
      for(int j : g.getSuccessors(i)) {
         nbCut += part[i] != part[j];
         if(part[i] < part[j]) nbCut0++;
      }
   int avgSize = g.getNbNodes() / nbParts;
   printf("Partition stats:\n  params: %d nodes\n", g.getNbNodes());
   printf("\n  nb parts: %d\n  average size: %d\n", nbParts, avgSize);
   printf("  nb cuts: %d\t  cuts/node: %lf\t cuts/part: %lf\n", nbCut, (double)nbCut/(double)g.getNbNodes(), (double)nbCut/(double)nbParts);
   printf(" %% cuts 0 -> 1: %lf\n", (double)nbCut0/(double)nbCut);
}

vector<int> convBisect(Graph g) {
   // we determine the best order for the two parts

   pair<vector<int>, int> p = patoh_bisect(g);
   vector<int>& part = p.first;
   //return part;
   int nbCuts[2] = {0, 0};
   for(int i = 0; i < g.getNbNodes(); i++) {
      int b = 0;
      for(int j : g.getSuccessors(i))
         if(part[i] != part[j])
            b = 1;
      nbCuts[part[i]] += b;
   }
   //if there are more cuts from 1 to 0, we switch parts
   if(nbCuts[1] > nbCuts[0])
      for(int i = 0; i < g.getNbNodes(); i++)
         part[i] = 1 - part[i];
   
   //printf("nb nodes : %d\n", g.getNbNodes());
   /*
   printf("Initial partition :\n");
   for(int i=0; i<g.getNbNodes(); i++)
      printf("%d : %d\n", i, part[i]);
      */
   
   vector<int> nbPred(g.getNbNodes(), 0);
   vector<int> nbSucc(g.getNbNodes(), 0);
   vector<int> convPart(g.getNbNodes(), 2);
   queue<int> ready;
   for(int i = 0 ; i < g.getNbNodes(); i++) {
      if(part[i] == 0 && g.getPredecessors(i).empty()) {
         ready.push(i);
         convPart[i] = 0;
      }
      if(part[i] == 1 && g.getSuccessors(i).empty()) {
         ready.push(i);
         convPart[i] = 1;
      }
      nbPred[i] = g.getPredecessors(i).size();
      nbSucc[i] = g.getSuccessors(i).size();
   }
   int nbOk[2] = {0, 0};
   while(!ready.empty()) {
      int nxt = ready.front(); 
      ready.pop();

      nbOk[part[nxt]]++;
      vector<int> adj = part[nxt] == 0 ? g.getSuccessors(nxt) : g.getPredecessors(nxt);
      vector<int>& nbAdj = part[nxt] == 0 ? nbPred : nbSucc;
      for(int j : adj)
         if(part[j] == part[nxt]) {
            nbAdj[j]--;
            if(convPart[j] == 2 && nbAdj[j] == 0){
               ready.push(j);
               convPart[j] = part[j];
            }
         }
   }

   //printf("nbOk : %d (%lf)   %d (%lf)\n", nbOk[0], (double)nbOk[0]/(double)g.getNbNodes(), nbOk[1], (double)nbOk[1]/(double)g.getNbNodes());
  // printf("rest : %d (%lf)\n", g.getNbNodes() - nbOk[0] - nbOk[1], (double)(g.getNbNodes() - nbOk[0] - nbOk[1])/(double)g.getNbNodes());
   if(nbOk[0] + nbOk[1] < g.getNbNodes()) {
      vector<int> newIds;
      vector<int> oldIds(g.getNbNodes());
      vector<vector<int>> newSucc;
      int curId = 0;
      for(int i = 0; i < g.getNbNodes(); i++)
         if(convPart[i] == 2) {
            newIds.push_back(i);
            oldIds[i] = curId++;
         }
      newSucc.resize(curId);
      for(int i = 0; i < curId; i++)
         for(int j : g.getSuccessors(newIds[i]))
            if(convPart[j] == 2) {
               newSucc[i].push_back(oldIds[j]); 
            }
      Graph newG(newSucc);
      vector<int> newConv = convBisect(newG);
      for(int i = 0; i < curId; i++)
         convPart[newIds[i]] = newConv[i];
   }

   /*
   printf("Convex partition :\n");
   for(int i=0; i<g.getNbNodes(); i++)
      printf("%d : %d\n", i, conv[i]);
      */
   return convPart;
}

void refineBisect(Graph& g, vector<int>& part)
{
   vector<int> delta(g.getNbNodes(), 0);
   vector<int> border(g.getNbNodes(), false);
   vector<int> nbSuccOther(g.getNbNodes(), 0);
   vector<int> nbPredOther(g.getNbNodes(), 0);
   vector<bool> moved(g.getNbNodes(), false);
   map<int, set<int> > toMove[2];
   vector<bool> waiting(g.getNbNodes(), false);

   int partSize[2] = {0, 0};

   for(int i = 0; i < g.getNbNodes(); i++) {
      partSize[part[i]]++;
      for(int j: g.getSuccessors(i)) 
         if(part[i] != part[j]) {
            nbSuccOther[i]++;
            border[i] = true;
         }
      for(int j: g.getPredecessors(i))
         if(part[i] != part[j]) {
            nbPredOther[i]++;
            border[i] = true;
         }
   }
   for(int i = 0; i < g.getNbNodes(); i++) 
      if(border[i]) {
         if(part[i] == 0) {
            delta[i] = -1;
            for(int j : g.getPredecessors(i))
               delta[i] += !border[j];
            if(nbSuccOther[i] == (int)g.getSuccessors(i).size() && delta[i] <= 0) {
               toMove[0][delta[i]].insert(i);
               waiting[i] = true;
            }
         }
         else {
            delta[i] = 1;
            for(int j : g.getPredecessors(i))
               if(nbSuccOther[j] == 1)
                  delta[i]--;
            if(nbPredOther[i] == (int)g.getPredecessors(i).size() && delta[i] <= 0) {
               toMove[1][delta[i]].insert(i);
               waiting[i] = true;
            }
         }
      }

   while(1) {
      pair<int, set<int>> min0, min1;
      min0.first = g.getNbNodes();
      min1.first = g.getNbNodes();
      if(toMove[0].begin() != toMove[0].end())
         min0 = *toMove[0].begin();
      if(toMove[1].begin() != toMove[1].end())
         min1 = *toMove[1].begin();
      int id;
      if(partSize[0] >= partSize[1] && min0.first < 0 && min0.first <= min1.first) {
         id = *min0.second.begin();
         min0.second.erase(id);
         if(min0.second.empty())
            toMove[0].erase(min0.first);
         part[id] = 1;
         for(int j : g.getSuccessors(id)) {
            toMove[1][delta[j]].erase(j);
            if(toMove[1][delta[j]].empty())
               toMove[1].erase(delta[j]);
            nbPredOther[j]--;
            if(nbPredOther[j] == 0)
               border[j] = false;
            if(g.getSuccessors(id).size() == 1)
               delta[j]--;
         }
         for(int j : g.getPredecessors(id)) {
            if(!border[j])
               for(int k : g.getSuccessors(j)) {
                  delta[k]++;
                  if(delta[k] <= 1) {
                     toMove[0][delta[k]-1].erase(k);
                     if(toMove[0][delta[k]-1].empty())
                        toMove[0].erase(delta[k]-1);
                  }
                  if(delta[k] <= 0)
                     toMove[0][delta[k]].insert(k);
               }
            border[j] = true;
            nbSuccOther[j]++;
            if(!moved[j] && nbSuccOther[j] == (int)g.getSuccessors(j).size())
               toMove[0][delta[j]].insert(j);
         }
         partSize[0]--;
         partSize[1]++;
      }
      else if(partSize[1] <= partSize[0] && min1.first < min0.first && min1.first < 0) {
         id = *min1.second.begin();
         min1.second.erase(id);
         if(min1.second.empty())
            toMove[1].erase(min1.first);
         for(int j : g.getSuccessors(id)) {
            toMove[0][delta[j]].erase(j);
            if(toMove[0][delta[j]].empty())
               toMove[0].erase(delta[j]);
            nbPredOther[j]++;
            if(!moved[j] && nbPredOther[j] == (int)g.getPredecessors(j).size())
               toMove[1][delta[j]].insert(j);
         }

         for(int j : g.getPredecessors(id)) {
            nbSuccOther[j]--;
            if(nbSuccOther[j] == 0) {
               border[j] = false;
               for(int k : g.getSuccessors(j)) {
                  delta[k]--;
                  if(delta[k] <= 1) {
                     toMove[0][delta[k]-1].erase(k);
                     if(toMove[0][delta[k]-1].empty())
                        toMove[0].erase(delta[k]-1);
                  }
                  if(delta[k] <= 0)
                     toMove[0][delta[k]].insert(k);
               }
            }
         }
         partSize[0]++;
         partSize[1]--;
      }
      else break;
      moved[id] = true;
   }
}

TreeExt* hierPart(Graph g, Graph* orig, vector<int> ids, int level, vector<int>* nbCutsLevel)
{
   if(g.getNbNodes() == 0)
      return NULL;
   if(g.getNbNodes() == 1)
      return new TreeExt(ids[0], 1, orig->getPredecessors(ids[0]).empty() ? 0 : 1);
   vector<int> bestBisect = convBisect(g);
   int bestCut = INT_MAX;
   if(g.getNbNodes() > 100)
      for(int t = 0; t < NB_TRY_BISECT; t++) {
         vector<int> bisect = convBisect(g);
         //refineBisect(g, bisect);
         int cut = 0;
         for(int i = 0; i < g.getNbNodes(); i++)
            for(int j : g.getSuccessors(i))
               cut += bisect[i] != bisect[j];
         if(cut < bestCut) {
            bestCut = cut;
            bestBisect = bisect;
         }
      }

   if((int)(*nbCutsLevel).size() <= level) {
      printf("begin level %d\n", level);
      (*nbCutsLevel).push_back(0);
   }

   vector<int> newIds[2];
   vector<int> oldIds(g.getNbNodes());
   for(int i = 0; i < g.getNbNodes(); i++) {
      oldIds[i] = newIds[bestBisect[i]].size();
      newIds[bestBisect[i]].push_back(ids[i]);
   }
   vector< vector<int> > succ[2];
   succ[0].resize(newIds[0].size());
   succ[1].resize(newIds[1].size());
   int nbCuts = 0;
   for(int i = 0; i < g.getNbNodes(); i++) {
      nbCuts += orig->getPredecessors(ids[i]).size() > g.getPredecessors(i).size();
      for(int j : g.getSuccessors(i)) {
         if(bestBisect[j] == bestBisect[i])
            succ[bestBisect[i]][oldIds[i]].push_back(oldIds[j]);
         else {
            (*nbCutsLevel)[level]++;
         }
      }
   }
   Graph g0(succ[0]);
   Graph g1(succ[1]);
   TreeExt* t0 = hierPart(g0, orig, newIds[0], level+1, nbCutsLevel);
   TreeExt* t1 = hierPart(g1, orig, newIds[1], level+1, nbCutsLevel);
   return new TreeExt(-1, g.getNbNodes(), nbCuts, t0, t1);
}

void write_ordering(FILE* f, TreeExt* t)
{
   if(t == NULL)
      return;
   if(t->id != -1)
      fprintf(f, "%d\n", t->id);
   write_ordering(f, t->succ[0]);
   write_ordering(f, t->succ[1]);
}

pair<int, int> calcCut(TreeExt* t, vector<int>& v)
{
   if(t == NULL)
      return make_pair(0, 0);
   pair<int, int> p0 = calcCut(t->succ[0], v);
   pair<int, int> p1 = calcCut(t->succ[1], v);
   int r = t->nbCuts - p0.first - p1.first;
   int m = max(t->nbCuts, max(p0.second, p1.second));
   //printf("%d  %d %d %d\n", t->id, t->nbCuts, r, m);
   v[m] += r;
   return make_pair(t->nbCuts, m);
}

int main(int argc, char* argv[])
{
   srand(time(NULL));
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

   vector<int> convPart = convBisect(graph);

   partStats(graph, convPart);
   vector<int> ids(graph.getNbNodes());
   for(int i = 0; i < graph.getNbNodes(); i++)
      ids[i] = i;
   vector<int> nbCutsLevel(0, 0);
   TreeExt* t = hierPart(graph, &graph, ids, 0, &nbCutsLevel);

   printf("nb cuts : \n");
   int tot = 0;
   int nbParts = 1;
   for(int l = 0; l < (int)nbCutsLevel.size(); l++) {
      nbParts *= 2;
      tot += nbCutsLevel[l];
      printf("  level %d : total %d (%lf)\n", l, tot, (double)tot / (double)(nbParts));
   }

   FILE* f_ord = fopen("order.txt", "w");
   write_ordering(f_ord, t);
   fclose(f_ord);
   //print_treee(t);
   //output_dot("nonconv.dot", graph, nonConv);
   //output_dot("conv.dot", graph, convPart);

   vector<int> diffCutS(graph.getNbNodes(), 0);
   calcCut(t, diffCutS);
   for(size_t i = 1; i < diffCutS.size(); i++) {
      diffCutS[i] += diffCutS[i-1];
   }
   FILE* f = fopen("cuts.txt", "w");
   for(int i = 1; i < graph.getNbNodes(); i++) {
      fprintf(f, "%d %d\n", i, diffCutS[i]);
      if(diffCutS[i] == 0)
          break;
   }
   fclose(f);

   return 0;
}

