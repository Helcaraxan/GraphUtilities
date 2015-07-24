#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <string>

#include <patoh.h>
//#include "metis_part.hpp"

using namespace std;

struct hgraph_t {
   int nvtxs, nhedges;
   int *eptr, *eind;
};

class Graph{
   private:
      int nbNodes;
      int nbEdges;
      vector<vector<int>> successors;
      vector<vector<int>> predecessors;

   public:
      Graph() : nbNodes(0), nbEdges(0) {}
      Graph(vector<vector<int>>& succ) {
         nbEdges = 0;
         nbNodes = succ.size();
         successors = succ;
         predecessors.resize(nbNodes);
         for(int i = 0; i < nbNodes; i++)
            for(int node : succ[i]) {
               predecessors[node].push_back(i);
               nbEdges++;
            }
      }

      int getNbNodes(){
         return nbNodes;
      }

      int getDeltaEdges(int i) {
         return predecessors[i].size() - successors[i].size();
      }

      vector<int> getSuccessors(int i){
         return successors[i];
      }

      vector<int> getPredecessors(int i){
         return predecessors[i];
      }

      /*
      graph_t* getMetisGraph()
      {
         mtmetis_adj_t* xadj = new mtmetis_adj_t[nbNodes + 1];
         mtmetis_vtx_t* adjncy = new mtmetis_vtx_t[2*nbEdges];
         unsigned iadj = 0;
         xadj[0] = 0;
         for(int node = 0; node < nbNodes; node ++) {
            for(int adj : successors[node])
               adjncy[iadj++] = adj;
            for(int adj : predecessors[node])
               adjncy[iadj++] = adj;
            xadj[node + 1] = iadj;
         }
         graph_t* g = new graph_t;
         g->nvtxs = nbNodes;
         g->xadj = xadj;
         g->adjncy = adjncy;
         return g;
      }
      */

      hgraph_t* getHyperGraph()
      {
         hgraph_t* hg = new hgraph_t;
         hg->nvtxs = nbNodes;
         hg->nhedges = nbNodes;
         int* eptr = new int[nbNodes + 2];
         int* eind = new int[nbEdges + nbNodes + 2];
         int ind = 0;
         eptr[0] = 0;
         for(int i = 0; i < nbNodes; i++) {
            eind[ind++] = i;
            for(int j : successors[i])
               eind[ind++] = j;
            eptr[i+1] = ind;
         }
         hg->eptr = eptr;
         hg->eind = eind;
         return hg;
      }
};

vector<vector<int>> parseGraFile(FILE* f);

void output_dot(const char* filename, Graph g, vector<int> part);

pair<vector<int>, int> patoh_bisect(Graph graph);
//vector<int> mtmetis_bisect(Graph graph);

#endif
