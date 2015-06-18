#include <cstdio>
#include <algorithm>
#include <iostream>
#include <string>

#include <mtmetis.h>

#include "metis-extra.hpp"

using namespace std;


// For debugging purpose
void print_tree_p(Tree * t, string p0, string p1) {
  if (t->id != -1) {
    cout << p0 << t->id << '\n' << p1 << '\n';
  } else {
    if (t->succ[0]) {
      string p00 = p0 + "---";
      string p01 = p0 + (t->succ[1] ? "|  " : "   "); 
      print_tree_p(t->succ[0], p00, p01);
    }

    if (t->succ[1]) {
      string p10 = p1 + "|--";
      string p11 = p1 + "   ";
      print_tree_p(t->succ[1], p10, p11);
    }
  }
}


void print_tree(Tree * t) {
  print_tree_p(t, "", "");
}


Tree* hier_partition(mtmetis_vtx_t nvtxs, mtmetis_vtx_t * ids, mtmetis_adj_t * xadj,
		mtmetis_vtx_t * adjncy, double * opts) {
  Tree * ret = new Tree();

  /*
     printf("DEBUG :\n nvxts : %d\n", nvtxs);
     for(int i=0; i<nvtxs; i++) {
     printf("%d : ", ids[i]);
     for(int j=xadj[i]; j<xadj[i+1]; j++)
     printf("%d ", ids[adjncy[j]]);
     printf("\n");
     }
     */


  if (nvtxs == 1)
    ret->id = ids[0];

  if (nvtxs == 2) {
    ret->succ[0] = new Tree(ids[0]);
    ret->succ[1] = new Tree(ids[1]);
  }

  if (nvtxs > 2) {
    mtmetis_pid_t * part = new mtmetis_pid_t[nvtxs];

    int err = mtmetis_partition_explicit(nvtxs, xadj, adjncy, NULL, NULL, opts, part, NULL);
    if(err != MTMETIS_SUCCESS) {
      printf("mt-metis error during partitioning: %s\n",
			  err == MTMETIS_ERROR_INVALIDINPUT ? "Invalid input" : "Not enough memory");
      exit(EXIT_FAILURE);
    }

    /*
       for(int i=0; i<nvtxs; i++)
       printf("%d : %d\n", ids[i], part[i]);
       */

    mtmetis_vtx_t nvtxs_s[2] = {0, 0};
    mtmetis_adj_t nedges_s[2] = {0, 0};
    mtmetis_vtx_t * ids_b = new mtmetis_vtx_t[nvtxs];

    // Finding new ids for both subgraphs;
    for (mtmetis_vtx_t i = 0; i < nvtxs; i++) {
      ids_b[i] = nvtxs_s[part[i]]++;
      for (mtmetis_adj_t j = xadj[i]; j < xadj[i+1]; j++)
        nedges_s[part[i]] += part[i] == part[adjncy[j]];
    }

    mtmetis_vtx_t * ids_s[2] = {new mtmetis_vtx_t[nvtxs_s[0]], new mtmetis_vtx_t[nvtxs_s[1]]};
    mtmetis_adj_t * xadj_s[2] = {new mtmetis_adj_t[nvtxs_s[0]+1], new mtmetis_adj_t[nvtxs_s[1]+1]};
    mtmetis_vtx_t * adjncy_s[2] = {new mtmetis_vtx_t[nedges_s[0]], new mtmetis_vtx_t[nedges_s[1]]};

    // Filling the new arrays
    int i_p[2] = {0, 0};
    int iadj_p[2] = {0, 0};
    xadj_s[0][0] = 0;
    xadj_s[1][0] = 0;
    for (mtmetis_vtx_t i = 0; i < nvtxs; i++) {
      ids_s[part[i]][i_p[part[i]]] = i;
      for (mtmetis_adj_t j = xadj[i]; j < xadj[i+1]; j++)
        if (part[i] == part[adjncy[j]]) 
          adjncy_s[part[i]][iadj_p[part[i]]++] = ids_b[adjncy[j]];

      xadj_s[part[i]][i_p[part[i]]+1] = iadj_p[part[i]];
      i_p[part[i]]++;
    }

    delete[] ids_b;
    delete[] part;

    // Recursive calls
    Tree * s0 = hier_partition(nvtxs_s[0], ids_s[0], xadj_s[0], adjncy_s[0], opts);
    Tree * s1 = hier_partition(nvtxs_s[1], ids_s[1], xadj_s[1], adjncy_s[1], opts);

    if (s0 == NULL) 
      swap(s0, s1);

    ret->succ[0] = s0;
    ret->succ[1] = s1;

    for (int i = 0; i < 2; i++) {
      delete[] ids_s[i];
      delete[] xadj_s[i];
      delete[] adjncy_s[i];
    }
  }

  return ret;
}

void set_t_path(Tree * tree, int depth, uint64_t path, uint64_t * t_path) {
  if (tree->id != -1) {
    t_path[tree->id] = path;
  }else {
    if (tree->succ[0])
      set_t_path(tree->succ[0], depth + 1, path, t_path);

    if (tree->succ[1])
      set_t_path(tree->succ[1], depth + 1, path + (UINT64_C(1) << (63 - depth)), t_path);
  }
}

uint64_t * paths_in_tree(graph_t &graph) {
  uint64_t * t_path = new uint64_t[graph.nvtxs];
  mtmetis_vtx_t * ids = new mtmetis_vtx_t[graph.nvtxs]; 

  double * opts = mtmetis_init_options();
  opts[MTMETIS_OPTION_NPARTS] = 2;
  opts[MTMETIS_OPTION_VERBOSITY] = 0;
  opts[MTMETIS_OPTION_NTHREADS] = 1;
  Tree * tree = hier_partition(graph.nvtxs, ids, graph.xadj, graph.adjncy, opts);

  set_t_path(tree, 0, 0, t_path);

  delete tree;
  delete[] ids;

  return t_path;
}
