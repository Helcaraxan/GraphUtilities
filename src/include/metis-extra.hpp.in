#ifndef METIS_EXTRA_HPP
#define METIS_EXTRA_HPP

#include <mtmetis.h>


// mt-METIS renaming macros
#define vtx_t mtmetis_vtx_t
#define adj_t mtmetis_adj_t


// Declarations for METIS compatibility
extern "C" {
  typedef struct graph_t {
    mtmetis_vtx_t nvtxs;	/* The # of vertices and edges in the graph */
    mtmetis_adj_t * xadj;		/* Pointers to the locally stored vertices */
    mtmetis_vtx_t * adjncy;        /* Array that stores the adjacency lists of nvtxs */
  } graph_internal_t;

  struct Tree{
    Tree * succ[2];
    int id;

    Tree(int id = -1, Tree* s1 = NULL, Tree* s2 = NULL) : id(id) {
      succ[0] = s1;
      succ[1] = s2;
    }

    ~Tree() {
      for(int i = 0; i < 2; i++)
        delete succ[i];
    }
  };
}

void print_tree(Tree * t);

/* partitions hierachically a graph with mt-metis */ 
Tree * hier_partition(mtmetis_vtx_t nvtxs, mtmetis_vtx_t * ids, mtmetis_adj_t * xadj,
    mtmetis_vtx_t * adjncy, double * opts);

/* builds the partition tree with mt-metis and returns a table with path in the tree
   encoded on 64 bits */
uint64_t * paths_in_tree(graph_t &graph);

#endif // METIS_EXTRA_HPP
