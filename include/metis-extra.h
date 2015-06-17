#include <mtmetis.h>

// Declarations for METIS compatibility
extern "C" {
  typedef struct graph_t {
    mtmetis_vtx_t nvtxs;	/* The # of vertices and edges in the graph */
    mtmetis_adj_t *xadj;		/* Pointers to the locally stored vertices */
    mtmetis_vtx_t *adjncy;        /* Array that stores the adjacency lists of nvtxs */
  } graph_t;
}
