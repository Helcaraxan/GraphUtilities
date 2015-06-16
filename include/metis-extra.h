#include <metis.h>

// Declarations for METIS compatibility
extern "C" {
  typedef int64_t idx_t;
  typedef double real_t;

  typedef struct graph_t {
    idx_t nvtxs;	/* The # of vertices and edges in the graph */
    idx_t *xadj;		/* Pointers to the locally stored vertices */
    idx_t *adjncy;        /* Array that stores the adjacency lists of nvtxs */
  } graph_t;
}
