#ifndef __DIJKSTRA_H__
#define __DIJKSTRA_H__

#ifdef __cplusplus
extern "C" {
#endif

typedef short vertex_t;
typedef short arc_t;
typedef float cost_t;
typedef short index_t;

typedef struct AdjacentList_{
  vertex_t V, *adjacent_vertices;
  arc_t A, *n_arcs;
  cost_t *costs;
} AdjacentList;


void malloc_adjl(AdjacentList *al);
void index_adjl(AdjacentList al, arc_t *index);
void free_adjl(AdjacentList *al);

void dijkstra ( AdjacentList adjl,		
		vertex_t u,
		char *vb,
		vertex_t nv,	
		vertex_t *heap,
		index_t *pos,
		cost_t *d,
		vertex_t *trace);

#ifdef __cplusplus
}
#endif

		
#endif
