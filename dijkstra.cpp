#include <stdio.h>
#include <stdlib.h>
#include "dijkstra.h"

#define int(i) (int)(i)
#define FOR(i, n) for((i)=0; (i)<(n); (i)++)
#define MALLOC(type, size) ((type*) malloc((size)*sizeof(type)))
#define FREE(var)	  \
	do { \
		if(var) free(var); \
		(var) = NULL; \
	} while(0)

void malloc_adjl(AdjacentList *al){
	(*al).n_arcs = MALLOC(arc_t, (*al).V + 1);
	(*al).adjacent_vertices = MALLOC(vertex_t, (*al).A);
	(*al).costs = MALLOC(cost_t, (*al).A);
}

void index_adjl(AdjacentList al, arc_t *index){
	vertex_t i, j, V = al.V;
	arc_t a, A = al.A;
	FOR(i, V) FOR(j, V) index[int(i)*V + j] = -1;
	FOR(i, V) for(a = al.n_arcs[i]; a<al.n_arcs[i+1]; a++)
		index[int(i)*V + al.adjacent_vertices[a]] = a;
}

void free_adjl(AdjacentList *al){
	FREE((*al).n_arcs);
	FREE((*al).adjacent_vertices);
	FREE((*al).costs);
}

void dijkstra ( AdjacentList adjl,		
                vertex_t u,
                char *vb,
                vertex_t nv,
		
                vertex_t *heap,
                index_t *pos,
                cost_t *d,
		
                vertex_t *trace ) 
{
	vertex_t i, dadvertex, childvertex, childvertex2;
	index_t dad, child, heapsize = 1;
	cost_t c, dadcost = 0.0, childcost = 0.0, childcost2 = 0.0;

	// Initialisation
	for(i = 0; i < adjl.V; i++) pos[i] = -1, trace[i] = -1;  
	d[u] = 0; heap[0] = u; pos[u] = 0;
  
	while(nv>0 && heapsize>0) {
		c = d[u = heap[0]];
		if(vb[u]) if((--nv) == 0) return;

		//////////////////////////////////////////////////////////////////
		// pop heap (upheap)
		//
		pos[ heap[dad = 0] = (dadvertex = heap[--heapsize]) ] = 0;
		dadcost = d[dadvertex];
		while( (child=dad*2+1) < heapsize-1){
			childcost  = d[childvertex  = heap[child]];
			childcost2 = d[childvertex2 = heap[child+1]];
			if(childcost2 < childcost) 
				child++, childvertex = childvertex2, childcost = childcost2;
			if(childcost > dadcost) break;
			heap[dad] = childvertex; pos[childvertex] = dad;
			dad = child;
		}
    
		if(child == heapsize-1)
			if(d[childvertex = heap[child]] < dadcost){
				heap[dad] = heap[child];
				pos[childvertex] = dad;
				dad = child;
			}    
		heap[dad] = dadvertex; pos[dadvertex] = dad;
		/////// end pop heap /////////////////////////////////////////////

		// loop over all adjacent vertices of the new vertex
		for(i = adjl.n_arcs[u]; i < adjl.n_arcs[u+1]; i++){
			childvertex = adjl.adjacent_vertices[i];
			childcost = c + adjl.costs[i];

			// check new vertex and/or cost reduction
			if(pos[childvertex] < 0){
				// this is a new vertex, insert to the end of the heap
				heap[child = heapsize] = childvertex;
				pos[childvertex] = heapsize++;
			} 
			else if (childcost >= d[childvertex]) continue; // no cost reduction
			else child = pos[childvertex]; // cost reduction
      
			// update cost and trace
			trace[childvertex] = u;
			d[childvertex] = childcost;
      
			//////////////////////////////////////////////////////////////////
			// update heap (upheap)
			//
			while(child > 0) {
				dadcost = d[dadvertex = heap[dad = (child-1)/2]];
				if(dadcost < childcost) break;
				heap[child] = dadvertex; pos[dadvertex] = child;
				child = dad;
			}
			heap[child] = childvertex; pos[childvertex] = child;
			/////// end upheap ///////////////////////////////////////////////
		}
	}

	//FREE(heap);
	//FREE(pos);
	//FREE(d);
}
