#include "stinger.h"
#include "stinger-atomics.h"
#include "xmalloc.h"

int64_t
bfs(struct stinger *S, int64_t nv, int64_t source, int64_t * distance) {
  int64_t * queue = xmalloc(sizeof(int64_t) * nv);
  int64_t * found = xcalloc(nv, sizeof(int64_t));


  int64_t q_start = 0;
  int64_t q_end = 1;
  int64_t q_next = 1;

  queue[0] = source;
  found[source] = 1;

  if(distance) {

    OMP("omp parallel for")
    for(int64_t v = 0; v < nv; v++)
      distance[v] = INT64_MAX;

    distance[source] = 0;

    while(q_start != q_end) {
      OMP("omp parallel for")
      for(int64_t q = q_start; q < q_end; q++) {
	int64_t s = queue[q];
	int64_t s_d = distance[s];
	STINGER_FORALL_EDGES_OF_VTX_BEGIN(S, s) {
	  if(!found[STINGER_EDGE_DEST]) {
	    if(0 == stinger_int64_fetch_add(found + STINGER_EDGE_DEST, 1)) {
	      int64_t d = stinger_int64_fetch_add(&q_next, 1);
	      queue[d] = STINGER_EDGE_DEST;
	      distance[STINGER_EDGE_DEST] = s_d + 1;
	    }
	  }
	} STINGER_FORALL_EDGES_OF_VTX_END();
      }
      q_start = q_end;
      q_end = q_next;
    }
  } else {
    while(q_start != q_end) {
      OMP("omp parallel for")
      for(int64_t q = q_start; q < q_end; q++) {
	int64_t s = queue[q];
	STINGER_FORALL_EDGES_OF_VTX_BEGIN(S, s) {
	  if(!found[STINGER_EDGE_DEST]) {
	    if(0 == stinger_int64_fetch_add(found + STINGER_EDGE_DEST, 1)) {
	      int64_t d = stinger_int64_fetch_add(&q_next, 1);
	      queue[d] = STINGER_EDGE_DEST;
	    }
	  }
	} STINGER_FORALL_EDGES_OF_VTX_END();
      }
      q_start = q_end;
      q_end = q_next;
    }
  }

  free(queue);
  free(found);
}

int64_t
bfs_in_set(struct stinger *S, int64_t nv, int64_t source, int64_t * distance, int64_t * members, int64_t target) {
  int64_t * queue = xmalloc(sizeof(int64_t) * nv);
  int64_t * found = xcalloc(nv, sizeof(int64_t));


  int64_t q_start = 0;
  int64_t q_end = 1;
  int64_t q_next = 1;

  queue[0] = source;
  found[source] = 1;

  if(distance) {

    OMP("omp parallel for")
    for(int64_t v = 0; v < nv; v++)
      distance[v] = INT64_MAX;

    distance[source] = 0;

    while(q_start != q_end) {
      OMP("omp parallel for")
      for(int64_t q = q_start; q < q_end; q++) {
	int64_t s = queue[q];
	int64_t s_d = distance[s];
	STINGER_FORALL_EDGES_OF_VTX_BEGIN(S, s) {
	  if(members[STINGER_EDGE_DEST] == target) {
	    if(!found[STINGER_EDGE_DEST]) {
	      if(0 == stinger_int64_fetch_add(found + STINGER_EDGE_DEST, 1)) {
		int64_t d = stinger_int64_fetch_add(&q_next, 1);
		queue[d] = STINGER_EDGE_DEST;
		distance[STINGER_EDGE_DEST] = s_d + 1;
	      }
	    }
	  }
	} STINGER_FORALL_EDGES_OF_VTX_END();
      }
      q_start = q_end;
      q_end = q_next;
    }
  } else {
    while(q_start != q_end) {
      OMP("omp parallel for")
      for(int64_t q = q_start; q < q_end; q++) {
	int64_t s = queue[q];
	STINGER_FORALL_EDGES_OF_VTX_BEGIN(S, s) {
	  if(!found[STINGER_EDGE_DEST]) {
	    if(0 == stinger_int64_fetch_add(found + STINGER_EDGE_DEST, 1)) {
	      int64_t d = stinger_int64_fetch_add(&q_next, 1);
	      queue[d] = STINGER_EDGE_DEST;
	    }
	  }
	} STINGER_FORALL_EDGES_OF_VTX_END();
      }
      q_start = q_end;
      q_end = q_next;
    }
  }

  free(queue);
  free(found);
}
