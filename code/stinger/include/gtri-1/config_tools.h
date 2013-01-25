#define EDGE_BUF_SIZE 250000

#include "glib.h"
#include "stinger-utils.h"
#include "stinger-atomics.h"

typedef struct {
  uint64_t weight;
  uint64_t timefirst;
  uint64_t timerecent;
} local_edge_t;

#define trees(X,Y) trees[X*tree_nv+Y]
GTree ** trees;
uint64_t tree_nv = 0;
local_edge_t * edgq;
uint64_t edgq_size = STINGER_MAX_LVERTICES * 4;
uint64_t edgq_top = 0;
uint8_t * dest_ptr = (uint8_t *)0x1;

gint
ptr_cmp(gconstpointer a, gconstpointer b) {
  return a - b;
}

void
plant_trees(int64_t nv) {
  trees = calloc(sizeof(GTree *), MAX_ETYPE * nv);
  edgq = malloc(sizeof(local_edge_t) * edgq_size);
  tree_nv = nv;
}

inline void insert(int64_t type, int64_t source, int64_t dest, int64_t weight, int64_t time) {
  local_edge_t * result = NULL;
  if(!trees(type,source)) {
    trees(type,source) = g_tree_new(ptr_cmp);
  } else {
    result = g_tree_lookup(trees(type,source), dest + dest_ptr);
  }
  if(result) {
    result->weight += weight;
    if(time > result->timerecent)
      result->timerecent = time;
    else if (time < result->timefirst)
      result->timefirst = time;
  } else {
    if(edgq_size <= edgq_top) {
      edgq_size += STINGER_MAX_LVERTICES;
      edgq = realloc(edgq, sizeof(local_edge_t) * edgq_size);
    }
    local_edge_t * my_edge = edgq + edgq_top++;
    my_edge->weight = weight;
    my_edge->timefirst = time;
    my_edge->timerecent = time;
    g_tree_insert(trees(type,source), dest + dest_ptr, my_edge);
  }
}

typedef struct {
  int64_t * off, * dest, * weight, * timerecent, * timefirst;
} csr_t;
typedef struct {
  csr_t * csr;
  int64_t v;
} workpace_t;

gboolean tree_to_csr(uint8_t * dest, local_edge_t * local_edge, workpace_t * workspace) {
  int64_t my_offset =(workspace->csr->off[workspace->v])++;
  workspace->csr->dest[my_offset] = dest - dest_ptr;
  workspace->csr->weight[my_offset] = local_edge->weight;
  workspace->csr->timefirst[my_offset] = local_edge->timefirst;
  workspace->csr->timerecent[my_offset] = local_edge->timerecent;

  return 0;
}

void trees_to_stinger(struct stinger * S, struct stinger * S_cluster, uint64_t nv) {
  tic();
  csr_t csr;

  csr.off = xmalloc(sizeof(int64_t) * (nv+2));
  for(uint64_t cur_type = 0; cur_type < MAX_ETYPE; cur_type++) {
    csr.off[0] = 0; csr.off[1] = 0;
    OMP("omp parallel for")
    for(int64_t v = 0; v < nv; v++) {
      if(trees(cur_type,v)) {
	csr.off[v+2] = g_tree_nnodes(trees(cur_type,v));
      } else {
	csr.off[v+2] = 0;
      }
    }

    prefix_sum(nv+2, csr.off);

    csr.dest		= malloc(sizeof(int64_t) * csr.off[nv+1]*4);
    csr.weight		= &csr.dest[csr.off[nv+1]];
    csr.timerecent	= &csr.dest[csr.off[nv+1]*2];
    csr.timefirst	= &csr.dest[csr.off[nv+1]*3];

    OMP("omp parallel for")
    for(int64_t v = 0; v < nv; v++) {
      if(trees(cur_type,v)) {
	workpace_t workspace; 
	workspace.csr = &csr;
	workspace.v = v+1;
	g_tree_foreach(trees(cur_type,v),(GTraverseFunc)tree_to_csr,&workspace);
      } 
    }

    {
      struct stinger_vb *restrict LVA = S->LVA;
      struct stinger_eb *restrict ebpool_priv = S->ebpool->ebpool;

      int64_t * blkoff = xmalloc(sizeof(int64_t) * (nv+1));
      blkoff[0] = 0;
      OMP("omp parallel for")
      for(int64_t v = 0; v < nv; v++) {
	blkoff[v+1] = csr.off[v+1] - csr.off[v];
	stinger_int64_fetch_add(&(S->LVA[v].outDegree), blkoff[v+1]);
	blkoff[v+1] = (blkoff[v+1] + STINGER_EDGEBLOCKSIZE - 1) / (STINGER_EDGEBLOCKSIZE);
      }
      prefix_sum(nv+1, blkoff);

      eb_index_t * block = xcalloc (blkoff[nv], sizeof (*block));

      new_blk_ebs (block, S, nv, blkoff, cur_type);

      OMP ("omp parallel for")
      for (int64_t v = 0; v < nv; ++v) {
	const size_t nextoff = csr.off[v + 1];
	size_t kgraph = csr.off[v];
	int64_t from;

	from = v;

	for (size_t kblk = blkoff[v]; kblk < blkoff[v + 1]; ++kblk) {
	  size_t n_to_copy, voff;
	  struct stinger_edge *restrict edge;
	  struct stinger_eb *restrict eb;
	  int64_t tslb = INT64_MAX, tsub = 0;

	  voff = kgraph;
	  kgraph += STINGER_EDGEBLOCKSIZE;

	  n_to_copy = STINGER_EDGEBLOCKSIZE;
	  if (voff + n_to_copy >= nextoff)
	    n_to_copy = nextoff - voff;

	  eb = ebpool_priv + block[kblk];
	  edge = &eb->edges[0];

	  for (size_t i = 0; i < n_to_copy; ++i) {
	    const int64_t to = csr.dest[voff + i];
	    stinger_int64_fetch_add (&LVA[to].inDegree, 1);
	    edge[i].neighbor = to;
	    edge[i].weight = csr.weight[voff + i];
	    edge[i].timeRecent = csr.timerecent[voff + i];
	    edge[i].timeFirst = csr.timefirst[voff + i];
	  }

	  for (size_t i = 0; i < n_to_copy; ++i) {
	    if (edge[i].timeFirst < tslb)
	      tslb = edge[i].timeFirst;
	    if (edge[i].timeRecent > tsub)
	      tsub = edge[i].timeRecent;
	  }

	  eb->smallStamp = tslb;
	  eb->largeStamp = tsub;
	  eb->numEdges = n_to_copy;
	  eb->high = n_to_copy;
	}

	if (blkoff[v] != blkoff[v + 1]) {
	  ebpool_priv[block[blkoff[v+1]-1]].next = LVA[from].edges;
	  LVA[from].edges = block[blkoff[v]];
	}
      }

      push_ebs (S, blkoff[nv], block);

      free (block);
      free (blkoff);
    }
    if(S_cluster) {
      struct stinger_vb *restrict LVA = S_cluster->LVA;
      struct stinger_eb *restrict ebpool_priv = S_cluster->ebpool->ebpool;

      int64_t * blkoff = xmalloc(sizeof(int64_t) * (nv+1));
      blkoff[0] = 0;
      OMP("omp parallel for")
      for(int64_t v = 0; v < nv; v++) {
	blkoff[v+1] = csr.off[v+1] - csr.off[v];
	stinger_int64_fetch_add(&(S_cluster->LVA[v].outDegree), blkoff[v+1]);
	blkoff[v+1] = (blkoff[v+1] + STINGER_EDGEBLOCKSIZE - 1) / (STINGER_EDGEBLOCKSIZE);
      }
      prefix_sum(nv+1, blkoff);

      eb_index_t * block = xcalloc (blkoff[nv], sizeof (*block));

      new_blk_ebs (block, S_cluster, nv, blkoff, cur_type);

      OMP ("omp parallel for")
      for (int64_t v = 0; v < nv; ++v) {
	const size_t nextoff = csr.off[v + 1];
	size_t kgraph = csr.off[v];
	int64_t from;

	from = v;

	for (size_t kblk = blkoff[v]; kblk < blkoff[v + 1]; ++kblk) {
	  size_t n_to_copy, voff;
	  struct stinger_edge *restrict edge;
	  struct stinger_eb *restrict eb;
	  int64_t tslb = INT64_MAX, tsub = 0;

	  voff = kgraph;
	  kgraph += STINGER_EDGEBLOCKSIZE;

	  n_to_copy = STINGER_EDGEBLOCKSIZE;
	  if (voff + n_to_copy >= nextoff)
	    n_to_copy = nextoff - voff;

	  eb = ebpool_priv + block[kblk];
	  edge = &eb->edges[0];

	  for (size_t i = 0; i < n_to_copy; ++i) {
	    const int64_t to = csr.dest[voff + i];
	    stinger_int64_fetch_add (&LVA[to].inDegree, 1);
	    edge[i].neighbor = to;
	    edge[i].weight = csr.weight[voff + i];
	    edge[i].timeRecent = csr.timerecent[voff + i];
	    edge[i].timeFirst = csr.timefirst[voff + i];
	  }

	  for (size_t i = 0; i < n_to_copy; ++i) {
	    if (edge[i].timeFirst < tslb)
	      tslb = edge[i].timeFirst;
	    if (edge[i].timeRecent > tsub)
	      tsub = edge[i].timeRecent;
	  }

	  eb->smallStamp = tslb;
	  eb->largeStamp = tsub;
	  eb->numEdges = n_to_copy;
	  eb->high = n_to_copy;
	}

	if (blkoff[v] != blkoff[v + 1]) {
	  ebpool_priv[block[blkoff[v+1]-1]].next = LVA[from].edges;
	  LVA[from].edges = block[blkoff[v]];
	}
      }

      push_ebs (S_cluster, blkoff[nv], block);

      free (block);
      free (blkoff);
    }

    
    free(csr.dest); 
  }

  free(csr.off);
  free(edgq);
  printf("\ttree to stinger %lf seconds\n", toc());
}

/* reverting to simpler time format (cuts parsing timestamps by > 10x) */
inline uint64_t parse_time(char * time) {
  uint64_t out = 0;
  out += 10000000000000 * (time[0] - '0')  + 
	 1000000000000  * (time[1] - '0')  + 
	 100000000000   * (time[2] - '0')  + 
         10000000000    * (time[3] - '0')  +
         1000000000     * (time[5] - '0')  + 
         100000000      * (time[6] - '0')  +
         10000000       * (time[8] - '0')  + 
         1000000        * (time[9] - '0')  +
         100000         * (time[11] - '0') + 
         10000          * (time[12] - '0') +
         1000           * (time[14] - '0') + 
         100            * (time[15] - '0') +
         10             * (time[17] - '0') + 
			  (time[18] - '0');
  return out;
}

#define PARSE void

#define BEGIN                             \
 (struct stinger * S, const char * dir, char delim) { \
  char filename[2048];                    \
                                          \
  char      *   buf = NULL;               \
  uint64_t      bufSize = 0;              \
  char      **  fields = NULL;            \
  uint64_t  *   lengths = NULL;           \
  uint64_t      fieldsSize = 0;           \
  uint64_t      count = 0;                \
                                          \
  sprintf(filename, "%s%s", dir,

#define FILENAME

#define FIELDS                            \
  ); enum {

#define ENDFIELDS                         \
  };                                      \
                                          \
  FILE * fp = fopen(filename, "r");       \
  if(NULL == fp) {                        \
    fprintf(stderr,"failed to open %s file\n", __func__); \
    return;                               \
  } tic();                                      \
  uint64_t results = 0;                   \
  while(!feof(fp)) {                      \
    readCSVLineDynamic(delim, fp, &buf, &bufSize, &fields, &lengths, &fieldsSize, &count); \
    if(count < 2 && lengths[0] < 2) break; \
    results++;

#define MAPVTX(X, F, T) \
      int64_t X;                                                                             \
      if(1 == stinger_new_physID_to_vtx_mapping(S, fields[F], lengths[F], &X)) { \
        stinger_set_vtype(S, X, T);                                           \
      }

#define PARSETIME(F,N)                             \
  uint64_t N = parse_time(fields[F]);

#define EDGE(TYPE,SRC,DEST,T)                              \
  insert(TYPE, SRC, DEST, 1, T);\
  insert(TYPE, DEST, SRC, 1, T);

#define END(X)                                       \
  }                                               \
  free(buf); free(fields); free(lengths); \
  printf("\tRead %ld %s rows in %lf seconds\n", results, #X, toc()); fflush(stdout); \
}
