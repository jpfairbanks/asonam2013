#include "stinger.h"
#include "stinger-shared.h"
#include "stinger-physmap-new.h"
#include "csv.h"
#include "xmalloc.h"
#include "timer.h"

/* BC */
#include "main.h"
#include "bcTreeDS.h"
#include "bcAlgorithms.h"
#include "dsUtils.h"
#include "dsAccum.h"
#include "timing_util.h"

/* GTRI-2 */
#include "clustering.h"

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
 * EDGE AND VERTEX TYPES
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#define ETYPES \
  TYPE(ETYPE_MENTION), \
  TYPE(MAX_ETYPE)

#define TYPE(X) X
typedef enum {
  ETYPES
} etype_t;
#undef TYPE

#define TYPE(X) #X
const char * etype_strings [] = {
  ETYPES
};
#undef TYPE

#define VTYPES \
  TYPE(VTYPE_NONE), \
  TYPE(VTYPE_USER), \
  TYPE(MAX_VTYPE)

#define TYPE(X) X
typedef enum {
  VTYPES
} vtype_t;
#undef TYPE

#define TYPE(X) #X
const char * vtype_strings [] = {
  VTYPES
};
#undef TYPE

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
 * Structs
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

typedef struct config {
  int64_t initial_size;
  int64_t batch_size;
  int64_t max_batches;
} config_t;

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
 * Function protos
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
void
parse_config(int argc, char ** argv, config_t * config);

void
histogram(struct stinger * S, uint64_t * labels, char * name, int64_t iteration);

void
histogram_float(struct stinger * S, float * scores, char * name, int64_t iteration);

void
REDUCTION(float_t** finalResultArray,float_t** parallelArray,int32_t threadCount);


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
 * MAIN
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
#define BC_ROOTS 256

int main(int argc, char *argv[]) {
  config_t config;
  parse_config(argc, argv, &config);

  struct stinger * S, * S_cluster;
  struct stinger_shared * S_shared, * S_cluster_shared;

  char * name = "/testing";
  S = stinger_shared_new(&S_shared, &name);
  S_cluster = stinger_shared_private(&S_cluster_shared, name);

  uint64_t threadCount = 1;
  #if defined(_OPENMP)
    OMP("omp parallel")
    {
      OMP("omp master")
      {
	threadCount = omp_get_num_threads();
      }
    }
  #endif

  /* set up data structures for algorithms */
  uint64_t * matches;
  uint64_t num_components;
  uint64_t * components = xmalloc(sizeof(uint64_t) * STINGER_MAX_LVERTICES);

  bcTree**    parallelForest    = createParallelForest(threadCount,STINGER_MAX_LVERTICES);
  float**     totalBCSS                 = createParallelBetweennessArray(threadCount,STINGER_MAX_LVERTICES);
  float**     finalBC           = createParallelBetweennessArray(1,STINGER_MAX_LVERTICES);
  uint32_t**  parallelSingleQueue = (uint32_t**)malloc(sizeof(uint32_t*)*threadCount);
  for(int i=0;i<threadCount;i++)
    parallelSingleQueue[i] = (uint32_t*)malloc(sizeof(uint32_t)*STINGER_MAX_LVERTICES);
  float timePerThread[threadCount];
  int rootsPerThread = BC_ROOTS/threadCount;
  uint32_t*** ppArray = createParallelParentArrayStinger(S ,STINGER_MAX_LVERTICES,threadCount);

  /* parsing structures */
  char      *   buf = NULL;
  uint64_t      bufSize = 0;
  char      **  fields = NULL;
  uint64_t  *   lengths = NULL;
  uint64_t      fieldsSize = 0;
  uint64_t      count = 0;
  int64_t	src, dest;


  /* begin parsing initial graph */
  {
    for(uint64_t e = 0; e < config.initial_size && !feof(stdin); e++) {
      readCSVLineDynamic(',', stdin, &buf, &bufSize, &fields, &lengths, &fieldsSize, &count);

      if(stinger_new_physID_to_vtx_mapping(S, fields[0], lengths[0], &src) == 1) {
	stinger_set_vtype(S, src, VTYPE_USER);
      }

      if(stinger_new_physID_to_vtx_mapping(S, fields[1], lengths[1], &dest) == 1) {
	stinger_set_vtype(S, dest, VTYPE_USER);
      }
      
      if(src != -1 && dest != -1) {
	stinger_insert_edge_pair(S, ETYPE_MENTION, src, dest, 1, e);
      }
    }
  }

  /* do first round of algorithms */
  printf("Clustering...\n");
  tic();
  static_multi_contract_clustering(&matches, STINGER_MAX_LVERTICES, S_cluster, S);
  printf("\tDone... %lf seconds \n", toc());
  histogram(S, matches, "communities", 0);
  free(matches);

  printf("Components...\n");
  tic();
  num_components = connected_components(S, components, STINGER_MAX_LVERTICES);
  printf("\tDone... %lf seconds \n", toc());
  histogram(S, components, "components", 0);

  {
    printf("\tPicking roots\n");
    tic();
    int64_t* selectedRoots=(int64_t*)malloc(sizeof(int64_t)*BC_ROOTS);
    int64_t r1 = 0;
    for(int64_t v = 0; v < STINGER_MAX_LVERTICES && r1 < BC_ROOTS; v++) {
      if(stinger_outdegree(S,v)) {
        selectedRoots[r1++] = v;
      }
    }
    printf("\tDone... %lf seconds \n", toc());

    printf("Starting BC...\n");
    tic();
    for(int v=0;v<STINGER_MAX_LVERTICES;v++)
      finalBC[0][v]=0;
    GO(bfsBrandes,STINGER,SINGLE,PA)(parallelForest,S,totalBCSS,(void**)parallelSingleQueue,
               ppArray,timePerThread,selectedRoots,rootsPerThread);
    REDUCTION(finalBC,totalBCSS,threadCount);
    printf("\tDone... %lf seconds \n", toc());
    histogram_float(S, finalBC[0], "bc", 0);
  }

  /* start batching */
  for(uint64_t b = 0; b < config.max_batches && !feof(stdin); b++) {
    printf("*** STARTING BATCH %ld ***\n", b+1);
    {
      for(uint64_t e = 0; e < config.batch_size && !feof(stdin); e++) {
	readCSVLineDynamic(',', stdin, &buf, &bufSize, &fields, &lengths, &fieldsSize, &count);

	if(stinger_new_physID_to_vtx_mapping(S, fields[0], lengths[0], &src) == 1) {
	  stinger_set_vtype(S, src, VTYPE_USER);
	}

	if(stinger_new_physID_to_vtx_mapping(S, fields[1], lengths[1], &dest) == 1) {
	  stinger_set_vtype(S, dest, VTYPE_USER);
	}
	
	if(src != -1 && dest != -1) {
	  stinger_insert_edge_pair(S, ETYPE_MENTION, src, dest, 1, e);
	}
      }
    }

    /* do algorithms */
    S_cluster = stinger_new();
    printf("Clustering...\n");
    tic();
    STINGER_PARALLEL_FORALL_EDGES_BEGIN(S, ETYPE_MENTION) {
      stinger_insert_edge(S_cluster, STINGER_EDGE_TYPE, STINGER_EDGE_SOURCE, STINGER_EDGE_DEST, STINGER_EDGE_WEIGHT, STINGER_EDGE_TIME_RECENT);
    } STINGER_PARALLEL_FORALL_EDGES_END();
    static_multi_contract_clustering(&matches, STINGER_MAX_LVERTICES, S_cluster, S);
    printf("\tDone... %lf seconds \n", toc());
    histogram(S, matches, "communities", b+1);
    free(matches);
    stinger_free_all(S_cluster);

    printf("Components...\n");
    tic();
    num_components = connected_components(S, components, STINGER_MAX_LVERTICES);
    printf("\tDone... %lf seconds \n", toc());
    histogram(S, components, "components", b+1);

    {
      printf("BC: Picking roots\n");
      tic();
      int64_t* selectedRoots=(int64_t*)malloc(sizeof(int64_t)*BC_ROOTS);
      int64_t r1 = 0;
      for(int64_t v = 0; v < STINGER_MAX_LVERTICES && r1 < BC_ROOTS; v++) {
	if(stinger_outdegree(S,v)) {
	  selectedRoots[r1++] = v;
	}
      }
      printf("\tDone... %lf seconds \n", toc());

      printf("BC: Starting ...\n");
      tic();
      for(int v=0;v<STINGER_MAX_LVERTICES;v++)
	finalBC[0][v]=0;
      GO(bfsBrandes,STINGER,SINGLE,PA)(parallelForest,S,totalBCSS,(void**)parallelSingleQueue,
		 ppArray,timePerThread,selectedRoots,rootsPerThread);
      REDUCTION(finalBC,totalBCSS,threadCount);
      printf("\tDone... %lf seconds \n", toc());
      histogram_float(S, finalBC[0], "bc", b+1);
    }
  }

}

void
parse_config(int argc, char ** argv, config_t * config) {
  config->initial_size = 1000;
  config->batch_size = 1000;
  config->max_batches = 1000;

  for(uint64_t i = 1; i < argc; i+=2) {
    char * str = argv[i];
    int len = strlen(argv[i]);
    if(!strncmp(str, "-", 1)) {
      str += 1; len -= 1;
      if(len) switch(*str) {
	case '-':
	  {
	    str++; len--;
	    if(len) switch(*str) {
	      case 'b':
		{
		  str++; len--;
		  if(!strncmp(str, "atchsize", len)) {
		    str += 8; len -= 8;
		    if(len == 0) {
		      /* --batchsize */
		      config->batch_size = atol(argv[i+1]);
		    }
		  }
		} break;
	      case 'i':
		{
		  str++; len--;
		  if(!strncmp(str, "nitialsize", len)) {
		    str += 10; len -= 10;
		    if(len == 0) {
		      /* --initialsize */
		      config->initial_size = atol(argv[i+1]);
		    }
		  }
		} break;
	      case 'n':
		{
		  str++; len--;
		  if(!strncmp(str, "umbatches", len)) {
		    str += 9; len -= 9;
		    if(len == 0) {
		      /* --numbatches */
		      config->max_batches = atol(argv[i+1]);
		    }
		  }
		} break;
	    }
	  } break;
	case 'b':
	  {
	    str++; len--;
	    if(len == 0) {
	      /* -b */
	      config->batch_size = atol(argv[i+1]);
	    }
	  } break;
	case 'i':
	  {
	    str++; len--;
	    if(len == 0) {
	      /* -i */
	      config->initial_size = atol(argv[i+1]);
	    }
	  } break;
	case 'n':
	  {
	    str++; len--;
	    if(len == 0) {
	      /* -n */
	      config->max_batches = atol(argv[i+1]);
	    }
	  } break;
      }
    }
  }
}

void
histogram(struct stinger * S, uint64_t * labels, char * name, int64_t iteration) {
  uint64_t * histogram = xcalloc(sizeof(uint64_t), STINGER_MAX_LVERTICES);
  uint64_t * sizes = xcalloc(sizeof(uint64_t), STINGER_MAX_LVERTICES);

  for(uint64_t v = 0; v < STINGER_MAX_LVERTICES; v++) {
    if(stinger_vtype(S, v) != VTYPE_NONE) {
      sizes[labels[v]]++;
    }
  }

  for(uint64_t v = 1; v < STINGER_MAX_LVERTICES; v++) {
    histogram[sizes[v]]++;
  }

  free(sizes);

  char filename[1024];
  sprintf(filename, "%s.%ld.csv", name, iteration);
  FILE * fp = fopen(filename, "w");
  for(uint64_t v = 1; v < STINGER_MAX_LVERTICES; v++) {
    if(histogram[v]) {
      fprintf(fp, "%ld, %ld\n", v, histogram[v]);
    }
  }
  fclose(fp);
  free(histogram);
}

void
histogram_float(struct stinger * S, float * scores, char * name, int64_t iteration) {
  int64_t max = 0;
  for(uint64_t v = 0; v < STINGER_MAX_LVERTICES; v++) {
    if(stinger_vtype(S, v) != VTYPE_NONE) {
      if((int64_t)(scores[v]) > max) {
	max = (int64_t) (scores[v]);
      }
    }
  }

  uint64_t * histogram = xcalloc(sizeof(uint64_t), (max+2));
  
  if(histogram) {
    for(uint64_t v = 0; v < STINGER_MAX_LVERTICES; v++) {
      histogram[(int64_t)(scores[v])]++;
    }

    char filename[1024];
    sprintf(filename, "%s.%ld.csv", name, iteration);
    FILE * fp = fopen(filename, "w");
    for(uint64_t v = 0; v < max+2; v++) {
      if(histogram[v]) {
	fprintf(fp, "%ld, %ld\n", v, histogram[v]);
      }
    }
    fclose(fp);
    free(histogram);
  } else {
    printf(" FAIL \n"); fflush(stdout);
  }
}

void REDUCTION(float_t** finalResultArray,float_t** parallelArray,int32_t threadCount) {
   #pragma parallel for
    for(int v=0;v<STINGER_MAX_LVERTICES;v++) {
        for(int t=0;t<threadCount;t++)
            finalResultArray[0][v]+=parallelArray[t][v];
    }
}
