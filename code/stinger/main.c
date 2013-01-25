#define _XOPEN_SOURCE

#include <time.h>
#include <stdint.h>
#include <glib.h>
#if defined(_OPENMP)
  #include <omp.h>
#endif
#include <unistd.h>

/* GTRI-1 */
#include "stinger.h"
#include "stinger-physmap-new.h"
#include "stinger-shared.h"
#include "csv.h"
#include "timer.h"
#include "xmalloc.h"
#include "alg_convert.h"

/* BC */
#include "main.h"
#include "bcTreeDS.h"
#include "bcAlgorithms.h"
#include "dsUtils.h"
#include "dsAccum.h"
#include "timing_util.h"

/* GTRI-2 */
#include "clustering.h"

/* GTRI-3 */
#include "seed_set_exp.h"

/* GTRI-6 */
#include "evicomb.h"
#include "evicomb-internal.h"

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
 * EDGE AND VERTEX TYPES
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#define ETYPES \
  TYPE(ETYPE_EMAIL_ALIAS), \
  TYPE(ETYPE_EMAIL_FROM), \
  TYPE(ETYPE_EMAIL_TO), \
  TYPE(ETYPE_ACCESSED), \
  TYPE(ETYPE_LOGON), \
  TYPE(ETYPE_LOGOFF), \
  TYPE(ETYPE_PRINT), \
  TYPE(ETYPE_IM_FROM), \
  TYPE(ETYPE_IM_TO), \
  TYPE(ETYPE_IM_ALIAS), \
  TYPE(ETYPE_CLIPBOARD), \
  TYPE(ETYPE_REGISTRY), \
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
  TYPE(VTYPE_DEVICE), \
  TYPE(VTYPE_WORKSTATION), \
  TYPE(VTYPE_FILE), \
  TYPE(VTYPE_EMAIL_ADDR), \
  TYPE(VTYPE_EMAIL), \
  TYPE(VTYPE_DOMAIN), \
  TYPE(VTYPE_URL), \
  TYPE(VTYPE_PATH), \
  TYPE(VTYPE_PRINTER), \
  TYPE(VTYPE_APPLICATION), \
  TYPE(VTYPE_IM), \
  TYPE(VTYPE_IMID), \
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
 * Algorithm and graph types
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#define ATYPES \
  TYPE(ATYPE_GTRI2), \
  TYPE(ATYPE_GTRI3), \
  TYPE(ATYPE_GTRIBC), \
  TYPE(ATYPE_GTRI6), \
  TYPE(MAX_ATYPE)

#define TYPE(X) X
typedef enum {
  ATYPES
} atype_t;
#undef TYPE

#define TYPE(X) #X
const char * atype_strings [] = {
  ATYPES
};
#undef TYPE

#define GTYPES \
  TYPE(GTYPE_EMAIL_EVENT), \
  TYPE(GTYPE_FILE_EVENT), \
  TYPE(GTYPE_IM_EVENT), \
  TYPE(GTYPE_LOGON), \
  TYPE(GTYPE_PRINTER_EVENT), \
  TYPE(GTYPE_PROCESS_EVENT), \
  TYPE(GTYPE_REGISTRY_EVENT), \
  TYPE(GTYPE_URL_EVENT), \
  TYPE(GTYPE_USER_EMAIL), \
  TYPE(GTYPE_USER_IM), \
  TYPE(MAX_GTYPE)

#define TYPE(X) X
typedef enum {
  GTYPES
} gtype_t;
#undef TYPE

#define TYPE(X) #X
const char * gtype_strings [] = {
  GTYPES
};
#undef TYPE

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
 * Global Settings
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
#define GLOBAL_DATE_FORMAT "%Y-%m/%d %H:%M:%S"

#define BC_ROOTS 256

#if STINGER_MAX_LVERTICES < 1ULL<<24
  #warning "You probably need more vertices"
#endif

#define HERE() printf("%s %d\n", __func__, __LINE__); fflush(stdout);

typedef struct {
  uint8_t is_client;
  uint8_t is_server;
  uint8_t algorithms[MAX_ATYPE];
  uint8_t graph[MAX_GTYPE];
} global_config;

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
 * Load config
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
#include "config_tools.h"
#include "config.h"

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
 * Function Protos
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
struct stinger * build_stinger_from_csv(char * csvDir, char delim, struct stinger ** S_cluster, global_config * conf);
struct stinger * build_shared_stinger_from_csv(char * csvDir, char delim, struct stinger_shared ** shared, char ** name, global_config * conf);
int parse_config(int argc, char ** argv, global_config * conf);
void print_config(global_config * conf, char ** argv, char * prefix);
void print_usage();

void REDUCTION(float_t** finalResultArray,float_t** parallelArray,int32_t threadCount);

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
 * MAIN
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
int
main(int argc, char *argv[])
{
  /* if so, go adjust stinger-defs.h */
  assert(MAX_ETYPE <= STINGER_NUMETYPES);

  if(argc < 3) {
    printf("usage: %s input_directory_path output_directory_path\n", argv[0]);
    print_usage();
    exit(-1);
  }

  global_config conf;
  parse_config(argc, argv, &conf);

  printf("Up and running. Config...\n");
  print_config(&conf, argv, "\t");

  /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
   * GTRI-1 Build graph, gather stats
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
  struct stinger * S, * S_cluster;
  struct stinger_shared * S_shared, * S_cluster_shared;
  char * name = argv[2];
  if(!conf.is_client && !conf.is_server) {
    printf("Building Graph...\n");
    S = build_stinger_from_csv(argv[1], ',', &S_cluster, &conf);
  } else if(conf.is_client) {
    printf("Mapping graph...\n");
    S = stinger_shared_map(&S_shared, argv[1]);
    S_cluster = stinger_shared_private(&S_cluster_shared, argv[1]);
  } else if (conf.is_server) {
    printf("Building and mapping graph...\n");
    S = build_shared_stinger_from_csv(argv[1], ',', &S_shared, &name, &conf);
  }

  printf("Graph(s) created. Running stats...");
  tic();
  printf("\n\tVertices: %ld\n\tEdges: %ld\n",
    stinger_num_active_vertices(S), stinger_total_edges(S));

  for(uint64_t t = 0; t < STINGER_NUMETYPES; t++) {
    printf("\tEBlocks[%s]: %ld\n", etype_strings[t], S->ETA[t].high);
  }

  printf("\tConsistency %ld\n", stinger_consistency_check(S, STINGER_MAX_LVERTICES));
  printf("\tDone. %lf seconds\n", toc());

  if(conf.is_server) {
    printf("Server pausing....\n");
    printf("\tgraph mapping:\n%s", name);
    fflush(stdout);
    pause();
    printf("\nServer unpaused. Shutting down...\n"); fflush(stdout);
    stinger_shared_free(S, S_shared, name);
    exit(0);
  }


  /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
   * Betweenness Centrality
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
  if(conf.algorithms[ATYPE_GTRIBC]) {
    printf("Running Betweenness Centrality...\n");
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

    printf("\tBuilding trees\n");
    bcTree**    parallelForest    = createParallelForest(threadCount,STINGER_MAX_LVERTICES);
    float**     totalBCSS                 = createParallelBetweennessArray(threadCount,STINGER_MAX_LVERTICES);
    float**     finalBC           = createParallelBetweennessArray(1,STINGER_MAX_LVERTICES);

    uint32_t**  parallelSingleQueue = (uint32_t**)malloc(sizeof(uint32_t*)*threadCount);
    for(int i=0;i<threadCount;i++)
      parallelSingleQueue[i] = (uint32_t*)malloc(sizeof(uint32_t)*STINGER_MAX_LVERTICES);


    printf("\tPicking roots\n");
    int64_t* selectedRoots=(int64_t*)malloc(sizeof(int64_t)*BC_ROOTS);
    int64_t r1 = 0;
    for(int64_t v = 0; v < STINGER_MAX_LVERTICES && r1 < BC_ROOTS; v++) {
      if(stinger_outdegree(S,v)) {
        selectedRoots[r1++] = v;
      }
    }

    printf("\tAllocating trees...\n");
    float timePerThread[threadCount];

    int rootsPerThread = BC_ROOTS/threadCount;
    uint32_t*** ppArray = createParallelParentArrayStinger(S ,STINGER_MAX_LVERTICES,threadCount);
    printf("\tStarting...\n");
    tic();
    for(int v=0;v<STINGER_MAX_LVERTICES;v++)
      finalBC[0][v]=0;
    GO(bfsBrandes,STINGER,SINGLE,PA)(parallelForest,S,totalBCSS,(void**)parallelSingleQueue,
               ppArray,timePerThread,selectedRoots,rootsPerThread);
    REDUCTION(finalBC,totalBCSS,threadCount);
    printf("\tDone... %lf seconds \n", toc());
    destroyParallelParentArray(ppArray,STINGER_MAX_LVERTICES,threadCount);

    char bcfilename[1024];
    sprintf(bcfilename, "%s/bc.csv", argv[2]);
    FILE * bcfile = fopen(bcfilename, "w");
    csvIfIDExistsfloat(bcfile, ',', S, vtype_strings, STINGER_MAX_LVERTICES, finalBC[0]);
    fclose(bcfile);

    free(selectedRoots);

    for(int i=0;i<threadCount;i++)
      free(parallelSingleQueue[i]);

    free(parallelSingleQueue);

    destroyParallelBetweennessArray(totalBCSS,threadCount);
    destroyParallelBetweennessArray(finalBC,1);
    destroyParallelForest(parallelForest,threadCount);
  }

  /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
   * GTRI-2 Community Detection
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
  if(conf.algorithms[ATYPE_GTRI2]) {
    printf("Beginning clustering...\n");
    uint64_t * matches;
    static_multi_contract_clustering(&matches, STINGER_MAX_LVERTICES, S_cluster, S);
    char communitiesfilename[1024];
    sprintf(communitiesfilename, "%s/communities.csv", argv[2]);
    FILE * communitiesfile = fopen(communitiesfilename, "w");
    csvIfIDExistsint64(communitiesfile, ',', S, vtype_strings, STINGER_MAX_LVERTICES, matches);
    fclose(communitiesfile);
    free(matches);

    printf("Beginning connected components...\n");
    uint64_t * components = xmalloc(sizeof(uint64_t) * STINGER_MAX_LVERTICES);
    tic();
    uint64_t num_components = connected_components(S, components, STINGER_MAX_LVERTICES);
    printf("\tDone. %lf seconds %ld components\n", toc(), num_components);
    char componentsfilename[1024];
    sprintf(componentsfilename, "%s/components.csv", argv[2]);
    FILE * componentsfile = fopen(componentsfilename, "w");
    csvIfIDExistsint64(componentsfile, ',', S, vtype_strings, STINGER_MAX_LVERTICES, components);
    fclose(componentsfile);
    free(components);
  }

  /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
   * GTRI-3 Seedset Expansion
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
  if(conf.algorithms[ATYPE_GTRI3]) {
    printf("Beginning seedset expansion...\n");
    int64_t num_seeds = 0;
    int64_t * seeds       = (int64_t*) xmalloc(STINGER_MAX_LVERTICES*sizeof(int64_t));

    printf("\tReading seeds...\n");
    char seedsfilename[1024];
    sprintf(seedsfilename, "%s/seeds.csv", argv[1]);
    FILE * seedsfile = fopen(seedsfilename, "r");
    char      *   buf = NULL;
    uint64_t      bufSize = 0;
    char      **  fields = NULL;
    uint64_t  *   lengths = NULL;
    uint64_t      fieldsSize = 0;
    uint64_t      count = 0;
    while(!feof(seedsfile)) {
      readCSVLineDynamic(',', seedsfile, &buf, &bufSize, &fields, &lengths, &fieldsSize, &count);
      if(lengths[0] < 1) break;
      int64_t v = stinger_physID_to_vtx(S, fields[0], lengths[0]);
      if(v >= 0) {
        seeds[num_seeds++] = v;
      }
    }
    fclose(seedsfile);
    free(buf); free(fields); free(lengths);
    int64_t * membership          = (int64_t*) xmalloc(STINGER_MAX_LVERTICES*sizeof(int64_t));
    int64_t * neighbors   = (int64_t*) xmalloc(STINGER_MAX_LVERTICES*sizeof(int64_t));
    int64_t * edge_between  = (int64_t*) xmalloc(STINGER_MAX_LVERTICES*sizeof(int64_t));
    int64_t * affected    = (int64_t*) xmalloc(STINGER_MAX_LVERTICES*sizeof(int64_t));

    double * neighbor_modval= (double*)xmalloc(STINGER_MAX_LVERTICES*sizeof(double));

    double history[HISTORY];
    int64_t history_idx = 0;
    int64_t total_comm_edges;
    int64_t internal_comm_edges;
    int64_t number_neighbors;

    printf("\tStarting...\n");
    tic();
    seed_set_expansion(S, seeds, num_seeds, STINGER_MAX_LVERTICES, stinger_total_edges(S), membership,
                       neighbors, edge_between, neighbor_modval, &total_comm_edges, &internal_comm_edges,
                       &number_neighbors,history, &history_idx);
    printf("\tDone... %lf seconds\n", toc());

    char setmembersfilename[1024];
    sprintf(setmembersfilename, "%s/setmembers.csv", argv[2]);
    FILE * setmembersfile = fopen(setmembersfilename, "w");
    csvIfIDExistsint64(setmembersfile, ',', S, vtype_strings, STINGER_MAX_LVERTICES, membership);
    fclose(setmembersfile);

    free(seeds); free(membership); free(neighbors); free(edge_between); free(affected); free(neighbor_modval);
  }

  /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
   * GTRI-6
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
  //this is disabled because it doesn't work
  if(conf.algorithms[ATYPE_GTRI6]) & 0 {

    printf("Beginning belief propogation...\n");

    printf("\tReading params..\n");
    char params[100][2048];
    FILE * paramsfile = fopen("files/params.txt", "r");

    uint64_t paramcount = 0;
    for(; paramcount < 100 && !feof(paramsfile); paramcount++) {
      fgets(params[paramcount], 2048, paramsfile);
      if(strlen(params[paramcount]) < 1) {
        paramcount--;
        continue;
      }
    }

    fclose(paramsfile);

    uint64_t algoIDs[] = {103,1150,1151,1154,1160,1161,1170,1171,1172,1173,1180,1181,1190};
    uint64_t sizes[]   = {100,500,1000,5000};

    int * anomalyIDs = malloc(sizeof(uint64_t) * 5000);

    for(uint64_t size = 0; size < 4; size++) {
      for(uint64_t p = 0; p < paramcount; p++) {
        char filename[2048];
        sprintf(filename, "/tmp/prodigalscores/avg/%ld/all", sizes[size]);
        printf("\tReading values from %s...\n", filename);

        int anomalycount = 0;
        FILE * list = fopen(filename, "r");
        for(; anomalycount < sizes[size] && !feof(list); anomalycount++) {
          char vertex[2048];
          fgets(vertex, 2048, list);
          if(strlen(vertex) < 1) {
            anomalycount--;
            continue;
          }
          anomalyIDs[anomalycount] = stinger_physID_to_vtx(S, vertex, strlen(vertex)-1);
          if(anomalyIDs[anomalycount] < 0)
            anomalycount--;
        }
        fclose(list);

        /* perform the belief propogation */
        printf("\tStarting... %ld %ld avg all\n", sizes[size], p);
        tic();
        evicomb(S, anomalycount, anomalyIDs, params[p], "avg_");
        printf("\tDone... %lf seconds\n", toc());
      }
    }

    for(uint64_t size = 0; size < 4; size++) {
      for(uint64_t p = 0; p < paramcount; p++) {
        char filename[2048];
        sprintf(filename, "/tmp/prodigalscores/max/%ld/all", sizes[size]);
        printf("\tReading values from %s...\n", filename);

        int anomalycount = 0;
        FILE * list = fopen(filename, "r");
        for(; anomalycount < sizes[size] && !feof(list); anomalycount++) {
          char vertex[2048];
          fgets(vertex, 2048, list);
          if(strlen(vertex) < 1) {
            anomalycount--;
            continue;
          }
          anomalyIDs[anomalycount] = stinger_physID_to_vtx(S, vertex, strlen(vertex));
          if(anomalyIDs[anomalycount] < 0)
            anomalycount--;
        }
        fclose(list);

        /* perform the belief propagation */
        printf("\tStarting... %ld %ld max all\n", sizes[size], p);
        tic();
        evicomb(S, anomalycount, anomalyIDs, params[p], "max_");
        printf("\tDone... %lf seconds\n", toc());
      }
    }
  }

  printf("Algorithms have completed. Closing.\n"); fflush(stdout);
  if(conf.is_client) {
    stinger_shared_unmap(S, S_shared, argv[1]);
  } else {
    stinger_free_all(S);
  }
}

int
parse_config(int argc, char ** argv, global_config * conf) {
  /* default config - not running client/server,
   * use all algorithms, all input tables*/
  conf->is_client = 0;
  conf->is_server = 0;
  for(uint64_t a = 0; a < MAX_ATYPE; a++) {
    conf->algorithms[a] = 1;
  }
  for(uint64_t g = 0; g < MAX_GTYPE; g++) {
    conf->graph[g] = 1;
  }

  argc -= 2;
  argv += 2;

  char      *   buf = NULL;
  uint64_t      bufSize = 0;
  char      **  fields = NULL;
  uint64_t  *   lengths = NULL;
  uint64_t      fieldsSize = 0;
  uint64_t      count = 0;

  int opt;
  while((opt = getopt(argc, argv, "a:g:cs?")) != -1) {
    switch (opt) {
      case 'a': {
        for(uint64_t a = 0; a < MAX_ATYPE; a++) {
          conf->algorithms[a] = 0;
        }
        splitLineCSVDynamic(':', optarg, strlen(optarg), &buf, &bufSize, &fields, &lengths, &fieldsSize, &count);
        for(uint64_t a = 0; a < count; a++) {
          if(lengths[a] > 0) {
            switch (*fields[a]) {
              case 'B':
              case 'b': {
                conf->algorithms[ATYPE_GTRIBC] = 1;
              } break;

              case 'c':
              case 'C':
              case '2': {
                conf->algorithms[ATYPE_GTRI2] = 1;
              } break;

              case 's':
              case 'S':
              case '3': {
                conf->algorithms[ATYPE_GTRI3] = 1;
              } break;

              case 'e':
              case 'E':
              case '6': {
                conf->algorithms[ATYPE_GTRI6] = 1;
              } break;

              default: {
                fprintf(stderr, "%s %d: Unrecognized algorithm: %s\n", __func__, __LINE__, fields[a]);
              } break;
            }
          }
        }
      } break;

      case 'g': {
        for(uint64_t g = 0; g < MAX_GTYPE; g++) {
          conf->graph[g] = 0;
        }
        splitLineCSVDynamic(':', optarg, strlen(optarg), &buf, &bufSize, &fields, &lengths, &fieldsSize, &count);
        for(uint64_t g = 0; g < count; g++) {
          if(lengths[g] > 0) {
            switch (*fields[g]) {
              case 'e':
              case 'E': {
                conf->graph[GTYPE_EMAIL_EVENT] = 1;
              } break;

              case 'f':
              case 'F': {
                conf->graph[GTYPE_FILE_EVENT] = 1;
              } break;

              case 'i':
              case 'I': {
                conf->graph[GTYPE_IM_EVENT] = 1;
              } break;

              case 'l':
              case 'L': {
                conf->graph[GTYPE_LOGON] = 1;
              } break;

              case 'p':
              case 'P': {
                if(fields[g][2] == 'i' || fields[g][1] == 'I') {
                  conf->graph[GTYPE_PRINTER_EVENT] = 1;
                } else if (fields[g][2] == 'o' || fields[g][2] == 'O') {
                  conf->graph[GTYPE_PROCESS_EVENT] = 1;
                } else {
                  fprintf(stderr, "%s %d: Unrecognized table: %s\n", __func__, __LINE__, fields[g]);
                }
              } break;

              case 'r':
              case 'R': {
                conf->graph[GTYPE_REGISTRY_EVENT] = 1;
              } break;

              case 'u':
              case 'U': {
                if(fields[g][1] == 'r' || fields[g][1] == 'R') {
                  conf->graph[GTYPE_URL_EVENT] = 1;
                } else if (fields[g][4] == 'e' || fields[g][4] == 'E') {
                  conf->graph[GTYPE_USER_EMAIL] = 1;

                } else if (fields[g][4] == 'i' || fields[g][4] == 'I') {
                  conf->graph[GTYPE_USER_IM] = 1;
                } else {
                  fprintf(stderr, "%s %d: Unrecognized table: %s\n", __func__, __LINE__, fields[g]);
                }
              } break;

              default: {
                fprintf(stderr, "%s %d: Unrecognized table: %s\n", __func__, __LINE__, fields[g]);
              } break;
            }
          }
        }
      } break;

      case 'c': {
        conf->is_client = 1;
        conf->is_server = 0;
      } break;

      case 's': {
        conf->is_client = 0;
        conf->is_server = 1;
      } break;

      default: {
        printf("Unrecognized option: %c\n", opt);
        print_usage();
        exit(-1);
      } break;
    }
  }

  if(buf)     free(buf);
  if(fields)  free(fields);
  if(lengths) free(lengths);
}

void
print_config(global_config * conf, char ** argv, char * prefix) {
  uint8_t algorithms[MAX_ATYPE];
  uint8_t graph[MAX_GTYPE];
  printf("%sreading-from: %s\n", prefix, argv[1]);
  printf("%swriting-to: %s\n", prefix, argv[2]);

  printf("%sclient-server: ", prefix);
  if(conf->is_client) {
    printf("client\n");
  } else if (conf->is_server) {
    printf("server\n");
  } else {
    printf("stand-alone\n");
  }

  printf("%sgraph-tables-enabled:\n", prefix);
  for(uint64_t g = 0; g < MAX_GTYPE; g++) {
    if(conf->graph[g]) {
      printf("%s\t%s\n", prefix, gtype_strings[g]);
    }
  }

  printf("%salgorithms-enabled:\n", prefix);
  for(uint64_t a = 0; a < MAX_ATYPE; a++) {
    if(conf->algorithms[a]) {
      printf("%s\t%s\n", prefix, atype_strings[a]);
    }
  }
}


struct stinger * build_stinger_from_csv(char * csvDir, char delim, struct stinger ** S_cluster, global_config * conf) {
  struct stinger * S = stinger_new();
  *S_cluster = stinger_new();

  plant_trees(STINGER_MAX_LVERTICES);

  if(conf->graph[GTYPE_EMAIL_EVENT])
    mentions(S, csvDir, delim);

  trees_to_stinger(S, *S_cluster, STINGER_MAX_LVERTICES);

  return S;
}

struct stinger * build_shared_stinger_from_csv(char * csvDir, char delim, struct stinger_shared ** shared, char ** name, global_config * conf) {
  struct stinger * S = stinger_shared_new (shared, name);
  printf("\n\tallocated %s\n", *name); fflush(stdout);

  plant_trees(STINGER_MAX_LVERTICES);

  if(conf->graph[GTYPE_EMAIL_EVENT])
    mentions (S, csvDir, delim);

  trees_to_stinger(S, NULL, STINGER_MAX_LVERTICES);

  return S;
}

void
print_usage() {
printf(
"STINGER for ADAMS\n"
"=================\n"
"\n"
"* Usage *\n"
"  Stand-alone usage:\n"
"  ./main /path/to/input/directory/ /path/to/output/ [-a algs:strings:here] [-g table:names:here]\n"
"\n"
"  Server usage:\n"
"  ./main /path/to/input/direcotry/ /servermappingname -s [-g table:names:here]\n"
"\n"
"  Client usage:\n"
"  ./main /servermappingname /path/to/output/directory/ -c [-a algs:strings:here]\n"
"\n"
"* Algorithms *\n"
"  GTRI-BC Betweenness Centrality\n"
"  Accepted strings: b B betweennessCentrality betweenness [any string starting with b or B]\n"
"\n"
"  GTRI-2 Community Detection / Clustering\n"
"  Accepted strings: 2 c C clustering community [any string starting with c or C]\n"
"\n"
"  GTRI-3 Seed Set Expansion\n"
"  Accepted strings: 3 s S seedset [any string starting with s or S]\n"
"\n"
"  GTRI-6 Evidence combination\n"
"  Accepted strings: 6 e E evidence [any string starting with e or E]\n"
"\n"
"* Accepted Table Names *\n"
"  emailEvent \n"
"  file \n"
"  imEvent \n"
"  logon \n"
"  printerEvent \n"
"  processEvent \n"
"  registryEvent \n"
"  urlEvent \n"
"  userEmail \n"
"  userIM \n"
);
}


void REDUCTION(float_t** finalResultArray,float_t** parallelArray,int32_t threadCount) {
   #pragma parallel for
    for(int v=0;v<STINGER_MAX_LVERTICES;v++) {
        for(int t=0;t<threadCount;t++)
            finalResultArray[0][v]+=parallelArray[t][v];
    }
}
