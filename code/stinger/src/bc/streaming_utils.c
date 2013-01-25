#include "streaming_utils.h"

#include <alloca.h>
#include <stdio.h>
#include <stdlib.h>

#include "stinger.h"
#include "stinger-traversal.h"
#include "stinger-utils.h"
#include "xmalloc.h"
#include "csr.h"

//#include "graph_defs.h"
//#include "graph_gen.h"
//#include "utils.h"

uint32_t ** squareMat(uint32_t NV) {
	uint32_t ** rtn = xmalloc(sizeof(uint32_t *) * NV);
	for(uint32_t i = 0; i < NV; i++) {
		rtn[i] = calloc(NV, sizeof(uint32_t));
	}
	return rtn;
}

void freeSquareMat(uint32_t ** sq, uint32_t NV) {
	for(uint32_t i = 0; i < NV; i++) {
		free(sq[i]);
	}
	free(sq);
}

bc_t ** squareMatFloat(uint32_t NV) {
	bc_t ** rtn = xmalloc(sizeof(bc_t *) * NV);
	for(uint32_t i = 0; i < NV; i++) {
		rtn[i] = calloc(NV, sizeof(bc_t));
	}
	return rtn;
}

void freeSquareMatFloat(bc_t ** sq, uint32_t NV) {
	for(uint32_t i = 0; i < NV; i++) {
		free(sq[i]);
	}
	free(sq);
}

void intersectAdjacencyWithLevel(uint32_t * outArray, uint32_t * numFound, uint32_t NV,
	uint32_t * adjacency, uint32_t * levelArray, uint32_t level) {
	for(uint32_t i = 0; i < NV; i++) {
		if(adjacency[i] > 0 && levelArray[i] == level) {
			adjacency[(*numFound)++] = i;
		}
	}
}

void moveUpTree(uint32_t * levels, uint32_t * pathsToRoot, uint32_t ** graph,
	uint32_t vertex, uint32_t ignoreVertex, uint32_t dist, uint32_t NV) {

	uint32_t touched[NV];
	uint32_t delta[NV];
	uint32_t BFSQueue[NV];
	uint32_t BFSQueueCur, BFSQueueNext;
	BFSQueueCur = BFSQueueNext = 0;

	BFSQueue[BFSQueueNext] = vertex;
	touched[vertex] = 1;
	pathsToRoot[vertex] = pathsToRoot[ignoreVertex];
	delta[BFSQueueNext++] = dist;
	touched[ignoreVertex] = NV+1;

	// while queue not empty
	while(BFSQueueCur != BFSQueueNext) {
		// traverse current queue head's adjacency list
		for(uint32_t i = 0; i < NV; i++) {
			// if adjacent and untouched
			if(graph[BFSQueue[BFSQueueCur]][i] > 0) {
				if(touched[i] < 0) {
					// compute distance the adjacent vertex should be moved
					int32_t computedDelta = levels[BFSQueue[BFSQueueCur]] - levels[i] - 1 + delta[BFSQueue[BFSQueueCur]];
					touched[i] = touched[BFSQueue[BFSQueueCur]] + 1;
					// if the adjacent vertex should be moved, put it in the queue
					if(computedDelta > 0) {
						pathsToRoot[i] = pathsToRoot[BFSQueue[BFSQueueCur]];
						delta[BFSQueueNext] = computedDelta;
						BFSQueue[BFSQueueNext++] = i;
					} else if(computedDelta == 0) {
						pathsToRoot[i] += pathsToRoot[BFSQueue[BFSQueueCur]];
					}
				} else if(touched[i] == touched[BFSQueue[BFSQueueCur]]) {
					pathsToRoot[i] += pathsToRoot[BFSQueue[BFSQueueCur]];
				}
			}
		}
		// move ourself and retire
		levels[BFSQueue[BFSQueueCur]] -= delta[BFSQueue[BFSQueueCur]];
		BFSQueueCur++;
	}
}

void addEdgeWithoutMovement(uint32_t * levels, uint32_t * pathsToRoot, uint32_t ** graph,
	uint32_t vertex, uint32_t deltaTop, uint32_t NV) {

	uint32_t touched[NV];
	uint32_t BFSQueue[NV];
	uint32_t delta[NV];
	uint32_t BFSQueueCur, BFSQueueNext;
	BFSQueueCur = BFSQueueNext = 0;

	BFSQueue[BFSQueueNext] = vertex;
	touched[vertex] = 1;
	delta[vertex] = deltaTop;
	pathsToRoot[vertex] += deltaTop;

	// while queue not empty
	while(BFSQueueCur != BFSQueueNext) {
	    	    printf("fuck");
		for(uint32_t i = 0; i < NV; i++) {
			if(graph[BFSQueue[BFSQueueCur]][i] > 0 && levels[BFSQueue[BFSQueueCur]] == levels[i] + 1) {
				if(touched[i] == 0) {
					BFSQueue[BFSQueueNext++] = i;
					touched[i] = touched[BFSQueue[BFSQueueCur]]+1;
					pathsToRoot[i] += delta[BFSQueue[BFSQueueCur]];
					delta[i] = delta[BFSQueue[BFSQueueCur]];
				} else if(touched[i] == touched[BFSQueue[BFSQueueCur]] + 1) {
				    printf("SHIT");
					pathsToRoot[i] += delta[BFSQueue[BFSQueueCur]];
					delta[i] += delta[BFSQueue[BFSQueueCur]];
				}
			}
		}
		BFSQueueCur++;
	}
}

void copyTreeIn(uint32_t * levelsDest, uint32_t * pathsToRootDest, uint32_t existingVertex, uint32_t newVertex, uint32_t *levelsSrc, uint32_t * pathsToRootSrc, uint32_t INF, uint32_t NV) {
	uint32_t pathsToRootMultiplier = pathsToRootDest[existingVertex];
	uint32_t levelsIncrease = levelsDest[existingVertex] + 1;
	for(uint32_t i = 0; i < NV; i++) {
		if(levelsSrc[i] != INF) {
			levelsDest[i] = levelsSrc[i] + levelsIncrease;
			pathsToRootDest[i] = pathsToRootSrc[i] * pathsToRootMultiplier;
		}
	}
}



void printGraph(uint32_t ** graph, uint32_t NV) {
	for(uint32_t i = 0; i < NV; i++) {
		for(uint32_t j = 0; j < NV; j++) {
			printf("%d", graph[i][j]);
		}
		printf("\n");
	}
}

void MatrixBCComput(uint32_t NV, uint32_t * pathsToRoot[],uint32_t * level[],bc_t* perTreeBC[],bc_t* BC)
{
    DB("Calculating BCs...\n")
	for(uint32_t i = 0; i < NV; i++) {
		VB("\tBC[%d]...\n",i)
		BC[i] = 0.0f;
		for(uint32_t j = 0; j < NV; j++) {
			for(uint32_t k = 0; k < NV; k++) {
				if(j == k || i == j || i == k)
					continue;
				if(d(j,i) + d(i,k) == d(j,k)) {
					perTreeBC[i][j] += (sigmaf(j,i) * sigmaf(i,k) / sigmaf(j,k));
				}
			}
			BC[i] += perTreeBC[i][j];
			VB("\t\tperTreeBC[%d][%d] = %f\n", i,j,perTreeBC[i][j]);
		}
		INFO("\tBC[%d] = %f\n", i, BC[i]);
	}
}


int64_t serial_shiloach_vishkin_componentsCSR (csrGraph* graph, int64_t * component_map)
{
  /* Initialize each vertex with its own component label in parallel */
  for (uint64_t i = 0; i < graph->NV; i++) {
    component_map[i] = i;
  }

  /* Iterate until no changes occur */
  while (1) {
    int changed = 0;

        for(int v=0; v<graph->NV;v++)
        {
            int src = v;
            for(int edge=graph->vertexPointerArray[src]; edge<graph->vertexPointerArray[src+1];edge++)
            {
                int dest = graph->edgeArray[edge];
                if (component_map[dest] < component_map[src])
                {
                    component_map[src] = component_map[dest];
                    changed++;
                }

            }
        }
    /* For all edges in the STINGER graph of type 0 in parallel, attempt to assign
       lesser component IDs to neighbors with greater component IDs */
/*    STINGER_FORALL_EDGES_BEGIN (S, 0) {
      if (component_map[STINGER_EDGE_DEST] <
          component_map[STINGER_EDGE_SOURCE]) {
        component_map[STINGER_EDGE_SOURCE] = component_map[STINGER_EDGE_DEST];
        changed++;
      }
    }
    STINGER_FORALL_EDGES_END ();
*/
    /* if nothing changed */
    if (!changed)
      break;

    /* Tree climbing with OpenMP parallel for */
    for (uint64_t i = 0; i < graph->NV; i++) {
      while (component_map[i] != component_map[component_map[i]])
        component_map[i] = component_map[component_map[i]];
    }
  }

  /* Count components */
  uint64_t components = 1;
  for (uint64_t i = 1; i < graph->NV; i++) {
    if (component_map[i] == i) {
      components++;
    }
  }

  return components;
}


int64_t serial_shiloach_vishkin_components (struct stinger *S, int64_t nv, int64_t * component_map)
{
  /* Initialize each vertex with its own component label in parallel */
  for (uint64_t i = 0; i < nv; i++) {
    component_map[i] = i;
  }

  /* Iterate until no changes occur */
  while (1) {
    int changed = 0;

    /* For all edges in the STINGER graph of type 0 in parallel, attempt to assign
       lesser component IDs to neighbors with greater component IDs */
    STINGER_FORALL_EDGES_BEGIN (S, 0) {
      if (component_map[STINGER_EDGE_DEST] <
          component_map[STINGER_EDGE_SOURCE]) {
        component_map[STINGER_EDGE_SOURCE] = component_map[STINGER_EDGE_DEST];
        changed++;
      }
    }
    STINGER_FORALL_EDGES_END ();

    /* if nothing changed */
    if (!changed)
      break;

    /* Tree climbing with OpenMP parallel for */
    for (uint64_t i = 0; i < nv; i++) {
      while (component_map[i] != component_map[component_map[i]])
        component_map[i] = component_map[component_map[i]];
    }
  }

  /* Count components */
  uint64_t components = 1;
  for (uint64_t i = 1; i < nv; i++) {
    if (component_map[i] == i) {
      components++;
    }
  }

  return components;
}


void selectRootsInMaxComponents(csrGraph* graph,int NV,int64_t* component_map, int64_t componentCount,
                                int numberOfRoots,int* countedEdges, int64_t* selectedRoots)
{
    int64_t* componentCounter=(int64_t*)malloc(sizeof(int64_t)*(graph->NV));
    int64_t* componentSelection=(int64_t*)malloc(sizeof(int64_t)*(graph->NV));

	// Initializing
    for (int v=0; v<NV; v++)
    {
         componentCounter[v]=0;
         componentSelection[v]=0;
    }

    // Creating a histogram for the connect components
    for (int v=0; v<NV; v++)
    {
        componentCounter[component_map[v]]++;

//		printf("(%d,", component_map[v]);
//		printf("%d),",  graph->vertexPointerArray[v+1]-graph->vertexPointerArray[v]);
	}

    int maxComp=0; int max=componentCounter[0];

    // Looking for component of the biggest size
    for (int v=0; v<NV; v++)
    {
        if (componentCounter[v]>max)
        {
            max=componentCounter[v];
            maxComp=0;
        }
    }
//printf("\nbiggest component size id : %d %d \n", max, maxComp) ;
    // Selecting the roots
    for(int r=0; r<numberOfRoots;r++)
    {
        int root = rand()%NV;
        if(componentSelection[root]==0 && component_map[root]==maxComp)
        {
            componentSelection[root]=1;
            selectedRoots[r]=root;
        }
        else
        {
            r--;
        }
    }

    *countedEdges=0;
    // Counting the edges in the graph
    for (int v=0; v<NV; v++)
    {
        if(component_map[v]==maxComp)
            *countedEdges+=graph->vertexPointerArray[v+1]-graph->vertexPointerArray[v];
    }
   free(componentCounter); free(componentSelection);
}


void readSnapGraph(struct stinger** G, char* filename)
{

	int NV = 5242 ,NE=28980;
    int count;


	FILE* inputFile = fopen(filename,"r");

	fscanf(inputFile,"%d\n%d\n",&NV,&NE);

    int64_t* edgeSrc = (int64_t*)xmalloc((NE+1)*sizeof(int64_t));
    int64_t* edgeDest = (int64_t*)xmalloc((NE+1)*sizeof(int64_t));
    int64_t* edgeWeight = (int64_t*)xmalloc((NE+1)*sizeof(int64_t));
	if(edgeSrc == NULL || edgeDest==NULL || edgeWeight==NULL)
		printf("NULL Pointer in random graph creation\n");

    for (int64_t i = 0; i < NE; i++)
    {
        edgeSrc[i] = 0;
        edgeDest[i] = 0;
        edgeWeight[i] = 0;
    }

	//printf("before loop %d %d\n",NV,NE);

	count=0;
	while(count<NE)
	{
		fscanf(inputFile,"%ld %ld\n", &edgeSrc[count], &edgeDest[count]);
//		printf("read %d %ld %ld %ld\n", count, edgeSrc[count], edgeDest[count], edgeWeight[count]);
		count++;
	}
//    printf("after loop %d\n",count);
//    printf("last value %ld %ld\n", edgeSrc[count-1], edgeDest[count-1]);


    fclose(inputFile);

//	printf("after close\n");

    *G = edge_list_to_stinger(NV,NE,edgeSrc,edgeDest,edgeWeight,NULL,NULL,0);
//    printf("Stinger created\n");


    free(edgeSrc); free(edgeDest); free(edgeWeight);

}

