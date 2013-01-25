#include "omp.h"

#include "timing_util.h"
#include "bcTreeDS.h"

/**
* @brief Computes betweenness centrality based on Brandes's algorithm.
*
* @param tree The tree that needs to be unallocated.
* @param someGraph The tree that needs to be unallocated.
* @param currRoot The vertex in the graph that will be used as the root of the BFS traversal.
* @param totalBC The BC score for each vertex (array) .
* @param someQueue The queue that will be used by the algorithm. It is possible to use a regular queue or a multi-level queue.
* @param parentParam The data structure that will be used to maintain the parents that are in the BFS traversal.
* @param someDataStructure
*
* @return None
*/

void PRE_NAME(bfsBrandesPerTree)(bcTree* tree, void* someGraph,uint32_t currRoot,float* totalBC,
                                     void* someQueue, void* parentParam, void* someDataStrucutre, void* someDataStrucutre2)
{
    #if ACCUM_TYPE==ACCUM_PARENT_ARRAY
    uint32_t** parentArray = (uint32_t **)parentParam;
    int32_t*   parentCounter = (int32_t*)someDataStrucutre2;
    for(uint32_t j = 0; j < tree->NV; j++)
    {
        parentCounter[j] = 0;
    }
    #elif ACCUM_TYPE==ACCUM_PARENT_LIST
    list_ptr* parentList = (list_ptr*)parentParam;
    for(uint32_t j = 0; j < tree->NV; j++)
    {
        emptyList(parentList[j]);
    }

    #endif

    #if FOR_TYPE==FOR_ADJ
    uint32_t** matGraph = (uint32_t**)someGraph;
    #elif FOR_TYPE==FOR_STINGER
    struct stinger* sStinger = (struct stinger*) someGraph;
    #elif FOR_TYPE==FOR_CSR
    csrGraph* graph = (csrGraph*)someGraph;
    #endif

    for(uint32_t j = 0; j < tree->NV; j++)
    {
        tree->level[j] = INFINITY_MY;
        tree->pathsToRoot[j] = INFINITY_MY;
        tree->delta[j] = 0;
    }
	tree->level[currRoot] = 0;
	tree->pathsToRoot[currRoot] = 1;

    #if QUEUE_TYPE==QUEUE_SINGLE
//    uint32_t Queue[tree->NV];
    uint32_t* Queue = (uint32_t*)someQueue;
    Queue[0] = currRoot;
    int32_t qStart=0,qEnd=1;

    // While queue is not empty
    while(qStart!=qEnd)
    {
        uint32_t currElement;
        currElement = Queue[qStart];
        qStart++;


    #elif QUEUE_TYPE==QUEUE_MULTI_ON_SINGLE
    uint32_t* Queue = (uint32_t*)someQueue;
    Queue[0] = currRoot;
    int32_t qStart=0,qEnd=1;
    int32_t* queueLevel = (int32_t*)someDataStrucutre;
    queueLevel[0]=0;
    int32_t deepestLevel = 0;
    int32_t currentLevel = 0;

    // While queue is not empty
    while(qStart!=qEnd)
    {
        uint32_t currElement;
        currElement = Queue[qStart];
        qStart++;



    #elif QUEUE_TYPE==QUEUE_MULTI

    uint32_t** queueBFSTREE = (uint32_t**)someQueue;
    int32_t levelCounter[tree->NV];

    for(uint32_t j = 0; j < tree->NV; j++)
    {
        levelCounter[j]=0;
    }
    int32_t deepestLevel = 0;
    int32_t currentLevel = 0;
    int32_t levelIterator = 0;
    queueBFSTREE[0][levelCounter[0]++]=currRoot;

    int stam;

//      printf("%d,%d,\n",levelCounter[0],levelIterator);
//    scanf("%d",&stam);

    while(levelCounter[currentLevel]!=0)
    {
        if(levelIterator==levelCounter[currentLevel])
        {
//            printf("\n");
            currentLevel++;
            levelIterator =0;
            continue;
        }
//    printf("%d,%d,%d,\n",currentLevel,levelCounter[currentLevel],levelIterator);
//    scanf("%d",&stam);


        uint32_t currElement = queueBFSTREE[currentLevel][levelIterator++];

    #endif


        #if FOR_TYPE==FOR_ADJ
        // Checking all the neighbors of the current elements
        for(uint32_t k = 0; k < tree->NV; k++)
        {
            if(matGraph[currElement][k] == 0)
                continue;
        #elif FOR_TYPE==FOR_STINGER
        STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger,currElement)
        {
            uint32_t k = STINGER_EDGE_DEST;
        #elif FOR_TYPE==FOR_CSR

        int startEdge = graph->vertexPointerArray[currElement];
        int stopEdge = graph->vertexPointerArray[currElement+1];

        for(int j=startEdge;startEdge<stopEdge;startEdge++)
        {
            uint32_t k = graph->edgeArray[startEdge];

        #endif

            // If this is a neighbor and has not been found
            if(tree->level[k] > tree->level[currElement])
            {
                // Checking if "k" has been found.
                if(tree->level[k]==INFINITY_MY)
                {
                    tree->level[k] = tree->level[currElement]+1;
                    tree->delta[k]=0;

                    #if QUEUE_TYPE==QUEUE_SINGLE
                    Queue[qEnd++] = k;

                    #elif QUEUE_TYPE==QUEUE_MULTI_ON_SINGLE
                    Queue[qEnd++] = k;

				    if(deepestLevel<tree->level[k])
				    {
                        deepestLevel=tree->level[k];
                        queueLevel[deepestLevel]=qStart;
				    }


                    #elif QUEUE_TYPE==QUEUE_MULTI

                    queueBFSTREE[tree->level[k]][levelCounter[tree->level[k]]++]=k;
				    // Checking if a "deeper level" has been reached.
				    if(deepestLevel<tree->level[k])
                        deepestLevel=tree->level[k];

                    #endif

                }

                if(tree->pathsToRoot[k] == INFINITY_MY)
                {
                    // k has not been found and therefore its paths to the roots are through its parent.
                    tree->pathsToRoot[k] = tree->pathsToRoot[currElement];

                    #if ACCUM_TYPE==ACCUM_PARENT_ARRAY
                    parentArray[k][parentCounter[k]++] = currElement;
                    #elif ACCUM_TYPE==ACCUM_PARENT_LIST
                    append(parentList[k],makeNode(currElement));
                    #endif
                }
                else
                {
                    // k has been found and has multiple paths to the root as it has multiple parents.
                    tree->pathsToRoot[k] += tree->pathsToRoot[currElement];

                    #if ACCUM_TYPE==ACCUM_PARENT_ARRAY
                    parentArray[k][parentCounter[k]++] = currElement;
                    #elif ACCUM_TYPE==ACCUM_PARENT_LIST
                    append(parentList[k],makeNode(currElement));
                    #endif

                }
            }

        #if FOR_TYPE==FOR_ADJ
        }
        #elif FOR_TYPE==FOR_STINGER
        }
        STINGER_FORALL_EDGES_OF_VTX_END();
        #elif FOR_TYPE==FOR_CSR
        }
        #endif
    }


    // Using Brandes algorithm to compute BC for a specific tree.
    // Essentially, use the stack which the elements are placed in depth-reverse order, to "climb" back
    // up the tree, all the way to the root.



    #if QUEUE_TYPE==QUEUE_SINGLE
    qEnd = qStart-1;

    while(qEnd>=0)
    {
        uint32_t currElement = Queue[qEnd];

    #elif QUEUE_TYPE==QUEUE_MULTI_ON_SINGLE

    queueLevel[deepestLevel+1]=qEnd;
    int32_t levelStart, levelEnd;
    levelStart=queueLevel[deepestLevel];
    levelEnd=queueLevel[deepestLevel+1];

    while(deepestLevel>=0)
    {
        if(levelStart==levelEnd)
        {
            deepestLevel--;
            levelStart=queueLevel[deepestLevel];
            levelEnd=queueLevel[deepestLevel+1];
            continue;
        }

        uint32_t currElement = Queue[levelStart];
        levelStart++;



    #elif QUEUE_TYPE==QUEUE_MULTI
//    printf("I got here %d %d\n",deepestLevel,levelCounter[deepestLevel]);

    while(deepestLevel >=0 && levelCounter[deepestLevel] > 0)
    {

        int32_t currQueue = 0;

        while(currQueue<levelCounter[deepestLevel])
        {
            // Removing last element from the queue
            uint32_t currElement = queueBFSTREE[deepestLevel][currQueue++];

    #endif

            #if ACCUM_TYPE==ACCUM_TRAVERSAL

                 #if FOR_TYPE==FOR_ADJ
                // Checking all the neighbors of the current elements
                for(uint32_t k = 0; k < tree->NV; k++)
                {
                    if(matGraph[currElement][k] == 0)
                        continue;

                #elif FOR_TYPE==FOR_STINGER

                STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger,currElement)
                {
                    uint32_t k = STINGER_EDGE_DEST;

                #elif FOR_TYPE==FOR_CSR

                int startEdge = graph->vertexPointerArray[currElement];
                int stopEdge = graph->vertexPointerArray[currElement+1];

                for(int j=startEdge;startEdge<stopEdge;startEdge++)
                {
                    uint32_t k = graph->edgeArray[startEdge];
                #endif

                    // If this is a neighbor and has not been found
                    if((tree->level[k] == (tree->level[currElement]-1)))
                    {
            #elif ACCUM_TYPE==ACCUM_PARENT_LIST
                    node_t* iterator = parentList[currElement]->head;
                    while(iterator!=NULL)
                    {
                        uint32_t k = iterator->id;
                        iterator = iterator->next;


            #elif ACCUM_TYPE==ACCUM_PARENT_ARRAY
                    for(int j=0; j<parentCounter[currElement];j++)
                    {
                          uint32_t k=parentArray[currElement][j];
            #endif
                           tree->delta[k] +=
                                ((bc_t)tree->pathsToRoot[k]/(bc_t)tree->pathsToRoot[currElement])*
                                (bc_t)(tree->delta[currElement]+1);


            #if ACCUM_TYPE==ACCUM_TRAVERSAL
                    }

                #if FOR_TYPE==FOR_ADJ
                }
                #elif FOR_TYPE==FOR_STINGER
                }
                STINGER_FORALL_EDGES_OF_VTX_END();
                #elif FOR_TYPE==FOR_CSR
                }
                #endif


            #elif ACCUM_TYPE==ACCUM_PARENT_LIST
                    }
            #elif ACCUM_TYPE==ACCUM_PARENT_ARRAY
                    }
            #endif

            if(currElement!=currRoot)
            {
                totalBC[currElement]+=tree->delta[currElement];
            }

    #if QUEUE_TYPE==QUEUE_SINGLE
        qEnd--;
    }

    #elif QUEUE_TYPE==QUEUE_MULTI_ON_SINGLE
    }
    #elif QUEUE_TYPE==QUEUE_MULTI
        }
        deepestLevel--;
    }
    #endif


    return;
}


void PRE_NAME(bfsBrandes)(bcTree** treeArray, void* someGraph,float** totalBCArray,
                          void** someQueueParam,uint32_t*** parentArrayParam,
                          float* timePerThread, int64_t* selectedRoots,int32_t rootsPerThread)
{
    tic_reset();
    #pragma omp parallel
    {
	#if defined(_OPENMP)
        int thread = omp_get_thread_num();
	#else
	int thread = 0;
	#endif	

/*        int numExtras = ROOTS%omp_get_num_threads();
        int segSize = ROOTS/omp_get_num_threads();
        int haveExtra = numExtras<thread? 1: 0;
*/
        int start,stop;
/*
        if(numExtras>thread)
        {
            start = (segSize+1)*thread;
            stop = start+segSize+1;

        }
        else
        {
            start = (segSize+1)*numExtras + (segSize)*(thread-numExtras);
            stop = start+segSize;
        }
*/
        start = (thread)*rootsPerThread;
        stop     = (thread+1)*rootsPerThread;

//        printf("%d, %d, %d\n", thread,start, stop);

        int32_t* extraPointer=NULL;
        int32_t* extraPointer2=NULL;
        #if QUEUE_TYPE==QUEUE_MULTI_ON_SINGLE
            extraPointer=(int32_t*)malloc(sizeof(uint32_t)*(*treeArray)->NV);
        #endif
        #if ACCUM_TYPE==ACCUM_PARENT_ARRAY
            extraPointer2=(int32_t*)malloc(sizeof(uint32_t)*(*treeArray)->NV);
        #endif

        float* totalBC = totalBCArray[thread];
        bcTree* tree = treeArray[thread];
        uint32_t** parentArray=NULL;
        if(parentArrayParam!=NULL)
            parentArray = parentArrayParam[thread];
//        else
//			printf("problem with parent array\n");


        void* someQueue = someQueueParam[thread];

        for(int i=0; i<tree->NV; i++)
        {
            totalBC[i]=0.0;
        }

        for(uint32_t i = start; i < stop; i++) {
           // printf("%d,",i); fflush(stdout);
            PRE_NAME(bfsBrandesPerTree)(tree,someGraph,selectedRoots[i],totalBC,someQueue,parentArray,extraPointer,extraPointer2);
        }


        #if QUEUE_TYPE==QUEUE_MULTI_ON_SINGLE
            free(extraPointer);
        #endif
        #if ACCUM_TYPE==ACCUM_PARENT_ARRAY
            free(extraPointer2);
        #endif


        timePerThread[thread] = tic_sinceReset();
    }
}

