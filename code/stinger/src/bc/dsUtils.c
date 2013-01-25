
#include "dsUtils.h"
#include "xmalloc.h"

// Creates the parent array needed for the computation of BC.
// The output of this function is an array of arrays. Each vertex is allocated an array based on it's adjaency.
uint32_t** createParentArray(csrGraph* graph,uint32_t NV)
{
    uint32_t** parentArray = (uint32_t**)malloc(NV*sizeof(uint32_t*));
    for(uint32_t v=0; v<NV;v++)
    {
        uint32_t edgeCount=graph->vertexPointerArray[v+1]-graph->vertexPointerArray[v];
        parentArray[v]=(uint32_t*)malloc(edgeCount*sizeof(uint32_t));
    }

    return parentArray;
}

// Destroys the parent array
void destroyParentArray(uint32_t** parentArray,uint32_t NV)
{

    for(uint32_t v=0; v<NV;v++)
    {
        free(parentArray[v]);
    }

    free(parentArray);
}

// Creates a parent array for each thread/core.
uint32_t*** createParallelParentArray(csrGraph* graph,uint32_t NV, uint32_t threadCount)
{
    uint32_t*** parallelParentArray = (uint32_t***)malloc(threadCount*sizeof(uint32_t**));

    for(uint32_t t=0; t<threadCount;t++)
    {
        parallelParentArray[t]=createParentArray(graph,NV);
        if(parallelParentArray[t]==NULL)
            printf("Failed to allocated memory for parallel parent array\n");
    }

    return parallelParentArray;
}

// Destroys the parent array of each thread/core.
void destroyParallelParentArray(uint32_t*** parallelParentArray,uint32_t NV,uint32_t threadCount)
{

    for(uint32_t t=0; t<threadCount;t++)
    {
        destroyParentArray(parallelParentArray[t],NV);
    }

    free(parallelParentArray);
}


 // Creates the parent array needed for the computation of BC.
// The output of this function is an array of arrays. Each vertex is allocated an array based on it's adjaency.
uint32_t** createParentArrayStinger(struct stinger* S,uint32_t NV)
{
    uint32_t** parentArray = (uint32_t**)malloc(NV*sizeof(uint32_t*));
    for(uint32_t v=0; v<NV;v++) {
	parentArray[v]=(uint32_t*)malloc(stinger_outdegree(S,v)*sizeof(uint32_t));
    }
    return parentArray;
}


// Creates a parent array for each thread/core.
uint32_t*** createParallelParentArrayStinger(struct stinger* S,uint32_t NV, uint32_t threadCount)
{
    uint32_t*** parallelParentArray = (uint32_t***)malloc(threadCount*sizeof(uint32_t**));

 //   	printf("thread count %d   %d \n", threadCount,NV);
   
	for(uint32_t t=0; t<threadCount;t++)
    {
//		printf("thread count %d\n", threadCount);
        parallelParentArray[t]=createParentArrayStinger(S,NV);
        if(parallelParentArray[t]==NULL)
            printf("Failed to allocated memory for parallel parent array\n");
    }

    return parallelParentArray;
}


 





// Creates the multi-level queue for the computation of BC. In this case O(NV^2) is allocated.
uint32_t** createMultiLevelQueue(uint32_t NV)
{
    uint32_t** multiLevelQueue = (uint32_t**)malloc(NV*sizeof(uint32_t*));

    for(uint32_t v=0; v<NV;v++)
    {
        multiLevelQueue[v]=(uint32_t*)malloc(NV*sizeof(uint32_t));
    }

    return multiLevelQueue;
}

// Destroys the multi level queue
void destroyMultiLevelQueue(uint32_t** multiLevelQueue,uint32_t NV)
{

    for(uint32_t v=0; v<NV;v++)
    {
        free(multiLevelQueue[v]);
    }

    free(multiLevelQueue);
}

// Creates a multi level queue for each thread/core.
uint32_t*** createParallelMultiLevelQueue(uint32_t NV, uint32_t threadCount)
{
    uint32_t*** parallelMultiLevelQueue = (uint32_t***)malloc(threadCount*sizeof(uint32_t**));

    for(uint32_t t=0; t<threadCount;t++)
    {
        parallelMultiLevelQueue[t]=createMultiLevelQueue(NV);
    }

    return parallelMultiLevelQueue;
}
// Destroys the multi level queue of each thread/core.
void destroyParallelMultiLevelQueue(uint32_t*** parallelMultiLevelQueue,uint32_t NV,uint32_t threadCount)
{

    for(uint32_t t=0; t<threadCount;t++)
    {
        destroyMultiLevelQueue(parallelMultiLevelQueue[t],NV);
    }

    free(parallelMultiLevelQueue);
}



// Creates the parent list needed for the computation of BC.
// The output of this function is an array of arrays. Each vertex is allocated an array based on it's adjaency.
bcTree** createParallelForest(int threadCount,int NV)
{
    bcTree** parallelForest = (bcTree**)xcalloc(threadCount,sizeof(bcTree*));

    for(int i=0;i<threadCount;i++)
    {
        parallelForest[i] = CreateTree(NV);
    }
    return parallelForest;

}
// Destroys the parent array
void destroyParallelForest(bcTree** parallelForest, int threadCount)
{
    for(int i=0;i<threadCount;i++)
    {
        DestroyTree(parallelForest[i]);
    }

    free(parallelForest);

}

// Creates a parent list for each thread/core.
list_ptr** createParallelList(int threadCount,int NV)
{
    list_ptr** parallelList = (list_ptr**)xcalloc(threadCount,sizeof(list_ptr**));
    for(int i=0; i<threadCount;i++)
    {
        makeArrayOfLists(&parallelList[i],NV);
    }

    return parallelList;
}

// Destroys the parent list of each thread/core.
void destroyParallelList(list_ptr** parallelList, int threadCount,int NV)
{
    for(int i=0; i<threadCount;i++)
    {
        destroyArrayOfLists(&parallelList[i],NV);
    }

    free(parallelList);
}


float** createParallelBetweennessArray(int threadCount,int NV)
{
    float** totalBC = (float**)malloc((threadCount)*sizeof(float*));

    for(int i=0;i<threadCount;i++)
    {
        totalBC[i] = malloc(sizeof(float)*NV);
    }
    return totalBC;
}

void destroyParallelBetweennessArray(float** parallelScore, int threadCount)
{
    for(int i=0; i<threadCount;i++)
    {
//        printf("destroying %d\n",i); fflush(stdout);
        free(parallelScore[i]);
    }
//        printf("destroying last\n"); fflush(stdout);
    free(parallelScore);
}
