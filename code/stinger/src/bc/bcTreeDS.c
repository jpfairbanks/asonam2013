

#include <alloca.h>
#include <stdio.h>
#include <stdlib.h>
#include "xmalloc.h"

#include "bcTreeDS.h"


/**
* @brief Creates the data structures needed for computing BC.
*       These include the level of each vertex in the BFS tree.
*       The number of shortest paths each vertex has to the root.
*       The delta value computed in the dependency accumulation stage/
*
* @param numVertices The number of vertices in the graph
*
* @return Returns the created data structure.
*/
bcTree* CreateTree(int numVertices)
{
    bcTree* newTree;
    newTree = (bcTree*)xmalloc(sizeof(bcTree));

    newTree->NV = numVertices;

    newTree->level = (int*)xcalloc(numVertices,sizeof(int));
    newTree->pathsToRoot = (int*)xcalloc(numVertices,sizeof(int));
    newTree->delta = (bc_t*)xcalloc(numVertices,sizeof(bc_t));

#if PARENT_INFO_ON
    makeArrayOfLists(&newTree->parentList,numVertices);
#endif



    return newTree;
}

/**
* @brief Destroys the allocated tree.
*
* @param deadTree The tree that needs to be unallocated.
*
* @return None
*/
void DestroyTree(bcTree* deadTree)
{
    free(deadTree->level);
    free(deadTree->pathsToRoot);
    free(deadTree->delta);

#if PARENT_INFO_ON
    destroyArrayOfLists(&deadTree->parentList, deadTree->NV);
#endif


    free(deadTree);
}

/**
* @brief Creates a tree for each vertex in the graph. Thus, a total of NV trees are created.
*
* @param newForest This is an OUT parameter that contains the entire forest.
* @param numVertices The number of vertices in the graph
*
* @return None (see parameter newForest)
*/
void CreateForest(bcForest** newForest, int numVertices)
{
    *newForest = (bcForest*)xmalloc(sizeof(bcForest));

    (*newForest)->NV = numVertices;
    (*newForest)->forest = (bcTree**)xcalloc(numVertices,sizeof(bcTree*));
    (*newForest)->totalBC = (bc_t*)xcalloc(numVertices,sizeof(bc_t));

    int i;
    for(i=0;i<numVertices;i++)
    {
        (*newForest)->forest[i] = CreateTree(numVertices);
    }
}



/**
* @brief Destroys the allocated forest.
*
* @param deadForest The forest that needs to be unallocated.
*
* @return None
*/
void DestroyForest(bcForest** deadForest)
{
    int i;
    for(i=0;i<(*deadForest)->NV ;i++)
    {
        DestroyTree((*deadForest)->forest[i]);
    }
    free((*deadForest)->totalBC);

    free((*deadForest)->forest);

    free(*deadForest);

    return;
}


