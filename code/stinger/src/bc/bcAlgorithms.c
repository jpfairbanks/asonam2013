
#include "bcAlgorithms.h"


void betCentCSR_TraverseNeig_BFSqueue(csrGraph* graph,uint32_t currRoot,float* totalBC,bcTree* tree)
{

    for(uint32_t j = 0; j < tree->NV; j++)
    {
        tree->level[j] = INFINITY_MY;
        tree->pathsToRoot[j] = INFINITY_MY;
        tree->delta[j] = 0;
    }
	tree->level[currRoot] = 0;
	tree->pathsToRoot[currRoot] = 1;

    uint32_t Stack[tree->NV];
    uint32_t Queue[tree->NV];

    Queue[0] = currRoot;
    int32_t qStart=0,qEnd=1;
    int32_t sStart=0;

    // While queue is not empty
    while(qStart!=qEnd)
    {
        uint32_t currElement = Queue[qStart];
        Stack[sStart] = currElement;
        sStart++;
        qStart++;

        int startEdge = graph->vertexPointerArray[currElement];
        int stopEdge = graph->vertexPointerArray[currElement+1];

        for(int j=startEdge;startEdge<stopEdge;startEdge++)
        {
            uint32_t k = graph->edgeArray[startEdge];

            // If this is a neighbor and has not been found
            if(tree->level[k] > tree->level[currElement])
            {
                // Checking if "k" has been found.
                if(tree->level[k]==INFINITY_MY)
                {
                    tree->level[k] = tree->level[currElement]+1;
                    Queue[qEnd++] = k;
                    tree->delta[k]=0;
                }

                if(tree->pathsToRoot[k] == INFINITY_MY)
                {
                    // k has not been found and therefore its paths to the roots are through its parent.
                    tree->pathsToRoot[k] = tree->pathsToRoot[currElement];
                }
                else
                {
                    // k has been found and has multiple paths to the root as it has multiple parents.
                    tree->pathsToRoot[k] += tree->pathsToRoot[currElement];
                }
            }
        }
    }

    // Using Brandes algorithm to compute BC for a specific tree.
    // Essentially, use the stack which the elements are placed in depth-reverse order, to "climb" back
    // up the tree, all the way to the root.

    int32_t sEnd = sStart-1;

    while(sEnd>=0)
    {
        uint32_t currElement = Stack[sEnd];

        int startEdge = graph->vertexPointerArray[currElement];
        int stopEdge = graph->vertexPointerArray[currElement+1];

        for(int j=startEdge;startEdge<stopEdge;startEdge++)
        {
            uint32_t k = graph->edgeArray[startEdge];

            // If this is a neighbor and has not been found
            if((tree->level[k] == (tree->level[currElement]-1)))
            {

                   tree->delta[k] +=
                        ((bc_t)tree->pathsToRoot[k]/(bc_t)tree->pathsToRoot[currElement])*
                        (bc_t)(tree->delta[currElement]+1);
            }
        }

        if(currElement!=currRoot)
        {
            totalBC[currElement]+=tree->delta[currElement];
        }

        sEnd--;
    }

}

void betCentStinger_TraverseNeig_BFSqueue(struct stinger* sStinger,uint32_t currRoot,float* totalBC,bcTree* tree)
{

    for(uint32_t j = 0; j < tree->NV; j++)
    {
        tree->level[j] = INFINITY_MY;
        tree->pathsToRoot[j] = INFINITY_MY;
        tree->delta[j] = 0;

    }
	tree->level[currRoot] = 0;
	tree->pathsToRoot[currRoot] = 1;

    uint32_t Stack[tree->NV];
    uint32_t Queue[tree->NV];

    Queue[0] = currRoot;
    int32_t qStart=0,qEnd=1;
    int32_t sStart=0;

    // While queue is not empty
    while(qStart!=qEnd)
    {
        uint32_t currElement = Queue[qStart];
        Stack[sStart] = currElement;
        sStart++;
        qStart++;

        STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger,currElement)
        {
            uint32_t k = STINGER_EDGE_DEST;

            // If this is a neighbor and has not been found
            if(tree->level[k] > tree->level[currElement])
            {
                // Checking if "k" has been found.
                if(tree->level[k]==INFINITY_MY)
                {
                    tree->level[k] = tree->level[currElement]+1;
                    Queue[qEnd++] = k;
                    tree->delta[k]=0;
                }

                if(tree->pathsToRoot[k] == INFINITY_MY)
                {
                    // k has not been found and therefore its paths to the roots are through its parent.
                    tree->pathsToRoot[k] = tree->pathsToRoot[currElement];

                }
                else
                {
                    // k has been found and has multiple paths to the root as it has multiple parents.
                    tree->pathsToRoot[k] += tree->pathsToRoot[currElement];
                }
            }
        }
        STINGER_FORALL_EDGES_OF_VTX_END();
    }

    // Using Brandes algorithm to compute BC for a specific tree.
    // Essentially, use the stack which the elements are placed in depth-reverse order, to "climb" back
    // up the tree, all the way to the root.

    int32_t sEnd = sStart-1;

    while(sEnd>=0)
    {
        uint32_t currElement = Stack[sEnd];

        STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger,currElement)
        {
            uint32_t k = STINGER_EDGE_DEST;

            // If this is a neighbor and has not been found
            if((tree->level[k] == (tree->level[currElement]-1)))
            {
                   tree->delta[k] +=
                        ((bc_t)tree->pathsToRoot[k]/(bc_t)tree->pathsToRoot[currElement])*
                        (bc_t)(tree->delta[currElement]+1);
            }

        }

        STINGER_FORALL_EDGES_OF_VTX_END();

        if(currElement!=currRoot)
        {
            totalBC[currElement]+=tree->delta[currElement];
        }

        sEnd--;
    }
}

/*
void betCentStinger_TraverseNeig_BFSqueue2(struct stinger* sStinger,uint32_t currRoot,float* totalBC,bcTree* tree)
{
    int NV = tree->NV;
    int level[NV], pathsToRoot[NV];
    float delta[NV];
    for(uint32_t j = 0; j < tree->NV; j++)
    {
        level[j] = INFINITY_MY;
        pathsToRoot[j] = INFINITY_MY;
        delta[j] = 0;

    }
	level[currRoot] = 0;
	pathsToRoot[currRoot] = 1;

    uint32_t Stack[NV];
    uint32_t Queue[NV];

    Queue[0] = currRoot;
    int32_t qStart=0,qEnd=1;
    int32_t sStart=0;

    // While queue is not empty
    while(qStart!=qEnd)
    {
        uint32_t currElement = Queue[qStart];
        Stack[sStart] = currElement;
        sStart++;
        qStart++;

        STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger,currElement)
        {
            uint32_t k = STINGER_EDGE_DEST;

            // If this is a neighbor and has not been found
            if(level[k] > level[currElement])
            {
                // Checking if "k" has been found.
                if(level[k]==INFINITY_MY)
                {
                    level[k] = level[currElement]+1;
                    Queue[qEnd++] = k;
                    delta[k]=0;
                }

                if(pathsToRoot[k] == INFINITY_MY)
                {
                    // k has not been found and therefore its paths to the roots are through its parent.
                    pathsToRoot[k] = pathsToRoot[currElement];

                }
                else
                {
                    // k has been found and has multiple paths to the root as it has multiple parents.
                    pathsToRoot[k] += pathsToRoot[currElement];
                }
            }
        }
        STINGER_FORALL_EDGES_OF_VTX_END();
    }

    // Using Brandes algorithm to compute BC for a specific tree.
    // Essentially, use the stack which the elements are placed in depth-reverse order, to "climb" back
    // up the tree, all the way to the root.

    int32_t sEnd = sStart-1;

    while(sEnd>=0)
    {
        uint32_t currElement = Stack[sEnd];

        STINGER_FORALL_EDGES_OF_VTX_BEGIN(sStinger,currElement)
        {
            uint32_t k = STINGER_EDGE_DEST;

            // If this is a neighbor and has not been found
            if((level[k] == (level[currElement]-1)))
            {
                   delta[k] +=
                        ((bc_t)pathsToRoot[k]/(bc_t)pathsToRoot[currElement])*
                        (bc_t)(delta[currElement]+1);
            }

        }

        STINGER_FORALL_EDGES_OF_VTX_END();

        if(currElement!=currRoot)
        {
            totalBC[currElement]+=delta[currElement];
        }

        sEnd--;
    }
}

*/
