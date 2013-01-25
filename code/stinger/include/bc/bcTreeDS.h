#pragma once

#include "streaming_utils.h"

#if PARENT_INFO_ON
#include "list.h"
#endif

typedef struct {
    int NV;
    int* level;
    int* pathsToRoot;
    bc_t* delta;

#if PARENT_INFO_ON
    list_ptr* parentList;
#endif

} bcTree;

typedef struct {
    bcTree** forest;
    bc_t* totalBC;
    int NV;
} bcForest;


bcTree* CreateTree(int numVertices);
void DestroyTree(bcTree* deadTree);

void CreateForest(bcForest** newForest, int numVertices);
void DestroyForest(bcForest** deadForest);


