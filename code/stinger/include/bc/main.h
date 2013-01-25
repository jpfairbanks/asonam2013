
#pragma once


typedef enum{
    TC_ER = 0,
    TC_RMAT,
    TC_REAL_DATA,       // DOES NOT USE ADJACENCY MATRIX
    TC_CC_ER,
    TC_SAME_GRAPH_ER,
    TC_SAME_GRAPH_RMAT,
} testCase;

typedef enum{
    UP_INSERT = 0,
    UP_DELETE,
} updateType;


extern uint32_t noComputeDiffComp;

extern uint32_t NV;
extern uint32_t NE;
extern uint32_t randomSeed;
extern uint32_t iterationCount;
extern char graphDataName[1024];
extern updateType opType;



extern testCase graphTestCase; // GraphType ER=0,RMAT=1,real graph=2, connec

typedef struct
{
    uint32_t movement;
    uint32_t adjacent;
    uint32_t sameLevel;
    uint32_t connectedComponents;
}StreamingExtraInfo;




