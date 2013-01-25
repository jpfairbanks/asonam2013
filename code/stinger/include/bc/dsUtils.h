#pragma once


#include "stinger.h"
#include "bcTreeDS.h"
#include "list.h"
#include "csr.h"

uint32_t** createParentArray(csrGraph* graph,uint32_t NV);
void destroyParentArray(uint32_t** parentArray,uint32_t NV);

uint32_t*** createParallelParentArray(csrGraph* graph,uint32_t NV, uint32_t threadCount);
void destroyParallelParentArray(uint32_t*** parallelParentArray,uint32_t NV,uint32_t threadCount);

uint32_t** createParentArrayStinger(struct stinger* s,uint32_t NV);
uint32_t*** createParallelParentArrayStinger(struct stinger* S,uint32_t NV, uint32_t threadCount); 


uint32_t** createMultiLevelQueue(uint32_t NV);
void destroyMultiLevelQueue(uint32_t** multiLevelQueue,uint32_t NV);
uint32_t*** createParallelMultiLevelQueue(uint32_t NV, uint32_t threadCount);
void destroyParallelMultiLevelQueue(uint32_t*** parallelMultiLevelQueue,uint32_t NV,uint32_t threadCount);

bcTree** createParallelForest(int threadCount,int NV);
void destroyParallelForest(bcTree** parallelForest, int threadCount);


list_ptr** createParallelList(int threadCount,int NV);
void destroyParallelList(list_ptr** parallelList, int threadCount,int NV);

float** createParallelBetweennessArray(int threadCount,int NV);
void destroyParallelBetweennessArray(float** parallelList, int threadCount);
