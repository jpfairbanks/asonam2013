#pragma once

#include <stdint.h>

#include "stinger.h"
#include "csr.h"

#define sigma(j,k) pathsToRoot[j][k]
#define sigmaf(j,k) ((bc_t)pathsToRoot[j][k])
#define d(j,k) level[j][k]

#define INFINITY_MY 1073741824

#define VB_ON 0
#define DB_ON 0
#define INFO_ON 0
#define STREAMING_INFO_ON 0
#define COUNTING_INFO_ON 0

#define PARENT_INFO_ON 0

#define BEFORE_AND_AFTER 1

#if VB_ON
#define VB(...) printf(__VA_ARGS__); fflush(stdout);
#else
#define VB(...)
#endif

#if DB_ON
#define DB(...) printf(__VA_ARGS__); fflush(stdout);
#else
#define DB(...)
#endif

#if INFO_ON
#define INFO(...) printf(__VA_ARGS__); fflush(stdout);
#else
#define INFO(...)
#endif

#if STREAMING_INFO_ON
#define STREAMING_INFO(...) printf(__VA_ARGS__); fflush(stdout);
#else
#define STREAMING_INFO(...)
#endif

#ifdef COUNTING_INFO_ON
extern uint64_t countFullVertex;
extern uint64_t countFullEdgeUp;
extern uint64_t countStreamVertexDown;
extern uint64_t countStreamVertexUp;
extern uint64_t countStreamEdgeDown;
extern uint64_t countStreamEdgeUp;
#endif


typedef double bc_t;

uint32_t ** squareMat(uint32_t NV);
void freeSquareMat(uint32_t ** sq, uint32_t NV);
bc_t ** squareMatFloat(uint32_t NV);
void freeSquareMatFloat(bc_t ** sq, uint32_t NV);
void printGraph(uint32_t ** graph, uint32_t NV);

void intersectAdjacencyWithLevel(uint32_t * outArray, uint32_t * numFound, uint32_t NV,
	uint32_t * adjacency, uint32_t * levelArray, uint32_t level);
void moveUpTree(uint32_t * levels, uint32_t * pathsToRoot, uint32_t ** graph,
	uint32_t vertex, uint32_t ignoreVertex, uint32_t dist, uint32_t NV);
void addEdgeWithoutMovement(uint32_t * levels, uint32_t * pathsToRoot, uint32_t ** graph,
	uint32_t vertex, uint32_t deltaTop, uint32_t NV);
void copyTreeIn(uint32_t * levelsDest, uint32_t * pathsToRootDest, uint32_t existingVertex, uint32_t newVertex,
	uint32_t *levelsSrc, uint32_t * pathsToRootSrc, uint32_t INF, uint32_t NV);


void MatrixBCComput(uint32_t NV, uint32_t * pathsToRoot[],uint32_t * level[],bc_t* perTreeBC[],bc_t* BC);


int64_t serial_shiloach_vishkin_components (struct stinger *S, int64_t nv, int64_t * component_map);
int64_t serial_shiloach_vishkin_componentsCSR (csrGraph* graph, int64_t * component_map);

void selectRootsInMaxComponents(csrGraph* graph,int NV,int64_t* component_map, int64_t componentCount,
                                int numberOfRoots,int* countedEdges, int64_t* selectedRoots);


void read_GML_graph(struct stinger** G, char* filename);
void readSnapGraph(struct stinger** G, char* filename);


void hostParseArgs(int argc, char** argv);
int updateEdge(int argc, char *argv[]);
