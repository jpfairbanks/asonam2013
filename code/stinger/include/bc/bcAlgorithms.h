#pragma once

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

#include "csr.h"
#include "bcTreeDS.h"

void betCentCSR_TraverseNeig_BFSqueue(csrGraph* graph,uint32_t currRoot,float* totalBC,bcTree* tree);
void betCentStinger_TraverseNeig_BFSqueue(struct stinger* sStinger,uint32_t currRoot,float* totalBC,bcTree* tree);
void betCentStinger_TraverseNeig_BFSqueue2(struct stinger* sStinger,uint32_t currRoot,float* totalBC,bcTree* tree);


