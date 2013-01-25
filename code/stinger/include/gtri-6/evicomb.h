#ifndef EVICOMB_H
#define EVICOMB_H

#include "stinger.h"
#include "evicomb-internal.h"
#include "alg_convert.h"

void evicomb(struct stinger *s, int anomaly_count, int *anomaly_ids, char *params, char * path);

#endif
