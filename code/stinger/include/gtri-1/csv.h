#ifndef  CSV_H
#define  CSV_H

#include <stdio.h>
#include "alg_convert.h"

void
readCSVLineDynamic(char delim, FILE * file, char ** buf, uint64_t * bufSize, char *** fields, uint64_t ** lengths, uint64_t * fieldsSize, uint64_t * count);

void
splitLineCSVDynamic(char delim, char * line, uint64_t lineSize, char ** buf, uint64_t * bufSize, char *** fields, uint64_t ** lengths, uint64_t * fieldsSize, uint64_t * count);

int
getIndex(char * string, char ** fields, uint64_t * lengths, uint64_t count);

void
printLine(char ** fields, uint64_t * lengths, uint64_t count);

void
csvIfIDExistsint64(FILE * fp, char delim, struct stinger * s, const char ** type_strings, uint64_t nv, int64_t * values);

void
csvIfIDExistsfloat(FILE * fp, char delim, struct stinger * s, const char ** type_strings, uint64_t nv, float * values);

void
csvIfIDExistsdouble(FILE * fp, char delim, struct stinger * s, const char ** type_strings, uint64_t nv, double * values);

void
csvIfIDExistsbeliefs(FILE * fp, char delim, struct stinger * S, const char ** type_strings, uint64_t nv, double * beliefs);

#endif  /*CSV_H*/
