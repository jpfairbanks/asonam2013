
#undef FOR_ADJ
#undef FOR_STINGER
#undef FOR_CSR
#define FOR_ADJ 1
#define FOR_STINGER 2
#define FOR_CSR 3

#define QUEUE_SINGLE 1
#define QUEUE_MULTI_ON_SINGLE 2
#define QUEUE_MULTI 3

#define ACCUM_TRAVERSAL 1
#define ACCUM_PARENT_LIST 2
#define ACCUM_PARENT_ARRAY 3

#define GO(x,y,z,w) x##y##z##w



#define ACCUM_TYPE ACCUM_PARENT_ARRAY

#undef FOR_TYPE
#undef PRE_NAME
#undef QUEUE_TYPE
#define FOR_TYPE FOR_STINGER
#define QUEUE_TYPE QUEUE_SINGLE
#define PRE_NAME(x) x##STINGER##SINGLE##PA
#include "genericBC.c"

#undef FOR_TYPE
#undef PRE_NAME
#undef QUEUE_TYPE
#define FOR_TYPE FOR_ADJ
#define QUEUE_TYPE QUEUE_SINGLE
#define PRE_NAME(x) x##ADJ##SINGLE##PA
#include "genericBC.c"


#undef FOR_TYPE
#undef PRE_NAME
#undef QUEUE_TYPE
#define FOR_TYPE FOR_CSR
#define QUEUE_TYPE QUEUE_SINGLE
#define PRE_NAME(x) x##CSR##SINGLE##PA
#include "genericBC.c"


#undef FOR_TYPE
#undef PRE_NAME
#undef QUEUE_TYPE
#define FOR_TYPE FOR_STINGER
#define QUEUE_TYPE QUEUE_MULTI
#define PRE_NAME(x) x##STINGER##MULTI##PA
#include "genericBC.c"

#undef FOR_TYPE
#undef PRE_NAME
#undef QUEUE_TYPE
#define FOR_TYPE FOR_ADJ
#define QUEUE_TYPE QUEUE_MULTI
#define PRE_NAME(x) x##ADJ##MULTI##PA
#include "genericBC.c"


#undef FOR_TYPE
#undef PRE_NAME
#undef QUEUE_TYPE
#define FOR_TYPE FOR_CSR
#define QUEUE_TYPE QUEUE_MULTI
#define PRE_NAME(x) x##CSR##MULTI##PA
#include "genericBC.c"



#undef FOR_TYPE
#undef PRE_NAME
#undef QUEUE_TYPE
#define FOR_TYPE FOR_STINGER
#define QUEUE_TYPE QUEUE_MULTI_ON_SINGLE
#define PRE_NAME(x) x##STINGER##M_ON_S##PA
#include "genericBC.c"

#undef FOR_TYPE
#undef PRE_NAME
#undef QUEUE_TYPE
#define FOR_TYPE FOR_ADJ
#define QUEUE_TYPE QUEUE_MULTI_ON_SINGLE
#define PRE_NAME(x) x##ADJ##M_ON_S##PA
#include "genericBC.c"


#undef FOR_TYPE
#undef PRE_NAME
#undef QUEUE_TYPE
#define FOR_TYPE FOR_CSR
#define QUEUE_TYPE QUEUE_MULTI_ON_SINGLE
#define PRE_NAME(x) x##CSR##M_ON_S##PA
#include "genericBC.c"









#undef ACCUM_TYPE
#define ACCUM_TYPE ACCUM_PARENT_LIST

#undef FOR_TYPE
#undef PRE_NAME
#undef QUEUE_TYPE
#define FOR_TYPE FOR_STINGER
#define QUEUE_TYPE QUEUE_SINGLE
#define PRE_NAME(x) x##STINGER##SINGLE##PL
#include "genericBC.c"

#undef FOR_TYPE
#undef PRE_NAME
#undef QUEUE_TYPE
#define FOR_TYPE FOR_ADJ
#define QUEUE_TYPE QUEUE_SINGLE
#define PRE_NAME(x) x##ADJ##SINGLE##PL
#include "genericBC.c"


#undef FOR_TYPE
#undef PRE_NAME
#undef QUEUE_TYPE
#define FOR_TYPE FOR_CSR
#define QUEUE_TYPE QUEUE_SINGLE
#define PRE_NAME(x) x##CSR##SINGLE##PL
#include "genericBC.c"


#undef FOR_TYPE
#undef PRE_NAME
#undef QUEUE_TYPE
#define FOR_TYPE FOR_STINGER
#define QUEUE_TYPE QUEUE_MULTI
#define PRE_NAME(x) x##STINGER##MULTI##PL
#include "genericBC.c"

#undef FOR_TYPE
#undef PRE_NAME
#undef QUEUE_TYPE
#define FOR_TYPE FOR_ADJ
#define QUEUE_TYPE QUEUE_MULTI
#define PRE_NAME(x) x##ADJ##MULTI##PL
#include "genericBC.c"


#undef FOR_TYPE
#undef PRE_NAME
#undef QUEUE_TYPE
#define FOR_TYPE FOR_CSR
#define QUEUE_TYPE QUEUE_MULTI
#define PRE_NAME(x) x##CSR##MULTI##PL
#include "genericBC.c"




#undef FOR_TYPE
#undef PRE_NAME
#undef QUEUE_TYPE
#define FOR_TYPE FOR_STINGER
#define QUEUE_TYPE QUEUE_MULTI_ON_SINGLE
#define PRE_NAME(x) x##STINGER##M_ON_S##PL
#include "genericBC.c"

#undef FOR_TYPE
#undef PRE_NAME
#undef QUEUE_TYPE
#define FOR_TYPE FOR_ADJ
#define QUEUE_TYPE QUEUE_MULTI_ON_SINGLE
#define PRE_NAME(x) x##ADJ##M_ON_S##PL
#include "genericBC.c"


#undef FOR_TYPE
#undef PRE_NAME
#undef QUEUE_TYPE
#define FOR_TYPE FOR_CSR
#define QUEUE_TYPE QUEUE_MULTI_ON_SINGLE
#define PRE_NAME(x) x##CSR##M_ON_S##PL
#include "genericBC.c"






#undef ACCUM_TYPE
#define ACCUM_TYPE ACCUM_TRAVERSAL

#undef FOR_TYPE
#undef PRE_NAME
#undef QUEUE_TYPE
#define FOR_TYPE FOR_STINGER
#define QUEUE_TYPE QUEUE_SINGLE
#define PRE_NAME(x) x##STINGER##SINGLE##PT
#include "genericBC.c"

#undef FOR_TYPE
#undef PRE_NAME
#undef QUEUE_TYPE
#define FOR_TYPE FOR_ADJ
#define QUEUE_TYPE QUEUE_SINGLE
#define PRE_NAME(x) x##ADJ##SINGLE##PT
#include "genericBC.c"


#undef FOR_TYPE
#undef PRE_NAME
#undef QUEUE_TYPE
#define FOR_TYPE FOR_CSR
#define QUEUE_TYPE QUEUE_SINGLE
#define PRE_NAME(x) x##CSR##SINGLE##PT
#include "genericBC.c"


#undef FOR_TYPE
#undef PRE_NAME
#undef QUEUE_TYPE
#define FOR_TYPE FOR_STINGER
#define QUEUE_TYPE QUEUE_MULTI
#define PRE_NAME(x) x##STINGER##MULTI##PT
#include "genericBC.c"

#undef FOR_TYPE
#undef PRE_NAME
#undef QUEUE_TYPE
#define FOR_TYPE FOR_ADJ
#define QUEUE_TYPE QUEUE_MULTI
#define PRE_NAME(x) x##ADJ##MULTI##PT
#include "genericBC.c"


#undef FOR_TYPE
#undef PRE_NAME
#undef QUEUE_TYPE
#define FOR_TYPE FOR_CSR
#define QUEUE_TYPE QUEUE_MULTI
#define PRE_NAME(x) x##CSR##MULTI##PT
#include "genericBC.c"



#undef FOR_TYPE
#undef PRE_NAME
#undef QUEUE_TYPE
#define FOR_TYPE FOR_STINGER
#define QUEUE_TYPE QUEUE_MULTI_ON_SINGLE
#define PRE_NAME(x) x##STINGER##M_ON_S##PT
#include "genericBC.c"

#undef FOR_TYPE
#undef PRE_NAME
#undef QUEUE_TYPE
#define FOR_TYPE FOR_ADJ
#define QUEUE_TYPE QUEUE_MULTI_ON_SINGLE
#define PRE_NAME(x) x##ADJ##M_ON_S##PT
#include "genericBC.c"


#undef FOR_TYPE
#undef PRE_NAME
#undef QUEUE_TYPE
#define FOR_TYPE FOR_CSR
#define QUEUE_TYPE QUEUE_MULTI_ON_SINGLE
#define PRE_NAME(x) x##CSR##M_ON_S##PT
#include "genericBC.c"





