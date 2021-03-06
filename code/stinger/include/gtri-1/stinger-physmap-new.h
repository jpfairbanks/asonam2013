#include "stinger.h"

int
stinger_new_physID_to_vtx_mapping(struct stinger * S, char * byte_string, uint64_t length, int64_t * vtx_out);

int64_t
stinger_physID_to_vtx(struct stinger * S, char * byte_string, uint64_t length);

int
stinger_vtx_to_physID(struct stinger * S, uint64_t vertexID, char ** outbuffer, uint64_t * outbufferlength);

int
stinger_vtx_to_physID_direct(struct stinger * S, uint64_t vertexID, char ** out_ptr, uint64_t * out_len);

void
stinger_save_physmap_to(struct stinger * S, char * name);

void
stinger_load_physmap_from(struct stinger * S, char * name);

void
stinger_save_physmap_to_file(struct stinger * S, FILE *);

void
stinger_load_physmap_from_file(struct stinger * S, FILE * fp);
