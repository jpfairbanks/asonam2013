#include "stinger.h"
#include "x86-full-empty.h"
#include "xmalloc.h"

#include <string.h>

uint64_t xor_hash(uint8_t * byte_string, uint64_t length) {
  uint64_t out = 0;
  uint64_t * byte64 = (uint64_t *)byte_string;
  while(length > 8) {
    length -= 8;
    out ^= *(byte64);
    byte64++;
  }
  if(length > 0) {
    uint64_t last = 0;
    uint8_t * cur = (uint8_t *)&last;
    uint8_t * byte8 = (uint8_t *)byte64;
    while(length > 0) {
      length--;
      *(cur++) = *(byte8++);
    }
    out ^= last;
  }

  /* bob jenkin's 64-bit mix */
  out = (~out) + (out << 21); // out = (out << 21) - out - 1;
  out = out ^ (out >> 24);
  out = (out + (out << 3)) + (out << 8); // out * 265
  out = out ^ (out >> 14);
  out = (out + (out << 2)) + (out << 4); // out * 21
  out = out ^ (out >> 28);
  out = out + (out << 31);

  return out;
}

/**
* @brief Gets or creates the mapping between the input byte-string and the vertex space 0 .. STINGER_MAX_LVERTICES-1
*
* Bytes string can contain null and other special characters.  It is treated as raw data. This function
* does make a copy of the input string.
*
* @param S Pointer to the stinger graph.
* @param byte_string An arbirary length string of bytes.
* @param length The length of the string. Must be correct.
* @param vtx_out The output vertex mapping.
* @return 0 on existing mapping found. 1 on new mapping created. -1 on failure (STINGER is full).
*/
int
stinger_new_physID_to_vtx_mapping(struct stinger * S, char * byte_string, uint64_t length, int64_t * vtx_out)
{
  uint64_t vtx = xor_hash(byte_string, length);
  vtx = vtx % STINGER_MAX_LVERTICES; 
  uint64_t init_index = vtx;
  while(1) {
    /* hoping compiler does strength reduction
     * MAX is usually power of 2 
     */
    if(0 == readff((uint64_t *)&S->LVA[vtx].physID)) {
      int64_t original = readfe((uint64_t *)&S->LVA[vtx].physID);
      if(original) {
	if(length == S->LVA[vtx].physIDlength &&
	  bcmp(byte_string,&(S->physMapBuffer[original]), length) == 0) {
	    writeef((uint64_t *)&S->LVA[vtx].physID, original);
	    *vtx_out = vtx;
	    return 0;
	}
      } else {
	int64_t place = stinger_int64_fetch_add(&(S->physMapTop), length);
	char * physID = &(S->physMapBuffer[place]);
	memcpy(physID, byte_string, length);
	S->LVA[vtx].physIDlength = length;
	writeef((uint64_t *)&S->LVA[vtx].physID, (uint64_t)place);
	*vtx_out = vtx;
	return 1;
      }
      writeef((uint64_t *)&S->LVA[vtx].physID, original);
    } else if(length == S->LVA[vtx].physIDlength &&
      bcmp(byte_string, &(S->physMapBuffer[readff((uint64_t *)&S->LVA[vtx].physID)]), length) == 0) {
	*vtx_out = vtx;
	return 0;
    }
    vtx++;
    vtx = vtx % STINGER_MAX_LVERTICES; 
    if(vtx == init_index) 
      return -1;
  }
}

/**
* @brief Get the mapping between the input byte-string and the vertex space 0 .. STINGER_MAX_LVERTICES-1 if it exists.
*
* @param S Pointer to the stinger graph.
* @param byte_string An arbirary length string of bytes.
* @param length The length of the string. Must be correct.
* @retrun The vertex ID or -1 on failure.
*/
int64_t
stinger_physID_to_vtx(struct stinger * S, char * byte_string, uint64_t length)
{
  uint64_t vtx = xor_hash(byte_string, length);
  vtx = vtx % STINGER_MAX_LVERTICES; 
  uint64_t init_index = vtx;
  while(1) {
    /* hoping compiler does strength reduction
     * MAX is usually power of 2 
     */
    if(0 == readff((uint64_t *)&S->LVA[vtx].physID)) {
      return -1;
    } else if(length == S->LVA[vtx].physIDlength &&
      bcmp(byte_string, &(S->physMapBuffer[readff((uint64_t *)&S->LVA[vtx].physID)]), length) == 0) {
	return vtx;
    }
    vtx++;
    vtx = vtx % STINGER_MAX_LVERTICES; 
    if(vtx == init_index) 
      return -1;
  }
}

/**
* @brief Get the mapping between the input vertex in the space space 0 .. STINGER_MAX_LVERTICES-1 and its representative
*	string if it exists.
*
* This function assumes that the input buffer might could already be allocated, but will allocate / reallocate as needed.
*
* @param S Pointer to the stinger graph.
* @param vertexID The vertex from which to get the string.
* @param outbuffer Output buffer where you would like the data to be placed.
* @param outbufferlength The length of the output buffer (will be set or increased if allocated / reallocated or if string is shorter).
* @retrun 0 on success, -1 on failure.
*/
int
stinger_vtx_to_physID(struct stinger * S, uint64_t vertexID, char ** outbuffer, uint64_t * outbufferlength)
{
    if(0 == readff((uint64_t *)&S->LVA[vertexID].physID)) {
      return -1;
    } else {
      if(*outbuffer == NULL || *outbufferlength == 0){
	*outbuffer = xmalloc(S->LVA[vertexID].physIDlength * sizeof(char));
	if(NULL == *outbuffer) {
	  return -1;
	}
      } else if(*outbufferlength < S->LVA[vertexID].physIDlength) {
	void * tmp = xrealloc(*outbuffer, S->LVA[vertexID].physIDlength);
	if(tmp) {
	  *outbuffer = tmp;
	} else {
	  return -1;
	}
      }
      memcpy(*outbuffer, &(S->physMapBuffer[readff((uint64_t *)&S->LVA[vertexID].physID)]), S->LVA[vertexID].physIDlength);
      *outbufferlength = S->LVA[vertexID].physIDlength;
      return 0;
    }
}

int
stinger_vtx_to_physID_direct(struct stinger * S, uint64_t vertexID, char ** out_ptr, uint64_t * out_len)
{
    *out_ptr = &S->physMapBuffer[readff((uint64_t *)&S->LVA[vertexID].physID)];
    *out_len = S->LVA[vertexID].physIDlength;
    if(out_ptr)
      return 0;

    *out_len = 0;
    return -1;
}
