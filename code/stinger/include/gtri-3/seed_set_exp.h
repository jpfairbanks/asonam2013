#ifndef  SEED_SET_EXP_H
#define  SEED_SET_EXP_H

#define HISTORY 200

void seed_set_expansion(const struct stinger *S, int64_t *seeds, int64_t num_seeds, int64_t nv, int64_t ne, int64_t *membership, int64_t *neighbors, int64_t *edge_between, double * neighbor_modval, int64_t *total_community_edges, int64_t *internal_community_edges, int64_t *number_neighbors,double *history, int64_t *history_idx);
int64_t contract(int64_t *neigh_array,int64_t* edge_between, int64_t length);
int64_t merge(int64_t *neighbors, int64_t *edge_between, int64_t *temp_neighbors, int64_t *temp_edge_between, int64_t length1, int64_t length2);
double max(double a, double b);
void update_community(const struct stinger *S, int64_t *seeds, int64_t num_seeds, int64_t nv, int64_t ne, int64_t *membership, int64_t *neighbors,  int64_t *total_community_edges, int64_t *internal_community_edges, int64_t *affected, int64_t num_affected, double *history, int64_t *history_idx);

#endif  /*SEED_SET_EXP_H*/
