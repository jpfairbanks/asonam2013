void
static_multi_contract_clustering (
    uint64_t ** matches,
    uint64_t nv,
    struct stinger * S,
    struct stinger * S_orig);

double
modularity_score (
    struct stinger * S, 
    uint64_t * cluster, 
    uint64_t nv, 
    uint64_t ne);

uint64_t connected_components (struct stinger * S, int64_t * components, uint64_t nv);
