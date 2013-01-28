#include <stdint.h>
#include <stdio.h>
#include <dirent.h>

#include "stinger.h"
#include "stinger-atomics.h"
#include "stinger-utils.h"
#include "xmalloc.h"
#include "timer.h"

#define HERE() { printf("%s %d reached\n", __func__, __LINE__); fflush(stdout); }

double
sum_all_edgeweights(
    struct stinger * S, 
    int64_t type)
{
  double sum = 0;
  struct stinger_eb * ebpool_priv = S->ebpool->ebpool;

  /* for each edge block */
  OMP("omp parallel for reduction(+:sum)")						
  MTA("mta assert parallel")						
  for(uint64_t eb_index = 0; eb_index < S->ETA[(type)].high; eb_index++) {	
    struct stinger_eb *  cur_eb = ebpool_priv + S->ETA[(type)].blocks[eb_index]; 
    uint64_t stop = stinger_eb_high(cur_eb);
    for(uint64_t e = 0; e < stop; e++) { 
      if(!stinger_eb_is_blank(cur_eb, e)) {                   
	sum += cur_eb->edges[e].weight;
      }								
    } 								
  }									
  return sum;
}

void
sequence(uint64_t * arr, uint64_t count) {
  OMP("omp parallel for")
  for(uint64_t i = 0; i < count; i++) {
    arr[i] = i;
  }
}

void
zero(double * arr, uint64_t count) {
  OMP("omp parallel for")
  for(uint64_t i = 0; i < count; i++) {
    arr[i] = 0;
  }
}

void
score_mean_variance_first_match(
    struct stinger * S, 
    uint64_t nv, 
    double volume, 
    double * scores, 
    uint64_t * matches, 
    double * sum, 
    double * sum_squares) 
{
  double local_sum = 0;
  double local_sum_squares = 0;

  /* precompute as much of the scores as possible */
  double half_volume = volume / 2;
  double two_div_vol_sqd = 2 / ((volume) * (volume));
  struct stinger_vb * LVA = S->LVA;

  /* for each vertex */
  OMP("omp parallel for reduction(+:local_sum) reduction(+:local_sum_squares)")
  for(uint64_t u = 0; u < nv; u++) {
    double best_found = 0;
    uint64_t best_match = matches[u];

    /* precompute more of the score (vtx specific part) */
    double wt_u_x2_div_vol_sqd = LVA[u].weight * two_div_vol_sqd;

    /* for all edge blocks of that vertex */
    STINGER_FORALL_EDGES_OF_VTX_BEGIN(S, u) {
      double wt_v = LVA[STINGER_EDGE_DEST].weight;
      double score = (((double)STINGER_EDGE_WEIGHT) / half_volume) - (wt_v * wt_u_x2_div_vol_sqd);
      /* check for the best score for this vertex*/
      if(score > best_found) {
	best_found = score;
	best_match = STINGER_EDGE_DEST;
      }
      /* sum the score and its square */
      local_sum += score;
      local_sum_squares += (score * score);
    } STINGER_FORALL_EDGES_OF_VTX_END();

    /* writeback the best score */
    scores[u] = best_found;
    matches[u] = best_match;
  }

  *sum = local_sum;
  *sum_squares = local_sum_squares;
}

void
filter_scores(
    uint64_t nv, 
    double * scores, 
    uint64_t * matches, 
    double sum, 
    double sum_squares) 
{
  /* compute the threshold */
  double mean = (sum) / nv;
  double mean_of_sq = (sum_squares) / nv;
  double variance = mean_of_sq - (mean * mean);
  double st_dev = sqrt(variance);
  double threshold = mean - 1.5 * st_dev;

  /* for all vertices */
  OMP("omp parallel for")
  for(uint64_t v = 0; v < nv; v++) {
    double score = scores[v];
    /* if score is below the threshold, don't match */
    if(score != 0 && score < threshold) {
      matches[v] = v;
    } else {
      /* if my match is trying to match on a lower-scoring edge, remove its match */
      uint64_t match = matches[v];
      if(scores[match] <= score) {
	matches[match] = match;
      }
    }
  }
}

void
tree_climb(
    uint64_t nv, 
    uint64_t * matches) 
{
  /* for all vertices */
  OMP("omp parallel for")
  for(uint64_t v = 0; v < nv; v++) {
    uint64_t older_match, old_match, match = v;
    old_match = v;
    /* climb the tree of matchings until we reach a root or cycle */
    do {
      older_match = old_match;
      old_match = match;
      match = matches[match];
      /* found a cycle - pick the lesser ID */
      if(match == older_match) {
	match = match > old_match ? old_match : match;
	break;
      }
    } while(old_match != match);
    matches[v] = match;
  }
}

int
multi_contract_root(
    struct stinger * S, 
    uint64_t nv,
    uint64_t * matches,
    int64_t timestamp)
{
  int work_remaining = 0;

  struct stinger_vb * LVA = S->LVA;
  struct stinger_eb * ebpool_priv = S->ebpool->ebpool;

  /* for all vertices */
  OMP("omp parallel for reduction(+:work_remaining)")
  for(uint64_t u = 0; u < nv; u++) {
    uint64_t match_u = matches[u];
    /* if it's a root */
    if(match_u == u) {
      /* for all edge blocks of the vertex */
      struct stinger_eb * currentBlock = ebpool_priv + LVA[u].edges;
      while(currentBlock != ebpool_priv) {
	struct stinger_edge * edge = currentBlock->edges;
	struct stinger_edge * last_edge = edge + currentBlock->high;
	/* for each edge in a block */
	for(; edge < last_edge; ++edge) {
	  int64_t v = edge->neighbor;
	  /* if the edge exists */
	  if(v >= 0) {
	    uint64_t match_v = matches[v];
	    /* and the edge is to be contracted, remove and increment vtx weight */
	    if(match_u == match_v) {
	      work_remaining = 1;
	      stinger_int64_fetch_add(&(LVA[match_u].weight), edge->weight);
	      edge->neighbor = ~v;
	      currentBlock->numEdges--;
	      stinger_int64_fetch_add(&(LVA[u].outDegree), -1); 
	      stinger_int64_fetch_add(&(LVA[v].inDegree), -1);
	    /* otherwise remove, remap, reinsert */
	    } else if(match_v != v) {
	      work_remaining = 1;
	      edge->neighbor = ~v;
	      currentBlock->numEdges--;
	      stinger_int64_fetch_add(&(LVA[u].outDegree), -1); 
	      stinger_int64_fetch_add(&(LVA[v].inDegree), -1);
	      stinger_incr_edge(S, currentBlock->etype, u, match_v, edge->weight, timestamp);
	    }
	  }
	}
	currentBlock = currentBlock->next + ebpool_priv;
      }
    }
  }

  return work_remaining;
}

int
multi_contract_tree(
    struct stinger * S, 
    uint64_t nv,
    uint64_t * matches,
    int64_t timestamp)
{
  int work_remaining = 0;

  struct stinger_vb * LVA = S->LVA;
  struct stinger_eb * ebpool_priv = S->ebpool->ebpool;

  /* for all vertices */
  OMP("omp parallel for reduction(+:work_remaining)")
  for(uint64_t u = 0; u < nv; u++) {
    uint64_t match_u = matches[u];
    struct stinger_eb * currentBlock = ebpool_priv + LVA[u].edges;
    if(match_u != u) {
      while(currentBlock != ebpool_priv) {
	struct stinger_edge * edge = currentBlock->edges;
	struct stinger_edge * last_edge = edge + currentBlock->high;
	/* for each edge in a block */
	for(; edge < last_edge; ++edge) {
	  int64_t v = edge->neighbor;
	  if(v >= 0) {
	    work_remaining = 1;
	    uint64_t match_v = matches[v];
	    if(match_u == match_v) {
	      stinger_int64_fetch_add(&(LVA[match_u].weight), edge->weight);
	    } else {
	      stinger_incr_edge(S, currentBlock->etype, match_u, match_v, edge->weight, timestamp);
	    }
	    edge->neighbor = ~v;
	    currentBlock->numEdges--;
	    stinger_int64_fetch_add(&(LVA[u].outDegree), -1); 
	    stinger_int64_fetch_add(&(LVA[v].inDegree), -1);
	  }
	}
	currentBlock = currentBlock->next + ebpool_priv;
      }
      stinger_int64_fetch_add(&S->LVA[match_u].weight, S->LVA[u].weight);
      S->LVA[u].weight = 0;
    }
  }

  return work_remaining;
}

double
modularity_score (
    struct stinger * S, 
    uint64_t * cluster, 
    uint64_t nv, 
    uint64_t ne)
{
  int64_t *Lss = (int64_t *)calloc(nv, sizeof(int64_t)), *Lsplus = (int64_t *)calloc(nv, sizeof(int64_t));
  double mod = 0.0, m_inv = 1.0/(double)ne;

  STINGER_FORALL_EDGES_BEGIN(S, 0) {
    if(STINGER_EDGE_SOURCE < STINGER_EDGE_DEST) {
      if(cluster[STINGER_EDGE_SOURCE] == cluster[STINGER_EDGE_DEST]) {
	Lss[cluster[STINGER_EDGE_SOURCE]]++;
      } else {
	Lsplus[cluster[STINGER_EDGE_SOURCE]]++;
      }
    }
  } STINGER_FORALL_EDGES_END();

  uint64_t total = 0;
  for(uint64_t i = 0; i < nv; i++) {
    if(cluster[i] >= nv)
      fprintf(stderr,"\nERROR: CLUSTER of %ld is %ld",i,cluster[i]);
    if(cluster[i] == i)
      total++;
    if(Lss[i] != 0.0 || Lsplus[i] != 0.0) {
      mod += ((double)Lss[i]) - (m_inv * ((double)Lsplus[i])) * (0.25 * ((double)Lsplus[i]));
    }
  }

  free(Lss); free(Lsplus);

  return mod * m_inv;
}

void
static_multi_contract_clustering (
    uint64_t ** matches,
    uint64_t nv,
    struct stinger * S,
    struct stinger * S_orig)
{

  double volume = 0;

  for(uint64_t t= 0; t < STINGER_NUMETYPES; t++) {
    volume += sum_all_edgeweights(S, t);
  }

  int work_remaining = 1;

  uint64_t iteration = 0;
  double sum, sum_squares;
  double * scores = xmalloc(nv * sizeof(double));
  *matches = xmalloc(nv * sizeof(uint64_t));

  OMP("omp parallel for")
  for(uint64_t v = 0; v < nv; v++) {
    stinger_set_vweight(S,v,0);
  }

  sequence(*matches, nv);

  int64_t max_iterations = 1000;
  if(getenv("MAX_ITERATIONS") && atoi(getenv("MAX_ITERATIONS")) > 0)
     max_iterations = atoi(getenv("MAX_ITERATIONS")) - 1;

  /* initialize the experiment */
  int ccexperiment = 0;
  int64_t scores_files = 0;
  const int64_t num_answers = 4;
  int64_t * num_in_answer, * comm_size, * distances, * inner_distances; 
  double * avg_distance, * avg_inner_distance;

  /* XXX hard code for now */
  char * red_team_node_names[] = {"1246400", "1657600", "1294800", "1683900"};
  int64_t red_team[4];

  if(getenv("CCEXPERIMENT")) {
    ccexperiment = 1;
    /* want to know: 
     * - how many in the same comm as the answer?
     * - what is the com size? 
     * - what is the distance from the answer, in comm, in graph? 
     */

    for(uint64_t i = 0; i < num_answers; i++) {
      red_team[i] = stinger_physID_to_vtx(S, red_team_node_names[i], strlen(red_team_node_names[i]));
    }

    if(getenv("SCORES") && getenv("SCORES_OUT")) {
      DIR * dp = opendir(getenv("SCORES"));
      struct dirent * ep;

      if(dp) {
	while(ep = readdir(dp)) {
	  if(ep->d_type == DT_REG)
	    scores_files++;
	}
      }

      num_in_answer = xcalloc(scores_files * num_answers * max_iterations, sizeof(int64_t));
      comm_size = xcalloc(max_iterations * num_answers , sizeof(int64_t));
      avg_distance = xcalloc(scores_files * num_answers * max_iterations, sizeof(double));
      avg_inner_distance = xcalloc(scores_files * num_answers * max_iterations, sizeof(double));
      distances = xcalloc(num_answers * nv, sizeof(int64_t));
      inner_distances = xcalloc(num_answers * nv, sizeof(int64_t));

      for(uint64_t i = 0; i < num_answers; i++) {
	if(red_team[i] > 0) {
	  bfs(S, nv, red_team[i], distances + i * nv);
	}
      }
    }
  }

  tic();
  while(work_remaining && iteration < max_iterations) {
    zero(scores, nv);
    work_remaining = 0;
    sum = 0;
    sum_squares = 0;

    score_mean_variance_first_match(S, nv, volume, scores, *matches, &sum, &sum_squares);

    filter_scores(nv, scores, *matches, sum, sum_squares);

    tree_climb(nv, *matches); 

    work_remaining += multi_contract_root(S, nv, *matches, 0);

    work_remaining += multi_contract_tree(S, nv, *matches, 0);

    /* update experiment data */
    if(ccexperiment) {
      OMP("omp parallel for")
      for(uint64_t v = 0; v < nv; v++) {
	for(uint64_t u = 0; u < num_answers; u++) {
	  if(red_team[u] >= 0 && (*matches)[v] == (*matches)[red_team[u]]) {
	    comm_size[iteration * num_answers + u]++;
	  }
	}
      }

      for(uint64_t i = 0; i < num_answers; i++) {
	bfs_in_set(S_orig, nv, red_team[i], inner_distances + i * nv, *matches, (*matches)[red_team[i]]);
      }

      if(getenv("SCORES") && getenv("SCORES_OUT")) {
	DIR * dp = opendir(getenv("SCORES"));
	FILE * fp = NULL;
	struct dirent * ep;
	
	if(dp) {
	  int cur_file = 0;

	  while(ep = readdir(dp)) {

	    if(ep->d_type == DT_REG) {
	      char filename[1024];
	      sprintf(filename, "%s%s", getenv("SCORES"), ep->d_name);
	      fp = fopen(filename, "r");

	      if(fp) {
		char * line = malloc(sizeof(char) * 20);
		int line_size = 20;
		int max = 100;
		while((getline(&line, &line_size, fp) > -1) && (max-- > 0)) {
		  int64_t vtx = stinger_physID_to_vtx(S, line, strlen(line)-1);

		  if(vtx != -1) {
		    for(uint64_t i = 0; i < num_answers; i++) {
		      if((red_team[i] > -1) && ((*matches)[vtx] == (*matches)[red_team[i]])) {
			num_in_answer[scores_files * num_answers * iteration + num_answers * cur_file + i]++;
			avg_distance[scores_files * num_answers * iteration + num_answers * cur_file + i] += distances[i * nv + vtx];
			//if(inner_distances[i * nv + vtx] != INT64_MAX)
			  avg_inner_distance[scores_files * num_answers * iteration + num_answers * cur_file + i] += inner_distances[i * nv + vtx];
		      }
		    }
		  }
		}

		for(uint64_t i = 0; i < 4; i++) {
		    if(num_in_answer[scores_files * num_answers * iteration + num_answers * cur_file + i]) {
		      avg_distance[scores_files * num_answers * iteration + num_answers * cur_file + i] /= num_in_answer[scores_files * num_answers * iteration + num_answers * cur_file + i];
		      avg_inner_distance[scores_files * num_answers * iteration + num_answers * cur_file + i] /= num_in_answer[scores_files * num_answers * iteration + num_answers * cur_file + i];
		    }
		}

		fclose(fp);
	      } else {
		fprintf(stderr, "Error opening scores output file %s\n", filename);
	      }
	      cur_file++;
	    }
	  }
	  closedir(dp);
	} else {
	  fprintf(stderr, "Error opening scores\n");
	}
      }
    }
    iteration++;
  }
  printf("\tcluster_time %lf seconds\n", toc());

  if(ccexperiment) {
    if(getenv("SCORES") && getenv("SCORES_OUT")) {
      DIR * dp = opendir(getenv("SCORES"));
      FILE * fp = NULL;
      struct dirent * ep;
      
      if(dp) {
	int cur_file = 0;
	while(ep = readdir(dp)) {

	  if(ep->d_type == DT_REG) {
	    char filename[1024];
	    sprintf(filename, "%s%s", getenv("SCORES_OUT"), ep->d_name);
	    fp = fopen(filename, "w");

	    if(fp) {
	      for(uint64_t i = 0; i < num_answers; i++) {
		fprintf(fp,"%s,num,", red_team_node_names[i]);
		for(uint64_t j = 0; j < iteration; j++) {
		  fprintf(fp,"%ld,", num_in_answer[scores_files * num_answers * j + num_answers * cur_file + i]);
		}
		fprintf(fp,"\n%s,size,", red_team_node_names[i]);
		for(uint64_t j = 0; j < iteration; j++) {
		  fprintf(fp,"%ld,", comm_size[num_answers * j + i]);
		}
		fprintf(fp,"\n%s,dist,", red_team_node_names[i]);
		for(uint64_t j = 0; j < iteration; j++) {
		  fprintf(fp,"%lf,",avg_distance[scores_files * num_answers * j + num_answers * cur_file + i]);
		}
		fprintf(fp,"\n%s,inner,", red_team_node_names[i]);
		for(uint64_t j = 0; j < iteration; j++) {
		  fprintf(fp,"%lf,",avg_inner_distance[scores_files * num_answers * j + num_answers * cur_file + i]);
		}
		fprintf(fp,"\n");
	      }

	      fclose(fp);
	    } else {
	      fprintf(stderr, "Error opening scores output file %s\n", filename);
	    }
	    cur_file++;
	  }
	}
	closedir(dp);
      } else {
	fprintf(stderr, "Error opening scores\n");
      }
    }
  }

  uint64_t count = 0;
  uint64_t edges = 0;
  for(uint64_t v = 0; v < nv; v++) {
    if(stinger_outdegree(S, v)) {
      edges += stinger_outdegree(S, v);
      if((*matches)[v] == v)
	count++;
    }
  }

  free(scores);

  printf("\tclusters %ld\n", count);
  printf("\tclusters_edges %ld\n", edges);
  printf("\tclustering_iterations %ld\n", iteration);
}

uint64_t connected_components (struct stinger * S, int64_t * components, uint64_t nv) {
  OMP("omp parallel for") 
  for(uint64_t v = 0; v < nv; v++) {
    components[v] = v;
  }

  while(1) {
    int done = 1;

    for(uint64_t t = 0; t < STINGER_NUMETYPES; t++) {
      STINGER_PARALLEL_FORALL_EDGES_BEGIN(S, t) {
	if(components[STINGER_EDGE_SOURCE] < components[STINGER_EDGE_DEST]) {
	  components[STINGER_EDGE_DEST] = components[STINGER_EDGE_SOURCE];
	  done = 0;
	} else if(components[STINGER_EDGE_DEST] < components[STINGER_EDGE_SOURCE]) {
	  components[STINGER_EDGE_SOURCE] = components[STINGER_EDGE_DEST];
	  done = 0;
	}
      } STINGER_PARALLEL_FORALL_EDGES_END();
    }

    if(done)
      break;

    OMP("omp parallel for")
    for(uint64_t v = 0; v < nv; v++) {
      while(components[v] != components[components[v]]) {
	components[v] = components[components[v]];
      }
    }
  }

  uint64_t num_components = 0;
  OMP("omp parallel for reduction(+:num_components)")
  for(uint64_t v = 0; v < nv; v++) {
    if(components[v] == v && stinger_outdegree(S,v))
      num_components++;
  }

  return num_components;
}

