#include <stdio.h>
#include <string.h>
#include "evicomb.h"
#include "evicomb-internal.h"
#include "stinger.h"
#include "stinger-physmap-new.h"

void sort_results (int i_start, int i_end, double *beliefs, int *ids) {
  if (i_start >= i_end)
    return;

  int left,
    right,
    temp_int;
  double pivot,
    temp_dbl;

  left = i_start;
  right = i_end;
  pivot = beliefs[(left << 1) + STATE_ANOMALOUS];
  
  while (left < right) {
    while ((beliefs[(right << 1) + STATE_ANOMALOUS] <= pivot) && (left < right))
      right--;
    while ((beliefs[(left << 1) + STATE_ANOMALOUS] >= pivot) && (left < right))
      left++;
    if (left != right) {
      temp_dbl = beliefs[(right << 1) + STATE_NORMAL];
      beliefs[(right << 1) + STATE_NORMAL] = beliefs[(left << 1) + STATE_NORMAL];
      beliefs[(left << 1) + STATE_NORMAL] = temp_dbl;
      temp_dbl = beliefs[(right << 1) + STATE_ANOMALOUS];
      beliefs[(right << 1) + STATE_ANOMALOUS] = beliefs[(left << 1) + STATE_ANOMALOUS];
        beliefs[(left << 1) + STATE_ANOMALOUS] = temp_dbl;
      temp_int = ids[left];
      ids[left] = ids[right];
      ids[right] = temp_int;
    }
    temp_dbl = beliefs[(i_start << 1) + STATE_NORMAL];
    beliefs[(i_start << 1) + STATE_NORMAL] = beliefs[(left << 1) + STATE_NORMAL];
    beliefs[(left << 1) + STATE_NORMAL] = temp_dbl;
    temp_dbl = beliefs[(i_start << 1) + STATE_ANOMALOUS];
    beliefs[(i_start << 1) + STATE_ANOMALOUS] = beliefs[(left << 1) + STATE_ANOMALOUS];
    beliefs[(left << 1) + STATE_ANOMALOUS] = temp_dbl;
    temp_int = ids[left];
    ids[left] = ids[i_start];
    ids[i_start] = temp_int;

    sort_results(i_start, left, beliefs, ids);
    sort_results(left + 1, i_end, beliefs, ids);
  }

  return;
}

int cmp_int (const void *p1, const void *p2) {
  return *((int *)p1) - *((int *)p2);
}

int bsrch (const int *arr, int size, int key) {
  int low,
    mid,
    high;

  low = 0;
  high = size - 1;
  
  while (low <= high) {
    mid = (low + high) >> 1;
    if (arr[mid] < key) {
      low = mid + 1;
    } else if (arr[mid] > key) {
      high = mid - 1;
    } else {

      return 1;
    }
  }

  return 0;
}

void evicomb (struct stinger *s, int anomaly_count, int *anomaly_ids, char *params, char *path) {
 
  int i,
    j,
    k,
    flag_conv,
    iter_count,
    ver_count,
    *sorted_ids;
  double temp_dbl,
    ini_bel_norm,
    ini_bel_anom,
    *m1,
    *m2,
    prop_mtx[STATE_COUNT][STATE_COUNT],
    *beliefs,
    *m_cur,
    *m_next;
  char output_file_name[1000],
    output_file_name_tmp[1000],
    output_file_name_ids[1000],
    output_file_name_rep[1000];
  FILE *output_file,
    *output_file_ids,
    *output_file_rep;
  int64_t buf[100];
  size_t neig_count;

  printf("params: %s anomaly_count: %d\n", params, anomaly_count);
  for(int i = 0; i < anomaly_count; i++) {
    printf("\tid[%d] = %d\n", i, anomaly_ids[i]);
  }

  /* parsing the parameters */
  sscanf(params, "%s %s %s %lf %lf %lf %lf ", output_file_name, output_file_name_ids, output_file_name_rep, &prop_mtx[STATE_NORMAL][STATE_ANOMALOUS], &prop_mtx[STATE_ANOMALOUS][STATE_ANOMALOUS], &ini_bel_norm, &ini_bel_anom);
  /* initialization */
  prop_mtx[STATE_NORMAL][STATE_NORMAL] = 1 - prop_mtx[STATE_NORMAL][STATE_ANOMALOUS];
  prop_mtx[STATE_ANOMALOUS][STATE_NORMAL] = 1 - prop_mtx[STATE_ANOMALOUS][STATE_ANOMALOUS];
  flag_conv = 0;
  iter_count = 0;
  ver_count = STINGER_MAX_LVERTICES; 
  /* allocate space */
  beliefs = (double *) calloc (ver_count, sizeof(double) << 1);
  m1 = (double *) calloc (ver_count, sizeof(double) << 1);
  m2 = (double *) calloc (ver_count, sizeof(double) << 1);
  sorted_ids = (int *) calloc (ver_count, sizeof(int));
  /* initialize beliefs */
  for (i = 0; i < ver_count; i++) {    
    beliefs[(i << 1) + STATE_ANOMALOUS] = ini_bel_norm;
    beliefs[(i << 1) + STATE_NORMAL] = 1 - ini_bel_norm;
  }
  for (i = 0; i < anomaly_count; i++) {
    beliefs[(anomaly_ids[i] << 1) + STATE_ANOMALOUS] = ini_bel_anom;
    beliefs[(anomaly_ids[i] << 1) + STATE_NORMAL] = 1 - ini_bel_anom;
  }
  /* initialize messages */
  m_cur = (double *)m1;
  m_next = (double *)m2;
  for (i = 0; i < ver_count; i++) {
    m_cur[(i << 1) + STATE_NORMAL] = 
      prop_mtx[STATE_NORMAL][STATE_NORMAL] * beliefs[(i << 1) + STATE_NORMAL] +
      prop_mtx[STATE_ANOMALOUS][STATE_NORMAL] * beliefs[(i << 1) + STATE_ANOMALOUS];
    m_cur[(i << 1) + STATE_ANOMALOUS] = 
      prop_mtx[STATE_NORMAL][STATE_ANOMALOUS] * beliefs[(i << 1) + STATE_NORMAL] +
      prop_mtx[STATE_ANOMALOUS][STATE_ANOMALOUS] * beliefs[(i << 1) +  STATE_ANOMALOUS];
  }
  /* propogate the beliefs until convergence achieved */
  while(!flag_conv) {
    /* calculate messages */
    for (i = 0; i < ver_count; i++) {
    m_cur[(i << 1) + STATE_NORMAL] = 
      prop_mtx[STATE_NORMAL][STATE_NORMAL] * beliefs[(i << 1) + STATE_NORMAL] +
      prop_mtx[STATE_ANOMALOUS][STATE_NORMAL] * beliefs[(i << 1) + STATE_ANOMALOUS];
    m_cur[(i << 1) + STATE_ANOMALOUS] = 
      prop_mtx[STATE_NORMAL][STATE_ANOMALOUS] * beliefs[(i << 1) + STATE_NORMAL] +
      prop_mtx[STATE_ANOMALOUS][STATE_ANOMALOUS] * beliefs[(i << 1) + STATE_ANOMALOUS];
    }
    /* perform belief updates */
    for (i = 0; i < ver_count; i++) {
      beliefs[(i << 1) + STATE_NORMAL] = 1;
      beliefs[(i << 1) + STATE_ANOMALOUS] = 1;
      STINGER_FORALL_EDGES_OF_VTX_BEGIN(s, i) {
        beliefs[(i << 1) + STATE_NORMAL] *= m_cur[(STINGER_EDGE_DEST << 1) + STATE_NORMAL];
        beliefs[(i << 1) + STATE_ANOMALOUS] *= m_cur[(STINGER_EDGE_DEST << 1) + STATE_ANOMALOUS];
      } STINGER_FORALL_EDGES_OF_VTX_END();
      temp_dbl = beliefs[(i << 1) + STATE_NORMAL] + beliefs[(i << 1) + STATE_ANOMALOUS];
      beliefs[(i << 1) + STATE_NORMAL] /= temp_dbl;
      beliefs[(i << 1) + STATE_ANOMALOUS] /= temp_dbl;
    }
    /* check if significant change occurs, if so update flag */
    if (iter_count >= MAX_ITERS)
      flag_conv = 1;
    else
      iter_count++;
  }

  /* hard code for now */
  uint64_t red_team_node_size = 4;
  char * red_team_node_names[] = {"1246400", "1657600", "1294800", "1683900"};
  uint8_t * friends_with_foe = calloc(STINGER_MAX_LVERTICES, sizeof(uint8_t));

  for(uint64_t i = 0; i < red_team_node_size; i++) {
    int64_t vtx = stinger_physID_to_vtx(s, red_team_node_names[i], strlen(red_team_node_names[i]));
    if(vtx >= 0) {
      friends_with_foe[vtx] = 2;
      STINGER_FORALL_EDGES_OF_VTX_BEGIN(s, vtx) {
	friends_with_foe[STINGER_EDGE_DEST] = 1;
      } STINGER_FORALL_EDGES_OF_VTX_END();
    }
  }

  /* output results */
  for (i = 0; i < ver_count; i++)
    sorted_ids[i] = i;
  sort_results(0, ver_count - 1, beliefs, sorted_ids);
  qsort(anomaly_ids, anomaly_count, sizeof(int), cmp_int);
  sprintf(output_file_name_tmp, "%s%s", path, output_file_name);
  output_file = fopen(output_file_name_tmp, "w");
  sprintf(output_file_name_tmp, "%s%s", path, output_file_name_ids);
  output_file_ids = fopen(output_file_name_tmp, "w");
  sprintf(output_file_name_tmp, "%s%s", path, output_file_name_rep);
  output_file_rep = fopen(output_file_name_tmp, "w");
  fprintf(output_file, "%s\n", params);
  fprintf(output_file_ids, "%s\n", params);
  fprintf(output_file_rep, "%s\n", params);
  for (i = 0; i < ver_count; i++) {
    fprintf(output_file, "%.3lf\n", beliefs[(i << 1) + STATE_ANOMALOUS]);   // score of ranked node i
    fprintf(output_file_ids, "%d\n", sorted_ids[i]);                        // id of ranked node i
    fprintf(output_file_rep, "%d %d\n", bsrch(anomaly_ids, anomaly_count, sorted_ids[i]), (int) friends_with_foe[i]); // whether ranked node i is a reported anomaly
  }
  fclose (output_file);
  fclose (output_file_ids);
  fclose (output_file_rep);
  /* deallocate space */
  free (beliefs);
  free (m1);
  free (m2);
  free (sorted_ids);
  
  return;
}

