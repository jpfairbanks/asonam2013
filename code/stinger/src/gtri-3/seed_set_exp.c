#include "stinger-atomics.h"
#include "stinger-utils.h"
#include "stinger.h"
#include "xmalloc.h"
#include "timer.h"
#include "seed_set_exp.h"
static int64_t nv, ne, naction;
static int64_t *restrict off;
static int64_t *restrict from;
static int64_t *restrict ind;
static int64_t *restrict weight;
static int64_t *restrict action;
static struct stinger *S;
/* handles for I/O memory */
static int64_t *restrict graphmem;
static int64_t *restrict actionmem;

static char *initial_graph_name = INITIAL_GRAPH_NAME_DEFAULT;
static char *action_stream_name = ACTION_STREAM_NAME_DEFAULT;

static long batch_size = BATCH_SIZE_DEFAULT;
static long nbatch = NBATCH_DEFAULT;

static double Modularity;
static int64_t Count;
static double updatetime;


#define ACTI(k) (action[2*(k)])
#define ACTJ(k) (action[2*(k)+1])
#define ACTI2(k) (act[2*(k)])
#define ACTJ2(k) (act[2*(k)+1])



#if 0
int main (const int argc, char *argv[])
{

 parse_args (argc, argv, &initial_graph_name, &action_stream_name,
              &batch_size, &nbatch);
  STATS_INIT ();

  load_graph_and_action_stream (initial_graph_name, &nv, &ne,
                                (int64_t **) & off, (int64_t **) & ind,
                                (int64_t **) & weight,
                                (int64_t **) & graphmem, action_stream_name,
                                &naction, (int64_t **) & action,
                                (int64_t **) & actionmem);

  print_initial_graph_stats (nv, ne, batch_size, nbatch, naction);
  BATCH_SIZE_CHECK ();

  int64_t num_seeds=4;
  //int64_t *seeds;
  int64_t seeds[4]={40,41,46,48};
 
   
  tic ();
  S=stinger_new();
  stinger_set_initial_edges (S, nv, 0, off, ind, weight, NULL, NULL, 0);
  PRINT_STAT_DOUBLE ("time_stinger", toc ());
  fflush (stdout);
    
    STINGER_FORALL_EDGES_OF_VTX_BEGIN(S,40){
       const int64_t dest= STINGER_EDGE_DEST;
   //     printf("%ld,\n",dest);
       }STINGER_FORALL_EDGES_OF_VTX_END();



  int64_t total_comm_edges=0;
  int64_t internal_comm_edges=0;
  int64_t number_neighbors=0;
  int64_t nedges=stinger_total_edges(S);
  int64_t *membership;
 int64_t *neighbors;
 int64_t *neighbors2;
  double *neighbor_modval;
  int64_t *edge_between;
  int64_t *affected;
  int64_t i;
  int64_t numaffected=0;
  Modularity=0;
  Count=0;
  updatetime=0;
  membership=(int64_t*)xmalloc(nv*sizeof(int64_t));

  neighbors=(int64_t*)xmalloc(nv*sizeof(int64_t));
  neighbors2=(int64_t*)xmalloc(nv*sizeof(int64_t));
  edge_between=(int64_t*)xmalloc(nv*sizeof(int64_t));
  neighbor_modval=(double*)xmalloc(nv*sizeof(double));
  affected=(int64_t*)xmalloc(nv*sizeof(int64_t));
  double history[HISTORY];
  int64_t history_idx=0;;
  tic();


  seed_set_expansion(S,seeds,num_seeds,nv,nedges,membership,neighbors, edge_between, neighbor_modval,&total_comm_edges, &internal_comm_edges,&number_neighbors,history,&history_idx);
  PRINT_STAT_DOUBLE("seed set time",toc());

//stream updates
/*
  for (i=0;i<nv;i++){
    neighbors2[i]=0;
  }
  for (i=0;i<number_neighbors;i++){
    neighbors2[neighbors[i]]=edge_between[i];
  }

 
int64_t ntrace = 0;
  if (1 == batch_size) {
    int actk;
    for (actk = 0; actk < nbatch; ++actk, ++ntrace) {
   
      int changed;

      int64_t i = ACTI (actk);
      int64_t j = ACTJ (actk);

      if (i < 0) {
        // Delete. 
        i = -i - 1;
        j = -j - 1;

        assert (i < nv);
        assert (j < nv);
        assert (i >= 0);
        assert (j >= 0);

        changed = stinger_remove_edge_pair (S, 0, i, j);
	tic();
        if (changed && i != j) {
	    if(membership[i]== membership[j]){
	    affected[0]=i;
	    affected[1]=j;
	    numaffected=2;
	    update_community(S,seeds,num_seeds,nv,nedges,membership,neighbors2,&total_comm_edges,&internal_comm_edges,affected,numaffected,history,&history_idx);

	    }

	  //  seed_set_expansion(S,seeds,num_seeds,nv,nedges,membership,neighbors, edge_between, neighbor_modval,&total_comm_edges, &internal_comm_edges,&number_neighbors,history,&history_idx);
	     updatetime+=toc();
        }
      } else {                    // Add an edge. 
        assert (i < nv);
        assert (j < nv);
        assert (i >= 0);
        assert (j >= 0);

        changed = stinger_insert_edge_pair (S, 0, i, j, 1, ntrace + 1);
	tic();
        if (changed && i != j) {
	   if (membership[i]!=membership[j]){
	    affected[0]=i;
	    affected[1]=j;
	    numaffected=2;
	    update_community(S,seeds,num_seeds,nv,nedges,membership,neighbors2,&total_comm_edges,&internal_comm_edges,affected,numaffected,history,&history_idx);
	    }
				
	   //   seed_set_expansion(S,seeds,num_seeds,nv,nedges,membership,neighbors, edge_between, neighbor_modval,&total_comm_edges, &internal_comm_edges,&number_neighbors,history,&history_idx);
	    updatetime+=toc();
        }
      }
    }
  }
else {
   
  int64_t * act = xmalloc (2 * 2 * batch_size * sizeof(*act));
  int64_t * insoff = xmalloc ((2 * 2 * batch_size + 1) * sizeof(*insoff));
  int64_t * deloff = xmalloc ((2 * 2 * batch_size + 1) * sizeof(*deloff));

    for (int actno = 0; actno < nbatch * batch_size; actno += batch_size) {
      const int endact = (actno + batch_size > naction ?
                          naction : actno + batch_size);
    
      numaffected=0;  
int64_t N = stinger_sort_actions (endact - actno, &action[2*actno], insoff, deloff, act);
 tic();
 
 int64_t k;
    for ( k = actno; k < endact; k++) {
      const int64_t i = ACTI(k);
      const int64_t j = ACTJ(k);
      
      
      if(i<0){
if (membership[-i-1]!=membership[-j-1]){
	    affected[numaffected++]=-i-1;
	    affected[numaffected++]=-j-1;
 
 }
      }
      else if (membership[i]==membership[j]){
	  affected[numaffected++]=i;
	  affected[numaffected++]=j;
	}
      
    }
 

     stinger_remove_and_insert_batch (S, 0, actno+1, N, insoff, deloff, act);
      update_community(S,seeds,num_seeds,nv,nedges,membership,neighbors2,&total_comm_edges,&internal_comm_edges,affected,numaffected,history,&history_idx);
      
       updatetime+=toc();

    }
    free(act);
    free(insoff);
    free(deloff);
 }
*/


 stinger_free_all (S);
 free (graphmem);
 free (actionmem); 
  STATS_END ();


}
#endif



/*
 * @params
 *
 *    S         - is the stinger structure
 *    seeds     - is an array of the seed vertices,
 *    num_seeds  - is the number of seed vertices,
 *    nv        - is the number of vertices,
 *    ne        - is the number of ed 
 *   membership        - unitialized array of size nv - will contain 1 if in commonuty, 0 otherwise
 *  --these are initialized outside of the function to allow updates in streaming--
 *   neighbors   - unitialized array of size nv 
 *   edge_between   - unitialized array of size nv
 *   neighbor_modval   - unitialized array of size nv
 *   total_comm_edges     -int64_t set to 0
 *   internal_comm_edges   -int64_t set to 0
 *   number_neighbors   -int64_t set to 0
 *   history         - array of size HISTORY;
 *   history_idx   - int64_t* set to 0
 */
void seed_set_expansion(const struct stinger *S, int64_t *seeds, int64_t num_seeds, int64_t nv, int64_t ne, int64_t *membership, int64_t *neighbors, int64_t *edge_between, double * neighbor_modval, int64_t *total_community_edges, int64_t *internal_community_edges, int64_t *number_neighbors,double *history, int64_t *history_idx){

 
  int64_t num_iterations=0;
  int64_t total_comm_edges=0;
  int64_t internal_comm_edges=0;
  int64_t num_neighbors=0;
  int64_t i,position;
  int64_t pos=0;
  double average_mod=0;
  double std_dev_mod=0;
  double modularity=0;
OMP("omp parallel for")
  for (i=0;i<nv;i++){
    membership[i]=0;
  }
OMP("omp parallel for")
  for (i=0;i<num_seeds;i++){
    membership[seeds[i]]=1;
    int64_t myinternaledges=0;
    int64_t d;
    int64_t deg=stinger_outdegree(S,seeds[i]);
    int64_t mypos=stinger_int64_fetch_add(&pos,deg);
    stinger_gather_typed_successors(S,0,seeds[i],&d,&neighbors[mypos],deg);
    
    num_neighbors+=deg;
    total_comm_edges+=deg;
    int64_t j;
    for (j=0;j<deg;j++){
      if (membership[neighbors[mypos+j]]==1){
	myinternaledges++;
      }
    }
    stinger_int64_fetch_add(&internal_comm_edges,myinternaledges);
  }
OMP("omp parallel for")
  for (i=0;i<num_neighbors;i++){
    edge_between[i]=1;
  }
OMP("omp parallel for")
  for (i=num_neighbors;i<nv;i++){
    edge_between[i]=0;
  }

  num_neighbors=contract(neighbors,edge_between,num_neighbors);

  
  int64_t total_num_merges=0;
  int64_t history_index=0;
  double history_average=0;
  double history_std_dev=0;
  int64_t change=1;
  int64_t maxnumedges=0;
 
 
  while (change==1){
    change=0;
     modularity=((double)internal_comm_edges-(((double)total_comm_edges)*(double)total_comm_edges)/((double)(4*ne)))/(double)ne;

   
    OMP("omp parallel for")
      for (i=0;i<num_neighbors;i++){
	neighbor_modval[i]=(double)edge_between[i]/(double)ne -((double)stinger_outdegree(S,neighbors[i])/(double)(2*ne))*((double)total_comm_edges/(double)(2*ne));
      }
    maxnumedges=0;
    for (i=0;i<num_neighbors;i++){
      average_mod=average_mod+neighbor_modval[i];
      if (edge_between[i]>maxnumedges){
	maxnumedges=edge_between[i];
      }
    }
    average_mod=average_mod/num_neighbors;
    for (i=0;i<num_neighbors;i++){
      std_dev_mod=std_dev_mod+(neighbor_modval[i]-average_mod)*(neighbor_modval[i]-average_mod);
    }
    std_dev_mod=std_dev_mod/num_neighbors;
    std_dev_mod=sqrt(std_dev_mod);
    
    if (total_num_merges>=HISTORY){
      history_average=0;
      history_std_dev=0;
      int64_t j;
      for (j=0;j<HISTORY;j++){
	history_average+=history[j];
      }
      history_average=history_average/(double)HISTORY;
      for (j=0;j<HISTORY;j++){
	history_std_dev=history_std_dev+(history[j]-history_average)*(history[j]-history_average);
      }
      history_std_dev=history_std_dev/(double)HISTORY;
    }
      
   
      int64_t *to_merge;
      int64_t *to_merge_edge_between;
      int64_t *new_neighbors;
      int64_t *new_neighbors_edge_between;
	
      to_merge=(int64_t*)xmalloc(num_neighbors*sizeof(int64_t));
      to_merge_edge_between=(int64_t*)xmalloc(num_neighbors*sizeof(int64_t));
      int64_t num_to_merge=0;
	
      int64_t j=num_neighbors-1;
       if (num_iterations<3 || total_num_merges<HISTORY){
	  if ( edge_between[j]==maxnumedges && neighbor_modval[j]>0){
	    while( edge_between[j]==maxnumedges && neighbor_modval[j]>0){
	      to_merge[num_to_merge]=neighbors[j];
	      to_merge_edge_between[num_to_merge]=edge_between[j];
	      int64_t curr_history_index=stinger_int64_fetch_add(&history_index,1)%HISTORY;
	      history[curr_history_index]=neighbor_modval[j];
	      j--;
	      num_to_merge++;
	      num_neighbors--;
	    }
	  }
	  int64_t lower,upper=j;
	  while(j>=0){
	    upper=j;
	    while( edge_between[j]<maxnumedges || neighbor_modval<=0){
	      j--;
	    }
	    lower=j+1;
	    while(edge_between[j]==maxnumedges && neighbor_modval>0){
	      to_merge[num_to_merge]=neighbors[j];
	      to_merge_edge_between[num_to_merge]=edge_between[j];
	      int64_t curr_history_index=stinger_int64_fetch_add(&history_index,1)%HISTORY;
	      history[curr_history_index]=neighbor_modval[j];
	      num_to_merge++;
	      num_neighbors--;
	      j--;
	    }
	    memmove(neighbors+j+1,neighbors+lower,(upper-lower+1)*sizeof(int64_t));
	    memmove(edge_between+j+1,edge_between+lower,(upper-lower+1)*sizeof(int64_t));
	    memmove(neighbor_modval+j+1,neighbor_modval+lower,(upper-lower+1)*sizeof(double));
	  }
       }
       else{
		
	 if (neighbor_modval[num_neighbors-1]>=max(0,max(average_mod+0*std_dev_mod,history_average-history_std_dev))){
	   while(neighbor_modval[j]>=max(0,max(average_mod+0*std_dev_mod,history_average-history_std_dev))){
	     to_merge[num_to_merge]=neighbors[j];
	     to_merge_edge_between[num_to_merge]=edge_between[j];
	     int64_t curr_history_index=stinger_int64_fetch_add(&history_index,1)%HISTORY;
	     history[curr_history_index]=neighbor_modval[j];
	     j--;
	     num_to_merge++;
	     num_neighbors--;
	   }
	 }
	 int64_t lower,upper=j;
	 while(j>=0){
	   upper=j;
	   while(j>=0 && neighbor_modval[j]<max(0,max(average_mod+0*std_dev_mod,history_average-history_std_dev))){
	     j--;
	   }
	   lower=j+1;
	   while(j>=0 && neighbor_modval[j]>=max(0,max(average_mod+0*std_dev_mod,history_average-history_std_dev))){
	     to_merge[num_to_merge]=neighbors[j];
	     to_merge_edge_between[num_to_merge]=edge_between[j];
	     int64_t curr_history_index=stinger_int64_fetch_add(&history_index,1)%HISTORY;
	     history[curr_history_index]=neighbor_modval[j];
	     num_to_merge++;
	     num_neighbors--;
	     j--;
	   }
	   memmove(neighbors+j+1,neighbors+lower,(upper-lower+1)*sizeof(int64_t));
	   memmove(edge_between+j+1,edge_between+lower,(upper-lower+1)*sizeof(int64_t));
	   memmove(neighbor_modval+j+1,neighbor_modval+lower,(upper-lower+1)*sizeof(double));
	 }
       }
	  
       total_num_merges+=num_to_merge;
       if (num_to_merge==0){
	 change=0;
       }
	  
       for (i=0;i<num_to_merge;i++){
	 int64_t currvertex=to_merge[i];
	 membership[currvertex]=1;
	 total_comm_edges+=stinger_outdegree(S,currvertex);
       }
       for (i=0;i<num_to_merge;i++){	    
	 internal_comm_edges+=to_merge_edge_between[i];
       }
       new_neighbors=(int64_t*)xmalloc(100*nv*sizeof(int64_t));
       position=0;

       OMP("omp parallel for")
       for (i=0;i<num_to_merge;i++){	    
	 int64_t outdeg=stinger_outdegree(S,to_merge[i]);
	 int64_t *tmp_neigh;
	 tmp_neigh=(int64_t*)xmalloc(outdeg*sizeof(int64_t));
	 int64_t count=0;
	 STINGER_FORALL_EDGES_OF_VTX_BEGIN(S,to_merge[i]){
	   const int64_t dest= STINGER_EDGE_DEST;
	   if (membership[dest]==0){		
	     tmp_neigh[count]=dest;
	     count++;
	   }
	 }STINGER_FORALL_EDGES_OF_VTX_END();
	    
	 int64_t y;
	 int64_t mypos=stinger_int64_fetch_add(&position,count);
	 for (y=0;y<count;y++){
	   new_neighbors[mypos+y]=tmp_neigh[y];
	 }
	 free(tmp_neigh);
       }
       new_neighbors_edge_between=(int64_t*)xmalloc(position*sizeof(int64_t));
       OMP("omp parallel for")
       for (i=0;i<position;i++){
	 new_neighbors_edge_between[i]=0;
       }
	
       int64_t num_new_neigh=contract(new_neighbors,new_neighbors_edge_between,position);
       num_neighbors=merge(neighbors,edge_between,new_neighbors,new_neighbors_edge_between,num_neighbors,num_new_neigh);
	  
       free(to_merge);
       free(to_merge_edge_between);
       free(new_neighbors);
       free(new_neighbors_edge_between);
       if (num_to_merge>0){
	 change=1;
       }
       num_iterations++;
  }
  *total_community_edges=total_comm_edges;
  *internal_community_edges=internal_comm_edges;
  *number_neighbors=num_neighbors;
  *history_idx=history_index;
  //Modularity+=modularity;
  //Count++;
}

int compare_func(const void *c1, const void *c2){
  const int64_t *v1,*v2;
  v1=(const int64_t *)c1;
  v2=(const int64_t *)c2;
  return (*v1-*v2);
}


int64_t contract(int64_t *neigh_array,int64_t* edge_between, int64_t length){
  qsort(neigh_array,length,sizeof(int64_t),compare_func);

  int64_t i=1,j=1,lastrepeat=-1,lastval=neigh_array[0],num_elements=length;
  int64_t curredges=1;
  while(j<length){
    if(neigh_array[i]==lastval){
      lastrepeat=i;
      while(neigh_array[i]==lastval){
	curredges++;
	i++;
	j++;
      }
      int64_t k;
      for (k=0;k<(length-j);k++){
	neigh_array[lastrepeat+k]=neigh_array[i+k];
      }
   
      num_elements=num_elements-(i-lastrepeat);
      i=lastrepeat;
      
    }
    lastval=neigh_array[i];
    edge_between[i-1]=curredges;
    curredges=1;
    i++;
    j++;

  }
  return num_elements;
}


int64_t merge(int64_t *neighbors, int64_t *edge_between, int64_t *temp_neighbors, int64_t *temp_edge_between, int64_t length1, int64_t length2){
  int64_t i=0,j=0,newlength=length1;
  while(i<length1 && j<length2){
      while(i<length1 && temp_neighbors[j]>neighbors[i]){
	i++;
      }
      if(temp_neighbors[j]==neighbors[i]){
	edge_between[i]+=temp_edge_between[j];
	i++;
	j++;
      }
      else {
	edge_between[newlength]=temp_edge_between[j];
	neighbors[newlength]=temp_neighbors[newlength];
	newlength++;
	j++;
      }
  }
  if(j<length2){
	edge_between[newlength]=temp_edge_between[j];
	neighbors[newlength]=temp_neighbors[newlength];
	newlength++;
	j++;
  }
  return newlength;
}

double max(double a, double b){
  if (a>b){
    return a;
  }
  else{
    return b;
  }
}



void update_community(const struct stinger *S, int64_t *seeds, int64_t num_seeds, int64_t nv, int64_t ne, int64_t *membership, int64_t *neighbors,  int64_t *total_community_edges, int64_t *internal_community_edges, int64_t *affected, int64_t num_affected, double *history, int64_t *history_idx){
  
  int64_t i,position, numtocheck=num_affected;
  int64_t *tocheck;
  int64_t *tochecktemp;
  int64_t *didchange;
  int64_t *added;
  int64_t total_comm_edges=*total_community_edges;
  int64_t internal_comm_edges=*internal_community_edges;
  int64_t history_index=*history_idx;
  double history_average=0,history_std_dev=0;
  tocheck=(int64_t*)xmalloc(nv*sizeof(int64_t));
  tochecktemp=(int64_t*)xmalloc(nv*sizeof(int64_t));
  didchange=(int64_t*)xmalloc(nv*sizeof(int64_t));
  added=(int64_t*)xmalloc(nv*sizeof(int64_t));
  for (i=0;i<num_affected;i++){
    tocheck[i]=affected[i];
    
  }
 
 check:
  
  OMP("omp parallel for")
  for (i=0;i<nv;i++){
    added[i]=0;
  }

  history_average=0;
  history_std_dev=0;
  for (i=0;i<HISTORY;i++){
    history_average+=history[i];
  }
  history_average=history_average/(double)HISTORY;
  for (i=0;i<HISTORY;i++){
    history_std_dev=history_std_dev+(history[i]-history_average)*(history[i]-history_average);
  }
  history_std_dev=history_std_dev/(double)HISTORY;
  

    OMP("omp parallel for")
  for (i=0;i<numtocheck;i++){

    int64_t v=tocheck[i];
    double modchange=0;
    int64_t totaledges=stinger_outdegree(S,v);
    int64_t tocomm_edges=0;
    int64_t j;
     STINGER_FORALL_EDGES_OF_VTX_BEGIN(S,v){
       const int64_t dest= STINGER_EDGE_DEST;
       if(membership[dest]==1){
	 tocomm_edges++;
       }
     }STINGER_FORALL_EDGES_OF_VTX_END();
     modchange=(double)tocomm_edges/(double)ne-((double)totaledges/(double)(2*ne))*((double)total_comm_edges/(double)(2*ne));
     if (membership[v]==0 && modchange>=history_average-history_std_dev){
       didchange[i]=1;
       stinger_int64_fetch_add(&internal_comm_edges,tocomm_edges);
       stinger_int64_fetch_add(&total_comm_edges,totaledges);
            int64_t h=stinger_int64_fetch_add(&history_index,1);
            history[h%HISTORY]=modchange;
         }
    
		else if(membership[v]==1 && modchange<=0){
	  
       didchange[i]=-1;
       stinger_int64_fetch_add(&internal_comm_edges,-1*tocomm_edges);
       stinger_int64_fetch_add(&total_comm_edges,-1*totaledges);
     }
     else{
       didchange[i]=0;
     }
     
  }
    OMP("omp parallel for")
  for (i=0;i<numtocheck;i++){
    if(didchange[i]==1){
      membership[tocheck[i]]=1;
    }
    else if(didchange[i]==-1){
      membership[tocheck[i]]=0;
    }
  }
  position=0;
      OMP("omp parallel for")
 for (i=0;i<numtocheck;i++){
   int64_t n=0;
   int64_t v=tocheck[i];
   int64_t *currneigh;
   currneigh=(int64_t*)xmalloc(stinger_outdegree(S,v)*sizeof(int64_t));
    if(didchange[i]==1){
     STINGER_FORALL_EDGES_OF_VTX_BEGIN(S,v){
       const int64_t dest= STINGER_EDGE_DEST;
       if(membership[dest]==0 && stinger_int64_fetch_add(added+dest,1)==0){
	
	 currneigh[n++]=dest;
       }
     }STINGER_FORALL_EDGES_OF_VTX_END();
     int64_t pos=stinger_int64_fetch_add(&position,n);

     memcpy(tochecktemp+pos,currneigh,n*sizeof(int64_t));
    }
    else if(didchange[i]==-1){
     STINGER_FORALL_EDGES_OF_VTX_BEGIN(S,v){
       const int64_t dest= STINGER_EDGE_DEST;
       if(membership[dest]==1 && stinger_int64_fetch_add(added+dest,1)==0){

	 currneigh[n++]=dest;
       }

     }STINGER_FORALL_EDGES_OF_VTX_END();
     int64_t pos=stinger_int64_fetch_add(&position,n);
     

     memcpy(tochecktemp+pos,currneigh,n*sizeof(int64_t));
    }
    free(currneigh);
 }

 memcpy(tocheck,tochecktemp,position*sizeof(int64_t));
 numtocheck=position;

 for (i=0;i<nv;i++){
   didchange[i]=0;
 }
 if (position){
   goto check;
 }

free(tocheck);
free(tochecktemp);
free(added);
free(didchange);
*total_community_edges=total_comm_edges;
*internal_community_edges=internal_comm_edges;
*history_idx=history_index;

double modularity=((double)internal_comm_edges-(((double)total_comm_edges)*(double)total_comm_edges)/((double)(4*ne)))/(double)ne;
 //Modularity+=modularity;
  //Count++;
}


