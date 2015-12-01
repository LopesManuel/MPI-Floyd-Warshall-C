#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mpi.h"
#include "Matrix_Treatment.h"

#define ROOT 0 

int main(int argc, char **argv) {

 
  int rank, M, j,i, *d_graph;
  int *local_matrix, *row_matrix, *col_matrix, *res_matrix, *rowIds, *colIds;
  int P, N, q, p_row, p_col;
  double start, finish;
  MPI_Status status;
 
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &P);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  //INPUT HANDLED BY THE ROOT PROCESSOR
  if (rank == ROOT){
    scanf("%d", &N);  
    q = check_fox_conditions(P,N);

    //Check's if the fox's conditions are met
    if(q == 0){
      MPI_Abort(MPI_COMM_WORLD, 0);
      return 1; //error
    }  

    d_graph = (int*)malloc((N*N) * sizeof(int));

    for(i=0; i < N; i++){
      for(j=0; j < N; j++){
	scanf("%d", &d_graph[GET_MTRX_POS(i,j,N)]);
	if (d_graph[GET_MTRX_POS(i,j,N)] == 0 && i != j) {
	  d_graph[GET_MTRX_POS(i,j,N)] = INF;
	}
      }
    }



    MPI_Bcast(&q, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if(q > 1)
      divide_matrix( d_graph, N, q); 
      
  }
  else{
    MPI_Bcast(&q, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
  }
  //---------------COMMON------------------
   
  int lngth = N / q;


  local_matrix = (int*)malloc((lngth*lngth) * sizeof(int));
  row_matrix   = (int*)malloc((lngth*lngth) * sizeof(int));
  col_matrix   = (int*)malloc((lngth*lngth) * sizeof(int));
  res_matrix   = (int*)malloc((lngth*lngth) * sizeof(int));
  
  if(q>1)
    chnkd_MPI_Recv(local_matrix, lngth*lngth, MPI_INT, 0);
  else
    local_matrix = d_graph;
    
  p_row = ( rank / q );
  p_col = ( rank % q );
    
  //CREATE COMMUNICATORS 
  MPI_Group MPI_GROUP_WORLD;
  MPI_Comm_group(MPI_COMM_WORLD, &MPI_GROUP_WORLD);
  MPI_Group row_group, col_group;
  MPI_Comm row_comm, col_comm, grid_comm;
  int tmp_row, tmp_col, proc;
  int row_process_ranks[q], col_process_ranks[q];
    
  for(proc = 0; proc < q; proc++){   
    row_process_ranks[proc] = (p_row * q) + proc;
    col_process_ranks[proc] = ((p_col + proc*q) %(q*q));
  }    
  radixsort(col_process_ranks, q);
  radixsort(row_process_ranks, q);

  MPI_Group_incl(MPI_GROUP_WORLD, q, row_process_ranks, &row_group);  
  MPI_Group_incl(MPI_GROUP_WORLD, q, col_process_ranks, &col_group);  
     
  MPI_Comm_create(MPI_COMM_WORLD, row_group, &row_comm);  
  MPI_Comm_create(MPI_COMM_WORLD, col_group, &col_comm);  

  if ((rank / q) == (rank % q)) {
      memcpy(row_matrix, local_matrix, (lngth*lngth) * sizeof(int));
  }
  int ln,d,flag;
  int step, rotation_src, rotation_dest, src;
  int count = 0;
  memcpy(res_matrix, local_matrix, (lngth*lngth) * sizeof(int));
  rotation_src = (p_row + 1) % q;
  rotation_dest = ((p_row - 1) + q) % q;
  ln = (lngth*q) << 1;
  start = MPI_Wtime();  

  for (d = 2; d < ln; d = d << 1) {
    memcpy(col_matrix, local_matrix, (lngth*lngth) * sizeof(int));
    for ( step = 0;  step < q;  step++) {
      src = (p_row +  step) % q;
      count++;
      if (src == p_col) {
	MPI_Bcast(local_matrix, lngth*lngth, MPI_INT, src, row_comm);
	floyd_warshall( local_matrix, col_matrix, res_matrix, lngth);
      } else {
	MPI_Bcast(row_matrix, lngth*lngth, MPI_INT, src, row_comm);
	floyd_warshall( row_matrix, col_matrix, res_matrix, lngth);
      }  
      if( step < q-1) 
        MPI_Sendrecv_replace(col_matrix, lngth*lngth, MPI_INT, rotation_dest, STD_TAG,rotation_src, STD_TAG, col_comm, MPI_STATUS_IGNORE);
  	
    }
    memcpy(local_matrix, res_matrix, (lngth*lngth) * sizeof(int));
  }
  
  
  int *sol;
  sol = malloc(N*N*sizeof(int));  
  
  MPI_Gather(res_matrix, lngth*lngth, MPI_INT, sol,  lngth*lngth, MPI_INT, 0, MPI_COMM_WORLD);
  
  if (rank == 0) {
    finish = MPI_Wtime();
    printf("Tempo de execução %f\n",finish - start);
  }
 
  if (rank == 0) {
    int row, col, pos_x, pos_y, pos, tmp_y, tmp_x;

    for (i = 0; i < P; i++) {
      pos_x = i / q;
      pos_y = i % q;
      pos = i * lngth*lngth;

      for (row = 0; row < lngth; row++) {
	for (col = 0; col < lngth; col++) {
          tmp_x = GET_MTRX_POS(pos_x,row,lngth);
          tmp_y = GET_MTRX_POS(pos_y,col,lngth);
          
	  if (sol[GET_MTRX_POS(row,col,lngth) + pos] == INF)
	    d_graph[GET_MTRX_POS(tmp_x,tmp_y,N)] = 0;
	  else
	    d_graph[GET_MTRX_POS(tmp_x,tmp_y,N)] = sol[GET_MTRX_POS(row,col,lngth) + pos];
	}
      }
    }
    prints_matrix(d_graph,N);
  }
  
  MPI_Finalize();
  return 0;
}
