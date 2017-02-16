#ifndef MATRIX_TREATMENT
#define MATRIX_TREATMENT

#define GET_MTRX_POS(ROW,COL,WIDTH) (WIDTH) * (ROW) + (COL)
#define BUFF_SIZE 999999
#define STD_TAG 1
#define INF  100000000
#define MPI_REDUCE_BLOCKSIZE 1000

/*******************************************************************************
* prints_matrix: Prints the passed matrix
*******************************************************************************/
void prints_matrix(int *d_graph, int N){
  int j,i;
  for (i = 0; i < N; i++){
    for(j=0; j < N-1; j++)
      printf("%d ",  d_graph[GET_MTRX_POS(i,j,N)]);
    printf("%d", d_graph[GET_MTRX_POS(i,j,N)]);
    printf("\n");
  }
}/*  prints_matrix */

/*******************************************************************************
* check_fox_conditions: Checks if it is possible to divide the matrix D1
                        in submatrices ( one submatrix to each processor  
                        running)
*******************************************************************************/
int check_fox_conditions(int p, int n){
  int temp;
  temp = sqrt(p);
  if(temp*temp == p){
    if(n%temp == 0){
      return temp;
    }
  }
  printf("It's not possible to apply the Fox algorithm\n");
  return 1;
}/* check_fox_conditions */


/*******************************************************************************
* divide_matrix: Divides the matrix accordingly to the number of 
                 processores running
*******************************************************************************/
int *divide_matrix(int *d_graph,int N , int q){
  int* temp_graph;
  int i,j,r,c, dest_count;
  int lngth = N / q;
  //aloca o espaço necessário para a matriz temporária
  temp_graph = (int*)malloc((lngth*lngth) * sizeof(int));

  dest_count = 0;
  for(r=0; r < q; r ++){
      for(c=0; c < q; c ++){
	// copia a parte da matriz correspondente ao processador (r,c)
	for (i = 0; i < lngth; i++) {
	  for (j = 0; j < lngth; j++) {
	    temp_graph[GET_MTRX_POS(i,j,lngth)] = d_graph[GET_MTRX_POS( i+(lngth*r),j+(lngth*c),N)] ;
	  }
	}
	chnkd_MPI_Send(temp_graph, lngth*lngth, MPI_INT,dest_count);
	dest_count++;
      }
  }
}/* divide_matrix */


/*******************************************************************************
* floyd_warshall: Calculates the minimum
*******************************************************************************/
void floyd_warshall(int *mtrx_a,int *mtrx_b,int *mtrx_res, int lngth){
  int i, j, k;
  for(i=0; i<lngth; i++){
    for(j=0; j<lngth; j++){
      for(k=0; k<lngth; k++){
	  if (mtrx_a[GET_MTRX_POS(i,k,lngth)] + mtrx_b[GET_MTRX_POS(k,j,lngth)] < mtrx_res[GET_MTRX_POS(i,j,lngth)])
	    mtrx_res[GET_MTRX_POS(i,j,lngth)] = mtrx_a[GET_MTRX_POS(i,k,lngth)] + mtrx_b[GET_MTRX_POS(k,j,lngth)];
	
      }
    }
  }
}/* floyd_warshall */


/*******************************************************************************
* radixsort: Orders the array in O(n) time using radixsort
*******************************************************************************/
void radixsort(int *vetor, int tamanho) {
    int i;
    int b[tamanho];
    int maior = vetor[0];
    int exp = 1;
 
    for (i = 0; i < tamanho; i++) 
	if (vetor[i] > maior)
    	    maior = vetor[i];

    while (maior/exp > 0) {
	int bucket[10] = { 0 };
    	for (i = 0; i < tamanho; i++)
    	    bucket[(vetor[i] / exp) % 10]++; 
    	for (i = 1; i < 10; i++)
    	    bucket[i] += bucket[i - 1];
    	for (i = tamanho - 1; i >= 0; i--)
    	    b[--bucket[(vetor[i] / exp) % 10]] = vetor[i];
    	for (i = 0; i < tamanho; i++)
    	    vetor[i] = b[i];
    	exp *= 10;
    }
} /* radixsort */

/*******************************************************************************
* chnkd_MPI_Send: Send array to MPI node by blocks to avoid buffer limit
*******************************************************************************/
int chnkd_MPI_Send(void *buf, long count, MPI_Datatype type,int dest){
  int size;
  long offset=0;
  int  chunck_size=MPI_REDUCE_BLOCKSIZE;
  
  if (!buf || count <= 0) return(MPI_ERR_COUNT);
  MPI_Type_size(type, &size);

  while (offset < count) {
    if (offset+chunck_size > count-1) 
      chunck_size=count-offset; 
    else 
      chunck_size=MPI_REDUCE_BLOCKSIZE;
    MPI_Send((void*)(buf+offset*size), chunck_size, type, dest, STD_TAG, MPI_COMM_WORLD);
    offset += chunck_size;
  }

  return MPI_SUCCESS;
} /* chnkd_MPI_Send */

/*******************************************************************************
* chnkd_MPI_Recv: Receives arrays from MPI nodes by blocks to avoid buffer limit
*                 the buffer must have been allocated previously.
*******************************************************************************/
int chnkd_MPI_Recv(void *buf,long count, MPI_Datatype type,int source)
{
  int size;
  long offset=0;
  int  chunck_size=MPI_REDUCE_BLOCKSIZE;
  
  if (!buf || count <= 0) return(MPI_ERR_COUNT);
  MPI_Type_size(type, &size);

  while (offset < count) {
    if (offset+chunck_size > count-1) 
      chunck_size=count-offset; 
    else 
      chunck_size=MPI_REDUCE_BLOCKSIZE;
    MPI_Recv((void*)(buf+offset*size), chunck_size, type, source, STD_TAG,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    offset += chunck_size;
  }

  return MPI_SUCCESS;
} /* mc_MPI_Recv */
#endif 
