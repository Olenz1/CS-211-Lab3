/*
*   Sieve of Eratosthenes
*
*   Programmed by Michael J. Quinn
*
*   Last modification: 7 September 2001
*/

#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define MIN(a, b)  ((a)<(b)?(a):(b))

int main(int argc, char *argv[]) {
   unsigned long int count;        /* Local prime count */
   double elapsed_time; /* Parallel execution time */
   unsigned long int first;        /* Index of first multiple */
   unsigned long int global_count = 0; /* Global prime count */
   unsigned long long int high_value;   /* Highest value on this proc */
   unsigned long int i;
   int id;           /* Process ID number */
   unsigned long int index;        /* Index of current prime */
   unsigned long long int low_value;    /* Lowest value on this proc */
   char *marked;       /* Portion of 2,...,'n' */
   unsigned long long int n;            /* Sieving from 2, ..., 'n' */
   int p;            /* Number of processes */
   unsigned long int proc0_size;   /* Size of proc 0's temparray */
   unsigned long int prime;        /* Current prime */
   unsigned long int size;         /* Elements in 'marked' */
   unsigned long int low_index; 
   unsigned long int high_index;
   char* temp_marked;
   unsigned long long int temp_low_value;
   unsigned long long int temp_high_value;
   unsigned long int temp_low_index; 
   unsigned long int temp_high_index;

   MPI_Init(&argc, &argv);

   /* Start the timer */

   MPI_Comm_rank(MPI_COMM_WORLD, &id);
   MPI_Comm_size(MPI_COMM_WORLD, &p);
   MPI_Barrier(MPI_COMM_WORLD);
   elapsed_time = -MPI_Wtime();

   if (argc != 2) {
      if (!id) printf("Command line: %s <m>\n", argv[0]);
      MPI_Finalize();
      exit(1);
   }

   n = atoll(argv[1]);

   /* Figure out this process's share of the array, as
      well as the integers represented by the first and
      last array elements */
   proc0_size = (n - 1) / p;
   if ((2 + proc0_size) < (int) sqrt((double) n)) {
      if (!id) printf("Too many processes\n");
      MPI_Finalize();
      exit(1);
   }
   
   int temp_n = (int)sqrt((double)n);
   int temp_N = (temp_n - 1) >> 1;

   temp_low_index = 0 * (temp_N / p) + MIN(0, temp_N % p); 
   temp_high_index = 1 * (temp_N / p) + MIN(1, temp_N % p) - 1; 
   temp_low_value = temp_low_index * 2 + 3; 
   temp_high_value = (temp_high_index + 1) * 2 + 1;

   temp_marked = (char*)malloc(temp_n);
   if(temp_marked == NULL) {
      printf("Cannot allocate enough memory \n");
      MPI_Finalize();
      exit(1);
   }

   for (i = 0; i < temp_n; i++) temp_marked[i] = 0;
   index = 0;
   prime = 3;
   do {
      first = (prime * prime - temp_low_value) >> 1;

      for (i = first; i < temp_n; i += prime) temp_marked[i] = 1;
      
      while(temp_marked[++index]);
      prime = 2 * index + 3;
      //if (p > 1) MPI_Bcast(&prime, 1, MPI_INT, 0, MPI_COMM_WORLD);
   } while (prime * prime <= temp_n);

   int N = (n - 1) >> 1;
   low_index = id * (N / p) + MIN(id, N % p);
   high_index = (id + 1) * (N / p) + MIN(id + 1, N % p) - 1;
   low_value = low_index * 2 + 3;
   high_value = (high_index + 1) * 2 + 1;
   size = (high_value - low_value) / 2 + 1;

   int LEVEL1_CACHE_size = 16384;      // default 16384
   int LEVEL2_CACHE_size = 2097152;     // default 2097152
   int LEVEL3_CACHE_size = 11582912;    // default 12582912

   int CACHE_size = atoi(argv[2]);

   int LEVEL1_CACHE_int = LEVEL1_CACHE_size / 4;
   int LEVEL2_CACHE_int = LEVEL2_CACHE_size / 4;
   int LEVEL3_CACHE_int = LEVEL3_CACHE_size / 4;
   int CACHE_int = CACHE_size / 4;

   int Block_size = CACHE_int / p ;
   int Block_num = size / Block_size;
   int Block_remain = size % Block_size;

   int Block_id = 0;
   int Block_N = Block_size - 1;
   int Block_low_index = Block_id * Block_N + MIN(Block_id, Block_remain) + low_index;
   int Block_high_index = (Block_id + 1) * Block_N + MIN(Block_id + 1, Block_remain) - 1 + low_index;
   int Block_low_value = Block_low_index * 2 + 3;
   int Block_high_value = (Block_high_index + 1) * 2 + 1;
   int Block_count;

   marked = (char *) malloc(Block_size);
   if (marked == NULL) {
      printf("Cannot allocate enough memory \n");
      MPI_Finalize();
      exit(1);
   }

   count = 0; 

   while(Block_id <= Block_num) {
      index = 0;
      prime = 3;
      Block_count = 0;
      
      for (i = 0; i < Block_size; i++) marked[i] = 0;

      do {
            if (prime * prime > Block_low_value) {
                  first = (prime * prime - Block_low_value) / 2;
            } else {
                  if (!(Block_low_value % prime)) first = 0;
                  else if (Block_low_value % prime % 2 == 0) first = prime - ((Block_low_value % prime) / 2);
                  else first = (prime - (Block_low_value % prime)) / 2;
            }

            for (i = first; i < Block_size; i += prime) {
                  marked[i] = 1;
            }

            while (temp_marked[++index]);

            prime = index * 2 + 3; 

         } while (prime * prime <= Block_high_value);

         for (i = 0; i < Block_size; i++) {
            if (marked[i] == 0) {
                  Block_count++;
            }
         }

         count += Block_count;

         Block_id++;
         Block_low_index = Block_id * Block_N + MIN(Block_id, Block_remain) + low_index;
         Block_high_index = (Block_id + 1) * Block_N + MIN(Block_id + 1, Block_remain) - 1 + low_index;
         Block_low_value = Block_low_index * 2 + 3;
         if (Block_id == Block_num) {
            Block_high_value = high_value;
            Block_high_index = high_index;
            Block_size = (Block_high_value - Block_low_value) / 2 + 1;
         } else {
            Block_high_value = (Block_high_index + 1) * 2 + 1;
         }

   }

   MPI_Reduce(&count, &global_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
   /* Stop the timer */
   elapsed_time += MPI_Wtime();


   /* Print the results */

   //global_count += 1;// 2 is even but also prime
   if (!id) {
      printf("The total number of prime: %ld, total time: %10.6f, total node %d\n", global_count, elapsed_time, p);
   }
   MPI_Finalize();
   return 0;

}


