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
   unsigned long int proc0_size;   /* Size of proc 0's subarray */
   unsigned long int prime;        /* Current prime */
   unsigned long int size;         /* Elements in 'marked' */
   unsigned long int low_index; 
   unsigned long int high_index;
   char* temp_marked;
   unsigned long long int temp_low_value;
   unsigned long long int temp_high_value;
   unsigned long int temp_low_index; 
   unsigned long int temp_high_index;
   unsigned long int temp_size;
   unsigned long int temp_index;


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

   /* Stop the timer */
   elapsed_time += MPI_Wtime();


   /* Print the results */

   global_count += 1;// 2 is even but also prime
   if (!id) {
      printf("The total number of prime: %ld, total time: %10.6f, total node %d\n", global_count, elapsed_time, p);
   }
   MPI_Finalize();
   return 0;

}


