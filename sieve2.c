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

   /* Figure out this process's share of the array, as
      well as the integers represented by the first and
      last array elements */


   int temp_n = (int)sqrt((double)n);
   int temp_N = (temp_n - 1) >> 1;

   temp_low_index = 0 * (temp_N / p) + MIN(0, temp_N % p); 
   temp_high_index = 1 * (temp_N / p) + MIN(1, temp_N % p) - 1; 
   temp_low_value = temp_low_index * 2 + 3; 
   temp_high_value = (temp_high_index + 1) * 2 + 1;

   /* Bail out if all the primes used for sieving are
      not all held by process 0 */

   proc0_size = (temp_n - 1) / p;

   if ((2 + proc0_size) < (int) sqrt((double) temp_n)) {
      if (!id) printf("Too many processes\n");
      MPI_Finalize();
      exit(1);
   }

   /* Allocate this process's share of the array. */

   temp_marked = (char *)malloc(temp_n);
   if (temp_marked == NULL) {
      printf("Cannot allocate enough memory\n");
      MPI_Finalize();
      exit(1);
   }

   for (i = 0; i < temp_n; i++) temp_marked[i] = 0;
   index = 0;
   prime = 3;
   do {
      // if (prime * prime > temp_low_value)
      //       first = (prime * prime - temp_low_value) >> 1;
      // else {
      //    if (!(low_value % prime)) first = 0;
      //    else if(low_value % prime % 2 == 0)
      //       first = prime - ((low_value % prime) >> 1);       // 此处在求局部first（数组中第一个可以被prime整除的数）的时候非常巧妙
      //    else
      //       first = (prime - (low_value % prime)) >> 1;
      // }
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

   marked = (char *) malloc(size);
   if (marked == NULL) {
      printf("Cannot allocate enough memory \n");
      MPI_Finalize();
      exit(1);
   }

   for (i = 0; i < size; i++) marked[i] = 0;

   index = 0;
   prime = 3;
   do {
      if (prime * prime > low_value)
         first = (prime * prime - low_value) >> 1;
      else {
         if (!(low_value % prime)) first = 0;
         else if (low_value % prime % 2 == 0) first = prime - ((low_value % prime) >> 1);
         else first = prime - ((low_value % prime) >> 1);
      }

      for (i = 0; i < size; i += prime) marked[i] = 1;

      while (temp_marked[++index]);
      prime = 2 * index + 3;
   } while (prime * prime <= n);

   count = 0; 
   for (i = 0; i < size; i++)
      if (!marked[i]) count++;
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


