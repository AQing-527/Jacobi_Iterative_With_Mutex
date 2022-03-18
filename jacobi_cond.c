#define _GNU_SOURCE

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>              /* For timing */
#include <sys/time.h>            /* For timing */
#include <sys/resource.h>
#include <string.h>
#include <pthread.h>
#include <stdbool.h>
#include <semaphore.h>


/****************Global****************************/

#define MAX(a,b) ((a)>(b)?(a):(b))
#define EPSILON 0.001            /* Termination condition */

char *filename;                  /* File name of output file */

/* Grid size */
int M = 200;                     /* Number of rows */
int N = 200;                     /* Number of cols */
long max_its = 1000000;          /* Maximum iterations, a safe bound to avoid infinite loop */
double final_diff;               /* Temperature difference between iterations at the end */

/* Thread count */
int thr_count = 2;

/* shared variables between threads */
/*************************************************************/
double** u;                   /* Previous temperatures */
double** w;                   /* New temperatures */


pthread_mutex_t lock; 
pthread_cond_t finished, startNext; /* Condition variables */
int count = 0; /* Number of threads that has finished the calculation */
int terminate = 0; /* Boolean to check whether terminate the process */
double its_diff = 0.0; /* Max temperature difference in each iteration */
long MICRO = 1000000;   /* Constant that is used to convert ps to sec */
typedef struct worker_thread{
   int id;
   struct rusage usage;
}worker_thread; /* A struct to store the information of each thread */




/**************************************************************/

int main (int argc, char *argv[])
{
   int      its;                 /* Iterations to converge */
   double   elapsed;             /* Execution time */
   struct timeval stime, etime;  /* Start and end times */
   struct rusage usage;

   void allocate_2d_array (int, int, double ***);
   void initialize_array (double ***);
   void print_solution (char *, double **);
   int  find_steady_state (void);

   /* For convenience of other problem size testing */
   if ((argc == 1) || (argc == 4)) {
      if (argc == 4) {
         M = atoi(argv[1]);
         N = atoi(argv[2]);
         thr_count = atoi(argv[3]);
      } // Otherwise use default grid and thread size
   } else {
     printf("Usage: %s [ <rows> <cols> <threads ]>\n", argv[0]);
     exit(-1);
   }
   printf("Problem size: M=%d, N=%d\nThread count: T=%d\n", M, N, thr_count);
   /* Create the output file */
   filename = argv[0];
   sprintf(filename, "%s.dat", filename);

   allocate_2d_array (M, N, &u);
   allocate_2d_array (M, N, &w);
   initialize_array (&u);
   initialize_array (&w);

   gettimeofday (&stime, NULL);
   its = find_steady_state();
   gettimeofday (&etime, NULL);

   elapsed = ((etime.tv_sec*1000000+etime.tv_usec)-(stime.tv_sec*1000000+stime.tv_usec))/1000000.0;

   printf("Converged after %d iterations with error: %8.6f.\n", its, final_diff);
   printf("Elapsed time = %8.4f sec.\n", elapsed);

   getrusage(RUSAGE_SELF, &usage);
   printf("Program completed - user: %.4f s, system: %.4f s\n",
      (usage.ru_utime.tv_sec + usage.ru_utime.tv_usec/1000000.0),
    (usage.ru_stime.tv_sec + usage.ru_stime.tv_usec/1000000.0));
   printf("no. of context switches: vol %ld, invol %ld\n\n",
  		  usage.ru_nvcsw, usage.ru_nivcsw);

   print_solution (filename, w);
}

/* Allocate two-dimensional array. */
void allocate_2d_array (int r, int c, double ***a)
{
   double *storage;
   int     i;
   storage = (double *) malloc (r * c * sizeof(double));
   *a = (double **) malloc (r * sizeof(double *));
   for (i = 0; i < r; i++)
      (*a)[i] = &storage[i * c];
}

/* Set initial and boundary conditions */
void initialize_array (double ***u)
{
   int i, j;

   /* Set initial values and boundary conditions */
   for (i = 0; i < M; i++) {
      for (j = 0; j < N; j++)
         (*u)[i][j] = 25.0;      /* Room temperature */
      (*u)[i][0] = 0.0;
      (*u)[i][N-1] = 0.0;
   }

   for (j = 0; j < N; j++) {
      (*u)[0][j] = 0.0;
      (*u)[M-1][j] = 1000.0;     /* Heat source */
   }
}

/* Print solution to standard output or a file */
void print_solution (char *filename, double **u)
{
   int i, j;
   char sep;
   FILE *outfile;

   if (!filename) { /* if no filename specified, print on screen */
      sep = '\t';   /* tab added for easier view */
      outfile = stdout;
   } else {
      sep = '\n';   /* for gnuplot format */
      outfile = fopen(filename,"w");
      if (outfile == NULL) {
         printf("Can't open output file.");
         exit(-1);
      }
   }

   /* Print the solution array */
   for (i = 0; i < M; i++) {
      for (j = 0; j < N; j++)
         fprintf (outfile, "%6.2f%c", u[i][j], sep);
      fprintf(outfile, "\n"); /* Empty line for gnuplot */
   }
   if (outfile != stdout)
      fclose(outfile);

}

/* Entry function of the worker threads */
void *thr_func(void *arg) {

   worker_thread * this_arg = (worker_thread *)arg;

   int worker_id = this_arg->id;
   int rows = M / thr_count;
   int start; /* The starting row of this worker thread */
   int end; /* The end row of this worker thread */
   int i,j;
   double inner_diff; /* Max temperature difference of this worker thread in this iteration */

   // Initialize staring row and end row
   if (worker_id == 0){
      start = 1;
      end = rows;
   }

   else if (worker_id == (thr_count-1)){
      start = rows*worker_id;
      end = N-1;
   }

   else{
      start = rows*worker_id;
      end = rows*(worker_id+1); 
   }
   
   while (1) {
      inner_diff = 0.0;
      // Do calculations
      for (i = start; i < end; i++) {
         for (j = 1; j < N-1; j++) {
            w[i][j] = 0.25 * (u[i-1][j] + u[i+1][j] + u[i][j-1] + u[i][j+1]);
            if (fabs(w[i][j] - u[i][j]) > inner_diff)
               inner_diff = fabs(w[i][j] - u[i][j]);
         }
      }

      // Get the lock to change the variables safyly
      pthread_mutex_lock(&lock);
      count = count + 1;
      if (inner_diff > its_diff){
         its_diff = inner_diff;
      }
      if (count == thr_count){
         pthread_cond_signal(&finished); /* The last one thread will send the signal of finish */
      }
      pthread_cond_wait(&startNext, &lock); /* Wait for the master thread to start next iteration */

      // Check whether the termination condition is satisfied
      if (terminate == 1){
         pthread_mutex_unlock(&lock);
         break;
      }
      pthread_mutex_unlock(&lock);
   }
   // Get the execution statistics of this working thread
   getrusage(RUSAGE_THREAD, &this_arg->usage);
   pthread_exit(0);
}


int find_steady_state (void)
{

   int its;             /* Iteration count */
   int i, j;
   struct rusage ru;

   pthread_t thread_id[thr_count];

   worker_thread thread_args[thr_count];

   for (int k=0; k<thr_count; k++){
      thread_args[k].id = k;
   }

   double user_time, sys_time;

   // Initialize mutex for syncronization
   pthread_mutex_init(&lock, NULL);
   // Initialize the condition varaibles
   pthread_cond_init(&finished, NULL);
   pthread_cond_init(&startNext, NULL);
   
   // Create threads
   for(i=0; i < thr_count; i++)
   {
      thread_args[i].id = i;
      pthread_create( &thread_id[i], NULL, thr_func, &thread_args[i]);
   }

   pthread_mutex_lock(&lock);
   for (its = 1; its <= max_its; its++) {
      while (count != thr_count){
         pthread_cond_wait(&finished, &lock); /* Wait for all worker threads finishing their work */
      }
      
      /* After all worker threads finishes their work */
      /* Swap matrix u, w */
      double **temp = u;
      u = w;
      w = temp;
      
      /* Reset the count */
      count = 0;

      /* Terminate if temperatures have converged or the max_its is reached */
      if (its_diff <= EPSILON || its == max_its){
         terminate = 1;
         final_diff = its_diff;
         pthread_cond_broadcast(&startNext);
         break;
      }
         
      its_diff=0.0;
      pthread_cond_broadcast(&startNext);
   }
   pthread_mutex_unlock(&lock);

   // Destroy the lock
   pthread_mutex_destroy(&lock);   
   pthread_cond_destroy(&finished);
   pthread_cond_destroy(&startNext);

   // Join threads
   for(j=0; j < thr_count; j++)
   {
      pthread_join(thread_id[j], NULL);
   }

   for (int k=0; k<thr_count; k++){
      user_time = thread_args[k].usage.ru_utime.tv_sec + (double)thread_args[k].usage.ru_utime.tv_usec / MICRO;
      sys_time = thread_args[k].usage.ru_stime.tv_sec + (double)thread_args[k].usage.ru_stime.tv_usec / MICRO;
      printf("Thread %d has completed - user: %.4f s, system: %.4f s\n", thread_args[k].id, user_time, sys_time);
   }

   getrusage(RUSAGE_SELF, &ru);
   user_time = ru.ru_utime.tv_sec + (double)ru.ru_utime.tv_usec / MICRO;
   sys_time = ru.ru_stime.tv_sec + (double)ru.ru_stime.tv_usec / MICRO;
   printf("find_steady_state - user: %.4f s, system: %.4f s\n", user_time, sys_time);

   return its;		//return the number of iterations


}
