/********************************************************** 
 * Name   : n-body.c
 * Author : Kevin Mooney
 * Date   : 14/04/16
 *
 * N-body simulation
 **********************************************************/


#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>

#include "integrator.h"

#define NDIMS   2
#define REORDER 0		/* Allow MPI to reorder tasks, off for now */

/* Natural Lennard-Jones Gradient */
REAL force( REAL x ) {
  return pow(x,-13) - pow(x,-7);
} 

/* Natural Lennard-Jones Potential */
REAL energy( REAL x ) {
  return ( pow(x,-12) - 2 * pow(x,-6) ) / 12.0;
}

/* Print a usage message */
void print_usage(char * progname) {
  fprintf(stderr,"Usage: %s [ -hnstw ]\n",progname);
  fprintf(stderr,"Performs gauss elimination\n");
  fprintf(stderr,"          -h  Set the height of cell \n");
  fprintf(stderr,"          -n  Set number of bodies  \n");
  fprintf(stderr,"          -s  Set number of steps   \n");
  fprintf(stderr,"          -t  Set step size         \n");
  fprintf(stderr,"          -w  Set the width of cell\n");
}

void print_header(GRID* g, double step_size, int n_steps) {
  printf("WIDTH      %lf\n",g->width  );
  printf("HEIGHT     %lf\n",g->height );
  printf("STEPS      %d \n",n_steps   );
  printf("TIME STEP  %lf\n",step_size );
  printf("NBODIES    %d \n",g->size   );
  printf("\n\n");
}

int main(int argc, char *argv[]) {
  
  MPI_Comm cart_comm;
  GRID      *  g;
  PARA_GRID * pg;
  struct timeval t0, t1;
  long elapsed;
  double step_size, width, height, delta;
  int option;
  int n_points, n_steps;
  int n_procs, rank;

  /* These are specific to a 2D problem, 
     will need to be rewritten if for some
     reason we want a higher dimension    */
  int dims[] = { 0, 0 };	/* dim array to be filled by MPI_Dims_create */
  int pbc [] = { 1, 1 };	/* Set periodic boundary conditions */

  /* defaults */
  n_points       = 2;
  step_size      = 0.00001;
  n_steps        = 100000;
  width = height = 2.0;
  delta          = 0.1; /* Need to check if this is a good value */

  while((option = getopt(argc,argv,"h:n:s:t:w:")) != -1) {
    switch(option) {
    case 'h':
      height    = atof(optarg);
      break;
    case 'n':
      n_points  = atoi(optarg);
      break;
    case 's':
      n_steps   = atoi(optarg);
      break;
    case 't':
      step_size = atof(optarg);
      break;
    case 'w':
      width     = atof(optarg);
      break;
    default :
      print_usage(argv[0]);
      return EXIT_FAILURE;
    }
  }
  
  /* Setup MPI */
  MPI_Init( &argc, &argv );
  MPI_Comm_size( MPI_COMM_WORLD,&n_procs );
  MPI_Comm_rank( MPI_COMM_WORLD,&rank );

  /* Create Cartesian Topology */
  MPI_Dims_create ( n_procs, NDIMS, dims );
  MPI_Cart_create ( MPI_COMM_WORLD, NDIMS, dims, pbc, REORDER, &cart_comm );

  int coord[2];
  MPI_Cart_coords(cart_comm,rank,2,coord);
  /* printf("P:%d My coordinates are %d %d\n",rank,coord[0],coord[1]); */

  /* Seed RNG with some variable plus rank*/
  srand48(3574235876*rank);

  /* INIT SIMULATION */
  g     = init_grid      ( width, height, n_points, force, energy );
  pg    = init_para_grid ( width, height, n_points, force, energy, delta, cart_comm );

  /* write_to_file(pg ); */

  /* START SIMULATION */
  /* print_header ( g, step_size, n_steps ) ; */
  /* leap_frog      ( g , step_size, n_steps ) ; */
  if ( rank == 0 )  gettimeofday(&t0, 0);

  para_leap_frog ( pg, step_size, n_steps ) ;

  if ( rank == 0 ) {
    gettimeofday(&t1, 0);
    elapsed = (t1.tv_sec-t0.tv_sec)*1000000 + t1.tv_usec-t0.tv_usec;  
    printf ("%d %d %ld\n", n_points, n_procs, elapsed );
  }

  /* FREE MEMORY */

  free_grid ( g );
  free_para_grid ( pg );
  
  MPI_Finalize();

  return EXIT_SUCCESS;
}
