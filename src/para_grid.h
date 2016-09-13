/********************************************************** 
 * Name   : para_grid.h
 * Author : Kevin Mooney
 * Date   : 22/04/16
 *
 * Struct for a parallel grid of particles
 *
 * REAL has just been typedefed as a double, 
 * kinda useless atm
 **********************************************************/

#pragma once

#include <mpi.h>

#include "point.h"

#define HALO_CAPACITY 10000

enum direction { UP, DOWN, LEFT, RIGHT, UPLEFT, UPRIGHT,DOWNLEFT, DOWNRIGHT };

typedef struct PARA_GRID_S {
  POINT **p;		        /* Array of points */
  REAL width, lwidth;		/* total width and local widht */
  REAL height, lheight;		/* total height and local height */
  REAL origin[2];
  REAL delta;			/* cut off distance */
  int size;			/* n_particles, should change the name.... */

  REAL (*f)(REAL);		/* force */
  REAL (*E)(REAL);		/* energy */

  MPI_Comm comm;		/* Cartesian Communicator */
  int rank;			/* This proccesses rank */
  int nprocs;

  int dir  [8];
  int ndims   ;
  int dims [2];
  int coord[2];			/* hard code in the dimensions, probably bad     */
  int pbc  [2];

  REAL halo_in[8][HALO_CAPACITY], halo_out[8][HALO_CAPACITY];     /* Array to hold incoming and outgoing particles */
  int  halo_size[8][2];
  int  halo_in_size[8][2];
  /* int  swap_size[8]; */
  int  halo_capacity;		/* How many particles the halo will hold, should give error if exceded */

} PARA_GRID ;

PARA_GRID * init_para_grid ( REAL _width, REAL _height, int _size , REAL (*_F)(REAL), REAL (*_E)(REAL),
		   REAL _delta, MPI_Comm _comm);

void free_para_grid ( PARA_GRID *g );

//void write_to_file ( PARA_GRID *g );
