/********************************************************** 
 * Name   : para_grid.c
 * Author : Kevin Mooney
 * Date   : 22/04/16
 *
 * Definitions of para_grid.h
 **********************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "para_grid.h"
#include "point.h"

PARA_GRID * init_para_grid ( REAL _width, REAL _height, int _nparticles, REAL (*_F)(REAL), REAL (*_E)(REAL),
			     REAL _delta, MPI_Comm _comm ) {
  PARA_GRID *g = malloc ( sizeof(PARA_GRID) );
  int   i;
  int   n[2];

  g->width   = _width;
  g->height  = _height;
  g->delta   = _delta;
  g->f       = _F;
  g->E       = _E;

  g->comm    = _comm;
  g->delta   = _delta;

  MPI_Comm_rank   ( _comm, &g->rank                     );
  MPI_Comm_size   ( _comm, &g->nprocs                   );
  MPI_Cartdim_get ( _comm, &g->ndims                    );
  MPI_Cart_get    ( _comm, 2, g->dims, g->pbc, g->coord );

  g->size    = _nparticles / g->nprocs; /* Everyone gens the same number of particles */
  g->lwidth  = _width  / g->dims[0];
  g->lheight = _height / g->dims[1];

  g->origin[0] = g->coord[0] * g->lwidth ;
  g->origin[1] = g->coord[1] * g->lheight;

  /* FIND ALL NEIGHBOURS */
  /* down left */
  n[0] = g->coord[0] - 1; n[1] = g->coord[1] - 1;
  MPI_Cart_rank   ( _comm, n, &g->dir[DOWNLEFT] );

  /* down */
  n[0] = g->coord[0];     n[1] = g->coord[1] - 1;
  MPI_Cart_rank   ( _comm, n, &g->dir[DOWN]     );

  /* down right */
  n[0] = g->coord[0] + 1; n[1] = g->coord[1] - 1;
  MPI_Cart_rank   ( _comm, n, &g->dir[DOWNRIGHT]   );

  /* left */
  n[0] = g->coord[0] - 1; n[1] = g->coord[1];
  MPI_Cart_rank   ( _comm, n, &g->dir[LEFT]        );

  /* right */
  n[0] = g->coord[0] + 1; n[1] = g->coord[1];
  MPI_Cart_rank   ( _comm, n, &g->dir[RIGHT]       );

  /* up left */
  n[0] = g->coord[0] - 1; n[1] = g->coord[1] + 1;
  MPI_Cart_rank   ( _comm, n, &g->dir[UPLEFT]      );

  /* up */
  n[0] = g->coord[0];     n[1] = g->coord[1] + 1;
  MPI_Cart_rank   ( _comm, n, &g->dir[UP]          );

  /* up right */
  n[0] = g->coord[0] + 1; n[1] = g->coord[1] + 1;
  MPI_Cart_rank   ( _comm, n, &g->dir[UPRIGHT]     );

  /* printf ( "%d (%f %f): ", g->rank, g->origin[0], g->origin[1] ); */
  /* for ( i=0; i<8; i++ ) */
  /*   printf ( "%d ", g->dir[i] ); */
  /* printf ("\n"); */

  /*************************************************/

  g->p = malloc ( _nparticles * sizeof( POINT * ) );
  for ( i=0; i < g->size; i++ )
    g->p[i] = init_point( g->origin[0] + drand48()*g->lwidth, g->origin[1] + drand48()*g->lheight, 0, 0 );

  return g;
}

void free_para_grid ( PARA_GRID *g ) {
  int i;
  for ( i=0; i < g->size; i++ )
    free_point ( g->p[i] );
  free(g->p);
  free(g);
}

/* void write_to_file ( PARA_GRID *g ) { */
/*   FILE *fp; */
/*   char filename[100]; */

/*   sprintf ( filename, "data_par/parafile_%d.dat", g->rank ); */

/*   fp = fopen ( filename, "w" ); */

/*   for ( int i=0; i < g->size; i++ ) */
/*     fprintf ( fp, "%10lf %10lf\n", g->p[i]->rx, g->p[i]->ry ); */
/*   fprintf ( fp, "\n" ); */

/* } */
