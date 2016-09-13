/********************************************************** 
 * Name   : grid.c
 * Author : Kevin Mooney
 * Date   : 14/04/16
 *
 * Definitions of grid.h
 **********************************************************/

#include <stdlib.h>

#include "grid.h"
#include "point.h"

GRID * init_grid ( REAL _width, REAL _height, int _nparticles, REAL (*_F)(REAL), REAL (*_E)(REAL)) {
  GRID *g = malloc ( sizeof(GRID) );
  int   i;

  g->width  = _width;
  g->height = _height;
  g->size   = _nparticles;
  g->f      = _F;
  g->E      = _E;

  g->p = malloc ( _nparticles * sizeof( POINT * ) );
  for ( i=0; i < _nparticles; i++ )
    g->p[i] = init_point( drand48()*_width, drand48()*_height, 0, 0 );

  return g;
}

void free_grid ( GRID *g ) {
  int i;
  for ( i=0; i < g->size; i++ )
    free_point ( g->p[i] );
  free(g->p);
  free(g);
}
