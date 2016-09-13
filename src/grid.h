/********************************************************** 
 * Name   : grid.h
 * Author : Kevin Mooney
 * Date   : 14/04/16
 *
 * Struct of a grid of particles
 **********************************************************/

#pragma once

#include "point.h"

typedef struct GRID_S {
  POINT **p;
  REAL width;
  REAL height;
  int size;
  REAL (*f)(REAL);
  REAL (*E)(REAL);
} GRID ;

GRID * init_grid ( REAL _width, REAL _height, int _size , REAL (*_F)(REAL), REAL (*_E)(REAL));

void free_grid ( GRID *g );
