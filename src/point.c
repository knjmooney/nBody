/********************************************************** 
 * Name   : point.c
 * Author : Kevin Mooney
 * Date   : 14/04/16
 *
 * Definitions of point.h
 **********************************************************/

#include <stdlib.h>

#include "point.h"

POINT * init_point ( REAL _rx, REAL _ry, REAL _vx, REAL _vy ) {
  POINT * p = malloc ( sizeof ( POINT ) );

  p->rx = _rx; p->ry = _ry;
  p->vx = _vx; p->vy = _vy;

  return p;
}

void copy_point ( POINT * p1, POINT *p2 ) {
  p1 -> rx = p2 -> rx;
  p1 -> ry = p2 -> ry;
  p1 -> vx = p2 -> vx;
  p1 -> vy = p2 -> vy;
  
}

void free_point ( POINT * p ) {
  free ( p ) ;
}
