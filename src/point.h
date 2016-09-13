
/********************************************************** 
 * Name   : point.h
 * Author : Kevin Mooney
 * Date   : 14/04/16
 *
 * Point struct to describe the current state of a point
 **********************************************************/

#pragma once

#include "types.h"

typedef struct POINT_S {
  REAL rx, ry;
  REAL vx, vy;
  REAL fx, fy;
} POINT;

POINT * init_point ( REAL _rx, REAL _ry, REAL _vx, REAL _vy );

void copy_point ( POINT * p1, POINT *p2 );

void free_point ( POINT * p);
