/********************************************************** 
 * Name   : integrator.h
 * Author : Kevin Mooney
 * Date   : 14/04/16
 *
 * Collection of integrators, currently only describes 
 * leap frog
 **********************************************************/

#pragma once

#include "grid.h"
#include "para_grid.h"

void leap_frog(GRID *g ,double h, int N );
void para_leap_frog(PARA_GRID *g ,double h, int N );

