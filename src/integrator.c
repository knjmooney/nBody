/********************************************************** 
 * Name   : integrator.c
 * Author : Kevin Mooney
 * Date   : 14/04/16
 *
 * Definitions of integrator.h
 *
 * Distance function doesn't calculate ghost distance
 *
 * Particles are only confined to a box when printing data
 **********************************************************/

#include <math.h>
#include <stdio.h>

#include "integrator.h"

/* Update a pair of forces */
/* To avoid calculating sqrts, you can deal wiht r^2 and rewrite E and F */
/* Can also avoid dividing by height and width all the time and just caluclating the inverse once */
inline void update_pairwise_force ( POINT *p1, POINT *p2 , REAL (*f)(REAL), double width, double height ) {
  REAL fx, fy, dx, dy, total_force, radial_seperation;

  /* Find closest dx, shamelessly stolen from wikipedia */
  dx  = p2->rx - p1->rx;
  dy  = p2->ry - p1->ry;
  dx -= width  * round ( dx / width  ); /* could precalculate */
  dy -= height * round ( dy / height );  

  radial_seperation = sqrt ( dx*dx + dy*dy );
  total_force       = - f ( radial_seperation ); /* force is in opposite direction */

  /* Could precalculate 1 / r */
  fx = total_force * ( dx ) / radial_seperation;
  fy = total_force * ( dy ) / radial_seperation;

  /* Append forces */
  p1->fx += fx;
  p1->fy += fy;
  p2->fx -= fx;
  p2->fy -= fy;
}

/* Update all the forces */
void update_forces ( GRID *g ) {
  POINT **p = g->p;
  int    i, j;

  /* reset forces */
  for ( i=0; i < g->size; i++ ) {
    p[i]->fx = 0;
    p[i]->fy = 0;
  }
    
  /* Recaluclate forces */
  for ( i=0; i < g->size; i++ ) {
    for ( j=i+1; j < g->size; j++ ) {
      update_pairwise_force ( p[i], p[j], g->f , g->width, g->height);
    }
  }
}

/* Calculates the total energy */
REAL total_energy ( GRID *g ) {
  REAL energy = 0.0;
  int i,j;
  POINT **p = g->p;

  REAL dx, dy;

  /* Add potentials */
  for ( i=0; i < g->size; i++ ) {
    for ( j=i+1; j < g->size; j++ ) {
      dx  = g->p[j]->rx - g->p[i]->rx;
      dy  = g->p[j]->ry - g->p[i]->ry;
      dx -= g->width  * round ( dx / g->width  ); /* could precalculate */
      dy -= g->height * round ( dy / g->height );

      energy += g->E( sqrt ( dx*dx + dy*dy ) );
    }
  }

  /* Add kinetic energy*/
  for ( i=0; i < g->size; i++ ) {
    energy += 0.5 * p[i]->vx*p[i]->vx + p[i]->vy*p[i]->vy;
  }
  return energy;
}

/* Confines all the particles into the one box */
void confine ( GRID *g ) {
  for ( int j=0; j < g->size; j++ ) {
    g->p[j]->rx -= g->width  * floor( g->p[j]->rx / g->width  );
    g->p[j]->ry -= g->height * floor( g->p[j]->ry / g->height );
  }
}

/* Print coordinates */
void print_coords(GRID *g, double energy) {
  for ( int j=0; j < g->size; j++ ) {
    printf ( "%10lf %10lf ", g->p[j]->rx, g->p[j]->ry );
  }
  printf("   %lf \n", ( energy - total_energy(g) ) / energy);
}

/* Leap frog integrator */
void leap_frog(GRID *g ,double h, int N ) {

  POINT **p = g->p;
  int i, j;
  /* REAL energy = total_energy(g); */

  update_forces ( g );

  /* Take half step */
   for ( i=0; i < g->size; i++ ) {
    p[i]->vx += 0.5 * h * p[i]->fx; 
    p[i]->vy += 0.5 * h * p[i]->fy; 
  }

  /* LEAP FROG */
  for ( i=0; i < N-1; i++ ) {

    update_forces ( g );

    for ( j=0; j < g->size; j++ ) {
      p[j]->rx += h * p[j]->vx;
      p[j]->ry += h * p[j]->vy;

      p[j]->vx += h * p[j]->fx; 
      p[j]->vy += h * p[j]->fy; 
    }    
    if ( i%1000 == 0 ) {
      /* confine the particles */
      confine ( g );
      /* print_coords(g,energy); */
    }
  }

  update_forces ( g );
  
  /* Last half step */
  for ( i=0; i < g->size; i++ ) {
    p[i]->rx += h * p[i]->vx;
    p[i]->ry += h * p[i]->vy;
    
    p[i]->vx += 0.5 * h * p[i]->fx; 
    p[i]->vy += 0.5 * h * p[i]->fy; 
  }
}

/* PARALLEL IMPLEMENTATION */

void add_to_halo ( PARA_GRID *g, POINT * p, enum direction way ) {
      g->halo_out [ way ][ g->halo_size[way][1] + g->halo_size[way][0]     ] = p->rx;
      g->halo_out [ way ][ g->halo_size[way][1] + g->halo_size[way][0] + 1 ] = p->ry;
      g->halo_size[ way ][0] += 2;
}

/* Must call init swap before init halo */
void init_halo ( PARA_GRID *g ) {
  POINT **p = g->p;
  int i;

  for ( i=0; i < 8; i++ ) {
    g->halo_size[i][0] = 0;
  }

  for ( i=0; i < g->size; i++ ) {
    REAL relx = p[i]->rx - g->origin[0];
    REAL rely = p[i]->ry - g->origin[1];
    
    if ( relx < g->delta ) {
      add_to_halo ( g, p[i], LEFT );
    }
    if ( rely < g->delta ) {
      add_to_halo ( g, p[i], DOWN );
    }
    if ( relx > g->lwidth  - g->delta ) {
      add_to_halo ( g, p[i], RIGHT );
    }
    if ( rely > g->lheight - g->delta ) {
      add_to_halo ( g, p[i], UP );
    }
  }  
}

/* Add the velocity and position of a point to the halo for swapping  */
void add_to_swap ( PARA_GRID *g, POINT * p, enum direction way ) {
      g->halo_out [ way ][ g->halo_size[way][1]     ] = p->rx;
      g->halo_out [ way ][ g->halo_size[way][1] + 1 ] = p->ry;
      g->halo_out [ way ][ g->halo_size[way][1] + 2 ] = p->vx;
      g->halo_out [ way ][ g->halo_size[way][1] + 3 ] = p->vy;
      g->halo_size[ way ][1] += 4;
}

/* Copies a point from the end of the list to the current index and frees the last element */
void delete_from_lgrid ( PARA_GRID *g, int index ) {
  g->size --;
  if ( index != g->size ) {
    copy_point ( g->p[index], g->p[g->size] );
  }
  free_point ( g->p[g->size] );
}

void add_to_lgrid ( PARA_GRID *g, REAL rx, REAL ry, REAL vx, REAL vy ) {
  g->p[g->size] = init_point ( rx, ry, vx, vy );
  g->size ++;
}

/* Appends the halo with particles to be swapped */
void init_swap ( PARA_GRID *g ) {
  POINT **p = g->p;
  int i;

  for ( i=0; i < 8; i++ ) {
    g->halo_size[i][1] = 0;
  }

  for ( i=0; i < g->size; i++ ) {
    REAL relx = p[i]->rx - g->origin[0];
    REAL rely = p[i]->ry - g->origin[1];
    
    if ( relx < 0.0) {
      add_to_swap ( g, p[i], LEFT );
      delete_from_lgrid ( g, i );
    }
    else if ( rely < 0.0 ) {
      add_to_swap ( g, p[i], DOWN );
      delete_from_lgrid ( g, i );
    }
    else if ( relx > g->lwidth ) {
      add_to_swap ( g, p[i], RIGHT );
      delete_from_lgrid ( g, i );
    }
    else if ( rely > g->lheight ) {
      add_to_swap ( g, p[i], UP );
      delete_from_lgrid ( g, i );
    }
  }  
}

void fix_ownership ( PARA_GRID *g ) {
  init_swap ( g );
  init_halo ( g );

  MPI_Request reqi[8], reqo[8];
  MPI_Status  stati[8], stato[8];

  for ( int i=0; i<8; i++ ) {
    MPI_Isend ( g->halo_size[i]   , 2, MPI_INT, g->dir[i], 0, g->comm, &reqo[i] );
    MPI_Irecv ( g->halo_in_size[i], 2, MPI_INT, g->dir[i], 0, g->comm, &reqi[i] );
  }

  for ( int i=0; i<8; i++ ) {
    MPI_Wait ( &reqo[i], &stato[i] );
    MPI_Wait ( &reqi[i], &stati[i] );
  }

  for ( int i=0; i<8; i++ ) {
    /* if ( g->halo_size[i][0] + g->halo_size[i][1] ) { */
    /*   printf ("%d->%d (%d,%d) : ",g->rank, g->dir[i], g->halo_size[i][0] , g->halo_size[i][1] ); */
    /*   for ( int j=0; j<g->halo_size[i][0] + g->halo_size[i][1]; j++ ) { */
    /* 	printf ( " %lf", g->halo_out[i][j] ); */
    /*   } */
    /*   printf ("\n"); */
    /* } */
    MPI_Isend ( g->halo_out[i], g->halo_size[i][0]    + g->halo_size[i][1]   , MPI_DOUBLE, 
		g->dir[i], 0, g->comm, &reqo[i] );

    MPI_Irecv ( g->halo_in[i] , g->halo_in_size[i][0] + g->halo_in_size[i][1], MPI_DOUBLE, 
		g->dir[i], 0, g->comm, &reqi[i] );
  }

  for ( int i=0; i<8; i++ ) {

    MPI_Wait ( &reqo[i], &stato[i] );
    MPI_Wait ( &reqi[i], &stati[i] );

    /* if ( g->halo_size[i][0] + g->halo_size[i][1] ) { */
    /*   printf ("%d->%d (%d,%d) : ",g->rank, g->dir[i], g->halo_size[i][0] , g->halo_size[i][1] ); */
    /*   for ( int j=0; j<g->halo_size[i][0] + g->halo_size[i][1]; j++ ) { */
    /* 	printf ( " %lf", g->halo_out[i][j] ); */
    /*   } */
    /*   printf ("\n"); */
    /* } */


    for ( int j=0; j < g->halo_in_size[i][1]; j+=4 ) {
      add_to_lgrid ( g, g->halo_in[i][j], g->halo_in[i][j+1], g->halo_in[i][j+2] ,g->halo_in[i][j+3] );
      /* printf ( "%d : %d\n", g->rank, g->dir[i] ); */
    }
  }


  return;
}

/* Update all the forces */
void para_update_forces ( PARA_GRID *g ) {
  POINT **p = g->p;
  int    i, j;

  fix_ownership ( g );

  /* reset forces */
  for ( i=0; i < g->size; i++ ) {
    p[i]->fx = 0;
    p[i]->fy = 0;
  }
    
  /* Recaluclate forces */
  for ( i=0; i < g->size; i++ ) {
    for ( j=i+1; j < g->size; j++ ) {
      update_pairwise_force ( p[i], p[j], g->f , g->width, g->height);
    }
  }
  for ( i=0; i < 8; i++ ) {
    for ( j=g->halo_in_size[i][1]; j < g->halo_in_size[i][1] + g->halo_in_size[i][0]; j+=2 ) {
      POINT *p = init_point ( g->halo_in[i][j], g->halo_in[i][j+1], 0, 0 );
      for ( int k = 0; k < g->size; k++ )
	update_pairwise_force ( g->p[k], p, g->f , g->width, g->height);
      free_point(p);
    }
  }

} 

/* Confines the particles to the binding box for printing */
void para_confine ( PARA_GRID *g ) {
  for ( int j=0; j < g->size; j++ ) {
    g->p[j]->rx -= g->width  * floor( g->p[j]->rx / g->width  );
    g->p[j]->ry -= g->height * floor( g->p[j]->ry / g->height );
  }
}

void write_to_file ( PARA_GRID *g, FILE *fp ) {
  fprintf ( fp ,"%10lf %10lf\n", g->origin[0], g->origin[1] );
  for ( int i=0; i < g->size; i++ )
    fprintf ( fp, "%10lf %10lf\n", g->p[i]->rx, g->p[i]->ry );
  fprintf ( fp, "\n\n\n" );
}

/* Leap frog integrator */
void para_leap_frog(PARA_GRID *g ,double h, int N ) {

  POINT **p = g->p;
  int i, j;
  FILE  *fp;
  char filename[100];

  sprintf ( filename, "data_par/parafile_%d.dat", g->rank );
  fp = fopen ( filename, "w" );
  write_to_file(g,fp);
      
  para_update_forces ( g );

  /* Take half step */
   for ( i=0; i < g->size; i++ ) {
    p[i]->vx += 0.5 * h * p[i]->fx; 
    p[i]->vy += 0.5 * h * p[i]->fy; 
  }

  /* LEAP FROG */
  for ( i=0; i < N-1; i++ ) {

    para_update_forces ( g );

    for ( j=0; j < g->size; j++ ) {
      p[j]->rx += h * p[j]->vx;
      p[j]->ry += h * p[j]->vy;

      p[j]->vx += h * p[j]->fx; 
      p[j]->vy += h * p[j]->fy; 
    }    
    /* if ( i%1000 == 0 ) { */
    /*   /\* confine the particles *\/ */
    /*   write_to_file(g,fp); */
    /* } */
    para_confine ( g );
  }

  /* Last half step */
  para_update_forces ( g );
  for ( i=0; i < g->size; i++ ) {
    p[i]->rx += h * p[i]->vx;
    p[i]->ry += h * p[i]->vy;
    
    p[i]->vx += 0.5 * h * p[i]->fx; 
    p[i]->vy += 0.5 * h * p[i]->fy; 
  }

  fclose(fp);
}
