/* 
   Copyright (C) 1997 1998 1999 2000 Dr. Preston B. Moore

Dr. Preston B. Moore
Associate Director, Center for Molecular Modeling (CMM)
University of Pennsylvania, Department of Chemistry, Box 188 
231 S. 34th St. Philadelphia, PA 19104-6323 USA
EMAIL: moore@cmm.chem.upenn.edu  
WWW: http://www.cmm.upenn.edu/~moore

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*/

/* subroutine to time processes */

#include <sys/types.h>
#include <time.h>
#include <stdio.h>
#define MAXTIME 2147.483648
#define MINTIME (-MAXTIME)

#ifdef PARA
#include <mpi.h>
#define WALL
#endif

#ifndef CLOCKS_PER_SEC
#define CLOCKS_PER_SEC  1000000
#endif

double cputime(void)
{
  static clock_t itime=0;
  static double to=0.,tn=0.,dtime;
#ifdef WALL
  time_t t_time=1;
#endif

#ifdef WALL
#ifdef PARA
  if(!itime){
    itime = 1;
    tn = to = (double)MPI_Wtime();
  } else {
    tn = (double)MPI_Wtime();
  }
  dtime = tn-to;
#else 
  if(!itime){
    itime = 1;
    tn = to = (double)time(&t_time);
  } else {
    tn = (double)time(&t_time);
  }
  dtime = tn-to;
#endif
#else
  itime = clock();
  tn = (double)((double)itime/(double)CLOCKS_PER_SEC);
  
  if(to>tn) dtime = (MAXTIME-to)+(tn-MINTIME);
  else dtime = tn-to;
#endif

  to = tn;
  return dtime;
}
/*-------------------------------------------------*/
double realtime(void)
{
  static double dtime;
  time_t t_time=1;

#ifdef PARA
  dtime = (double)MPI_Wtime();
#else 
  dtime = (double)time(&t_time);
#endif

  return dtime;
}
