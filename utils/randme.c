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

/* code that returns a random number between 0-1 */

#include <stdlib.h>
#include <math.h>

/* proto_types */
double ran1(long *);
double randme(void);
void srandme(int);

/* #define RAND  */
/* #define RANDOM  */
/* #define DRAND48 */
#define NR 

#ifdef RAND
#define MAXRND 32767 /* 2^15 -1 */
double randme(void){ return (rand()/(double)MAXRND);}
void srandme(int iseed){  srand(iseed);}
#endif
/*------------------------------------------------------*/
#ifdef RANDOM
#define MAXRND 2147483647  /* 2^31 - 1 */
double randme(void){ return ((double)random()/(double)MAXRND);}
void srandme(int iseed) { srandom(iseed);}
#endif
/*------------------------------------------------------*/
#ifdef DRAND48
double randme(void){  return drand48();}
void srandme(int iseed){ srand48((long) iseed);}
#endif
/*------------------------------------------------------*/
#ifdef NR

static long seed;
double randme(void)
{return ran1(&seed);}

void srandme(int iseed)
{
  if(iseed < 0 ) seed = (long) iseed;
  else seed = (long) -iseed;
  ran1(&seed);
}

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-12
#define RNMX (1.0-EPS)

double ran1(long *idum)
{
  int j;
  long k;
  static long iy=0;
  static long iv[NTAB];
  double temp;	

  if(*idum <=0 || !iy){
    if(-(*idum)<1) *idum=1;
    else *idum = -(*idum);
    for(j=NTAB+7;j>=0;j--){
      k=(*idum)/IQ;
      *idum=IA*(*idum-k*IQ)-IR*k;
      if(*idum < 0) *idum += IM;
      if(j <NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ;
  *idum=IA*(*idum-k*IQ)-IR*k;
  if(*idum < 0) *idum += IM;
  j=iy/NDIV;
  iy=iv[j];
  iv[j] = *idum;
  if((temp=AM*iy) > RNMX) return RNMX;
  else return temp;
}
#endif
    


