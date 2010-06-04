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


/* zero a double vector */

#include <math.h>
#include <stdlib.h>

void md_error(char *);

void azero(double x[],const int n)
{

	azero_new(x, n);

}
/*-----------------------------------------------------------*/
/* take projected array and expand them back into the full space */
void expandztox(double *x,double *y,double *z,double *q,double **d2v,
		int m,int n)
{

	expandztox_new(x, y, z, q, d2v, m, n);

}

/*-----------------------------------------------------------*/
/* take unprojected array and project them back into the subspace */
void projectxtoz(double *x,double *y,double *z,double *q,double **d2v,
		 int m,int n)
{

	projectxtoz_new(x, y, z, q, d2v, m, n);

}
/*-------------------------------------------------------------*/
void mass_weight(int n,double *x,double *y,double *z,double *amass,int iflag)
{

	mass_weight_new(n, x, y, z, amass, iflag);

}
/*----------------------------------------------------------------*/
