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


#include "md.h"

/*
#define PRINT_DIPOLE
#define PRINT_DEFIELD
#define PRINT_TRIATIC
#define PRINT_T_DIPOLE
*/

#define GETDIPDER


void getdid(SIMPARMS *simparms,COORDS *coords,INTER *inter,double **dmumat)
{

	getdid_new(simparms, coords, inter, dmumat);

}
/*----------------------------------------------------------------*/

void getddid(SIMPARMS *simparms,COORDS *coords,double **amat)
{

	getddid_new(simparms, coords, amat);

}

/*----------------------------------------------------------------*/

void getdiatic(SIMPARMS *simparms,COORDS *coords,
	       int i_now,int j_now,double **diatic)
{

	getdiatic_new(simparms, coords, i_now, j_now, diatic);

}
/*----------------------------------------------------------------*/

void test_diatic(SIMPARMS *simparms,COORDS *coords)
{

	test_diatic_new(simparms, coords);

}  

/*----------------------------------------------------------------*/

void getefield(SIMPARMS *simparms,COORDS *coords,double *bfield)
{

	getefield_new(simparms, coords, bfield);

}
/*---------------------------------------------------------------*/
void getdefield(SIMPARMS *simparms,COORDS *coords,INTER *inter,double **amat)
{

	getdefield_new(simparms, coords, inter, amat);

}
/*---------------------------------------------------------------*/
void getmu(SIMPARMS *simparms,COORDS *coords,double *mu)
{

	getmu_new(simparms, coords, mu);

}
/*----------------------------------------------------------------*/

void getdadq(SIMPARMS *simparms,COORDS *coords,double **amat)
{

	getdadq_new(simparms, coords, amat);

}

/*---------------------------------------------------------------*/
