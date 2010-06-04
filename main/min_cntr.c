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

/* routine to minimize a function in one dimension */

#include "md.h"

void mal_verify(int);
#define DXSTEP .01
#define IRSTEP 100
#define PRINT_FREQ

/* static varibales to go over the head of the minimize routines */
static NMODES nmodes;
static ENERGIES energies;
static SIMPARMS *simparms_s;
static COORDS *coords_s;
static INTER *inter_s;
static NGBR *ngbr_s;
static int ndim;
static double *svec;

double func(double *x);
double func_ir(double *x);
void dfunc(double *x,double *dx);
void print_potential(SIMPARMS *,COORDS *coords);
void load_fix(COORDS *coords,int *ndim, double **p);

/*-----------------------------------------------------------*/

void min_cntr(char *command,FILENAMES *filenames,SIMPARMS *simparms,
	      WRITE_STEP *write_step,COORDS *coords,INTER *inter,NGBR *ngbr)
{

	min_cntr_new(command, filenames, simparms, write_step, coords, inter, ngbr);

}
/*-------------------------------------------------------------------*/

#define KCOM 10
double func(double *x)
{

	return(func_new(x));

}
/*-------------------------------------------------------------------*/

void dfunc(double *x,double *y)
{

	dfunc_new(x, y);

}
  
/*-------------------------------------------------------------------*/
double func_ir(double *x)
{

	return(func_ir_new(x));

}

/*-----------------------------------------------------------*/

void print_potential(SIMPARMS *simparms,COORDS *coords)
{

	print_potential_new(simparms, coords);

}

/*-----------------------------------------------------------*/
void load_fix(COORDS *coords,int *ndim, double **p)
{

	load_fix_new(coords, ndim, p);

}
/*-----------------------------------------------------------*/
