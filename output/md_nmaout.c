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

/* subroutines to deal with output of nma */

/* #define MOL_EIGEN */
#include "md.h"
#define ZERO /* remove zero frequencies from freq store */


void write_nma_start(int natoms,int npt,double dt,double frmax,
		     double frmin,int iunits,int nsave,NMODES *nmodes)
{

	write_nma_start_new(natoms, npt, dt, frmax, frmin, iunits, nsave, nmodes);

}

/*------------------------------------------------------------------*/
void save_nma(SIMPARMS *simparms,FILENAMES *filenames,int npoints,
	      double freq_min,NMODES *nmodes)
{

	save_nma_new(simparms, filenames, npoints, freq_min, nmodes);

}
/*-----------------------------------------------------------------------*/
void print_mat(int nr,int nc,double **matrix)
{

	print_mat_new(nr, nc, matrix);

}
/*------------------------------------------------------------------------*/
void psysvec(FILENAMES *filenames,int nconf,int natoms,int iunits,int npoints,
	     double *px,double *py,double *pz,double *amass,
	     double *hmat,double **fcmat,double *dn,double **molmat,
	     int nspec,int *nmol,int *napm)
{  

	psysvec_new(filenames, nconf, natoms, iunits, npoints, px, py, pz, amass, hmat, fcmat, dn, molmat, nspec, nmol, napm);

}
/*---------------------------------------------------------------------*/
void store_freq(SIMPARMS *simparms,int nmodes,double *dn,double time)
{
	store_freq_new(simparms, nmodes, dn, time);

}

/*---------------------------------------------------------------------*/
/* treat imaginary modes as real modes at the magnitude of their imaginary frequency */
void imag2real(int npoints,double df,double frmin,double a1[],double a2[])
{

	imag2real_new(npoints, df, frmin, a1, a2);

}

/*-------------------------------------------------------------------*/
/* calculate the participation ratio
   per Space,Rabitz,Askar JCP vol 99, 1993, pp 9070-9079  */

void partratio(double **fcmat,double *dn,double *part,int natoms)
{

	partratio_new(fcmat, dn, part, natoms);

}
                                         
                                         
