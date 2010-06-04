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

/* routines that zero total momentum, scale temperature, zero COM, etc.. */

#include "md.h"

/* #define DEBUG */
/*----------------------------------------------------------------------*/
void samvel(SIMPARMS *simparms,COORDS *coords,int irandom)
{
	samvel_new(simparms, coords, irandom);

}
/*--------------------------------------------------------------------*/
void zerotm(int natoms,double *vx,double *vy,double *vz,double *amass)
{
	int i, ii;
	double *v;

	v = calloc(3*natoms, sizeof(double));
	/* interleave */
	for(i = 0; i < natoms; i++) {

		ii = i*3;
		v[ii+0] = vx[i];
		v[ii+1] = vy[i];
		v[ii+2] = vz[i];

	}

	zerotm_new(natoms, v, amass);

	free(v);

}
/*---------------------------------------------------------------------*/
double get_temp(int natoms,double *vx,double *vy,double *vz,
		double *amass,int ndof)
{


	int i, ii;
	double *v;

	v = calloc(3*natoms, sizeof(double));
	/* interleave */
	for(i = 0; i < natoms; i++) {

		ii = i*3;
		v[ii+0] = vx[i];
		v[ii+1] = vy[i];
		v[ii+2] = vz[i];

	}

	get_temp_new(natoms, v, amass, ndof);

	free(v);

}
/*---------------------------------------------------------------------*/
void scale(int natoms,double *vx,double *vy,double *vz,
	   double *amass,double temp,int ndof)
{

	int i, ii;
	double *v;

	v = calloc(3*natoms, sizeof(double));
	/* interleave */
	for(i = 0; i < natoms; i++) {

		ii = i*3;
		v[ii+0] = vx[i];
		v[ii+1] = vy[i];
		v[ii+2] = vz[i];

	}

	scale_new(natoms, v, amass, temp, ndof);

	free(v);

}

/*-----------------------------------------------------------*/
void radial(int natoms,double *px,double *py,double *pz,
	    double *vx,double *vy,double *vz,double *amass)
{


	int i, ii;
	double *v, *p;

	v = calloc(3*natoms, sizeof(double));
	p = calloc(3*natoms, sizeof(double));
	/* interleave */
	for(i = 0; i < natoms; i++) {

		ii = i*3;
		v[ii+0] = vx[i];
		v[ii+1] = vy[i];
		v[ii+2] = vz[i];
		p[ii+0] = px[i];
		p[ii+1] = py[i];
		p[ii+2] = pz[i];

	}

	radial_new(natoms, p, v, amass);

	free(v);
	free(p);

}
/*-----------------------------------------------------------*/
void zerocm(int natoms,double *px,double *py,double *pz,double *amass)
{
	int i, ii;
	double *p;

	p = calloc(3*natoms, sizeof(double));
	/* interleave */
	for(i = 0; i < natoms; i++) {

		ii = i*3;
		p[ii+0] = px[i];
		p[ii+1] = py[i];
		p[ii+2] = pz[i];

	}

	zerocm_new(natoms, p, amass);

	free(p);

}
/*--------------------------------------------------------------*/
void zero_angm(int natoms,double *px,double *py,double *pz,
	       double *vx,double *vy,double *vz,double *amass)
{


	int i, ii;
	double *v, *p;

	v = calloc(3*natoms, sizeof(double));
	p = calloc(3*natoms, sizeof(double));
	/* interleave */
	for(i = 0; i < natoms; i++) {

		ii = i*3;
		v[ii+0] = vx[i];
		v[ii+1] = vy[i];
		v[ii+2] = vz[i];
		p[ii+0] = px[i];
		p[ii+1] = py[i];
		p[ii+2] = pz[i];

	}

	zero_angm_new(natoms, p, v, amass);

	free(v);
	free(p);

}
/*------------------------------------------------------------------*/
void rotate(int natoms,double *px,double *py,double *pz,double *amass)
{

	int i, ii;
	double *p;

	p = calloc(3*natoms, sizeof(double));
	/* interleave */
	for(i = 0; i < natoms; i++) {

		ii = i*3;
		p[ii+0] = px[i];
		p[ii+1] = py[i];
		p[ii+2] = pz[i];

	}

	rotate_new(natoms, p, amass);

	free(p);

}
/*------------------------------------------------------------------*/
