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

/* subroutines to handle io from md program */

#include "md.h"

/*------------------------------------------------------------------*/
/* determine if i and j are on the same molecule */
int same_mole(int iatom,int jatom,int *itype_species,int *itype_molecules)
{

	return(same_mole_new(iatom, jatom, itype_species, itype_molecules));

}
/*--------------------------------------------------------------*/
/* open file for writing and check to make sure that the 
   file was opened correctly! */
FILE *cfopenw(char *name)
{

	return(cfopenw_new(name));

}

/*--------------------------------------------------------------*/
#ifndef cmalloc
void *cmalloc(size_t mem)
{

	cmalloc_new(mem);

}
#endif
/*--------------------------------------------------------------*/
#ifndef crealloc
void *crealloc(void *ptr,size_t mem)
{

	crealloc_new(ptr, mem);

}
#endif
/*--------------------------------------------------------------*/
void allocate_coords(SIMPARMS *simparms,COORDS *coords)
{

	allocate_coords_new(simparms, coords);

}
void free_coords(COORDS *coords)
{

	free_coords_new(coords);

}
/*--------------------------------------------------------------*/
void allocate_extend(SIMPARMS *simparms,COORDS *coords)
{
	allocate_extend_new(simparms, coords);

}
void free_extend(SIMPARMS *simparms,COORDS *coords)
{

	free_extend_new(simparms, coords);

}
/*--------------------------------------------------------------*/
void allocate_exttab(SIMPARMS *simparms,COORDS *coords,int ntable)
{

	allocate_exttab_new(simparms, coords, ntable);

}
void free_exttab(COORDS *coords)
{

	free_exttab_new(coords);

}

/*--------------------------------------------------------------*/
void rmass(int natoms,double *fx,double *fy,double *fz,double *amass)
{

	rmass_new(natoms, fx, fy, fz, amass);

}
/*---------------------------------------------------------------------*/
/* routine to make sure that anint give the correct numbers!!!!*/
void test_anint(void)
{

	test_anint_new();

}
/*------------------------------------------------------------------------*/

void free_all_pointers(SIMPARMS *simparms,COORDS *coords,
		       INTER *inter,NGBR *ngbr)
{

	free_all_pointers_new(simparms, coords, inter, ngbr);

}

