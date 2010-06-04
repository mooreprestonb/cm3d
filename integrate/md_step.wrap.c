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

/*------------------------------------------------------*/

void step_init(SIMPARMS *simparms,COORDS *coords,ENERGIES *energies,
	       INTER *inter,NGBR *ngbr)
{

	v1fv3(simparms->natoms, coords->p, coords->px, coords->py, coords->pz);
	v1fv3(simparms->natoms, coords->v, coords->vx, coords->vy, coords->vz);
	v1fv3(simparms->natoms, coords->fs, coords->fxs, coords->fys, coords->fzs);
	v1fv3(simparms->natoms, coords->fl, coords->fxl, coords->fyl, coords->fzl);
	v1fv3(simparms->natoms, coords->fa, coords->fxa, coords->fya, coords->fza);
	v1fv3(simparms->natoms, coords->fr, coords->fxr, coords->fyr, coords->fzr);
	v1fv3(simparms->natoms, coords->ft, coords->fxt, coords->fyt, coords->fzt);

	step_init_new(simparms, coords, energies, inter, ngbr);

	v3fv1(simparms->natoms, coords->p, coords->px, coords->py, coords->pz);
	v3fv1(simparms->natoms, coords->v, coords->vx, coords->vy, coords->vz);
	v3fv1(simparms->natoms, coords->fs, coords->fxs, coords->fys, coords->fzs);
	v3fv1(simparms->natoms, coords->fl, coords->fxl, coords->fyl, coords->fzl);
	v3fv1(simparms->natoms, coords->fa, coords->fxa, coords->fya, coords->fza);
	v3fv1(simparms->natoms, coords->fr, coords->fxr, coords->fyr, coords->fzr);
	v3fv1(simparms->natoms, coords->ft, coords->fxt, coords->fyt, coords->fzt);


}
/*------------------------------------------------------*/

void step(SIMPARMS *simparms, COORDS *coords, ENERGIES *energies, INTER *inter, NGBR *ngbr) {


	v1fv3(simparms->natoms, coords->p, coords->px, coords->py, coords->pz);
	v1fv3(simparms->natoms, coords->v, coords->vx, coords->vy, coords->vz);
	v1fv3(simparms->natoms, coords->fs, coords->fxs, coords->fys, coords->fzs);
	v1fv3(simparms->natoms, coords->fl, coords->fxl, coords->fyl, coords->fzl);
	v1fv3(simparms->natoms, coords->fa, coords->fxa, coords->fya, coords->fza);
	v1fv3(simparms->natoms, coords->fr, coords->fxr, coords->fyr, coords->fzr);
	v1fv3(simparms->natoms, coords->ft, coords->fxt, coords->fyt, coords->fzt);

	step_new(simparms, coords, energies, inter, ngbr);

	v3fv1(simparms->natoms, coords->p, coords->px, coords->py, coords->pz);
	v3fv1(simparms->natoms, coords->v, coords->vx, coords->vy, coords->vz);
	v3fv1(simparms->natoms, coords->fs, coords->fxs, coords->fys, coords->fzs);
	v3fv1(simparms->natoms, coords->fl, coords->fxl, coords->fyl, coords->fzl);
	v3fv1(simparms->natoms, coords->fa, coords->fxa, coords->fya, coords->fza);
	v3fv1(simparms->natoms, coords->fr, coords->fxr, coords->fyr, coords->fzr);
	v3fv1(simparms->natoms, coords->ft, coords->fxt, coords->fyt, coords->fzt);


}
/*-----------------------------------------------------------*/
