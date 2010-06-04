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


void md_startup(FILENAMES *,SIMPARMS *,SUBSPACE *,
		WRITE_STEP *,INTER *,NGBR *,COORDS *);

void md_mdcntr(char *,FILENAMES*,SIMPARMS *,WRITE_STEP *,COORDS *,
	       INTER *,NGBR *);

void colvarcntr(char *,FILENAMES*,SIMPARMS *,WRITE_STEP *,COORDS *,
		INTER *,NGBR *);

void bcast_system(SIMPARMS *simparms,COORDS *coords);
void gettngbr(NGBR *ngbr);

void md_subcntr(char *,FILENAMES *,SIMPARMS *,WRITE_STEP *,COORDS *,
		SUBSPACE *,INTER *,NGBR *);

void nma_control(char *command,FILENAMES *,SIMPARMS *simparms,COORDS *coords,
		 INTER *,WRITE_STEP *write_step,SUBSPACE *subspace);
void nma_allocate(SIMPARMS *simparms,NMODES *nmodes);

void min_cntr(char*,FILENAMES*,SIMPARMS*,WRITE_STEP*,COORDS*,INTER*,NGBR *);


/* --------------------------------sub_proj.c ----------------*/
void init_subproj(SIMPARMS *simparms);
void sub_proj(SIMPARMS *simparms,COORDS *coords,INTER *inter,NMODES *);




void md_startup_new(FILENAMES *,SIMPARMS *,SUBSPACE *,
		WRITE_STEP *,INTER *,NGBR *,COORDS *);

void md_mdcntr_new(char *,FILENAMES*,SIMPARMS *,WRITE_STEP *,COORDS *,
	       INTER *,NGBR *);

void colvarcntr_new(char *,FILENAMES*,SIMPARMS *,WRITE_STEP *,COORDS *,
		INTER *,NGBR *);

void bcast_system_new(SIMPARMS *simparms,COORDS *coords);
void gettngbr_new(NGBR *ngbr);

void md_subcntr_new(char *,FILENAMES *,SIMPARMS *,WRITE_STEP *,COORDS *,
		SUBSPACE *,INTER *,NGBR *);

void nma_control_new(char *command,FILENAMES *,SIMPARMS *simparms,COORDS *coords,
		 INTER *,WRITE_STEP *write_step,SUBSPACE *subspace);
void nma_allocate_new(SIMPARMS *simparms,NMODES *nmodes);

void min_cntr_new(char*,FILENAMES*,SIMPARMS*,WRITE_STEP*,COORDS*,INTER*,NGBR *);


/* --------------------------------sub_proj.c ----------------*/
void init_subproj_new(SIMPARMS *simparms);
void sub_proj_new(SIMPARMS *simparms,COORDS *coords,INTER *inter,NMODES *);

