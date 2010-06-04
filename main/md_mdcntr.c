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

/* controling subroutine for molecular dynamcis 
   using nose thermostats (NVT) and extened system variables
   for the volume
   */

#include "md.h"

/* #define USER_FUNC */

/* might get different numbers from the eigensolver */
/* #define NO_ROTATE */
#define FREEZEFORCE
#ifdef FREEZEFORCE
			   void save_freezeforce(int,COORDS *coords);
#endif

/*------------------------------------------------------------------*/
void md_mdcntr(char *command,FILENAMES *filenames,SIMPARMS *simparms,
	       WRITE_STEP *write_step,COORDS *coords,INTER *inter,
	       NGBR *ngbr)
{

	md_mdcntr_new(command, filenames, simparms, write_step, coords, inter, ngbr);

}

/*---------------------------------------------------------------------*/
#ifdef PARA
void bcast_system(SIMPARMS *simparms,COORDS *coords)
{

	bcast_system_new(simparms, coords);

}
#endif
/*---------------------------------------------------------------------*/

void gettngbr(NGBR *ngbr)
{

	gettngbr_new(ngbr);

}
/*---------------------------------------------------------------------*/
