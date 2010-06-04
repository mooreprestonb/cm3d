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

/* routine have to do with pressure and volume */

#include "md.h"
/* #define ISO */
/* #define ZERO_TM */
/*------------------------------------------------------------------*/
void set_baro(SIMPARMS *simparms,COORDS *coords,double tauv)
{

	set_baro_new(simparms, coords, tauv);

}
/*----------------------------------------------------------------*/
void forceV(SIMPARMS *simparms,COORDS *coords,ENERGIES *energies)
{

	forceV_new(simparms, coords, energies);

}
/*----------------------------------------------------------------*/
void samvbar(SIMPARMS *simparms,COORDS *coords)
{

	samvbar_new(simparms, coords);

}

/*----------------------------------------------------------------*/
void getengv(SIMPARMS *simparms,COORDS *coords,double *zkev,double *potv)
{

	getengv_new(simparms, coords, zkev, potv);

}

/*--------------------------------------------------------------------*/
void getengv_cheese(SIMPARMS *simparms,COORDS *coords,double *zkev,
		    double *potv,double *zken)
{

	getengv_cheese_new(simparms, coords, zkev, potv, zken);

}

/*--------------------------------------------------------------------*/
void intbar(SIMPARMS *simparms,COORDS *coords,double dta2,int iflag)
{

	intbar_new(simparms, coords, dta2, iflag);

}
/*--------------------------------------------------------------------*/
void app_ux2(int natoms,double dta2,COORDS *coords)
{

	app_ux2_new(natoms, dta2, coords);

}

/*----------------------------------------------------------------*/
void app_up2(int natoms,double dta2,double dnf,COORDS *coords)
{

	app_up2_new(natoms, dta2, dnf, coords);

}

/*----------------------------------------------------------------*/
void app_uh(double dta2,COORDS *coords)
{

	app_uh_new(dta2, coords);

}
