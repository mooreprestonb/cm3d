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

/* routine to preform periodic boundry conditions */

#include "md.h"

void period_new(int natoms,double dxs[],double dys[],double dzs[],
	    double hmat[],double hmati[],int iperd,int iensemble)
{
  int jpart;
  double sx,sy,sz;

  switch(iperd){
    /*=================================================================*/
     case 3:        /* Three dimensions */
#include "vectorize.h"
       for(jpart=0;jpart < natoms; ++jpart) {
         sx = dxs[jpart]*hmati[0]+dys[jpart]*hmati[1]+dzs[jpart]*hmati[2];
         sy = dxs[jpart]*hmati[3]+dys[jpart]*hmati[4]+dzs[jpart]*hmati[5];
         sz = dxs[jpart]*hmati[6]+dys[jpart]*hmati[7]+dzs[jpart]*hmati[8];
         sx -= anint(sx);
         sy -= anint(sy);
         sz -= anint(sz);
         dxs[jpart] = sx*hmat[0]+sy*hmat[1]+sz*hmat[2];
         dys[jpart] = sx*hmat[3]+sy*hmat[4]+sz*hmat[5];
         dzs[jpart] = sx*hmat[6]+sy*hmat[7]+sz*hmat[8];
       }
       break;
       /*==============================================================*/
     case 1:   /* One dimensions */
#include "vectorize.h"
       for(jpart=0;jpart < natoms; ++jpart) {
         sx  = dxs[jpart]*hmati[0];
         sx -= anint(sx);
         dxs[jpart] = sx*hmat[0];
       }
       break;
       /*============================================================*/
     case 2: /* Two dimensions */
#include "vectorize.h"
       for(jpart=0;jpart < natoms; ++jpart) {
         sx = dxs[jpart]*hmati[0]+dys[jpart]*hmati[1];
         sy = dxs[jpart]*hmati[3]+dys[jpart]*hmati[4];
         sx -= anint(sx);
         sy -= anint(sy);
         dxs[jpart] = sx*hmat[0]+sy*hmat[1];
         dys[jpart] = sx*hmat[3]+sy*hmat[4];
       }/* endfor */
       break;
  }/* end switch */
}
/*---------------------------------------------------------------------*/
#ifdef HALF
       if (iensemble < 3) {
         double hmxx,hmyy,hmzz;
         double hmxxi,hmyyi,hmzzi,hmxxhalf,hmyyhalf,hmzzhalf;
         hmxx=hmat[0]; hmyy=hmat[4]; hmzz=hmat[8];
         hmxxhalf=0.5*hmat[0]; hmyyhalf=0.5*hmat[4]; hmzzhalf=0.5*hmat[8];
         
         for(jpart=0;jpart < natoms; ++jpart) {
           
           while(dxs[jpart] > hmxxhalf) dxs[jpart] -= hmxx ;
           while(dxs[jpart] < -hmxxhalf) dxs[jpart] += hmxx ;   
           
           while(dys[jpart] > hmyyhalf) dys[jpart] -= hmyy ;
           while(dys[jpart] < -hmyyhalf ) dys[jpart] += hmyy ;
           
           while(dzs[jpart] > hmzzhalf) dzs[jpart] -= hmzz ;
           while(dzs[jpart] < -hmzzhalf) dzs[jpart] += hmzz ;
         }
       } else {
#endif 

