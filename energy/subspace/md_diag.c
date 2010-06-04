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

/* routines dealing with the diagonalization of a matrix */

#include "md.h"

/*--------------------------------------------------------------*/
void diagonalize(double **fcmat,int nmol3,double *dn,
   double freq_min,double freq_max,int *num_real,int *num_imag,int kdiag)
{
  if(kdiag==0){
    rs_me(nmol3,dn,fcmat,1);
    ceigsrtv(dn,fcmat,nmol3);
  }
  if(kdiag==1){
    printf("freq_min = %g, freq_max = %g, num_real = %d, num_imag = %d\n",
	   freq_min,freq_max,*num_real,*num_imag);
    fprintf(stderr,"LAN NOT IMPLIMENTED!!\n");
    exit(1);
  }  
}
/*--------------------------------------------------------------*/

