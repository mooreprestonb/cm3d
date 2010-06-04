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

/* output things that the user might want to know */

#include "md.h"
      
void output_simparm(SIMPARMS *simparms,WRITE_STEP *write_step,int istep);

/*-------------------------------------------------------------------*/
void output_subparm(SIMPARMS *simparms,WRITE_STEP *write_step,
		    SUBSPACE *subspace,int istep)
{
  int i,nr,ni,nt;
  LINE line;

  /* output stuff that the user might want to know */
  output_simparm(simparms,write_step,istep);

  nr = ni = nt = 0;
  sprintf(line,"----------- %d subvectors --------------",subspace->num_subst);
  md_stdout(line);
  sprintf(line,"     #  Vectors  Real  Imaginary\n");
  md_stdout(line);
  for(i=0;i<subspace->num_subst;i++){
    sprintf(line,"%6d %6d %6d %6d",i+1,subspace->num_vecsub_max[i],
	   subspace->num_real_max[i],subspace->num_imag_max[i]);
    md_stdout(line);
    nr += subspace->num_real_max[i];
    ni += subspace->num_imag_max[i];
    nt += subspace->num_vecsub_max[i];
  }
  sprintf(line,"----------------------------------------------");
  md_stdout(line);
  sprintf(line,"total  %6d %6d %6d",nt,nr,ni);
  md_stdout(line);
}
/*---------------------------------------------------------*/ 
