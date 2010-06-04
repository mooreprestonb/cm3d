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

void write_vectors(int natoms,double *x0,double *y0,double *z0,double *amass,
		   SUBSPACE *subspace,double **d2v,double *dn,double cell)
{
  int i,j,l;
  WORD name;
  double factor,mass;
  FILE *fp;

  /* calculate some normalization factor */
  factor = 0.;
  for(i=0;i<natoms;i++){
    factor += sqrt(amass[i]);
  }
  factor /= cbrt((double)natoms);

  for(i=0;i<subspace->nstate;i++){
    sprintf(name,"%s_%d.%.2f",subspace->vecfile,i,dn[i]);
    fp=cfopenw(name);
    fprintf(fp,"%d 1 1\n",natoms);
    for(l=j=0;j<natoms;j++,l+=3){
      mass = factor/sqrt(amass[j]);
      fprintf(fp,"%9.5g %9.5g %9.5g %9.5g %9.5g %9.5g\n",x0[j],y0[j],z0[j],
	      d2v[l  ][i]*mass,d2v[l+1][i]*mass,d2v[l+2][i]*mass);
    }
    fclose(fp);
  }
  /* print header */

  sprintf(name,"%s.init",subspace->vecfile);
  fp=cfopenw(name);

  fprintf(fp,"%d %d\n",natoms,1);

  for(i=0;i<natoms;i++){
    fprintf(fp,"%g %g %g\n",x0[i],y0[i],z0[i]);
  }
  fprintf(fp,"%g 0. 0.\n0. %g 0.\n0. 0. %g\n",cell,cell,cell);
  fclose(fp);
}
