/*   Copyright (C) 1997 1998 1999 2000 Dr. Preston B. Moore

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


void getddid(SIMPARMS *simparms,COORDS *coords,double **mat2);
void getamatrix(SIMPARMS *simparms,COORDS *coords,double **amat,int nef,INTER *inter,double *bvec);
void adjustAlphap_der(COORDS *coords,double **fmat,int n,int imode,double *alphap_der);
void getdadq(SIMPARMS *simparms,COORDS *coords,double **amat,double **fmat,int imode);
void getpolder(SIMPARMS *simparms,COORDS *coords,INTER *inter,NMODES *nmodes);
void getdRkdqt(COORDS *coords, double **fmat, int k, int imode, double *r);
void getSnkBond(COORDS *coords, int n, int k, double *s);
void ludcmp(double **amat,int n3,int *indx,double *);
void lubksb(double **amat,int n3,int *indx,double *mu);
void getefield(SIMPARMS *simparms,COORDS *coords,double *bfield);
void getewef(SIMPARMS *simparms,COORDS *coords,double *bvec);
void getmu(SIMPARMS *simparms,COORDS *coords,double *bfield);
void getdefield(SIMPARMS *simparms,COORDS *coords,INTER *,double **amat);
void test_diatic(SIMPARMS *simparms,COORDS *coords);

void getdid(SIMPARMS *simparms,COORDS *coords,INTER *inter,double **dmumat);
