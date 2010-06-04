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
#ifdef JUNK
int main(int argc,char *argv[])
{
   int i,nitem,nproc;
   nitem = atoi(argv[1]);
   nproc = atoi(argv[2]);

   for ( i=0; i<nproc; i++ ) {
      int d,n;
      decomp1d(nitem,nproc,i,&d,&n);
      printf("%3d:  %6d %6d\n",i,n,d);
   }

   return(0);
}
#endif

/*
 * MPI datatypes to be created automatically during initialization
 */
MPI_Datatype MPI_DCOMPLEX;

void My_MPI_Init(void)
/******************************************************************/
{
   MPI_Type_contiguous(2,MPI_DOUBLE,&MPI_DCOMPLEX);
   MPI_Type_commit(&MPI_DCOMPLEX);
}


void decomp1d(int nitem,int nproc,int rank,int *this_disp,
                     int *this_n)
/******************************************************************/
/* This function calculates the displacement into a list of nitem
 * objects for process rank, 0<=rank<nproc, in a group of nproc
 * processes. The number of items on this rank are also calculated
 * and returned in this_n. The decomposition places nitems/nproc
 * on all, and distributes the remainder across the first nitem%nproc.
 */
{
   int nb,nr;
   nb = nitem/nproc;
   nr = nitem%nproc;
   if ( rank<nr ) {
      *this_disp = rank*(nb+1);
      *this_n    = nb+1;
   } else {
      *this_disp = nr*(nb+1)+(rank-nr)*nb;
      *this_n    = nb;
   }
}

