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

#ifdef SERIAL
#  error "SERIAL not implemented yet!"
#else
#  include <mpi.h>
#endif

/*
 * MPI datatype definitions that I set up automatically
 * with initialization.
 */
extern MPI_Datatype MPI_DCOMPLEX;

/*
 * Macro to insure proper initialization of my MPI datatypes
 */
void My_MPI_Init(void);
#define MPI_Init(argc,argv) \
   ( MPI_Init(argc,argv), My_MPI_Init())


void decomp1d(int nitem,int nproc,int rank,int *this_disp,
                     int *this_n);


#endif
