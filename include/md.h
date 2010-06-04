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


/* file md.h which has common definitions for the md code mike */

#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>

/* #define DEBUG */

#include "md_defines.h"
#include "md_numbers.h"
#include "md_macros.h"
#include "md_typedef.h"
#include "md_proto.h"

/* #define PARA */

#ifdef PARA
#include "mpi.h"
#endif

/* #define MEMDEBUG */
#ifdef MEMDEBUG
/* #include "malloc.h" */
#include "memdebug.h"
#endif

#ifdef DMALLOC
#include "/home/moore/lib/include/dmalloc.h"
#endif

#ifdef SGI
#include <scsl_blas.h>
#endif
