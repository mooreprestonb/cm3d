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

#define IEEE
#ifdef CRAY
#undef IEEE
#endif
#ifdef LINUX
#undef IEEE
#endif

#ifndef MAX
#define MAX(A,B) ((A>B)?(A):(B))
#endif

#ifndef MIN
#define MIN(A,B) ((A<B)?(A):(B))
#endif

#ifndef TZ
#define TZ(A) ((A==0.0)?(1):(A))
#endif

#define MAGIC 6755399441055744.0  /* 3 * 2^51*/
#define MNINT(X) (((X)+MAGIC)-MAGIC)

#ifndef anint
#ifndef IEEE
#define anint(A) ((int)((A)+((A)>=0?.5:-.5)))
#else
#define anint(X) MNINT(X)
#endif
#endif

#define COPY_SIGN(X,Y) ((Y>=0)?fabs(X):-fabs(X))

#ifndef cbrt
#define cbrt(A) (pow((A),1./3.))
#endif
