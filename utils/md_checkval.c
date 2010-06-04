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


/* routine to loop over arguments and check for infinities and NaN's */
 
#include <stdlib.h>
#include <limits.h>
#include <math.h>

#ifdef NOCHECKVAL
int check_vals()
{
  return 0;
}
#else
/* #define _SUN_ */
#ifdef _SUN_
#include <varargs.h>
#else
#include <stdarg.h>
#endif

void md_warning(char *);

#ifndef DBL_MAX /* if limits.h does not define this, include it explicitly */
#define DBL_MAX  1.7976931348623157E+308  /* max decimal */
#endif

#ifdef _SUN_
int check_vals(va_alist) va_dcl
#else
int check_vals(int num,...)
#endif
{
#ifdef _SUN_
  int num;
  double t;
  va_list argptr;

  va_start(argptr);
  num = va_arg(argptr,int);

#else
  double t;
  va_list argptr;

  /* initialize argptr */
  va_start(argptr,num);
#endif

#ifdef IBM
  /* ibm doesn't conform */
  return 0;
#endif
  while(num--){
    t = va_arg(argptr,double);
    
    /* check for inf's */
    if(fabs(t) > DBL_MAX){
      md_warning("There are infinite values from the initial values!");
      return 1;
    }
    
    /* check for NaNs */
    if(t != t){
      md_warning("There are NaNs from the calculation of the initial values!");
      return 1;
    }
  }
  va_end(argptr);
  return 0;
}

#endif /* NOCHECKVAL */

/*---------------------------------------------------------*/
