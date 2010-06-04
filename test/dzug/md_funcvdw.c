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

/* subroutines for functions */
/* divides derivative by x at end for efficency */

#include <math.h>
#include <stdio.h>
#ifndef M_PI
#define M_PI     3.14159265358979323846
#endif

#include "md.h"
char line[MAXLINELEN];

/*-----------------------------------------------------------------------*/
void func_lj(double x,double *f,double *df,double *d2f,double eps,double sig)
{
  *f = 4.*eps*(pow(sig/x,12.)-pow(sig/x,6.));
  *df = 4.*eps*(6.*pow(sig,6.)*pow(x,-7.)-12.*pow(sig,12.)*pow(x,-13.));
  *d2f = 4.*eps*(156.*pow(sig,12.)*pow(x,-14.)-42.*pow(sig,6.)*pow(x,-8.));
  *df /= x;
}
void func_tab_pot(FILE *ftab, double x,double *f,double *df,double *d2f)
{
  double offset;

  fscanf(ftab,"%le %le %le %le \n",&offset,f,df,d2f);

  if(fabs(offset-x) > 0.0001) {
    sprintf(line,"expect distance (%e) and tabulated distance (%e) differ\n",x,offset  );
    sprintf(line,"%s for tabulated potential. Check max_dist in paramfile.\n",line);
    sprintf(line,"%s should be (max_dist^2+2 * min_dist^2/997)*997/999.\n",line);
    md_error(line);
  }
  *df /= x;
}

/*-----------------------------------------------------------------------*/
void func_lj96(double x,double *f,double *df,double *d2f,double eps,double sig)
{
  *f = (27./4.)*eps*(pow(sig/x,9.)-pow(sig/x,6.));
  *df = (27./4.)*eps*(6.*pow(sig,6.)*pow(x,-7.)-9.*pow(sig,9.)*pow(x,-10.));
  *d2f = (27./4.)*eps*(90.*pow(sig,9.)*pow(x,-11.)-42.*pow(sig,6.)*pow(x,-8.));
  *df /= x;
}

void func_lj64(double x,double *f,double *df,double *d2f,double eps,double sig)
{
  *f = 6.75*eps*(pow(sig/x,6.)-pow(sig/x,4.));
  *df = 6.75*eps*(-6.*pow(sig,6.)*pow(x,-7.)+4.*pow(sig,4.)*pow(x,-5.));
  *d2f = 6.75*eps*(-20.*pow(sig,4.)*pow(x,-6.)+42.*pow(sig,6.)*pow(x,-8.));
  *df /= x;
}

/*-----------------------------------------------------------------------*/
// The Dzugutov equation can be divided into two parts:
//
// V(r) = f1(r) + f2(r)
//
// Whether or not each part is included depends on the theta function, which
// is 0 or 1.
/*-----------------------------------------------------------------------*/

void func_dzug(double x,double *f,double *df,double *d2f,double eps, double sig)
{
  double r = x/sig; // x is in "reduced" units... this puts it in "real" units
  
  // A,B,a,b,c,d,m are constant for now  
  float adzug,bdzug,a,b,c,d,m;
  adzug = 5.82;
  bdzug = 1.280;
  a     = /*0.982585;*/1.87;
  b     = 1.94;
  c     = 1.10;
  d     = 0.27;
  m     = 16;

  // calculate theta (Heaviside) function's values for f1 and f2
  // each piece of the function is only included in 'f' when it's theta is 1
  int t1 = (r<=a);
  int t2 = (r<=b);

  double f1, df1, d2f1;
  double f2, df2, d2f2;
  
  if(t1)
  {
    double Aexp1 = adzug*exp(c/(r-a));  // these two terms are used frequently
    double rmB   = ( pow(r,-m)-bdzug );
    
    f1   = Aexp1*rmB;
    df1  = ( -Aexp1*m*pow(r,(-1-m))*( -(a*a*m)+(2*a*m*r) - r*(c +m*r - (bdzug*c*pow(r,m))) ) ) /
              ((a-r)*(a-r));

    // part1 & part2 just help d2f1 look cleaner... change if it is too inefficient
    double part1 = ( bdzug*c*(-2*a + c + 2*r) ) / pow((a-r),4);
    double part2 = pow(r,-m)*
                   ( ((c*c)/pow((a-r),4)) + ((m*(1+m))/(r*r)) +
                   (2*c*(r-(m*a)+(m*r)))/(r*pow((r-a),3)) );
    d2f1         = Aexp1*(-part1 + part2);
  } else {f1 = 0; df1 = 0; d2f1 = 0;}
  
  if(t2)
  {
    f2   = bdzug*exp( d/(r-b) );
    df2  = (-f2*d)/((b-r)*(b-r));
    d2f2 = (-df2*(-2*b + d + 2*r) ) / ((b-r)*(b-r));
  } else {f2 = 0; df2 = 0; d2f2 = 0;}
  
  *f   = eps*(f1 + f2);
  *df  = (eps/sig)*(df1  + df2);
  *d2f = (eps/(sig*sig))*(d2f1 + d2f2);
  *df /= x;
}

/*-----------------------------------------------------------------------*/
void func_will(double x,double *f,double *df,double *d2f,
	       double awill,double bwill,double c6,double c8,double c10)
{
  *f = awill*exp(-bwill*x) -c10*pow(x,-10.)-c8*pow(x,-8.)-c6*pow(x,-6.);
  *df = (-awill*bwill*exp(-bwill*x)
	 + 10.*c10*pow(x,-11.)+8.*c8*pow(x,-9.)+6.*c6*pow(x,-7.));
  *d2f = (awill*bwill*bwill*exp(-bwill*x) 
	  - 110.*c10*pow(x,-12.) - 72.*c8*pow(x,-10.) - 42.*c6*pow(x,-8.));

  *df /= x;
}
/*-----------------------------------------------------------------------*/
void func_aziz(double x,double *f,double *df,double *d2f,
	       double awill,double bwill,double c6,double c8,double c9,
	       double c10,double yamma,double drmexp,double gamma)
{
  double r2,r6,r8,r9,r10,arg,sw,disp,v1,v2,dv2,ddisp,dv1;
  double d2v2,d2v1,sw1,sw2,d2disp;

  r2 = x*x;
  r10 = r2*(r8 = r2*(r6 = r2*r2*r2));
  r9 = r10*x;
  if(x>drmexp){
    sw = 1.;
    sw1 = sw2 = 0.;
  } else {
    arg = drmexp/x - 1.;
    sw = exp(-arg*arg);
    sw1 = 2.*arg*drmexp*sw/r2;
    sw2 = 2.*drmexp*sw*(2.*arg*arg*drmexp-drmexp-2.*arg*x)/(r2*r2);
  }
  disp = -c6/r6-c8/r8-c9/r9-c10/r10;
  ddisp = (6.*c6/r6 + 8.*c8/r8 + 9.*c9/r9 + 10.*c10/r10)/x;
  d2disp = (-42.*c6/r6-72.*c8/r8-90.*c9/r9-110.*c10/r10)/r2;
  v1 = sw*disp;
  v2 = awill*pow(x,yamma)*exp(-bwill*x-gamma*r2);

  dv1 = sw*ddisp + sw1*disp;
  dv2 = (yamma/x-bwill-2.*gamma*x)*v2;

  d2v1 = 2.*sw1*ddisp+sw2*disp+sw*d2disp;
  d2v2 = v2*(bwill*bwill-2.*gamma+4.*bwill*gamma*x+4.*gamma*gamma*r2-
	     4.*gamma*yamma-gamma/r2-2.*bwill*yamma/x+yamma*yamma/r2);
  
  *f   = (v1+v2);
  *df  = (dv1+dv2);
  *d2f = (d2v1+d2v2);

  *df /= x;
}
/*-----------------------------------------------------------------------*/
void func_hydbnd(double x,double *f,double *df,double *d2f,
		 double c12,double d10)
{
  *f   = c12*pow(x,-12.)-d10*pow(x,-10.);
  *df  = 10.*d10*pow(x,-11.)-12.*c12*pow(x,-13.);
  *d2f = 156.*c12*pow(x,-14.)-110.*d10*pow(x,-12.);
  *df /= x;
}
/*-----------------------------------------------------------------------*/
void func_kerf(double x,double *f,double *df,double *d2f,double alpha)
{
  *f = erfc(alpha*x)/x;
  *df = -erfc(alpha*x)/(x*x)-2.*alpha*exp(-alpha*alpha*x*x)/(sqrt(M_PI)*x);
  *d2f = (2.*(2.*alpha*x+2.*alpha*alpha*alpha*x*x*x+
	      exp(alpha*alpha*x*x)*sqrt(M_PI)*erfc(alpha*x))/
	  (exp(alpha*alpha*x*x)*sqrt(M_PI)*x*x*x));
  *df /= x;
}
/*-----------------------------------------------------------------------*/
void func_erf(double x,double *f,double *df,double *d2f,double alpha)
{
  double ax,ax2,spi;
  
  spi = 1./sqrt(M_PI);
  ax = alpha*x;
  ax2 = ax*ax;
  *f = erf(ax)/x;
  *df = 2.*alpha*spi*exp(-ax2)/x - erf(ax)/(x*x);
  *d2f = (-4.*spi*pow(alpha,3.)*exp(-ax2)-4.*alpha*spi*exp(-ax2)/(x*x)
	  +2.*erf(ax)/(x*x*x));
  *df /= x;
}
/*-----------------------------------------------------------------------*/
void func_coul(double x,double *f,double *df,double *d2f)
{
  *f = 1./x;
  *df = -1./(x*x);
  *d2f = 2./(x*x*x);
  *df /= x;
}
/*-----------------------------------------------------------------------*/
void func_hautman(double x,double *f,double *df,double *d2f,
		  double zmin,double c12,double c3)
{
  double dz;
  
  dz = x - zmin;
  if(dz==0.0) dz = 1.e-10; /* a small number */

  *f = c12*pow(dz,-12.)-c3*pow(dz,-3.);
  *df = -12.*c12*pow(dz,-13.)+3.*c3*pow(dz,-4.);
  *d2f = 156.*c12*pow(dz,-14.)-12.*c3*pow(dz,-5.);
  /* WARNING this is a surface pot and we don't want to do a 
   *df /= x; which is usefull for other potentials */
}

