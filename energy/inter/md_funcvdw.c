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

void func_lj124(double x,double *f,double *df,double *d2f,double eps,double sig)
{
  *f = (1.5*sqrt(3.))*eps*(pow(sig/x,12.)-pow(sig/x,4.));
  *df =(1.5*sqrt(3.))*eps*(4.*pow(sig,4.)*pow(x,-5.)-12.*pow(sig,12.)*pow(x,-13.));
  *d2f = (1.5*sqrt(3.))*eps*(156.*pow(sig,12.)*pow(x,-14.)-20.*pow(sig,4.)*pow(x,-6.));
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
// The Dzugutov equation can be divided into two parts:
//
// V(r) = f1(r) + f2(r)
//
// Whether or not each part is included depends on the theta function, which
// is 0 or 1.
/*-----------------------------------------------------------------------*/

void func_dzug(double x,double *f,double *df,double *d2f,double eps, double sig)
{
  int t1,t2;
  
  double r = x/sig; // x is in "reduced" units... this puts it in "real" units
  
  // A,B,a,b,c,d,m are constant for now  
  float adzug,bdzug,a,b,c,d,m;
  double f1, df1, d2f1;
  double f2, df2, d2f2;
  double part1,part2;
  
  adzug = 5.82;
  bdzug = 1.280;
  a     = /*0.982585;*/1.87;
  b     = 1.94;
  c     = 1.10;
  d     = 0.27;
  m     = 16;
  t1 = (r<=a);
  t2 = (r<=b);

  // calculate theta (Heaviside) function's values for f1 and f2
  // each piece of the function is only included in 'f' when it's theta is 1
  
  if(t1) {
    /* terms used frequently */
    double Aexp1 = adzug*exp(c/(r-a));  
    double rmB   = ( pow(r,-m)-bdzug );
    
    f1   = Aexp1*rmB;
    df1  = ( -Aexp1*m*pow(r,(-1-m))*
	     ( -(a*a*m)+(2*a*m*r) - r*(c +m*r - (bdzug*c*pow(r,m))) ) ) /
      ((a-r)*(a-r));
    
    /* part1 & part2 just help d2f1 look cleaner... 
       change if it is too inefficient 
    */
    part1 = ( bdzug*c*(-2*a + c + 2*r) ) / pow((a-r),4);
    part2 = pow(r,-m)*
                   ( ((c*c)/pow((a-r),4)) + ((m*(1+m))/(r*r)) +
		     (2*c*(r-(m*a)+(m*r)))/(r*pow((r-a),3)) );
    d2f1         = Aexp1*(-part1 + part2);
  } else {
    f1 = 0; df1 = 0; d2f1 = 0;
  }
  
  if(t2) {
    f2   = bdzug*exp( d/(r-b) );
    df2  = (-f2*d)/((b-r)*(b-r));
    d2f2 = (-df2*(-2*b + d + 2*r) ) / ((b-r)*(b-r));
  } else {
    f2 = 0; df2 = 0; d2f2 = 0;
  }
  
  *f   = eps*(f1 + f2);
  *df  = (eps/sig)*(df1  + df2);
  *d2f = (eps/(sig*sig))*(d2f1 + d2f2);
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

/*-------------------------------------------------------------------------*/
void func_sphere_atom(double d0,double r,double *f,double *df,double *d2f,double *dfr,double eps,double sig)
{
  int cgwater = 1; /* 1 = true, 0 = false */
  double t1,t2,t3,t4,t5,t6,t7,t8,t9;
  double a,b;
  double d;
  d= d0-r; /*atom to the sphere surface*/ 
  if (d<0.0){
    fprintf(stderr,"Error: particle inside sphere\n");
  }
  if(cgwater) {
    a=0.0338;
    b=0.118;
    t1 = 4.0*a*eps*pow(sig,9.0)*pow(r,3.0)/5.0/pow(d,6.0)/(r+d)/pow((2.0*r+d),6.0); 
    t2 = 35.0*pow(d,4.0)+140.0*r*pow(d,3.0)+252.0*r*r*d*d+224.0*pow(r,3.0)*d+80.0*pow(r,4.0); 
    t3 = 8.0*b*eps*pow(sig,6.0)*pow(r,3.0)/pow(d,3.0)/pow((2.0*r+d),3.0); 
    t4 = 105.0*pow(d,6.0) + 630.0*pow(d,5.0)*r + 1764.0*pow(d,4.0)*r*r + 2856.0*pow(d,3.0)*pow(r,3.0) + 2736.0*d*d*pow(r,4.0) + 1440.0*d*pow(r,5.0) + 320.0*pow(r,6.0); 
    t5 = 7.0*pow(d0,6.0) + 35.0*pow(d0,4.0)*r*r + 21.0*d0*d0*pow(r,4.0) + pow(r,6); 
    t6 = d0*d0-r*r; 

    
    *f = t1*t2 - t3; 
    *df = (-3.0)*t1*t4/(d*(d0)*(d0+r)) + 6.0*t3*d0/(d*(d0+r)); 
    *dfr = 12.0*a*eps*pow(sig,9.0)*pow(r,2.0)/(d0*pow(t6,7.0))*t5 - 24.0*b*eps*pow(sig,6.0)*r*r*(d0*d0+r*r)/pow(t6,4.0); 

    } else {
  double rho = 0.113;
  t1 = 16.0*3.14159265358979323*rho*eps*pow(r,3.0);
  t2 = 2.0*r+d;
  t3 = 15.0*pow(d,6.0)+90.0*r*pow(d,5.0)+288.0*r*r*pow(d,4.0)+552.0*pow(r,3.0)*pow(d,3.0)+648.0*pow(r,4.0)*d*d+432.0*pow(r,5.0)*d+128.0*pow(r,6.0);
  t4 = 5.0*pow(d,7.0) + 35.0*pow(d,6.0)*r + 132.0*pow(d,5.0)*r*r + 310.0*pow(d,4.0)*pow(r,3.0) + 472.0*pow(d,3.0)*pow(r,4.0) + 456.0*d*d*pow(r,5.0) + 256.0*d*pow(r,6.0) + 64.0*pow(r,7.0);
  t5 = 5.0*pow(r,8.0) + 60.0*pow(r,6.0)*d0*d0 + 126.0*pow(r,4.0)*pow(d0,4.0) + 60.0*r*r*pow(d0,6.0) + 5.0*pow(d0,8.0);
  t6 = d0*d0-r*r;

  *f =  t1*pow(sig,12.0)*t3/45.0/pow(d,9.0)/pow(t2,9.0) - t1*pow(sig,6.0)/3.0/pow(d,3.0)/pow(t2,3.0);
  *df = -4.0*pow(sig,12.0)*t1*t4/(5*pow(d*t2,10.0)) + 2.0*pow(sig,6.0)*t1*d0/pow(d*t2,4.0);
  *dfr =   t1*pow(sig,12.0)*t5/(5*r*pow(t6,10.0)) - t1*pow(sig,6.0)*(r*r+d0*d0)/(r*pow(t6,4.0));
  }

  *d2f = 0.0; /*not used */
  *df /= d0;
}
