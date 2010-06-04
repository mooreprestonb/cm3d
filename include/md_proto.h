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

/* protofiles for all routines in md */

#include "../utils/proto_utils.h"
#include "../input/proto_input.h"
#include "../output/proto_output.h"
#include "../integrate/proto_integrate.h"
#include "../main/proto_main.h"
#include "../energy/proto_energy.h"
#include "../minimize/minimize.h"
#include "../energy/polar/proto_polar.h"

void get_species(FILE *,SPECIES *);
void get_atom(FILE *,char *,int ,ATOM_TOPOL *);
void get_bond(FILE *,char *,int,int *,int *,WORD *);
void get_bend(FILE *,char *,int,int *,int *,int *,int *,WORD *);
void get_bendquart(FILE *,char *,int,int *,int *,int *,WORD *);
void get_tors(FILE *,char *,int,int *,int *,int *,int *,WORD *tto);
void get_onefour(FILE *,char *,int,int *,int *);
void get_bondx(FILE *,char *,int,int *,int *,int *,WORD *);

void keep(char *);

void save_eigvec(char *,int ,int ,double *,double **);
void read_eigvec(char *,int ,int ,double *,double **);

void ngradv2(double *,double *,double *,int,int, double (*)(double *,double *,double *),double **);

void get_fcmat(SIMPARMS *simparms,COORDS *coords,INTER *,double **fcmat);
void massweight_fcmat(int,double **,double *);
void reflect_fcmat(int ,double **);
void azero(double [],int);
void expandztox(double *,double *,double *,double *,double **,int ,int);
void projectxtoz(double *,double *,double *,double *,double **,int,int);
void mass_weight(int ,double *,double *,double *,double *,int);

int exclij(int i,int j,int *nexcl,int **excl);

void user_initial(int,int,double,double*,double*,double*,
		  double*,double*,double*,double*,double*,double*);
void user_function(int,int,double,double*,double*,double*,
		   double*,double*,double*,double*,double*,double*);

void mal_verify(int);
