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

/* prototypes for intergrate routines */

/* -----------------------------md_averages----------------------------*/
void set_econi(ENERGIES *energies);
void accumulate(double *,ENERGIES *);
void zero_averages(SIMPARMS *simparms,ENERGIES *energies);

void substep_init(SIMPARMS*,COORDS*,ENERGIES *,SUBSPACE *,INTER *,NGBR *);
void update_subspace(SIMPARMS *,COORDS *,INTER *,SUBSPACE *,NGBR *);
void sub_step(SIMPARMS *,COORDS *,ENERGIES *,INTER *,NGBR *,SUBSPACE *);

void step_init(SIMPARMS *,COORDS *,ENERGIES *,INTER *,NGBR *);
void step(SIMPARMS *,COORDS *,ENERGIES *,INTER *,NGBR *);

void set_baro(SIMPARMS *simparms,COORDS *coords,double tauv);
void forceV(SIMPARMS *simparms,COORDS *coords,ENERGIES *energies);
void samvbar(SIMPARMS *simparms,COORDS *coords);
void getengv(SIMPARMS *simparms,COORDS *coords,double *zkev,double *potv);
void getengv_cheese(SIMPARMS *simparms,COORDS *coords,double *zkev,
		    double *potv,double *zken);
void intbar(SIMPARMS *simparms,COORDS *coords,double dta2,int iflag);

void set_therm(SIMPARMS *simparms,COORDS *coords,
	       int nspec,SPECIES spec_root,double tau_nhc);
void samveta(SIMPARMS *,COORDS *);
void forcen(SIMPARMS *,COORDS *);
void getengn(SIMPARMS *,COORDS *,double *,double *);
void getengn_cheese(SIMPARMS *,COORDS *,double *,double *);
void intnhc(SIMPARMS *,COORDS *,double);
void period(int natoms,double dxs[],double dys[],double dzs[],
	    double hmat[],double hmati[],int iperd,int iensemble);

/* -----------------------md_zero.c---------------------------*/
void samvel(SIMPARMS *,COORDS *,int);
void zerotm(int ,double *,double *,double *,double *);
void rotate(int natoms,double *px,double *py,double *pz,double *amass);
void scale(int ,double *,double *,double *,double *,double,int);
void radial(int ,double *,double *,double *,
	    double *,double *,double *,double *);
void zerocm(int,double *,double *,double *,double *);
void zero_angm(int,double*,double*,double*,double*,double*,double*,double*);
double get_temp(int ,double *,double *,double *,double *,int );


/*----------------------------------freqbin.c----------------------------*/
void freqbin(SIMPARMS *simparms,int npoints,double freq_min,NMODES *nmodes);

/* --------------------------------------baro.c--------------------------*/
void app_ux2(int natoms,double dta2,COORDS *coords);
void app_up2(int natoms,double dta2,double dnf,COORDS *coords);
void app_uh(double dta2,COORDS *coords);









/* -----------------------------md_averages----------------------------*/
void set_econi_new(ENERGIES *energies);
void accumulate_new(double *,ENERGIES *);
void zero_averages_new(SIMPARMS *simparms,ENERGIES *energies);

void substep_init_new(SIMPARMS*,COORDS*,ENERGIES *,SUBSPACE *,INTER *,NGBR *);
void update_subspace_new(SIMPARMS *,COORDS *,INTER *,SUBSPACE *,NGBR *);
void sub_step_new(SIMPARMS *,COORDS *,ENERGIES *,INTER *,NGBR *,SUBSPACE *);

void step_init_new(SIMPARMS *,COORDS *,ENERGIES *,INTER *,NGBR *);
void step_new(SIMPARMS *,COORDS *,ENERGIES *,INTER *,NGBR *);

void set_baro_new(SIMPARMS *simparms,COORDS *coords,double tauv);
void forceV_new(SIMPARMS *simparms,COORDS *coords,ENERGIES *energies);
void samvbar_new(SIMPARMS *simparms,COORDS *coords);
void getengv_new(SIMPARMS *simparms,COORDS *coords,double *zkev,double *potv);
void getengv_cheese_new(SIMPARMS *simparms,COORDS *coords,double *zkev, double *potv,double *zken);
void intbar_new(SIMPARMS *simparms,COORDS *coords,double dta2,int iflag);

void set_therm_new(SIMPARMS *simparms,COORDS *coords, int nspec,SPECIES spec_root,double tau_nhc);
void samveta_new(SIMPARMS *,COORDS *);
void forcen_new(SIMPARMS *,COORDS *);
void getengn_new(SIMPARMS *,COORDS *,double *,double *);
void getengn_cheese_new(SIMPARMS *,COORDS *,double *,double *);
void intnhc_new(SIMPARMS *,COORDS *,double);
void period_new(int natoms,double dxs[],double dys[],double dzs[], double hmat[],double hmati[],int iperd,int iensemble);

/* -----------------------md_zero.c---------------------------*/

void samvel_new(SIMPARMS *,COORDS *,int);
void zerotm_new(int ,double *,double *);
void rotate_new(int ,double *,double *);
void scale_new(int ,double *,double *,double,int);
void radial_new(int ,double *, double *, double *);
void zerocm_new(int,double *,double *);
void zero_angm_new(int,double*,double*,double*);
double get_temp_new(int ,double *,double *,int );

/*----------------------------------freqbin.c----------------------------*/
void freqbin_new(SIMPARMS *simparms,int npoints,double freq_min,NMODES *nmodes);

/* --------------------------------------baro.c--------------------------*/
void app_ux2_new(int natoms,double dta2,COORDS *coords);
void app_up2_new(int natoms,double dta2,double dnf,COORDS *coords);
void app_uh_new(double dta2,COORDS *coords);


