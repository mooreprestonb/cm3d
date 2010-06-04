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


/* functional prototype for energy routines */

void force(SIMPARMS *,COORDS *,INTER *,NGBR *,int);
void getvireal(SIMPARMS *,COORDS*,ENERGIES *,INTER *,NGBR *,int);
void geteng(SIMPARMS *,COORDS *,ENERGIES *);
void force_extern(SIMPARMS *,COORDS *,int);
void getvextern(SIMPARMS *,COORDS *,double *,double *,double *,double *,int);
void freeze_atom(SIMPARMS *simparms,COORDS *coords);
void get_p2mt(int natoms,double [],double [],double [],double [],double [9]);

/*--------------------- md_finter.c -------------------------------*/
void set_exclude(SIMPARMS *simparms,COORDS *coords,INTER *inter);
void force_ter(SIMPARMS *,COORDS *,INTER *inter,NGBR *ngbr,int iflg);
void getvirter(SIMPARMS *,COORDS *,INTER *,NGBR *,double*,double*,double*,
                double [9]);
void fcinter(SIMPARMS *,COORDS *,INTER *,double **);
void fcintern(SIMPARMS *,COORDS *,INTER *,double **);
void get_ngbr(SIMPARMS *,COORDS *,INTER *,NGBR *);
void check_distance(SIMPARMS *,COORDS *,INTER *);
void getnpair(NGBR *);
void nrefsh(SIMPARMS *,COORDS *,INTER *,NGBR *);
double pot_inter(double *px,double *py,double *pz);
void allocate_lists(int natoms,int nlen,NGBR *ngbr,COORDS *coords,int *mem);
void free_lists(NGBR *ngbr,COORDS *coords);

/*-------------------------md_fnopol.c--------------------------*/
void force_npol(SIMPARMS *,COORDS *,int ,int *,int *,int ,int *,INTER *,int);
void vinter_npol(int,int *,int *,int *,double *,double *,double *,double *,
		 int ,int **,double *,double *,double *,int ,
		 double **,double **,double *,double *,double *,
                 double [9],int);
void getsvirter(double *,double *,double *,double [9]);
void zero_potl(void);

/*--------------------- md_ewald.c -------------------------------*/
void check_distance_e(SIMPARMS *simparms,COORDS *coords,INTER *inter);
void force_ter_e(SIMPARMS *,COORDS *,INTER *,NGBR *,int iflg);
void getvirter_e(SIMPARMS *,COORDS *,INTER *,NGBR *,double*,double*,
		 double[],int);
void ewald_setup(SIMPARMS *,COUL *coul,double *qch,int kmax,double alpha);
void k_ewald(COORDS *coords,double *vkspace);
void fk_ewald(SIMPARMS *,COORDS *);
void getvewald(SIMPARMS *,COORDS *,double *,double *,double [9]);
void ecorr(SIMPARMS *,COORDS *);
void r_ewald(SIMPARMS *simparms,COORDS *coords,double *vrspace);
void fr_ewald(SIMPARMS *simparms,COORDS *coords);

void allocate_lists_e(int ncharge,NGBR *ngbr,int *mem);

/*---------------------finter_e.c------------------------------------*/
void zero_pote(void);
void force_e(int ,int ,int *,int *,COORDS *,INTER *,int ,int ,int,int );
void vinter_e(int ,int *,int *,COORDS *,INTER *,int iperd,
	      double *spelec,double *spelech,double *swelec,
              double swelectensor[9],int,int,int);
void getsvirter_e(double *velec,double *welec,double welectensor[9]);
void force_dipole(int , double *,double *,double *,double *,double *,double *,
		  double *,double *);

double pot_dipole(int nlen,double *dx,double *dy,double *dz,double *qch,
		  double *alpha);

/*---------------------elec_ngbr.c------------------------------------*/
void count_charge(SIMPARMS *,COORDS *,INTER *,NGBR *,int *,int *);
void set_charges(SIMPARMS *,COORDS *,INTER *,NGBR *,STARTUP *);
void set_ecorr(SIMPARMS *,COORDS *,int *,int **);
void get_ngbr_e(SIMPARMS *simparms,COORDS *coords,INTER *inter,NGBR *ngbr);
void ngbr_liste(int,int *,int *,int,COORDS *,INTER *,NGBR *,int);
void cnt_ngbr_liste(int,int *,int *,int,COORDS *,INTER *,int *,int *,int);

/* ----------------------set_link.c--------------------------*/
void get_lnklist(int natoms,int ncell,int *head,int *llist,int iperd,
		 double *hmat,double *x,double *y,double *z);
void set_link(SIMPARMS *,NGBR *,double *,double,double);
void get_lnklist_e(int natoms,int ncell,int *head,int *llist,int iperd,
		   double *hmat,double *x,double *y,double *z,double *q);
void count_link_pair(int ncell,int iperd,double hmat[],int ibegin,int iend,
		     int *nvdw,int *nelec,double rmax,double rmaxe);

/*--------------------- intra molecular energies --------------------------*/
void exclude_bonds(int *,int **,BONDS *);
void fbond(SIMPARMS *,COORDS *);
void getvbond(SIMPARMS *,COORDS *,double *,double *,double [9]);
void fcbond(COORDS *,double **);
void fcbondn(COORDS *,double **);
double potbond(double *,double *,double *);

void exclude_bends(int *,int **,BENDS *);
void fbend(SIMPARMS *,COORDS *);
void getvbend(SIMPARMS *,COORDS *,double *,double *,double [9]);
void fcbend(COORDS *,double **);
void fcbendn(COORDS *,double **);
double potbend(double *,double *,double *);

void exclude_tors(int *,int **,TORS *);
void ftors(SIMPARMS *,COORDS *);
void getvtors(SIMPARMS *,COORDS *,double *,double *,double [9]);
void fctors(COORDS *,double **);
void fctorsn(COORDS *,double **);
double pottors(double *,double *,double *);

void exclude_onfo(int *,int **,ONFO *);
void fonfo(SIMPARMS *,COORDS *);
void getvonfo(SIMPARMS *,COORDS *,double*,double*,double*,double*,
               double [9],double [9]);
void fconefour(COORDS *,double **);
void fconefourn(COORDS *,double **);
double pot_onefour(double *,double *,double *);

void fbondx(SIMPARMS *,COORDS *);
void getvbondx(SIMPARMS *,COORDS *,double *,double *,double [9]);
void fcbondx(COORDS *,double **);
void fcbondxn(COORDS *,double **);
double potbondx(double *,double *,double *);

/*---------------------- md_funcvdw.c -----------------------------*/
void func_lj(double,double *,double *,double *,double,double);
void func_will(double,double*,double*,double*,
	       double,double,double,double,double);
void func_aziz(double ,double *,double *,double *,double ,double ,
	       double ,double ,double ,double ,double ,double ,double);
void func_hydbnd(double ,double *,double *,double *,double ,double );
void func_kerf(double ,double *,double *,double *,double );
void func_coul(double ,double *,double *,double *);
void func_hautman(double x,double *f,double *df,double *d2f,
		  double zmin,double c12,double c3);
void func_lj64(double,double *,double *,double *,double,double);
void func_lj96(double,double *,double *,double *,double,double);
void func_lj124(double,double *,double *,double *,double,double);
void func_tab_pot(FILE *, double ,double *,double *,double *);
void func_dzug(double ,double *,double *,double *,double,double);

/*------------------------ md_subspace.c ---------------------------*/
void diagonalize(double **,int,double*,double,double,int*,int*,int);
void find_zeros(int ,double *,int *,int *,int *);
void pack(SUBSPACE *,int, int,double *,double **);
void pack_d2v(int,SUBSPACE*,int,double**,double*,double**,double*,int*,int);

void get_subspace(SIMPARMS *simparms,COORDS *coords,INTER *,
		  double *x2,double *y2,double *z2,SUBSPACE *subspace);
void get_subvec(SIMPARMS *,COORDS *,INTER *,SUBSPACE *);
void get_grp_fcmat(int,int*,int,int*,double**,SUBSPACE*,int*,double***);

/*---------------------- md_nmafcntr.c -----------------------------*/
void fmat(SIMPARMS *simparms,COORDS *coords,INTER *inter,NMODES *nmodes);
void copy(int ,int ,double **,int ,int,double **,int ,int );
void get_axis(int,int,int*,int*,COORDS *,double**);
void overper(int ,int,int*,int*,double **,double **,double **);

/*----------------------- friciton.c ----------------------------------*/
void friction(SIMPARMS *simparms,COORDS *coords,INTER *inter,NMODES *nmodes);

/*----------------------- md_spectrum ---------------------------------*/
void avg_spec(int natoms,NMODES *nmodes);
void get_spectrum(int ,double *,double *,double **,double *,double **);
void get_specq(int ,double *,double *,double **,double *);
void get_spect_full(int ,double *,double *,double **,double *,double **);

/* ------------------------hills---------------------------------------*/

void fcolvar_cntr(SIMPARMS *,COORDS *,int force_flag);
void hills(SIMPARMS *simparms,COORDS *coords);
void add_hills(COLVAR *colvar);
