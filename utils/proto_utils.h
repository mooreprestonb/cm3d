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

/* function prototypes for mathlib directory */

double gethinv9(double hmat[], double hmati[]);
double get_deth(double []);
double ddot1(int n,double *a,double *b);
double dsum1(int n,double *a);
double ddot1a(int n,double *a,int astep,double *b,int bstep);
double dsum1a(int n,double *a,int astep);
int locate(int *,int,int);
int insert(int *ia,int *n,int jpoint);
void splint(double [], double [], double [], int , double , double *);
void spline(double [], double [], int , double , double , double []);
void polint(double xa[], double ya[], int n, double x, double *y, double *dy);

double *dvector(int,int);
double **dmatrix(int,int,int,int);
double **drematrix(double **,int,int,int,int,int,int,int,int);
double ***d3tensor(int nrl,int nrh,int ncl,int nch,int ndl,int ndh);
double cputime(void);
double realtime(void);
void free_dvector(double *,int,int);
void free_dmatrix(double **,int,int,int,int);
void free_d3tensor(double ***,int,int,int,int,int,int);

void srandme(int);
void ggauss(int ,double *);
void *cmalloc(size_t);
void *crealloc(void *,size_t);
void decomp1d(int,int,int,long *,long *);
void decomp_trig(int,int,int,long *,long *);

void rs_me(int n,double *,double **,int);
void ctred2(double **,int,double [],double []);
void ctred2v(double **,int,double [],double []);
void ctqli(double [],double [],int,double **);
void ctqliv(double [],double [],int,double **);
void ceigsrt(double [],double **,int);
void ceigsrtv(double [],double **,int);
void invres(int,double **,double **);
void vsrtasnd(double d[],double **v,int m,int n);

/*-----------------------matopt.c-------------------------------*/
void matmul(int n,double **a,double **b,double **c);
void vecmat(int n,double **a,double *x,double *y);

void matmulv(int n,double a[],double b[],double c[]);
void transmatv(int n,double a[],double b[]);
void sym33(double mat[9]);
void test_anint(void);

/*-------------------------util.c------------------------------*/
int same_mole(int iatom,int jatom,int *itype_species,int *itype_molecules);
FILE *cfopenw(char *name);
void *cmalloc(size_t mem);
void allocate_coords(SIMPARMS *simparms,COORDS *coords);
void allocate_extend(SIMPARMS *simparms,COORDS *coords);
void allocate_exttab(SIMPARMS *simparms,COORDS *coords,int);
void rmass(int natoms,double *fx,double *fy,double *fz,double *amass);
void free_exttab(COORDS *coords);
void free_coords(COORDS *coords);
void free_extend(SIMPARMS *,COORDS *coords);
void free_all_pointers(SIMPARMS *simparms,COORDS *coords,INTER *,NGBR *);

/*-------------------------checkvals.c------------------------------*/
int check_vals(int,...);
/*-------------------------checkvals.c------------------------------*/

/*------------------------md_char.c----------------------------*/
int get_dict_num(int,WORD [],char *);
void set_sim_keyword(int *,WORD **dict);
void get_meta_key(FILE *,char *,char *,int *,int *,char *,KEY *,int *);
void get_keys(FILE *,char *,char *,int *,int *,int *,KEY *);
void free_key(int,KEY *);
void syntax_error(const char *,FILE *,int ,char *,int);
void print_dict(WORD,int,WORD []);
void err_found(WORD,int,char *,int);
void write_key(KEY *);
void get_sim_keys(const char *,const char *,int *,KEY *);
int getnxtint(char *);
/*------------------------md_char.c----------------------------*/

void v1fv3(int ,double *, double *, double *, double *);
void v3fv1(int ,double *, double *, double *, double *);






double gethinv9_new(double hmat[], double hmati[]);
double get_deth_new(double []);
double ddot1_new(int n,double *a,double *b);
double dsum1_new(int n,double *a);
double ddot1a_new(int n,double *a,int astep,double *b,int bstep);
double dsum1a_new(int n,double *a,int astep);
int locate_new(int *,int,int);
int insert_new(int *ia,int *n,int jpoint);
void splint_new(double [], double [], double [], int , double , double *);
void spline_new(double [], double [], int , double , double , double []);
void polint_new(double xa[], double ya[], int n, double x, double *y, double *dy);

void ggauss_new(int ,double *);
void *cmalloc_new(size_t);
void *crealloc_new(void *,size_t);
void decomp1d_new(int,int,int,long *,long *);
void decomp_trig_new(int,int,int,long *,long *);

void rs_me_new(int n,double *,double **,int);
void ctred2_new(double **,int,double [],double []);
void ctred2v_new(double **,int,double [],double []);
void ctqli_new(double [],double [],int,double **);
void ctqliv_new(double [],double [],int,double **);
void ceigsrt_new(double [],double **,int);
void ceigsrtv_new(double [],double **,int);
void invres_new(int,double **,double **);
void vsrtasnd_new(double d[],double **v,int m,int n);


/*-------------------------util.c------------------------------*/
int same_mole_new(int iatom,int jatom,int *itype_species,int *itype_molecules);
FILE *cfopenw_new(char *name);
void *cmalloc_new(size_t mem);
void allocate_coords_new(SIMPARMS *simparms,COORDS *coords);
void allocate_extend_new(SIMPARMS *simparms,COORDS *coords);
void allocate_exttab_new(SIMPARMS *simparms,COORDS *coords,int);
void rmass_new(int natoms,double *fx,double *fy,double *fz,double *amass);
void free_exttab_new(COORDS *coords);
void free_coords_new(COORDS *coords);
void free_extend_new(SIMPARMS *,COORDS *coords);
void free_all_pointers_new(SIMPARMS *simparms,COORDS *coords,INTER *,NGBR *);

