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


/* routines that perform output */

/*-------------- md_output.c ---------------------*/
void usage(char *);
void md_error(char *);
void md_warning(char *);
void md_stdout(char *);
void set_stderrnow(char *s);
void set_stdoutnow(char *s);

void writestart(FILENAMES *,int,SIMPARMS *,ENERGIES *,COORDS *);
void writescreen(SIMPARMS *,NGBR *,ENERGIES *,double *);
void finish(FILENAMES *,SIMPARMS *,ENERGIES *,COORDS *);
void write_coord(char *name,int natoms,COORDS *coords);
void write_all_coords(SIMPARMS *simparms,COORDS *coords);

void output_initval(SIMPARMS *,COORDS *,NGBR *,ENERGIES *);
void output_subparm(SIMPARMS *,WRITE_STEP *,SUBSPACE *,int);
void output_simparm(SIMPARMS *,WRITE_STEP *,int);

void openfiles(SIMPARMS *,WRITE_STEP *,FILENAMES *);
void write_vectors(int,double*,double*,double*,double*,
		   SUBSPACE*,double**,double*,double);
void save_conf(FILE *,int,COORDS *);
void save_vel(FILE *,int ,COORDS *);
void save_force(FILE *,int ,COORDS *);
void save_inst(FILE *,FILE *,FILE *,int,double,ENERGIES *,double*);
void save_hills(char *hillname,int istep,double dt,COLVAR *colvar);
void save_colv(FILE *fcolv,int istep,double dt,COLVAR *colvar);
void read_hills(char *hillfile,int istep,double dt,COLVAR *colvar);

void print_mat(int,int,double **);
void psysvec(FILENAMES *,int,int,int,int,double *,double *,double *,double *,
	     double *,double **,double*,double **,int,int *,int *);
void partratio(double **fcmat,double *dn,double *part,int natoms);
void save_nma(SIMPARMS *,FILENAMES *,int ,double ,NMODES *);
void write_nma_start(int,int,double,double,double,int,int,NMODES *);
void store_freq(SIMPARMS *,int nmodes,double *dn,double time);





void print_mat_new(int,int,double **);
void psysvec_new(FILENAMES *,int,int,int,int,double *,double *,double *,double *,
	     double *,double **,double*,double **,int,int *,int *);
void partratio_new(double **fcmat,double *dn,double *part,int natoms);
void save_nma_new(SIMPARMS *,FILENAMES *,int ,double ,NMODES *);
void write_nma_start_new(int,int,double,double,double,int,int,NMODES *);
void store_freq_new(SIMPARMS *,int nmodes,double *dn,double time);


void imag2real(int,double,double,double[],double[]);
void imag2real_new(int,double,double,double[],double[]);



