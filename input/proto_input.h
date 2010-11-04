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


/* prototype files in input directory */

void readpos(char *command,char *initfile,SIMPARMS *simparms,
	     int *estep,ENERGIES *energies,COORDS *coords);

void read_sim_inputs(SIMPARMS *simparms,SUBSPACE *subspace,
		     WRITE_STEP *write_step,FILENAMES *filenames,
		     INTER *inter,NGBR *ngbr,COORDS *coords,STARTUP *startup);

void set_sim_default(SIMPARMS *simparms,WRITE_STEP *write_step,
		     SUBSPACE *subspace,FILENAMES *filenames,
		     INTER *inter,NGBR *ngbr,COORDS *coords,STARTUP *startup);


void read_setfile(SIMPARMS *,WORD,int *,SPECIES *,COORDS *);

void read_parmfile(SIMPARMS *,COORDS *,STARTUP *,INTER *,NGBR *);

void read_top_file(SIMPARMS *,STARTUP *,ATOM_TOPOL **,
		   int **,int **, WORD **,
		   int **,int **,int **,int **, WORD **,
		   int **,int **,int **,int **, WORD **,
		   int **,int **,
		   int **,int **,int **, WORD **);

/* ----------------------- md_nmain.c ------------------------------*/
void open_nma_in(char *);
void read_header_conf(int,int *,double *);
int read_config(int ,double *,double *,double *,double *);
void read_cont(char [],int ,double *,double *,int ,NMODES *);


/* ----------------------- md_searchbase.c ------------------------------*/
void search_ter_base(SIMPARMS *,STARTUP *,INTER *,NGBR *);
void search_bond_base(SIMPARMS *,char *,int **,int **,WORD **,
		      ATOM_TOPOL**,SPECIES ,int ,BONDS *);
void search_bend_base(SIMPARMS *,char*,int**,int**,int**,int**,WORD**,
		      ATOM_TOPOL**,SPECIES,int,BENDS *);
void search_tors_base(SIMPARMS *,char*,int**,int**,int**,int**,WORD**,
		      ATOM_TOPOL**,SPECIES,int,TORS *);
void search_onfo_base(SIMPARMS *,char*,int,int**,int**,ATOM_TOPOL**,SPECIES,
		      int,double,double,int,ONFO *);
void search_bocross_base(SIMPARMS *,char*,int**,int**,int**,WORD**,
			 ATOM_TOPOL**,SPECIES,int,BONDXS *);
void search_bendquart_base(SIMPARMS *,char*,int**,int**,int**,WORD**,
		      ATOM_TOPOL**,SPECIES,int,BENDQUART *);

/* ------------------------- md_getparm.c ------------------------------*/
void get_ter_parm(int,META_KEY *,char *,WORD,WORD,
   double*,double*,double*,double*,double*,double*,INTER *,int ,int *,double);
void get_bondparm(int,META_KEY *,char *, WORD, WORD, WORD,int *,
		  double *,double *,double *,int *,double *);
void get_bondxparm(FILE *,char *,WORD,WORD,WORD,WORD,
		   double*,double*,double*,int *);
void get_bendparm(FILE *,char *,WORD,WORD,WORD,int *, WORD ,double *,double *,int *,
		  int,int,int,double *);
void get_torsparm(FILE *,char *,WORD,WORD,WORD,WORD,WORD,int *,
		  double *,double *,VPHI ,int *);
void get_onfo_parm(FILE *,char *,WORD,WORD,double *,double *,double *,
		   double *,double *,double *,int);

void get_bendquartparm(FILE *,char *,WORD,WORD,WORD,WORD,double*,double*,int *,
		       int,int,int);

/* ------------------------- read_freqfile ------------------------------*/

void read_freqfile(char *,int,double *);

/* ------------------------- read_extern_file ----------------------------*/

void read_extern_pot(char *filename,int ntypes,WORD types[],
		     int ntable,COORDS *coords);

void read_table(char *file,int ntable,double zmin,double zmax,
		double *vtab,double *dvtab,double *dv2tab);

/* prototypes for parse_mike */
void parse_file(char *filename,int *num_meta,META_KEY *meta_root);
void free_meta_keys(int nmeta,META_KEY **meta_key);

/* ------------------------- read colvar file ----------------------------*/

void colvar_init(int ncolvar,COLVAR *);
void read_colvarfile(SIMPARMS *simparms,char *colvarfile,COORDS *coords);
