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

typedef char WORD[MAXWORDLEN];
typedef char LINE[MAXLINELEN];
typedef double VPHI[MAXPOWER_TORS];

/* ----------------system.h defs ------------------*/

typedef struct KEYS {
  char *name,*arg;
  int itype;
  struct KEYS *next;
} KEYS;

typedef struct META_KEY {
  char *meta;
  int nkeys;
  KEYS *keys;
  struct META_KEY *next;
} META_KEY;

/* ----------------system.h defs ------------------*/

typedef struct keys {
  WORD keyword,keyarg;
  struct keys *next; 
} KEY;

typedef struct species {
  WORD filename,thermopt,name,group;
  int nmol,molidx,napm,nbond,nbend,ntors,n14,npath,nbondx,nbendquart;
  struct species *next; 
} SPECIES;

typedef struct coul {
  int ncharge,*icharge;
  int *nexclude,**exclude;
  int *ka,*kb,*kc,*kup,ktot;
  int necorr,*iecorr,*jecorr;
  double kappa,kappa2,sqrt_pi,ewald_self,ewald_bkg,falp2;
  double vk_ewald,vp_ewald,pot_ecorr,prs_ecorr;
  double prs_ecorr_tensor[9];
  double vp_ewald_tensor[9];
  double *cossc,*sinsc,*helr,*heli;
} COUL;

typedef struct bonds {
  int itypes,nbonds,*itypbond,*idxbond,*ibond,*jbond,*ifix;
  double *eqbond,*fkbond,*bmorse,*aq;
} BONDS;

typedef struct bends {
  int itypes,nbends,*itypbend,*idxbend,*ibend,*jbend,*kbend,*ifix;
  double *eqbend,*fkbend,*aq;
} BENDS;

typedef struct tors {
  int itypes,ntorss,*idxtors,*itors,*jtors,*ktors,*ltors,*ifix;
  int *ind_pot,*ind_pot_num;
  double *eqtors,*fktors;
  VPHI *vphi;
} TORS;

typedef struct onfo {
  int itypes,n14s,*i14,*j14,*i14dx,n14table;
  double **v14tab,**dv14tab,**d2v14tab;
  double *r14min2,*r14max2,*dx214tab;
  double scale_onefour,scale_onefour_e;
} ONFO;

typedef struct bondxs {
  int itypes,nbondxs,*idxbox,*ibondx,*jbondx,*kbondx,*ifix;
  double *eq1bondx,*eq2bondx,*fkbondx;
} BONDXS;

/* added June 2000, JCS */
typedef struct bendquart {

  int itypes,nbendquarts,*idxbendquart,*ibendquart,*jbendquart,*kbendquart;
  int *ifix;
  double *eqbendquart,*fkbendquart;
} BENDQUART;
/* end of added June 2000, JCS */

typedef struct subspace {
  int nstate,num_subst,stemp_update;
  int *num_vecsub,*num_real,*num_imag;
  int *num_vecsub_max,*num_real_max,*num_imag_max;
  int igetvec,nstate_max;
  double freq_min,freq_max;
  WORD vecfile,vecconf;
  double *q0,*q1,*q2;
  double **d2v,*d2vn;
} SUBSPACE;

typedef struct simparms {
  int icalc_type,istart,nma_type,fstep,istep,nstep;
  int natoms,nspec,ntypes,iensemble;
  int ninter,ninnra,ninrra,nbar,ivol,iperd,imin_type,nocell_all;
  int iunits,ndof,nlen,npoints,iextern,nhc;
  int ntherm,nchain,max_exclude,nfreeze;
  int rank,size,mem_bytes,nyoshida,ipolar;
  int *itype,*iatom,*imolecule,*ispecies,*igroup;
  double temp,dt,pext,vlrc,wlrc,scaleeps,scalecharge,scaletemp;
  double min_tol,wallclock,wyosh[3];
  WORD *atom,*mole,*group;
} SIMPARMS;

typedef struct sfilename {
  FILE *fconf,*fvel,*fham,*feng,*fext,*fforce,*fcolv;
  WORD restart,initfile,configs,velfile,instham,insteng,setfile,forcefile;
  WORD nmfile,specfile,molvec,sysvec,extfile,instext,partrat,colvarcntrfile;
  WORD command,infile,outfile,errfile,colvfile,hillfile;
} FILENAMES;

typedef struct senergies {
  int estep;
  double avham,avpotra,avpoter,avzke,avpotn,avzken,avpotv,avzkev;
  double avtemp,avprs,avvol,acpu;
  double avprs_tensor[9],prs_tensor[9];
  double prsi,tiout,ham;
  double poth,potra,poter,zke,zken,potn,zkev,potv;
  double pot_inter,pot_bond,pot_bend,pot_bendquart,pot_tors,pot_extern,pot_extern_e;
  double pot_onfo,pot_onfo_e,pot_elec,pot_recip;
  double avpot_inter,avpot_bond,avpot_bend,avpot_bendquart,avpot_tors;
  double avpot_onfo,avpot_onfo_e,avpot_elec,avpot_recip;
  double avpot_extern,avpot_extern_e;
  double econv,econ2;
  double econi,deltae2;
  double vinter,vhinter,velec,vewald,vbond,vbend,vbendquart,vtors,vonfo,vonfo_e;
  double vextern,vextern_e,vbondx;
  double winter,welec,wewald,wbond,wbend,wbendquart,wtors,wonfo,wonfo_e,wbondx;
  double wintertensor[9],welectensor[9],wbondtensor[9],wbondxtensor[9];
  double wbendtensor[9],wbendquarttensor[9],wonfotensor[9],wonfo_etensor[9];
  double wtorstensor[9],wewaldtensor[9],WItensor[9],Wtensor[9];
  double wextern,wextern_e;
  double US,U,UI,W,WI;
  double colvar;
} ENERGIES;

typedef struct satom_topol {
  int atm_idx;
  WORD type,group;
  double mass,charge,alpha;
} ATOM_TOPOL;

typedef struct swrite_step {
  int nscrn,ndump,nvel,nconf,ninst,nupdss,nforce,ncolv;
  int resamvel,psysvec,rescalevel,peigval;
} WRITE_STEP;

typedef struct colvar {
  int ncolvar,lrestart,nhills,ltunehills;
  int maxstephill,minstephill,lasthill;/* max and min of step before we place a hill */
  int *lconst_vel,*lhills,*type;
  int *ndim,**natom,***atom; /* ndim[ncolvar] # of centers in your constraint */
  /* natom[ncolvar,ndim] # of atoms in each center */
  /* atom[ncolvar,ndim,natom] array of atom indexes */
  double ecolvar;
  double *amass,*fk;
  double *pcolvar,*vcolvar,*fcolvar,*pistcolvar;
  double *pcolvarlast,*pcolvarmin,*pcolvarmax,*vcolvarlast,*pcolvarmean,*pcolvarstd;
  double *f_hills,*f_spring,*f_spring_cum;
  double *hillwidth, hilldepth;
  double **t_pcolvar;  /* t_pcolvar[ncolvar,nhills] */
  double **t_hillwidth, *t_hilldepth;
  /* t_hillwidth[ncolvar,nhills] t_hilldepth[nhills] */
  double *min_val, *max_val;   /* [ncolvar] coll.var. boundary wall parameters */
  double *mindwidth,mindtol2; /* min distance relative to hillwidth before we place a hill */
  int *istep;
  /* Steve adds */
  double meanforce_inst, meanforce_cum;
  double *temp, *thrmmass, *f_boundry;
} COLVAR;


typedef struct scoords {
  int nbreak,*ibreak; /* nlength scratch */
  int *ithm,*ifreeze;
  int *indx,*jndx,*ityp,*jtyp,*intr;
  double *amass, *qch, *alpha;
  double cmx, cmy, cmz;

/* XXX new ones */
	double *p, *v;		/* position and velocity */
	double *fs, *fl;
	double *fr, *fa, *ft;
/* end new ones */

/* TBC */
  double *px, *py, *pz;
  double *vx, *vy, *vz;
  double *fxs,*fys,*fzs,*fxl,*fyl,*fzl,*fxr,*fyr,*fzr;
  double *fxa,*fya,*fza,*fxt,*fyt,*fzt;
/* end TBC */

  double *ux,*uy,*uz;
  double *xold,*yold,*zold;
  double *xdis,*ydis,*zdis;  /* nlength scratch */
  double **eta,**veta,**feta,**mnh,**gkt,*p2mt;
  double *bar,*vbar,*mbar;
  double zmin_extern,zmax_extern,dxtab_extern;
  double **dvtab_extern,**vtab_extern,*dvtab_extern_e,*vtab_extern_e;
  double pvol,vvol,fvol,mvol;
  double vvol9[9],fvol9[9];
  double hmat[9],hmati[9];
  double cell_old;
  double *scr_buf,*scr_rec;
  double *ffreeze;
  COUL coul;
  BONDS bonds;
  BENDS bends;
  TORS tors;
  ONFO onfo;
  BONDXS bondxs;
  BENDQUART bendquarts;
  COLVAR colvar;
} COORDS;

typedef struct sngbr {
  int update,ilist;
  int nprs,nprt,nprse,nprte;
  int nprst,nprtt,nprset,nprtet;
  int *inps,*jnps,*inpt,*jnpt,*offs,*offl,*offt;
  int *inpse,*jnpse,*inpte,*jnpte,*offse,*offte;

  int ncells,*lnklist,*hlist;
  int npaircell,*icell,*jcell,npaircelle,*icelle,*jcelle;
} NGBR;

typedef struct sinter {
  int ntable,ntype;
  int *itype,**map,*nexclude,**exclude;
  double **vtab,**dvtab,**d2vtab,**switchs;
  double *rminl_cut2,*rmaxl_cut2,*dx2tab;
  double *rmins_cut2,*rmaxs_cut2;
  double *rmaxs_skin2,*rmaxl_skin2,*rminl_skin2;
  double *dxn,*dyn,*dzn;
  double skin_ter,rheal;
  /* electrostatic */
  double *coultab,*dcoultab,*switchse,dx2coul;
  double rmins_coul2,rmaxs_coul2,rminl_coul2,rmaxl_coul2;
  double rmaxs_skin2_e,rmaxl_skin2_e,rminl_skin2_e;
} INTER;

typedef struct sstartup {
  WORD inter_file,bond_file,bend_file,tors_file,onfo_file;
  SPECIES spec_root;
  int i,nspec,kmax,ishift,ivlrc;
  double tau_nhc,tau_vol,alp_ewald;
  double scale_onfo,scale_onfo_e;
  double rcute_max,rcute_min,rcute_resp;
  LINE nvecstore,nvecreal,nvecimag;
} STARTUP;

typedef struct nmodes {
  int nconf,nconf_off,flo,ier,imodes;
  int nspec,*napm,*nmol;
  double maxfreq,minfreq,acpu,df;
  double **fcmat,*dn,*en,*en2,*freq,*fric,*fpos,*ftot;
  double **molmat,**decmod,**permod,*fricmod,*fricmodq,*fricmodcq;
  double **decmodd, **decmodq, **decmoddq,**decmodrr,**decmod_isorr;
  double *aniso_raman,*iso_raman,*raman,*redram,*oke,*iso;
  double *suspx,*suspy,*suspz;
  double *avg_freq,*avg_strn;
  double *part,*fpart;
  double *A_xx1,*A_zz1,*A_xz1,*A_zx1;
  double *A_xx2,*A_zz2,*A_xz2,*A_zx2;
  double *mux1,*muy1,*muz1,*mux2,*muy2,*muz2;
} NMODES;
