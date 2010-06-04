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

typedef struct keys {
  WORD keyword,keyarg;
  struct keys *next; 
} KEY;

typedef struct species {
  WORD filename,thermopt,name;
  int nmol,molidx,napm,nbond,nbend,ntors,n14,nbondx;
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
  int itypes,nbonds,*idxbond,*ibond,*jbond,*ifix;
  double *eqbond,*fkbond;
} BONDS;

typedef struct bends {
  int itypes,nbends,*idxbend,*ibend,*jbend,*kbend,*ifix;
  double *eqbend,*fkbend;
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

typedef struct atom_map {
  int type,atom,molecule,species,group;
} ATOM_MAP;

typedef struct simparms {
  int icalc_type,istart,fstep,istep,nstep,natoms,nspec,iensemble;
  int ninner,ninnra,ninrra,nbar,ivol,iperd,imin_type,nocell_all;
  int iunits,ndof,nlen,npoints,iextern,nhc;
  int ntherm,nchain,max_exclude,nfreeze;
  int rank,size,mem_bytes,nyoshida;
  double temp,dt,pext,vlrc,wlrc,scaleeps,scalecharge;
  double min_tol,wallclock,wyosh[3];
  ATOM_MAP *atom_map;
  WORD *atom,*mole,*group;
  LINE yoshidas;
} SIMPARMS;

typedef struct senergies {
  int estep;
  double avham,avpotra,avpoter,avzke,avpotn,avzken,avpotv,avzkev;
  double avtemp,avprs,avvol,acpu;
  double avprs_tensor[9],prs_tensor[9];
  double prsi,tiout,ham;
  double poth,potra,poter,zke,zken,potn,zkev,potv;
  double pot_inter,pot_bond,pot_bend,pot_tors,pot_extern,pot_extern_e;
  double pot_onfo,pot_onfo_e,pot_elec,pot_recip;
  double avpot_inter,avpot_bond,avpot_bend,avpot_tors;
  double avpot_onfo,avpot_onfo_e,avpot_elec,avpot_recip;
  double avpot_extern,avpot_extern_e;
  double econv,econ2;
  double econi,deltae2;
  double vinter,vhinter,velec,vewald,vbond,vbend,vtors,vonfo,vonfo_e;
  double vextern,vextern_e,vbondx;
  double winter,welec,wewald,wbond,wbend,wtors,wonfo,wonfo_e,wbondx;
  double wintertensor[9],welectensor[9],wbondtensor[9],wbondxtensor[9];
  double wbendtensor[9],wonfotensor[9],wonfo_etensor[9];
  double wtorstensor[9],wewaldtensor[9],WItensor[9],Wtensor[9];
  double wextern,wextern_e;
  double US,U,UI,W,WI;
} ENERGIES;

typedef struct sfilename {
  WORD restart,initfile,configs,velfile,instham,insteng,setfile,forcefile;
  WORD nmfile,specfile,molvec,sysvec,command,infile,extfile,instext;
  FILE *fconf,*fvel,*fham,*feng,*fext,*fforce;
} FILENAMES;

typedef struct satom_topol {
  int atm_idx;
  WORD type,group;
  double mass,charge,alpha;
} ATOM_TOPOL;

typedef struct swrite_step {
  int nscrn,ndump,nvel,nconf,ninst,nupdss,nforce;
  int resamvel,psysvec,rescalevel;
} WRITE_STEP;

typedef struct scoords {
  int nbreak,*ibreak; /* nlength scratch */
  int *ithm,*ifreeze;
  int *indx,*jndx,*ityp,*jtyp,*intr;
  double *px,*py,*pz,*vx,*vy,*vz,*amass,*qch,*alpha;
  double *ux,*uy,*uz;
  double **eta,**veta,**feta,**mnh,**gkt,*p2mt;
  double *bar,*vbar,*mbar;
  double *fxs,*fys,*fzs,*fxl,*fyl,*fzl,*fxr,*fyr,*fzr;
  double *fxa,*fya,*fza,*fxt,*fyt,*fzt;
  double rmin_extern,rmax_extern,dxtab_extern;
  double *dvtab_extern,*vtab_extern,*dvtab_extern_e,*vtab_extern_e;
  double pvol,vvol,fvol,mvol;
  double vvol9[9],fvol9[9];
  double hmat[9],hmati[9];
  double *scr_buf,*scr_rec;
  double *xold,*yold,*zold,cell_old;
  double *xdis,*ydis,*zdis;  /* nlength scratch */
  COUL coul;
  BONDS bonds;
  BENDS bends;
  TORS tors;
  ONFO onfo;
  BONDXS bondxs;
} COORDS;

typedef struct sngbr {
  int update,ilist;
  int nprs,nprl,nprt,nprse,nprle,nprte;
  int nprst,nprlt,nprtt,nprset,nprlet,nprtet;
  int *inps,*jnps,*inpl,*jnpl,*inpt,*jnpt,*offs,*offl,*offt;
  int *inpse,*jnpse,*inple,*jnple,*inpte,*jnpte,*offse,*offle,*offte;

  int ncells,*lnklist,*hlist;
  int npaircell,*icell,*jcell,npaircelle,*icelle,*jcelle;
} NGBR;

typedef struct sinter {
  int ntable,ntype;
  int *itype,**map,*nexclude,**exclude;
  double **vtab,**dvtab,**d2vtab,**switchl,**switchs;
  double *rminl_cut2,*rmaxl_cut2,*dx2tab;
  double *rmins_cut2,*rmaxs_cut2;
  double *rmaxs_skin2,*rmaxl_skin2,*rminl_skin2;
  double *dxn,*dyn,*dzn;
  double skin_ter,rheal;
  /* electrostatic */
  double *coultab,*dcoultab,*switchle,*switchse,dx2coul;
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
  double **fcmat,*dn,*en,*en2,*freq,*fric,*fpos;
  double **molmat,**decmod,**permod,*fricmod,*fricmodq;
  double *suspx,*suspy,*suspz;
  double *avg_freq,*avg_strn;
} NMODES;
