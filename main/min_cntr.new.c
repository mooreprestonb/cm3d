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

/* routine to minimize a function in one dimension */

#include "md.h"

void mal_verify(int);
#define DXSTEP .01
#define IRSTEP 100
#define PRINT_FREQ

/* static varibales to go over the head of the minimize routines */
static NMODES nmodes;
static ENERGIES energies;
static SIMPARMS *simparms_s;
static COORDS *coords_s;
static INTER *inter_s;
static NGBR *ngbr_s;
static int ndim;
static double *svec;

double func(double *x);
double func_ir(double *x);
void dfunc(double *x,double *dx);
void print_potential(SIMPARMS *,COORDS *coords);
void load_fix(COORDS *coords,int *ndim, double **p);

/*-----------------------------------------------------------*/

void min_cntr_new(char *command,FILENAMES *filenames,SIMPARMS *simparms,
	      WRITE_STEP *write_step,COORDS *coords,INTER *inter,NGBR *ngbr)
{
  int i,j,iter,istep;
  double *p,**xmat;
  double tol,fret;
  double *freq_start;
  LINE line;

  simparms->iensemble = 0; /* set to NVE no matter what */
  simparms->ntherm = 0;
  simparms->ivol = 0;

  if(simparms->rank==0){
    readpos(command,filenames->initfile,simparms,
	    &energies.estep,&energies,coords);
    if(simparms->iperd==0){
      rotate(simparms->natoms,coords->px,coords->py,coords->pz,
	     coords->amass);
      zero_angm(simparms->natoms,coords->px,coords->py,coords->pz,
         coords->vx,coords->vy,coords->vz,coords->amass);
    }
    zero_averages(simparms,&energies);
    openfiles(simparms,write_step,filenames);
  }
  /* broadcast possitions and velocities, masses and charge 
     use the fact that they are contiguous in memory */
#ifdef PARA
  bcast_system(simparms,coords);
#endif

  gethinv9(coords->hmat,coords->hmati);

  ngbr->ilist = 0; /* nolist */
  step_init(simparms,coords,&energies,inter,ngbr);
  geteng(simparms,coords,&energies);

  save_conf(filenames->fconf,simparms->natoms,coords);

  /* output stuff the user might want to know */
  if(simparms->rank==0){
    output_simparm(simparms,write_step,simparms->istep);

    output_initval(simparms,coords,ngbr,&energies);
    
    check_vals(9,energies.potra,energies.poter,energies.zke,
	       energies.potn,energies.zken,energies.potv,energies.zkev,
	       energies.tiout,energies.prsi);
    
    set_econi(&energies);
  }

  ndim = 3*simparms->natoms;
  p = cmalloc(ndim*sizeof(double));
  for(i=0;i<ndim;i++) p[i] = coords->px[i];
  tol = simparms->min_tol;
  
  simparms_s = simparms;
  coords_s = coords;
  ngbr_s = ngbr;
  inter_s = inter;

  energies.acpu = cputime();
  energies.acpu = 0.;
  sprintf(line,"Bytes of memory allocated so far is %d",simparms->mem_bytes);
  md_stdout(line);
  md_stdout("Starting Minimization\n");

  switch(simparms->imin_type){
  case 0: /* powell method */
    printf("Initial Energy = %g  (Powell method)\n",func(p));
    xmat = dmatrix(0,ndim-1,0,ndim-1);
    for(i=0;i<ndim*ndim;i++) xmat[0][i] = 0.;
    for(i=0;i<ndim;i++) xmat[i][i] = DXSTEP;
    iter = ITMAX;
    for(istep=1;istep<=simparms->nstep && iter==ITMAX;istep++){
      powell(p,xmat,ndim,tol,&iter,&fret,func,ITMAX);
      printf("Energy at step %d = %g K (iter = %d)\n",istep,func(p),--iter);
      energies.acpu += cputime();
      if(write_step->nconf && simparms->istep%(write_step->nconf) == 0){
	save_conf(filenames->fconf,simparms->natoms,coords);
      }
    }
    printf("Final Energy = %g (iter = %d)\n",func(p),iter);
    break;
  case 1: /* conjugate gradient methods */
    printf("Initial Energy = %g (Conjugate gradient)\n",func(p));
    iter = ITMAX;
    for(istep=1;istep<=simparms->nstep && iter==ITMAX;istep++){
      frprmn(p,ndim,tol,&iter,&fret,func,dfunc,ITMAX);
      printf("Energy at step %d = %g K (iter = %d)\n",istep,func(p),iter);
      energies.acpu += cputime();
      if(write_step->nconf && simparms->istep%(write_step->nconf) == 0){
	save_conf(filenames->fconf,simparms->natoms,coords);
      }
    }
    printf("Final Energy = %g (iter = %d)\n",func(p),iter);
    break;
  case 2: /* simplex method */
    printf("Initial Energy = %g (Simplex Method)\n",func(p));
    xmat = dmatrix(0,ndim,0,ndim-1);
    for(i=0;i<(ndim+1);i++) {
      for(j=0;j<ndim;j++) xmat[i][j] = p[i];
    }
    for(i=0;i<ndim;i++) xmat[i][i] += DXSTEP;
    for(i=0;i<ndim+1;i++){
      p[i] = func(xmat[i]);
    }
    iter = ITMAX;
    for(istep=1;istep<=simparms->nstep && iter==ITMAX;istep++){
      amoeba(xmat,p,ndim,tol,func,&iter,ITMAX*ndim);
      energies.acpu += cputime();
      if(iter<ITMAX) break;
      if(write_step->nconf && simparms->istep%(write_step->nconf) == 0){
	save_conf(filenames->fconf,simparms->natoms,coords);
      }
    }
    for(i=0;i<ndim;i++) coords_s->px[i] = xmat[0][i];
    printf("Final Energy = %g (iter = %d)\n",func(xmat[0]),iter/ndim);
    break;

  case 3: /* minimize IR */
    printf("IR calculation first minimizing structure\n");

#ifdef MIN_FIRST
    printf("Minimizing structure\n");
    printf("Initial Energy = %g  (Powell method)\n",func(p));
    xmat = dmatrix(0,ndim-1,0,ndim-1);
    for(i=0;i<ndim*ndim;i++) xmat[0][i] = 0.;
    for(i=0;i<ndim;i++) xmat[i][i] = DXSTEP;
    for(istep=1;istep<=simparms->nstep;istep++){
      powell(p,xmat,ndim,tol,&iter,&fret,func,ITMAX);
      printf("Energy at step %d = %g K (iter = %d)\n",istep,func(p),iter);
      energies.acpu += cputime();
      if(iter<ITMAX) break;
    }
    printf("Final Energy = %g (iter = %d)\n",func(p),iter);
    free_dmatrix(xmat,0,ndim-1,0,ndim-1);
#endif

    free(p);
    /* allocate NM memory */
    nmodes.flo = 0;
    nma_allocate(simparms,&nmodes);
    sprintf(line,"Bytes of memory allocated so far is %d",simparms->mem_bytes);
    md_stdout(line);

    svec = cmalloc(simparms_s->natoms*3*sizeof(double));
    read_freqfile(filenames->nmfile,simparms->natoms*3,svec);
    
    /* inital frequencies */

    printf("Calculating initial frequencies \n");
    fmat(simparms_s,coords_s,inter_s,&nmodes);    
    /* diagonalize system to get eigenvalues */
    rs_me(simparms->natoms*3,nmodes.dn,nmodes.fcmat,0);    
    ceigsrtv(nmodes.dn,nmodes.fcmat,simparms->natoms*3);

#ifdef PRINT_FREQ
    freq_start = cmalloc(simparms->natoms*3*sizeof(double));
    for(i=0;i<simparms->natoms*3;i++){
      freq_start[i] = nmodes.dn[i];
      printf("%d %g %g\n",i,nmodes.dn[i],svec[i]);
    }
#endif

    printf("Finding adjustable parameters\n");
    load_fix(coords,&ndim,&p);

    printf("Initial Difference = %g  (Powell method)\n",func_ir(p));
    print_potential(simparms,coords);
    xmat = dmatrix(0,ndim-1,0,ndim-1);
    for(i=0;i<ndim*ndim;i++) xmat[0][i] = 0.;
    for(i=0;i<ndim;i++) xmat[i][i] = IRSTEP;
    for(istep=1;istep<=simparms->nstep;istep++){
      powell(p,xmat,ndim,tol*tol,&iter,&fret,func_ir,ITMAX);
      printf("Difference at step %d = %g K (iter = %d)\n",
	     istep,func_ir(p),iter);
      energies.acpu += cputime();
      if(iter<ITMAX) break;
    }
    printf("Final Difference = %g (iter = %d)\n",func_ir(p),iter);

    print_potential(simparms,coords);

#ifdef JUNK
    fmat(simparms_s,coords_s,inter_s,&nmodes);
    /* diagonalize system to get eigenvalues */
    rs_me(simparms->natoms*3,nmodes.dn,nmodes.fcmat,0);
    ceigsrtv(nmodes.dn,nmodes.fcmat,simparms->natoms*3);
    for(i=0;i<simparms->natoms*3;i++){
      printf("%d %g %g %g\n",i,nmodes.dn[i],svec[i],freq_start[i]);
    }
#endif

    break;
  }

  energies.acpu += cputime();
  sprintf(line,"CPU time = %g (%g min per CPU sec)\n",energies.acpu,
	  energies.acpu/(istep*ITMAX+iter));
  md_stdout(line);  

  write_coord(filenames->initfile,simparms->natoms,coords);
  save_conf(filenames->fconf,simparms->natoms,coords);

  exit(0);
}
/*-------------------------------------------------------------------*/

#define KCOM 10
double func_new(double *x)
{
  int i,n;
  double am,fun,ax,ay,az;

  for(i=0;i<ndim;i++) {
    coords_s->px[i] = x[i];
  }
  n = ndim/3;

  /* add a harmonic term to the center of mass */
  fun = 0.;
  ax = ay = az = 0.;
  for(i=0;i<n;i++){
    am = coords_s->amass[i];
    ax += am*coords_s->px[i];
    ay += am*coords_s->py[i];
    az += am*coords_s->pz[i];
  }
  am = KCOM*(ax*ax+ay*ay+az*az);

  getvireal(simparms_s,coords_s,&energies,inter_s,ngbr_s,0);

  fun = energies.U + energies.UI + am;
  return(fun);
}
/*-------------------------------------------------------------------*/

void dfunc_new(double *x,double *y)
{
  int i,n,n2;
  double ax,ay,az,am;

  n = ndim/3;
  n2 = n+n;
  ax = ay = az = 0.;
  for(i=0;i<n;i++){
    am = coords_s->amass[i];
    ax += am*coords_s->px[i];
    ay += am*coords_s->py[i];
    az += am*coords_s->pz[i];
  }
  ax *= 2.*KCOM;
  ay *= 2.*KCOM;
  az *= 2.*KCOM;

  for(i=0;i<ndim;i++) {
    coords_s->px[i] = x[i];
  }

  force(simparms_s,coords_s,inter_s,ngbr_s,1);
  force(simparms_s,coords_s,inter_s,ngbr_s,2);
  force(simparms_s,coords_s,inter_s,ngbr_s,3);
  getvireal(simparms_s,coords_s,&energies,inter_s,ngbr_s,1);

  for(i=0;i<ndim;i++){
    y[i]=coords_s->fxa[i]+coords_s->fxr[i]+coords_s->fxl[i];
  }

  /* add harmonic term of the center of mass */
  for(i=0;i<n;i++){
    am = coords_s->amass[i];
    y[i   ] += ax*am;
    y[i+n ] += ay*am;
    y[i+n2] += az*am;
  }
}
  
/*-------------------------------------------------------------------*/
double func_ir_new(double *x)
{
  int i,ii,nfix,i_tors;
  double dq,fun;

  nfix = 0;

  for(i=0;i<coords_s->bonds.itypes;i++){
    if(coords_s->bonds.ifix[i]==1) coords_s->bonds.fkbond[i] = x[nfix++];
  }
  for(i=0;i<coords_s->bends.itypes;i++){
    if(coords_s->bends.ifix[i]==1) coords_s->bends.fkbend[i] = x[nfix++];
  }
  for(i=0;i<coords_s->bondxs.itypes;i++){
    if(coords_s->bondxs.ifix[i]==1) coords_s->bondxs.fkbondx[i] = x[nfix++];
  }
  for(i=0;i<coords_s->tors.itypes;i++){
    if(coords_s->tors.ifix[i]==1){
      i_tors = coords_s->tors.ind_pot_num[coords_s->tors.ind_pot[i]];
      switch(coords_s->tors.ind_pot[i]){
      case 0:
	coords_s->tors.fktors[i_tors] = x[nfix++];
	break;
      case 1:
	for(ii=0;ii<MAXPOWER_TORS;ii++){
	  coords_s->tors.vphi[i_tors][ii] = x[nfix++];
	}
	break;
      }
    }
  }
  fmat(simparms_s,coords_s,inter_s,&nmodes);
  
  /* diagonalize system to get eigenvalues */
  rs_me(simparms_s->natoms*3,nmodes.dn,nmodes.fcmat,0);    
  ceigsrt(nmodes.dn,nmodes.fcmat,simparms_s->natoms*3);

  fun = 0;
  for(i=0;i<simparms_s->natoms*3;i++){
    dq = svec[i] - nmodes.dn[i];
    fun += dq*dq;
  }
  return(fun);
}

/*-----------------------------------------------------------*/

void print_potential_new(SIMPARMS *simparms,COORDS *coords)
{
  int i,ii,ioff,nb,ia,ib,ic,id,i_tors;

  for(i=0;i<coords->bonds.itypes;i++){
    for(nb=0;nb<coords->bonds.nbonds;nb++){
      if(i==coords->bonds.idxbond[nb]){
	ia = simparms->itype[coords->bonds.ibond[nb]];
	ib = simparms->itype[coords->bonds.jbond[nb]];
	printf("Bond %s %s -- fk = %g K\n",simparms->atom[ia],
	       simparms->atom[ib],coords->bonds.fkbond[i]);
	break;
      }
    }
  }
  ioff = i;
  for(i=0;i<coords->bends.itypes;i++){
    for(nb=0;nb<coords->bends.nbends;nb++){
      if(i==coords->bends.idxbend[nb]){
	ia = simparms->itype[coords->bends.ibend[nb]];
	ib = simparms->itype[coords->bends.jbend[nb]];
	ic = simparms->itype[coords->bends.kbend[nb]];
	printf("Bend %s %s %s -- fk = %g K\n",
	       simparms->atom[ia],simparms->atom[ib],
	       simparms->atom[ic],coords->bends.fkbend[i]);
	break;
      }
    }
  }
  ioff += i;
  for(i=0;i<coords->bondxs.itypes;i++){
    for(nb=0;nb<coords->bondxs.nbondxs;nb++){
      if(i==coords->bondxs.idxbox[nb]){
	ia = simparms->itype[coords->bondxs.ibondx[nb]];
	ib = simparms->itype[coords->bondxs.jbondx[nb]];
	ic = simparms->itype[coords->bondxs.kbondx[nb]];
	printf("Bond Cross %s %s %s -- fk = %g K\n",
	       simparms->atom[ia],simparms->atom[ib],
	       simparms->atom[ic],coords->bondxs.fkbondx[i]);
	break;
      }
    }
  }
  ioff += i;
  for(i=0;i<coords_s->tors.itypes;i++){
    for(nb=0;nb<coords->tors.ntorss;nb++){
      if(i==coords->tors.idxtors[nb]){
	ia = simparms->itype[coords->tors.itors[nb]];
	ib = simparms->itype[coords->tors.jtors[nb]];
	ic = simparms->itype[coords->tors.ktors[nb]];
	id = simparms->itype[coords->tors.ltors[nb]];
	printf("Torsion %s %s %s %s -- ",simparms->atom[ia],
	       simparms->atom[ib],simparms->atom[ic],simparms->atom[id]);
	break;
      }
    }
    i_tors = coords->tors.ind_pot_num[coords_s->tors.ind_pot[i]];
    switch(coords_s->tors.ind_pot[i]){
    case 0:
      printf("fk = %g\n",coords_s->tors.fktors[i_tors]);
      ioff++;
      break;
    case 1:
      printf("\n\tseries :");
      for(ii=0;ii<MAXPOWER_TORS;ii++){
	printf("%4g ",coords_s->tors.vphi[i_tors][ii]);
	ioff++;
      }
      printf("\n");
      break;
    }
  }
}

/*-----------------------------------------------------------*/
void load_fix_new(COORDS *coords,int *ndim, double **p)
{
  /* count and load up variables to minimize */

  int i,ii,maxdim,i_tors;

  maxdim  = coords->bonds.itypes;
  maxdim += coords->bondxs.itypes;
  maxdim += coords->bends.itypes;
  maxdim += 7*coords->tors.itypes;

  (*p) = cmalloc(maxdim*sizeof(double));
  
  *ndim = 0;
  for(i=0;i<coords->bonds.itypes;i++){
    if(coords->bonds.ifix[i]==1) (*p)[(*ndim)++] = coords->bonds.fkbond[i];
  }

  for(i=0;i<coords->bends.itypes;i++){
    if(coords->bends.ifix[i]==1) (*p)[(*ndim)++] = coords->bends.fkbend[i];
  }

  for(i=0;i<coords->bondxs.itypes;i++){
    if(coords->bondxs.ifix[i]==1) (*p)[(*ndim)++] = coords->bondxs.fkbondx[i];
  }

  for(i=0;i<coords->tors.itypes;i++){
    if(coords->tors.ifix[i]==1){
      i_tors = coords->tors.ind_pot_num[coords->tors.ind_pot[i]];
      switch(coords->tors.ind_pot[i]){
      case 0:
	(*p)[(*ndim)++] = coords->tors.fktors[i_tors];
	break;
      case 1:
	for(ii=0;ii<MAXPOWER_TORS;ii++){
	  (*p)[(*ndim)++] = coords->tors.vphi[i_tors][ii];
	}
	break;
      }
    }
  }
 
  if(*ndim==0){
    md_error("NO adjustable paramaters to fit spectra");
  }
  fprintf(stdout,"Number of adjustable parameters = %d\n",*ndim);

  (*p) = realloc((*p),(*ndim)*sizeof(double));

}
/*-----------------------------------------------------------*/
