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

/* subroutines to handle io from md program */

#include "md.h"

/*------------------------------------------------------------------*/
/* determine if i and j are on the same molecule */
int same_mole_new(int iatom,int jatom,int *itype_species,int *itype_molecules)
{
  if(itype_species[iatom] != itype_species[jatom]) return 0;
  if(itype_molecules[iatom] != itype_molecules[jatom]) return 0;

  return 1;
}
/*--------------------------------------------------------------*/
/* open file for writing and check to make sure that the 
   file was opened correctly! */
FILE *cfopenw_new(char *name)
{
  LINE line;
  FILE *fp;

  if((fp=fopen(name,"r"))!=NULL){
    sprintf(line,"File %s exists! overwriting!",name);
    md_warning(line);
    fclose(fp);
  }

  if((fp=fopen(name,"w"))==NULL){
    sprintf(line,"can't open %s for writing",name);
    md_error(line);
  }
  return fp;
}

/*--------------------------------------------------------------*/
#ifndef cmalloc
void *cmalloc_new(size_t mem)
{
  void *p;
  LINE line;

#ifdef DMALLOC
  dmalloc_verify(NULL);
#endif
  if(mem==0) return NULL;
  
  p = malloc(mem);
  if(p == NULL){
    sprintf(line,"can't allocate %ld bytes of memory",(long)mem);
    md_error(line);
  }
#ifdef DMALLOC
  dmalloc_verify(NULL);
#endif
  return p;
}
#endif
/*--------------------------------------------------------------*/
#ifndef crealloc
void *crealloc_new(void *ptr,size_t mem)
{
  void *p;
  LINE line;

  if(ptr==NULL) md_warning("Reallocating a NULL pointer");

  if(mem==0){
    free(ptr);
    return NULL;
  }
  if((p = realloc(ptr,mem))==NULL){
    sprintf(line,"can't reallocate %ld bytes of memory",(long)mem);
    md_error(line);
  }
  return p;
}
#endif
/*--------------------------------------------------------------*/
void allocate_coords_new(SIMPARMS *simparms,COORDS *coords)
{
  int n;
  double *p;
  LINE line;

  n = simparms->natoms;
  /* allocate memory for positions and velocities and force variables*/
  p = (double *)cmalloc(9*n*sizeof(double));
  coords->px = p;
  coords->py = p + n;
  coords->pz = p + 2*n;
  coords->vx = p + 3*n;
  coords->vy = p + 4*n;
  coords->vz = p + 5*n;
  coords->amass = p + 6*n;
  coords->qch   = p + 7*n;
  coords->alpha = p + 8*n;

  simparms->mem_bytes += 9*n*sizeof(double);

  /* allocate forces */
  coords->fxs = (double *)cmalloc(15*n*sizeof(double));
  coords->fys = coords->fxs+ 1*n;
  coords->fzs = coords->fxs+ 2*n;
  coords->fxa = coords->fxs+ 3*n;
  coords->fya = coords->fxs+ 4*n;
  coords->fza = coords->fxs+ 5*n;
  coords->fxl = coords->fxs+ 6*n;
  coords->fyl = coords->fxs+ 7*n;
  coords->fzl = coords->fxs+ 8*n;
  coords->fxr = coords->fxs+ 9*n;
  coords->fyr = coords->fxs+10*n;
  coords->fzr = coords->fxs+11*n;
  coords->fxt = coords->fxs+12*n;
  coords->fyt = coords->fxs+13*n;
  coords->fzt = coords->fxs+14*n;

  simparms->mem_bytes += 15*n*sizeof(double);

  if(simparms->rank==0){
    sprintf(line,"Bytes of memory allocated after coordinates  = %d",
	    simparms->mem_bytes);
    md_stdout(line);
  }

  /* XXX JB */

        coords->p = calloc(simparms->natoms*3, sizeof(double));
        coords->v = calloc(simparms->natoms*3, sizeof(double));
        coords->fs = calloc(simparms->natoms*3, sizeof(double));
        coords->fl = calloc(simparms->natoms*3, sizeof(double));
        coords->fa = calloc(simparms->natoms*3, sizeof(double));
        coords->fr = calloc(simparms->natoms*3, sizeof(double));
        coords->ft = calloc(simparms->natoms*3, sizeof(double));


}
void free_coords_new(COORDS *coords)
{
  free(coords->px);
  free(coords->fxs);
}
/*--------------------------------------------------------------*/
void allocate_extend_new(SIMPARMS *simparms,COORDS *coords)
{
  if(simparms->nbar>0) {
    coords->bar  = (double *)cmalloc((simparms->nbar)*sizeof(double));
    coords->vbar = (double *)cmalloc((simparms->nbar)*sizeof(double));
  }
  coords->eta  = dmatrix(0,simparms->ntherm,0,simparms->nchain);
  coords->veta = dmatrix(-1,simparms->ntherm,0,simparms->nchain);
  simparms->mem_bytes += 2*simparms->nbar*sizeof(double);
  simparms->mem_bytes += 2*(simparms->ntherm+2)*(simparms->nchain+1)*
    sizeof(double);
}
void free_extend_new(SIMPARMS *simparms,COORDS *coords)
{
  if(simparms->nbar>0){
    free(coords->bar);
    free(coords->vbar);
    free(coords->mbar);
  }
  if(simparms->ntherm>0){
    free_dmatrix(coords->eta,0,simparms->ntherm,0,simparms->nchain);
    free_dmatrix(coords->veta,-1,simparms->ntherm,0,simparms->nchain);
    free(coords->p2mt);
    free_dmatrix(coords->feta,0,simparms->ntherm-1,0,simparms->nchain-1);
    free_dmatrix(coords->gkt,0,simparms->ntherm-1,0,simparms->nchain-1);
    free_dmatrix(coords->mnh,0,simparms->ntherm-1,0,simparms->nchain-1);
  }
}
/*--------------------------------------------------------------*/
void allocate_exttab_new(SIMPARMS *simparms,COORDS *coords,int ntable)
{
  int ntypes = simparms->ntypes;

  coords->vtab_extern  = dmatrix(0,ntypes,0,ntable-1);
  coords->dvtab_extern = dmatrix(0,ntypes,0,ntable-1);
  /* coords->dv2tab_extern = dmatrix(0,ntypes,0,ntable-1); */

  coords->vtab_extern_e  = (double *)cmalloc(2*ntable*sizeof(double));
  coords->dvtab_extern_e = coords->vtab_extern_e+ntable;;
  /*  coords->dv2tab_extern_e = (double *)cmalloc(2*ntable*sizeof(double));*/

  simparms->mem_bytes += 4*ntable*sizeof(double);
}
void free_exttab_new(COORDS *coords)
{
  free(coords->dvtab_extern);
  free(coords->dvtab_extern_e);
}

/*--------------------------------------------------------------*/
void rmass_new(int natoms,double *fx,double *fy,double *fz,double *amass)
{
  int i;
  double rmass;
  for(i=0;i<natoms;i++){
    rmass = 1./amass[i];
    fx[i] *= rmass;
    fy[i] *= rmass;
    fz[i] *= rmass;
  }
}
/*---------------------------------------------------------------------*/
/* routine to make sure that anint give the correct numbers!!!!*/
void test_anint_new(void)
{
  int ierr=0;
  double a;

  a = anint(-1.49999999);
  if(a != -1.0) ierr=1;
  a = anint(-0.50000001);
  if(a != -1.0) ierr=1;
  a = anint( 1.49999999);
  if(a !=  1.0) ierr=1;
  a = anint( 0.50000001);
  if(a !=  1.0) ierr=1;

  if(ierr==1) md_error("anint doesn't work!\n\tYou need to recompile with the correct defines!");
}
/*------------------------------------------------------------------------*/

void free_all_pointers_new(SIMPARMS *simparms,COORDS *coords,
		       INTER *inter,NGBR *ngbr)
{
  int i;

  free_coords(coords);
  if(simparms->nfreeze>0) free(coords->ifreeze);    
  free(simparms->itype);  free(simparms->iatom);
  free(simparms->imolecule);  free(simparms->ispecies);
  free(simparms->igroup);  free(simparms->group);  free(simparms->mole);

  for(i=0;i<simparms->natoms;i++){
    free(inter->exclude[i]); free(coords->coul.exclude[i]);
  }
  free(inter->exclude); free(inter->nexclude);
  free(coords->coul.nexclude); free(coords->coul.exclude);
  
  free(inter->itype);
  free(inter->map[0]);  free(inter->map);
  free(inter->rmaxs_cut2);  free(inter->rmins_cut2);
  free(inter->rmaxl_cut2);  free(inter->rminl_cut2);
  free(inter->rmaxs_skin2);  free(inter->rmaxl_skin2);free(inter->rminl_skin2);
  free(inter->vtab[0]);    free(inter->vtab);  
  free(inter->dvtab[0]);   free(inter->dvtab);   
  free(inter->switchs[0]); free(inter->switchs);
  free(inter->dx2tab);   free(inter->dxn);
  free(inter->coultab);

  free(coords->bonds.idxbond);  
  free(coords->bonds.ibond); free(coords->bonds.jbond);
  if(coords->bonds.nbonds>0){
    free(coords->bonds.ifix);
    free(coords->bonds.eqbond);free(coords->bonds.fkbond);
  }
  free(coords->bends.idxbend); 
  free(coords->bends.ibend);   free(coords->bends.jbend);
  free(coords->bends.kbend);
  if(coords->bends.nbends>0){
    free(coords->bends.ifix);
    free(coords->bends.eqbend);free(coords->bends.fkbend);
  }
  free(coords->tors.idxtors);
  free(coords->tors.itors);  free(coords->tors.jtors);
  free(coords->tors.ktors);  free(coords->tors.ltors);
  if(coords->tors.ntorss>0){
    free(coords->tors.eqtors); free(coords->tors.fktors);
    free(coords->tors.vphi);   free(coords->tors.ifix);
    free(coords->tors.ind_pot);free(coords->tors.ind_pot_num);
  }
  free(coords->onfo.i14dx);
  free(coords->onfo.i14); free(coords->onfo.j14);
  if(coords->onfo.n14s>0){
    free(coords->onfo.r14max2); free(coords->onfo.r14min2);
    free(coords->onfo.dx214tab);
    free_dmatrix(coords->onfo.v14tab,
		 0,coords->onfo.itypes,0,coords->onfo.n14table-1);
    free_dmatrix(coords->onfo.dv14tab,
		 0,coords->onfo.itypes,0,coords->onfo.n14table-1);
  }

  free(coords->bondxs.idxbox); 
  free(coords->bondxs.ibondx); free(coords->bondxs.jbondx);
  free(coords->bondxs.kbondx);
  if(coords->bondxs.nbondxs>0){
    free(coords->bondxs.ifix);free(coords->bondxs.fkbondx);
    free(coords->bondxs.eq1bondx); free(coords->bondxs.eq2bondx);
  }

  if(coords->coul.ncharge>0)  {
     free(coords->coul.icharge);
     free(coords->coul.cossc);
     free(coords->coul.iecorr);
     free(coords->coul.jecorr);
     free(coords->coul.ka);
     free(coords->coul.kb);
     free(coords->coul.kc);
     free(coords->coul.kup);
  }
   
  if(simparms->nbar>0) free(coords->mbar); 
  if(simparms->ntherm>0){
    free_dmatrix(coords->feta,0,simparms->ntherm-1,0,simparms->nchain-1);
    free_dmatrix(coords->gkt,0,simparms->ntherm-1,0,simparms->nchain-1);
    free_dmatrix(coords->mnh,0,simparms->ntherm-1,0,simparms->nchain-1);
    free(coords->p2mt); 
  }

  free(coords->ithm);
  free_extend(simparms,coords);

  if(simparms->iextern == 1) free_exttab(coords);

  free_lists(ngbr,coords);
}

