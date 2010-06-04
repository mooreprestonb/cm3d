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

#include "md.h"

/* #define ANALYTIC */
/* #define WRITE_TABLE */

#ifdef ANALYTIC
static double alpha_ewald;
#endif

void set_charges(SIMPARMS *simparms,COORDS *coords,INTER *inter,
		 NGBR *ngbr,STARTUP *startup)
{
  int i;
  double xt,f,df,d2f,br,sm,vshift;
  LINE line;
#ifdef WRITE_TABLE
  FILE *fp;
#endif

#ifdef ANALYTIC
  alpha_ewald = alp_ewald;
#endif

  /* set min and max cutoff to use for special coulumn cutoff */
  inter->rmins_coul2 = startup->rcute_min*startup->rcute_min;
  inter->rmaxs_coul2 = startup->rcute_resp*startup->rcute_resp;
  inter->rminl_coul2 = ((startup->rcute_resp-inter->rheal)*
			(startup->rcute_resp-inter->rheal));
  inter->rmaxl_coul2 = startup->rcute_max*startup->rcute_max;

  inter->rmaxs_skin2_e = ((startup->rcute_resp+inter->skin_ter)*
			  (startup->rcute_resp+inter->skin_ter));
  inter->rmaxl_skin2_e = ((startup->rcute_max+inter->skin_ter)*
			  (startup->rcute_max+inter->skin_ter));
  inter->rminl_skin2_e = ((startup->rcute_resp-inter->rheal-inter->skin_ter)*
			  (startup->rcute_resp-inter->rheal-inter->skin_ter));
  
  inter->dx2coul = (((double)inter->ntable-3.)/
		    (inter->rmaxl_coul2-inter->rmins_coul2));
  
  inter->coultab  = (double *)cmalloc(3*inter->ntable*sizeof(double));
  inter->dcoultab = inter->coultab + inter->ntable;
  inter->switchse = inter->coultab + inter->ntable*2;  
  simparms->mem_bytes += 3*inter->ntable*sizeof(double);

  for(i=0;i<inter->ntable;i++){
    xt = i/inter->dx2coul + inter->rmins_coul2;
    if(xt<0){
      sprintf(line,"Real space charge table values are negative! xt = %g",xt);
      md_error(line);
    }
    if(simparms->iperd == 3){
      func_kerf(sqrt(xt),&f,&df,&d2f,startup->alp_ewald);
    } else {
      func_coul(sqrt(xt),&f,&df,&d2f);
    }
    inter->coultab[i] = f;
    inter->dcoultab[i] = df;
    if(xt<=inter->rminl_coul2){
      inter->switchse[i] = 1.; 
    }
    if(xt>inter->rminl_coul2 && xt < inter->rmaxs_coul2){
      br = (sqrt(xt)-sqrt(inter->rminl_coul2))/inter->rheal; 
      sm = (1.+ br*br*(2.*br-3.));
      inter->switchse[i] = sm; 
    }
    if(xt >= inter->rmaxs_coul2){
      inter->switchse[i] = 0.; 
    }
  }

  if(startup->ishift){
    vshift=inter->coultab[inter->ntable-3];
    for(i=0;i<inter->ntable;i++){
      inter->coultab[i] -= vshift; 
    }
  }
#ifdef WRITE_TABLE
  fp = fopen("table_coul.dat","w");
  for(i= 0;i<inter->ntable;i++){
    xt = i/inter->dx2coul + inter->rmins_coul2;
    fprintf(fp,"%d %g %g %g %g\n",i,sqrt(xt),inter->coultab[i],
	    inter->dcoultab[i],inter->switchse[i]);
  }
  fclose(fp);
#endif

  /* set up vector of charged particles */
  coords->coul.ncharge = 0;
  coords->coul.icharge = (int *)cmalloc(simparms->natoms*sizeof(int));
  for(i=0;i<simparms->natoms;i++){
    if(coords->qch[i] != 0.){
      coords->coul.icharge[coords->coul.ncharge++] = i;
    }
  }
  if(coords->coul.ncharge==0)
    free(coords->coul.icharge);
  else
    coords->coul.icharge = (int *)realloc(coords->coul.icharge,
					  coords->coul.ncharge*sizeof(int));
  
  simparms->mem_bytes += coords->coul.ncharge*sizeof(int);
  sprintf(line,"Bytes of memory allocated so far = %d",simparms->mem_bytes);
  md_stdout(line);

  if(simparms->iperd == 3){
    ewald_setup(simparms,&coords->coul,coords->qch,startup->kmax,
		startup->alp_ewald);
  }
}

/*-----------------------------------------------------------------------*/
void set_ecorr(SIMPARMS *simparms,COORDS *coords,int *nexcl,int **excl)
{
  int i,j,nreal;
  LINE line;

  if(simparms->rank==0){ md_stdout("Setting up electrostatic exclusion list");}
  nreal = 0;
  for(i=0;i<simparms->natoms;i++) if(coords->qch[i] != 0.0) nreal++;
  if(simparms->rank==0){ 
    sprintf(line,"Setting up electrostatic exclusion list (ncharge = %d)",
	    nreal);
    md_stdout(line);
  }
  
  coords->coul.nexclude = (int *)cmalloc(simparms->natoms*sizeof(int));
  coords->coul.exclude  = (int **)cmalloc(simparms->natoms*sizeof(int *));
  simparms->mem_bytes += simparms->natoms*sizeof(int);
  simparms->mem_bytes += simparms->natoms*sizeof(int *);
  for(i=0;i<simparms->natoms;i++) {
    coords->coul.nexclude[i] = 0;
    coords->coul.exclude[i] = (int *)cmalloc((nexcl[i]+1)*sizeof(int));
  }

  /* count all possible exclusions (allocate at least one :-) */
  nreal = 0.;
  for(i=0;i<simparms->natoms;i++){ nreal += nexcl[i]; }
  if(nreal==0) nreal=1;
  coords->coul.iecorr=(int *)cmalloc(nreal*sizeof(int));
  coords->coul.jecorr=(int *)cmalloc(nreal*sizeof(int));
  
  nreal = coords->coul.necorr = 0;
  for(i=0;i<simparms->natoms;i++){
    if(coords->qch[i] != 0.0){
      for(j=0;j<nexcl[i];j++){
        if(coords->qch[excl[i][j]] != 0.0){
          coords->coul.exclude[i][coords->coul.nexclude[i]] = excl[i][j];
          coords->coul.nexclude[i]++;
        }
      }
      coords->coul.exclude[i][coords->coul.nexclude[i]] = -1;
    }
    if(simparms->iperd==3){
      for(j=0;j<coords->coul.nexclude[i];j++){
        if(coords->qch[i]*coords->qch[coords->coul.exclude[i][j]] != 0.){
          coords->coul.iecorr[coords->coul.necorr]=i;
          coords->coul.jecorr[coords->coul.necorr]=coords->coul.exclude[i][j];
          coords->coul.necorr++;
        }
      }
    }
    nreal += coords->coul.nexclude[i];
  }
  for(i=0;i<simparms->natoms;i++){
    coords->coul.exclude[i]=(int *)realloc(coords->coul.exclude[i],
       (coords->coul.nexclude[i]+1)
       *sizeof(int));
    simparms->mem_bytes += (coords->coul.nexclude[i]+1)*sizeof(int);
  }
  
  if(coords->coul.necorr != 0){
    coords->coul.iecorr = (int *)realloc(coords->coul.iecorr,
       coords->coul.necorr*sizeof(int));
    coords->coul.jecorr = (int *)realloc(coords->coul.jecorr,
       coords->coul.necorr*sizeof(int));
    simparms->mem_bytes += 2*coords->coul.necorr*sizeof(int);
  } else {
    free(coords->coul.iecorr);
    free(coords->coul.jecorr);
  }

  if(simparms->rank==0){
    sprintf(line,"There are %d real space electrostatic corrections",nreal);
    md_stdout(line);
    sprintf(line,"There are %d coulomb space electrostatic corrections",
	    coords->coul.necorr);
    md_stdout(line);
    sprintf(line,"Bytes of memory allocated so far = %d",simparms->mem_bytes);
    md_stdout(line);
  }
}
/*-----------------------------------------------------------------------*/
void get_ngbr_e(SIMPARMS *simparms,COORDS *coords,INTER *inter,NGBR *ngbr)
{
  int i,j,k,ncharge,ii,jj,*excl,nexcl,mem_change,mem_add;
  int icell,jcell,*inp,*jnp,iperd,ivol;
  long ibegin,iend;
  static int npso,npto;
  LINE line;

  gethinv9(coords->hmat,coords->hmati);
  ncharge=coords->coul.ncharge;
  mem_add = MAX(simparms->nlen,ncharge);
  ivol = simparms->ivol;
  iperd = simparms->iperd;

    mem_change = simparms->mem_bytes;
  if(ngbr->update == 0){
    if(ngbr->ilist==1){
      npso = npto = 0;
      count_charge(simparms,coords,inter,ngbr,&npso,&npto);
      npso = (int)(double)(npso*1.2+1.);
      npto = (int)(double)(npto*1.2+1.); /* add safety margin */
      ngbr->inpse = (int *)cmalloc(npso*sizeof(int));
      ngbr->jnpse = (int *)cmalloc(npso*sizeof(int));
      ngbr->inpte = (int *)cmalloc(npto*sizeof(int));
      ngbr->jnpte = (int *)cmalloc(npto*sizeof(int));
      simparms->mem_bytes += 2*npso*sizeof(int);
      simparms->mem_bytes += 2*npto*sizeof(int);
    }
    allocate_lists_e(ncharge,ngbr,&i);
    simparms->mem_bytes += i;
  }
  ngbr->nprse = ngbr->nprte = 0;
  ngbr->offse[0] = ngbr->offte[0] = 0;
  
  switch(ngbr->ilist){
     case 0:    /* nolist */
       break;
     case 1:    /* verlist */
       simparms->mem_bytes -= 2*npso*sizeof(int);
       simparms->mem_bytes -= 2*npto*sizeof(int);
       
       inp = cmalloc(2*simparms->nlen*sizeof(int));
       jnp = inp + simparms->nlen;
       if(simparms->nocell_all==0){
         /* decompose i so that they is roughly the same number of
            i-j pairs on each processors */
         decomp_trig(ncharge,simparms->size,simparms->rank,&ibegin,&iend);
         k = 0;
         for(i=ibegin;i<iend;i++){
           ii = coords->coul.icharge[i];
           excl = &(coords->coul.exclude[ii][0]);
           nexcl = *(excl);
           if(npso < ngbr->nprse + ncharge){
             npso += ncharge;
             ngbr->inpse = (int *)realloc(ngbr->inpse,npso*sizeof(int));
             ngbr->jnpse = (int *)realloc(ngbr->jnpse,npso*sizeof(int));
           }
           if(npto < ngbr->nprte + ncharge){
             npto += ncharge;
             ngbr->inpte = (int *)realloc(ngbr->inpte,npto*sizeof(int));
             ngbr->jnpte = (int *)realloc(ngbr->jnpte,npto*sizeof(int));
           }
           
           for(j=i+1;j<ncharge;j++){
             jj = coords->coul.icharge[j];
             if(nexcl == jj){
               nexcl = *(++excl);
             } else {
               inp[k] = ii;
               jnp[k] = jj;
               k++;
               if(k==simparms->nlen){
                 ngbr_liste(k,inp,jnp,iperd,coords,inter,ngbr,ivol);
                 k=0;
               }
             }
           }
           ngbr_liste(k,inp,jnp,iperd,coords,inter,ngbr,ivol);
           k=0;
           ngbr->offse[i+1] = ngbr->nprse;  
           ngbr->offte[i+1] = ngbr->nprte;
         }
       } else {
         get_lnklist_e(simparms->natoms,ngbr->ncells,ngbr->hlist,ngbr->lnklist,
            iperd,coords->hmat,coords->px,coords->py,coords->pz,coords->qch);
         k = 0;
         for(i=0;i<ngbr->npaircelle;i++){
           icell = ngbr->icelle[i];
           jcell = ngbr->jcelle[i];
           
           ii = ngbr->hlist[icell];
           
           while(ii!=-1){
             if(icell==jcell){
               jj = ngbr->lnklist[ii];
             } else {
               jj = ngbr->hlist[jcell];
             }
             if(npso < ngbr->nprse + mem_add){
               npso += mem_add;
               ngbr->inpse = (int *)realloc(ngbr->inpse,npso*sizeof(int));
               ngbr->jnpse = (int *)realloc(ngbr->jnpse,npso*sizeof(int));
             }
             if(npto < ngbr->nprte + mem_add){
               npto += mem_add;
               ngbr->inpte = (int *)realloc(ngbr->inpte,npto*sizeof(int));
               ngbr->jnpte = (int *)realloc(ngbr->jnpte,npto*sizeof(int));
             }
             while(jj!=-1){
               if(!exclij(ii,jj,coords->coul.nexclude,coords->coul.exclude)){ 
                 inp[k] = ii;
                 jnp[k] = jj;
                 k++;
                 
                 if(k==simparms->nlen){
                   ngbr_liste(k,inp,jnp,iperd,coords,inter,ngbr,ivol);
                   k = 0;
                 }
               }
               jj = ngbr->lnklist[jj];
             } /* end while(jj!=-1) */
             ii = ngbr->lnklist[ii];
           } 
         }
         ngbr_liste(k,inp,jnp,iperd,coords,inter,ngbr,ivol);
       }
       simparms->mem_bytes += 2*npso*sizeof(int);
       simparms->mem_bytes += 2*npto*sizeof(int);
       free(inp);
     case 2:   /* linklist */
       break;
  } /* end switch */
  
  if(simparms->rank==0 && (mem_change-simparms->mem_bytes)!= 0 ){
    md_stdout("Electrostatic Neighbor list memory changed!");
    sprintf(line,"Bytes of memory allocated so far = %d",simparms->mem_bytes);
    md_stdout(line);
  }
}
/*---------------------------------------------------------*/
void allocate_lists_e(int ncharge,NGBR *ngbr,int *mem)
{
  *mem = 0;
  ncharge++; /* add one for last ncharge */
  ngbr->offse = (int *)cmalloc(ncharge*sizeof(int));
  ngbr->offte = (int *)cmalloc(ncharge*sizeof(int));
  *mem += 2*ncharge*sizeof(int);
}

/*---------------------------------------------------------*/
void count_charge(SIMPARMS *simparms,COORDS *coords,INTER *inter,NGBR *ngbr,
		  int *npso,int *npto)
{
  int i,j,k,ncharge,ii,jj,*excl,nexcl;
  int icell,jcell,*inp,*jnp;
  long ibegin,iend;

  ncharge = coords->coul.ncharge;
  inp = cmalloc(2*simparms->nlen*sizeof(int));
  jnp = inp + simparms->nlen;

  if(simparms->nocell_all==0){
    decomp_trig((long)ncharge,simparms->size,simparms->rank,&ibegin,&iend);
    
    k = 0;
    for(i=ibegin;i<iend;i++){
      ii = coords->coul.icharge[i];
      excl = &(coords->coul.exclude[ii][0]);
      nexcl = *(excl);
      for(j=i+1;j<ncharge;j++){
        if(nexcl == j){
          nexcl = *(++excl);
        } else {
          inp[k] = i;
          jnp[k] = j;
          k++;
          if(k==simparms->nlen){
            cnt_ngbr_liste(k,inp,jnp,simparms->iperd,coords,inter,
               npso,npto,simparms->ivol);
            k=0;
          }
        }
      }
      cnt_ngbr_liste(k,inp,jnp,simparms->iperd,coords,inter,npso,npto,
         simparms->ivol);
      k=0;
    }
  } else {
    get_lnklist_e(simparms->natoms,ngbr->ncells,ngbr->hlist,ngbr->lnklist,
       simparms->iperd,coords->hmat,
       coords->px,coords->py,coords->pz,coords->qch);
    
    k = 0;
    for(i=0;i<ngbr->npaircelle;i++){
      icell = ngbr->icelle[i];
      jcell = ngbr->jcelle[i];
      ii = ngbr->hlist[icell];
      while(ii!=-1){
        if(icell==jcell){
          jj = ngbr->lnklist[ii];
        } else {
          jj = ngbr->hlist[jcell];
        }
        while(jj!=-1){
          if(!exclij(ii,jj,coords->coul.nexclude,coords->coul.exclude)){ 
            inp[k] = ii;
            jnp[k] = jj;
            k++;
            
            if(k==simparms->nlen){
              cnt_ngbr_liste(k,inp,jnp,simparms->iperd,coords,inter,
                 npso,npto,simparms->ivol);
              k = 0;
            }
          }
          jj = ngbr->lnklist[jj];
        } /* end while(jj!=-1) */
        ii = ngbr->lnklist[ii];
      } 
    }
    cnt_ngbr_liste(k,inp,jnp,simparms->iperd,coords,inter,npso,npto,
       simparms->ivol);
  }
  free(inp);
}




