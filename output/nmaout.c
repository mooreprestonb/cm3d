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

/* subroutines to deal with output of nma */

#include "md.h"
#define ZERO /* remove zero frequencies from freq store */

static void qweighti(int,double,double,double[],double[],double[],double[],double[],double[],double,int);

void write_nma_start(int natoms,int npt,double dt,double frmax,
		     double frmin,int iunits,int nsave,NMODES *nmodes)
{
  int i;

  fprintf(stdout,"\n****Getting Spectrum****\n");
  fprintf(stdout,"# of atoms  = %d\n",natoms);
  fprintf(stdout,"# of species = %d\n",nmodes->nspec);
  for(i=0;i<nmodes->nspec;i++){
    fprintf(stdout,"Species %d has %d molecules and %d atoms/molecule\n",
	    i,nmodes->nmol[i],nmodes->napm[i]);
  }
  fprintf(stdout,"# of points = %d\n",npt);
  fprintf(stdout,"configurations were saved every %g ps\n",dt);
  if(iunits ==3){
    fprintf(stdout,"max frequency = %g (cm^-1)\n",frmax);
    fprintf(stdout,"min frequency = %g (cm^-1)\n",frmin);
  } else {
    fprintf(stdout,"max frequency = %g (ps^-1)\n",frmax);
    fprintf(stdout,"min frequency = %g (ps^-1)\n",frmin);
  }
  fprintf(stdout,"Saving spectra every %d configurations\n",nsave);
  fprintf(stdout,"%d configurations already done\n",nmodes->nconf_off);
}

/*------------------------------------------------------------------*/
void save_nma(SIMPARMS *simparms,FILENAMES *filenames,int npoints,
	      double freq_min,NMODES *nmodes)
{
  int i,j,nconf;
  double fac,*a2,*a3;
  double *a4,*a5;
  FILE *fpout=NULL;
  LINE line;
  
  nconf = nmodes->nconf+nmodes->nconf_off;
  if(nconf==0){
    md_error("No configurations in config file!");
  }
  /* open output file */
  fpout=cfopenw(filenames->nmfile);

  fprintf(stdout,"Saving modes file %s with %d configurations averaged\n",
	  filenames->nmfile,nconf);
  
  fac = 1./(simparms->natoms*3*nconf*nmodes->df);
  fprintf(fpout,"# %d %d %d DOS and projected DOS along axis\n",
	  simparms->natoms,nmodes->nspec,nmodes->nconf_off);
  for(i=0;i<npoints;i++){
    fprintf(fpout,"%.10g %.10g ",i*nmodes->df+freq_min,nmodes->freq[i]*fac);
    for(j=0;j<nmodes->imodes*nmodes->nspec;j++) 
      fprintf(fpout,"%.10g ",nmodes->decmod[j][i]*fac);
    fprintf(fpout,"\n");
  }
  fclose(fpout);

  /* open spectrum file */
  fac = 1./(nconf);
  fac *= 1./nmodes->df;
  if(simparms->iunits==3) fac *= 1.e8/CSPEED;

  a2 = cmalloc(4*npoints*sizeof(double));
  a3 = a2 + npoints;
  a4 = a3 + npoints;
  a5 = a4 + npoints;
/* added by HABS HA and BS 3/30/98    */
/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/*   qweighti(npoints,nmodes->df,freq_min,nmodes->fric,a2,a3,
	   simparms->temp,simparms->iunits); */
/* modified by HABS HA and BS 3/30/98    */
/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
  qweighti(npoints,nmodes->df,freq_min,nmodes->fric,nmodes->ftot,a2,a3,a4,a5,
     simparms->temp,simparms->iunits);

  fpout=cfopenw(filenames->specfile);
  fprintf(fpout,"# %d %d %d %.10g (all units are in (A^2/system)\n",
	  simparms->natoms,nmodes->nspec,nconf,nmodes->df);
/* modified by HABS HA and BS 3/30/98    */
/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
  fprintf(fpout,
     "# freq, QID, QID_R + QID_I, QID_R + Weight*QID_I, Charge, Tot, Tot R+I, Tot R+ Weight*I \n");
  for(i=0;i<npoints;i++){
    fprintf(fpout,"%14.10g %14.10g %14.10g %14.10g %14.10g %14.10g %14.10g %14.10g\n",
	    i*nmodes->df+freq_min,nmodes->fric[i]*fac,a3[i]*fac,
	    a2[i]*fac,nmodes->fpos[i]*fac,nmodes->ftot[i]*fac,a4[i]*fac,a5[i]*fac);
  }
/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
  fclose(fpout);
  free(a2);

  /* open spectrum average file */
  fac = 1./(nconf);
  if(simparms->iunits==3) fac *= 1.e8/CSPEED;
  sprintf(line,"%s_a",filenames->specfile);
  fpout=cfopenw(line);
  fprintf(fpout,"# %d %d %d average_w average_str\n",
	  simparms->natoms,nmodes->nspec,nconf);
  for(i=0;i<simparms->natoms*3;i++){
    fprintf(fpout,"%d %.10g %.10g\n",
	    i,nmodes->avg_freq[i]*fac,nmodes->avg_strn[i]*fac);
  }
  fclose(fpout);
  
  if(simparms->iunits == 3){
    fprintf(stdout,"Frequencies found between %g and %g (cm^1)\n",
	    nmodes->minfreq,nmodes->maxfreq);
  } else {
    fprintf(stdout,"Frequencies found between %g and %g (ps^-1)\n",
	    nmodes->minfreq,nmodes->maxfreq);
  }
  fprintf(stdout,"cputime/config = %g seconds\n",nmodes->acpu/nconf);

  fprintf(stdout,"\n");
  fflush(stdout); fflush(stderr);
}
/*-----------------------------------------------------------------------*/
void print_mat(int nr,int nc,double **matrix)
{
  int i,j;
  
  printf("{");
  for(i=0;i<nr-1;i++){
    printf("{");
    for(j=0;j<nc-1;j++)
      printf("%10g,",matrix[i][j]);
    printf("%10g},\n",matrix[i][j]);
  }
  printf("{");
  for(j=0;j<nc;j++)
    printf("%10g,",matrix[i][j]);
  printf("%10g}}\n",matrix[i][j]);
}
/*------------------------------------------------------------------------*/
void psysvec(FILENAMES *filenames,int nconf,int natoms,int iunits,
	     double *px,double *py,double *pz,double *amass,
	     double *hmat,double **fcmat,double *dn,double **molmat,
	     int nspec,int *nmol,int *napm)
{  
  int i,j,l;
  double mass,mmass,fr;
  WORD name;
  FILE *fp;
#ifdef MOL_EIGEN
  int namol,namol3,ispec,ioff;
  double *dx,*dy,*dz,*r2,**axis,*d,*e;
  LINE line;
#endif
  
  mmass = 0.;
  for(i=0;i<natoms;i++){ mmass += amass[i];}

  for(i=0;i<3*natoms;i++){
#ifdef DEBUG
    printf("%d %.10g\n",i,dn[i]);
#endif
    fr=dn[i];
    if (fr<0.0) fr= -sqrt(-fr);
    else fr=sqrt(fr);
    dn[i]=fr;
  }
  /*  convert eigenvalues in ps^-1 to cm-1 if you want before recording */
  if(iunits==3) for (i=0;i<3*natoms;i++) dn[i] *= FCONV;

  for(i=0;i<natoms*3;i++){
    sprintf(name,"%s%.3d",filenames->sysvec,i);
    if((fp=fopen(name,"a+"))==NULL){
      fprintf(stderr,"ERROR: can't open %s for writing\n",name);
      exit(1);
    }
    fprintf(fp,"# %d %g\n",natoms,dn[i]);
    for(l=j=0;j<natoms;j++,l+=3){
      mass = (double)natoms/sqrt(amass[j]);
      fprintf(fp,"%9.5g %9.5g %9.5g %9.5g %9.5g %9.5g\n",px[j],py[j],pz[j],
	      fcmat[l  ][i]*mass,fcmat[l+1][i]*mass,fcmat[l+2][i]*mass);
    }
    fclose(fp);
  }

#ifdef MOL_EIGEN
  int namol,namol3,ispec,ioff;
  double *dx,*dy,*dz,*r2,**axis,*d,*e;

  ioff = 0;
  for(ispec=0;ispec<nspec;ispec++){
    namol = napm[ispec];
    namol3 = namol*3;

    sprintf(line,"spec %d %d napm = %d\n",ispec,nmol[ispec],namol);
    md_stdout(line);

    axis = dmatrix(0,(int)DIM-1,0,(int)DIM-1);
    dx = (double *)cmalloc(namol*sizeof(double));
    dy = (double *)cmalloc(namol*sizeof(double));
    dz = (double *)cmalloc(namol*sizeof(double));
    r2 = (double *)cmalloc(namol*sizeof(double));
    d = (double *)cmalloc(namol3*sizeof(double));
    e = (double *)cmalloc(namol3*sizeof(double));
    
    for(imol=0;imol<nmol[ispec];imol++){      
      cx = cy = cz = 0.;
      mmass = 0.;
      for(j=ioff,i=0;i<namol;i++,j++){
        mmass += mass = amass[j];
        cx += px[j]*mass; cy += py[j]*mass; cz += pz[j]*mass;
      }
      cx /= mmass; cy /= mmass; cz /= mmass;
      
      for(i=0,j=ioff;i<namol;i++,j++){
        dx[i] = px[j] - cx; dy[i] = py[j] - cy; dz[i] = pz[j] - cz;
        r2[i] = dx[i]*dx[i] + dy[i]*dy[i] + dz[i]*dz[i];
      }

      for(i=0;i<(int)DIM;i++) for(j=0;j<(int)DIM;j++) axis[i][j] = 0.;
      
      for(i=0;i<namol;i++){
        mass = amass[i+ioff];
        axis[0][0] += mass*(r2[i] - dx[i]*dx[i]);
        axis[1][1] += mass*(r2[i] - dy[i]*dy[i]);
        axis[2][2] += mass*(r2[i] - dz[i]*dz[i]);
        
        axis[0][1] -= mass*dx[i]*dy[i];
        axis[0][2] -= mass*dx[i]*dz[i];
        axis[1][2] -= mass*dy[i]*dz[i];
      }
      
      axis[1][0] = axis[0][1];
      axis[2][0] = axis[0][2];
      axis[2][1] = axis[1][2];
      
      /* solve for eigenvectors which are the priciple axis*/
      
      if(namol>1){
        ctred2v(axis,(int)DIM,d,e);
        ctqliv(d,e,(int)DIM,axis);
        ceigsrtv(d,axis,(int)DIM);
      }

      sprintf(name,"%s%d.%d.axis",filenames->molvec,ispec,imol);
      fp=cfopenw(name);
      fprintf(fp,"3 1 1\n");
      fprintf(fp,"%8.5g %8.5g %8.5g %8.5g %8.5g %8.5g\n",cx,cy,cz,
	      axis[0][0],axis[1][0],axis[2][0]);
      fprintf(fp,"%8.5g %8.5g %8.5g %8.5g %8.5g %8.5g\n",cx,cy,cz,
	      axis[0][1],axis[1][1],axis[2][1]);
      fprintf(fp,"%8.5g %8.5g %8.5g %8.5g %8.5g %8.5g\n",cx,cy,cz,
	      axis[0][2],axis[1][2],axis[2][2]);
      fclose(fp);

      if(namol>1){
	/* get eigenvectors of isolated molecule */
        k = ioff*DIM;
        ctred2v(&(molmat[k]),namol3,d,e);
        ctqliv(d,e,namol3,&(molmat[k]));
        ceigsrtv(d,&(molmat[k]),namol3);
        
        for(i=0;i<namol3;i++){
          fr=d[i];
          if (fr<0.0) fr= -sqrt(-fr);
          else fr=sqrt(fr);
          d[i]=fr;
        }
	/*  convert eigenvalues in ps^-1 to cm-1 if you want */
        if(iunits==3) for (i=0;i<namol3;i++) d[i] *= FCONV;
        
        for(i=0;i<namol3;i++){
          sprintf(name,"%s%d_%d_%d_%g",filenames->molvec,ispec,imol,i,d[i]);
          fp=cfopenw(name);
          fprintf(fp,"%d 1 1\n",namol);
          for(l=j=0;j<namol;j++,l+=3){
            mass = (double)namol/sqrt(amass[ioff+j]);
            fprintf(fp,"%9.5g %9.5g %9.5g %9.5g %9.5g %9.5g\n",
               px[ioff+j],py[ioff+j],pz[ioff+j],
               molmat[k+l  ][i]*mass,
               molmat[k+l+1][i]*mass,
               molmat[k+l+2][i]*mass);
          }
          fclose(fp);
        }
      }
      ioff += namol; /* update the atom offset */
    }	/* loop over next molecule */
    /* free up temparary vectors */
    free_dmatrix(axis,0,(int)DIM-1,0,(int)DIM-1);
    free(dx); free(dy); free(dz); free(r2); free(d); free(e);
  } /* loop over next species */
#endif
}
/*---------------------------------------------------------------------*/
void store_freq(int nmodes,double *dn,double time)
{
  int i,zero[3],np;
  char name[50];
  double fq,fr;
  FILE *fp;
  LINE line;


  zero[0] = zero[1] = zero[2] = -1;
#ifdef ZERO
  fr = DBL_MAX;
  for(i=0;i<nmodes;i++){
    if(fabs(dn[i]) < fr){zero[0] = i;  fr = fabs(dn[i]); }
  }
  fr = DBL_MAX;
  for(i=0;i<nmodes;i++){
    if(fabs(dn[i]) < fr && i != zero[0]){
      zero[1] = i;  
      fr = fabs(dn[i]);
    }
  }
  fr = DBL_MAX;
  for(i=0;i<nmodes;i++){
    if(fabs(dn[i])<fr&&i!=zero[0]&&i!=zero[1]){
      zero[2] = i;
      fr=fabs(dn[i]);
    }
  }
  if(dn[zero[0]]>ERRMAX || dn[zero[1]] > ERRMAX || 
     dn[zero[2]] > ERRMAX){
    sprintf(line,"zero frequency exceeds tolerence of %g",ERRMAX);
    md_warning(line);
    sprintf(line,"%d %d %d %g %g %g",zero[0],zero[1],zero[2],
	    dn[zero[0]],dn[zero[1]],dn[zero[2]]);
    md_warning(line);
  }
#endif
  
  np = 0;
  for(i=0;i<nmodes;i++){
    if(i != zero[0] && i != zero[1] && i != zero[2]){
      sprintf(name,"freq/w.%d",np++);
      if((fp = fopen(name,"a+"))==NULL){
	fprintf(stderr,"ERROR: can't append file %s\n",name);
	exit(1);
      }
      if(dn[i]<0) fq = -sqrt(-dn[i]);
      else fq = sqrt(dn[i]);
      
      fprintf(fp,"%g %.15g\n",time,fq);
      fclose(fp);
    }
  }
  if(np != nmodes-3) fprintf(stderr,"MODES Don't match %d %d\n",np,nmodes-3);
}

/*---------------------------------------------------------------------*/
/* quantum weight the imaginary modes so that detailed balence is obeyed 
   c_q(w) = exp(beta hbar w/2) c_c(w)
 */
/* void qweighti(int npoints,double df,double frmin,double a1[],double a2[],
   double a3[],double temp,int iunits) */
/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* added by HABS HA and BS 3/30/98  Heather Ahlborn  */
/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
void qweighti(int npoints,double df,double frmin,double a1[],double b1[],
   double a2[],double a3[],double a4[],double a5[],double temp,int iunits)
{
  int i,j,k,n0;
  double fr,arg,fac,ufac=1.;

  if(iunits==3) ufac /= FCONV;
/*  changed to hbar from h b/c freq is omega in ps-1*  HABS 4/29/98 */
/*  arg = HBAR*2.*M_PI/(BOLTS*temp); /* because cm^-1 is nu not omega */
  arg = HBAR/(BOLTS*temp);
  
  for(i=0;i<npoints;i++){
    fr = i*df+frmin;
    if(fr < 0.) a2[i] = -fr;
    if(fr>=0.) {
      n0 = i;
      break;
    }
  }
  if(n0!=i) md_error("All negative frequencies! Can't get spectrum");
  
  /* add correct quantum weight from negative
     to postive spectum and zero negative portion */

  for(i=0;i<n0 && i<(npoints-n0);i++){
    j = n0-i-1;
    k = n0+i;
    fr = fabs(a2[j]*ufac);
/* factor of two was wrong in line below on quantum correction, fixed 4/29/98 HABS */
/*    fac = sinh(fr*arg)/(arg*fr); */
    fac = sinh(fr*arg/2.)/(arg*fr/2.);
    a2[k] = fac*a1[j]+a1[k];
/* added by HABS HA and BS 3/30/98    */
/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
    a5[k] = fac*b1[j]+b1[k];
    a3[k] = a1[j]+a1[k];
    a4[k] = b1[j]+b1[k];
/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
    a2[j] = a3[j] = 0.;
/* added by HABS HA and BS 3/30/98    */
/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
    a4[j] = a5[j] = 0.;
  }
  for(i=2*n0;i<npoints;i++) a3[i] = a2[i] = a1[i];
/* added by HABS HA and BS 3/30/98    */
/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
  for(i=2*n0;i<npoints;i++) a4[i] = a5[i] = b1[i];
}

/*
      do i=1,450
         read(22,*)r1,r2,r3,r4,r5,r6,r7
ccc         if(r1.lt.0.)then
            temper=236.
            arg=6.626e-22/2./1.38e-23/temper
            fac=sinh(abs(r1/5.3)*arg)/arg/abs(r1/5.3)
c qid imaginary is total qid - real qid
            r8=(r5-r2)
            r9=fac*(r5-r2)
            r10=r2+r9
            write(23,78)r1,r2,r8,r5,r9,r10
            write(24,*)r1,fac
ccc         endif
      enddo
      */
