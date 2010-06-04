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

/* routine that are used in subspace dynamics */

#include "md.h"

/* #define NUMERIC  */
/* #define WRITE_FCMAT */
/* #define PRINT_TRAN */
/* #define IDENTITY */
/* #define DEBUG */

/*-----------------------------------------------------------------*/
void get_subspace(SIMPARMS *simparms,COORDS *coords,INTER *inter,
		  double *x2,double *y2,double *z2,SUBSPACE *subspace)
{
  int natoms;
  double *vx,*vy,*vz,*amass;
#ifdef DEBUG
  int i,j,k;
  double xke,yke,zke;
#endif  

  natoms = simparms->natoms;
  amass = coords->amass;
  vx = coords->vx; vy = coords->vy; vz = coords->vz;

#ifdef IDENTITY
  /* allocate memory for subspace */
  subspace->nstate = natoms*3;
  subspace->d2v = dmatrix(0,natoms*3-1,0,subspace->nstate-1);
  for(i=0;i<natoms*3;i++){
    for(j=i+1;j<natoms*3;j++){
      subspace->d2v[i][j] = subspace->d2v[j][i] = 0.;
    }
    subspace->d2v[i][i] = 1.;
  }
#else
  get_subvec(simparms,coords,inter,subspace);
#endif

#ifdef DEBUG
#ifdef PRINT_TRAN
  for(i=0;i<subspace->nstate;i++){
    printf("%4d",i);
    for(j=0;j<natoms*3;j++){
      printf(" %9g",(subspace->d2v)[j][i]);
    }
    printf("\n");
  }
#endif
  /* calculate pT p */
  yke = zke = 0.;
  for(i=0;i<subspace->nstate;i++){
    for(j=0;j<subspace->nstate;j++){
      xke = 0.;
      for(k=0;k<natoms*3;k++){
	xke += (subspace->d2v)[k][i]*(subspace->d2v)[k][j];
      }
      if(i == j){
	if(fabs(xke-1.) > ERRMAX) printf("%d %d %g != 1\n",i,j,xke);
	yke += xke;
      } else  {
	if(fabs(xke)>ERRMAX) printf("%d %d %g != 0\n",i,j,xke);
	zke += fabs(xke);
      }
    }
  }
  printf("trace = %g  sum off = %g\n",yke,zke);
#endif

  /* Set up the reduced coord velocities */
  mass_weight(natoms,vx,vy,vz,amass,1);
  projectxtoz(vx,vy,vz,subspace->q1,subspace->d2v,subspace->nstate,natoms);
  expandztox(vx,vy,vz,subspace->q1,subspace->d2v,subspace->nstate,natoms);
  mass_weight(natoms,vx,vy,vz,amass,-1);

  /* Set up the reduced coord accelerations */
  mass_weight(natoms,x2,y2,z2,coords->amass,1);
  projectxtoz(x2,y2,z2,subspace->q2,subspace->d2v,subspace->nstate,natoms);
  expandztox(x2,y2,z2,subspace->q2,subspace->d2v,subspace->nstate,natoms);
  mass_weight(natoms,x2,y2,z2,coords->amass,-1);

  /*  Set up the reduced coords positions
      !!!!Defintion makes the intial q0's=0!!!! */

  /* q0=0 is the initial condition of the subspace positions */
  azero(subspace->q0,subspace->nstate);
  /* unproject positions (which just zero's position displacements) 
     this is silly since we get them back in SO DON'T DO IT */

#ifdef PRINT_TRAN
  for(i=0;i<subspace->nstate;i++){
    printf("transform %d %g %g %g\n",i,
	   subspace->q0[i],subspace->q1[i],subspace->q2[i]);
  }
#endif
}
/*-----------------------------------------------------------------*/
void get_subvec(SIMPARMS *simparms,COORDS *coords,INTER *inter,
		SUBSPACE *subspace)
{
  int i,isubst,nr,ni,nt;
  int kdiag,*indx,nagp,nsubvec,ngrp_min,ngrp_max,ngroup;
  double *dn,*en,**fcmat,**fcgp;
  double freq_max,freq_min;
  LINE line;
  
  nsubvec = 0;
  ngroup = simparms->igroup[simparms->natoms-1] + 1;

  /* count subspaces */
  for(isubst=0;isubst<subspace->num_subst;isubst++){
    ngrp_min = isubst*ngroup/subspace->num_subst;
    ngrp_max = (isubst+1)*ngroup/subspace->num_subst;
    if(subspace->num_vecsub_max[isubst] == -1){
      sprintf(line,"Changing the number of vectors in substructure %d",isubst);
      md_warning(line);
      subspace->num_vecsub[isubst] = 0;
      for(i=0;i<simparms->natoms;i++){
	if(simparms->igroup[i] >= ngrp_min && simparms->igroup[i] < ngrp_max){
	  subspace->num_vecsub[isubst] += 3;  /* three dof/atom */
	  nsubvec += 3;
	}
      }
    } else {
      nsubvec += subspace->num_vecsub_max[isubst];
      subspace->num_vecsub[isubst] = subspace->num_vecsub_max[isubst];
    }
  }
  if(subspace->nstate_max == 0) {
    sprintf(line,"Changing nstate to be sum of nsubvec = %d",nsubvec);
    md_warning(line);
    subspace->nstate = nsubvec;
  }

  /* allocate memory for subspace */
  subspace->d2v = dmatrix(0,simparms->natoms*3-1,0,subspace->nstate-1);
  subspace->d2vn = (double *)cmalloc(simparms->natoms*3*sizeof(double));
  simparms->mem_bytes += simparms->natoms*3*subspace->nstate + simparms->natoms*3;

  if(subspace->igetvec == -1){ /* get subspace do not read them in */

    /* temp set up for diagonalize variables */
    freq_max = subspace->freq_max;
    freq_min = subspace->freq_min;
    nsubvec = kdiag = 0;
    
    /* allocate temporary memory */
    indx = (int *)cmalloc(simparms->natoms*sizeof(int));
    en = (double *)cmalloc(simparms->natoms*3*sizeof(double));
    dn = (double *)cmalloc(simparms->natoms*3*sizeof(double));
    fcmat = dmatrix(0,simparms->natoms*3-1,0,simparms->natoms*3-1);
    
    get_fcmat(simparms,coords,inter,fcmat);
    
#ifdef WRITE_FCMAT
    md_stdout("****   FCMAT  ******");
    for(i=0;i<simparms->natoms*3;i++){
      for(j=0;j<simparms->natoms*3;j++){
	printf("%d %d %g\n",j,i,fcmat[j][i]);
      }
    }
#endif

    /* diagonalize each subspace */
    for(isubst=0;isubst<subspace->num_subst;isubst++){
      /* pick out matrix of forces */

      get_grp_fcmat(isubst,indx,simparms->natoms*3,simparms->igroup,
		    fcmat,subspace,&nagp,&fcgp);

#ifdef WRITE_FCMAT
      md_stdout("****   FCGP  ******");
      for(i=0;i<nagp*3;i++){
	for(j=0;j<nagp*3;j++){
	  printf("%d %d %g\n",j,i,fcgp[j][i]);
	}
      }
#endif

      /* diagonalize */
      sprintf(line,"Diagonalizeing sub block of force constant matrix %dx%d",
	      nagp*3,nagp*3);
      md_stdout(line);
      diagonalize(fcgp,nagp*3,dn,freq_min,freq_max,
         subspace->num_real,subspace->num_imag,kdiag);
      
#ifdef WRITE_FCMAT
      md_stdout("****   EIGENVALUES  ******");
      for(i=0;i<nagp*3;i++){
        printf("%d %g\n",i,dn[i]);
      }
      md_stdout("****   EIGENVECTORS  ******");
      for(i=0;i<nagp*3;i++){
        for(j=0;j<nagp*3;j++){
          printf("%d %d %g\n",j,i,fcgp[j][i]);
        }
      }
#endif

      pack(subspace,isubst,nagp*3,dn,fcgp);
      
      pack_d2v(isubst,subspace,nagp,subspace->d2v,subspace->d2vn,
         fcgp,dn,indx,nsubvec);
      
      free_dmatrix(fcgp,0,nagp*3-1,0,nagp*3-1);
      nsubvec += subspace->num_vecsub[isubst];
      nagp = 0;
    }    
    free_dmatrix(fcmat,0,simparms->natoms*3-1,0,simparms->natoms*3-1);
    free(dn);  free(en);  free(indx);
  } else if(subspace->igetvec == 1){
    dn = (double *)cmalloc(simparms->natoms*3*sizeof(double));
    read_eigvec(subspace->vecfile,subspace->nstate,simparms->natoms*3,dn,
       subspace->d2v);
    free(dn);
  } else {
    md_error("ERROR:iflgvec != 1 or -1");
  }
  /* save_eigvec(vec_file,subspace->nstate,simparms->natoms*3,dn,
     subspace->d2v); */
  
  if(nsubvec != subspace->nstate){
    sprintf(line,"Number of spacespaces %d != %d number of states!",
	    nsubvec,subspace->nstate);
    sprintf(line,"%s\n\tChanging number of states to %d",line,nsubvec);
    md_warning(line);
    subspace->d2v = drematrix(subspace->d2v,0,simparms->natoms*3-1,0,
			      subspace->nstate-1,0,simparms->natoms*3-1,0,
			      nsubvec-1);
    subspace->nstate = nsubvec;
  }


#ifdef WRITE_FCMAT
  md_stdout("****   TRANSFORM  ******");
  for(j=0;j<subspace->nstate;j++){
    for(i=0;i<simparms->natoms*3;i++){
      printf("%d %d %g\n",i,j,subspace->d2v[i][j]);
    }
  }
#endif
  vsrtasnd(subspace->d2vn,subspace->d2v,subspace->nstate,simparms->natoms*3);

  nr = ni = nt = 0;
  sprintf(line,"----------- %d subvectors --------------",subspace->num_subst);
  md_stdout(line);
  md_stdout("     #  Vectors  Real  Imaginary");
  for(i=0;i<subspace->num_subst;i++){
    sprintf(line,"%6d %6d %6d %6d",i+1,subspace->num_vecsub[i],
	    subspace->num_real[i],subspace->num_imag[i]);
    md_stdout(line);
    nr += subspace->num_real[i];
    ni += subspace->num_imag[i];
    nt += subspace->num_vecsub[i];
  }
  md_stdout("-------------------------------------------");
  sprintf(line,"total  %6d %6d %6d",nt,nr,ni);
  md_stdout(line);
}
/*--------------------------------------------------------------*/
void get_grp_fcmat(int isubr,int *indx,int natoms3,int *igroup,double **fcmat,
		   SUBSPACE *subspace,int *nagp,double ***fcgp)
{
  int i,j,ngroup,ngrp_min,ngrp_max,natoms;
  LINE line;

  natoms = natoms3/3;
  ngroup = igroup[natoms-1] + 1;
  
#ifdef DEBUG
  printf("ngroup = %d\n",ngroup);
#endif
  if(ngroup % subspace->num_subst !=0){
    sprintf(line,"%d substructures not divisable by %d groups",
	    subspace->num_subst,ngroup);
    md_error(line);
  }

  /* count atoms in this substruct */
  *nagp = 0;
  ngrp_min = isubr*ngroup/subspace->num_subst;
  ngrp_max = (isubr+1)*ngroup/subspace->num_subst;

  for(i=0;i<natoms;i++){
    if(igroup[i] >= ngrp_min && igroup[i] < ngrp_max){
      indx[*nagp] = i;
      ++(*nagp);
    }
  }

#ifdef DEBUG
  printf("Number of atoms in substructure %d is %d\n",isubr,*nagp);
#endif
  *fcgp = dmatrix(0,*nagp*3-1,0,*nagp*3-1);

  for(i=0;i<*nagp;i++){
    for(j=0;j<*nagp;j++){
      (*fcgp)[i*3  ][j*3  ] = fcmat[indx[i]*3  ][indx[j]*3  ];
      (*fcgp)[i*3  ][j*3+1] = fcmat[indx[i]*3  ][indx[j]*3+1];
      (*fcgp)[i*3  ][j*3+2] = fcmat[indx[i]*3  ][indx[j]*3+2];
      (*fcgp)[i*3+1][j*3  ] = fcmat[indx[i]*3+1][indx[j]*3  ];
      (*fcgp)[i*3+1][j*3+1] = fcmat[indx[i]*3+1][indx[j]*3+1];
      (*fcgp)[i*3+1][j*3+2] = fcmat[indx[i]*3+1][indx[j]*3+2];
      (*fcgp)[i*3+2][j*3  ] = fcmat[indx[i]*3+2][indx[j]*3  ];
      (*fcgp)[i*3+2][j*3+1] = fcmat[indx[i]*3+2][indx[j]*3+1];
      (*fcgp)[i*3+2][j*3+2] = fcmat[indx[i]*3+2][indx[j]*3+2];
    }
  }
}
/*----------------------------------------------------------------*/
/* pack the vector that you want into the lowest nvectors */
void pack(SUBSPACE *subspace,int isubst,int natoms3,double *dn,double **fcmat)
{
  int i,j,k,izero;
  double dtmp;
  LINE line;

  /* dn eigenvalues are sorted in decending order */
  i = 0;
  while(dn[i]>0. && i < natoms3){ /* find zero */
    i++;
  }
  izero = i;

  /* ------ set up num_real + num_imag  = num_subvec--------*/
  /* allow for maximum */
  if(subspace->num_vecsub_max[isubst] == -1){
    subspace->num_vecsub[isubst] = natoms3;
  }

  /* count up how many imaginary we want */
  if(subspace->num_imag_max[isubst] == -1){
    subspace->num_imag[isubst] = MIN(natoms3-izero,
				     subspace->num_vecsub[isubst]);
  } else {
    subspace->num_imag[isubst] = subspace->num_imag_max[isubst];
  }

  /* count up how many real we want */
  if(subspace->num_real_max[isubst] == -1){
    subspace->num_real[isubst] = MIN(izero,subspace->num_vecsub[isubst]
				     - subspace->num_imag[isubst]);
  } else {
    subspace->num_real[isubst] = subspace->num_real_max[isubst];
  }
 
  /* redo num_vecsub so that we have only the ones we want */
  subspace->num_vecsub[isubst] = subspace->num_real[isubst] + 
    subspace->num_imag[isubst];

  /* ------ consistancy checks ---------*/
  if(subspace->num_real[isubst]+subspace->num_imag[isubst] > natoms3){
    sprintf(line,"Number of eigenvectors (%d+%d) requested is > then the",
	    subspace->num_real[isubst],subspace->num_imag[isubst]);
    sprintf(line,"%s number of eigenvectors (%d)",line,natoms3);
    md_error(line);
  }
  if(subspace->num_real[isubst]>izero){
    sprintf(line,"Number of real eigenvectors (%d) requested is > then",
	    subspace->num_real[isubst]);
    sprintf(line,"%s number of real eigenvectors found (%d)",line,izero);
    md_error(line);
  }
  if(subspace->num_imag[isubst]>natoms3-izero){
    sprintf(line,"Number of imaginary eigenvectors (%d) requested is > then",
	    subspace->num_imag[isubst]);
    sprintf(line,"%s number of imaginary eigenvectors found (%d)",
	    line,natoms3-izero);
    md_error(line);
  }

  /* ------------- pack positive eigenvalues ------------*/ 
  for(i=0;i<MIN(subspace->num_real[isubst],izero/2);i++){
    k = izero - i-1;
    dtmp = dn[i];
    dn[i] = dn[k];
    dn[k]=dtmp;
    for (j=0;j<natoms3;j++){
      dtmp = fcmat[j][i];
      fcmat[j][i]=fcmat[j][k];
      fcmat[j][k]=dtmp;
    }
  }
  /* --------------- pack negative eigenvalues --------------------*/
  for(i=0;i<subspace->num_imag[isubst];i++){
    k = subspace->num_real[isubst]+i;
    dtmp  = dn[i+izero];
    dn[i+izero] = dn[k];
    dn[k]=dtmp;
    for (j=0;j<natoms3;j++){
      dtmp = fcmat[j][i+izero];
      fcmat[j][i+izero]=fcmat[j][k];
      fcmat[j][k]=dtmp;
    }
  }
  
#ifdef DEBUG
  printf("subspace group %d is %d %d %d\n",isubst,subspace->num_vecsub[isubst],
	 subspace->num_real[isubst],subspace->num_imag[isubst]);
#endif

}
/*------------------------------------------------------------------*/
/* unpack the subspace vectors into d2v */
void pack_d2v(int isubst,SUBSPACE *subspace,int nagp,double **d2v,
	      double *eigval,double **fcgp,double *dn,int *indx,int nsubvec)
{
  int i,j,k;
  double freq;
  LINE line;

  if(subspace->num_vecsub_max[isubst] == -1){
    subspace->num_vecsub[isubst] = subspace->num_real[isubst] + 
      subspace->num_imag[isubst];
  } else {
    subspace->num_vecsub[isubst] = subspace->num_vecsub_max[isubst];
  }
  
  if(nagp*3<subspace->num_vecsub[isubst]){
    sprintf(line,"Number of subspace vectors specified (%d) is greater then",
	    subspace->num_vecsub[isubst]);
    sprintf(line,"%s the number of vectors (%d) in substructure %d",
	    line,nagp*3,isubst);
    md_error(line);
  }
  for(i=0;i<subspace->num_vecsub[isubst];i++){
    freq = (dn[i]>=0)?sqrt(dn[i]):-sqrt(-dn[i]);
    eigval[i] = freq;
#ifdef DEBUG
    printf("%5d %5d %10g %5d %5d %5d\n",i,nsubvec,freq,nagp,
	   isubst,subspace->num_vecsub[isubst]);
#endif
    k = i+nsubvec;
    for(j=0;j<nagp;j++){
      if(freq < subspace->freq_min || freq > subspace->freq_max){
	md_error("Frequency out of range in subspace");
      }
      d2v[indx[j]*3  ][k] = fcgp[j*3  ][i];
      d2v[indx[j]*3+1][k] = fcgp[j*3+1][i];
      d2v[indx[j]*3+2][k] = fcgp[j*3+2][i];
    }
  }
}
/*----------------------------------------------------------------*/
