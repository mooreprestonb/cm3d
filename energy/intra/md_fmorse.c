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

From space@space1.chemistry.duq.edu  Tue Oct 17 17:08:48 1995
Posted-Date: Tue, 17 Oct 1995 17:08:46 -0400
Received-Date: Tue, 17 Oct 1995 17:08:46 -0400
Received: (from space@localhost) by space1.chemistry.duq.edu (940816.SGI.8.6.9/8.6.9) id RAA02581 for moore@sg10.chem.upenn.edu; Tue, 17 Oct 1995 17:08:54 -0400
From: "Brian Space" <space@space1.chemistry.duq.edu>
Message-Id: <9510171708.ZM2579@space1.chemistry.duq.edu>
Date: Tue, 17 Oct 1995 17:08:54 -0400
In-Reply-To: "Preston B. Moore" <moore@sg10.chem.upenn.edu>
        "Re: code" (Oct 17,  3:42pm)
References: <9510111248.ZM12602@sg10.chem.upenn.edu> 
	<9510111321.ZM2549@space1.chemistry.duq.edu> 
	<9510111418.ZM13095@sg10.chem.upenn.edu> 
	<9510111412.ZM2771@space1.chemistry.duq.edu> 
	<9510111612.ZM13803@sg10.chem.upenn.edu> 
	<9510111617.ZM3276@space1.chemistry.duq.edu> 
	<9510111630.ZM14362@sg10.chem.upenn.edu> 
	<9510121121.ZM7435@space1.chemistry.duq.edu> 
	<9510121220.ZM26208@sg10.chem.upenn.edu> 
	<9510161154.ZM19568@space1.chemistry.duq.edu> 
	<9510161210.ZM25668@sg10.chem.upenn.edu> 
	<9510171539.ZM2258@space1.chemistry.duq.edu> 
	<9510171542.ZM10528@sg10.chem.upenn.edu>
X-Mailer: Z-Mail (3.2.0 26oct94 MediaMail)
To: "Preston B. Moore" <moore@sg10.chem.upenn.edu>
Subject: Re: code
Mime-Version: 1.0
Content-Type: multipart/mixed;
	boundary="PART-BOUNDARY=.19510171708.ZM2579.chemistry.duq.edu"
Status: ORr

--
--PART-BOUNDARY=.19510171708.ZM2579.chemistry.duq.edu
Content-Type: text/plain; charset=us-ascii

Hi P i am sending you the new md_fbond

 I recreated all the routines for morse expcept the numerical routine,
I comment all the relecvant stuff where i made changes
with bspace6 flags and I put important info at the  top !!!

It complies except for not knowing the 2 new aprameters
for the morse potential !


GO MORSE
Brian


-- 
Brian Space
Professor Brian Space
Department of Chemistry and Biochemistry
Duquesne University
Pittsburgh, PA 15282-1503

space@space1.chemistry.duq.edu
phone:(412)396-4732
fax:  (412)396-5683

--PART-BOUNDARY=.19510171708.ZM2579.chemistry.duq.edu
X-Zm-Content-Name: md_fbond.c
Content-Description: Text
Content-Type: text/plain ; name="md_fbond.c" ; charset=us-ascii

/* bspace6 for our morse changes we can rember that our harmonic frequency is 2.*D_morse*rho_morse*rho_morse and i wrote no numerical routineeeeeeeee and i don't have the params (D_morse and rho_morse)
passed in */
/* subroutines dealing with bonded interactions
   search_bond_base - searches the data base and sets up initial vectors
   fbond - gets the force of the bond
   getvbond - gets the potential energy of the bond
   bonded   - check to see if to atoms are bonded
*/

#include "md.h" 

int nbonds,*idxbond,*ibond,*jbond;
double *eqbond,*fkbond;

int ibon,jbon;
double eqbon,fkbon;

/*-----------------------------------------------------------------------*/
int bond_pair(int iatom,int jatom)
{
  int i,j,itemp,nb;

  if(iatom>jatom){ itemp = iatom;iatom = jatom; jatom = itemp;}

  for(nb=0;nb<nbonds;nb++){
    i = ibond[nb]; j = jbond[nb];
    if(i>j){itemp = i; i = j; j = itemp; }
    if(iatom == i && jatom == j) return 1;
  }
  return 0;
}
    
/*-----------------------------------------------------------------------*/
void exclude_bonds(int *nexclude,int **exclude)
{
  int i,j,nb,ne,nt;

  for(nb=0;nb<nbonds;nb++){
    if(ibond[nb]<jbond[nb]){i=ibond[nb];j=jbond[nb];
    } else {j=ibond[nb];i=jbond[nb];}
    ne=0;
    while(j>exclude[i][ne] && exclude[i][ne] != -1) ne++;
    if(exclude[i][ne] != j){
      nexclude[i]++;
      for(nt=nexclude[i]-1;nt>ne;nt--){
	exclude[i][nt] = exclude[i][nt-1];
      }
      exclude[i][ne] = j;
    }
  }
}

/*-----------------------------------------------------------------------*/

void fbond(double *px,double *py,double *pz,double *fx,double *fy,double *fz)
{
  int i,j,itype,nb;
  double dx,dy,dz,r,fr,dr;

  for(nb=0;nb<nbonds;nb++){
    i = ibond[nb]; j = jbond[nb]; itype = idxbond[nb];
    dx = px[i] - px[j];  dy = py[i] - py[j];  dz = pz[i] - pz[j];
    r = sqrt(dx*dx+dy*dy+dz*dz);
    dr = r-eqbond[itype];
    fr = fkbond[itype]*dr/r;

    dx *= fr; dy *= fr; dz *= fr;
    fx[i] -= dx;  fy[i] -= dy;  fz[i] -= dz;
    fx[j] += dx;  fy[j] += dy;  fz[j] += dz;
  }
}
/*-----------------------------------------------------------------------*/
/* bspace6 */
void morsefbond(double *px,double *py,double *pz,double *fx,double *fy,double *fz)
{
  int i,j,itype,nb;
  double dx,dy,dz,r,fr,dr,mpot;;

  for(nb=0;nb<nbonds;nb++){
    i = ibond[nb]; j = jbond[nb]; itype = idxbond[nb];
    dx = px[i] - px[j];  dy = py[i] - py[j];  dz = pz[i] - pz[j];
    r = sqrt(dx*dx+dy*dy+dz*dz);
    dr = r-eqbond[itype];
    mpot= 1.-exp(-rho_morse*dr);
    /* old harm     fr = fkbond[itype]*dr/r; */
    /* bspace6 fr is 1/r (dv/dr) for morse */
    fr = (2.*D_morse*rho_morse*mpot*(1.-mpot))/r;
    dx *= fr; dy *= fr; dz *= fr;
    fx[i] -= dx;  fy[i] -= dy;  fz[i] -= dz;
    fx[j] += dx;  fy[j] += dy;  fz[j] += dz;
  }
}
/*-------------------------------------------------------------------*/
void morsegetvbond(double *px,double *py,double *pz,double *UI,double *WI)
{
  int i,j,itype,nb;
  double dx,dy,dz,r,dr,spot=0.,sprs=0.,mpot;

  for(nb=0;nb<nbonds;nb++){
    i = ibond[nb];   j = jbond[nb];   itype = idxbond[nb];
    dx = px[i] - px[j];   dy = py[i] - py[j];   dz = pz[i] - pz[j];
    r = sqrt(dx*dx+dy*dy+dz*dz);
    dr = r - eqbond[itype];
    mpot= 1.-exp(-rho_morse*dr);
    spot+=D_morse*mpot*mpot;
    /*    spot += .5*fkbond[itype]*dr*dr;      sprs -= fkbond[itype]*dr*r; */
    /* bspace6 i am tabulating r * (dv/dr) */
    sprs -= r*(2.*D_morse*rho_morse*mpot*(1.-mpot));
  }
  *UI += spot;  *WI += sprs;
}
/*---------------------------------------------------------------*/
void getvbond(double *px,double *py,double *pz,double *UI,double *WI)
{
  int i,j,itype,nb;
  double dx,dy,dz,r,dr,spot=0.,sprs=0.;

  for(nb=0;nb<nbonds;nb++){
    i = ibond[nb];   j = jbond[nb];   itype = idxbond[nb];
    dx = px[i] - px[j];   dy = py[i] - py[j];   dz = pz[i] - pz[j];
    r = sqrt(dx*dx+dy*dy+dz*dz);
    dr = r - eqbond[itype];
    spot += .5*fkbond[itype]*dr*dr;
    sprs -= fkbond[itype]*dr*r;
  }
  *UI += spot;  *WI += sprs;
}
/*---------------------------------------------------------------*/
void fcbond(double *px,double *py,double *pz,double **fcmat)
{
  int i,j,itype,nb,ii,jj,kk;
  double dx,dy,dz,r,r3,dr,drx,dry,drz;
  double dxx,dyy,dzz,dxy,dxz,dyz,dudr,du2dr2;

  for(nb=0;nb<nbonds;nb++){
    i = ibond[nb];
    j = jbond[nb];
    itype = idxbond[nb];

    dx = px[i] - px[j];
    dy = py[i] - py[j];
    dz = pz[i] - pz[j];
    r = sqrt(dx*dx+dy*dy+dz*dz);

    dr = r-eqbond[itype];
    dudr = fkbond[itype]*dr;
    du2dr2 = fkbond[itype];

    r = 1./r;
    r3 = (r*r*r);
    drx=dx*r;
    dry=dy*r;
    drz=dz*r;
    dxx = dx*dx*r3-r;
    dyy = dy*dy*r3-r;
    dzz = dz*dz*r3-r;
    dxy = dx*dy*r3;
    dxz = dx*dz*r3;
    dyz = dy*dz*r3;

    dxx = du2dr2*drx*drx - dudr*dxx;
    dxy = du2dr2*drx*dry - dudr*dxy;
    dxz = du2dr2*drx*drz - dudr*dxz;
    dyy = du2dr2*dry*dry - dudr*dyy;
    dyz = du2dr2*dry*drz - dudr*dyz;
    dzz = du2dr2*drz*drz - dudr*dzz;

    ii=i*3-1;
    jj=j*3-1;

    /* fill in diagonal elements of matrix  */

    fcmat[ii+1][ii+1] += dxx;
    fcmat[ii+1][ii+2] += dxy;
    fcmat[ii+1][ii+3] += dxz;
    fcmat[ii+2][ii+1] += dxy;
    fcmat[ii+2][ii+2] += dyy;
    fcmat[ii+2][ii+3] += dyz;
    fcmat[ii+3][ii+1] += dxz;
    fcmat[ii+3][ii+2] += dyz;
    fcmat[ii+3][ii+3] += dzz;
    
    fcmat[jj+1][jj+1] += dxx;
    fcmat[jj+1][jj+2] += dxy;
    fcmat[jj+1][jj+3] += dxz;
    fcmat[jj+2][jj+1] += dxy;
    fcmat[jj+2][jj+2] += dyy;
    fcmat[jj+2][jj+3] += dyz;
    fcmat[jj+3][jj+1] += dxz;
    fcmat[jj+3][jj+2] += dyz;
    fcmat[jj+3][jj+3] += dzz;

      /*  Fill in off-diagonal elements  */

    if(ii>jj) {kk=ii;ii=jj;jj=kk;}
    fcmat[ii+1][jj+1] -= dxx;
    fcmat[ii+1][jj+2] -= dxy;
    fcmat[ii+1][jj+3] -= dxz;
    fcmat[ii+2][jj+1] -= dxy;
    fcmat[ii+2][jj+2] -= dyy;
    fcmat[ii+2][jj+3] -= dyz;
    fcmat[ii+3][jj+1] -= dxz;
    fcmat[ii+3][jj+2] -= dyz;
    fcmat[ii+3][jj+3] -= dzz;
  }
}

/*-----------------------------------------------------------------------*/
/*---------------------------------------------------------------*/
void morsefcbond(double *px,double *py,double *pz,double **fcmat)
{
  int i,j,itype,nb,ii,jj,kk;
  double dx,dy,dz,r,r3,dr,drx,dry,drz;
  double dxx,dyy,dzz,dxy,dxz,dyz,dudr,du2dr2,mpot;

  for(nb=0;nb<nbonds;nb++){
    i = ibond[nb];
    j = jbond[nb];
    itype = idxbond[nb];

    dx = px[i] - px[j];
    dy = py[i] - py[j];
    dz = pz[i] - pz[j];
    r = sqrt(dx*dx+dy*dy+dz*dz);

    dr = r-eqbond[itype];
    /*bspace6*/
    mpot= 1.-exp(-rho_morse*dr);
    dudr = (2.*D_morse*rho_morse*mpot*(1.-mpot));
    /*     dudr = fkbond[itype]*dr; */
    du2dr2 = 2.*D_morse*rho_morse*rho_morse*(1.-mpot)*(1.-2.*mpot);
    /*    du2dr2 = fkbond[itype]; */

    r = 1./r;
    r3 = (r*r*r);
    drx=dx*r;
    dry=dy*r;
    drz=dz*r;
    dxx = dx*dx*r3-r;
    dyy = dy*dy*r3-r;
    dzz = dz*dz*r3-r;
    dxy = dx*dy*r3;
    dxz = dx*dz*r3;
    dyz = dy*dz*r3;

    dxx = du2dr2*drx*drx - dudr*dxx;
    dxy = du2dr2*drx*dry - dudr*dxy;
    dxz = du2dr2*drx*drz - dudr*dxz;
    dyy = du2dr2*dry*dry - dudr*dyy;
    dyz = du2dr2*dry*drz - dudr*dyz;
    dzz = du2dr2*drz*drz - dudr*dzz;

    ii=i*3-1;
    jj=j*3-1;

    /* fill in diagonal elements of matrix  */

    fcmat[ii+1][ii+1] += dxx;
    fcmat[ii+1][ii+2] += dxy;
    fcmat[ii+1][ii+3] += dxz;
    fcmat[ii+2][ii+1] += dxy;
    fcmat[ii+2][ii+2] += dyy;
    fcmat[ii+2][ii+3] += dyz;
    fcmat[ii+3][ii+1] += dxz;
    fcmat[ii+3][ii+2] += dyz;
    fcmat[ii+3][ii+3] += dzz;
    
    fcmat[jj+1][jj+1] += dxx;
    fcmat[jj+1][jj+2] += dxy;
    fcmat[jj+1][jj+3] += dxz;
    fcmat[jj+2][jj+1] += dxy;
    fcmat[jj+2][jj+2] += dyy;
    fcmat[jj+2][jj+3] += dyz;
    fcmat[jj+3][jj+1] += dxz;
    fcmat[jj+3][jj+2] += dyz;
    fcmat[jj+3][jj+3] += dzz;

      /*  Fill in off-diagonal elements  */

    if(ii>jj) {kk=ii;ii=jj;jj=kk;}
    fcmat[ii+1][jj+1] -= dxx;
    fcmat[ii+1][jj+2] -= dxy;
    fcmat[ii+1][jj+3] -= dxz;
    fcmat[ii+2][jj+1] -= dxy;
    fcmat[ii+2][jj+2] -= dyy;
    fcmat[ii+2][jj+3] -= dyz;
    fcmat[ii+3][jj+1] -= dxz;
    fcmat[ii+3][jj+2] -= dyz;
    fcmat[ii+3][jj+3] -= dzz;
  }
}

/*-----------------------------------------------------------------------*/
void fcbondn(double *px,double *py,double *pz,double **fcmat)
{
  int nb,itype;

  for(nb=0;nb<nbonds;nb++){

    ibon = ibond[nb];
    jbon = jbond[nb];
    itype = idxbond[nb];

    eqbon = eqbond[itype];
    fkbon = fkbond[itype];

    ngradv2(px,py,pz,ibon,ibon,potbond,fcmat);
    ngradv2(px,py,pz,ibon,jbon,potbond,fcmat);
    ngradv2(px,py,pz,jbon,ibon,potbond,fcmat);
    ngradv2(px,py,pz,jbon,jbon,potbond,fcmat);
  }
}
/*----------------------------------------------------------------*/
double potbond(double *px,double *py,double *pz)
{
  double dx,dy,dz,r,dr,spot;

  dx = px[ibon] - px[jbon];
  dy = py[ibon] - py[jbon];
  dz = pz[ibon] - pz[jbon];
  r = sqrt(dx*dx+dy*dy+dz*dz);

  dr = r-eqbon;
  spot = .5*fkbon*dr*dr;

  return spot;
}
/*------------------------------------------------------------------*/
double morsepotbond(double *px,double *py,double *pz)
{
  double dx,dy,dz,r,dr,spot,mpot;

  dx = px[ibon] - px[jbon];
  dy = py[ibon] - py[jbon];
  dz = pz[ibon] - pz[jbon];
  r = sqrt(dx*dx+dy*dy+dz*dz);

  dr = r-eqbon;
    mpot= 1.-exp(-rho_morse*dr);
    spot+=D_morse*mpot*mpot;

  return spot;
}

--PART-BOUNDARY=.19510171708.ZM2579.chemistry.duq.edu--

