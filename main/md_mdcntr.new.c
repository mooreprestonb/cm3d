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

/* controling subroutine for molecular dynamcis 
   using nose thermostats (NVT) and extened system variables
   for the volume
*/

#include "md.h"

/* #define USER_FUNC */

/* might get different numbers from the eigensolver */
/* #define NO_ROTATE */
#define FREEZEFORCE
#ifdef FREEZEFORCE
void save_freezeforce(int,COORDS *coords);
#endif

/*! subroutine to see if we should change the variables of CM3D while running
 */
int checksigfile(void){
  WORD sigfile;
  FILE *fp;
  int pid = (int) getpid();
  sprintf(sigfile,".stop.%d",pid);

  if((fp = fopen(sigfile,"r"))!=NULL){
    md_warning("sigfile detected, stopping! (and deleting sigfile)");
    fclose(fp);
    remove(sigfile);
    return 1;
  }
  return 0;
}

/*! subroutine to check for existant of .print.pid
 */
int checksigprint(void){
  WORD sigfile;
  FILE *fp;
  int pid = (int) getpid();

  sprintf (sigfile,".print.%d",pid);

  if((fp = fopen(sigfile,"r"))!=NULL){
    md_warning("print sigfile detected: printing values"); 
    fclose(fp);
    remove(sigfile);
    return 1;
  }
  return 0;
}

/*------------------------------------------------------------------*/
void md_mdcntr_new(char *command,FILENAMES *filenames,SIMPARMS *simparms,
		   WRITE_STEP *write_step,COORDS *coords,INTER *inter,
		   NGBR *ngbr)
{
  int exitall=0,exitsum=0;
  double acttime,oldtime,difftime,starttime;
  ENERGIES energies;

  if(simparms->rank==0){
    readpos(command,filenames->initfile,simparms,
	    &energies.estep,&energies,coords);
    if(simparms->istart == 1){
      samvel(simparms,coords,1);
      samveta(simparms,coords);
      if(simparms->ivol) samvbar(simparms,coords);
    }
    if(simparms->iperd==0 && simparms->istart!=3){
#ifndef NO_ROTATE 
      rotate(simparms->natoms,coords->px,coords->py,coords->pz,
	     coords->amass);
#endif
      zero_angm(simparms->natoms,coords->px,coords->py,coords->pz,
         coords->vx,coords->vy,coords->vz,coords->amass);

      /* rescale velocities */
      if(simparms->istart==1) samvel(simparms,coords,0);
    }
    
    zero_averages(simparms,&energies);
  }
  /* broadcast possitions and velocities, masses and charge 
     use the fact that they are contiguous in memory */
#ifdef PARA
  bcast_system(simparms,coords);
#endif
  
  gethinv9(coords->hmat,coords->hmati);

  if(!simparms->ivol){
    if(simparms->nbar > 1){
      md_error("you can't barostat a constant volume!\n \tChange either ivol or nbar in the control file.\n");
      exit(1);
    }
  }

  /* setup integrator and get energies */
  energies.acpu = 0.;
  energies.acpu = cputime();
  step_init(simparms,coords,&energies,inter,ngbr);

  if(simparms->iunits==1){
    getengn_cheese(simparms,coords,&energies.zken,&energies.potn);
    getengv_cheese(simparms,coords,&energies.zkev,&energies.potv,
		   &energies.zken);
  } else {
    getengn(simparms,coords,&energies.zken,&energies.potn);
    getengv(simparms,coords,&energies.zkev,&energies.potv);
  }

  /* output stuff the user might want to know */
  if(simparms->rank==0){
    output_simparm(simparms,write_step,simparms->istep);
    output_initval(simparms,coords,ngbr,&energies);
        
    set_econi(&energies);
  }

#ifdef USER_FUNC
  md_stdout("User Function defined");
  user_initial(simparms->istep,simparms->natoms,simparms->dt,
	       coords->px,coords->py,coords->pz,
	       coords->vx,coords->vy,coords->vz,
	       coords->amass,coords->qch,coords->hmat);
#endif

  /* initial energies */
  energies.acpu = cputime();
  simparms->fstep=simparms->istep;
  if(simparms->istart == 3) ++simparms->istep;
  if(simparms->rank==0){
    writescreen(simparms,ngbr,&energies,coords->hmat);

    if(check_vals(9,energies.potra,energies.poter,energies.zke,
		  energies.potn,energies.zken,energies.potv,energies.zkev,
		  energies.tiout,energies.prsi)){
      free_all_pointers(simparms,coords,inter,ngbr);
      exit(1);
    }

    /* **************  Commence MD Steps ******************** */
    md_stdout("\n ***** Begining Molecular Dynamics Steps *****");
  }
  energies.acpu += cputime();
  oldtime=starttime=realtime();

  for(;simparms->istep<=simparms->nstep;++simparms->istep){

    step(simparms,coords,&energies,inter,ngbr);
    geteng(simparms,coords,&energies);
    gettngbr(ngbr);

    if(simparms->iunits==1){
      getengn_cheese(simparms,coords,&energies.zken,&energies.potn);
      getengv_cheese(simparms,coords,&energies.zkev,&energies.potv,
		     &energies.zken);
    } else {
      getengn(simparms,coords,&energies.zken,&energies.potn);
      getengv(simparms,coords,&energies.zkev,&energies.potv);
    }
    
    if(simparms->rank==0){
      accumulate(coords->hmat,&energies);
      
      /* if .print.pid exist write to screen */
      if(checksigprint()){writescreen(simparms,ngbr,&energies,coords->hmat);}
  

      /* dump restart info and output averages every ndump passes */
      if(write_step->nscrn && !(simparms->istep%(write_step->nscrn))){
        energies.acpu += cputime();
        writescreen(simparms,ngbr,&energies,coords->hmat);
      }
      if(write_step->ndump && !(simparms->istep%(write_step->ndump))){
        writestart(filenames,simparms->istep,simparms,&energies,coords);
      }
      if(write_step->rescalevel&&!(simparms->istep%(write_step->rescalevel))){
        samvel(simparms,coords,0);
        if(simparms->iperd==0){
          zero_angm(simparms->natoms,coords->px,coords->py,coords->pz,
             coords->vx,coords->vy,coords->vz,coords->amass);
        }
      }
      if(write_step->resamvel && !(simparms->istep%(write_step->resamvel))){
        samvel(simparms,coords,1);
        if(simparms->iperd==0){
          zero_angm(simparms->natoms,coords->px,coords->py,coords->pz,
             coords->vx,coords->vy,coords->vz,coords->amass);
        }
        samveta(simparms,coords);
        if(simparms->ivol) samvbar(simparms,coords); 
      }
      /* -------  record instantaneous properties ------- */
      
      if(write_step->nconf && simparms->istep%(write_step->nconf) == 0){
        save_conf(filenames->fconf,simparms->natoms,coords);
      }
      if(write_step->nvel && simparms->istep%(write_step->nvel) == 0){
        save_vel(filenames->fvel,simparms->natoms,coords);
      }
      
      if(write_step->nforce && simparms->istep%(write_step->nforce) == 0){
        save_force(filenames->fforce,simparms->natoms,coords);
      }
      
      if(write_step->ninst && simparms->istep%(write_step->ninst) == 0){
        save_inst(filenames->fham,filenames->feng,filenames->fext,
           simparms->istep,simparms->dt,&energies,coords->hmat);
      }
#ifdef FREEZEFORCE
      save_freezeforce(simparms->nfreeze,coords);
#endif
      
#ifdef USER_FUNC
      user_function(simparms->istep,simparms->natoms,simparms->dt,
         coords->px,coords->py,coords->pz,
         coords->vx,coords->vy,coords->vz,
         coords->amass,coords->qch,coords->hmat);
#endif 
    }
    /* calculate how much real time has passed and see if we need to exit */
    acttime=realtime();
    difftime=acttime-oldtime; if (difftime<0) difftime= -difftime;
    oldtime=acttime;
#ifdef DEBUG
    printf("here comes wallclock and acttime %g,%g,%g %g %g \n",
       simparms->wallclock,difftime,acttime-starttime,starttime,acttime);
#endif
    exitsum += checksigfile();
    if ((acttime-starttime+difftime) > simparms->wallclock) ++exitall;
#ifdef PARA
    if(exitall) fprintf(stdout,"Rank %d has exceeded the time limit %g\n",
       simparms->rank,simparms->wallclock);
    MPI_Allreduce(&exitall,&exitsum,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    bcast_system(simparms,coords);
#else
    exitsum += exitall;
#endif
    if (exitsum) break;
  } /* for istep loop */
  
  if(simparms->rank==0){
    if (exitsum) md_warning("Kicked out by timelimit or signal file!");
    md_stdout("\n *** Finnished Molecular Dynamics Steps ***\n");
    if (!exitall) simparms->istep -= 1;
    finish(filenames,simparms,&energies,coords);
  }
}

/*---------------------------------------------------------------------*/
#ifdef PARA
void bcast_system_new(SIMPARMS *simparms,COORDS *coords)
{
  int n;

  MPI_Bcast(&simparms->istep,1,MPI_INT,0,MPI_COMM_WORLD);

  MPI_Bcast(coords->px,8*simparms->natoms,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(coords->hmat,9,MPI_DOUBLE,0,MPI_COMM_WORLD);

  n = (simparms->ntherm+1)*(simparms->nchain+1);
  MPI_Bcast(&coords->eta[0][0],n,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&coords->veta[0][0],n,MPI_DOUBLE,0,MPI_COMM_WORLD);

  MPI_Bcast(&coords->pvol,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&coords->vvol,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(coords->vvol9,9,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&coords->mvol,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

  if(simparms->nbar!=0){
    MPI_Bcast(coords->bar,simparms->nbar,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(coords->vbar,simparms->nbar,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(coords->mbar,simparms->nbar,MPI_DOUBLE,0,MPI_COMM_WORLD);
  }
}
#endif
/*---------------------------------------------------------------------*/

void gettngbr_new(NGBR *ngbr)
{
#ifdef PARA
  MPI_Reduce(&ngbr->nprs,&ngbr->nprst,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&ngbr->nprse,&ngbr->nprset,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&ngbr->nprt,&ngbr->nprtt,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&ngbr->nprte,&ngbr->nprtet,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
#else
  ngbr->nprst  = ngbr->nprs;
  ngbr->nprset = ngbr->nprse;
  ngbr->nprtt  = ngbr->nprt;
  ngbr->nprtet = ngbr->nprte;
#endif
}
/*---------------------------------------------------------------------*/
