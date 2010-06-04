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

/* Main C program for simulations */

#include "md.h"

#include <sys/types.h>
#include <unistd.h>

void args_filename_new(FILENAMES *filenames,int argc,char *argv[])
{
  int i;

  strcpy(filenames->errfile,"");
  strcpy(filenames->outfile,"");
  switch(argc){
     case 4:
       strcpy(filenames->errfile,argv[3]);
       set_stderrnow(argv[3]);       
     case 3:
       strcpy(filenames->outfile,argv[2]);       
       set_stdoutnow(argv[2]);
     case 2:
       strcpy(filenames->infile,argv[1]);
       strcpy(filenames->command,argv[0]);
       break;
     default:
       fprintf(stderr,"ERROR: in command line arguments... \n");
       fprintf(stderr,"you must specify at least the sim_parm file\n");
       fprintf(stderr,"argc = %d, \n",argc);
       for(i=0;i<argc;i++) fprintf(stderr,"argv[%d] = \"%s\"\n",i,argv[i]);
       exit(1);
       break;
  }
}

int main_new(int argc,char *argv[]) {

  char name[MAXWORDLEN];
  int pid;
  SIMPARMS simparms;
  FILENAMES filenames;
  COORDS coords;
  SUBSPACE subspace;
  WRITE_STEP write_step;
  NGBR ngbr;
  INTER inter;
  LINE line;
  
  if(argc == 1) usage(argv[0]);
  simparms.mem_bytes = 0;

  gethostname(name,MAXWORDLEN);
  pid = (int) getpid();

#ifdef PARA
  MPI_Init( &argc, &argv );  /* initialize MPI */
  MPI_Comm_rank( MPI_COMM_WORLD, &simparms.rank );
  MPI_Comm_size( MPI_COMM_WORLD, &simparms.size );
#else
  simparms.rank = 0;
  simparms.size = 1;
#endif

  args_filename(&filenames,argc,argv); /* MPI must be initialized first */
  
#ifdef PARA
  sprintf(line,"Running MPI_md %s on %s (id=%d) of rank %d and size %d",
	  argv[0],name,pid,simparms.rank,simparms.size);
  md_stdout(line);
#else
  sprintf(line,"Running serial job %s on %s (id=%d)",argv[0],name,pid);
  md_stdout(line);
#endif

  /* read in the simulation data */
  md_startup(&filenames,&simparms,&subspace,&write_step,&inter,&ngbr,&coords);

  switch(simparms.icalc_type){
  case 0:  /* md simulation */
    md_mdcntr(argv[0],&filenames,&simparms,&write_step,&coords,&inter,&ngbr);
    break;
  case 1:  /* subspace simulation */
    md_subcntr(argv[0],&filenames,&simparms,&write_step,&coords,
	       &subspace,&inter,&ngbr);
    break;
  case 2:   /* minimization */
    min_cntr(argv[0],&filenames,&simparms,&write_step,&coords,&inter,&ngbr);
    break;
  case 3:             /* Normal mode anaylsis */
    nma_control(argv[0],&filenames,&simparms,&coords,&inter,&write_step,
		&subspace);
    break;
  case 4: /* colvar simulation */
    colvarcntr(argv[0],&filenames,&simparms,&write_step,&coords,&inter,&ngbr);
    break;
  default:
    md_error("unknown calculation type\n");
    break;
  }
#ifdef PARA
  MPI_Finalize();
#endif

#ifdef DMALLOC
  dmalloc_verify(NULL);
  dmalloc_shutdown();
#endif

  sprintf(line,"%s on %s has completed.\n:-)",argv[0],name);
  md_stdout(line);

  exit(0);
  return 0;

}

/*-----------------------------------------------------------------------*/
