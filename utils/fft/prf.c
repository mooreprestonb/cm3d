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
   $Author: moore $
   $Revision: 1.1.1.1 $
   Copyright (c) 1997 P. Jeffrey Ungar and the University of Pennsylvania
   */
/**********************************************************************/

static char rcsid[]="$Header: /CVS/cm3d/utils/fft/prf.c,v 1.1.1.1 2008-01-29 00:35:26 moore Exp $";

#include "prf.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>

/*
 * Timer structure with my_advance member to allow both
 * inclusive and exclusive timing measurements
 */
typedef struct timer_s {
   double start;           /* time when last started */
   double my_advance;      /* global timer advance when last started */
   double acc;             /* accumulated active time */
} Timer;

/*
 * Profiler type that tracks wall and cpu times.  Maybe more
 * data to come!
 */
#ifndef PRF_PROFILER_T
#define PRF_PROFILER_T
typedef struct profiler_s *Profiler;
#endif
struct profiler_s {
   char *name;             /* name of this profiler -- may be "" */
   int id;                 /* unique id of this profiler */ 
   Timer wall;             /* wall clock timer */
   Timer cpu;              /* cpu timer */ 
   int ntimes;             /* number of times measured */
   int active;             /* 1 if running, 0 if not started */
};

/*
 * Type to build linked list of all allocated profilers.
 * This way we can dump results later, etc.
 */
typedef struct prf_list_n_s prf_list_n;
struct prf_list_n_s {
   Profiler P;
   prf_list_n *next;
};
static prf_list_n prf_listhead = {0,0};
static prf_list_n *pprf_listtail = &prf_listhead;
static void prfAddProfilerToList(Profiler P);

/*
 * Stack to keep track of profilers as they are started and
 * stopped -- we want to warn when profiled regions cross over
 * each other.
 */
#define PRF_STACKSIZE  100
static int prfStack[PRF_STACKSIZE];        /* stack of profilers */
static int prfStackPos=0;                  /* top of stack */
static void prfPushCheck(const Profiler P);
static void prfPopCheck(const Profiler P);

/*
 * Other functions to assert the state of the profiler when
 * starting or stopping
 */
static void prfAssertActive(const Profiler P);
static void prfAssertInactive(const Profiler P);



static double prfGlobalCpuTimerAdvance=0;  /* to allow exclusive cpu timings */
static double prfGlobalWallTimerAdvance=0; /* to allow exclusive wall timings */
static int prfCount=0;                     /* number of created profilers */ 

static double cpu_time(void);
static double wall_time(void);



void NewProfiler(const char *name,Profiler *pP)
/**************************************************************************/
/* Create a new profiler on the heap and record it in the list of all
 * profilers. If profiler already exists, then do nothing.
 */
{
   Profiler P;
   int nlen;

   if ( !*pP ) {
      if ( !(P = (Profiler) malloc(sizeof(struct profiler_s))) ) {
         fprintf(stderr,"Failed to allocate profiler in NewProfiler()\n");
         exit(1);
      }
      *pP = P;
   } else {
      return;
   }

   if ( name ) {
      nlen = strlen(name);
   } else {
      nlen = 0;
   }
   if ( !(P->name = (char *) malloc(nlen+1)) ) {
      fprintf(stderr,"Failed to allocate profiler name in NewProfiler()\n");
      exit(1);
   }
   if ( name ) {
      strcpy(P->name,name);
   } else {
      strcpy(P->name,"");
   }
   P->wall.start      = 0;
   P->wall.my_advance = 0;
   P->wall.acc        = 0;
   P->cpu.start       = 0;
   P->cpu.my_advance  = 0;
   P->cpu.acc         = 0;
   P->ntimes          = 0;
   P->active          = 0;
   P->id              = prfCount++;
   
   prfAddProfilerToList(P);
}


void ProfileStart(Profiler P)
/**************************************************************************/
/* Record starting times as well as the global timer advances (for
 * timings exclusive of any interpolated profiling starts+stops).
 */
{
   prfAssertInactive(P);
   P->wall.start      = wall_time();
   P->wall.my_advance = prfGlobalWallTimerAdvance;
   P->cpu.start       = cpu_time();
   P->cpu.my_advance  = prfGlobalCpuTimerAdvance;
   P->active = 1;
   prfPushCheck(P);
}


void ProfileStop(Profiler P)
/*************************************************************************/
/* Stop this profiler and record the net accumulated active time
 * exclusive of any other profilers that were active in the same
 * time period.
 */
{
   double wdel,cdel;
   double net_wadv,net_cadv;

   prfPopCheck(P);
   prfAssertActive(P);

   net_cadv = prfGlobalCpuTimerAdvance  - P->cpu.my_advance;
   net_wadv = prfGlobalWallTimerAdvance - P->wall.my_advance;
   cdel = cpu_time() - (P->cpu.start)  - net_cadv;
   wdel = wall_time()- (P->wall.start) - net_wadv;
   
   P->wall.acc += wdel;
   P->cpu.acc  += cdel;
   prfGlobalWallTimerAdvance += wdel;
   prfGlobalCpuTimerAdvance  += cdel;
   
   P->active = 0;
   P->ntimes++;
}


void ProfileStopI(Profiler P)
/*************************************************************************/
/* Stop this profiler and accumulate simple time difference since
 * started.
 */
{
   double wdel,cdel;

   prfPopCheck(P);
   prfAssertActive(P);
      
   cdel  = cpu_time()-(P->cpu.start);
   wdel  = wall_time()-(P->wall.start);
   P->wall.acc += wdel;
   P->cpu.acc  += cdel;
   P->active = 0;
   P->ntimes++;
}


double ProfileGetWallTime(Profiler P)
/*************************************************************************/
{
   return ( P->wall.acc );
}

double ProfileGetCpuTime(Profiler P)
/*************************************************************************/
{
   return ( P->cpu.acc );
}

int ProfileGetCount(Profiler P)
/*************************************************************************/
{
   return ( P->ntimes );
}


void DumpProfileResults(FILE *fp)
/*************************************************************************/
{
   prf_list_n *pprfl = prf_listhead.next;
   
   fprintf(fp,"\n\nProfiling Results\n\n");

   fprintf(fp," ID        Profiler Name  CPU Tot   Wall Tot  calls"
      "  CPU/call  Wall/call\n");
   fprintf(fp,"---------------------------------------------------"
      "---------------------\n");
   while ( pprfl ) {
      Profiler P = pprfl->P;
      char *name  = P->name;
      double wacc = ProfileGetWallTime(P);
      double cacc = ProfileGetCpuTime(P);
      int nt      = ProfileGetCount(P);
      int id      = P->id;
      fprintf(fp,"%3d %20s %8.2f  %9.3f   %4d  %8.2f  %9.3f\n",
         id,name,cacc,wacc,nt,cacc/nt,wacc/nt);
      pprfl = pprfl->next;
   }
   fprintf(fp,"\n");
}


static double cpu_time()
/*************************************************************************/
{
   return ((double) clock())/CLOCKS_PER_SEC;
}

static double wall_time()
/*************************************************************************/
{
   struct timeval t;
   gettimeofday(&t,NULL);
   return ((t.tv_sec-80000000)+1e-6*(t.tv_usec));
}

static void DeleteProfiler(Profiler P)
/**************************************************************************/
/* It doesn't make any sense to use this externally (for now), or internally
 * for that matter! It is a profiler destructor but it would hose the
 * linked list of all profilers, etc.
 */
{
   if ( P ) {
      if ( P->name ) {
         free(P->name);
      }
      free(P);
   }
}



static void prfAddProfilerToList(Profiler P)
/**************************************************************************/
/* adds this profiler to the linked list this module uses to keep track
 * of all existing profilers
 */
{
   pprf_listtail =
      pprf_listtail->next = (prf_list_n *) malloc(sizeof(prf_list_n));
   if ( !pprf_listtail ) {
      fprintf(stderr,"Allocation error adding to Profiler list\n");
      exit(1);
   }
   pprf_listtail->P   = P;
   pprf_listtail->next = 0;
}

static void prfPushCheck(const Profiler P)
/**************************************************************************/
/* Push a profiler onto the stack, checking it first to see that there
 * is room.
 */
{
   if ( prfStackPos>=PRF_STACKSIZE ) {
      fprintf(stderr,"PRF ERROR:   Profiler stack error starting [%3d] (%s)\n",
         P->id,P->name);
      DumpProfileResults(stderr);
      exit(1);
   } else {
      prfStack[prfStackPos++] = P->id;
   }
}

static void prfPopCheck(const Profiler P)
/**************************************************************************/
/* Pop a profiler off the stack (if valid) and check to make sure it is the
 * same as the profiler being stopped.
 */
{
   int isp = prfStackPos-1;
   
   if ( isp<0 || prfStack[isp]!=P->id ) {
      fprintf(stderr,"PRF ERROR:   Profiler stack error stopping [%3d] (%s)\n",
         P->id,P->name);
      DumpProfileResults(stderr);
      exit(1);
   } else {
      prfStackPos=isp;
   }
}


static void prfAssertActive(const Profiler P)
/**************************************************************************/
/* Aborts if the given profiler is not active (used when stopping, since it
 * makes no sense to stop a profiler that isn't running).
 */
{
   if ( !P->active ) {
      fprintf(stderr,"PRF ERROR:   Profiler [%3d] (%s) stopped before starting\n",
         P->id,P->name);
      DumpProfileResults(stderr);
      exit(1);
   }
}

static void prfAssertInactive(const Profiler P)
/**************************************************************************/
/* Aborts if the givne profiler is active (used when starting, since it
 * makes no sense to start a profiler that is already running).
 */
{
   if ( P->active ) {
      fprintf(stderr,"PRF ERROR:   Profiler [%3d] (%s) restarted before stopped\n",
         P->id,P->name);
      DumpProfileResults(stderr);
      exit(1);
   }
}

