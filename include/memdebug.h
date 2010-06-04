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
/********************************************************************************/
/********************************************************************************/
/*                               */
/* Copyright (C) 1992                        */
/* Centre de Recherche Public Henri Tudor (CRP-HT)            */
/* 6, rue Coudenhove-Kalergi                     */
/* L1359 Luxembourg-Kirchberg                     */
/*                               */
/* Authors : Schmit Rene                      */
/*                               */
/* This software may be copied, distributed, ported and modified in source or   */
/* object format as long as :                     */
/*                               */
/*    1) No distribution for commercial purposes is made.         */
/*    2) No third-party copyrights (such as runtime licenses) are involved   */
/*    3) This copyright notice is not removed or changed.         */
/*                               */
/* No responsibility is assumed for any damages that may result       */
/* from any defect in this software.                  */
/*                               */
/********************************************************************************/

/*
**   If no memory checking is required (i.e., memdebug is undefined at compilation time)
**      the memory checking functions are wiped out by the following macro definitions.
*/

#ifndef __memdebug__
#define __memdebug__

#include <stdlib.h>

#ifndef MEMDEBUG

#define VarInfo(p1)
#define generate_MemdebugError()
#define   print_MemdebugStatistics()
#define   check_MemdebugError()      0

#define set_MemdebugOptions(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11)

/* replaced by above expression, as the preprocessor may have difficulties with the continuation character

#define set_MemdebugOptions(p_GeneralStatistics,p_AlphabeticalList,\
p_NotFreeList,p_CallSequenceList,p_SpuriousFreeList,\
p_PrintContents,p_DestroyContents,p_GenerateErrorCount,p_MaximalMemoryAvailable,\
p_StatisticsFileName,p_ErrorFileName)
*/

#endif

/*
   We need the following definitions only when memory checking is required at compilation
      time. This is indicated by defining the memdebug pre-processor symbol.
       
   The VAXC preprocessor doesn't implement stringization: "???"
*/

#ifdef MEMDEBUG
#ifndef VAX

#define malloc(p_Size)    __check_malloc   (p_Size,#p_Size,__FILE__,__LINE__)
#define calloc(p_Number,p_Size)   __check_calloc   (p_Number,p_Size,#p_Size,__FILE__,__LINE__)
#define realloc(p_Name,p_Size)   __check_realloc   (p_Name,p_Size,#p_Name,#p_Size,__FILE__,__LINE__)
#define free(p_Name)       __check_free   (p_Name,#p_Name,__FILE__,__LINE__)

#else

#define malloc(p_Size)       __check_malloc   (p_Size,"p_Size",__FILE__,__LINE__)
#define calloc(p_Number,p_Size)   __check_calloc   (p_Number,p_Size,"p_Size",__FILE__,__LINE__)
#define realloc(p_Name,p_Size)   __check_realloc   (p_Name,p_Size,   "p_Name","p_Size",__FILE__,__LINE__)
#define free(p_Name)      __check_free   (p_Name,"p_Name",__FILE__,__LINE__)

#endif
#endif

/**************************************************************************/
/******************************* Data Types *******************************/
/**************************************************************************/

#ifndef __t_biState__
#define __t_biState__

/*
   This is kind of a BOOLEAN type, but with more and thus clearer
      names for the different states
*/

enum t_biState 
{
   c_Off   = 0,
   c_On,
   
   c_False = 0,
   c_True,
   
   c_No   = 0,
   c_Yes
};

typedef enum t_biState t_biState;

#endif

/**************************************************************************/

#ifdef _MSDOS
#define __MSDOS__
#endif

#ifndef __MSDOS__
#ifndef macintosh
#define macintosh
#endif
#endif

#ifdef macintosh
#endif

#ifdef __MSDOS__
#endif

#ifdef MEMDEBUG

#ifdef __cplusplus
extern "C" {
#endif

/*************************************************************************/
/************************* Function Prototypes ***************************/
/*************************************************************************/

void  * __check_malloc     (   size_t      p_Size,
                     char     * p_SizeExpression,
                     char     * p_FileName,
                     long      p_LineNumber);
                  
void  * __check_calloc     (   size_t      p_Number,
                     size_t      p_Size,
                     char     * p_SizeExpression,
                     char     * p_FileName,
                     long      p_LineNumber);
                  
void  *   __check_realloc     ( void     * p_Pointer,
                     size_t      p_Size,
                     char     *   p_NameString,
                     char     *   p_SizeExpression,
                     char     * p_FileName,
                     long      p_LineNumber);
   
void   __check_free     ( void     * p_Pointer,
                     char     * p_NameString,
                     char     * p_FileName,
                     long      p_LineNumber);

void generate_MemdebugError        ( void );
void print_MemdebugStatistics     ( void );
int    check_MemdebugError        ( void );
void set_MemdebugOptions        ( t_biState      p_GeneralStatistics,
                           t_biState      p_AlphabeticalList,
                           t_biState      p_NotFreeList,
                           t_biState      p_CallSequenceList,
      
                           t_biState      p_SpuriousFreeList,
                           
                           t_biState      p_PrintContents,
                           t_biState      p_DestroyContents,
                           
                           long         p_GenerateErrorCount,
                           unsigned long   p_MaximalMemoryAvailable,
                           
                           char        * p_StatisticsFileName,
                           char        * p_ErrorFileName
                          );


#ifdef __cplusplus

}   // close the extern "C" declaration

/********************************************************************************/


#define VarInfo(p_VarName)   __set_LineInfo(#p_VarName,__FILE__,__LINE__);

void __set_LineInfo (   char  * p_VariableName,
                  char  * p_FileName,
                  int      p_LineNumber);
                  
#endif

#endif
#endif

