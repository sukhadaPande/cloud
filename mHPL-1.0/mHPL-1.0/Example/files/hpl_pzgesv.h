/* 
 * -- High Performance Computing Linpack Benchmark (HPL)                
 *    HPL - 1.0a - January 20, 2004                          
 *    Antoine P. Petitet                                                
 *    University of Tennessee, Knoxville                                
 *    Innovative Computing Laboratories                                 
 *    (C) Copyright 2000-2004 All Rights Reserved                       
 *                                                                      
 * -- Copyright notice and Licensing terms:                             
 *                                                                      
 * Redistribution  and  use in  source and binary forms, with or without
 * modification, are  permitted provided  that the following  conditions
 * are met:                                                             
 *                                                                      
 * 1. Redistributions  of  source  code  must retain the above copyright
 * notice, this list of conditions and the following disclaimer.        
 *                                                                      
 * 2. Redistributions in binary form must reproduce  the above copyright
 * notice, this list of conditions,  and the following disclaimer in the
 * documentation and/or other materials provided with the distribution. 
 *                                                                      
 * 3. All  advertising  materials  mentioning  features  or  use of this
 * software must display the following acknowledgement:                 
 * This  product  includes  software  developed  at  the  University  of
 * Tennessee, Knoxville, Innovative Computing Laboratories.             
 *                                                                      
 * 4. The name of the  University,  the name of the  Laboratory,  or the
 * names  of  its  contributors  may  not  be used to endorse or promote
 * products  derived   from   this  software  without  specific  written
 * permission.                                                          
 *                                                                      
 * -- Disclaimer:                                                       
 *                                                                      
 * THIS  SOFTWARE  IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES,  INCLUDING,  BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE UNIVERSITY
 * OR  CONTRIBUTORS  BE  LIABLE FOR ANY  DIRECT,  INDIRECT,  INCIDENTAL,
 * SPECIAL,  EXEMPLARY,  OR  CONSEQUENTIAL DAMAGES  (INCLUDING,  BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA OR PROFITS; OR BUSINESS INTERRUPTION)  HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT,  STRICT LIABILITY,  OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
 */ 
#ifndef HPL_PGESV_H
#define HPL_PGESV_H
/*
 * ---------------------------------------------------------------------
 * Include files
 * ---------------------------------------------------------------------
 */
#include "hpl_zmisc.h"
#include "hpl_zblas.h"
#include "hpl_zauxil.h"

#include "hpl_pzmisc.h"
#include "hpl_zgrid.h"
#include "hpl_zcomm.h"
#include "hpl_pzauxil.h"
#include "hpl_zpanel.h"
#include "hpl_pzfact.h"
/*
 * ---------------------------------------------------------------------
 * #typedefs and data structures
 * ---------------------------------------------------------------------
 */
typedef enum
{
   HPL_SWAP00        = 451,                      /* Use HPL_pzlaswp00 */
   HPL_SWAP01        = 452,                      /* Use HPL_pzlaswp01 */
   HPL_SW_MIX        = 453, /* Use HPL_pzlaswp00_ for small number of */
                            /* columns, and HPL_pzlaswp01_ otherwise. */
   HPL_NO_SWP        = 499
} HPL_T_SWAP;

typedef struct HPL_ZS_palg
{
   HPL_T_TOP           btopo;               /* row broadcast topology */
   int                 depth;                     /* look-ahead depth */
   int                 nbdiv;            /* recursive division factor */
   int                 nbmin;         /* recursion stopping criterium */
   HPL_T_FACT          pfact;                   /* panel fact variant */
   HPL_T_FACT          rfact;               /* recursive fact variant */
   HPL_ZT_PFA_FUN       pffun;              /* panel fact function ptr */
   HPL_ZT_RFA_FUN       rffun;          /* recursive fact function ptr */
   HPL_ZT_UPD_FUN       upfun;                      /* update function */
   HPL_T_SWAP          fswap;                   /* Swapping algorithm */
   int                 fsthr;                   /* Swapping threshold */
   int                 equil;                        /* Equilibration */
   int                 align;              /* data alignment constant */
} HPL_ZT_palg;

typedef struct HPL_ZS_pmat
{
#ifdef HPL_CALL_VSIPL
   vsip_block_d        * block;
#endif
   zcomplex              * A;            /* pointer to local piece of A */
   zcomplex              * X;             /* pointer to solution vector */
   int                 n;                      /* global problem size */
   int                 nb;                         /* blocking factor */
   int                 ld;                 /* local leading dimension */
   int                 mp;                    /* local number of rows */
   int                 nq;                 /* local number of columns */
   int                 info;                    /* computational flag */
   int                 * global_ipiv;      /* pointer to global vector */
} HPL_ZT_pmat;
/*
 * ---------------------------------------------------------------------
 * #define macro constants
 * ---------------------------------------------------------------------
 */
#define    MSGID_BEGIN_PFACT   1001              /* message id ranges */
#define    MSGID_END_PFACT     2000
#define    MSGID_BEGIN_FACT    2001
#define    MSGID_END_FACT      3000
#define    MSGID_BEGIN_PTRSV   3001
#define    MSGID_END_PTRSV     4000
 
#define    MSGID_BEGIN_COLL    9001
#define    MSGID_END_COLL     10000
/*
 * ---------------------------------------------------------------------
 * #define macros definitions
 * ---------------------------------------------------------------------
 */
#define    MNxtMgid( id_, beg_, end_ ) \
                             (( (id_)+1 > (end_) ?  (beg_) : (id_)+1 ))
/*
 * ---------------------------------------------------------------------
 * Function prototypes
 * ---------------------------------------------------------------------
 */
void                             HPL_pzipid
STDC_ARGS( (
   HPL_ZT_panel *,
   int *,
   int *
) );
void                             HPL_pzlindx0
STDC_ARGS( (
   HPL_ZT_panel *,
   const int,
   int *,
   int *,
   int *,
   int *
) );
void                             HPL_pzlaswp00N
STDC_ARGS( (
   HPL_ZT_panel *,
   int *,
   HPL_ZT_panel *,
   const int
) );
void                             HPL_pzlaswp00T
STDC_ARGS( (
   HPL_ZT_panel *,
   int *,
   HPL_ZT_panel *,
   const int
) );

void                             HPL_perm
STDC_ARGS( (
   const int,
   int *,
   int *,
   int *
) );
void                             HPL_logsort
STDC_ARGS( (
   const int,
   const int,
   int *,
   int *,
   int *
) );
void                             HPL_pzlindx10
STDC_ARGS( (
   HPL_ZT_panel *,
   const int,
   const int *,
   int *,
   int *,
   int *
) );
void                             HPL_pzlindx1
STDC_ARGS( (
   HPL_ZT_panel *,
   const int,
   const int *,
   int *,
   int *,
   int *,
   int *,
   int *,
   int *,
   int *,
   int *
) );
void                             HPL_zspreadN
STDC_ARGS( (
   HPL_ZT_panel *,
   int *,
   HPL_ZT_panel *,
   const enum HPL_SIDE,
   const int,
   zcomplex *,
   const int,
   const int,
   const int *,
   const int *,
   const int *
) );
void                             HPL_zspreadT
STDC_ARGS( (
   HPL_ZT_panel *,
   int *,
   HPL_ZT_panel *,
   const enum HPL_SIDE,
   const int,
   zcomplex *,
   const int,
   const int,
   const int *,
   const int *,
   const int *
) );
void                             HPL_zequil
STDC_ARGS( (
   HPL_ZT_panel *,
   int *,
   HPL_ZT_panel *,
   const enum HPL_TRANS,
   const int,
   zcomplex *,
   const int,
   int *,
   const int *,
   const int *,
   int *
) );
void                             HPL_zrollN
STDC_ARGS( (
   HPL_ZT_panel *,
   int *,
   HPL_ZT_panel *,
   const int,
   zcomplex *,
   const int,
   const int *,
   const int *,
   const int *
) );
void                             HPL_zrollT
STDC_ARGS( (
   HPL_ZT_panel *,
   int *,
   HPL_ZT_panel *,
   const int,
   zcomplex *,
   const int,
   const int *,
   const int *,
   const int *
) );
void                             HPL_pzlaswp01N
STDC_ARGS( (
   HPL_ZT_panel *,
   int *,
   HPL_ZT_panel *,
   const int
) );
void                             HPL_pzlaswp01T
STDC_ARGS( (
   HPL_ZT_panel *,
   int *,
   HPL_ZT_panel *,
   const int
) );

void                             HPL_pzupdateNN
STDC_ARGS( (
   HPL_ZT_panel *,
   int *,
   HPL_ZT_panel *,
   const int
) );
void                             HPL_pzupdateNT
STDC_ARGS( (
   HPL_ZT_panel *,
   int *,
   HPL_ZT_panel *,
   const int
) );
void                             HPL_pzupdateTN
STDC_ARGS( (
   HPL_ZT_panel *,
   int *,
   HPL_ZT_panel *,
   const int
) );
void                             HPL_pzupdateTT
STDC_ARGS( (
   HPL_ZT_panel *,
   int *,
   HPL_ZT_panel *,
   const int
) );

void                             HPL_pzgesv0
STDC_ARGS( (
   HPL_T_grid *,
   HPL_ZT_palg *,
   HPL_ZT_pmat *
) );
void                             HPL_pzgesvK1
STDC_ARGS( (
   HPL_T_grid *,
   HPL_ZT_palg *,
   HPL_ZT_pmat *
) );
void                             HPL_pzgesvK2
STDC_ARGS( (
   HPL_T_grid *,
   HPL_ZT_palg *,
   HPL_ZT_pmat *
) );
void                             HPL_pzgesv
STDC_ARGS( (
   HPL_T_grid *,
   HPL_ZT_palg *,
   HPL_ZT_pmat *
) );
 
void                             HPL_pztrsv
STDC_ARGS( (
   HPL_T_grid *,
   HPL_ZT_pmat *
) );

#endif
/*
 * End of hpl_pzgesv.h
 */ 
