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
#ifndef HPL_COMM_H
#define HPL_COMM_H
/*
 * ---------------------------------------------------------------------
 * Include files
 * ---------------------------------------------------------------------
 */
#include "hpl_pcmisc.h"
#include "hpl_cpanel.h"
/*
 * ---------------------------------------------------------------------
 * #typedefs and data structures
 * ---------------------------------------------------------------------
 */
typedef enum
{
   HPL_1RING         = 401,                        /* Increasing ring */
   HPL_1RING_M       = 402,             /* Increasing ring (modified) */
   HPL_2RING         = 403,                      /* Increasing 2-ring */
   HPL_2RING_M       = 404,           /* Increasing 2-ring (modified) */
   HPL_BLONG         = 405,                         /* long broadcast */
   HPL_BLONG_M       = 406               /* long broadcast (modified) */
} HPL_T_TOP;
/*
 * ---------------------------------------------------------------------
 * #define macro constants
 * ---------------------------------------------------------------------
 */
#define    HPL_FAILURE            0
#define    HPL_SUCCESS            1
#define    HPL_KEEP_TESTING       2
/*
 * ---------------------------------------------------------------------
 * comm function prototypes
 * ---------------------------------------------------------------------
 */
int                              HPL_csend
STDC_ARGS( (
   ccomplex *,
   int,
   int,
   int,
   MPI_Comm
) );
int                              HPL_crecv
STDC_ARGS( (
   ccomplex *,
   int,
   int,
   int,
   MPI_Comm
) );
int                              HPL_csdrv
STDC_ARGS( (
   ccomplex *,
   int,
   int,
   ccomplex *,
   int,
   int,
   int,
   MPI_Comm
) );
int                              HPL_cbinit
STDC_ARGS( (
   HPL_CT_panel *
) );
int                              HPL_cbcast
STDC_ARGS( (
   HPL_CT_panel *,
   int *
) );
int                              HPL_cbwait
STDC_ARGS( (
   HPL_CT_panel *
) );
int                              HPL_cpackL
STDC_ARGS( (
   HPL_CT_panel *,
   const int,
   const int,
   const int
) );
void                             HPL_ccopyL
STDC_ARGS( (
   HPL_CT_panel *
) );
 
int HPL_cbinit_1ring STDC_ARGS( ( HPL_CT_panel *        ) );
int HPL_cbcast_1ring STDC_ARGS( ( HPL_CT_panel *, int * ) );
int HPL_cbwait_1ring STDC_ARGS( ( HPL_CT_panel *        ) );
 
int HPL_cbinit_1rinM STDC_ARGS( ( HPL_CT_panel *        ) );
int HPL_cbcast_1rinM STDC_ARGS( ( HPL_CT_panel *, int * ) );
int HPL_cbwait_1rinM STDC_ARGS( ( HPL_CT_panel *        ) );
 
int HPL_cbinit_2ring STDC_ARGS( ( HPL_CT_panel *        ) );
int HPL_cbcast_2ring STDC_ARGS( ( HPL_CT_panel *, int * ) );
int HPL_cbwait_2ring STDC_ARGS( ( HPL_CT_panel *        ) );
 
int HPL_cbinit_2rinM STDC_ARGS( ( HPL_CT_panel *        ) );
int HPL_cbcast_2rinM STDC_ARGS( ( HPL_CT_panel *, int * ) );
int HPL_cbwait_2rinM STDC_ARGS( ( HPL_CT_panel *        ) );
 
int HPL_cbinit_blong STDC_ARGS( ( HPL_CT_panel *        ) );
int HPL_cbcast_blong STDC_ARGS( ( HPL_CT_panel *, int * ) );
int HPL_cbwait_blong STDC_ARGS( ( HPL_CT_panel *        ) );
 
int HPL_cbinit_blonM STDC_ARGS( ( HPL_CT_panel *        ) );
int HPL_cbcast_blonM STDC_ARGS( ( HPL_CT_panel *, int * ) );
int HPL_cbwait_blonM STDC_ARGS( ( HPL_CT_panel *        ) );

#endif
/*
 * End of hpl_ccomm.h
 */
