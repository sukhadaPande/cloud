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
#ifndef HPL_PFACT_H
#define HPL_PFACT_H
/*
 * ---------------------------------------------------------------------
 * Include files
 * ---------------------------------------------------------------------
 */
#include "hpl_cmisc.h"
#include "hpl_cblas.h"
#include "hpl_cgesv.h"

#include "hpl_pcmisc.h"
#include "hpl_pcauxil.h"
#include "hpl_cpanel.h"
/*
 * ---------------------------------------------------------------------
 * #typedefs and data structures
 * ---------------------------------------------------------------------
 */
typedef void (*HPL_CT_PFA_FUN)
(  HPL_CT_panel *,   const int,       const int,       const int,
   ccomplex * );
typedef void (*HPL_CT_RFA_FUN)
(  HPL_CT_panel *,   const int,       const int,       const int,
   ccomplex * );
typedef void (*HPL_CT_UPD_FUN)
(  HPL_CT_panel *,   int *,           HPL_CT_panel *,   const int ); 
/*
 * ---------------------------------------------------------------------
 * Function prototypes
 * ---------------------------------------------------------------------
 */
void                             HPL_clocmax
STDC_ARGS( (
   HPL_CT_panel *,
   const int,
   const int,
   const int,
   ccomplex *
) );

void                             HPL_clocswpN
STDC_ARGS( (
   HPL_CT_panel *,
   const int,
   const int,
   ccomplex *
) );
void                             HPL_clocswpT
STDC_ARGS( (
   HPL_CT_panel *,
   const int,
   const int,
   ccomplex *
) );
void                             HPL_pcmxswp
STDC_ARGS( (
   HPL_CT_panel *,
   const int,
   const int,
   const int,
   ccomplex *
) );

void                             HPL_pcpancrN
STDC_ARGS( (
   HPL_CT_panel *,
   const int,
   const int,
   const int,
   ccomplex *
) );
void                             HPL_pcpancrT
STDC_ARGS( (
   HPL_CT_panel *,
   const int,
   const int,
   const int,
   ccomplex *
) );
void                             HPL_pcpanllN
STDC_ARGS( (
   HPL_CT_panel *,
   const int,
   const int,
   const int,
   ccomplex *
) );
void                             HPL_pcpanllT
STDC_ARGS( (
   HPL_CT_panel *,
   const int,
   const int,
   const int,
   ccomplex *
) );
void                             HPL_pcpanrlN
STDC_ARGS( (
   HPL_CT_panel *,
   const int,
   const int,
   const int,
   ccomplex *
) );
void                             HPL_pcpanrlT
STDC_ARGS( (
   HPL_CT_panel *,
   const int,
   const int,
   const int,
   ccomplex *
) );

void                             HPL_pcrpancrN
STDC_ARGS( (
   HPL_CT_panel *,
   const int,
   const int,
   const int,
   ccomplex *
) );
void                             HPL_pcrpancrT
STDC_ARGS( (
   HPL_CT_panel *,
   const int,
   const int,
   const int,
   ccomplex *
) );
void                             HPL_pcrpanllN
STDC_ARGS( (
   HPL_CT_panel *,
   const int,
   const int,
   const int,
   ccomplex *
) );
void                             HPL_pcrpanllT
STDC_ARGS( (
   HPL_CT_panel *,
   const int,
   const int,
   const int,
   ccomplex *
) );
void                             HPL_pcrpanrlN
STDC_ARGS( (
   HPL_CT_panel *,
   const int,
   const int,
   const int,
   ccomplex *
) );
void                             HPL_pcrpanrlT
STDC_ARGS( (
   HPL_CT_panel *,
   const int,
   const int,
   const int,
   ccomplex *
) );

void                             HPL_pcfact
STDC_ARGS( (
   HPL_CT_panel *
) );
 
#endif
/*
 * End of hpl_pcfact.h
 */
