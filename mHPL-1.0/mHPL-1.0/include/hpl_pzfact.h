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
#include "hpl_zmisc.h"
#include "hpl_zblas.h"
#include "hpl_zgesv.h"

#include "hpl_pzmisc.h"
#include "hpl_pzauxil.h"
#include "hpl_zpanel.h"
/*
 * ---------------------------------------------------------------------
 * #typedefs and data structures
 * ---------------------------------------------------------------------
 */
typedef void (*HPL_ZT_PFA_FUN)
(  HPL_ZT_panel *,   const int,       const int,       const int,
   zcomplex * );
typedef void (*HPL_ZT_RFA_FUN)
(  HPL_ZT_panel *,   const int,       const int,       const int,
   zcomplex * );
typedef void (*HPL_ZT_UPD_FUN)
(  HPL_ZT_panel *,   int *,           HPL_ZT_panel *,   const int ); 
/*
 * ---------------------------------------------------------------------
 * Function prototypes
 * ---------------------------------------------------------------------
 */
void                             HPL_zlocmax
STDC_ARGS( (
   HPL_ZT_panel *,
   const int,
   const int,
   const int,
   zcomplex *
) );

void                             HPL_zlocswpN
STDC_ARGS( (
   HPL_ZT_panel *,
   const int,
   const int,
   zcomplex *
) );
void                             HPL_zlocswpT
STDC_ARGS( (
   HPL_ZT_panel *,
   const int,
   const int,
   zcomplex *
) );
void                             HPL_pzmxswp
STDC_ARGS( (
   HPL_ZT_panel *,
   const int,
   const int,
   const int,
   zcomplex *
) );

void                             HPL_pzpancrN
STDC_ARGS( (
   HPL_ZT_panel *,
   const int,
   const int,
   const int,
   zcomplex *
) );
void                             HPL_pzpancrT
STDC_ARGS( (
   HPL_ZT_panel *,
   const int,
   const int,
   const int,
   zcomplex *
) );
void                             HPL_pzpanllN
STDC_ARGS( (
   HPL_ZT_panel *,
   const int,
   const int,
   const int,
   zcomplex *
) );
void                             HPL_pzpanllT
STDC_ARGS( (
   HPL_ZT_panel *,
   const int,
   const int,
   const int,
   zcomplex *
) );
void                             HPL_pzpanrlN
STDC_ARGS( (
   HPL_ZT_panel *,
   const int,
   const int,
   const int,
   zcomplex *
) );
void                             HPL_pzpanrlT
STDC_ARGS( (
   HPL_ZT_panel *,
   const int,
   const int,
   const int,
   zcomplex *
) );

void                             HPL_pzrpancrN
STDC_ARGS( (
   HPL_ZT_panel *,
   const int,
   const int,
   const int,
   zcomplex *
) );
void                             HPL_pzrpancrT
STDC_ARGS( (
   HPL_ZT_panel *,
   const int,
   const int,
   const int,
   zcomplex *
) );
void                             HPL_pzrpanllN
STDC_ARGS( (
   HPL_ZT_panel *,
   const int,
   const int,
   const int,
   zcomplex *
) );
void                             HPL_pzrpanllT
STDC_ARGS( (
   HPL_ZT_panel *,
   const int,
   const int,
   const int,
   zcomplex *
) );
void                             HPL_pzrpanrlN
STDC_ARGS( (
   HPL_ZT_panel *,
   const int,
   const int,
   const int,
   zcomplex *
) );
void                             HPL_pzrpanrlT
STDC_ARGS( (
   HPL_ZT_panel *,
   const int,
   const int,
   const int,
   zcomplex *
) );

void                             HPL_pzfact
STDC_ARGS( (
   HPL_ZT_panel *
) );
 
#endif
/*
 * End of hpl_pzfact.h
 */
