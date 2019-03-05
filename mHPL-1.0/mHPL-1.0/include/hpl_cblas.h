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
#ifndef HPL_BLAS_H
#define HPL_BLAS_H
/*
 * ---------------------------------------------------------------------
 * Include files
 * ---------------------------------------------------------------------
 */
#include "hpl_cmisc.h"
/*
 * ---------------------------------------------------------------------
 * typedef definitions
 * ---------------------------------------------------------------------
 */
enum HPL_ORDER
{  HplRowMajor = 101,  HplColumnMajor  = 102 };
enum HPL_TRANS
{  HplNoTrans  = 111,  HplTrans        = 112,  HplConjTrans    = 113 };
enum HPL_UPLO
{  HplUpper    = 121,  HplLower        = 122 };
enum HPL_DIAG
{  HplNonUnit  = 131,  HplUnit         = 132 };
enum HPL_SIDE
{  HplLeft     = 141,  HplRight        = 142 }; 

#ifdef HPL_CALL_CBLAS
/*
 * ---------------------------------------------------------------------
 * The C interface of the BLAS is available ...
 * ---------------------------------------------------------------------
 * #define macro constants
 * ---------------------------------------------------------------------
 */
#define    CBLAS_INDEX         int
 
#define    CBLAS_ORDER         HPL_ORDER
#define    CblasRowMajor       HplRowMajor
#define    CblasColMajor       HplColMajor
 
#define    CBLAS_TRANSPOSE     HPL_TRANS
#define    CblasNoTrans        HplNoTrans
#define    CblasTrans          HplTrans
#define    CblasConjTrans      HplConjTrans
 
#define    CBLAS_UPLO          HPL_UPLO
#define    CblasUpper          HplUpper
#define    CblasLower          HplLower
 
#define    CBLAS_DIAG          HPL_DIAG
#define    CblasNonUnit        HplNonUnit
#define    CblasUnit           HplUnit
 
#define    CBLAS_SIDE          HPL_SIDE
#define    CblasLeft           HplLeft
#define    CblasRight          HplRight
/*
 * ---------------------------------------------------------------------
 * CBLAS Function prototypes
 * ---------------------------------------------------------------------
 */
CBLAS_INDEX       cblas_icamax
STDC_ARGS(
(  const int,       const ccomplex *,  const int ) );
void              cblas_cswap
STDC_ARGS(
(  const int,       ccomplex *,        const int,       ccomplex *,
   const int ) );
void              cblas_ccopy
STDC_ARGS(
(  const int,       const ccomplex *,  const int,       ccomplex *,
   const int ) );
void              cblas_caxpy
STDC_ARGS(
(  const int,       const ccomplex,    const ccomplex *,  const int,
   ccomplex *,        const int ) );
void              cblas_cscal
STDC_ARGS(
(  const int,       const ccomplex,    ccomplex *,        const int ) );

void              cblas_cgemv
STDC_ARGS(
(  const enum CBLAS_ORDER,           const enum CBLAS_TRANSPOSE,
   const int,       const int,       const ccomplex,    const ccomplex *,
   const int,       const ccomplex *,  const int,       const ccomplex,
   ccomplex *,        const int ) );

void              cblas_cger
STDC_ARGS(
(  const enum CBLAS_ORDER,           const int,       const int,
   const ccomplex,    const ccomplex *,  const int,       const ccomplex *,
   const int,       ccomplex *,        const int ) );
void              cblas_ctrsv
STDC_ARGS(
(  const enum CBLAS_ORDER,           const enum CBLAS_UPLO,
   const enum CBLAS_TRANSPOSE,       const enum CBLAS_DIAG,
   const int,       const ccomplex *,  const int,       ccomplex *,
   const int ) );

void              cblas_cgemm
STDC_ARGS(
(  const enum CBLAS_ORDER,           const enum CBLAS_TRANSPOSE,
   const enum CBLAS_TRANSPOSE,       const int,       const int,
   const int,       const ccomplex,    const ccomplex *,  const int,
   const ccomplex *,  const int,       const ccomplex,    ccomplex *,
   const int ) );
void              cblas_ctrsm
STDC_ARGS(
(  const enum CBLAS_ORDER,           const enum CBLAS_SIDE,
   const enum CBLAS_UPLO,            const enum CBLAS_TRANSPOSE,
   const enum CBLAS_DIAG,            const int,       const int,
   const ccomplex,    const ccomplex *,  const int,       ccomplex *,
   const int ) );
/*
 * ---------------------------------------------------------------------
 * HPL C BLAS macro definition
 * ---------------------------------------------------------------------
 */
#define    HPL_cswap           cblas_cswap
#define    HPL_ccopy           cblas_ccopy
#define    HPL_caxpy           cblas_caxpy
#define    HPL_cscal           cblas_cscal
#define    HPL_icamax          cblas_icamax

#define    HPL_cgemv           cblas_cgemv
#define    HPL_ctrsv           cblas_ctrsv
#define    HPL_cger            cblas_cger

#define    HPL_cgemm           cblas_cgemm
#define    HPL_ctrsm           cblas_ctrsm

#endif

#ifdef HPL_CALL_FBLAS
/*
 * ---------------------------------------------------------------------
 * Use the Fortran 77 interface of the BLAS ...
 * ---------------------------------------------------------------------
 * Defaults: Add_, F77_INTEGER=int, StringSunStyle
 * ---------------------------------------------------------------------
 */
#ifndef NoChange
#ifndef UpCase
#ifndef Add__
#ifndef Add_

#define Add_

#endif
#endif
#endif
#endif

#ifndef F77_INTEGER
#define    F77_INTEGER         int
#else
#define    HPL_USE_F77_INTEGER_DEF
#endif

#ifndef StringCrayStyle
#ifndef StringStructVal
#ifndef StringStructPtr
#ifndef StringSunStyle

#define StringSunStyle

#endif
#endif
#endif
#endif
/*
 * ---------------------------------------------------------------------
 * Fortran 77 <-> C interface
 * ---------------------------------------------------------------------
 *
 * These macros identifies how Fortran routines will be called.
 *
 * Add_     : the Fortran compiler expects the name of C functions to be
 * in all lower case and to have an underscore postfixed it (Suns, Intel
 * compilers expect this).
 *
 * NoChange : the Fortran compiler expects the name of C functions to be
 * in all lower case (IBM RS6K compilers do this).
 *
 * UpCase   : the Fortran compiler expects the name of C functions to be
 * in all upcase. (Cray compilers expect this).
 *
 * Add__    : the Fortran compiler in use is f2c, a Fortran to C conver-
 * ter.
 */
#ifdef NoChange
/*
 * These defines  set  up  the  naming scheme required to have a FORTRAN
 * routine called by a C routine with the following  FORTRAN to C inter-
 * face:
 *
 *          FORTRAN DECLARATION            C CALL
 *          SUBROUTINE CGEMM(...)          cgemm(...)
 */
#define    F77cswap               cswap
#define    F77cscal               cscal
#define    F77ccopy               ccopy
#define    F77caxpy               caxpy
#define    F77icamax              icamax

#define    F77cgemv               cgemv
#define    F77ctrsv               ctrsv
#define    F77cger                cgeru

#define    F77cgemm               cgemm
#define    F77ctrsm               ctrsm

#endif

#ifdef UpCase
/*
 * These defines  set  up  the  naming scheme required to have a FORTRAN
 * routine called by a C routine with the following  FORTRAN to C inter-
 * face:
 *
 *          FORTRAN DECLARATION            C CALL
 *          SUBROUTINE CGEMM(...)          CGEMM(...)
 */
#ifdef CRAY_BLAS
                                                                                
#define    F77cswap               CSWAP
#define    F77cscal               CSCAL
#define    F77ccopy               CCOPY
#define    F77caxpy               CAXPY
#define    F77icamax              ICAMAX
                                                                                
#define    F77cgemv               CGEMV
#define    F77ctrsv               CTRSV
#define    F77cger                CGERU
                                                                                
#define    F77cgemm               CGEMM
#define    F77ctrsm               CTRSM
                                                                                
#else

#define    F77cswap               CSWAP
#define    F77cscal               CSCAL
#define    F77ccopy               CCOPY
#define    F77caxpy               CAXPY
#define    F77icamax              ICAMAX

#define    F77cgemv               CGEMV
#define    F77ctrsv               CTRSV
#define    F77cger                CGERU

#define    F77cgemm               CGEMM
#define    F77ctrsm               CTRSM

#endif

#endif

#ifdef Add_
/*
 * These defines  set  up  the  naming scheme required to have a FORTRAN
 * routine called by a C routine  with the following  FORTRAN to C inter-
 * face:
 *
 *          FORTRAN DECLARATION            C CALL
 *          SUBROUTINE CGEMM(...)          cgemm_(...)
 */
#define    F77cswap               cswap_
#define    F77cscal               cscal_
#define    F77ccopy               ccopy_
#define    F77caxpy               caxpy_
#define    F77icamax              icamax_

#define    F77cgemv               cgemv_
#define    F77ctrsv               ctrsv_
#define    F77cger                cgeru_

#define    F77cgemm               cgemm_
#define    F77ctrsm               ctrsm_

#endif

#ifdef Add__
/*
 * These defines  set  up  the  naming scheme required to have a FORTRAN
 * routine called by a C routine  with the following  FORTRAN to C inter-
 * face:
 *
 *          FORTRAN DECLARATION            C CALL
 *          SUBROUTINE CGEMM(...)          cgemm_(...)
 */
#define    F77cswap               cswap_
#define    F77cscal               cscal_
#define    F77ccopy               ccopy_
#define    F77caxpy               caxpy_
#define    F77icamax              icamax_
 
#define    F77cgemv               cgemv_
#define    F77ctrsv               ctrsv_
#define    F77cger                cgeru_
 
#define    F77cgemm               cgemm_
#define    F77ctrsm               ctrsm_
 
#endif
/*
 * ---------------------------------------------------------------------
 * Typedef definitions and conversion utilities
 * ---------------------------------------------------------------------
 */
#ifdef StringCrayStyle

#include <fortran.h>
                      /* Type of character argument in a FORTRAN call */
#define    F77_CHAR            _fcd
                                    /* Character conversion utilities */
#define    HPL_F2C_CHAR(c)     (*(_fcdtocp(c) ))
#define    HPL_C2F_CHAR(c)     (_cptofcd(&(c), 1))

#define    F77_CHAR_DECL       F77_CHAR          /* input CHARACTER*1 */

#endif
/* ------------------------------------------------------------------ */
#ifdef StringStructVal
                      /* Type of character argument in a FORTRAN call */
typedef struct { char *cp; F77_INTEGER len; } F77_CHAR;
                                    /* Character conversion utilities */
#define    HPL_F2C_CHAR(c)     (*(c.cp))

#define    F77_CHAR_DECL       F77_CHAR          /* input CHARACTER*1 */

#endif
/* ------------------------------------------------------------------ */
#ifdef StringStructPtr
                      /* Type of character argument in a FORTRAN call */
typedef struct { char *cp; F77_INTEGER len; } F77_CHAR;
                                    /* Character conversion utilities */
#define    HPL_F2C_CHAR(c)     (*(c->cp))

#define    F77_CHAR_DECL       F77_CHAR *        /* input CHARACTER*1 */

#endif
/* ------------------------------------------------------------------ */
#ifdef StringSunStyle
                      /* Type of character argument in a FORTRAN call */
#define    F77_CHAR            char *
                                    /* Character conversion utilities */
#define    HPL_F2C_CHAR(c)     (*(c))
#define    HPL_C2F_CHAR(c)     (&(c))

#define    F77_CHAR_DECL       F77_CHAR          /* input CHARACTER*1 */
#define    F77_1_CHAR          , F77_INTEGER
#define    F77_2_CHAR          F77_1_CHAR F77_1_CHAR
#define    F77_3_CHAR          F77_2_CHAR F77_1_CHAR
#define    F77_4_CHAR          F77_3_CHAR F77_1_CHAR

#endif
/* ------------------------------------------------------------------ */

#ifndef F77_1_CHAR
#define    F77_1_CHAR
#define    F77_2_CHAR
#define    F77_3_CHAR
#define    F77_4_CHAR
#endif

#define    F77_INT_DECL        const F77_INTEGER *   /* input integer */
#define    F77_SIN_DECL        const ccomplex *         /* input scalar */
#define    F77_VIN_DECL        const ccomplex *         /* input vector */
#define    F77_VINOUT_DECL     ccomplex *        /* input/output matrix */
#define    F77_MIN_DECL        const ccomplex *         /* input matrix */
#define    F77_MINOUT_DECL     ccomplex *        /* input/output matrix */
 
#ifdef CRAY_PVP_ENV                      /* Type of FORTRAN functions */
#define    F77_VOID_FUN        extern fortran void      /* subroutine */
#define    F77_INT_FUN         extern fortran int /* integer function */
#else
#define    F77_VOID_FUN        extern void              /* subroutine */
#define    F77_INT_FUN         extern int         /* integer function */
#endif
/*
 * ---------------------------------------------------------------------
 * Fortran 77 BLAS function prototypes
 * ---------------------------------------------------------------------
 */
F77_VOID_FUN    F77cswap
STDC_ARGS(
(  F77_INT_DECL,    F77_VINOUT_DECL, F77_INT_DECL,    F77_VINOUT_DECL,
   F77_INT_DECL ) );
F77_VOID_FUN    F77cscal
STDC_ARGS(
(  F77_INT_DECL,    F77_SIN_DECL,    F77_VINOUT_DECL, F77_INT_DECL ) );
F77_VOID_FUN    F77ccopy
STDC_ARGS(
(  F77_INT_DECL,    F77_VIN_DECL,    F77_INT_DECL,    F77_VINOUT_DECL,
   F77_INT_DECL ) );
F77_VOID_FUN    F77caxpy
STDC_ARGS(
(  F77_INT_DECL,    F77_SIN_DECL,    F77_VIN_DECL,    F77_INT_DECL,
   F77_VINOUT_DECL, F77_INT_DECL ) );
F77_INT_FUN     F77icamax
STDC_ARGS(
(  F77_INT_DECL,    F77_VIN_DECL,    F77_INT_DECL ) );

F77_VOID_FUN    F77cgemv
STDC_ARGS(
(  F77_CHAR_DECL,   F77_INT_DECL,    F77_INT_DECL,    F77_SIN_DECL,
   F77_MIN_DECL,    F77_INT_DECL,    F77_VIN_DECL,    F77_INT_DECL,
   F77_SIN_DECL,    F77_VINOUT_DECL, F77_INT_DECL     F77_1_CHAR ) );
F77_VOID_FUN    F77cger
STDC_ARGS(
(  F77_INT_DECL,    F77_INT_DECL,    F77_SIN_DECL,    F77_VIN_DECL,
   F77_INT_DECL,    F77_VIN_DECL,    F77_INT_DECL,    F77_MINOUT_DECL,
   F77_INT_DECL ) );
F77_VOID_FUN    F77ctrsv
STDC_ARGS(
(  F77_CHAR_DECL,   F77_CHAR_DECL,   F77_CHAR_DECL,   F77_INT_DECL,
   F77_MIN_DECL,    F77_INT_DECL,    F77_VINOUT_DECL, F77_INT_DECL
   F77_3_CHAR ) );

F77_VOID_FUN    F77cgemm
STDC_ARGS(
(  F77_CHAR_DECL,   F77_CHAR_DECL,   F77_INT_DECL,    F77_INT_DECL,
   F77_INT_DECL,    F77_SIN_DECL,    F77_MIN_DECL,    F77_INT_DECL,
   F77_MIN_DECL,    F77_INT_DECL,    F77_SIN_DECL,    F77_MINOUT_DECL,
   F77_INT_DECL     F77_2_CHAR ) );
F77_VOID_FUN    F77ctrsm
STDC_ARGS(
(  F77_CHAR_DECL,   F77_CHAR_DECL,   F77_CHAR_DECL,   F77_CHAR_DECL,
   F77_INT_DECL,    F77_INT_DECL,    F77_SIN_DECL,    F77_MIN_DECL,
   F77_INT_DECL,    F77_MINOUT_DECL, F77_INT_DECL     F77_4_CHAR ) );

#endif
/*
 * ---------------------------------------------------------------------
 * HPL BLAS Function prototypes
 * ---------------------------------------------------------------------
 */
#ifndef HPL_CALL_CBLAS

int                              HPL_icamax
STDC_ARGS( (
   const int,
   const ccomplex *,
   const int
) );
void                             HPL_caxpy
STDC_ARGS( (
   const int,
   const ccomplex,
   const ccomplex *,
   const int,
   ccomplex *,
   const int
) );
void                             HPL_ccopy
STDC_ARGS( (
   const int,
   const ccomplex *,
   const int,
   ccomplex *,
   const int
) );
void                             HPL_cscal
STDC_ARGS( (
   const int,
   const ccomplex,
   ccomplex *,
   const int
) );
void                             HPL_cswap
STDC_ARGS( (
   const int,
   ccomplex *,
   const int,
   ccomplex *,
   const int
) );
void                             HPL_cgemv
STDC_ARGS( (
   const enum HPL_ORDER,
   const enum HPL_TRANS,
   const int,
   const int,
   const ccomplex,
   const ccomplex *,
   const int,
   const ccomplex *,
   const int,
   const ccomplex,
   ccomplex *,
   const int
) );
void                             HPL_cger
STDC_ARGS( (
   const enum HPL_ORDER,
   const int,
   const int,
   const ccomplex,
   const ccomplex *,
   const int,
   ccomplex *,
   const int,
   ccomplex *,
   const int
) );
void                             HPL_ctrsv
STDC_ARGS( (
   const enum HPL_ORDER,
   const enum HPL_UPLO,
   const enum HPL_TRANS,
   const enum HPL_DIAG,
   const int,
   const ccomplex *,
   const int,
   ccomplex *,
   const int
) );
void                             HPL_cgemm
STDC_ARGS( (
   const enum HPL_ORDER,
   const enum HPL_TRANS,
   const enum HPL_TRANS,
   const int,
   const int,
   const int,
   const ccomplex,
   const ccomplex *,
   const int,
   const ccomplex *,
   const int,
   const ccomplex,
   ccomplex *,
   const int
) );
void                             HPL_ctrsm
STDC_ARGS( (
   const enum HPL_ORDER,
   const enum HPL_SIDE,
   const enum HPL_UPLO,
   const enum HPL_TRANS,
   const enum HPL_DIAG,
   const int,
   const int,
   const ccomplex,
   const ccomplex *,
   const int,
   ccomplex *,
   const int
) );

#endif

#endif
/*
 * hpl_cblas.h
 */
