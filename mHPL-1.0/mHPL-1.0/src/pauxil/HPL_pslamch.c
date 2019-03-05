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
 * ---------------------------------------------------------------------
 */ 
/*
 * Include files
 */
#include "shpl.h"

#ifdef STDC_HEADERS
float HPL_pslamch
(
   MPI_Comm                         COMM,
   const HPL_T_MACH                 CMACH
)
#else
float HPL_pslamch
( COMM, CMACH )
   MPI_Comm                         COMM;
   const HPL_T_MACH                 CMACH;
#endif
{
/* 
 * Purpose
 * =======
 *
 * HPL_pslamch determines  machine-specific  arithmetic  constants  such  as
 * the relative machine precision (eps),  the safe minimum(sfmin) such that
 * 1/sfmin does not overflow, the base of the machine (base), the precision
 * (prec),  the  number  of  (base)  digits in the  mantissa  (t),  whether
 * rounding occurs in addition (rnd = 1.0 and 0.0 otherwise),  the  minimum
 * exponent before  (gradual)  underflow (emin),  the  underflow  threshold
 * (rmin)- base**(emin-1), the largest exponent before overflow (emax), the
 * overflow threshold (rmax)  - (base**emax)*(1-eps).
 *
 * Arguments
 * =========
 *
 * COMM    (global/local input)          MPI_Comm
 *         The MPI communicator identifying the process collection.
 *
 * CMACH   (global input)                const HPL_T_MACH
 *         Specifies the value to be returned by HPL_pslamch            
 *            = HPL_MACH_EPS,   HPL_pslamch := eps (default)            
 *            = HPL_MACH_SFMIN, HPL_pslamch := sfmin                    
 *            = HPL_MACH_BASE,  HPL_pslamch := base                     
 *            = HPL_MACH_PREC,  HPL_pslamch := eps*base                 
 *            = HPL_MACH_MLEN,  HPL_pslamch := t                        
 *            = HPL_MACH_RND,   HPL_pslamch := rnd                      
 *            = HPL_MACH_EMIN,  HPL_pslamch := emin                     
 *            = HPL_MACH_RMIN,  HPL_pslamch := rmin                     
 *            = HPL_MACH_EMAX,  HPL_pslamch := emax                     
 *            = HPL_MACH_RMAX,  HPL_pslamch := rmax                     
 *          
 *         where                                                        
 *          
 *            eps   = relative machine precision,                       
 *            sfmin = safe minimum,                                     
 *            base  = base of the machine,                              
 *            prec  = eps*base,                                         
 *            t     = number of digits in the mantissa,                 
 *            rnd   = 1.0 if rounding occurs in addition,               
 *            emin  = minimum exponent before underflow,                
 *            rmin  = underflow threshold,                              
 *            emax  = largest exponent before overflow,                 
 *            rmax  = overflow threshold.
 *
 * ---------------------------------------------------------------------
 */ 
/*
 * .. Local Variables ..
 */
   float                     param;
/* ..
 * .. Executable Statements ..
 */
   param = HPL_slamch( CMACH );

   switch( CMACH )
   {
      case HPL_MACH_EPS   :
      case HPL_MACH_SFMIN :
      case HPL_MACH_EMIN  :
      case HPL_MACH_RMIN  :
         (void) HPL_sall_reduce( (void *)(&param), 1, HPL_FLOAT,
                                HPL_smax, COMM );
         break;
      case HPL_MACH_EMAX  :
      case HPL_MACH_RMAX  :
         (void) HPL_sall_reduce( (void *)(&param), 1, HPL_FLOAT,
                                HPL_smin, COMM );
         break;
      default             :
         break;
   } 

   return( param );
/*
 * End of HPL_pslamch
 */
}
