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
#include "zhpl.h"
 
#ifdef STDC_HEADERS
void HPL_pzgesvK1
(
   HPL_T_grid *                     GRID,
   HPL_ZT_palg *                     ALGO,
   HPL_ZT_pmat *                     A
)
#else
void HPL_pzgesvK1
( GRID, ALGO, A )
   HPL_T_grid *                     GRID;
   HPL_ZT_palg *                     ALGO;
   HPL_ZT_pmat *                     A;
#endif
{
/* 
 * Purpose
 * =======
 *
 * HPL_pzgesvK1 factors a N+1-by-N matrix using LU factorization with row
 * partial pivoting.  The main algorithm  is the "right looking" variant
 * with look-ahead.  The  lower  triangular factor is left unpivoted and
 * the pivots are not returned. The right hand side is the N+1 column of
 * the coefficient matrix.
 *
 * Arguments
 * =========
 *
 * GRID    (local input)                 HPL_T_grid *
 *         On entry,  GRID  points  to the data structure containing the
 *         process grid information.
 *
 * ALGO    (global input)                HPL_ZT_palg *
 *         On entry,  ALGO  points to  the data structure containing the
 *         algorithmic parameters.
 *
 * A       (local input/output)          HPL_ZT_pmat *
 *         On entry, A points to the data structure containing the local
 *         array information.
 *
 * ---------------------------------------------------------------------
 */ 
/*
 * .. Local Variables ..
 */
   HPL_ZT_panel                * * panel = NULL;
   HPL_ZT_UPD_FUN              HPL_pzupdate;
   int                        N, depth, icurcol=0, j, jb, jj=0, jstart,
                              k, mycol, n, nb, nn, npcol, nq,
                              tag=MSGID_BEGIN_FACT, test=HPL_KEEP_TESTING;
/* ..
 * .. Executable Statements ..
 */
   mycol = GRID->mycol; npcol        = GRID->npcol;
   depth = ALGO->depth; HPL_pzupdate = ALGO->upfun;
   N     = A->n;        nb           = A->nb; 

   if( N <= 0 ) return;
/*
 * Allocate a panel list of length depth + 1 (depth >= 1)
 */
   panel = (HPL_ZT_panel **)malloc( (depth+1)*sizeof( HPL_ZT_panel *) );
   if( panel == NULL )
   { HPL_pabort( __LINE__, "HPL_pzgesvK1", "Memory allocation failed" ); }
/*
 * Create and initialize the first depth panels
 */
   nq = HPL_numroc( N+1, nb, nb, mycol, 0, npcol ); nn = N; jstart = 0;

   for( k = 0; k < depth; k++ )
   {
      jb = Mmin( nn, nb );
      HPL_pzpanel_new( GRID, ALGO, nn, nn+1, jb, A, jstart, jstart,
                       tag, &panel[k] );
      nn -= jb; jstart += jb;
      if( mycol == icurcol ) { jj += jb; nq -= jb; }
      icurcol = MModAdd1( icurcol, npcol );
      tag     = MNxtMgid( tag, MSGID_BEGIN_FACT, MSGID_END_FACT );
   }
/*
 * Initialize the lookahead - Factor jstart columns: panel[0..depth-1]
 */
   for( k = 0, j = 0; k < depth; k++ )
   {
      jb = jstart - j; jb = Mmin( jb, nb ); j += jb;
/*
 * Factor and broadcast k-th panel - use long topology for those
 */
      HPL_pzfact(         panel[k] );
      (void) HPL_zbinit(   panel[k] );
      do
      { (void) HPL_zbcast( panel[k], &test ); }
      while( test != HPL_SUCCESS );
      (void) HPL_zbwait(   panel[k] );
/*
 * Partial update of the depth-1-k panels in front of me
 */
      if( k < depth - 1 )
      {
         nn = HPL_numrocI( jstart-j, j, nb, nb, mycol, 0, npcol );
         HPL_pzupdate( NULL, NULL, panel[k], nn );
      }
   }
/*
 * Main loop over the remaining columns of A
 */
   for( j = jstart; j < N; j += nb )
   {
      n = N - j; jb = Mmin( n, nb );
/*
 * Allocate current panel resources - Finish latest update - Factor and
 * broadcast current panel
 */
      HPL_pzpanel_new( GRID, ALGO, n, n+1, jb, A, j, j, tag, &panel[depth] );
 
      if( mycol == icurcol )
      {
         nn = HPL_numrocI( jb, j, nb, nb, mycol, 0, npcol );
         for( k = 0; k < depth; k++ )   /* partial updates 0..depth-1 */
            HPL_pzupdate( NULL, NULL, panel[k], nn );
         HPL_pzfact(    panel[depth] );       /* factor current panel */
      }
      else { nn = 0; }
          /* Finish the latest update and broadcast the current panel */
      (void) HPL_zbinit( panel[depth] );
      HPL_pzupdate(     panel[depth], &test, panel[0], nq-nn );
      (void) HPL_zbwait( panel[depth] );
/*
 * Release latest panel resources - circular  of the panel pointers
 * Go to the next process row and column -  update  the message ids  for
 * broadcast
 */
      (void) HPL_pzpanel_disp( &panel[0] );
      for( k = 0; k < depth; k++ ) panel[k] = panel[k+1];
 
      if( mycol == icurcol ) { jj += jb; nq -= jb; }
      icurcol = MModAdd1( icurcol, npcol );
      tag     = MNxtMgid( tag, MSGID_BEGIN_FACT, MSGID_END_FACT );
   }
/*
 * Clean-up: Finish updates - release panels and panel list
 */
   nn = HPL_numrocI( 1, N, nb, nb, mycol, 0, npcol );
   for( k = 0; k < depth; k++ )
   {
      HPL_pzupdate( NULL, NULL, panel[k], nn );
      (void) HPL_pzpanel_disp( &panel[k] );
   }
 
   if( panel ) free( panel );
/*
 * End of HPL_pzgesvK1
 */
}
