/*
 Fortran callable interface
*/
#define DLEN_ 9
#define CTXT_ 1
#define M_ 2
#define N_ 3
#define NB_ 5
#define MB_ 4
#define LLD_ 8
#define RSRC_ 6
#define ALIGNMENT 8

#include "zhpl.h"
#include "mpi.h"
#include <assert.h>


#ifdef Add__

#define BLACS_GET blacs_get__
#define BLACS_GRIDINFO blacs_gridinfo__

#define HPL_BLACSINIT_REAL hpl_zblacsinit__
#define HPL_BLACSINIT_LAYER hpl_zblacsinit_ 

#define HPL_PZGESV_REAL hpl_pzgesv__
#define HPL_PZGESV_LAYER hpl_pzgesv_

#define HPL_MATINIT_REAL hpl_zmatinit__
#define HPL_MATINIT_LAYER hpl_zmatinit_

#else

#define BLACS_GET blacs_get_
#define BLACS_GRIDINFO blacs_gridinfo_


#define HPL_BLACSINIT_REAL hpl_zblacsinit_ 
#define HPL_BLACSINIT_LAYER hpl_zblacsinit__ 

#define HPL_PZGESV_REAL hpl_pzgesv_ 
#define HPL_PZGESV_LAYER hpl_pzgesv__

#define HPL_MATINIT_REAL hpl_zmatinit_ 
#define HPL_MATINIT_LAYER hpl_zmatinit__
#endif


static HPL_T_grid grid ;
static HPL_ZT_pmat mat ;



extern void pzswap_( int *n, 
             zcomplex *x, int *ix, int *jx, int *descX, int *incx,
             zcomplex *y, int *iy, int *jy, int *descY, int *incy );

extern void pzcopy_( int *n, 
             zcomplex *x, int *ix, int *jx, int *descX, int *incx,
             zcomplex *y, int *iy, int *jy, int *descY, int *incy );






extern void BLACS_GET (int *ConTxt, int *what, int *val);

extern void BLACS_GRIDINFO ( int *ConTxt, int *nprow, int *npcol, 
                        int *myprow, int *mypcol );

extern void hpl_blacsinit_(int *blacs_context );
extern void hpl_blacsinit__(int *blacs_context );

extern void hpl_pzgesv_(  int *pn, zcomplex *A, int *descA, int *ipiv, int *info );
extern void hpl_pzgesv__(  int *pn, zcomplex *A, int *descA, int *ipiv, int *info );






void HPL_BLACSINIT_REAL( int *blacs_context) 
{

  int what = 10;
  int ivalue = 0;
  MPI_Comm mpi_comm = MPI_COMM_WORLD;
  MPI_Comm new_mpi_comm = MPI_COMM_WORLD;

  int nprow = 1;
  int npcol = 1;
  int myprow = 0;
  int mypcol = 0;

  int nprow_hpl = 1;
  int npcol_hpl = 1;
  int myprow_hpl = 0;
  int mypcol_hpl = 0;

  int irank = 0;
  int isvalid = 0;
  int use_column_order = 0;

  int color = -1;
  int key = 0;

  int mpierr = MPI_SUCCESS;
  int ierr = MPI_SUCCESS;

  HPL_T_ORDER order = HPL_COLUMN_MAJOR;

  BLACS_GRIDINFO ( blacs_context, &nprow, &npcol, &myprow, &mypcol );
  assert( nprow >= 1 );
  assert( npcol >= 1 );
  assert( (0 <= myprow) && (myprow < nprow) );
  assert( (0 <= mypcol) && (mypcol < npcol) );

  what = 10;
  BLACS_GET ( blacs_context, &what, &ivalue );

  mpi_comm = (MPI_Comm) ivalue;


  use_column_order =  (order == HPL_COLUMN_MAJOR);
  if (use_column_order) {
    irank = myprow + mypcol*nprow;
    }
  else {
    irank = mypcol + myprow*npcol;
    };

  color = MPI_UNDEFINED;
  isvalid = (0 <= myprow) && (myprow < nprow) && (nprow >= 1) &&
            (0 <= mypcol) && (mypcol < npcol) && (npcol >= 1);
  if (isvalid) {
      color = 0;
      };

  /*
   setup new mpi communicator with the exact mapping of 
   processors
  */


  key = irank;
  mpierr =  MPI_Comm_split(  mpi_comm, color, key, &new_mpi_comm );
  assert( mpierr == MPI_SUCCESS );


  ierr = HPL_grid_init( new_mpi_comm, order, nprow, npcol, &grid );
  assert( ierr == MPI_SUCCESS );



  /* 
   zcomplex check
  */

  ierr = HPL_grid_info( &grid, &nprow_hpl, &npcol_hpl, 
                               &myprow_hpl, &mypcol_hpl );
  assert( ierr == MPI_SUCCESS );

  assert( nprow_hpl == nprow );
  assert( npcol_hpl == npcol );
  assert( myprow_hpl == myprow );
  assert( mypcol_hpl == mypcol );

}


/*
  compute ld and amount of memory to match what HPL wants
  code taken from HPL_pztest.c
*/

void HPL_MATINIT_REAL( int *pn, int *pnb,    int *ld, int *ineed )
{
   HPL_ZT_palg  ALGO_t;
   HPL_ZT_palg *ALGO = &ALGO_t;
   HPL_ZT_pmat mat;

   int nq = 0;
   int ip2 = 0;
   int ii = 0;
   int N = *pn;
   int NB = *pnb;

   int ierr = 0;
   int myrow = 0;
   int mycol = 0;
   int nprow = 1;
   int npcol = 1;

  ierr = HPL_grid_info( &grid, &nprow, &npcol, 
                               &myrow, &mycol );
  assert( ierr == MPI_SUCCESS );



   ALGO->align = ALIGNMENT;

   mat.n  = N; mat.nb = NB; mat.info = 0;
   mat.mp = HPL_numroc( N, NB, NB, myrow, 0, nprow );
   nq     = HPL_numroc( N, NB, NB, mycol, 0, npcol );
   mat.nq = nq + 1;

/*
 * Allocate matrix, right-hand-side, and vector solution x. [ A | b ] is
 * N by N+1.  One column is added in every process column for the solve.
 * The  result  however  is stored in a 1 x N vector replicated in every
 * process row. In every process, A is lda * (nq+1), x is 1 * nq and the
 * workspace is mp. 
 *
 * Ensure that lda is a multiple of ALIGN and not a power of 2
 */
   mat.ld = ( ( Mmax( 1, mat.mp ) - 1 ) / ALGO->align ) * ALGO->align;
   do
   {
      ii = ( mat.ld += ALGO->align ); ip2 = 1;
      while( ii > 1 ) { ii >>= 1; ip2 <<= 1; }
   }
   while( mat.ld == ip2 );


   *ineed = (ALGO->align + (mat.ld+1)*(mat.nq)) + N;
   *ld = mat.ld;

}

void HPL_MATINIT_LAYER( int *pn, int *pnb,  int *ld, int *ineed ) 
{
  HPL_MATINIT_REAL( pn, pnb, ld, ineed );
}

  
void HPL_BLACSINIT_LAYER( int *blacs_context )
{
  HPL_BLACSINIT_REAL( blacs_context );
}


void HPL_PZGESV_REAL(  int *pn, zcomplex *A, int *descA, int *ipiv, int *info )
{


   int                        nval  [HPL_MAX_PARAM],
                              nbval [HPL_MAX_PARAM],
                              pval  [HPL_MAX_PARAM],
                              qval  [HPL_MAX_PARAM],
                              nbmval[HPL_MAX_PARAM],
                              ndvval[HPL_MAX_PARAM],
                              ndhval[HPL_MAX_PARAM];

   HPL_T_FACT                 pfaval[HPL_MAX_PARAM],
                              rfaval[HPL_MAX_PARAM];

   HPL_T_TOP                  topval[HPL_MAX_PARAM];

   HPL_ZT_palg                 algo;
   HPL_T_test                 test;
   int                        L1notran, Unotran, align, equil,  
                              inbm, indh, indv, ipfa,  irfa, itop,
                              ns, nbs, nbms, ndhs, ndvs,
                              npfs, npqs, nrfs, ntps, 
                              rank, size, tswap;
   HPL_T_ORDER                pmapping;
   HPL_T_FACT                 rpfa;
   HPL_T_SWAP                 fswap;


int mpierr = MPI_SUCCESS;
int ierr = MPI_SUCCESS;


int nq = 0;
int n = *pn;


int nb = descA[NB_];
int ld = descA[LLD_];


int nprow = 0;
int npcol = 0;
int myrow = 0;
int mycol = 0;

int isroot = 0;
int j = 0;



ierr = HPL_grid_info( &grid, &nprow, &npcol, &myrow, &mycol );
assert( ierr == MPI_SUCCESS );
assert( (0 <= myrow) && (myrow < nprow) );
assert( (0 <= mycol) && (mycol < npcol) );



isroot = (myrow == 0) && (myrow == 0);


mat.info = 0;
mat.n = n;
mat.nb = nb;
mat.ld = ld;
mat.mp = HPL_numroc( n, nb, nb, myrow, 0, nprow );
nq = HPL_numroc( n, nb, nb, mycol, 0, npcol );
mat.nq = nq + 1;
mat.A = A;
mat.X = mat.A + (mat.ld * mat.nq);

/*
 *Debug: allocate the pivot vector
 */

 // mat.global_ipiv = (int *)malloc( n * sizeof( int ) );
 mat.global_ipiv = (int *) (A + (ALIGNMENT + (mat.ld+1)*(mat.nq)));
 {
	int k;
	for (k = 0; k < n; k++)
		mat.global_ipiv[k] = -1;
 }
 assert( (1 <= mat.ld) && (mat.mp <= mat.ld) );



/*
 * Set up the algorithm parameters
 */


   mpierr = MPI_Comm_rank( MPI_COMM_WORLD, &rank );
   assert( mpierr == MPI_SUCCESS );

   mpierr = MPI_Comm_size( MPI_COMM_WORLD, &size );
   assert( mpierr == MPI_SUCCESS );

/*
HPLinpack benchmark input file
Innovative Computing Laboratory, University of Tennessee
HPL.out      output file name (if any)
8            device out (6=stdout,7=stderr,file)
1            # of problems sizes (N)
1064280       Ns
1            # of NBs
60           NBs
0            PMAP process mapping (0=Row-,1=Column-major)
1            # of process grids (P x Q)
100            Ps
104             Qs
16.0         threshold
1            # of panel fact
2            PFACTs (0=left, 1=Crout, 2=Right)
1            # of recursive stopping criterium
4            NBMINs (>= 1)
1            # of panels in recursion
2            NDIVs
1            # of recursive panel fact.
2            RFACTs (0=left, 1=Crout, 2=Right)
1            # of broadcast
1            BCASTs (0=1rg,1=1rM,2=2rg,3=2rM,4=Lng,5=LnM)
1            # of lookahead depth
1            DEPTHs (>=0)
0            SWAP (0=bin-exch,1=long,2=mix)
64           swapping threshold
1            L1 in (0=transposed,1=no-transposed) form
1            U  in (0=transposed,1=no-transposed) form
0            Equilibration (0=no,1=yes)
8            memory alignment in zcomplex (> 0)
##### This line (no. 32) is ignored (it serves as a separator). ######
0                      		Number of additional problem sizes for PTRANS
1200 10000 30000        	values of N
10                       	number of additional blocking sizes for PTRANS
30 31 32 33 63 64 65 66 67 68                       values of NB
 */
 ipfa = 0;
 inbm = 0;
 indv = 0;

 irfa = 0;
 itop = 0;
 indh = 0;

#undef USE_HPL_PZINFO

#ifdef USE_HPL_PZINFO
   printf("before HPL_pzinfo\n");
   HPL_pzinfo( &test, &ns, nval, &nbs, nbval, &pmapping, &npqs, pval, qval,
               &npfs, pfaval, &nbms, nbmval, &ndvs, ndvval, &nrfs, rfaval,
               &ntps, topval, &ndhs, ndhval, &fswap, &tswap, &L1notran,
               &Unotran, &equil, &align );



#else
/*
pfaval  == PF  panel factorization algorithm
nbmval  == NBM recursive stopping critera
ndvval  == NDV number of panels in recursion
rfaval  == RF  recursive factorization algorithms 
topval  == TP  broadcast (along rows) topologies
ndhval  == DH  lookahead depths 0 is no-lookahead
*/


/*
1            # of panel fact
2            PFACTs (0=left, 1=Crout, 2=Right)
*/
         j = 2;
         if(      j == 0 ) pfaval[ ipfa ] = HPL_LEFT_LOOKING;
         else if( j == 1 ) pfaval[ ipfa ] = HPL_CROUT;
         else if( j == 2 ) pfaval[ ipfa ] = HPL_RIGHT_LOOKING;
         else              pfaval[ ipfa ] = HPL_RIGHT_LOOKING;


/*
1            # of recursive stopping criterium
4            NBMINs (>= 1)
*/
nbmval[inbm] = 4;



/*
1            # of panels in recursion
2            NDIVs
*/
ndvval[indv] = 2;


/*
1            # of recursive panel fact.
2            RFACTs (0=left, 1=Crout, 2=Right)
*/
         j = 2;
         if(      j == 0 ) rfaval[ irfa ] = HPL_LEFT_LOOKING;
         else if( j == 1 ) rfaval[ irfa ] = HPL_CROUT;
         else if( j == 2 ) rfaval[ irfa ] = HPL_RIGHT_LOOKING;
         else              rfaval[ irfa ] = HPL_RIGHT_LOOKING;


/*
1            # of broadcast
1            BCASTs (0=1rg,1=1rM,2=2rg,3=2rM,4=Lng,5=LnM)
*/

         j = 1;
         if(      j == 0 ) topval[ itop ] = HPL_1RING;
         else if( j == 1 ) topval[ itop ] = HPL_1RING_M;
         else if( j == 2 ) topval[ itop ] = HPL_2RING;
         else if( j == 3 ) topval[ itop ] = HPL_2RING_M;
         else if( j == 4 ) topval[ itop ] = HPL_BLONG;
         else if( j == 5 ) topval[ itop ] = HPL_BLONG_M;
         else              topval[ itop ] = HPL_1RING_M;

/*
1            # of lookahead depth
1            DEPTHs (>=0)
*/
ndhval[indh] = 1;



/*
0            SWAP (0=bin-exch,1=long,2=mix)
*/
      j = 0;
      if(      j == 0 ) fswap = HPL_SWAP00;
      else if( j == 1 ) fswap = HPL_SWAP01;
      else if( j == 2 ) fswap = HPL_SW_MIX;
      else              fswap = HPL_SWAP01;


/*
64           swapping threshold
1            L1 in (0=transposed,1=no-transposed) form
1            U  in (0=transposed,1=no-transposed) form
0            Equilibration (0=no,1=yes)
8            memory alignment in zcomplex (> 0)
*/

tswap = 64;
L1notran =  1;
Unotran = 1;
equil = 0;

align = ALIGNMENT;
#endif


#ifdef DEBUG
   if (isroot) {
   printf("pfaval[ipfa] %d \n", pfaval[ipfa]);
   printf("nbmval[inbm] %d \n", nbmval[inbm]);
   printf("ndvval[indv] %d \n", ndvval[indv]);

   printf("rfaval[irfa] %d \n", rfaval[irfa]);
   printf("topval[itop] %d \n", topval[itop]);
   printf("ndhval[indh] %d \n", ndhval[indh]);


   printf("fswap %d \n", fswap);
   printf("tswap %d \n", tswap);
   printf("L1notran %d \n", L1notran);
   printf("Unotran %d \n", Unotran);
   printf("equil %d \n", equil);
   };
#endif




              algo.btopo = topval[itop]; algo.depth = ndhval[indh];
              algo.nbmin = nbmval[inbm]; algo.nbdiv = ndvval[indv];

              algo.pfact = rpfa = pfaval[ipfa];

              if( L1notran != 0 )
              {
                 if( rpfa == HPL_LEFT_LOOKING ) algo.pffun = HPL_pzpanllN;
                 else if( rpfa == HPL_CROUT   ) algo.pffun = HPL_pzpancrN;
                 else                           algo.pffun = HPL_pzpanrlN;

                 algo.rfact = rpfa = rfaval[irfa];
                 if( rpfa == HPL_LEFT_LOOKING ) algo.rffun = HPL_pzrpanllN;
                 else if( rpfa == HPL_CROUT   ) algo.rffun = HPL_pzrpancrN;
                 else                           algo.rffun = HPL_pzrpanrlN;

                 if( Unotran != 0 ) algo.upfun = HPL_pzupdateNN;
                 else               algo.upfun = HPL_pzupdateNT;
              }
              else
              {
                 if( rpfa == HPL_LEFT_LOOKING ) algo.pffun = HPL_pzpanllT;
                 else if( rpfa == HPL_CROUT   ) algo.pffun = HPL_pzpancrT;
                 else                           algo.pffun = HPL_pzpanrlT;

                 algo.rfact = rpfa = rfaval[irfa];
                 if( rpfa == HPL_LEFT_LOOKING ) algo.rffun = HPL_pzrpanllT;
                 else if( rpfa == HPL_CROUT   ) algo.rffun = HPL_pzrpancrT;
                 else                           algo.rffun = HPL_pzrpanrlT;

                 if( Unotran != 0 ) algo.upfun = HPL_pzupdateTN;
                 else               algo.upfun = HPL_pzupdateTT;
              }

              algo.fswap = fswap; algo.fsthr = tswap;
              algo.equil = equil; algo.align = align;
 

/*
 avoid complications with alignment
 set alignment to 1 
*/
  algo.align = ALIGNMENT;


#ifdef DEBUG
/*
 print out matrix
*/

  { 
    int ia = 0;
    int ja = 0;
    zcomplex aij = 0;

    if (isroot) {
      printf("n %d ld %d nb %d \n",
        mat.n, mat.ld, mat.nb );
      printf("mat.mp %d mat.nq %d \n",
        mat.mp, mat.nq );
      };

    for (ja=0; ja < mat.n; ja++) {
    for (ia=0; ia < mat.n; ia++ ) {
      HPL_pzelget( &grid, mat.A, ia,ja, mat.ld, mat.nb, &aij );
      if (isroot) {
        printf("ia %d ja %d aij %lf \n", 
                ia,ja,aij );
        };
      };
      };

    for (ia=0; ia < mat.n; ia++) {
      ja = mat.n ;
      HPL_pzelget( &grid, mat.A, ia, ja, mat.ld, mat.nb, &aij );
      if (isroot) {
        printf("ia %d ja %d Xi %lf  \n", 
                ia, ja, aij );
        };
      };
  }
 
#endif
  
  HPL_pzgesv( &grid, &algo, &mat );


/*
 * Debug: make the global pivot vector and do row swap
 */
{
  int k = 0;
  int nn = 0, L = 0;
  ( void ) HPL_zall_reduce( (int *)(mat.global_ipiv), 
                 mat.n, HPL_INT, HPL_zmax, grid.all_comm );
  
 /* debug
  for (k = 0; k < mat.n; k++) {
	  printf("global_ipiv[%d] = %d\n", k, mat.global_ipiv[k]);
          };
  */
  
  for (k = 0; k < mat.n; k++) {
	  nn = ( (int) k / mat.nb ) * mat.nb;
	  L = mat.global_ipiv[k];
	  if ( L != k ) {
		 int Lp1 = L + 1, kk = k + 1, one = 1;

                 // printf("Lp1 %d kk %d \n", Lp1, kk );

		 pzswap_( &nn, mat.A, &kk, &one, descA, 
                    &(descA[M_]), mat.A, &Lp1, &one, descA, &(descA[M_]) );
  	  }
  }
 

}


	
/*
 * Debug: Take the local pivot vector from global_ipiv
 */

//Method 1
{
	int i = 0;
	int row_processor = 0;
	int is_mine = 0;
	int myrow = grid.myrow;
	int nprow = grid.nprow;
	int nprocs = grid.nprocs;
	int nb = mat.nb;
	int *global_ipiv = mat.global_ipiv;
	int n = mat.n;
	int lrow;

	int srcproc = descA[RSRC_];
        int mb = descA[MB_]; 
        int ifirst = n+1;
        int istart = n+1;
        int iend = -1;
	int roc = HPL_numroc(mat.n, mb, mb, grid.myrow, 0, grid.nprow);




        /* -------------------------------------------------  */
        /* find ifirst, the first row index on this processor */
        /* -------------------------------------------------  */
        ifirst = n+1;
        for(i=0; i < n; i += mb) {
	   HPL_indxg2lp( &lrow, &row_processor, i, mb, mb, srcproc, nprow );
           is_mine = (row_processor == myrow);
           if (is_mine) {
              ifirst = i; 
              break;
              };
           };
        /* ----------------------------------------------- */
        /* take advantage of the 2D block cyclic structure */
        /* ----------------------------------------------- */
        for(istart=ifirst; istart < n; istart += mb*nprow) {
           int loff = 0;

	   HPL_indxg2lp(&loff,&row_processor,istart,mb,mb,srcproc,nprow);
           iend = istart + mb;
           if (iend > n) {
               iend = n;
               };
           for(i=istart; i < iend; i++) {
#ifdef DEBUG
	       HPL_indxg2lp(&lrow,&row_processor,i,mb,mb,srcproc,nprow);
               assert( row_processor == myrow );
               assert( lrow == (loff + (i-istart)) );
#endif
               ipiv[ loff + (i-istart) ] = global_ipiv[i] + 1;
               };
            };
              
               
              
          

          

#if (0)
	for (i = 0; i < mat.n; i++) {
		HPL_indxg2lp( &lrow, &row_processor, i, nb, nb, srcproc, nprow );
		
		is_mine = ( row_processor == myrow );
		if ( is_mine ) { 
//			printf("local row = %d; global row = %d; row_processor = %d\n", lrow, i, row_processor);
			ipiv[lrow] = global_ipiv[i] + 1;
			}
	}
#endif
//	for (i = 0; i < mat.n; i++) 
//		printf("global_ipiv[%d] = %d\n", i, global_ipiv[i]);
//	for (i = 0; i < roc; i++)
//		printf("ipiv[%d] = %d\n", i, ipiv[i]);
	
}

/*  Method 3 (Problem with the pzcopy_)
{
	int roc = HPL_numroc(mat.n, mat.nb, mat.nb, grid.myrow, 0, grid.nprow);
	printf("roc = %d\n", roc);
	
	int desc_global_ipiv[DLEN_];
	int desc_ipiv[DLEN_];
	zcomplex ipiv[roc];
	zcomplex global[mat.n];
	int m, n, mb, nb, rsrc, csrc,lld, ictxt;

	//setup descriptor for global_ipiv

	m = descA[M_];
	mb = m;
	n = 1;
	nb = 1;
	rsrc = -1;
	csrc = -1;
	lld = m;
	ictxt =descA[CTXT_];

	descset_(desc_global_ipiv, &m, &n, &mb, &nb, &rsrc, &csrc, &ictxt, &lld);
	
	//setip descriptor for ipiv

	m = descA[M_];
	mb = descA[MB_];
	n = 1;
	nb = 1;
	rsrc = descA[RSRC_];
	csrc = -1;
	lld = descA[LLD_];
	ictxt = descA[CTXT_];
	
	descset_(desc_ipiv, &m, &n, &mb, &nb, &rsrc, &csrc, &ictxt, &lld);

//	int kk;
//	for (kk = 0; kk < mat.n; kk++) {
//		global[kk] = (zcomplex)mat.global_ipiv[kk];
//		printf("global[%d] = %f\n", kk, global[kk]);
//	}
	
	int one = 1;
	pzcopy_(&one, global, &one, &one, desc_global_ipiv, &one, ipiv, &one, &one, desc_ipiv, &one);

//	for (kk = 0; kk < roc; kk++)
//		printf("ipiv[%d] = %f\n", kk, ipiv[kk]);

}
*/

/*
 * Debug: free the memories
 */
  // free( mat.global_ipiv );
  
#ifdef DEBUG
/*
 print content of X
 */
  {
   int ix = 0; 
   int jx = 0;
   int ia = 0;
   int ja = 0;
   int i = 0;
   int ld = 1;
   zcomplex xij = 0;
   zcomplex aij = 0;

   for (i=0; i < mat.n; i++) {
       ix = 0; jx = i;
       HPL_pzelget( &grid, mat.X, ix,jx, ld, mat.nb, &xij );

       ia = i; ja = mat.n;
       HPL_pzelget( &grid, mat.A, ia,ja, mat.ld, mat.nb, &aij );
       if (isroot) {
          printf("jx %d xij %lf aij %lf \n", 
              jx, xij, aij );
          };
       /* HPL_pzelset( &grid, mat.A, ia,ja, mat.ld, mat.nb, xij ); */
       };

  }
#endif

/*
  copy X to A
*/
  {
  int i = 0;
  int ldx = 1;
  int inc1 = ldx;
  int ix = 1;
  int jx = 1;

  int ia = 1;
  int ja = descA[M_]+1;
  int inc2 = 1;

  int descX[DLEN_];

  for(i=0; i < DLEN_; i++) {
      descX[i] = descA[i];
      };
  descX[M_] = 1;
  descX[N_] = mat.n;
  descX[LLD_] = ldx;

  pzcopy_( &(mat.n), mat.X, &ix,&jx, descX, &inc1, 
                    mat.A, &ia,&ja, descA, &inc2);
  }


  *info = mat.info; 
  if (mat.info == MPI_SUCCESS)  {
     *info = 0;
     };

}

void HPL_PZGESV_LAYER( int *pn, zcomplex *A, int *descA, int *ipiv, int *info )
{
  HPL_PZGESV_REAL( pn, A, descA, ipiv, info );
}
