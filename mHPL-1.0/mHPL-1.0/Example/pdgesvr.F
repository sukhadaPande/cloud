        subroutine pdgesvr( n, nrhs, A, ia, ja, descA,                     &
     &          ipiv, B, ib, jb, descB, info )
        integer, intent(in) :: n,nrhs, ia,ja, ib,jb
        integer, intent(in) :: descA(DLEN_), descB(DLEN_)
        integer, intent(inout) :: info

        integer, intent(inout) :: ipiv(*)
        real(kind=dp), intent(inout) :: A(*)
        real(kind=dp), intent(inout) :: B(*)
*       ---------------
*       local variables
*       ---------------
        logical :: isok, isdone, isroot

        integer :: descAr(DLEN_), descX(DLEN_)
        integer :: ictxt, myrow,mycol, nprow,npcol
        integer :: Locp, Locq, lrindx, lcindx, lld, rsrc,csrc
        integer :: Locp_B, Locq_B


        integer :: iter
        integer, parameter :: maxiter = 30

        integer :: Xsize, Asize, ierr
        real(kind=sp), dimension(:), allocatable :: Ar, Xr
        real(kind=dp), dimension(:), allocatable :: X, R
        real(kind=dp)  :: anorm, eps, bwdcte, xnorm
        real(kind=dp)  :: rnorm,rnorm_old,bwdmax

        real(kind=sp), dimension(:), allocatable :: rwork
        real(kind=dp), dimension(:), allocatable :: dwork
        integer :: ldwork 
        integer :: m,Mp0,Np0

        real(kind=dp), parameter :: one = 1.0d0
        real(kind=dp) :: alpha,beta

#ifdef USE_HPL
        integer :: hpl_lld, hpl_ineed
#endif

*       ----------------------------------------------
*       norm = 'F','f','E', or 'e' for Frobenius norm
*       norm = 'M','m' for max norm
*       norm = '1','O','o' for 1-norm
*       norm = 'I','i', for infinity norm
*       ----------------------------------------------
        character :: normtype = 'I'   

        ictxt = descA(CTXT_)
        call blacs_gridinfo( ictxt, nprow,npcol,myrow,mycol )
        isroot = (myrow.eq.0).and.(mycol.eq.0)

        m = n
        Mp0 = iceil( m, descA(MB_) )*descA(MB_)
        Np0 = iceil( n, descA(NB_) )*descA(NB_)
        ldwork = max(Mp0,Np0) + 1
        allocate( dwork(ldwork),rwork(ldwork),stat=ierr)
        if (ierr.ne.0) then
           write(*,*) 'allocate(dwork) return ierr=',ierr
        endif

         anorm = pdlange(normtype, n,n, A,ia,ja,descA, dwork )
         eps = dlamch('P')
         bwdmax = 1.0d0
*        bwdcte = anorm * eps * min(4.0d0,sqrt(dble(n))/6.0d0)
*        bwdcte = anorm * eps * max(4.0d0,sqrt(dble(n))/6.0d0)
         bwdcte = anorm * eps * sqrt( dble(n) ) * bwdmax

      

*       ----------------------------------
*       make a copy with reduced precision
*       ----------------------------------

*       ----------------------------------------------------
*       set rsrc,csrc so data is exactly aligned
*       on the same processor, no communication is required
*       ----------------------------------------------------
        call infog2l( ia,ja, descA, nprow,npcol, myrow,mycol,               &
     &         lrindx, lcindx, rsrc, csrc )

*       ---------------------------------------
*       for simplicity, use (rsrc,csrc) = (0,0)
*       ---------------------------------------
        rsrc = 0
        csrc = 0

        Locp = numroc( n, descA(MB_), myrow, rsrc, nprow)
        Locq = numroc( n+1, descA(NB_), mycol, csrc, npcol)
        lld = max(1,Locp)
        Asize = max(1,Locp*Locq)


#ifdef USE_HPL
        call hpl_dmatinit( n, descA(MB_), hpl_lld, hpl_ineed )

        if (idebug .ge. 2) then
           write(*,*) 'myrow,mycol,hpl_lld,hpl_ineed ',                     &
     &              myrow,mycol,hpl_lld,hpl_ineed
           write(*,*) 'n,mb,nb ', n, descA(MB_),descA(NB_)
        endif

        call assert(hpl_lld.ge.lld, 'hpl_lld < lld ', hpl_lld )
        call assert(hpl_ineed.ge.Asize,'hpl_ineed < Asize ',hpl_ineed)

        lld = hpl_lld
        Asize = hpl_ineed
#endif

        allocate( Ar( Asize ), stat=ierr)
        if (ierr.ne.0) then
          write(*,*) 'allocate(Ar) return ierr=',ierr
        endif
        Ar(:) = 0
     

        call descinit( descAr, n,n+1, descA(MB_), descA(NB_),                 &
     &       rsrc, csrc, descA(CTXT_), lld, info )
        isok = (info.eq.0)
        if (.not.isok) then
                write(*,*) 'decinit for Ar return info=',info
                stop '** error in pdgesvr '
        endif

#ifdef USE_PROF
        call blacs_barrier( descAr(CTXT_),'A')
        call profstart('pzgesvr: pconvertd2s A')
#endif

        
        call pconvertd2s( n,n, A,ia,ja,descA, Ar,1,1,descAr )

#ifdef USE_PROF
        call blacs_barrier( descAr(CTXT_),'A')
        call profend('pzgesvr: pconvertd2s A')
#endif


*       ------------------
*       make a copy of rhs
*       ------------------
        call infog2l( ib,jb, descB, nprow,npcol, myrow,mycol,              &
     &          lrindx, lcindx, rsrc, csrc )
        Locp_B = numroc( n, descB(MB_), myrow, rsrc, nprow)
        Locq_B = numroc( nrhs, descB(NB_), mycol, csrc, npcol)
        Xsize = max(1, Locp_B*Locq_B)
        allocate( X( Xsize ), R(Xsize), Xr(Xsize), stat=ierr )
        if (ierr.ne.0) then
          write(*,*) 'allocate(X) return ierr=',ierr
        endif

*       ----------------
*       setup descriptor
*       ----------------
        lld = max(1, Locp_B )
        call descinit( descX, n, nrhs,  descB(MB_), descB(NB_),            &
     &            rsrc, csrc, descB(CTXT_), lld, info )


        isok = (info.eq.0)
        if (.not.isok) then
                write(*,*) 'descinit for X return info=',info
                stop '** error in pdgesvr '
        endif

        call pdgecopy( n,nrhs, B,ib,jb,descB, X,1,1,descX )

        call pconvertd2s( n,nrhs, X,1,1,descX, Xr,1,1,descX )

        if ((idebug.ge.2) .and. isroot) then
*       ------------
*       double check
*       ------------
        rnorm = pslange(normtype,n,nrhs,Xr,1,1,descX,rwork)
        write(*,*) 'norm(Xr) after c2z ', rnorm


        rnorm = pdlange(normtype, n,nrhs, B,ib,jb,descB, dwork )
        write(*,*) 'norm(B) ', rnorm
        call pconverts2d( n,nrhs, Xr,1,1,descX, R,1,1,descX )


        rnorm = pdlange(normtype,n,nrhs,R,1,1,descX,dwork)
        write(*,*) 'norm(R) after c2z ', rnorm

        alpha = -1.0d0
        beta = 1.0d0
        call pdgeadd( 'N', n,nrhs, alpha, X,1,1,descX,                     &
     &                 beta,R,1,1,descX)
        rnorm = pdlange( normtype, n,nrhs, R,1,1,descX,dwork)
        write(*,*) 'rnorm from conversion ',rnorm

        endif
*       ---------------------------------------------
*       perform LU factorization in reduced precision
*       ---------------------------------------------

#ifdef USE_HPL

#ifdef USE_PROF
        call blacs_barrier(ictxt,'A')
        call profstart('hpl_psgesv')
#endif

        call hpl_psgesv(n, Ar, descAr,ipiv,info)
        call assert(info.eq.0,'hpl_psgesv return info',info)
        
#ifdef USE_PROF
        call blacs_barrier(ictxt,'A')
        call profend('hpl_psgesv')
#endif

#else




#ifdef USE_PROF
        call blacs_barrier(ictxt,'A')
        call profstart('psgetrf')
#endif
        call psgetrf( n,n, Ar, 1,1, descAr, ipiv, info )
#ifdef USE_PROF
        call blacs_barrier(ictxt,'A')
        call profend('psgetrf')
#endif


#endif

        isok = (info.eq.0)
        if (.not.isok) then
*          ---------------------------------
*          system may be too ill conditioned
*          ---------------------------------
           info  = -(n+1)
           goto 999
        endif
*       -----------------------
*       generate  initial guess
*       -----------------------
        call pconvertd2s( n,nrhs, B,ib,jb,descB, Xr,1,1,descX )

#ifdef USE_PROF
        call blacs_barrier(ictxt,'A')
        call profstart('psgetrs')
#endif
        call psgetrs( 'N', n, nrhs, Ar, 1,1, descAr, ipiv,                 &
     &        Xr, 1,1, descX, info )
#ifdef USE_PROF
        call blacs_barrier(ictxt,'A')
        call profend('psgetrs')
#endif


        call pconverts2d( n,nrhs, Xr,1,1,descX,X,1,1,descX)

        if ((idebug.ge.2).and.isroot) then
                write(*,*) 'initial B(1:5) ', B(1:5)
                write(*,*) 'initial Xr(1:5) ', Xr(1:5)
                write(*,*) 'initial X(1:5) ', X(1:5)
        endif

*       ----------------        
*       compute residual
*       ----------------        
        call pdgecopy( n,nrhs, B, ib,jb, descB, R, 1,1,descX )

        alpha = -one
        beta = one

#ifdef USE_PROF
        call profstart('pzgesvr:pzgemm')
#endif
        call pdgemm('N','N',n,nrhs,n, alpha, A,ia,ja,descA,                 &
     &        X,1,1,descX,  beta, R, 1,1, descX )
#ifdef USE_PROF
        call profend('pzgesvr:pzgemm')
#endif

        xnorm = pdlange(normtype,n,nrhs,B,ib,jb,descB,dwork)
        rnorm_old = xnorm
        do iter=1,maxiter
*          ---------------------
*          check for convergence
*          ---------------------
           rnorm =  pdlange( normtype, n,nrhs, R,1,1,descX, dwork )
           isdone = (rnorm.le.xnorm*bwdcte) .and. (iter.ge.6)

           if ((idebug.ge.1).and.isroot) then
             write(*,*) 'iter, rnorm ', iter,rnorm
           endif

           if (isdone) then
                   exit
           endif

           if ((iter.ge.7).and.(rnorm.gt.rnorm_old)) then
              exit
           endif
*          ------------------
*          compute correction
*          ------------------
           call pconvertd2s( n,nrhs, R,1,1, descX, Xr,1,1,descX )
#ifdef USE_PROF
           call profstart('psgetrs')
#endif
           call psgetrs( 'N', n,nrhs, Ar, 1,1, descAr, ipiv,               &
     &              Xr,1,1,descX,info )
#ifdef USE_PROF
           call profend('psgetrs')
#endif
*          ----------------
*          update solution
*          ----------------
           call pconverts2d( n,nrhs, Xr, 1,1,descX, R,1,1,descX )
           if ((idebug.ge.2) .and.(isroot)) then
             rnorm = pslange(normtype,n,nrhs, Xr,1,1,descX,rwork)
             write(*,*) 'norm(Xr correction) ', rnorm

             rnorm = pdlange(normtype,n,nrhs, R,1,1,descX,dwork)
             write(*,*) 'norm( correction ) ',  rnorm
             write(*,*) 'R(1:5) ', R(1:5)
             write(*,*) 'Xr(1:5) ', Xr(1:5)
             write(*,*) 'X(1:5) ', X(1:5)
           endif

           alpha =   one
           beta = one
           call pdgeadd( 'N', n,nrhs, alpha, R,1,1,descX,                  &
     &                    beta, X,1,1,descX )

*          ----------------
*          compute residual
*          ----------------
           call pdgecopy( n,nrhs, B,ib,jb,descB, R, 1,1, descX )
           alpha = -one
           beta = one
#ifdef USE_PROF
           call profstart('pzgesvr:pzgemm')
#endif
           call pdgemm('N','N',n,nrhs,n,                                    &
     &          alpha, A,ia,ja,descA, X,1,1,descX,                          &
     &          beta, R,1,1,descX )
#ifdef USE_PROF
           call profend('pzgesvr:pzgemm')
#endif
         rnorm_old = rnorm
        enddo

*       ---------------------
*       finally copy solution
*       ---------------------
        call pdgecopy( n,nrhs, X,1,1,descX, B,ib,jb,descB )
        info = 0
        if (.not.isdone) then
                info = -(n+1)
        endif

 999    continue

        deallocate( Ar, X, R, Xr, dwork,rwork,stat=ierr)
        if (ierr.ne.0) then
         write(*,*) 'deallocate(Ar) return ierr=',ierr
        endif

        return

        end subroutine pdgesvr
