C NODEBUG 
C NODEBUG 
C
Cdcputime00() is part of BLACS.
C
      double precision function dcputime00()
      double precision MPI_Wtime
      external MPI_Wtime
      dcputime00 = MPI_Wtime()
      return
      end
      subroutine profinit()
      double precision dictstart(256 ), dicttotal(256 )
      integer dictcount(256 )
      integer nroutine, nlevels
      character*80 dictname(256 )
      character*80 lastroutine(256 )
      common /profilecom/ dictstart,dicttotal,dictcount,nroutine,       
     &nlevels, dictname, lastroutine
      save /profilecom/
      double precision dcputime00 
      external dcputime00 
      integer i 
      nroutine = 0
      do 23000 i = 1 , 256 
      dictname(i) = ' '
      dictstart(i) = 0.0
      dictcount(i) = 0
      dicttotal(i) = 0.0
23000 continue
23001 continue
      nlevels = 0
      do 23002 i = 1 , 256 
      lastroutine(i) = ' '
23002 continue
23003 continue
      return
      end
      subroutine profstart(rname)
      character*(*) rname
      double precision dictstart(256 ), dicttotal(256 )
      integer dictcount(256 )
      integer nroutine, nlevels
      character*80 dictname(256 )
      character*80 lastroutine(256 )
      common /profilecom/ dictstart,dicttotal,dictcount,nroutine,       
     &nlevels, dictname, lastroutine
      save /profilecom/
      double precision dcputime00 
      external dcputime00 
      integer j, i, ipos
      logical found
      character*80 name
C ======= start execution =  
      name = rname
      nlevels = nlevels + 1
      lastroutine(nlevels) = name
      found = .false. 
      do 23004 j = 1 , nroutine 
      i = nroutine - j + 1
C count down loop  
C heuristic filter for faster execution  
      if(dictname(i) (1: 1) .ne. name(1:1))then
      goto 23004
      endif
      found = ( dictname(i) .eq. name ) 
      if(found)then
      ipos = i
      goto 23005
      endif
23004 continue
23005 continue
      if(.not.found)then
      nroutine = nroutine + 1
      ipos = nroutine
      dictname(ipos) = name 
      dictcount(ipos) = 0
      dicttotal(ipos) = 0.0
      endif
      dictstart(ipos) = dcputime00 ()
      dictcount(ipos) = dictcount(ipos) + 1
      return
      end
      subroutine profend(rname)
      character*(*) rname
      double precision dictstart(256 ), dicttotal(256 )
      integer dictcount(256 )
      integer nroutine, nlevels
      character*80 dictname(256 )
      character*80 lastroutine(256 )
      common /profilecom/ dictstart,dicttotal,dictcount,nroutine,       
     &nlevels, dictname, lastroutine
      save /profilecom/
      double precision dcputime00 
      external dcputime00 
      integer j, i, ipos
      logical found
      character*80 name
      double precision tend
C ======= start execution =  
      tend = dcputime00 ()
      name = rname
      if(.not. ( name .eq. lastroutine(nlevels) ) )then
      print *, ' ** profend: name != lastroutine ', name
      stop ' ** ERROR ** '
      endif
      found = .false. 
      do 23014 j = 1 , nroutine 
      i = nroutine - j + 1
C count down loop  
C heuristic filter for faster execution  
      if(dictname(i) (1: 1) .ne. name(1:1))then
      goto 23014
      endif
      found = ( dictname(i) .eq. name ) 
      if(found)then
      ipos = i
      goto 23015
      endif
23014 continue
23015 continue
      if(.not.found)then
      print *, ' ** profend: routine name not found', name
      stop ' ** ERROR ** '
      endif
      dicttotal(ipos) = dicttotal(ipos) + (tend - dictstart(ipos))
      nlevels = nlevels - 1
      return
      end
      subroutine profstat()
      double precision dictstart(256 ), dicttotal(256 )
      integer dictcount(256 )
      integer nroutine, nlevels
      character*80 dictname(256 )
      character*80 lastroutine(256 )
      common /profilecom/ dictstart,dicttotal,dictcount,nroutine,       
     &nlevels, dictname, lastroutine
      save /profilecom/
      double precision dcputime00 
      external dcputime00 
      integer i 
      character*80 fname
      fname = 'profstat.dat'
      open(17 , file = fname, access = 'sequential', form = 'formatted' 
     &)
      rewind(17 )
      do 23022 i = 1 , nroutine 
      write(17 , 999) dictname(i), dictcount(i), dicttotal(i)
      write(6, 999) dictname(i), dictcount(i), dicttotal(i)
999   format(a20, ' was called ', i10, ' times, total ', f10 .2,        
     &' secs')
23022 continue
23023 continue
      close(17 )
      return
      end
