*
* $Id$
*
c************************************************************************
c Complex 3D FFT interface.
c Based on a blockwise distribution of slabs dimensioned id1xid2.
c The blockwise (and contiguous) distribution *must* be performed
c outside of this routine.
c 
c GENERAL compilation:
c       uses the portable "EMSL" FFT,  pfft1.0. 
c IBM SP2:
c       uses the faster ESSL FFT.  
c 
c Non-destructive interface (see below).
c ------------------------------------------------------------------------
      SUBROUTINE nwfft3d(nd1,nd2,nd3,idir,
     .                   x1,x2,map1,map2,doMap,
     .                   iwrklen,iwrk,rwrklen,rwrk,jstatus,ffttype)
      implicit none
#include "errquit.fh"
      INTEGER nd1,nd2,nd3, isign,jstatus,idir
      INTEGER map1(*),map2(*)
      INTEGER doMap
      DOUBLE PRECISION x1(*),x2(*)
      INTEGER rwrklen
      DOUBLE PRECISION rwrk(*)
      INTEGER iwrklen
      INTEGER iwrk(*)
      INTEGER ffttype
#include "mafdecls.fh"
c *** input parameters described
c
c  nd1,nd2,nd3: dimension of x1 and x2
c  isign: 1 => forward fft;  -1 => reverse fft; -2 => reverse, no scaling
c  x1: data to be fft-ed. input(x1) = output(x1), i.e., non-destructive.
c  x2: output fft(x1). Note that EMSL uses a destructive interface, while
c      TYPE: double precsion  => x1(dim), dim  = nd1*nd2*nd3*2 
c            double complex   => x1(dim), dim  = nd1*nd2*nd3 
c  map1: (EMSL) id of processor owning ith plane: Computed here according
c        to a blockwise decomposition.
c  map2: (EMSL) id of processor owning ith plane after transpose: computed here
c        according to a blockwise decomposition.
c        NOTE1: map1 and map2 dimensions are greater than nd3-1
c        NOTE2: map1 and map2 are NOT used for ESSL call (may dimension to 1)
c  doMap: (EMSL) forces recomputation of map; map computed on first
c              call and saved.  Can also be use 
c  rwrklen,iwrklen: if <> 0, use for EMSL routine scratch memory 
c  ffttype:  1=>EMSL fft;  2=>ESSL fft
c
c *** local variables
      integer d1,d3,owner
      double precision scale
      integer blksz,computedMap
      integer n123
#ifdef ESSLFFT
      integer icontext,ip(40) 
      integer initESSL,iESSL
#endif
#ifdef EMSLFFT
      integer  nwords, nrealscr,ma_nrealscr, ma_nintscr
      logical  bstatus
      integer  dblh,inth, idbl,iint 
      integer  initSCRATCH, mamem 
      integer nvecs1, nvecs2, n1max, n2max, nele, nelebk
      integer lstcnt, ii, m2, n3
#endif
c       integer isize
c       parameter (isize=128*128*128)
c       integer iwork(isize)
c	double precision work(isize)
       integer nodeid,nnodes,me,nproc
#ifdef ESSLFFT
      data iESSL /0/
      save iESSL
#endif
#ifdef EMSLFFT
      data initSCRATCH,ma_nrealscr,ma_nintscr /0,0,0/
      save initSCRATCH, ma_nrealscr,ma_nintscr, dblh,inth,idbl,iint
#endif
*
      data computedMap /0/
      save computedMap
*
      me = nodeid()
      nproc = nnodes()
      isign=idir
      if(isign.lt.-1) isign=-1

c      if (me.EQ.0) print*,' ENTER FFT ',nd1,nd2,nd3,' isign=',isign

c --- calculate constants ----
      n123 = nd1*nd2*nd3
c scale for forward (isign==1) and reverse (isign==-1) fft 
      if (idir.EQ.-1) then
         scale = 1.D0/dble(n123)
      else
         scale = 1.D0
      endif

#ifdef ESSLFFT
      ffttype = 2
#else
      ffttype = 1
#endif
     
#ifdef EMSLFFT
      blksz = nd3/ nproc
#endif

#ifdef ESSLFFT
c     ** do not use defaults
      ip(1) = 1
c     ** return normal form, i.e., x2 = Transpose(FFT(x1)), 
c       where FFT transposes x1.
      ip(2) = 1
c     ** fft x1(nd1,nd2,*) 
      ip(20) = nd1
      ip(21) = nd2

      if (iESSL.EQ.0) then
         iESSL = 1
         jstatus = initESSL(icontext)
         if (jstatus.LT.0) then
           call blacs_abort(icontext,1)
           call errquit('nwfft3d: ESSL initialization failed',0,
     &       UNKNOWN_ERR)
         endif
      endif
#endif



c determine map if necessary (required once) for EMSL routine
#ifdef EMSLFFT
      if (blksz.LE.0) then
       call errquit('nwfft3d:: internal error: bad blocksize',0,
     &       UNKNOWN_ERR)
      endif
      if (doMap.EQ.1 .OR. computedMap.EQ.0) then
        owner = 0
        do d3 = 1, nd3
          if (d3.GT.blksz*(owner+1)) owner = owner+1
          map1(d3) = owner
          map2(d3) = owner
        enddo
        computedMap = 1
      endif
c     
      
      NELEBK = 2 * ND1           ! complex to complex
      NELEBK = 2 * ( ND1 / 2 + 1 ) ! real-> complex
c
      if ( me.eq.0) then
         do ii=1, nd3
            write(*,*) ' map1, map2 ', ii, map1(ii), map2(ii)
         enddo
         call ffflush()
      endif
c
      NELEBK = max(2*nd1, 2*(nd1/2+1))
      NVECS1 = LSTCNT( ME, MAP1, nd3 )
      NVECS2 = LSTCNT( ME, MAP2, nd3 )
      N1MAX = 0
      N2MAX = 0
      do iI = 0, NPROC-1
         N1MAX = MAX( N1MAX, LSTCNT( iI, MAP1, N3 ) )
         N2MAX = MAX( N2MAX, LSTCNT( iI, MAP2, M2 ) )
      enddo
      NELE = NELEBK * MAX( M2 * N1MAX, N3 * N2MAX )
      nrealscr = nele
c
      if (rwrklen.NE.0.AND.iwrklen.NE.0) then
        mamem = 0 
      else
        mamem = 1 
      endif
      if(mamem.EQ.1.AND.
     .         (nrealscr.GT.ma_nrealscr.OR.initSCRATCH.EQ.0))then
c free previous space if present 
        if (ma_nrealscr.GT.0) then
           bstatus = MA_free_heap(dblh) 
        endif
        nwords = MA_inquire_heap(MT_DBL)
        if (nwords.LT.nrealscr) then
          call errquit('nwfft3d: insufficient heap for dscratch',0,
     &       MA_ERR)
        else 
           if(.NOT.MA_alloc_get(MT_DBL,nrealscr,'fftdbl',dblh,idbl))then
            call errquit('nwfft3d: insufficient heap for dscratch',0,
     &       MA_ERR)
           endif
        endif
        ma_nrealscr = nrealscr
      endif

c integer scratch depends on the number of processes. Assume, therefore,
c it will not change during the course of a run. (see pfft1.0/pfft/chk3d.F)
      if (mamem.EQ.1.AND.initSCRATCH.EQ.0) then
        ma_nintscr = 5*nproc+100
        nwords = MA_inquire_heap(MT_INT)
        if (nwords.LT.ma_nintscr) then
           call errquit('nwfft3d: insufficient heap for iscratch',0,
     &       MA_ERR)
        else 
          if (.NOT.MA_alloc_get(MT_INT,ma_nintscr,
     .                         'fftint',inth,iint)) then
           call errquit('nwfft3d: insufficient heap for iscratch',0,
     &       MA_ERR)
          endif
        endif
        initSCRATCH = 1
      endif
#endif



#ifdef ESSLFFT
      call pdcft3(x1,x2,nd1,nd2,nd3,isign,scale,icontext,ip)
#endif
#ifdef EMSLFFT
c     * make the EMSL destructive interface, nondestructive
c
      do d1 = 1, nd1*nd2*blksz*2
         x2(d1) = x1(d1)
      enddo
c	print*,'-->',isign,nd1,nd2,nd3,x2(1),map1(1),map2(1),blksz
c     call pfft3d(isign,nd1,nd2,nd3,x2,map1,map2,
c    .            isize,work,isize,iwork,jstatus)

      
      if (mamem.EQ.1) then
        call pfft3d(isign,nd1,nd2,nd3,x2,map1,map2,
     .       ma_nrealscr,dbl_mb(idbl),ma_nintscr,int_mb(iint),jstatus)
      else
        call pfft3d(isign,nd1,nd2,nd3,x2,map1,map2,
     .       rwrklen,rwrk,iwrklen,iwrk,jstatus)
      endif
      if (isign.EQ.-1) then
        do d1 = 1, nd1*nd2*blksz*2
          x2(d1) = x2(d1)*scale
        enddo
      endif
#endif

c status check
      if (jstatus.NE.0) then
#ifdef EMSLFFT
        call errquit('nwfft3d: pfft3d() error return',0, UNKNOWN_ERR)
#endif 
#ifdef FFTDEBUG
        print*,me,'WARNING: forward fft status = ',jstatus,nd1,nd2,nd3
#endif 
      endif

      return
      end

c************************************************************************

      double precision FUNCTION inorm(x,y,n)
      implicit none
      integer n
      double precision x(*),y(*)
      double precision maxval, newval
      integer i, me, nodeid

      me = nodeid()
      maxval = 0.D0
      do i = 1, n
        newval = ABS(x(i)-y(i)) 
        if (newval.GT.maxval) then
c        write(6,*) me,' i,x1,x2: ',i,x(i),y(i)
          maxval=newval
        endif
      enddo 
      inorm = maxval
      return
      end

c************************************************************************
