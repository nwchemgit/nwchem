
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
      SUBROUTINE nwfft3d(nd1,nd2,nd3,isign,
     .                   x1,x2,map1,map2,doMap,jstatus)
      implicit none
      INTEGER nd1,nd2,nd3, isign, jstatus
      INTEGER map1(*), map2(*)
      INTEGER doMap
      DOUBLE PRECISION x1(*),x2(*)
#include "mafdecls.fh"
c *** input parameters described
c
c  nd1,nd2,nd3: dimension of x1 and x2
c  isign: 1 => forward fft;  -1 => reverse fft
c  x1: data to be fft'ed. input(x1) = output(x1), i.e., non-destructive.
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
c
c *** local variables
      integer d1,d2,d3, owner, n123
      double precision scale
      integer blksz, computedMap
      data computedMap /0/
      save computedMap
#ifdef ESSLFFT
      integer icontext, ip(40) 
      integer initESSL, iESSL
      data iESSL /0/
      save iESSL
#endif
#ifdef EMSLFFT
      integer  nwords, nscratch
      logical  bstatus
      integer  dblh,inth, idbl,iint
      integer  initSCRATCH, n123ma
      data initSCRATCH,n123ma /0,0/
      save initSCRATCH, n123ma, dblh,inth,idbl,iint
#endif
        integer isize
c       parameter (isize=128*128*128)
c       integer iwork(isize)
c	double precision work(isize)
       integer nodeid,nnodes,me,nproc

      me = nodeid()
      nproc = nnodes()

c     if (me.EQ.0) print*,' ENTER FFT ',nd1,nd2,nd3,' isign=',isign

c --- calculate constants ----
      n123 = nd1*nd2*nd3
c scale for forward (isign==1) and reverse (isign==-1) fft 
      if (isign.EQ.-1) then
         scale = 1.D0/REAL(n123)
      else
         scale = 1.D0
      endif
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
           call errquit('nwfft3d: ESSL initialization failed',me)
         endif
      endif
#endif

#ifdef EMSLFFT
c allocate scratch space: number of (complex) elements x 2 
      if (n123.GT.n123ma.OR.initSCRATCH.EQ.0) then
        initSCRATCH = 1
c free previous space if present 
        if (n123ma.GT.0) then
           bstatus = MA_free_heap(dblh) 
           bstatus = MA_free_heap(inth) 
        endif
        nscratch = n123*2*2 
        nwords = MA_inquire_heap(MT_DBL)
        if (nwords.LT.nscratch) then
          call errquit('nwfft3d: insufficient heap for dscratch',me)
        else 
           if(.NOT.MA_alloc_get(MT_DBL,nscratch,'fftdbl',dblh,idbl))then
            call errquit('nwfft3d: insufficient heap for dscratch',me)
           endif
        endif

        nwords = MA_inquire_heap(MT_INT)
        if (nwords.LT.nscratch) then
           call errquit('nwfft3d: insufficient heap for iscratch',me)
        else 
           if(.NOT.MA_alloc_get(MT_INT,nscratch,'fftint',inth,iint))then
            call errquit('nwfft3d: insufficient heap for iscratch',me)
           endif
        endif
        n123ma = n123
      endif
#endif


c determine map if necessary (required once) for EMSL routine
#ifdef EMSLFFT
      if (blksz.LE.0) then
        call errquit('nwfft3d:: internal error: bad blocksize',blksz)
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
#endif


#ifdef ESSLFFT
      call pdcft3(x1,x2,nd1,nd2,nd3,isign,scale,icontext,ip)
#endif
#ifdef EMSLFFT
c     * make the EMSL destructive interface, nondestructive
      do d1 = 1, nd1*nd2*blksz*2
         x2(d1) = x1(d1)
      enddo
c	print*,'-->',isign,nd1,nd2,nd3,x2(1),map1(1),map2(1),blksz
c     call pfft3d(isign,nd1,nd2,nd3,x2,map1,map2,
c    .            isize,work,isize,iwork,jstatus)
      call pfft3d(isign,nd1,nd2,nd3,x2,map1,map2,
     .            nscratch,dbl_mb(idbl),nscratch,int_mb(iint),jstatus)
      if (isign.EQ.-1) then
        do d1 = 1, nd1*nd2*blksz*2
          x2(d1) = x2(d1)*scale
        enddo
      endif
#endif

c status check
      if (jstatus.NE.0) then
#ifdef EMSLFFT
        call errquit('nwfft3d: pfft3d() error return',me)
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
