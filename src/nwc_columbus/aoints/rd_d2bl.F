c
c***********************************************************************
c
*deck indast
      subroutine rd_d2bl(
     & vals,   labs,   nbuf,   buffer,
     & last,   sfile,  info)
c
c
      implicit none
c
c  ##  paramter & common block section
c
      integer nipv,iretbv
      parameter (nipv=4,iretbv=0)
c
      integer   ninfmx
      parameter(ninfmx=10)
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
c  ##  integer section
c
      integer info(ninfmx),itypea,itypeb,ibvtyp,ierr,ibitv(1)
      integer last,labs(4,*)
      integer nbuf
      integer sfile
c
c  ##  real*8
c
      real*8 buffer(*)
      real*8 vals(*)
c
c
c-----------------------------------------------------
c
cgk debug
       !write(6,*)'gk: before sifrd2'
       !write(6,*)'gk:  sfile = ', sfile
       !write(6,*)'gk:  info(1):', info(1)
       !write(6,*)'gk:  info(2):', info(2)
       !write(6,*)'gk:  info(3):', info(3)
       !write(6,*)'gk:  info(4):', info(4)
       !write(6,*)'gk:  info(5):', info(5)
       !write(6,*)'gk:  last = ', last
cgk end
      call sifrd2(sfile,   info,   nipv,    iretbv,
     &            buffer,  nbuf,   last,    itypea,
     &            itypeb,  ibvtyp, vals,
     &            labs,    ibitv,  ierr)
cgk debug
       !write(6,*)'gk: after sifrd2'
       !write(6,*)'in rd_d2bl, itypea = ', itypea
       !write(6,*)'in rd_d2bl, itypeb = ', itypeb
       !write(6,*)'in rd_d2bl, last = ', last
       !write(6,*)'in rd_d2bl, labs(1,1) = ', labs(1,1)
       !write(6,*)'in rd_d2bl, labs(2,1) = ', labs(2,1)
       !write(6,*)'in rd_d2bl, labs(3,1) = ', labs(3,1)
       !write(6,*)'in rd_d2bl, labs(4,1) = ', labs(4,1)
       !write(6,*)'            vals(1)   = ', vals(1)
cgk end
c
      if (ierr.gt.0) call bummer("error in sifrd2",0,faterr)
      return
      end
