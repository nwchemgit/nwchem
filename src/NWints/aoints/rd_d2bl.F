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
      integer*4 nipv,iretbv
      parameter (nipv=4,iretbv=0)
c
      integer*4   ninfmx
      parameter(ninfmx=10)
c
c     # bummer error types.
      integer*4   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
c  ##  integer*4 section
c
      integer*4 info(ninfmx),itypea,itypeb,ifmt,ibvtyp,ierr,ibitv(1)
      integer*4 last,labs(4,*)
      integer*4 nbuf
      integer*4 sfile
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
!      write(*,*)'gk: before sifrd2'
!      write(*,*)'gk:  sfile = ', sfile
!      write(*,*)'gk:  info(1):', info(1)
!      write(*,*)'gk:  info(2):', info(2)
!      write(*,*)'gk:  info(3):', info(3)
!      write(*,*)'gk:  info(4):', info(4)
!      write(*,*)'gk:  info(5):', info(5)
!      write(*,*)'gk:  last = ', last
cgk end
      call sifrd2(sfile,   info,   nipv,    iretbv,
     &            buffer,  nbuf,   last,    itypea,
     &            itypeb,  ifmt,   ibvtyp,  vals,
     &            labs,    ibitv,  ierr)
cgk debug
!      write(*,*)'gk: after sifrd2'
!      write(*,*)'in rd_d2bl, itypea = ', itypea
!      write(*,*)'in rd_d2bl, itypeb = ', itypeb
!      write(*,*)'in rd_d2bl, last = ', last
!      write(*,*)'in rd_d2bl, labs(1,1) = ', labs(1,1)
!      write(*,*)'in rd_d2bl, labs(2,1) = ', labs(2,1)
!      write(*,*)'in rd_d2bl, labs(3,1) = ', labs(3,1)
!      write(*,*)'in rd_d2bl, labs(4,1) = ', labs(4,1)
!      write(*,*)'            vals(1)   = ', vals(1)
cgk end
c
      if (ierr.gt.0) call bummer("error in sifrd2",0,faterr)
      return
      end
