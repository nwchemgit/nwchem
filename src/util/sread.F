      subroutine sread(itape,a,n)
      implicit none
#include "errquit.fh"
*
* $Id$
*
      integer n,itape
      double precision a(n)
c
      integer left,nbuf,ibuf,lenbuf,m
      parameter (lenbuf = 524287) ! 4*1024*1024/8 - 1
c
      if (n.le.0) return
      left = n
      nbuf = (n-1)/lenbuf + 1
      do 10 ibuf = 1,nbuf
        m = min(lenbuf, left)
        call ssread(itape, a(1 + (ibuf-1)*lenbuf), m)
        left = left - m
10    continue
      if (left.ne.0) call errquit('sread: left .ne. 0',left, INPUT_ERR)
c
      end
      subroutine ssread(itape,a,n)
      implicit none
#include "errquit.fh"
      integer n,itape
      double precision a(n)
c
      read(itape,end=10,err=20) a
      return
 10   call errquit('ssread: end of file reading from unit ', itape,
     &       DISK_ERR)
 20   call errquit('ssread: error reading from unit ', itape,
     &       DISK_ERR)
      return
      end
c
      subroutine sread_ft(itape,a,n,ierr)
      implicit none
#include "errquit.fh"
      integer n,itape,ierr
      double precision a(n)
c
      integer left,nbuf,ibuf,lenbuf,m
      parameter (lenbuf = 524287) ! 4*1024*1024/8 - 1
c
      if (n.le.0) return
      left = n
      nbuf = (n-1)/lenbuf + 1
      do 10 ibuf = 1,nbuf
        m = min(lenbuf, left)
        call ssread_ft(itape, a(1 + (ibuf-1)*lenbuf), m,ierr)
        if(ierr.ne.0) return
        left = left - m
10    continue
      if (left.ne.0) then
         ierr=2046
         write(0,*)'sread: left .ne. 0',left
      endif
      return
c
      end
      subroutine ssread_ft(itape,a,n,ierr)
      implicit none
#include "errquit.fh"
      integer n,itape,ierr
      double precision a(n)
c
      read(itape,end=10,err=20,iostat=ierr) a
      return
 10   ierr=-20
      write(0,100) 'ssread: end of file reading from unit ', itape,
     .     ' IOSTAT= ',ierr
      return
 20   ierr=20
      write(0,100) 'ssread: error reading from unit ', itape,
     .     ' IOSTAT= ',ierr
      return
 100  format(A,I4,A,I4)
c
      end
