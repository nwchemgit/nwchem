      subroutine et_det(op,handle,n,detval)
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
#include "apiP.fh"
#include "bas.fh"
#include "geom.fh"
#include "cscf.fh"
c
c returns the determinant of the square matrix with file handle 'handle'.
c n is the dimension of the array, detval is the returned value of the
c determinant. If op='N', a=a. If op='T', a=a'.
c
      character*1 op
      integer handle
      integer n,i,j,ipvt(n),info
      double precision detval,a(n,n),work(n)
c
c
c
      do i = 1,n
      do j = 1,n
       if(op.eq.'N') call ga_get(handle, i, i, j, j, a(i,j), 1)  
       if(op.eq.'T') call ga_get(handle, j, j, i, i, a(i,j), 1)  
      enddo
      enddo
c
      call dgefa(a,n,n,ipvt,info)
      call dgedi (a,n,n,ipvt,detval,work,10)
c
      return
      end

