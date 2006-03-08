      subroutine esp_rms(x,qf,g,v)
c
c $Id: esp_rms.F,v 1.4 2006-03-08 05:13:00 d3j191 Exp $
c
      implicit none
c
#include "esp_params.fh"
#include "esp_common.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      real*8 x(3,natoms),qf(ndim,6)
      real*8 g(3,mxgrid),v(mxgrid)
c
      integer i,j,k
      real*8 eq
c
      do 1 k=1,nf
      rms(1,k)=zero
      rms(2,k)=zero
      do 2 i=ngfr,ngto
      eq=zero
      do 3 j=1,natoms
      eq=eq+qf(j,k)/dsqrt((g(1,i)-x(1,j))*(g(1,i)-x(1,j))+
     + (g(2,i)-x(2,j))*(g(2,i)-x(2,j))+(g(3,i)-x(3,j))*(g(3,i)-x(3,j)))
    3 continue
      rms(1,k)=rms(1,k)+(v(i)-eq)*(v(i)-eq)
      rms(2,k)=rms(2,k)+dabs(hundrd*(v(i)-eq)/v(i))
    2 continue
    1 continue
      if(np.gt.1) then
      call ga_dgop(mre_004,rms,2*nf,'+')
      endif
      do 4 k=1,nf
      rms(1,k)=cau2kj*dsqrt(rms(1,k))/dble(ngrid)
      rms(2,k)=rms(2,k)/dble(ngrid)
    4 continue
c
      return
      end
