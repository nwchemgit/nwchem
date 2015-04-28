      subroutine util_sillysort(value,
     R     kaux,klist,n,found)
      implicit none
#include "mafdecls.fh"
C
C     $Id: ccsd_t_neword.F 26884 2015-02-25 02:06:49Z edo $
c very simple sorting
c
      integer value    ! [in]
      integer n        ! [in]
      integer found    ! [in/out]
      integer kaux(*)  ! [in] values to sort
      integer klist(*) ! [out] index in orginal sorting
c
      integer i,j,m
c
      do m=1,n
         if(kaux(m).gt.value) then
            found=found+1
            klist(found)=m
c     invalidate
            kaux(m)=-99
         endif
      enddo
      return
      end
