      subroutine ana_esppt(isel,xs,xw,qs,qw)
c
      implicit none
c
#include "ana_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "rtdb.fh"
c
      integer isel(msa)
      real*8 xs(msa,3),xw(mwm,mwa,3),qs(msa),qw(mwa)
      integer i,j
      real*8 espval,point(3)
c
      open(unit=lfnpin,file=filpin(1:index(filpin,' ')-1),
     + form='formatted',status='old')
      open(unit=lfnpin,file=filpou(1:index(filpou,' ')-1),
     + form='formatted',status='unknown')
      rewind(lfnpin)
      rewind(lfnpou)
      espval=0.0d0
      do 1 j=1,nsa
      if(isel(j).gt.0) then
      espval=espval+138.9354*qs(i)/sqrt((point(1)-xs(i,1))**2+
     + (point(2)-xs(i,2))**2+(point(3)-xs(i,3))**2)
      endif
    1 continue
      write(*,'(3f12.6,e12.5)') point,espval
      close(unit=lfnpin,status='keep')
      close(unit=lfnpou,status='keep')
c
      return
      end
