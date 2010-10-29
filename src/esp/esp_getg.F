      subroutine esp_getg(g,v)
c
c $Id$
c
      implicit none
c
#include "esp_params.fh"
#include "esp_common.fh"
#include "mafdecls.fh"
c
      real*8 g(3,mxgrid),v(mxgrid)
      integer i,j
c
      call util_file_name('grid',.false.,.false.,grdfil)
      open(unit=lfngrd,file=grdfil,form='formatted',status='old')
      rewind(lfngrd)
      read(lfngrd,'(i10,f20.10)') ngrid,charge
      do 1 i=1,ngrid
      read(lfngrd,'(4f20.10)') (g(j,i),j=1,3),v(i)
    1 continue
      close(unit=lfngrd)
c
      return
      end
