      subroutine ga_create_list(nxtot,nytot,name,nx,ny,nz,handle)
*
* $Id: ga_list.F,v 1.5 2004-01-31 22:44:56 d3j191 Exp $
*
      implicit none
c
      integer nxtot,nytot,nx,ny,nz
      integer handle(*)
      character*2 name
      character*4 nami
c
      integer i
c
#include "global.fh"
#include "mafdecls.fh"
c
      do 1 i=1,nz
      write(nami,'(a,i2.2)') name,i
      if(.not.ga_create(mt_dbl,nxtot,nytot,nami,nx,ny,handle(i))) then
      call md_abort('Failed to create global array list',i)
      endif
    1 continue
c
      return
      end
      subroutine ga_destroy_list(nz,handle)
c
      implicit none
c
      integer nz
      integer handle(*)
c
      integer i
c
#include "global.fh"
c
      do 1 i=1,nz
      if(.not.ga_destroy(handle(i)))
     + call md_abort('Failed to destroy global array list',i)
    1 continue
c
      return
      end
      subroutine ga_zero_list(nz,handle)
c
      implicit none
c
      integer nz
      integer handle(*)
c
      integer i
c
#include "global.fh"
c
      do 1 i=1,nz
      call ga_zero(handle(i))
    1 continue
c
      return
      end
