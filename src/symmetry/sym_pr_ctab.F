      subroutine sym_print_char_table(geom)
C$Id$
      implicit none
#include "errquit.fh"
#include "sym.fh"
      integer geom              ! [input]
c
c     Print the character table for the group to Fortran unit 6
c
      integer maxireps
      parameter (maxireps=20)
      integer nir, nop, iclass_dim(maxireps)
      character*8 zir(maxireps), zclass(maxireps)
      double precision chars(maxireps*maxireps)
      character*8 group_name
      integer i, j
c
      call sym_group_name(geom,group_name)
      if (.not. sym_char_table(group_name, nop, nir, iclass_dim,
     $     zir, zclass, chars)) 
     $     call errquit('sym_print_char_table: no char table',0,
     &       GEOM_ERR)
      write(6,111) group_name
      write(6,112) (iclass_dim(i),zclass(i),i=1,nir)
      do j = 1, nir
         write(6,113) zir(j), (chars((j-1)*nir+i),i=1,nir)
      enddo
      write(6,*)
 111  format(//8x,10x,'Character table for group ', a/)
 112  format(8x,8x,1x,20(i2,a6))
 113  format(8x,a8,20f8.3)

c
      end
