      subroutine sym_irrepname(geom,irrep,name)
c
c $Id$
c
c Returns the name of irrep.  If irrep is zero, returns
c '?' as the name.
c
      implicit none
#include "geom.fh"
#include "sym.fh"
      integer geom       ! [input] Geometry handle
      integer irrep      ! [input] Irrep number
      character*(*) name ! [output] name of the irrep
      character*8 zname  ! name of symmetry group
      integer maxop
      integer maxireps
      parameter (maxop=120)
      parameter (maxireps=20)
      integer nop
      integer nir
      integer class_dim(maxireps)
      character*8 zir(maxireps)
      character*8 zclass(maxireps)
      double precision chars(maxireps*maxireps)
c
      call sym_group_name(geom,zname)
      if (.not.sym_char_table(zname,nop,nir,class_dim,
     1  zir,zclass,chars)) call errquit
     2  ('sym_irrepname: failed to get character table',0,0)
      if (irrep.eq.0) then
        name='?'
      else
        name=zir(irrep)
      endif
      return
      end
