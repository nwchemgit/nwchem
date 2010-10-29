      subroutine tce_nud1(dipole)
c
c $Id$
c
c     Computes nuclear contributions of dipole moments
c
      implicit none
#include "rtdb.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "sym.fh"
#include "sf.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      double precision dipole(3)       ! Total dipole moments
      integer natoms                   ! Number of atoms
      double precision center(3,10000) ! Nuclear positions
      double precision charge(10000)   ! Nuclear charges
      character*16     tag(10000)      ! Nuclear tags
      integer iatom
      logical  geom_ncent
      external geom_ncent
      logical  geom_cart_get
      external geom_cart_get
      logical nodezero
c
      nodezero = (ga_nodeid().eq.0)
      if (.not.geom_ncent(geom,natoms))
     1  call errquit('tce_nud1: geom_ncent failed',0,GEOM_ERR)
      if (natoms.gt.10000)
     1  call errquit('tce_nud1: too many atoms',natoms,GEOM_ERR)
      if (.not.geom_cart_get(geom,natoms,tag,center,charge))
     1  call errquit('tce_nud1: geom_cart_get failed',0,GEOM_ERR)
      dipole(1) = 0.0d0
      dipole(2) = 0.0d0
      dipole(3) = 0.0d0
      do iatom = 1,natoms
        if (nodezero.and.util_print('nud1',print_debug)) then
          write(LuOut,*) iatom,center(1,iatom),center(2,iatom),
     1                         center(3,iatom),charge(iatom)
        endif
        dipole(1) = dipole(1) + center(1,iatom)*charge(iatom)
        dipole(2) = dipole(2) + center(2,iatom)*charge(iatom)
        dipole(3) = dipole(3) + center(3,iatom)*charge(iatom)
      enddo
      return
      end
