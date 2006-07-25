c
c$Id:
c
      logical function atom_gfac(znuc,gfactor,isotope)
c
c     Function will extract nuclear g-factor for a given atom
c
      implicit none
#include "atomdata.fh"
c
      double precision znuc       ! [in]  nuclear charge
      double precision gfactor    ! [out] nuclear g-factor
c
      integer iatom, niso, i, iloc, isotope
      double precision abun, spin, gfac
c
      iatom = atomlocator(2,nint(znuc))
      niso = atomlocator(1,nint(znuc))
c
c     Determine isotope with spin and highest abundance
c
      abun = 0.0d0
      spin = 0.0d0
      gfac = 0.0d0
      isotope = 0
      do i = 1, niso
         iloc = iatom+i-1
         if (atomprop(3,iloc).ne.0.0d0) then
            if (atomprop(2,iloc).ge.abun) then
               isotope = nint(atomprop(1,iloc))
               abun = atomprop(2,iloc)
               spin = atomprop(3,iloc)
               gfac = atomprop(4,iloc)
            endif
         endif
      enddo
      gfactor = 0.0d0
      if (spin.gt.0.0d0) gfactor = gfac/spin
c
      atom_gfac = .true.
c
      return
      end      
