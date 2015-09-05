      subroutine dimqm_energy(nAtoms, muind, eqme, eqmn)
      implicit none
#include "stdio.fh"
#include "util.fh"
#include "dimqm_constants.fh"
#include "global.fh"
      integer nAtoms
      double precision muind(3,nAtoms)
      double precision eqme(3, nAtoms)
      double precision eqmn(3, nAtoms)
c    Local Variables
      double precision en_dipole_el
      double precision en_dipole_nuc
      double precision energy
      character*50 d
      character*50 dd
      
      if(ga_nodeid().ne.0) return
c
      dd =
     $ '=================================================='
      d =
     $ '--------------------------------------------------'
c
c    Calculate the DIM energy
      en_dipole_el = -HALF * SUM(muind * eqme)
      en_dipole_nuc= -HALF * SUM(muind * eqmn)
      energy = en_dipole_el + en_dipole_nuc
c
c
      write(luout,*) dd
      write(LuOut,*) 'DIM/QM Energy'
      write(LuOut,*) d
      write(LuOut,111) en_dipole_el
      write(LuOut,222) en_dipole_nuc
      write(LuOut,333) energy
      write(luout,*) dd
      call util_flush(LuOut)
 111  format(' Polarization (Dipole, el.) = ',f20.10)
 222  format(' Polarization (Dipole, nuc) = ',f20.10)
 333  format('              Total Energy  = ',f20.10)

      end subroutine dimqm_energy
