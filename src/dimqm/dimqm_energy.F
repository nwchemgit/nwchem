      subroutine dimqm_energy(nAtoms, muind, eqme, eqmn)
      implicit none
#include "stdio.fh"
#include "util.fh"
#include "dimqm_constants.fh"
#include "global.fh"
c  jbecca: START
#include "dimqm.fh"
c  jbecca: END
      integer nAtoms
      double precision muind(3,nAtoms)
      double precision eqme(3, nAtoms)
      double precision eqmn(3, nAtoms)
c    Local Variables
      double precision en_dipole_el
      double precision en_dipole_nuc
      double precision efld_int
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
      efld_int  =  ZERO
      if (lefield) then
         efld_int  =  SUM(muind(1,:) * finiteFieldX)
     $             +  SUM(muind(2,:) * finiteFieldY)
     $             +  SUM(muind(3,:) * finiteFieldZ)
         efld_int  =  -HALF * efld_int
      endif
      energy = en_dipole_el + en_dipole_nuc
c
c
      write(luout,*) dd
      write(LuOut,*) 'DIM/QM Energy'
      write(LuOut,*) d
      write(LuOut,111) en_dipole_el
      write(LuOut,222) en_dipole_nuc
c      write(LuOut,444) efld_int
      write(LuOut,333) energy
      write(luout,*) dd
      write(luout,*) ''
      if (lefield) then 
         write(luout,*) dd
         write(luout,*)'DIM System Energy'
         write(LuOut,*) d
         write(luout,555) efld_int
         write(luout,*) dd
      endif

      call util_flush(LuOut)
 111  format(' Polarization (Dipole, el.) = ',f20.10)
 222  format(' Polarization (Dipole, nuc) = ',f20.10)
 333  format('              Total Energy  = ',f20.10)
 444  format(' Interaction with Efield    = ',f20.10)
 555  format('                      Total = ',f20.10)

      end subroutine dimqm_energy
