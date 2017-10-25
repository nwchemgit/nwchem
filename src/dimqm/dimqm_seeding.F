      subroutine dimqm_seed_init_tolerance(request)
      implicit none
#include "dimqm.fh"
#include "stdio.fh"
#include "global.fh"

      double precision request
c   Initialize finding the seed
      dimqm_seeded = .false.
c   Return if seeding is off
      if(dimqm_noseed) return
c   If the user set tolerance is less than the pre-seed tolerance, save the user's
c   setting and change to the pre-seed tolerance.  This will be reversed once seeding
c   has started.
      if(dimtol < request) then
        dimtol0 = dimtol
        dimtol = request
        if(ga_nodeid().eq.0 .and. ldebug) then
          write(luout,*) "User requested tolerance below", request
          write(luout,*) "Setting pre-seed tolerance to", request
        end if
      end if
      end subroutine dimqm_seed_init_tolerance

      subroutine dimqm_check_dipoles(tolerance, error)
      implicit none
#include "dimqm.fh"
#include "stdio.fh"
#include "global.fh"

      double precision tolerance, error
c   Return if seeding is off
      if(dimqm_noseed) return
c   Return if seed already selected
      if(dimqm_seeded) return
c   Return if this is the first cycle
      if(lfirst) return
c   Check to see if the error is below our seeding tolerance,
      if(error < tolerance) then
c       Return if we've already set the seed
c       Otherwise, set the seed and revert the tolerance back to user requested level
        dimqm_seeded = .true.
        if(ga_nodeid().eq.0 .and. ldebug) then
          write(luout,*) "Dipole error below ", tolerance
          write(luout,*) "Setting current system as seed"
          write(luout,*) "Reverting tolerance back to", dimtol0
        end if
        dimtol = dimtol0
      end if
      end subroutine dimqm_check_dipoles

      subroutine dimqm_check_scf(de, rms, derr, ipol, tolerance)
      implicit none
#include "dimqm.fh"
#include "stdio.fh"
#include "global.fh"

      double precision de, rms(2), tolerance, derr(2)
      integer ipol
      logical e, d1, d2, g1, g2
c   Return if dimqm is off
      if (.not.ldimqm) return
c   Return if seeding is off
      if(dimqm_noseed) return
c   Return if the seed hasn't been set yet
      if(.not.dimqm_seeded) return
      e = .false.
      d1 = .false.
      d2 = .false.
      g1 = .false.
      g2 = .false.
c   Check to see if the SCF RMS-density and delta E are below a given tolerance
c   If not, then this is not a good seed system and we reset the seed flag
      if(de .gt. tolerance) e = .true.
      if(dsqrt(rms(1)) .gt. tolerance) d1 = .true.
      if(derr(1) .gt. tolerance) g1 = .true.
      if(ipol .gt. 1) then
        if(dsqrt(rms(2)) .gt. tolerance) d2 = .true.
        if(derr(2) .gt. tolerance) g2 = .true.
      end if
      if(e .or. d1 .or. d2 .or. g1 .or. g2) then
        if(ga_nodeid().eq.0 .and. ldebug) then
          write(luout,*) "Seed was set, but some SCF error"
          write(luout,*) " is above ", tolerance
          write(luout,*) "Unsetting seed and reverting tolerance back"
          write(luout,*) "to 1.0d-4"
        end if
        dimqm_seeded = .false.
        dimtol = 1.0d-4
      end if
      end subroutine dimqm_check_scf
