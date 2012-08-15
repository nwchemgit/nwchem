c===============================================================================================
ccase...Coulomb attenuation utilities
c===============================================================================================
c     need to initialize doscreen since case_setflags might not be called
      block data case_data
#include "case.fh"
      data doscreen /.false./
      data cam_srhf /.false./
      end
c
      subroutine case_setflags(screenflag)

c     Turn on/off the attenuation based on the screenflag

      implicit none

#include "stdio.fh"
#include "mafdecls.fh"
#include "case.fh"

      logical screenflag

      doscreen = .false.
      doscreen = screenflag
      call ga_sync()   ! synchronise over all processors

      return
      end
c
c===============================================================================================
ccase...Calculate the Coulomb attenuation (CASE) factors
c===============================================================================================
c
c================================================================================================
ccase...For the Texas 2e4c integrals
c================================================================================================
c
      subroutine case_txs(rho, rys, const)

      implicit none

#include "stdio.fh"
#include "mafdecls.fh"
#include "case.fh"

      double precision rho, rys, const
      double precision cam_omega2,mu

c     Calculate the screening factors
      mu = 1.d0
      cam_omega2 = cam_omega*cam_omega
      if (doscreen) mu = cam_omega2/(cam_omega2 + rho)

c     Debug print
c     write(luout,*) "In case_txs: cam_omega:",cam_omega
c     write(luout,*) "In case_txs: cam_exch:",cam_exch
c     write(luout,*) "In case_txs: doscreen:",doscreen

c     Multiply in the screening factors
      rys = mu*rys
      const = const*dsqrt(mu)
c 
      return
      end
c
c=================================================================================================
ccase...For the McMurchie-Davidson 2e4c integrals
c=================================================================================================
c
      subroutine case_md(rho)

      implicit none

#include "stdio.fh"
#include "mafdecls.fh"
#include "case.fh"

      double precision rho
      double precision cam_omega2,mu

c     Calculate the screening factors
      mu = 1.d0
      cam_omega2 = cam_omega*cam_omega
      if (doscreen) mu = cam_omega2/(cam_omega2 + rho)
c
c     Debug print
c     write(luout,*) "In case_md: cam_omega:",cam_omega
c     write(luout,*) "In case_md: cam_exch:",cam_exch
c     write(luout,*) "In case_md: doscreen:",doscreen

c     Multiply in the screening factors
      rho = mu*rho
c
      return
      end
c $Id$
