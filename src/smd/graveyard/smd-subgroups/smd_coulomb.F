      subroutine smd_coulomb_compute()
      implicit none

#include "inp.fh"
#include "mafdecls.fh"
#include "smd_system.fh"
#include "stdio.fh"
#include "errquit.fh"
c
      character*30 pname
      logical operiodic
c
      pname = "smd_coulomb"
      write(*,*) "in "//pname
c
      operiodic = smd_system_periodic()

c
c     NOTE THAT UPDATED VERLET LISTS
c     ARE REQUIRED
c     -----------------------------
      if(operiodic) then
        call smd_ewald_self()
        call smd_ewald_recip()
        call smd_ewald_excl()
        call smd_ewald_real()
      else
        call errquit("non periodic conditions not implemented",0,0)
      end if
      return
      END
c $Id$
