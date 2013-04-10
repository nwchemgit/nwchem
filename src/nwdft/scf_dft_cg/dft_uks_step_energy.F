      subroutine dft_uks_step_energy(rtdb, step, g_search,
     $     eone, etwo, exc, enrep, energy, edisp, g_grad, nexc)
      implicit none
#include "errquit.fh"
#include "cuhf.fh"
#include "cscfps.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "cscf.fh"
c
      integer rtdb
      double precision step
      integer g_search
      double precision eone, etwo, enrep, energy, edisp
      double precision exc(2)
      integer nexc
      integer g_grad
c
c     Compute the energy and gradient at the given step
c
      integer g_tmp(2)
      if (oscfps) call pstat_on(ps_search)
c
c     Apply rotation
c
      if (.not. ga_create(MT_DBL, nbf, nmo,
     $                    'dft_uks_step_energy: alpha MOs', 32, 32,
     $                    g_tmp(1)))
     $          call errquit('dft_uks_step_energy: MOs', 0, GA_ERR)
      if (.not. ga_create(MT_DBL, nbf, nmo,
     $                    'dft_uks_step_energy: beta MOs', 32, 32,
     $                    g_tmp(2)))
     $          call errquit('dft_uks_step_energy: MOs', 0, GA_ERR)
c
      call uhf_rotate(step, g_search, g_movecs, g_tmp)
c
c     Compute energy and gradient
c     
      call dft_uks_energy(rtdb, g_tmp, eone, etwo, exc, enrep,
     $                    energy, edisp, g_grad, nexc)
c     
c     Tidy up
c
      if (.not. ga_destroy(g_tmp(1))) call errquit
     $     ('dft_uks_step_energy: tmp?', 0, GA_ERR)
      if (.not. ga_destroy(g_tmp(2))) call errquit
     $     ('dft_uks_step_energy: tmp?', 0, GA_ERR)
c
      if (oscfps) call pstat_off(ps_search)
c
      end
c $Id$
