      subroutine rohf_step_energy(rtdb, step, g_search,
     $     eone, etwo, enrep, energy, g_grad)
*
* $Id: rohf_step.F,v 1.5 2000-05-23 05:39:07 mdupuis Exp $
*
      implicit none
#include "crohf.fh"
#include "cscfps.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "cscf.fh"
c
      integer rtdb
      double precision step
      integer g_search
      double precision eone, etwo, enrep, energy
      integer g_grad
c
      double precision scale
c
c     Compute the energy and gradient at the given step
c
      integer g_tmp
      if (oscfps) call pstat_on(ps_search)
c
c     Generate MOs * Exp(-step*K) in temporary array
c     
      if (.not. ga_create(MT_DBL, nbf, nmo, 'rohf_s_e: MOs',
     $     32, 32, g_tmp)) call errquit('rohf_s_e: MOs', 0)
c      
      if (step .ne. 0.0d0) then
         call ga_dscal(g_search, step)
         call rohf_k2cf(basis,
     $        nbf, nmo, nclosed, nopen,
     $        g_search, g_movecs, g_tmp)
         scale = 1.0d0 / step
         call ga_dscal(g_search, scale)
      else
         call ga_copy(g_movecs, g_tmp)
      endif
c
c     Compute energy and gradient
c     
      call rohf_energy( g_tmp, eone, etwo, enrep, energy, g_grad,
     $                  rtdb)
c     
c     Tidy up
c
      if (.not. ga_destroy(g_tmp)) call errquit
     $     ('rohf_nr_solve: tmp?', 0)
c
      if (oscfps) call pstat_off(ps_search)
c
      end
