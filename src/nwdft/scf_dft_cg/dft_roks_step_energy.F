      subroutine dft_roks_step_energy(rtdb, step, g_search,
     $     eone, etwo, exc, enrep, energy, edisp, g_grad, nexc)
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "crohf.fh"
#include "cscfps.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "cscf.fh"
c
      integer rtdb
      double precision step
      integer g_search
      double precision eone, etwo, enrep, energy, exc(2), edisp
      integer g_grad
      integer nexc
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
      if (.not. ga_create(MT_DBL, nbf, nmo, 'dft_roks_step_energy: MOs',
     $     32, 32, g_tmp))
     $  call errquit('dft_roks_step_energy: MOs', 0, GA_ERR)
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
      call dft_roks_energy(g_tmp, eone, etwo, exc, enrep, energy,
     $                     edisp, g_grad, rtdb, nexc)
c     
c     Tidy up
c
      if (.not. ga_destroy(g_tmp)) call errquit
     $     ('dft_roks_step_energy: tmp?', 0, GA_ERR)
c
      if (oscfps) call pstat_off(ps_search)
c
      end
