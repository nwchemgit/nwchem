      subroutine rohf_canon(oaufbau, oprint)
C$Id: rohf_canon.F,v 1.2 1995-10-21 19:39:33 d3g681 Exp $
      implicit none
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "cscfps.fh"
#include "cscf.fh"
#include "crohf.fh"
      logical oaufbau
      logical oprint
c
      integer g_u, g_fock, g_tmp
      double precision one, zero
      data one, zero/1.d0, 0.d0/
c
c     This routine assumes that rohf_energy/rohf_fock have been called
c     so that the contents of /crohf/ are current.
c
c     Diagonalize the ROHF 'Fock' matrix
c     
c     If (oaufbau) 
c        diagonalize the whole thing and allow mixing of closed-open-virt
c     else
c        diagonalize separately the cloed-closed, open-open, and 
c        virt-virt parts
c
c     Transform Fock matrices and MO coefficients into the new canonical basis
c
      if (oscfps) call pstat_on(ps_diag)
c
      if (.not. ga_create(MT_DBL, nbf, nbf, 'rohf_canon: u',
     $     1, 1, g_u)) call errquit('rohf_canon: ga failed for u', 0)
      call ga_zero(g_u)
      if (.not. ga_create(MT_DBL, nbf, nbf, 'rohf_canon: fock',
     $     1, 1, g_fock)) call errquit
     $     ('rohf_canon: ga failed for fock', 0)
c
      call rohf_get_fock(g_fock)
c
      if (oaufbau) then
#if defined(PARALLEL_DIAG)
         call ga_diag_std(g_fock, g_u, dbl_mb(k_eval))
#else 
         call ga_diag_std_seq(g_fock, g_u, dbl_mb(k_eval))
#endif
      else
c
c     closed-closed piece
c
         call rohf_canon_subspace(g_fock, g_u, dbl_mb(k_eval),
     $        1, nclosed)
c
         if (nopen .gt. 0) then
c
c     open-open piece
c
            call rohf_canon_subspace(g_fock, g_u,
     $           dbl_mb(k_eval+nclosed),
     $           nclosed+1, nclosed+nopen)
         endif
c     
c     virt-virt piece
c     
         call rohf_canon_subspace(g_fock, g_u,
     $        dbl_mb(k_eval+nclosed+nopen),
     $        nclosed+nopen+1, nbf)

      endif
c
      if (oscfps) call pstat_off(ps_diag)
c
      call movecs_fix_phase(g_u)
c
c     Apply rotation to orbitals and fock matrix
c
      if (.not. ga_create(MT_DBL, nbf, nbf, 'rohf_canon: tmp',
     $     1, 1, g_tmp)) call errquit
     $     ('rohf_canon: ga failed for tmp', 0)
c
      call ga_copy(g_movecs, g_tmp)
      call ga_dgemm('n', 'n', nbf, nbf, nbf, one, g_tmp, g_u,
     $     zero, g_movecs)
c
      call ga_dgemm('n', 'n', nbf, nbf, nbf, one, crohf_g_fcv, g_u,
     $     zero, g_tmp)
      call ga_dgemm('t', 'n', nbf, nbf, nbf, one, g_u, g_tmp, 
     $     zero, crohf_g_fcv)
c
      if (nopen .gt. 0) then
         call ga_dgemm('n', 'n', nbf, nbf, nbf, one, crohf_g_fcp, g_u,
     $        zero, g_tmp)
         call ga_dgemm('t', 'n', nbf, nbf, nbf, one, g_u, g_tmp, 
     $        zero, crohf_g_fcp)
         call ga_dgemm('n', 'n', nbf, nbf, nbf, one, crohf_g_fpv, g_u,
     $        zero, g_tmp)
         call ga_dgemm('t', 'n', nbf, nbf, nbf, one, g_u, g_tmp, 
     $        zero, crohf_g_fpv)
      endif
c
      if (oprint .and. ga_nodeid().eq.0) then
         write(6,*)
         write(6,*)
         call util_print_centered(6, 'Eigenvalues', 20, .true.)
         call output(dbl_mb(k_eval), 1, min(nclosed+nopen+5,nbf), 
     $        1, 1, nbf, 1, 1)
         call util_flush(6)
      endif
c
      if (.not. ga_destroy(g_u))
     $     call errquit('rohf_canon: destroy', 0)
      if (.not. ga_destroy(g_fock))
     $     call errquit('rohf_canon: destroy', 0)
      if (.not. ga_destroy(g_tmp))
     $     call errquit('rohf_canon: destroy', 0)
c
      end
      subroutine rohf_canon_subspace(g_fock, g_u, evals, ilo, ihi)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
      integer g_fock, g_u, ilo, ihi
      double precision evals(*)
c
      integer n
      integer g_v
      integer g_tmp
c     
      n = ihi - ilo + 1
c
      if (.not. ga_create(MT_DBL, n, n, 'rohf_canon: tmp',
     $     1, 1, g_tmp)) call errquit
     $     ('rohf_canon: ga failed for tmp', 0)
      if (.not. ga_create(MT_DBL, n, n, 'rohf_canon: v',
     $     1, 1, g_v)) call errquit
     $     ('rohf_canon: ga failed for v', 0)
c
      call ga_copy_patch('n', 
     $     g_fock, ilo, ihi, ilo, ihi, 
     $     g_tmp, 1, n, 1, n)
c
#if defined(PARALLEL_DIAG)
      call ga_diag_std(g_tmp, g_v, evals)
#else 
      call ga_diag_std_seq(g_tmp, g_v, evals)
#endif
      call ga_copy_patch('n', 
     $     g_v, 1, n, 1, n,
     $     g_u, ilo, ihi, ilo, ihi)
      if (.not. (ga_destroy(g_tmp) .and. ga_destroy(g_v)))
     $     call errquit('rohf_canon: ga_destroy ?', 0)
c
      end

