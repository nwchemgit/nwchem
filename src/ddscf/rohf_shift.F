      subroutine rohf_make_shift(shift)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "crohf.fh"
#include "cscf.fh"
c
c     Return minimum shift necessary for diagonal approximation
c     to 1-e Hessian to be positive definite
c
      double precision shift
c
      integer l_diag, k_diag
      double precision shift_cv, shift_pv, shift_cp
c
      if (.not. ma_push_get(mt_dbl, nbf, 'rohf_m_s: diag', 
     $     l_diag, k_diag)) call errquit
     $     ('rohf_make_shift: insufficient memory', nbf)
c
      call ga_sync()
c
      shift_cv = 0.0d0
      shift_pv = 0.0d0
      shift_cp = 0.0d0
c
      call rohf_make_shift_2(crohf_g_fcv, dbl_mb(k_diag),
     $     1, nclosed, nclosed+nopen+1, nbf,
     $     4.0d0, shift_cv)
      if (nopen .gt. 0) then
         call rohf_make_shift_2(crohf_g_fpv, dbl_mb(k_diag),
     $        nclosed+1, nclosed+nopen, nclosed+nopen+1, nbf, 
     $        2.0d0, shift_pv)
         call rohf_make_shift_2(crohf_g_fcv, dbl_mb(k_diag),
     $        1, nclosed, nclosed+1, nclosed+nopen, 
     $        2.0d0, shift_cp)
      endif
c
      shift = max(shift_cv, shift_pv, shift_cp)
c
      if (.not. ma_pop_stack(l_diag)) call errquit('rohf_n_s: e',0)
c
      end
      subroutine rohf_make_shift_2(g_fock, diag,
     $     ilo, ihi, alo, ahi, scale, shift)
      implicit none
c
c     return scale * min(i,a) (F(a) - F(i))
c     for i/a in given ranges
c
      integer g_fock
      integer ilo, ihi, alo, ahi
      double precision diag(*), scale, shift
c
      double precision minocc, minvirt
      integer i, a
c
      call ga_get_diagonal(g_fock, diag)
c
      minocc = -diag(ilo)
      do i = ilo, ihi
         minocc = min(-diag(i),minocc)
      enddo
      minvirt = diag(alo)
      do a = alo, ahi
         minvirt = min(diag(a),minvirt)
      enddo
c
      shift = scale * abs(min(0.0d0,minocc + minvirt))
c
      end
