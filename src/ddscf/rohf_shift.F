      subroutine rohf_make_shift(shift)
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "crohf.fh"
#include "cscf.fh"
c
c     Return minimum shift necessary for diagonal approximation
c     to 1-e Hessian to be positive definite
c
c     Now modified to exclude symmetry forbidden rotations if
c     symmetry adaption is enabled
c
      double precision shift
c
      integer l_diag, k_diag
      double precision shift_cv, shift_pv, shift_cp
c
      if (.not. ma_push_get(mt_dbl, nmo, 'rohf_m_s: diag', 
     $     l_diag, k_diag)) call errquit
     $     ('rohf_make_shift: insufficient memory', nmo, MA_ERR)
c
      call ga_sync()
c
      shift_cv = 0.0d0
      shift_pv = 0.0d0
      shift_cp = 0.0d0
c
      call rohf_make_shift_2(crohf_g_fcv, dbl_mb(k_diag),
     $     1, nclosed, nclosed+nopen+1, nmo, 
     $     oadapt, int_mb(k_irs), 4.0d0, shift_cv)
      if (nopen .gt. 0) then
         call rohf_make_shift_2(crohf_g_fpv, dbl_mb(k_diag),
     $        nclosed+1, nclosed+nopen, nclosed+nopen+1, nmo, 
     $        oadapt, int_mb(k_irs), 2.0d0, shift_pv)
*     I think that this is not necessary here and results in
*     open shell d/f systems being shifted uncessarily
*         call rohf_make_shift_2(crohf_g_fcv, dbl_mb(k_diag),
*     $        1, nclosed, nclosed+1, nclosed+nopen, 
*     $        oadapt, int_mb(k_irs), 2.0d0, shift_cp)
      endif
c
      shift = max(shift_cv, shift_pv, shift_cp)
c
      if (.not. ma_pop_stack(l_diag)) call errquit('rohf_n_s: e',0,
     &       MA_ERR)
c
      end
      subroutine rohf_make_shift_2(g_fock, diag,
     $     ilo, ihi, alo, ahi, oadapt, irs, scale, shift)
      implicit none
*
c     return scale * min(i,a) (F(a) - F(i))
c     for i/a in given ranges and of same symmetry
c
      integer g_fock
      integer ilo, ihi, alo, ahi
      logical oadapt
      integer irs(*)
      double precision diag(*), scale, shift
c
      integer i, a
c
      call ga_get_diagonal(g_fock, diag)
c
      shift = 0.0d0
      if (oadapt) then
        do i = ilo, ihi
            do a = alo, ahi
                if (irs(i).eq.irs(a))
     $             shift = min(shift, diag(a)-diag(i))
            enddo
        enddo
      else
        do i = ilo, ihi
            do a = alo, ahi
                shift = min(shift, diag(a)-diag(i))
            enddo
        enddo
      endif
c
      shift = - scale * shift
c
      end
