c-----------------------------------------------------------------------
c
      subroutine wfn1_gen_mat(wfn1_inst,wfn1_mat)
      implicit none
!>
!> \brief calculate the 1-electron matrices
!>
!> Calculate the overlap, the 1-electron hamiltonian, and the 
!> Cauchy-Schwarz matrices.
!>
#include "wfn1_param.fh"
#include "wfn1_matcnst.fh"
#include "bas.fh"
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
!>
      type(wfn1_prminst), intent(in)    :: wfn1_inst
      type(wfn1_matcnst), intent(inout) :: wfn1_mat
c
      integer max1e ! The integral buffer size.
      integer mem1  ! The scratch memory size for the integral code.
c
      integer l_buf, k_buf ! The buffer memory segment
      integer l_scr, k_scr ! The scratch memory segment
c
      integer iproc   ! the processor rank
      integer nproc   ! the number of processors
      integer nshell  ! the number of shells
      integer ishell  ! the i-shell number
      integer jshell  ! the j-shell number
c
      integer ilo, ihi ! the basis function limits
      integer jlo, jhi ! the basis function limits
      integer klo, khi ! the local matrix block limits
      integer llo, lhi ! the local matrix block limits
c
      character*13 pname
      pname = "wfn1_gen_mat:"
c
      iproc = ga_nodeid()
      nproc = ga_nnodes()
c
      call ga_zero(wfn1_mat%g_s)
      call ga_zero(wfn1_mat%g_tv)
      call ga_zero(wfn1_mat%g_cs)
      call ga_distribution(wfn1_mat%g_s,iproc,klo,khi,llo,lhi)
c
      if (.not.bas_numcont(wfn1_inst%basis,nshell)) call errquit
     &  (pname//" bas_numcont failed for basis",wfn1_inst%basis,
     &   BASIS_ERR)
c
      call int_mem_1e(max1e,mem1)
      if (.not. MA_push_get(MT_DBL,max1e,pname//'buf',l_buf,k_buf))
     &  call errquit(pname//' ma buf failed', max1e, MA_ERR)
      if (.not. MA_push_get(MT_DBL,mem1,pname//'scr',l_scr,k_scr))
     &  call errquit(pname//' ma scr failed', mem1, MA_ERR)
c
      do jshell = 1, nshell
        if (.not.bas_cn2bfr(wfn1_inst%basis,jshell,jlo,jhi)) then
          call errquit(pname//' bas_cn2bfr',jshell,BASIS_ERR)
        endif
        if (klo.le.jlo.and.jlo.le.khi) then
          do ishell = 1, nshell
            if (.not.bas_cn2bfr(wfn1_inst%basis,ishell,ilo,ihi)) then
              call errquit(pname//' bas_cn2bfr',ishell,BASIS_ERR)
            endif
            if (llo.le.ilo.and.ilo.le.lhi) then
c
c             Overlap
c
              call int_1eov(wfn1_inst%basis,jshell,
     &                      wfn1_inst%basis,ishell,
     &                      mem1,dbl_mb(k_scr),
     &                      max1e,dbl_mb(k_buf))
              call ga_acc(wfn1_mat%g_s,ilo,ihi,jlo,jhi,dbl_mb(k_buf),
     &                    ihi-ilo+1,1.0d0)
c
c             Kinetic energy
c
              call int_1eke(wfn1_inst%basis,jshell,
     &                      wfn1_inst%basis,ishell,
     &                      mem1,dbl_mb(k_scr),
     &                      max1e,dbl_mb(k_buf))
              call ga_acc(wfn1_mat%g_tv,ilo,ihi,jlo,jhi,dbl_mb(k_buf),
     &                    ihi-ilo+1,1.0d0)
c
c             Nuclear attraction integrals
c
              call int_1epe(wfn1_inst%basis,jshell,
     &                      wfn1_inst%basis,ishell,
     &                      mem1,dbl_mb(k_scr),
     &                      max1e,dbl_mb(k_buf))
              call ga_acc(wfn1_mat%g_tv,ilo,ihi,jlo,jhi,dbl_mb(k_buf),
     &                    ihi-ilo+1,1.0d0)
c
c             Cauchy-Schwarz integrals (skipped for now)
c
            endif
          enddo
        endif
      enddo
c
      if (.not. MA_pop_stack(l_scr)) call errquit
     &  (pname//' pop scr failed', 0, MA_ERR)
      if (.not. MA_pop_stack(l_buf)) call errquit
     &  (pname//' pop buf failed', 0, MA_ERR)
c
      end
c
c-----------------------------------------------------------------------
