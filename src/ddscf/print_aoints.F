      subroutine print_aoints(basis)
C$Id$
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "util.fh"
c
      integer basis
c
      integer ish, jsh, ksh, lsh, mem2, max2e, l_g, k_g, l_scr, k_scr,
     $     nsh, ilo, ihi, jlo, jhi, klo, khi, llo, lhi
      logical status
c
      status = bas_numcont(basis, nsh)
      if (.not. status) call errquit('print_aoints: info failed', 0,
     &       BASIS_ERR)
      call int_mem_2e4c(max2e, mem2)
c
c     scratch space needed for 2-electron routines
c     
      status = .true.
      status = status .and. 
     $     ma_push_get(MT_DBL, max2e, 'print_aoints: buf', l_g, k_g)
      status = status .and.
     $     ma_push_get(MT_DBL, mem2, 'print_aoints: scr', l_scr, k_scr)
      if (.not. status)
     $     call errquit('print_aoints: insufficient memory', mem2,
     &       MA_ERR)
c
c     integrals evaluated to full precision
c
      call int_acc_set(1d-30)
c     
      do ish = 1, nsh
         if (.not. bas_cn2bfr(basis, ish, ilo, ihi))
     $        call errquit('print_aoints: bas_cn2bfr', ish, BASIS_ERR)
         do jsh = 1, nsh
            if (.not. bas_cn2bfr(basis, jsh, jlo, jhi))
     $           call errquit('print_aoints: bas_cn2bfr', jsh,
     &       BASIS_ERR)
            do ksh = 1, nsh
               if (.not. bas_cn2bfr(basis, ksh, klo, khi))
     $              call errquit('print_aoints: bas_cn2bfr', ksh,
     &       BASIS_ERR)
               do lsh = 1, nsh
                  if (.not. bas_cn2bfr(basis, lsh, llo, lhi))
     $                 call errquit('print_aoints: bas_cn2bfr', lsh,
     &       BASIS_ERR)
c
                  call int_2e4c(basis, ish, jsh, basis, ksh, lsh,
     $                 mem2, dbl_mb(k_scr), max2e, dbl_mb(k_g))
c
                  call print_ao_buf(ilo, ihi, jlo, jhi, klo, khi, 
     $                 llo, lhi, dbl_mb(k_g))
c
               enddo
            enddo
         enddo
      enddo
c
c
      if (.not. ma_pop_stack(l_scr)) call errquit('print_aoints: pop',0,
     &       MA_ERR)
      if (.not. ma_pop_stack(l_g)) call errquit('print_aoints: pop', 1,
     &       MA_ERR)
c
      call int_acc_std()
c
      end
      subroutine print_ao_buf(ilo, ihi, jlo, jhi, klo, khi, 
     $     llo, lhi, g)
      implicit none
      integer ilo, ihi, jlo, jhi, klo, khi, llo, lhi
      double precision g(llo:lhi,klo:khi,jlo:jhi,ilo:ihi)
c
      integer i, j, k, l
c
      do i = ilo, ihi
         do j = jlo, jhi
            do k = klo, khi
               do l = llo, lhi
                  if (abs(g(l,k,j,i)).gt.1d-8) then
                     write(6,1) i,j,k,l,g(l,k,j,i)
 1                   format(' ao ',4i5,f20.8)
                  endif
               enddo
            enddo
         enddo
      enddo
c
      end
