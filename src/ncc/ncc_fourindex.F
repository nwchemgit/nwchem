c===============================================================c
c                                                               c
c       NCC - a New Coupled-Cluster code for NWChem             c
c                                                               c
c       Developed by:                                           c
c                                                               c
c               Jeff R. Hammond                                 c
c               Leadership Computing Facility                   c
c               Argonne National Laboratory                     c
c               jhammond@mcs.anl.gov                            c
c                                                               c
c               Karol Kowalski                                  c
c               Environmental Molecular Sciences Laboratory     c
c               Pacific Northwest National Laboratory           c
c               karol.kowalski@pnl.gov                          c
c                                                               c
c               Marta Włoch                                     c
c               Department of Chemistry                         c
c               Michigan Technological University               c
c               wloch@mtu.edu                                   c
c                                                               c
c===============================================================c
      subroutine ncc_fourindex(rtdb,basis,
     1                         noa,nob,nva,nvb,nmoa,nmob,nbf,
     2                         v2oooo,v2ooov,v2oovv,v2ovvv,v2vvvv,
     3                         moveca,movecb)
c
c $Id$
c
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "inp.fh"
#include "bas.fh"
#include "geom.fh"
#include "sym.fh"
#include "schwarz.fh"
c
c     object handles
c
      integer rtdb             ! RTDB handle
      integer basis            ! BASIS handle
c
      integer i,j,k,l,offset
      integer nbf,nshells      ! Number of AO basis functions and shells
      integer noa,nob,nva,nvb,nmoa,nmob
c
      integer k_ao2e,l_ao2e
      integer k_qtr1,l_qtr1
      integer k_half,l_half
      integer k_qtr3,l_qtr3
      integer k_mo2e,l_mo2e
c
      double precision moveca(nbf,nmoa)
      double precision movecb(nbf,nmob)
c
      double precision v2oooo(noa,noa,noa,noa)
      double precision v2ooov(noa,noa,noa,nva)
      double precision v2oovv(noa,noa,nva,nva)
      double precision v2ovvv(noa,nva,nva,nva)
      double precision v2vvvv(nva,nva,nva,nva)
c
#ifdef DEBUG_PRINT
      print*,'top of ncc_fourindex'
#endif
c
c     ==========
c     Allocate
c     ==========
c
      if (.not.ma_push_get(mt_dbl,nmoa*nmoa*nmoa*nmoa,
     1                     'mo2e',l_mo2e,k_mo2e)) then
          call errquit ('ncc_fourindex: ma_push_get mo2e',
     1                  noa*nmoa*nmoa*nmoa,MA_ERR)
      endif
c
      if (.not.ma_push_get(mt_dbl,nbf*nmoa*nmoa*nmoa,
     1                     'qtr3',l_qtr3,k_qtr3)) then
          call errquit ('ncc_fourindex: ma_push_get qtr3',
     1                  nbf*nmoa*nmoa*nmoa,MA_ERR)
      endif
c
      if (.not.ma_push_get(mt_dbl,nbf*nbf*nmoa*nmoa,
     1                     'half',l_half,k_half)) then
          call errquit ('ncc_fourindex: ma_push_get half',
     1                  nbf*nbf*nmoa*nmoa,MA_ERR)
      endif
c
      if (.not.ma_push_get(mt_dbl,nbf*nbf*nbf*nmoa,
     1                     'qtr1',l_qtr1,k_qtr1)) then
          call errquit ('ncc_fourindex: ma_push_get qtr1',
     1                  nbf*nbf*nbf*nmoa,MA_ERR)
      endif
c
      if (.not.ma_push_get(mt_dbl,nbf*nbf*nbf*nbf,
     1                     'ao2e',l_ao2e,k_ao2e)) then
          call errquit ('ncc_fourindex: ma_push_get ao2e',
     1                  nbf*nbf*nbf*nbf,MA_ERR)
      endif
c
c     ==========
c     Fill ao2e array
c     ==========
c
      call ncc_make_ao2e(rtdb,basis,nbf,
     1                   dbl_mb(k_ao2e))
c
c     ==========
c     Do first half-transformation
c     ==========
c
      call ncc_half_trans1(nmoa,nmob,nbf,
     1                     dbl_mb(k_ao2e),
     2                     dbl_mb(k_qtr1),
     3                     dbl_mb(k_half),
     4                     moveca,movecb)
c
c     ==========
c     Do second half-transformation
c     ==========
c
      call ncc_half_trans2(nmoa,nmob,nbf,
     1                     dbl_mb(k_half),
     2                     dbl_mb(k_qtr3),
     3                     dbl_mb(k_mo2e),
     4                     moveca,movecb)
c
c      call ncc_print_mo2e(nmoa,nmob,dbl_mb(k_mo2e))
c
c     ==========
c     Deallocate
c     ==========
c
      if (.not.ma_pop_stack(l_ao2e)) then
        call errquit('ncc_fourindex: MA problem ao2e ',0,MA_ERR)
      endif
c
      if (.not.ma_pop_stack(l_qtr1)) then
        call errquit('ncc_fourindex: MA problem qtr1 ',0,MA_ERR)
      endif
c
      if (.not.ma_pop_stack(l_half)) then
        call errquit('ncc_fourindex: MA problem half ',0,MA_ERR)
      endif
c
      if (.not.ma_pop_stack(l_qtr3)) then
        call errquit('ncc_fourindex: MA problem qtr3 ',0,MA_ERR)
      endif
c
      if (.not.ma_pop_stack(l_mo2e)) then
        call errquit('ncc_fourindex: MA problem mo2e ',0,MA_ERR)
      endif
c
#ifdef DEBUG_PRINT
      print*,'end of ncc_fourindex'
#endif
c
      return
c
      end





      subroutine ncc_half_trans1(nmoa,nmob,nbf,
     1                           ao2e,qtr1,half,Ua,Ub)
      implicit none
c
      integer i,j,k,l ! MO indices
      integer a,b,c,d ! AO indices
      integer nbf,nmoa,nmob
c
      double precision ao2e(nbf,nbf,nbf,nbf)
      double precision qtr1(nmoa,nbf,nbf,nbf)
      double precision half(nmoa,nmoa,nbf,nbf)
c
      double precision Ua(nbf,nmoa) ! moveca
      double precision Ub(nbf,nmob) ! movecb
c
#ifdef DEBUG_PRINT
      print*,'top of ncc_half_trans1'
#endif
c
      do a=1,nbf
        do b=1,nbf
          do c=1,nbf
            do i=1,nmoa
              qtr1(i,c,b,a) = 0.0d0
              do d=1,nbf
                qtr1(i,c,b,a) = qtr1(i,c,b,a) + Ua(d,i)*ao2e(d,c,b,a)
              enddo
            enddo
          enddo
        enddo
      enddo
c
      do a=1,nbf
        do b=1,nbf
          do j=1,nmoa
            do i=1,nmoa
              half(i,j,b,a) = 0.0d0
              do c=1,nbf
                half(i,j,b,a) = half(i,j,b,a) + Ua(c,j)*qtr1(i,c,b,a)
              enddo
            enddo
          enddo
        enddo
      enddo
c
#ifdef DEBUG_PRINT
      print*,'end of ncc_half_trans1'
#endif
c
      return
c
      end



      subroutine ncc_half_trans2(nmoa,nmob,nbf,
     1                           half,qtr3,mo2e,Ua,Ub)
      implicit none
c
      integer i,j,k,l ! MO indices
      integer a,b,c,d ! AO indices
      integer nbf
      integer noa,nob,nva,nvb,nmoa,nmob
c
      double precision half(nmoa,nmoa,nbf,nbf)
      double precision qtr3(nmoa,nmoa,nmoa,nbf)
      double precision mo2e(nmoa,nmoa,nmoa,nmoa)
c
      double precision Ua(nbf,nmoa) ! moveca
      double precision Ub(nbf,nmob) ! movecb
c
#ifdef DEBUG_PRINT
      print*,'top of ncc_half_trans2'
#endif
c
      do a=1,nbf
        do k=1,nmoa
          do j=1,nmoa
            do i=1,nmoa
              qtr3(i,j,k,a) = 0.0d0
              do b=1,nbf
                qtr3(i,j,k,a) = qtr3(i,j,k,a) + Ua(b,k)*half(i,j,b,a)
              enddo
            enddo
          enddo
        enddo
      enddo
c
      do l=1,nmoa
        do k=1,nmoa
          do j=1,nmoa
            do i=1,nmoa
              mo2e(i,j,k,l) = 0.0d0
              do a=1,nbf
                mo2e(i,j,k,l) = mo2e(i,j,k,l) + Ua(a,l)*qtr3(i,j,k,a)
              enddo
            enddo
          enddo
        enddo
      enddo
c
#ifdef DEBUG_PRINT
      print*,'end of ncc_half_trans2'
#endif
c
      return
c
      end









      subroutine ncc_make_ao2e(rtdb,basis,nbf,ao2e)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "bas.fh"
#include "schwarz.fh"
c
c     object handles
c
      integer rtdb             ! RTDB handle
      integer basis            ! BASIS handle
c
      integer i,j,k,l
      integer nbf,nshells      ! Number of AO basis functions and shells
c
      integer s_int2e,k_int2e,l_int2e
      integer s_scr2e,k_scr2e,l_scr2e
      integer ish,ilo,ihi,irng
      integer jsh,jlo,jhi,jrng
      integer ksh,klo,khi,krng
      integer lsh,llo,lhi,lrng
c
      double precision tol2e
      double precision schwarz_ij,schwarz_kl
c
      double precision ao2e(nbf,nbf,nbf,nbf)
c
#ifdef DEBUG_PRINT
      print*,'top of ncc_make_ao2e'
#endif
c
c     =======================
c     basis information
c     =======================
c
      if (.not.rtdb_get(rtdb,'ncc:tol2e',mt_dbl,1,tol2e)) then
        tol2e = 1.0d-10
      endif
      call int_acc_set(tol2e)
c
      if (.not.bas_numcont(basis,nshells)) then
        call errquit('ncc_make_ao2e: failed to read nshells',0,
     1               BASIS_ERR)
      endif
c
c     =======================
c     Initialize integral memory
c     =======================
c
      call int_mem_2e4c(s_int2e,s_scr2e)
      if (.not.ma_push_get(mt_dbl,s_int2e,'int2e',l_int2e,k_int2e)) then
        call errquit('ncc_make_ao2e: MA problem int2e',s_int2e,MA_ERR)
      endif
      if (.not.ma_push_get(mt_dbl,s_scr2e,'scr2e',l_scr2e,k_scr2e)) then
        call errquit('ncc_make_ao2e: MA problem scr2e',s_scr2e,MA_ERR)
      endif
c
c     ================
c     Loop over shells
c     ================
c
      do ish=1,nshells
        if (.not.bas_cn2bfr(basis,ish,ilo,ihi)) then
            call errquit('ncc_make_ao2e: bas_cn2bfr i',i,BASIS_ERR)
        endif
        irng = ihi - ilo + 1
        do jsh=1,nshells
          if (.not.bas_cn2bfr(basis,jsh,jlo,jhi)) then
            call errquit('ncc_make_ao2e: bas_cn2bfr j',j,BASIS_ERR)
          endif
          jrng = jhi - jlo + 1
          schwarz_ij = schwarz_shell(ish,jsh)
          do ksh=1,nshells
            if (.not.bas_cn2bfr(basis,ksh,klo,khi)) then
              call errquit('ncc_make_ao2e: bas_cn2bfr k',k,BASIS_ERR)
            endif
            krng = khi - klo + 1
            do lsh=1,nshells
              if (.not.bas_cn2bfr(basis,lsh,llo,lhi)) then
                call errquit('ncc_make_ao2e: bas_cn2bfr l',l,BASIS_ERR)
              endif
              lrng = lhi - llo + 1
              schwarz_kl = schwarz_shell(ksh,lsh)
              if ((schwarz_ij*schwarz_kl).ge.tol2e) then
                call int_2e4c(basis,ish,jsh,basis,ksh,lsh,
     1                        s_scr2e,dbl_mb(k_scr2e),
     2                        s_int2e,dbl_mb(k_int2e))
                call put_ao2e_to_buf(ilo,ihi,jlo,jhi,
     1                               klo,khi,llo,lhi,nbf,
     2                               dbl_mb(k_int2e),ao2e)

! #ifdef DEBUG_PRINT
!                 write(6,*) '========================================='
!                 write(6,1) ish,jsh,ksh,lsh
!                 write(6,2) irng,jrng,krng,lrng
!  1              format('(ish,jsh,ksh,lsh) = ',4i5)
!  2              format('(irng,jrng,krng,lrng) = ',4i5)
!                 call ncc_print_ao2e(ilo,ihi,jlo,jhi,klo,khi,llo,lhi,
!      1                              dbl_mb(k_int2e),0.0d0)
! #endif
              endif
            enddo
          enddo
        enddo
      enddo
c
c
c
! #ifdef DEBUG_PRINT
!       do i=1,nbf
!         do j=1,nbf
!           do k=1,nbf
!             do l=1,nbf
!                 write(6,3) i,j,k,l,ao2e(l,k,j,i)
!  3              format('ao2e(',i5,',',i5,',',i5,',',i5,') = ',f24.14)
!             enddo
!           enddo
!         enddo
!       enddo
! #endif
c
c     ==========
c     Deallocate
c     ==========
c
      if (.not.ma_pop_stack(l_scr2e)) then
        call errquit('ncc_make_ao2e: MA problem scr2e ',0,MA_ERR)
      endif
      if (.not.ma_pop_stack(l_int2e)) then
        call errquit('ncc_make_ao2e: MA problem int2e ',0,MA_ERR)
      endif
c
#ifdef DEBUG_PRINT
      print*,'end of ncc_make_ao2e'
#endif
c
      return
c
      end
