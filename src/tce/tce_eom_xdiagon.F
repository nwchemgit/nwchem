      subroutine tce_eom_xdiagon(needx1,needx2,needx3,needx4,
     1  size_x1,size_x2,size_x3,size_x4,
     2  k_x1_offset,k_x2_offset,k_x3_offset,k_x4_offset,
     3  d_r1,d_r2,d_r3,d_r4,omega,residual)
c
c $Id: tce_eom_xdiagon.F,v 1.1 2003-07-25 17:56:56 sohirata Exp $
c
c     Form a subspace matrix and diagonalize
c
      implicit none
#include "stdio.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "tce.fh"
#include "util.fh"
      logical needx1
      logical needx2
      logical needx3
      logical needx4
      integer size_x1
      integer size_x2
      integer size_x3
      integer size_x4
      integer k_x1_offset
      integer k_x2_offset
      integer k_x3_offset
      integer k_x4_offset
      integer d_r1
      integer d_r2
      integer d_r3
      integer d_r4
      double precision omega(*)
      double precision residual(*)
      double precision overlap
      logical nodezero
      integer ivec, jvec
      integer l_hbar, k_hbar
      integer l_vr, k_vr
      integer l_vl, k_vl
      integer l_ei, k_ei
      integer l_wk, k_wk
      integer info
      integer newnxtrials
      double precision ddotfile
      external ddotfile
      double precision au2ev   ! Conversion factor from a.u. to eV
      parameter (au2ev=27.2113961d0)
      character*255 filename
c
      nodezero = (ga_nodeid().eq.0)
c
c     Allocate
c
      if (.not.ma_push_get(mt_dbl,nxtrials*nxtrials,'hbar',
     1  l_hbar,k_hbar)) call errquit('tce_eom_xdiagon: MA problem',0)
      if (.not.ma_push_get(mt_dbl,nxtrials*nxtrials,'vr',
     1  l_vr,k_vr)) call errquit('tce_eom_xdiagon: MA problem',1)
      if (.not.ma_push_get(mt_dbl,nxtrials,'vl',
     1  l_vl,k_vl)) call errquit('tce_eom_xdiagon: MA problem',2)
      if (.not.ma_push_get(mt_dbl,nxtrials,'ei',
     1  l_ei,k_ei)) call errquit('tce_eom_xdiagon: MA problem',4)
      if (.not.ma_push_get(mt_dbl,4*nxtrials,'work',
     1  l_wk,k_wk)) call errquit('tce_eom_xdiagon: MA problem',5)
c
c     Form similarity transformed Hamiltonian subspace matrix
c
      do ivec = 1, nxtrials
        do jvec = 1, nxtrials
          dbl_mb(k_hbar+(ivec-1)*nxtrials+jvec-1) = 0.0d0
          if (needx1) then
            if (.not.x1_exist(ivec))
     1        call errquit('tce_eom_xdiagon: x1 not found',ivec)
            if (.not.xp1_exist(jvec))
     1        call errquit('tce_eom_xdiagon: xp1 not found',jvec)
            dbl_mb(k_hbar+(ivec-1)*nxtrials+jvec-1) = 
     1      dbl_mb(k_hbar+(ivec-1)*nxtrials+jvec-1) + 
     2      ddotfile(x1(jvec),xp1(ivec),size_x1)
          endif
          if (needx2) then
            if (.not.x2_exist(ivec))
     1        call errquit('tce_eom_xdiagon: x2 not found',ivec)
            if (.not.xp2_exist(jvec))
     1        call errquit('tce_eom_xdiagon: xp2 not found',jvec)
            dbl_mb(k_hbar+(ivec-1)*nxtrials+jvec-1) = 
     1      dbl_mb(k_hbar+(ivec-1)*nxtrials+jvec-1) + 
     2      ddotfile(x2(jvec),xp2(ivec),size_x2)
          endif
          if (needx3) then
            if (.not.x3_exist(ivec))
     1        call errquit('tce_eom_xdiagon: x3 not found',ivec)
            if (.not.xp3_exist(jvec))
     1        call errquit('tce_eom_xdiagon: xp3 not found',jvec)
            dbl_mb(k_hbar+(ivec-1)*nxtrials+jvec-1) = 
     1      dbl_mb(k_hbar+(ivec-1)*nxtrials+jvec-1) + 
     2      ddotfile(x3(jvec),xp3(ivec),size_x3)
          endif
          if (needx4) then
            if (.not.x4_exist(ivec))
     1        call errquit('tce_eom_xdiagon: x4 not found',ivec)
            if (.not.xp4_exist(jvec))
     1        call errquit('tce_eom_xdiagon: xp4 not found',jvec)
            dbl_mb(k_hbar+(ivec-1)*nxtrials+jvec-1) = 
     1      dbl_mb(k_hbar+(ivec-1)*nxtrials+jvec-1) + 
     2      ddotfile(x4(jvec),xp4(ivec),size_x4)
          endif
        enddo
      enddo
      if (nodezero.and.util_print('hbar',print_debug)) then
        write(LuOut,*)
        call ma_print(dbl_mb(k_hbar),nxtrials,nxtrials,
     1    "Similarity Transformed Hamiltonian")
      endif
c     call dcopy(nxtrials**2,dbl_mb(k_hbar),1,dbl_mb(k_gbar),1)
c
c     Diagonalize the non-Hermitian matrix
c
      call tddft_dgeev('n','v',nxtrials,dbl_mb(k_hbar),nxtrials,
     1  omega,dbl_mb(k_ei),dbl_mb(k_vl),1,dbl_mb(k_vr),nxtrials,
     2  dbl_mb(k_wk),4*nxtrials,info)
      if (info.ne.0) call errquit('tce_eom_xdiagon: dgeev failed',info)
c     call tddft_schmidt(nxtrials,dbl_mb(k_vr))
      call tddft_sort(nxtrials,omega,dbl_mb(k_vr),'A')
      if (nodezero.and.util_print('omega',print_debug)) then
        do ivec = 1, nxtrials
          write(LuOut,9000) omega(ivec),omega(ivec)*au2ev
        enddo
      endif
      if (nodezero.and.util_print('vectors',print_debug)) then
        call ma_print(dbl_mb(k_vr),nxtrials,nxtrials,'Eigenvectors')
      endif
c
c     Compute residuals
c
      newnxtrials = nxtrials
      do ivec = 1, nroots_reduced
        residual(ivec) = 0.0d0
        if (needx1) then
          call tce_zero(d_r1,size_x1)
          call daxpyfile(nxtrials,dbl_mb(k_vr+(ivec-1)*nxtrials),
     1      xp1,d_r1,size_x1)
          call dfill(nxtrials,0.0d0,dbl_mb(k_ei),1)    
          call daxpy(nxtrials,-omega(ivec),
     1      dbl_mb(k_vr+(ivec-1)*nxtrials),1,dbl_mb(k_ei),1)
          call daxpyfile(nxtrials,dbl_mb(k_ei),
     1      x1,d_r1,size_x1)
          residual(ivec) = residual(ivec) 
     1                   + ddotfile(d_r1,d_r1,size_x1)
          if (nodezero.and.util_print('residual',print_debug)) then
            write(LuOut,*) ' X1: vector',ivec
            write(LuOut,*) ' X1: residual',residual(ivec)
            call sf_print(d_r1,size_x1)
          endif
        endif
        if (needx2) then
          call tce_zero(d_r2,size_x2)
          call daxpyfile(nxtrials,dbl_mb(k_vr+(ivec-1)*nxtrials),
     1      xp2,d_r2,size_x2)
          call dfill(nxtrials,0.0d0,dbl_mb(k_ei),1)    
          call daxpy(nxtrials,-omega(ivec),
     1      dbl_mb(k_vr+(ivec-1)*nxtrials),1,dbl_mb(k_ei),1)
          call daxpyfile(nxtrials,dbl_mb(k_ei),
     1      x2,d_r2,size_x2)
          residual(ivec) = residual(ivec) 
     1                   + ddotfile(d_r2,d_r2,size_x2)
          if (nodezero.and.util_print('residual',print_debug)) then
            write(LuOut,*) ' X2: vector',ivec
            write(LuOut,*) ' X2: residual',residual(ivec)
            call sf_print(d_r2,size_x2)
          endif
        endif
        if (needx3) then
          call tce_zero(d_r3,size_x3)
          call daxpyfile(nxtrials,dbl_mb(k_vr+(ivec-1)*nxtrials),
     1      xp3,d_r3,size_x3)
          call dfill(nxtrials,0.0d0,dbl_mb(k_ei),1)    
          call daxpy(nxtrials,-omega(ivec),
     1      dbl_mb(k_vr+(ivec-1)*nxtrials),1,dbl_mb(k_ei),1)
          call daxpyfile(nxtrials,dbl_mb(k_ei),
     1      x3,d_r3,size_x3)
          residual(ivec) = residual(ivec) 
     1                   + ddotfile(d_r3,d_r3,size_x3)
          if (nodezero.and.util_print('residual',print_debug)) then
            write(LuOut,*) ' X3: vector',ivec
            write(LuOut,*) ' X3: residual',residual(ivec)
            call sf_print(d_r3,size_x3)
          endif
        endif
        if (needx4) then
          call tce_zero(d_r4,size_x4)
          call daxpyfile(nxtrials,dbl_mb(k_vr+(ivec-1)*nxtrials),
     1      xp4,d_r4,size_x4)
          call dfill(nxtrials,0.0d0,dbl_mb(k_ei),1)    
          call daxpy(nxtrials,-omega(ivec),
     1      dbl_mb(k_vr+(ivec-1)*nxtrials),1,dbl_mb(k_ei),1)
          call daxpyfile(nxtrials,dbl_mb(k_ei),
     1      x4,d_r4,size_x4)
          residual(ivec) = residual(ivec) 
     1                   + ddotfile(d_r4,d_r4,size_x4)
          if (nodezero.and.util_print('residual',print_debug)) then
            write(LuOut,*) ' X4: vector',ivec
            write(LuOut,*) ' X4: residual',residual(ivec)
            call sf_print(d_r4,size_x4)
          endif
        endif
        residual(ivec) = dsqrt(residual(ivec))
        if (nodezero.and.util_print('residual',print_debug)) then
          write(LuOut,*) ivec,' residual',residual(ivec)
        endif
c
c       New trial vectors
c
        if (residual(ivec).gt.thresh) then
          if (newnxtrials.ge.maxtrials) call errquit
     1      ('tce_eom_xdiagon: nxtrials exceeded maxtrials',newnxtrials)
          if (needx1) call tce_jacobi_x1(d_r1,k_x1_offset)
          if (needx2) call tce_jacobi_x2(d_r2,k_x2_offset)
          if (needx3) call tce_jacobi_x3(d_r3,k_x3_offset)
          if (needx4) call tce_jacobi_x4(d_r4,k_x4_offset)
          if (needx1) call dscalfile(1.0d0/residual(ivec),d_r1,size_x1)
          if (needx2) call dscalfile(1.0d0/residual(ivec),d_r2,size_x2)
          if (needx3) call dscalfile(1.0d0/residual(ivec),d_r3,size_x3)
          if (needx4) call dscalfile(1.0d0/residual(ivec),d_r4,size_x4)
          do jvec = 1, newnxtrials
            overlap = 0.0d0
            if (needx1) overlap=overlap+ddotfile(d_r1,x1(jvec),size_x1)
            if (needx2) overlap=overlap+ddotfile(d_r2,x2(jvec),size_x2)
            if (needx3) overlap=overlap+ddotfile(d_r3,x3(jvec),size_x3)
            if (needx4) overlap=overlap+ddotfile(d_r4,x4(jvec),size_x4)
            if (needx1) call daxpyfile(1,-overlap,x1(jvec),d_r1,size_x1)
            if (needx2) call daxpyfile(1,-overlap,x2(jvec),d_r2,size_x2)
            if (needx3) call daxpyfile(1,-overlap,x3(jvec),d_r3,size_x3)
            if (needx4) call daxpyfile(1,-overlap,x4(jvec),d_r4,size_x4)
          enddo
          overlap = 0.0d0
          if (needx1) overlap=overlap+ddotfile(d_r1,d_r1,size_x1)
          if (needx2) overlap=overlap+ddotfile(d_r2,d_r2,size_x2)
          if (needx3) overlap=overlap+ddotfile(d_r3,d_r3,size_x3)
          if (needx4) overlap=overlap+ddotfile(d_r4,d_r4,size_x4)
          overlap = dsqrt(overlap)
          if (needx1) call dscalfile(1.0d0/overlap,d_r1,size_x1)
          if (needx2) call dscalfile(1.0d0/overlap,d_r2,size_x2)
          if (needx3) call dscalfile(1.0d0/overlap,d_r3,size_x3)
          if (needx4) call dscalfile(1.0d0/overlap,d_r4,size_x4)
          if (overlap.gt.thresh) then
            newnxtrials = newnxtrials + 1
            if (needx1) then
              call tce_filenameindexed(newnxtrials,'x1',filename)
              call createfile(filename,x1(newnxtrials),size_x1)
              call copyfile(d_r1,x1(newnxtrials),size_x1)
              x1_exist(newnxtrials) = .true.
            endif
            if (needx2) then
              call tce_filenameindexed(newnxtrials,'x2',filename)
              call createfile(filename,x2(newnxtrials),size_x2)
              call copyfile(d_r2,x2(newnxtrials),size_x2)
              x2_exist(newnxtrials) = .true.
            endif
            if (needx3) then
              call tce_filenameindexed(newnxtrials,'x3',filename)
              call createfile(filename,x3(newnxtrials),size_x3)
              call copyfile(d_r3,x3(newnxtrials),size_x3)
              x3_exist(newnxtrials) = .true.
            endif
            if (needx4) then
              call tce_filenameindexed(newnxtrials,'x4',filename)
              call createfile(filename,x4(newnxtrials),size_x4)
              call copyfile(d_r4,x4(newnxtrials),size_x4)
              x4_exist(newnxtrials) = .true.
            endif
          endif
        endif
      enddo
      if (newnxtrials.eq.nxtrials) then
        do ivec = 1,nroots_reduced
          if (residual(ivec).gt.thresh)
     1      call errquit('tce_eom_xdiagon: failed to extend subspace',0)
        enddo
      endif
      nxtrials = newnxtrials
c
c     Deallocate
c
      if (.not.ma_pop_stack(l_wk))
     1  call errquit('tce_eom_xdiagon: MA problem',6)
      if (.not.ma_pop_stack(l_ei))
     1  call errquit('tce_eom_xdiagon: MA problem',7)
      if (.not.ma_pop_stack(l_vl))
     1  call errquit('tce_eom_xdiagon: MA problem',9)
      if (.not.ma_pop_stack(l_vr))
     1  call errquit('tce_eom_xdiagon: MA problem',10)
      if (.not.ma_pop_stack(l_hbar))
     1  call errquit('tce_eom_xdiagon: MA problem',11)
 9000 format(1x,f20.10,' hartree ',f20.10,' eV')
      return
      end
