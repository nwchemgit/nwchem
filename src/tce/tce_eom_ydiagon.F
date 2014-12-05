      subroutine tce_eom_ydiagon(needy1,needy2,needy3,needy4,
     1  size_y1,size_y2,size_y3,size_y4,
     2  k_y1_offset,k_y2_offset,k_y3_offset,k_y4_offset,
     3  d_r1,d_r2,d_r3,d_r4,omega,residual)
c
c $Id$
c
c     Form a subspace matrix and diagonalize
c
      implicit none
#include "stdio.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "tce.fh"
#include "tce_main.fh"
#include "tce_diis.fh"
#include "errquit.fh"
#include "util.fh"
      logical needy1
      logical needy2
      logical needy3
      logical needy4
      integer size_y1
      integer size_y2
      integer size_y3
      integer size_y4
      integer k_y1_offset
      integer k_y2_offset
      integer k_y3_offset
      integer k_y4_offset
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
      integer newnytrials
      double precision ddotfile,ddotfile_1,ddotfile_2
      external ddotfile,ddotfile_1,ddotfile_2
      double precision au2ev   ! Conversion factor from a.u. to eV
      parameter (au2ev=27.2113961d0)
      character*255 filename
c     logical hbar_restore
c     external hbar_restore
c     double precision hbar
c
      nodezero = (ga_nodeid().eq.0)
c
c     Allocate
c
      if (.not.ma_push_get(mt_dbl,nytrials*nytrials,'hbar',
     1  l_hbar,k_hbar)) call errquit('tce_eom_ydiagon: MA problem',0,
     2  UNKNOWN_ERR)
      if (.not.ma_push_get(mt_dbl,nytrials*nytrials,'vr',
     1  l_vr,k_vr)) call errquit('tce_eom_ydiagon: MA problem',1,
     2  UNKNOWN_ERR)
      if (.not.ma_push_get(mt_dbl,nytrials,'vl',
     1  l_vl,k_vl)) call errquit('tce_eom_ydiagon: MA problem',2,
     2  UNKNOWN_ERR)
      if (.not.ma_push_get(mt_dbl,nytrials,'ei',
     1  l_ei,k_ei)) call errquit('tce_eom_ydiagon: MA problem',4,
     2  UNKNOWN_ERR)
      if (.not.ma_push_get(mt_dbl,4*nytrials,'work',
     1  l_wk,k_wk)) call errquit('tce_eom_ydiagon: MA problem',5,
     2  UNKNOWN_ERR)
c
c     Form similarity transformed Hamiltonian subspace matrix
c
      do ivec = 1, nytrials
        do jvec = 1, nytrials
c         if (hbar_restore(ivec,jvec,hbar)) then
c           dbl_mb(k_hbar+(ivec-1)*nytrials+jvec-1) = hbar
c         else
            dbl_mb(k_hbar+(ivec-1)*nytrials+jvec-1) = 0.0d0
            if (needy1) then
              if (.not.y1_exist(ivec))
     1          call errquit('tce_eom_ydiagon: y1 not found',ivec,
     2          UNKNOWN_ERR)
              if (.not.yp1_exist(jvec))
     1          call errquit('tce_eom_ydiagon: yp1 not found',jvec,
     2          UNKNOWN_ERR)
              dbl_mb(k_hbar+(ivec-1)*nytrials+jvec-1) = 
     1        dbl_mb(k_hbar+(ivec-1)*nytrials+jvec-1) + 
     2        ddotfile_1(y1(jvec),yp1(ivec),size_y1)
            endif
            if (needy2) then
              if (.not.y2_exist(ivec))
     1          call errquit('tce_eom_ydiagon: y2 not found',ivec,
     2          UNKNOWN_ERR)
              if (.not.yp2_exist(jvec))
     1          call errquit('tce_eom_ydiagon: yp2 not found',jvec,
     2          UNKNOWN_ERR)
              dbl_mb(k_hbar+(ivec-1)*nytrials+jvec-1) = 
     1        dbl_mb(k_hbar+(ivec-1)*nytrials+jvec-1) + 
     2        ddotfile_2(y2(jvec),yp2(ivec),size_y2)
            endif
            if (needy3) then
              if (.not.y3_exist(ivec))
     1          call errquit('tce_eom_ydiagon: y3 not found',ivec,
     2          UNKNOWN_ERR)
              if (.not.yp3_exist(jvec))
     1          call errquit('tce_eom_ydiagon: yp3 not found',jvec,
     2          UNKNOWN_ERR)
              dbl_mb(k_hbar+(ivec-1)*nytrials+jvec-1) = 
     1        dbl_mb(k_hbar+(ivec-1)*nytrials+jvec-1) + 
     2        ddotfile(y3(jvec),yp3(ivec),size_y3)
            endif
            if (needy4) then
              if (.not.y4_exist(ivec))
     1          call errquit('tce_eom_ydiagon: y4 not found',ivec,
     2          UNKNOWN_ERR)
              if (.not.yp4_exist(jvec))
     1          call errquit('tce_eom_ydiagon: yp4 not found',jvec,
     2          UNKNOWN_ERR)
              dbl_mb(k_hbar+(ivec-1)*nytrials+jvec-1) = 
     1        dbl_mb(k_hbar+(ivec-1)*nytrials+jvec-1) + 
     2        ddotfile(y4(jvec),yp4(ivec),size_y4)
            endif
c           call hbar_store(ivec,jvec,
c    1        dbl_mb(k_hbar+(ivec-1)*nytrials+jvec-1))
c         endif
        enddo
      enddo
      if (nodezero.and.util_print('hbar',print_debug)) then
        write(LuOut,*)
        call ma_print(dbl_mb(k_hbar),nytrials,nytrials,
     1    "Similarity Transformed Hamiltonian")
      endif
c     call dcopy(nytrials**2,dbl_mb(k_hbar),1,dbl_mb(k_gbar),1)
c
c     Diagonalize the non-Hermitian matrix
c
cjeff      call tce_dgeev('n','v',nytrials,dbl_mb(k_hbar),nytrials,
      call dgeev('n','v',nytrials,dbl_mb(k_hbar),nytrials,
     1  omega,dbl_mb(k_ei),dbl_mb(k_vl),1,dbl_mb(k_vr),nytrials,
     2  dbl_mb(k_wk),4*nytrials,info)
      if (info.ne.0) call errquit('tce_eom_ydiagon: dgeev failed',info,
     1  UNKNOWN_ERR)
c     call tce_schmidt(nytrials,dbl_mb(k_vr))
      call tce_sort(nytrials,omega,dbl_mb(k_vr),'A')
      if (nodezero.and.util_print('omega',print_debug)) then
        do ivec = 1, nytrials
          write(LuOut,9000) omega(ivec),omega(ivec)*au2ev
        enddo
      endif
      if (nodezero.and.util_print('vectors',print_debug)) then
        call ma_print(dbl_mb(k_vr),nytrials,nytrials,'Eigenvectors')
      endif
c
c     Compute residuals
c
      newnytrials = nytrials
      do ivec = 1, nroots_reduced
        residual(ivec) = 0.0d0
        if (needy1) then
          call tce_zero(d_r1,size_y1)
          call daxpyfile(nytrials,dbl_mb(k_vr+(ivec-1)*nytrials),
     1      yp1,d_r1,size_y1)
          call dfill(nytrials,0.0d0,dbl_mb(k_ei),1)    
          call daxpy(nytrials,-omega(ivec),
     1      dbl_mb(k_vr+(ivec-1)*nytrials),1,dbl_mb(k_ei),1)
          call daxpyfile(nytrials,dbl_mb(k_ei),
     1      y1,d_r1,size_y1)
          residual(ivec) = residual(ivec) 
     1                   + ddotfile_1(d_r1,d_r1,size_y1)
          if (nodezero.and.util_print('residual',print_debug)) then
            write(LuOut,*) ' Y1: vector',ivec
            write(LuOut,*) ' Y1: residual',residual(ivec)
            call sf_print(d_r1,size_y1)
          endif
        endif
        if (needy2) then
          call tce_zero(d_r2,size_y2)
          call daxpyfile(nytrials,dbl_mb(k_vr+(ivec-1)*nytrials),
     1      yp2,d_r2,size_y2)
          call dfill(nytrials,0.0d0,dbl_mb(k_ei),1)    
          call daxpy(nytrials,-omega(ivec),
     1      dbl_mb(k_vr+(ivec-1)*nytrials),1,dbl_mb(k_ei),1)
          call daxpyfile(nytrials,dbl_mb(k_ei),
     1      y2,d_r2,size_y2)
          residual(ivec) = residual(ivec) 
     1                   + ddotfile_2(d_r2,d_r2,size_y2)
          if (nodezero.and.util_print('residual',print_debug)) then
            write(LuOut,*) ' Y2: vector',ivec
            write(LuOut,*) ' Y2: residual',residual(ivec)
            call sf_print(d_r2,size_y2)
          endif
        endif
        if (needy3) then
          call tce_zero(d_r3,size_y3)
          call daxpyfile(nytrials,dbl_mb(k_vr+(ivec-1)*nytrials),
     1      yp3,d_r3,size_y3)
          call dfill(nytrials,0.0d0,dbl_mb(k_ei),1)    
          call daxpy(nytrials,-omega(ivec),
     1      dbl_mb(k_vr+(ivec-1)*nytrials),1,dbl_mb(k_ei),1)
          call daxpyfile(nytrials,dbl_mb(k_ei),
     1      y3,d_r3,size_y3)
          residual(ivec) = residual(ivec) 
     1                   + ddotfile(d_r3,d_r3,size_y3)
          if (nodezero.and.util_print('residual',print_debug)) then
            write(LuOut,*) ' Y3: vector',ivec
            write(LuOut,*) ' Y3: residual',residual(ivec)
            call sf_print(d_r3,size_y3)
          endif
        endif
        if (needy4) then
          call tce_zero(d_r4,size_y4)
          call daxpyfile(nytrials,dbl_mb(k_vr+(ivec-1)*nytrials),
     1      yp4,d_r4,size_y4)
          call dfill(nytrials,0.0d0,dbl_mb(k_ei),1)    
          call daxpy(nytrials,-omega(ivec),
     1      dbl_mb(k_vr+(ivec-1)*nytrials),1,dbl_mb(k_ei),1)
          call daxpyfile(nytrials,dbl_mb(k_ei),
     1      y4,d_r4,size_y4)
          residual(ivec) = residual(ivec) 
     1                   + ddotfile(d_r4,d_r4,size_y4)
          if (nodezero.and.util_print('residual',print_debug)) then
            write(LuOut,*) ' Y4: vector',ivec
            write(LuOut,*) ' Y4: residual',residual(ivec)
            call sf_print(d_r4,size_y4)
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
          if (newnytrials.ge.maxtrials) call errquit
     1      ('tce_eom_ydiagon: nytrials exceeded maxtrials',newnytrials,
     2      CALC_ERR)
          if (needy1) call tce_jacobi_y1(d_r1,k_y1_offset)
          if (needy2) call tce_jacobi_y2(d_r2,k_y2_offset)
          if (needy3) call tce_jacobi_y3(d_r3,k_y3_offset)
          if (needy4) call tce_jacobi_y4(d_r4,k_y4_offset)
          if (needy1) call dscalfile(1.0d0/residual(ivec),d_r1,size_y1)
          if (needy2) call dscalfile(1.0d0/residual(ivec),d_r2,size_y2)
          if (needy3) call dscalfile(1.0d0/residual(ivec),d_r3,size_y3)
          if (needy4) call dscalfile(1.0d0/residual(ivec),d_r4,size_y4)
          do jvec = 1, newnytrials
            overlap = 0.0d0
            if(needy1) overlap=overlap+ddotfile_1(d_r1,y1(jvec),size_y1)
            if(needy2) overlap=overlap+ddotfile_2(d_r2,y2(jvec),size_y2)
            if (needy3) overlap=overlap+ddotfile(d_r3,y3(jvec),size_y3)
            if (needy4) overlap=overlap+ddotfile(d_r4,y4(jvec),size_y4)
            if (needy1) call daxpyfile(1,-overlap,y1(jvec),d_r1,size_y1)
            if (needy2) call daxpyfile(1,-overlap,y2(jvec),d_r2,size_y2)
            if (needy3) call daxpyfile(1,-overlap,y3(jvec),d_r3,size_y3)
            if (needy4) call daxpyfile(1,-overlap,y4(jvec),d_r4,size_y4)
          enddo
          overlap = 0.0d0
          if (needy1) overlap=overlap+ddotfile_1(d_r1,d_r1,size_y1)
          if (needy2) overlap=overlap+ddotfile_2(d_r2,d_r2,size_y2)
          if (needy3) overlap=overlap+ddotfile(d_r3,d_r3,size_y3)
          if (needy4) overlap=overlap+ddotfile(d_r4,d_r4,size_y4)
          overlap = dsqrt(overlap)
          if (needy1) call dscalfile(1.0d0/overlap,d_r1,size_y1)
          if (needy2) call dscalfile(1.0d0/overlap,d_r2,size_y2)
          if (needy3) call dscalfile(1.0d0/overlap,d_r3,size_y3)
          if (needy4) call dscalfile(1.0d0/overlap,d_r4,size_y4)
          if (overlap.gt.thresh) then
            newnytrials = newnytrials + 1
            if (needy1) then
              call tce_filenameindexed(newnytrials,'y1',filename)
              call createfile(filename,y1(newnytrials),size_y1)
              call gatoeaf(y1(newnytrials))
              call copyfile(d_r1,y1(newnytrials),size_y1)
              call gatodra(y1(newnytrials))
              y1_exist(newnytrials) = .true.
            endif
            if (needy2) then
              call tce_filenameindexed(newnytrials,'y2',filename)
              call createfile(filename,y2(newnytrials),size_y2)
              call gatoeaf(y2(newnytrials))
              call copyfile(d_r2,y2(newnytrials),size_y2)
              call gatodra(y2(newnytrials))
              y2_exist(newnytrials) = .true.
            endif
            if (needy3) then
              call tce_filenameindexed(newnytrials,'y3',filename)
              call createfile(filename,y3(newnytrials),size_y3)
              call gatoeaf(y3(newnytrials))
              call copyfile(d_r3,y3(newnytrials),size_y3)
              call gatodra(y3(newnytrials))
              y3_exist(newnytrials) = .true.
            endif
            if (needy4) then
              call tce_filenameindexed(newnytrials,'y4',filename)
              call createfile(filename,y4(newnytrials),size_y4)
              call gatoeaf(y4(newnytrials))
              call copyfile(d_r4,y4(newnytrials),size_y4)
              call gatodra(y4(newnytrials))
              y4_exist(newnytrials) = .true.
            endif
          endif
        else
c         if (nodezero) write(LuOut,9010) ivec
          if (needy1) then
            call tce_zero(d_r1,size_y1)
            call daxpyfile(nytrials,dbl_mb(k_vr+(ivec-1)*nytrials),
     1        y1,d_r1,size_y1)
            call tce_filenameindexed(ivec,'yc1',filename)
            call createfile(filename,yc1(ivec),size_y1)
            call gatoeaf(yc1(ivec))
            call copyfile(d_r1,yc1(ivec),size_y1)
            call gatodra(yc1(ivec))
            yc1_exist(ivec) = .true.
          endif
          if (needy2) then
            call tce_zero(d_r2,size_y2)
            call daxpyfile(nytrials,dbl_mb(k_vr+(ivec-1)*nytrials),
     1        y2,d_r2,size_y2)
            call tce_filenameindexed(ivec,'yc2',filename)
            call createfile(filename,yc2(ivec),size_y2)
            call gatoeaf(yc2(ivec))
            call copyfile(d_r2,yc2(ivec),size_y2)
            call gatodra(yc2(ivec))
            yc2_exist(ivec) = .true.
          endif
          if (needy3) then
            call tce_zero(d_r3,size_y3)
            call daxpyfile(nytrials,dbl_mb(k_vr+(ivec-1)*nytrials),
     1        y3,d_r3,size_y3)
            call tce_filenameindexed(ivec,'yc3',filename)
            call createfile(filename,yc3(ivec),size_y3)
            call gatoeaf(yc3(ivec))
            call copyfile(d_r3,yc3(ivec),size_y3)
            call gatodra(yc3(ivec))
            yc3_exist(ivec) = .true.
          endif
          if (needy4) then
            call tce_zero(d_r4,size_y4)
            call daxpyfile(nytrials,dbl_mb(k_vr+(ivec-1)*nytrials),
     1        y4,d_r4,size_y4)
            call tce_filenameindexed(ivec,'yc4',filename)
            call createfile(filename,yc4(ivec),size_y4)
            call gatoeaf(yc4(ivec))
            call copyfile(d_r4,yc4(ivec),size_y4)
            call gatodra(yc4(ivec))
            yc4_exist(ivec) = .true.
          endif
        endif
      enddo
      if (newnytrials.eq.nytrials) then
        do ivec = 1,nroots_reduced
          if (residual(ivec).gt.thresh)
     1      call errquit('tce_eom_ydiagon: failed to extend subspace',0,
     2      CALC_ERR)
        enddo
      endif
      nytrials = newnytrials
c
c     Deallocate
c
      if (.not.ma_pop_stack(l_wk))
     1  call errquit('tce_eom_ydiagon: MA problem',6,MA_ERR)
      if (.not.ma_pop_stack(l_ei))
     1  call errquit('tce_eom_ydiagon: MA problem',7,MA_ERR)
      if (.not.ma_pop_stack(l_vl))
     1  call errquit('tce_eom_ydiagon: MA problem',9,MA_ERR)
      if (.not.ma_pop_stack(l_vr))
     1  call errquit('tce_eom_ydiagon: MA problem',10,MA_ERR)
      if (.not.ma_pop_stack(l_hbar))
     1  call errquit('tce_eom_ydiagon: MA problem',11,MA_ERR)
 9000 format(1x,f20.10,' hartree ',f20.10,' eV')
 9010 format(1x,'Root ',i3,' converged')
      return
      end
