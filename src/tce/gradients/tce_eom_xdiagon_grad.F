      subroutine tce_eom_xdiagon_grad(needx1,needx2,needx3,needx4,
     1  size_x1,size_x2,size_x3,size_x4,
     2  k_x1_offset,k_x2_offset,k_x3_offset,k_x4_offset,
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
      double precision cpu, wall
      double precision cpu1, wall1
      double precision cpu2, wall2
      logical hbar_restore
      external hbar_restore
      double precision hbar
c
      nodezero = (ga_nodeid().eq.0)
c
c     Allocate
c
c
      if (.not.ma_push_get(mt_dbl,nxtrials*nxtrials,'hbar',
     1  l_hbar,k_hbar)) call errquit('tce_eom_xdiagon: MA problem',0,
     2  MA_ERR)
      if (.not.ma_push_get(mt_dbl,nxtrials*nxtrials,'vr',
     1  l_vr,k_vr)) call errquit('tce_eom_xdiagon: MA problem',1,
     2  MA_ERR)
      if (.not.ma_push_get(mt_dbl,nxtrials,'vl',
     1  l_vl,k_vl)) call errquit('tce_eom_xdiagon: MA problem',2,
     2  MA_ERR)
      if (.not.ma_push_get(mt_dbl,nxtrials,'ei',
     1  l_ei,k_ei)) call errquit('tce_eom_xdiagon: MA problem',4,
     2  MA_ERR)
      if (.not.ma_push_get(mt_dbl,4*nxtrials,'work',
     1  l_wk,k_wk)) call errquit('tce_eom_xdiagon: MA problem',5,
     2  MA_ERR)
c
c     Form similarity transformed Hamiltonian subspace matrix
c
c     cpu = - util_cpusec()
c     wall = - util_wallsec()
      cpu1 = 0.0d0
      wall1 = 0.0d0
      cpu2 = 0.0d0
      wall2 = 0.0d0
      do ivec = 1, nxtrials
        do jvec = 1, nxtrials
c         if (hbar_restore(ivec,jvec,hbar)) then
c           cpu = - util_cpusec()
c           wall = - util_wallsec()
c           dbl_mb(k_hbar+(ivec-1)*nxtrials+jvec-1) = hbar
c           cpu1 = cpu1 + cpu + util_cpusec()
c           wall1 = wall1 + wall + util_wallsec()
c         else
            cpu = - util_cpusec()
            wall = - util_wallsec()
            dbl_mb(k_hbar+(ivec-1)*nxtrials+jvec-1) = 0.0d0
            if (needx1) then
              if (.not.x1_exist(ivec))
     1          call errquit('tce_eom_xdiagon: x1 not found',ivec,
     2          UNKNOWN_ERR)
              if (.not.xp1_exist(jvec))
     1          call errquit('tce_eom_xdiagon: xp1 not found',jvec,
     2          UNKNOWN_ERR)
              dbl_mb(k_hbar+(ivec-1)*nxtrials+jvec-1) = 
     1        dbl_mb(k_hbar+(ivec-1)*nxtrials+jvec-1) + 
     2        ddotfile(x1(jvec),xp1(ivec),size_x1)
            endif
            if (needx2) then
              if (.not.x2_exist(ivec))
     1          call errquit('tce_eom_xdiagon: x2 not found',ivec,
     2          UNKNOWN_ERR)
              if (.not.xp2_exist(jvec))
     1          call errquit('tce_eom_xdiagon: xp2 not found',jvec,
     2          UNKNOWN_ERR)
              dbl_mb(k_hbar+(ivec-1)*nxtrials+jvec-1) = 
     1        dbl_mb(k_hbar+(ivec-1)*nxtrials+jvec-1) + 
     2        ddotfile(x2(jvec),xp2(ivec),size_x2)
            endif
            if (needx3) then
              if (.not.x3_exist(ivec))
     1          call errquit('tce_eom_xdiagon: x3 not found',ivec,
     2          UNKNOWN_ERR)
              if (.not.xp3_exist(jvec))
     1          call errquit('tce_eom_xdiagon: xp3 not found',jvec,
     2          UNKNOWN_ERR)
              dbl_mb(k_hbar+(ivec-1)*nxtrials+jvec-1) = 
     1        dbl_mb(k_hbar+(ivec-1)*nxtrials+jvec-1) + 
     2        ddotfile(x3(jvec),xp3(ivec),size_x3)
            endif
            if (needx4) then
              if (.not.x4_exist(ivec))
     1          call errquit('tce_eom_xdiagon: x4 not found',ivec,
     2          UNKNOWN_ERR)
              if (.not.xp4_exist(jvec))
     1          call errquit('tce_eom_xdiagon: xp4 not found',jvec,
     2          UNKNOWN_ERR)
              dbl_mb(k_hbar+(ivec-1)*nxtrials+jvec-1) = 
     1        dbl_mb(k_hbar+(ivec-1)*nxtrials+jvec-1) + 
     2        ddotfile(x4(jvec),xp4(ivec),size_x4)
            endif
c           call hbar_store(ivec,jvec,
c    1        dbl_mb(k_hbar+(ivec-1)*nxtrials+jvec-1))
c           cpu2 = cpu2 + cpu + util_cpusec()
c           wall2 = wall2 + wall + util_wallsec()
c         endif
        enddo
      enddo
c     if (nodezero.and.util_print('eomtime',print_debug)) then
c       write(LuOut,9020) 'HBAR_RESTORE',cpu1, wall1
c       write(LuOut,9020) 'HBAR_STORE',cpu2, wall2
c       call util_flush(LuOut)
c     endif
c     cpu = cpu + util_cpusec()
c     wall = wall + util_wallsec()
c     if (nodezero.and.util_print('eomtime',print_debug)) then
c       write(LuOut,9020) 'HBAR',cpu, wall
c       call util_flush(LuOut)
c     endif
      if (nodezero.and.util_print('hbar',print_debug)) then
        write(LuOut,*)
        call ma_print(dbl_mb(k_hbar),nxtrials,nxtrials,
     1    "Similarity Transformed Hamiltonian")
      endif
c     call dcopy(nxtrials**2,dbl_mb(k_hbar),1,dbl_mb(k_gbar),1)
c
c     Diagonalize the non-Hermitian matrix
c
      cpu = - util_cpusec()
      wall = - util_wallsec()
      call util_dgeev('n','v',nxtrials,dbl_mb(k_hbar),nxtrials,
     1  omega,dbl_mb(k_ei),dbl_mb(k_vl),1,dbl_mb(k_vr),nxtrials,
     2  dbl_mb(k_wk),4*nxtrials,info)
      if (info.ne.0) call errquit('tce_eom_xdiagon: dgeev failed',info,
     1  UNKNOWN_ERR)
c     call tddft_schmidt(nxtrials,dbl_mb(k_vr))
      call util_nwsort(nxtrials,omega,dbl_mb(k_vr),'A')
      if (nodezero.and.util_print('omega',print_debug)) then
        do ivec = 1, nxtrials
          write(LuOut,9000) omega(ivec),omega(ivec)*au2ev
        enddo
      endif
      cpu = cpu + util_cpusec()
      wall = wall + util_wallsec()
      if (nodezero.and.util_print('eomtime',print_debug)) then
        write(LuOut,9020) 'DGEEV & SORT',cpu, wall
        call util_flush(LuOut)
      endif
      if (nodezero.and.util_print('vectors',print_debug)) then
        call ma_print(dbl_mb(k_vr),nxtrials,nxtrials,'Eigenvectors')
      endif
c
c     Compute residuals
c
      newnxtrials = nxtrials
      do ivec = 1, nroots_reduced
        if (nodezero.and.util_print('eomtime',print_debug)) then
          write(LuOut,9030) ivec
          call util_flush(LuOut)
        endif
        cpu = - util_cpusec()
        wall = - util_wallsec()
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
        cpu = cpu + util_cpusec()
        wall = wall + util_wallsec()
        if (nodezero.and.util_print('eomtime',print_debug)) then
          write(LuOut,9020) 'RESIDUAL',cpu, wall
          call util_flush(LuOut)
        endif
c
c       New trial vectors
c
        cpu = - util_cpusec()
        wall = - util_wallsec()
        if (residual(ivec).gt.thresh) then
          if (newnxtrials.ge.maxtrials) call errquit
     1      ('tce_eom_xdiagon: nxtrials exceeded maxtrials',newnxtrials,
     2      UNKNOWN_ERR)
          if (needx1) call tce_jacobi_x1(d_r1,k_x1_offset)
          if (needx2) call tce_jacobi_x2(d_r2,k_x2_offset)
          if (needx3) then
           if(activecalc) then 
            call tce_jacobi_x3a(d_r3,k_x3_offset)
           else 
            call tce_jacobi_x3(d_r3,k_x3_offset)
           end if
          end if
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
              call gatoeaf(x1(newnxtrials))
              call copyfile(d_r1,x1(newnxtrials),size_x1)
              call gatodra(x1(newnxtrials))
              x1_exist(newnxtrials) = .true.
            endif
            if (needx2) then
              call tce_filenameindexed(newnxtrials,'x2',filename)
              call createfile(filename,x2(newnxtrials),size_x2)
              call gatoeaf(x2(newnxtrials))
              call copyfile(d_r2,x2(newnxtrials),size_x2)
              call gatodra(x2(newnxtrials))
              x2_exist(newnxtrials) = .true.
            endif
            if (needx3) then
              call tce_filenameindexed(newnxtrials,'x3',filename)
              call createfile(filename,x3(newnxtrials),size_x3)
              call gatoeaf(x3(newnxtrials))
              call copyfile(d_r3,x3(newnxtrials),size_x3)
              call gatodra(x3(newnxtrials))
              x3_exist(newnxtrials) = .true.
            endif
            if (needx4) then
              call tce_filenameindexed(newnxtrials,'x4',filename)
              call createfile(filename,x4(newnxtrials),size_x4)
              call gatoeaf(x4(newnxtrials))
              call copyfile(d_r4,x4(newnxtrials),size_x4)
              call gatodra(x4(newnxtrials))
              x4_exist(newnxtrials) = .true.
            endif
          endif
        else
c         if (nodezero) write(LuOut,9010) ivec
          if (needx1) then
            call tce_zero(d_r1,size_x1)
            call daxpyfile(nxtrials,dbl_mb(k_vr+(ivec-1)*nxtrials),
     1        x1,d_r1,size_x1)
            if (.not.xc1_exist(ivec)) then
               call tce_filenameindexed(ivec,'xc1',filename)
               call createfile(filename,xc1(ivec),size_x1)
            endif
            call gatoeaf(xc1(ivec))
            call copyfile(d_r1,xc1(ivec),size_x1)
            call gatodra(xc1(ivec))
            xc1_exist(ivec) = .true.
          endif
          if (needx2) then
            call tce_zero(d_r2,size_x2)
            call daxpyfile(nxtrials,dbl_mb(k_vr+(ivec-1)*nxtrials),
     1        x2,d_r2,size_x2)
            if (.not.xc2_exist(ivec)) then
               call tce_filenameindexed(ivec,'xc2',filename)
               call createfile(filename,xc2(ivec),size_x2)
            endif
            call gatoeaf(xc2(ivec))
            call copyfile(d_r2,xc2(ivec),size_x2)
            call gatodra(xc2(ivec))
            xc2_exist(ivec) = .true.
          endif
          if (needx3) then
            call tce_zero(d_r3,size_x3)
            call daxpyfile(nxtrials,dbl_mb(k_vr+(ivec-1)*nxtrials),
     1        x3,d_r3,size_x3)
            if (.not.xc3_exist(ivec)) then
               call tce_filenameindexed(ivec,'xc3',filename)
               call createfile(filename,xc3(ivec),size_x3)
            endif
            call gatoeaf(xc3(ivec))
            call copyfile(d_r3,xc3(ivec),size_x3)
            call gatodra(xc3(ivec))
            xc3_exist(ivec) = .true.
          endif
          if (needx4) then
            call tce_zero(d_r4,size_x4)
            call daxpyfile(nxtrials,dbl_mb(k_vr+(ivec-1)*nxtrials),
     1        x4,d_r4,size_x4)
            if (.not.xc4_exist(ivec)) then
               call tce_filenameindexed(ivec,'xc4',filename)
               call createfile(filename,xc4(ivec),size_x4)
            endif
            call gatoeaf(xc4(ivec))
            call copyfile(d_r4,xc4(ivec),size_x4)
            call gatodra(xc4(ivec))
            xc4_exist(ivec) = .true.
          endif
        endif
        cpu = cpu + util_cpusec()
        wall = wall + util_wallsec()
        if (nodezero.and.util_print('eomtime',print_debug)) then
          write(LuOut,9020) 'NEW TRIAL',cpu, wall
          call util_flush(LuOut)
        endif
      enddo
      if (newnxtrials.eq.nxtrials) then
        do ivec = 1,nroots_reduced
          if (residual(ivec).gt.thresh)
     1      call errquit('tce_eom_xdiagon: failed to extend subspace',0,
     2      CALC_ERR)
        enddo
      endif
      nxtrials = newnxtrials
c
c     Deallocate
c
      if (.not.ma_pop_stack(l_wk))
     1  call errquit('tce_eom_xdiagon: MA problem',6,MA_ERR)
      if (.not.ma_pop_stack(l_ei))
     1  call errquit('tce_eom_xdiagon: MA problem',7,MA_ERR)
      if (.not.ma_pop_stack(l_vl))
     1  call errquit('tce_eom_xdiagon: MA problem',9,MA_ERR)
      if (.not.ma_pop_stack(l_vr))
     1  call errquit('tce_eom_xdiagon: MA problem',10,MA_ERR)
      if (.not.ma_pop_stack(l_hbar))
     1  call errquit('tce_eom_xdiagon: MA problem',11,MA_ERR)
 9000 format(1x,f20.10,' hartree ',f20.10,' eV')
 9010 format(1x,'Root ',i3,' converged')
 9020 format(1x,A,' cpu:',f10.2,' wall:',f10.2)
 9030 format(1x,'IVEC = ',i3)
      return
      end
