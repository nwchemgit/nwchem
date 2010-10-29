      subroutine cis_xdiagon(size_x1,k_x1_offset,d_r1,omega,residual)
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
      integer size_x1
      integer k_x1_offset
      integer d_r1
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
c      logical hbar_restore
c      external hbar_restore
c      double precision hbar
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
        if(.not.x1_exist(ivec))
     1    call errquit('cis_xdiagon: x1 not found',ivec,
     2      unknown_err)
        do jvec = 1, nxtrials
          cpu = - util_cpusec()
          wall = - util_wallsec()
          dbl_mb(k_hbar+(ivec-1)*nxtrials+jvec-1) = 0.0d0
          if (.not.xp1_exist(jvec))
     1      call errquit('ci_xdiagon: xp1 not found',jvec,
     2          UNKNOWN_ERR)
          dbl_mb(k_hbar+(ivec-1)*nxtrials+jvec-1) = 
     1    dbl_mb(k_hbar+(ivec-1)*nxtrials+jvec-1) + 
     2    ddotfile(x1(jvec),xp1(ivec),size_x1)
        enddo
      enddo
      if (nodezero.and.util_print('hbar',print_debug)) then
        write(LuOut,*)
        call ma_print(dbl_mb(k_hbar),nxtrials,nxtrials,
     1    "Similarity Transformed Hamiltonian")
      endif
c
c     Diagonalize the non-Hermitian matrix
c
      cpu = - util_cpusec()
      wall = - util_wallsec()
      call tce_dgeev('n','v',nxtrials,dbl_mb(k_hbar),nxtrials,
     1  omega,dbl_mb(k_ei),dbl_mb(k_vl),1,dbl_mb(k_vr),nxtrials,
     2  dbl_mb(k_wk),4*nxtrials,info)
      if (info.ne.0) call errquit('tce_eom_xdiagon: dgeev failed',info,
     1  UNKNOWN_ERR)
      call tce_sort(nxtrials,omega,dbl_mb(k_vr),'A')
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
        call tce_zero(d_r1,size_x1)
        call daxpyfile(nxtrials,dbl_mb(k_vr+(ivec-1)*nxtrials),
     1     xp1,d_r1,size_x1)
        call dfill(nxtrials,0.0d0,dbl_mb(k_ei),1)
        call daxpy(nxtrials,-omega(ivec),
     1      dbl_mb(k_vr+(ivec-1)*nxtrials),1,dbl_mb(k_ei),1)
        call daxpyfile(nxtrials,dbl_mb(k_ei),
     1     x1,d_r1,size_x1)
        residual(ivec) = residual(ivec) 
     1                 + ddotfile(d_r1,d_r1,size_x1)
        if (nodezero.and.util_print('residual',print_debug)) then
          write(LuOut,*) ' X1: vector',ivec
          write(LuOut,*) ' X1: residual',residual(ivec)
          call sf_print(d_r1,size_x1)
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
     1      ('cis_xdiagon: nxtrials exceeded maxtrials',newnxtrials,
     2      UNKNOWN_ERR)
          call tce_jacobi_x1(d_r1,k_x1_offset,0.0d0)
          call dscalfile(1.0d0/residual(ivec),d_r1,size_x1)
          do jvec = 1, newnxtrials
            overlap = 0.0d0
            overlap=overlap+ddotfile(d_r1,x1(jvec),size_x1)
            call daxpyfile(1,-overlap,x1(jvec),d_r1,size_x1)
          enddo
          overlap = 0.0d0
          overlap=overlap+ddotfile(d_r1,d_r1,size_x1)
          overlap = dsqrt(overlap)
          call dscalfile(1.0d0/overlap,d_r1,size_x1)
          if (overlap.gt.thresh) then
            newnxtrials = newnxtrials + 1
            call tce_filenameindexed(newnxtrials,'x1',filename)
            call createfile(filename,x1(newnxtrials),size_x1)
            call gatoeaf(x1(newnxtrials))
            call copyfile(d_r1,x1(newnxtrials),size_x1)
            call gatodra(x1(newnxtrials))
            x1_exist(newnxtrials) = .true.
          endif
        else
          call tce_zero(d_r1,size_x1)
          call daxpyfile(nxtrials,dbl_mb(k_vr+(ivec-1)*nxtrials),
     1        x1,d_r1,size_x1)
          call tce_filenameindexed(ivec,'x1_ini',filename)
          call createfile(filename,x1_ini(ivec),size_x1)
          call gatoeaf(x1_ini(ivec))
          call copyfile(d_r1,x1_ini(ivec),size_x1)
          call gatodra(x1_ini(ivec))
          x1_ini_exist(ivec) = .true.
c
c
c          call tce_filenameindexed(ivec,'x2_ini',filename)
c          call createfile(filename,x2_ini(ivec),size_x2)
c          call gatodra(x2_ini(ivec))
c          call gatoeaf(x2_ini(ivec))
c          x2_ini_exist(ivec) = .true.
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
     1      call errquit('cis_xdiagon: failed to extend subspace',0,
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
