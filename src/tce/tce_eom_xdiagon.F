      subroutine tce_eom_xdiagon(needx1,needx2,needx3,needx4,
     1  size_x1,size_x2,size_x3,size_x4,
     2  k_x1_offset,k_x2_offset,k_x3_offset,k_x4_offset,
     3  d_r1,d_r2,d_r3,d_r4,omega,residual,k_hbar,iteration,
     +  eaccsd,ipccsd)
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
      integer l_shbar,k_shbar
      integer l_vr, k_vr
      integer l_vl, k_vl
      integer l_ei, k_ei
      integer l_wk, k_wk
      integer info
      integer newnxtrials
      double precision ddotfile,ddotfile_1,ddotfile_2
      external ddotfile,ddotfile_1,ddotfile_2
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
      integer iteration
ckbn --- EA/IPCCSD
c#ifdef EACCSD
      logical eaccsd
c#endif
c#ifdef IPCCSD
      logical ipccsd
c#endif
c
      nodezero = (ga_nodeid().eq.0)
c
c     Allocate
c
c Little change here - little matrxi l_shbar is introduced
      if (.not.ma_push_get(mt_dbl,nxtrials*nxtrials,'hbar',
     1  l_shbar,k_shbar)) call errquit('tce_eom_xdiagon: MA problem',0,
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
c --------NEW SOLVER----------------------------------
        if(nxtrials.gt.hbard) 
     &    call errquit('tce_eom_xdiagon: nxtrails too big',4,MA_ERR)
c ----------------------------------------------------
c
c     Form similarity transformed Hamiltonian subspace matrix
c
c     cpu = - util_cpusec()
c     wall = - util_wallsec()
      cpu1 = 0.0d0
      wall1 = 0.0d0
      cpu2 = 0.0d0
      wall2 = 0.0d0
c
c declare smaller matrix nxtrials*nxtrials and copy 
c whatever need to be copied from k_hbar
c
      if(iteration.eq.1) then
         do ivec = 1, nxtrials
         do jvec = 1, nxtrials
            cpu = - util_cpusec()
            wall = - util_wallsec()
            dbl_mb(k_hbar+(ivec-1)*hbard+jvec-1) = 0.0d0
            if (needx1) then
              if (.not.x1_exist(ivec))
     1          call errquit('tce_eom_xdiagon: x1 not found',ivec,
     2          UNKNOWN_ERR)
              if (.not.xp1_exist(jvec))
     1          call errquit('tce_eom_xdiagon: xp1 not found',jvec,
     2          UNKNOWN_ERR)
              dbl_mb(k_hbar+(ivec-1)*hbard+jvec-1) =
     1        dbl_mb(k_hbar+(ivec-1)*hbard+jvec-1) +
     2        ddotfile_1(x1(jvec),xp1(ivec),size_x1)
            endif
            if (needx2) then
              if (.not.x2_exist(ivec))
     1          call errquit('tce_eom_xdiagon: x2 not found',ivec,
     2          UNKNOWN_ERR)
              if (.not.xp2_exist(jvec))
     1          call errquit('tce_eom_xdiagon: xp2 not found',jvec,
     2          UNKNOWN_ERR)
              dbl_mb(k_hbar+(ivec-1)*hbard+jvec-1) =
     1        dbl_mb(k_hbar+(ivec-1)*hbard+jvec-1) +
     2        ddotfile_2(x2(jvec),xp2(ivec),size_x2)
            endif
            if (needx3) then
              if (.not.x3_exist(ivec))
     1          call errquit('tce_eom_xdiagon: x3 not found',ivec,
     2          UNKNOWN_ERR)
              if (.not.xp3_exist(jvec))
     1          call errquit('tce_eom_xdiagon: xp3 not found',jvec,
     2          UNKNOWN_ERR)
              dbl_mb(k_hbar+(ivec-1)*hbard+jvec-1) =
     1        dbl_mb(k_hbar+(ivec-1)*hbard+jvec-1) +
     2        ddotfile(x3(jvec),xp3(ivec),size_x3)
            endif
            if (needx4) then
              if (.not.x4_exist(ivec))
     1          call errquit('tce_eom_xdiagon: x4 not found',ivec,
     2          UNKNOWN_ERR)
              if (.not.xp4_exist(jvec))
     1          call errquit('tce_eom_xdiagon: xp4 not found',jvec,
     2          UNKNOWN_ERR)
              dbl_mb(k_hbar+(ivec-1)*hbard+jvec-1) =
     1        dbl_mb(k_hbar+(ivec-1)*hbard+jvec-1) +
     2        ddotfile(x4(jvec),xp4(ivec),size_x4)
            endif
         enddo
         enddo
      else
c
         do ivec = nxtrials-nroots_reduced+1,nxtrials !-- ivec column
         do jvec = 1, nxtrials !--------------- jvec row
            cpu = - util_cpusec()
            wall = - util_wallsec()
            dbl_mb(k_hbar+(ivec-1)*hbard+jvec-1) = 0.0d0
            if (needx1) then
              if (.not.x1_exist(ivec))
     1          call errquit('tce_eom_xdiagon: x1 not found',ivec,
     2          UNKNOWN_ERR)
              if (.not.xp1_exist(jvec))
     1          call errquit('tce_eom_xdiagon: xp1 not found',jvec,
     2          UNKNOWN_ERR)
              dbl_mb(k_hbar+(ivec-1)*hbard+jvec-1) = 
     1        dbl_mb(k_hbar+(ivec-1)*hbard+jvec-1) + 
     2        ddotfile_1(x1(jvec),xp1(ivec),size_x1)
            endif
            if (needx2) then
              if (.not.x2_exist(ivec))
     1          call errquit('tce_eom_xdiagon: x2 not found',ivec,
     2          UNKNOWN_ERR)
              if (.not.xp2_exist(jvec))
     1          call errquit('tce_eom_xdiagon: xp2 not found',jvec,
     2          UNKNOWN_ERR)
              dbl_mb(k_hbar+(ivec-1)*hbard+jvec-1) = 
     1        dbl_mb(k_hbar+(ivec-1)*hbard+jvec-1) + 
     2        ddotfile_2(x2(jvec),xp2(ivec),size_x2)
            endif
            if (needx3) then
              if (.not.x3_exist(ivec))
     1          call errquit('tce_eom_xdiagon: x3 not found',ivec,
     2          UNKNOWN_ERR)
              if (.not.xp3_exist(jvec))
     1          call errquit('tce_eom_xdiagon: xp3 not found',jvec,
     2          UNKNOWN_ERR)
              dbl_mb(k_hbar+(ivec-1)*hbard+jvec-1) = 
     1        dbl_mb(k_hbar+(ivec-1)*hbard+jvec-1) + 
     2        ddotfile(x3(jvec),xp3(ivec),size_x3)
            endif
            if (needx4) then
              if (.not.x4_exist(ivec))
     1          call errquit('tce_eom_xdiagon: x4 not found',ivec,
     2          UNKNOWN_ERR)
              if (.not.xp4_exist(jvec))
     1          call errquit('tce_eom_xdiagon: xp4 not found',jvec,
     2          UNKNOWN_ERR)
              dbl_mb(k_hbar+(ivec-1)*hbard+jvec-1) = 
     1        dbl_mb(k_hbar+(ivec-1)*hbard+jvec-1) + 
     2        ddotfile(x4(jvec),xp4(ivec),size_x4)
            endif
         enddo ! -------------------------- jvec - row
         enddo ! -------------------------- ivec - column
c
         do ivec = 1, nxtrials-nroots_reduced !---------- ivec column
         do jvec = nxtrials-nroots_reduced+1,nxtrials
            cpu = - util_cpusec()
            wall = - util_wallsec()
            dbl_mb(k_hbar+(ivec-1)*hbard+jvec-1) = 0.0d0
            if (needx1) then
              if (.not.x1_exist(ivec))
     1          call errquit('tce_eom_xdiagon: x1 not found',ivec,
     2          UNKNOWN_ERR)
              if (.not.xp1_exist(jvec))
     1          call errquit('tce_eom_xdiagon: xp1 not found',jvec,
     2          UNKNOWN_ERR)
              dbl_mb(k_hbar+(ivec-1)*hbard+jvec-1) = 
     1        dbl_mb(k_hbar+(ivec-1)*hbard+jvec-1) + 
     2        ddotfile_1(x1(jvec),xp1(ivec),size_x1)
            endif
            if (needx2) then
              if (.not.x2_exist(ivec))
     1          call errquit('tce_eom_xdiagon: x2 not found',ivec,
     2          UNKNOWN_ERR)
              if (.not.xp2_exist(jvec))
     1          call errquit('tce_eom_xdiagon: xp2 not found',jvec,
     2          UNKNOWN_ERR)
              dbl_mb(k_hbar+(ivec-1)*hbard+jvec-1) = 
     1        dbl_mb(k_hbar+(ivec-1)*hbard+jvec-1) + 
     2        ddotfile_2(x2(jvec),xp2(ivec),size_x2)
            endif
            if (needx3) then
              if (.not.x3_exist(ivec))
     1          call errquit('tce_eom_xdiagon: x3 not found',ivec,
     2          UNKNOWN_ERR)
              if (.not.xp3_exist(jvec))
     1          call errquit('tce_eom_xdiagon: xp3 not found',jvec,
     2          UNKNOWN_ERR)
              dbl_mb(k_hbar+(ivec-1)*hbard+jvec-1) = 
     1        dbl_mb(k_hbar+(ivec-1)*hbard+jvec-1) + 
     2        ddotfile(x3(jvec),xp3(ivec),size_x3)
            endif
            if (needx4) then
              if (.not.x4_exist(ivec))
     1          call errquit('tce_eom_xdiagon: x4 not found',ivec,
     2          UNKNOWN_ERR)
              if (.not.xp4_exist(jvec))
     1          call errquit('tce_eom_xdiagon: xp4 not found',jvec,
     2          UNKNOWN_ERR)
              dbl_mb(k_hbar+(ivec-1)*hbard+jvec-1) = 
     1        dbl_mb(k_hbar+(ivec-1)*hbard+jvec-1) + 
     2        ddotfile(x4(jvec),xp4(ivec),size_x4)
            endif
         enddo ! -------------------------- jvec - row
         enddo ! -------------------------- ivec - column
      endif
c
c
c
c  copy from l_hbar to l_shbar
c
        do ivec=1,nxtrials !row
        do jvec=1,nxtrials !column
           dbl_mb(k_shbar+(ivec-1)*nxtrials+jvec-1)=
     &     dbl_mb(k_hbar+(ivec-1)*hbard+jvec-1)
        enddo
        enddo
c
c
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
c
c
c  DIAGONALIZE l_shbar
c
c
      call util_dgeev('n','v',nxtrials,dbl_mb(k_shbar),nxtrials,
     1  omega,dbl_mb(k_ei),dbl_mb(k_vl),1,dbl_mb(k_vr),nxtrials,
     2  dbl_mb(k_wk),4*nxtrials,info)
      if (info.ne.0) call errquit('tce_eom_xdiagon: dgeev failed',info,
     1  UNKNOWN_ERR)
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
     1                   + ddotfile_1(d_r1,d_r1,size_x1)
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
     1                   + ddotfile_2(d_r2,d_r2,size_x2)
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
          if((.not. ipccsd) .and. (.not. eaccsd)) then
           if (needx1) call tce_jacobi_x1(d_r1,k_x1_offset)
ckbn --- IPCCSD
#ifdef IPCCSD
          else if(ipccsd) then
           if (needx1) call tce_jacobi_ipx1(d_r1,k_x1_offset)
#endif
ckbn --- EACCSD
#ifdef EACCSD
          else if(eaccsd) then
           if (needx1) call tce_jacobi_eax1(d_r1,k_x1_offset)
#endif
          endif
          if((.not. ipccsd) .and. (.not. eaccsd)) then
           if (needx2) call tce_jacobi_x2(d_r2,k_x2_offset)
ckbn --- IPCCSD
#ifdef IPCCSD
          else if(ipccsd) then
           if (needx2) call tce_jacobi_ipx2(d_r2,k_x2_offset)
#endif
ckbn --- EACCSD
#ifdef EACCSD
          else if(eaccsd) then
           if (needx2) call tce_jacobi_eax2(d_r2,k_x2_offset)
#endif
          endif
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
          if (needx1) overlap=overlap+ddotfile_1(d_r1,x1(jvec),size_x1)
          if (needx2) overlap=overlap+ddotfile_2(d_r2,x2(jvec),size_x2)
            if (needx3) overlap=overlap+ddotfile(d_r3,x3(jvec),size_x3)
            if (needx4) overlap=overlap+ddotfile(d_r4,x4(jvec),size_x4)
            if (needx1) call daxpyfile(1,-overlap,x1(jvec),d_r1,size_x1)
            if (needx2) call daxpyfile(1,-overlap,x2(jvec),d_r2,size_x2)
            if (needx3) call daxpyfile(1,-overlap,x3(jvec),d_r3,size_x3)
            if (needx4) call daxpyfile(1,-overlap,x4(jvec),d_r4,size_x4)
          enddo
          overlap = 0.0d0
          if (needx1) overlap=overlap+ddotfile_1(d_r1,d_r1,size_x1)
          if (needx2) overlap=overlap+ddotfile_2(d_r2,d_r2,size_x2)
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
            call copyfile(d_r1,xc1(ivec),size_x1)
            call gatodra(xc1(ivec))
            xc1_exist(ivec) = .true.
          endif
          if (needx2) then
            call tce_zero(d_r2,size_x2)
            call daxpyfile(nxtrials,dbl_mb(k_vr+(ivec-1)*nxtrials),
     1        x2,d_r2,size_x2)
            call copyfile(d_r2,xc2(ivec),size_x2)
            call gatodra(xc2(ivec))
            xc2_exist(ivec) = .true.
          endif
          if (needx3) then
            call tce_zero(d_r3,size_x3)
            call daxpyfile(nxtrials,dbl_mb(k_vr+(ivec-1)*nxtrials),
     1        x3,d_r3,size_x3)
            call copyfile(d_r3,xc3(ivec),size_x3)
            call gatodra(xc3(ivec))
            xc3_exist(ivec) = .true.
          endif
          if (needx4) then
            call tce_zero(d_r4,size_x4)
            call daxpyfile(nxtrials,dbl_mb(k_vr+(ivec-1)*nxtrials),
     1        x4,d_r4,size_x4)
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
      if (.not.ma_pop_stack(l_shbar))
     1  call errquit('tce_eom_xdiagon: MA problem',11,MA_ERR)
 9000 format(1x,f20.10,' hartree ',f20.10,' eV')
 9010 format(1x,'Root ',i3,' converged')
 9020 format(1x,A,' cpu:',f10.2,' wall:',f10.2)
 9030 format(1x,'IVEC = ',i3)
      return
      end
c
c
c
c
c
c
      subroutine deltax_ort_norm(k_x1_offset,k_x2_offset,
     1                           size_x1,size_x2,
     1                           imicro,xdiis)
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
      integer imicro,xdiis
      integer  i,j
      integer d_ex1,size_ex1
      integer l_ex1_offset,k_ex1_offset
      integer k_x1_offset,k_x2_offset
      integer size_x1,size_x2
      double precision over,over1,over2
      character*255 filename
c
      irrep_c=irrep_x
      irrep_d=irrep_x
c
      call tce_e_offset(l_ex1_offset,k_ex1_offset,size_ex1)
      call createfile(filename,d_ex1,size_ex1)
      call tce_zero(d_ex1,size_ex1)
      call tce_filename('ex1',filename)
c
      do i=1,imicro-1
       over=0.0d0
       over1=0.0d0
       over2=0.0d0
       call tce_zero(d_ex1,size_ex1)
       call c1_d1(x1(imicro),x1(i),d_ex1,k_x1_offset,k_x1_offset,
     1            k_ex1_offset)
       call reconcilefile(d_ex1,1)
       call get_block(d_ex1,over1,1,0)
       call tce_zero(d_ex1,size_ex1)
       call c2_d2(x2(imicro),x2(i),d_ex1,k_x2_offset,k_x2_offset,
     1            k_ex1_offset)
       call reconcilefile(d_ex1,1)
       call get_block(d_ex1,over2,1,0)
       over=over1+over2
c fix this and change header
       call daxpyfile(1,(-1.0d0)*over,x1(i),
     1                x1(imicro),size_x1)
       call daxpyfile(1,(-1.0d0)*over,x2(i),
     1                x2(imicro),size_x2)
      enddo
c normalization of x(imicro)
       over=0.0d0
       over1=0.0d0
       over2=0.0d0
       call tce_zero(d_ex1,size_ex1)
       call c1_c1(x1(imicro),d_ex1,k_x1_offset,k_ex1_offset)
       call reconcilefile(d_ex1,1)
       call get_block(d_ex1,over1,1,0)
       call tce_zero(d_ex1,size_ex1)
       call c2_c2(x2(imicro),d_ex1,k_x2_offset,k_ex1_offset)
       call reconcilefile(d_ex1,1)
       call get_block(d_ex1,over2,1,0)
       over=over1+over2
       over=dsqrt(over)
       call dscalfile(1.0d0/over,x1(imicro),size_x1)
       call dscalfile(1.0d0/over,x2(imicro),size_x2)
c
c
      call deletefile(d_ex1)
      if (.not.ma_pop_stack(l_ex1_offset))
     1    call errquit("tce_energy: deall1",24,MA_ERR)
c
      return 
      end 
c
c
c
c
c
      subroutine form_best_guess(k_x1_offset,k_x2_offset,
     &                     size_x1,size_x2,
     &                     c_eff,imicro_x,xdiis,iroot)
c
c $Id$
c
c  xc(iroot) = sum(i=1,...,imicro) x(i)
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
      integer imicro_x,xdiis,residual
      integer i,j
      integer d_ex1,size_ex1
      integer l_ex1_offset,k_ex1_offset
      integer k_x1_offset,k_x2_offset
      integer size_x1,size_x2
      integer iroot
      double precision over,over1,over2
      double precision c_eff(xdiis)
c
      irrep_c=irrep_x
      irrep_d=irrep_x
c
c
      call tce_zero(xc1(iroot),size_x1)
      call tce_zero(xc2(iroot),size_x2)
c
      do i=1,imicro_x
       call daxpyfile(1,c_eff(i),x1(i),
     1                xc1(iroot),size_x1)
       call daxpyfile(1,c_eff(i),x2(i),
     1                xc2(iroot),size_x2)
      enddo
c
c
      return 
      end 
c
c
c
c
c
      subroutine form_correction(k_x1_offset,k_x2_offset,
     &                     size_x1,size_x2,
     &                     c_eff,omega,imicro,xdiis,residual)
c
c $Id$
c
c  delta_x = (omega - H-bar) xc
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
      integer imicro,xdiis
      integer i,j
      integer d_ex1,size_ex1
      integer l_ex1_offset,k_ex1_offset
      integer k_x1_offset,k_x2_offset
      integer size_x1,size_x2
      double precision over,over1,over2
      double precision c_eff(xdiis),omega,residual
      character*255 filename
      logical nodezero         ! True if node 0
c
       nodezero=(ga_nodeid().eq.0)
c
      irrep_c=irrep_x
      irrep_d=irrep_x
c
c
      do i=1,imicro-1
       call daxpyfile(1,(-1.0d0)*omega*c_eff(i),x1(i),
     1                x1(imicro),size_x1)
       call daxpyfile(1,c_eff(i),xp1(i),
     1                x1(imicro),size_x1) 
       call daxpyfile(1,(-1.0d0)*omega*c_eff(i),x2(i),
     1                x2(imicro),size_x2)
       call daxpyfile(1,c_eff(i),xp2(i),
     1                x2(imicro),size_x2)
c  
       if(dabs(zlshift).gt.0.0d0) then !-------
         call daxpyfile(1,zlshift*c_eff(i),x1(i),
     1                  x1(imicro),size_x1)
         call daxpyfile(1,zlshift*c_eff(i),x2(i),
     1                  x2(imicro),size_x2)
       end if                          !-------
c
      enddo
c 
c calculating residual
c
      over1=0.0d0
      over2=0.0d0
c
      call tce_e_offset(l_ex1_offset,k_ex1_offset,size_ex1)
      call createfile(filename,d_ex1,size_ex1)
      call tce_zero(d_ex1,size_ex1)
      call tce_filename('ex1',filename)
c
      call c1_c1(x1(imicro),d_ex1,k_x1_offset,k_ex1_offset)
      call reconcilefile(d_ex1,1)
      call get_block(d_ex1,over1,1,0)
c
      call tce_zero(d_ex1,size_ex1)
      call c2_c2(x2(imicro),d_ex1,k_x2_offset,k_ex1_offset)
      call reconcilefile(d_ex1,1)
      call get_block(d_ex1,over2,1,0)
      residual=dsqrt(over1+over2)
c
c
      call deletefile(d_ex1)
      if (.not.ma_pop_stack(l_ex1_offset))
     1    call errquit("tce_energy: deall2",24,MA_ERR)
c
c JACOBI
c
      call tce_jacobi_x1(x1(imicro),k_x1_offset,omega)
      call tce_jacobi_x2(x2(imicro),k_x2_offset,omega)
c
c
c
      return 
      end 
c
c
c
c
c
c
      subroutine form_heff(k_x1_offset,k_x2_offset,size_x1,size_x2,
     &                     c_eff,h_eff,omega,imicro,xdiis,iter_eomcc)
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
      integer imicro,xdiis,residual
      integer  i,j,iter_eomcc
      integer d_ex1,size_ex1
      integer l_ex1_offset,k_ex1_offset
      integer k_x1_offset,k_x2_offset
      integer size_x1,size_x2
      double precision over,over1,over2
      double precision h_eff(xdiis,xdiis),c_eff(xdiis) 
      double precision h_eff_copy(xdiis,xdiis)
      double precision omega,xmax
c
      character*1 jobvl,jobvr
      integer lwork,info_x
      double precision wr(xdiis),wi(xdiis)
      double precision vl(xdiis,xdiis),vr(xdiis,xdiis)
      double precision work(4*xdiis)
c 
      character*255 filename
c
      logical nodezero         ! True if node 0
c
       nodezero=(ga_nodeid().eq.0)
c
      jobvl='N'
      jobvr='V'
      lwork=4*xdiis
c
      irrep_c=irrep_x
      irrep_d=irrep_x
c
      call tce_e_offset(l_ex1_offset,k_ex1_offset,size_ex1)
      call createfile(filename,d_ex1,size_ex1)
      call tce_zero(d_ex1,size_ex1)
      call tce_filename('ex1',filename)
c
      do i=1,imicro-1 !-------- i-loop
       over=0.0d0
       over1=0.0d0
       over2=0.0d0
c imicro - row
       call tce_zero(d_ex1,size_ex1)
       call c1_d1(x1(imicro),xp1(i),d_ex1,k_x1_offset,k_x1_offset,
     1            k_ex1_offset)
       call reconcilefile(d_ex1,1)
       call get_block(d_ex1,over1,1,0)
       call tce_zero(d_ex1,size_ex1)
       call c2_d2(x2(imicro),xp2(i),d_ex1,k_x2_offset,k_x2_offset,
     1            k_ex1_offset)
       call reconcilefile(d_ex1,1)
       call get_block(d_ex1,over2,1,0)
       over=over1+over2
c fix this and change header
       h_eff(imicro,i)=over
c normalization of x(imicro)
       over=0.0d0
       over1=0.0d0
       over2=0.0d0
c imicro - column
       call tce_zero(d_ex1,size_ex1)
       call c1_d1(x1(i),xp1(imicro),d_ex1,k_x1_offset,k_x1_offset,
     1            k_ex1_offset)
       call reconcilefile(d_ex1,1)
       call get_block(d_ex1,over1,1,0)
       call tce_zero(d_ex1,size_ex1)
       call c2_d2(x2(i),xp2(imicro),d_ex1,k_x2_offset,k_x2_offset,
     1            k_ex1_offset)
       call reconcilefile(d_ex1,1)
       call get_block(d_ex1,over2,1,0)
       over=over1+over2
c fix this and change header
       h_eff(i,imicro)=over
c normalization of x(imicro)
       enddo !-------------------- i-loop
c h_eff(imicro,imicro)
       over=0.0d0
       over1=0.0d0
       over2=0.0d0
c nxrtials - column
       call tce_zero(d_ex1,size_ex1)
       call c1_d1(x1(imicro),xp1(imicro),d_ex1,
     1            k_x1_offset,k_x1_offset,
     1            k_ex1_offset)
       call reconcilefile(d_ex1,1)
       call get_block(d_ex1,over1,1,0)
       call tce_zero(d_ex1,size_ex1)
       call c2_d2(x2(imicro),xp2(imicro),d_ex1,
     1            k_x2_offset,k_x2_offset,
     1            k_ex1_offset)
       call reconcilefile(d_ex1,1)
       call get_block(d_ex1,over2,1,0)
       over=over1+over2
c fix this and change header
       h_eff(imicro,imicro)=over
c
      call deletefile(d_ex1)
      if (.not.ma_pop_stack(l_ex1_offset))
     1    call errquit("tce_energy: deall3",24,MA_ERR)
c
       do i=1,xdiis
       do j=1,xdiis
        h_eff_copy(i,j)=h_eff(i,j)
       enddo
       enddo
c big number on diagonal beyond current space
       do i=imicro+1,xdiis
         h_eff(i,i)=10000000000.0d0
       enddo
c
c
c
        call util_dgeev(jobvl,jobvr,xdiis,h_eff,xdiis,wr,wi,
     &             vl,xdiis,vr,xdiis,work,lwork,info_x)
c
c selection
c
       xmax=0.0d0
       j=0
       do i=1,xdiis
        if((dabs(vr(1,i)).gt.xmax).and.(dabs(wi(i)).lt.0.000010d0)) then
         xmax=dabs(vr(1,i))
         j=i
        end if
       enddo
       if(j.eq.0) then
        call errquit('Eigen-solver - imaginary roots',xdiis,CALC_ERR)
       end if       
c
       xmax=0.0d0
       do i=1,xdiis
        if(i.ne.j) then
         if(dabs(vr(1,i)).gt.xmax) xmax=dabs(vr(1,i))
        end if
       enddo
       if((xmax.gt.0.10d0).and.(iter_eomcc.gt.xdiis)) then
        zlshift=0.20d0
        if(nodezero) then
         write(6,*)'   level shift is applied ............'
         call util_flush(6)
        end if 
       end if
c
       xmax=0.0d0
       omega=wr(j)
       do i=1,xdiis
        c_eff(i)=vr(i,j)
        xmax=xmax+vr(i,j)*vr(i,j)
       enddo
       do i=1,xdiis
        c_eff(i)=(1.0d0/dsqrt(xmax))*c_eff(i)
       enddo
c
c
       do i=1,xdiis
       do j=1,xdiis
        h_eff(i,j)=h_eff_copy(i,j)
       enddo
       enddo
c 
      return 
      end 
c
c
c
c
c
      subroutine norm_xc_vector(k_x1_offset,k_x2_offset,
     1                           size_x1,size_x2,iroot)
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
      integer iroot
      integer i,j
      integer d_ex1,size_ex1
      integer l_ex1_offset,k_ex1_offset
      integer k_x1_offset,k_x2_offset
      integer size_x1,size_x2
      double precision over,over1,over2
      character*255 filename
c
      irrep_c=irrep_x
      irrep_d=irrep_x
c
      call tce_e_offset(l_ex1_offset,k_ex1_offset,size_ex1)
      call createfile(filename,d_ex1,size_ex1)
      call tce_zero(d_ex1,size_ex1)
      call tce_filename('ex1',filename)
c normalization of x(imicro)
       over=0.0d0
       over1=0.0d0
       over2=0.0d0
       call tce_zero(d_ex1,size_ex1)
       call c1_c1(xc1(iroot),d_ex1,k_x1_offset,k_ex1_offset)
       call reconcilefile(d_ex1,1)
       call get_block(d_ex1,over1,1,0)
       call tce_zero(d_ex1,size_ex1)
       call c2_c2(xc2(iroot),d_ex1,k_x2_offset,k_ex1_offset)
       call reconcilefile(d_ex1,1)
       call get_block(d_ex1,over2,1,0)
       over=over1+over2
       over=dsqrt(over)
       call dscalfile(1.0d0/over,xc1(iroot),size_x1)
       call dscalfile(1.0d0/over,xc2(iroot),size_x2)
c
c
      call deletefile(d_ex1)
      if (.not.ma_pop_stack(l_ex1_offset))
     1    call errquit("tce_energy: deall4",24,MA_ERR)
c
      return 
      end 
c
c
c
      subroutine eom_starts(size_x1,k_x1_offset,d_f1,k_f1_offset,d_v2,
     1 k_v2_offset,k_omegax)
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
c
      integer size_x1
      integer k_x1_offset
      integer d_f1
      integer k_f1_offset
      integer d_v2
      integer k_v2_offset
      integer k_omegax,l_omegax
      integer k_residual,l_residual
      logical nodezero
      integer ivec,i,j
      integer iopen_f
c
      character*1 jobvl,jobvr
      integer lwork,info_x
      integer k_matrix,l_matrix
      integer k_wr,l_wr,k_wi,l_wi,k_vl,l_vl,k_vr,l_vr
      integer k_work,l_work
c
      integer d_ex1,size_ex1
      integer l_ex1_offset,k_ex1_offset
      double precision over,xmin
c
      double precision cpu, wall
      double precision au2ev
      parameter (au2ev=27.2113961d0)
      character*255 filename
c
      nodezero=(ga_nodeid().eq.0)
c
      irrep_c=irrep_x
      irrep_d=irrep_x
c
      iopen_f=0
      call cis_xguess(size_x1,k_x1_offset,iopen_f) !nxtrials x1 vectors created
        do ivec=1,nroots_reduced
          call tce_filenameindexed(ivec,'x1_ini',filename)
          call createfile(filename,x1_ini(ivec),size_x1)
          call tce_zero(x1_ini(ivec),size_x1)
          x1_ini_exist(ivec) = .true.
        enddo
      iopen_f=1
      call cis_xguess(size_x1,k_x1_offset,iopen_f) !nxtrials x1 vectors created
      do ivec=1,nxtrials !nxtrials x1 vectors formed here
        call tce_filenameindexed(ivec,'xp1',filename)
        call createfile(filename,xp1(ivec),size_x1)
        call tce_zero(xp1(ivec),size_x1)
        xp1_exist(ivec) = .true.
        call cis_x1(d_f1,xp1(ivec),d_v2,x1(ivec),
     1              k_f1_offset,k_x1_offset,k_v2_offset,k_x1_offset)
        call reconcilefile(xp1(ivec),size_x1)
      enddo 
c
      jobvl='N'
      jobvr='V'
      lwork=4*nxtrials
c
c  opening dgeev matrices
       if (.not.MA_PUSH_GET(mt_dbl,nxtrials*nxtrials,'matrix',
     1      l_matrix,k_matrix)) call errquit('k_matrix',1,MA_ERR)
       if (.not.MA_PUSH_GET(mt_dbl,nxtrials,'wr',
     1      l_wr,k_wr)) call errquit('k_wr',1,MA_ERR)
       if (.not.MA_PUSH_GET(mt_dbl,nxtrials,'wi',
     1      l_wi,k_wi)) call errquit('k_wi',1,MA_ERR)
       if (.not.MA_PUSH_GET(mt_dbl,nxtrials*nxtrials,'l_vl',
     1      l_vl,k_vl)) call errquit('k_vl',1,MA_ERR)
       if (.not.MA_PUSH_GET(mt_dbl,nxtrials*nxtrials,'l_vr',
     1      l_vr,k_vr)) call errquit('k_vr',1,MA_ERR)
       if (.not.MA_PUSH_GET(mt_dbl,4*nxtrials,'l_work',
     1      l_work,k_work)) call errquit('k_work',1,MA_ERR)
c form matrix
       do i=1,nxtrials*nxtrials
        dbl_mb(k_matrix+i-1)=0.0d0
       enddo 
c
       call tce_e_offset(l_ex1_offset,k_ex1_offset,size_ex1)
       call createfile(filename,d_ex1,size_ex1)
       call tce_zero(d_ex1,size_ex1)
       call tce_filename('ex1',filename)
c
       do i=1,nxtrials
       do j=1,nxtrials
        over=0.0d0
        call tce_zero(d_ex1,size_ex1)
        call c1_d1(x1(i),xp1(j),d_ex1,
     1            k_x1_offset,k_x1_offset,
     1            k_ex1_offset)
        call reconcilefile(d_ex1,1)
        call get_block(d_ex1,over,1,0)
        dbl_mb(k_matrix+(j-1)*nxtrials+i-1)=over
       enddo
       enddo
c
       call deletefile(d_ex1)
      if (.not.ma_pop_stack(l_ex1_offset))
     1    call errquit("tce_energy: deall5",24,MA_ERR)
c
c  diagonalize matrix
        call util_dgeev(jobvl,jobvr,nxtrials,dbl_mb(k_matrix),
     &   nxtrials,dbl_mb(k_wr),dbl_mb(k_wi),
     &   dbl_mb(k_vl),nxtrials,dbl_mb(k_vr),
     &   nxtrials,dbl_mb(k_work),lwork,info_x)
c nroots_reduced lowest energies
       do ivec=1,nroots_reduced
          xmin=100000.0d0 
          j=0
          do i=1,nxtrials
           if(dbl_mb(k_wr+i-1).lt.xmin) then
            xmin=dbl_mb(k_wr+i-1)
            j=i
           end if
          enddo
          dbl_mb(k_omegax+ivec-1)=dbl_mb(k_wr+j-1)
          dbl_mb(k_wr+j-1)=1000000.0d0
          do i=1,nxtrials 
           call daxpyfile(1,dbl_mb(k_vr+(j-1)*nxtrials+i-1),x1(i),
     1                    x1_ini(ivec),size_x1)
          enddo
       enddo
c
       if(.not.MA_POP_STACK(l_work))
     &    call errquit('l_work_e',4,MA_ERR)
       if(.not.MA_POP_STACK(l_vr))
     &    call errquit('l_vr_e',4,MA_ERR)
       if(.not.MA_POP_STACK(l_vl))
     &    call errquit('l_vl_e',4,MA_ERR)
       if(.not.MA_POP_STACK(l_wi))
     &    call errquit('l_wi_e',4,MA_ERR)
       if(.not.MA_POP_STACK(l_wr))
     &    call errquit('l_wr_e',4,MA_ERR)
       if(.not.MA_POP_STACK(l_matrix))
     &    call errquit('l_matrix_e',4,MA_ERR)
c
         do ivec=1,maxtrials
           if(x1_exist(ivec))  then
            call deletefile(x1(ivec))
            x1_exist(ivec) = .false.
           end if
           if(xp1_exist(ivec))  then
            call deletefile(xp1(ivec))
            xp1_exist(ivec) = .false.
           end if
         enddo
c
       return
       end
c
c
c
c
c
c
      subroutine deltax_ort_norm_act(k_x1_offset,k_x2_offset,
     1                           size_x1,size_x2,
     1                           imicro,xdiis)
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
      integer imicro,xdiis
      integer  i,j
      integer d_ex1,size_ex1
      integer l_ex1_offset,k_ex1_offset
      integer k_x1_offset,k_x2_offset
      integer size_x1,size_x2
      double precision over,over1,over2
      character*255 filename
c
      irrep_c=irrep_x
      irrep_d=irrep_x
c
      call tce_e_offset(l_ex1_offset,k_ex1_offset,size_ex1)
      call createfile(filename,d_ex1,size_ex1)
      call tce_zero(d_ex1,size_ex1)
      call tce_filename('ex1',filename)
c
      do i=1,imicro-1
       over=0.0d0
       over1=0.0d0
       over2=0.0d0
       call tce_zero(d_ex1,size_ex1)
       call c1_d1(x1(imicro),x1(i),d_ex1,k_x1_offset,k_x1_offset,
     1            k_ex1_offset)
       call reconcilefile(d_ex1,1)
       call get_block(d_ex1,over1,1,0)
       call tce_zero(d_ex1,size_ex1)
       call c2_d2_act(x2(imicro),x2(i),d_ex1,k_x2_offset,k_x2_offset,
     1            k_ex1_offset)
       call reconcilefile(d_ex1,1)
       call get_block(d_ex1,over2,1,0)
       over=over1+over2
c fix this and change header
       call daxpyfile(1,(-1.0d0)*over,x1(i),
     1                x1(imicro),size_x1)
       call daxpyfile(1,(-1.0d0)*over,x2(i),
     1                x2(imicro),size_x2)
      enddo
c normalization of x(imicro)
       over=0.0d0
       over1=0.0d0
       over2=0.0d0
       call tce_zero(d_ex1,size_ex1)
       call c1_c1(x1(imicro),d_ex1,k_x1_offset,k_ex1_offset)
       call reconcilefile(d_ex1,1)
       call get_block(d_ex1,over1,1,0)
       call tce_zero(d_ex1,size_ex1)
       call c2_c2_act(x2(imicro),d_ex1,k_x2_offset,k_ex1_offset)
       call reconcilefile(d_ex1,1)
       call get_block(d_ex1,over2,1,0)
       over=over1+over2
       over=dsqrt(over)
       call dscalfile(1.0d0/over,x1(imicro),size_x1)
       call dscalfile(1.0d0/over,x2(imicro),size_x2)
c
c
      call deletefile(d_ex1)
      if (.not.ma_pop_stack(l_ex1_offset))
     1    call errquit("tce_energy: deall1",24,MA_ERR)
c
      return 
      end 
c
c
c
      subroutine form_correction_act(k_x1_offset,k_x2_offset,
     &                     size_x1,size_x2,
     &                     c_eff,omega,imicro,xdiis,residual)
c
c $Id$
c
c  delta_x = (omega - H-bar) xc
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
      integer imicro,xdiis
      integer i,j
      integer d_ex1,size_ex1
      integer l_ex1_offset,k_ex1_offset
      integer k_x1_offset,k_x2_offset
      integer size_x1,size_x2
      double precision over,over1,over2
      double precision c_eff(xdiis),omega,residual
      character*255 filename
      logical nodezero         ! True if node 0
c
       nodezero=(ga_nodeid().eq.0)
c
      irrep_c=irrep_x
      irrep_d=irrep_x
c
c
c
      do i=1,imicro-1
       call daxpyfile(1,(-1.0d0)*omega*c_eff(i),x1(i),
     1                x1(imicro),size_x1)
       call daxpyfile(1,c_eff(i),xp1(i),
     1                x1(imicro),size_x1) 
       call daxpyfile(1,(-1.0d0)*omega*c_eff(i),x2(i),
     1                x2(imicro),size_x2)
       call daxpyfile(1,c_eff(i),xp2(i),
     1                x2(imicro),size_x2)
      enddo
c 
c calculating residual
c
      over1=0.0d0
      over2=0.0d0
c
      call tce_e_offset(l_ex1_offset,k_ex1_offset,size_ex1)
      call createfile(filename,d_ex1,size_ex1)
      call tce_zero(d_ex1,size_ex1)
      call tce_filename('ex1',filename)
c
      call c1_c1(x1(imicro),d_ex1,k_x1_offset,k_ex1_offset)
      call reconcilefile(d_ex1,1)
      call get_block(d_ex1,over1,1,0)
c
      call tce_zero(d_ex1,size_ex1)
      call c2_c2_act(x2(imicro),d_ex1,k_x2_offset,k_ex1_offset)
      call reconcilefile(d_ex1,1)
      call get_block(d_ex1,over2,1,0)
      residual=dsqrt(over1+over2)
c
c
      call deletefile(d_ex1)
      if (.not.ma_pop_stack(l_ex1_offset))
     1    call errquit("tce_energy: deall2",24,MA_ERR)
c
c JACOBI
c
      call tce_jacobi_x1(x1(imicro),k_x1_offset,omega)
      call tce_jacobi_x2a(x2(imicro),k_x2_offset,omega)
c
c
c
      return 
      end 
c
c
c
c
c
c
      subroutine form_heff_act(k_x1_offset,k_x2_offset,size_x1,size_x2,
     &                     c_eff,h_eff,omega,imicro,xdiis)
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
      integer imicro,xdiis,residual
      integer  i,j
      integer d_ex1,size_ex1
      integer l_ex1_offset,k_ex1_offset
      integer k_x1_offset,k_x2_offset
      integer size_x1,size_x2
      double precision over,over1,over2
      double precision h_eff(xdiis,xdiis),c_eff(xdiis) 
      double precision h_eff_copy(xdiis,xdiis)
      double precision omega,xmax
c
      character*1 jobvl,jobvr
      integer lwork,info_x
      double precision wr(xdiis),wi(xdiis)
      double precision vl(xdiis,xdiis),vr(xdiis,xdiis)
      double precision work(4*xdiis)
c
      character*255 filename
c
      logical nodezero         ! True if node 0
c
       nodezero=(ga_nodeid().eq.0)
c
      jobvl='N'
      jobvr='V'
      lwork=4*xdiis
c
      irrep_c=irrep_x
      irrep_d=irrep_x
c
      call tce_e_offset(l_ex1_offset,k_ex1_offset,size_ex1)
      call createfile(filename,d_ex1,size_ex1)
      call tce_zero(d_ex1,size_ex1)
      call tce_filename('ex1',filename)
c
      do i=1,imicro-1 !-------- i-loop
       over=0.0d0
       over1=0.0d0
       over2=0.0d0
c imicro - row
       call tce_zero(d_ex1,size_ex1)
       call c1_d1(x1(imicro),xp1(i),d_ex1,k_x1_offset,k_x1_offset,
     1            k_ex1_offset)
       call reconcilefile(d_ex1,1)
       call get_block(d_ex1,over1,1,0)
       call tce_zero(d_ex1,size_ex1)
       call c2_d2_act(x2(imicro),xp2(i),d_ex1,k_x2_offset,k_x2_offset,
     1            k_ex1_offset)
       call reconcilefile(d_ex1,1)
       call get_block(d_ex1,over2,1,0)
       over=over1+over2
c fix this and change header
       h_eff(imicro,i)=over
c normalization of x(imicro)
       over=0.0d0
       over1=0.0d0
       over2=0.0d0
c imicro - column
       call tce_zero(d_ex1,size_ex1)
       call c1_d1(x1(i),xp1(imicro),d_ex1,k_x1_offset,k_x1_offset,
     1            k_ex1_offset)
       call reconcilefile(d_ex1,1)
       call get_block(d_ex1,over1,1,0)
       call tce_zero(d_ex1,size_ex1)
       call c2_d2_act(x2(i),xp2(imicro),d_ex1,k_x2_offset,k_x2_offset,
     1            k_ex1_offset)
       call reconcilefile(d_ex1,1)
       call get_block(d_ex1,over2,1,0)
       over=over1+over2
c fix this and change header
       h_eff(i,imicro)=over
c normalization of x(imicro)
       enddo !-------------------- i-loop
c h_eff(imicro,imicro)
       over=0.0d0
       over1=0.0d0
       over2=0.0d0
c nxrtials - column
       call tce_zero(d_ex1,size_ex1)
       call c1_d1(x1(imicro),xp1(imicro),d_ex1,
     1            k_x1_offset,k_x1_offset,
     1            k_ex1_offset)
       call reconcilefile(d_ex1,1)
       call get_block(d_ex1,over1,1,0)
       call tce_zero(d_ex1,size_ex1)
       call c2_d2_act(x2(imicro),xp2(imicro),d_ex1,
     1            k_x2_offset,k_x2_offset,
     1            k_ex1_offset)
       call reconcilefile(d_ex1,1)
       call get_block(d_ex1,over2,1,0)
       over=over1+over2
c fix this and change header
       h_eff(imicro,imicro)=over
c
      call deletefile(d_ex1)
      if (.not.ma_pop_stack(l_ex1_offset))
     1    call errquit("tce_energy: deall3",24,MA_ERR)
c
       do i=1,xdiis
       do j=1,xdiis
        h_eff_copy(i,j)=h_eff(i,j)
       enddo
       enddo
c big number on diagonal beyond current space
       do i=imicro+1,xdiis
         h_eff(i,i)=10000000000.0d0
       enddo
c
c
c
        call util_dgeev(jobvl,jobvr,xdiis,h_eff,xdiis,wr,wi,
     &             vl,xdiis,vr,xdiis,work,lwork,info_x)
c
c selection
c
       xmax=0.0d0
       j=0
       do i=1,xdiis
        if((dabs(vr(1,i)).gt.xmax).and.(dabs(wi(i)).lt.0.000010d0)) then
         xmax=dabs(vr(1,i))
         j=i
        end if
       enddo
       if(j.eq.0) then
        call errquit('Eigen-solver - imaginary roots',xdiis,CALC_ERR)
       end if       
c
       xmax=0.0d0
       omega=wr(j)
       do i=1,xdiis
        c_eff(i)=vr(i,j)
        xmax=xmax+vr(i,j)*vr(i,j)
       enddo
       do i=1,xdiis
        c_eff(i)=(1.0d0/dsqrt(xmax))*c_eff(i)
       enddo
c
       do i=1,xdiis
       do j=1,xdiis
        h_eff(i,j)=h_eff_copy(i,j)
       enddo
       enddo
c 
      return 
      end 
c
c
c
c
c
      subroutine norm_xc_vector_act(k_x1_offset,k_x2_offset,
     1                           size_x1,size_x2,iroot)
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
      integer iroot
      integer i,j
      integer d_ex1,size_ex1
      integer l_ex1_offset,k_ex1_offset
      integer k_x1_offset,k_x2_offset
      integer size_x1,size_x2
      double precision over,over1,over2
      character*255 filename
c
      irrep_c=irrep_x
      irrep_d=irrep_x
c
      call tce_e_offset(l_ex1_offset,k_ex1_offset,size_ex1)
      call createfile(filename,d_ex1,size_ex1)
      call tce_zero(d_ex1,size_ex1)
      call tce_filename('ex1',filename)
c normalization of x(imicro)
       over=0.0d0
       over1=0.0d0
       over2=0.0d0
       call tce_zero(d_ex1,size_ex1)
       call c1_c1(xc1(iroot),d_ex1,k_x1_offset,k_ex1_offset)
       call reconcilefile(d_ex1,1)
       call get_block(d_ex1,over1,1,0)
       call tce_zero(d_ex1,size_ex1)
       call c2_c2_act(xc2(iroot),d_ex1,k_x2_offset,k_ex1_offset)
       call reconcilefile(d_ex1,1)
       call get_block(d_ex1,over2,1,0)
       over=over1+over2
       over=dsqrt(over)
       call dscalfile(1.0d0/over,xc1(iroot),size_x1)
       call dscalfile(1.0d0/over,xc2(iroot),size_x2)
c
c
      call deletefile(d_ex1)
      if (.not.ma_pop_stack(l_ex1_offset))
     1    call errquit("tce_energy: deall4",24,MA_ERR)
c
      return 
      end 
c
c
c
c
c
      subroutine norm_x_vector_mr(k_x1_offset,k_x2_offset,
     1                           size_x1,size_x2,totr)
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
      integer iroot
      integer i,j,totr
      integer d_ex1,size_ex1
      integer l_ex1_offset,k_ex1_offset
      integer k_x1_offset,k_x2_offset
      integer size_x1,size_x2
      double precision over,over1,over2
      character*255 filename
c
      irrep_c=irrep_x
      irrep_d=irrep_x
c
      call tce_e_offset(l_ex1_offset,k_ex1_offset,size_ex1)
      call createfile(filename,d_ex1,size_ex1)
      call tce_zero(d_ex1,size_ex1)
      call tce_filename('ex1',filename)
c normalization of initial  space
      do i=1,totr
c
       do j=1,i-1
        over=0.0d0
        over1=0.0d0
        over2=0.0d0
        call tce_zero(d_ex1,size_ex1)
        call c1_d1(x1(i),x1(j),d_ex1,k_x1_offset,k_x1_offset,
     1             k_ex1_offset)
        call reconcilefile(d_ex1,1)
        call get_block(d_ex1,over1,1,0)
        call tce_zero(d_ex1,size_ex1)
        call c2_d2(x2(i),x2(j),d_ex1,k_x2_offset,k_x2_offset,
     1             k_ex1_offset)
        call reconcilefile(d_ex1,1)
        call get_block(d_ex1,over2,1,0)
        over=over1+over2
c fix this and change header
        call daxpyfile(1,(-1.0d0)*over,x1(j),
     1                 x1(i),size_x1)
        call daxpyfile(1,(-1.0d0)*over,x2(j),
     1                 x2(i),size_x2)
       enddo ! ---- over j
c
c normalization of x(i)
       over=0.0d0
       over1=0.0d0
       over2=0.0d0
       call tce_zero(d_ex1,size_ex1)
       call c1_c1(x1(i),d_ex1,k_x1_offset,k_ex1_offset)
       call reconcilefile(d_ex1,1)
       call get_block(d_ex1,over1,1,0)
       call tce_zero(d_ex1,size_ex1)
       call c2_c2(x2(i),d_ex1,k_x2_offset,k_ex1_offset)
       call reconcilefile(d_ex1,1)
       call get_block(d_ex1,over2,1,0)
       over=over1+over2
       over=dsqrt(over)
       call dscalfile(1.0d0/over,x1(i),size_x1)
       call dscalfile(1.0d0/over,x2(i),size_x2)
c
      enddo !over i ------
c
      call deletefile(d_ex1)
      if (.not.ma_pop_stack(l_ex1_offset))
     1    call errquit("tce_energy: deall1",24,MA_ERR)
c
      return 
      end 
c
c
c
c
      subroutine form_heff_mr(k_x1_offset,k_x2_offset,size_x1,size_x2,
     &                      c_eff,h_eff,omega,imicro,xdiis,iter_eomcc)
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
      integer imicro,xdiis,residual
      integer  i,j,iter_eomcc,ipos,ipos1,ipos2,imin
      integer d_ex1,size_ex1
      integer l_ex1_offset,k_ex1_offset
      integer k_x1_offset,k_x2_offset
      integer size_x1,size_x2
      double precision over,over1,over2
      double precision h_eff(xdiis*nroots_reduced,xdiis*nroots_reduced)
      double precision c_eff(xdiis*nroots_reduced,nroots_reduced) 
      double precision h_copy(xdiis*nroots_reduced,xdiis*nroots_reduced)
      double precision omega(nroots_reduced),xmax,xmin
c
      character*1 jobvl,jobvr
      integer lwork,info_x
      double precision wr(xdiis*nroots_reduced),wi(xdiis*nroots_reduced)
      double precision vl(xdiis*nroots_reduced,xdiis*nroots_reduced)
      double precision vr(xdiis*nroots_reduced,xdiis*nroots_reduced)
      double precision work(4*xdiis*nroots_reduced)
c 
      character*255 filename
c
      logical nodezero         ! True if node 0
c
       nodezero=(ga_nodeid().eq.0)
c
      jobvl='N'
      jobvr='V'
      lwork=4*xdiis*nroots_reduced
c
      irrep_c=irrep_x
      irrep_d=irrep_x
c 
      call tce_e_offset(l_ex1_offset,k_ex1_offset,size_ex1)
      call createfile(filename,d_ex1,size_ex1)
      call tce_zero(d_ex1,size_ex1)
      call tce_filename('ex1',filename)
c
      do i=1,(imicro-1)*nroots_reduced !-------- i-loop
c
       do j=1,nroots_reduced
        over=0.0d0
        over1=0.0d0
        over2=0.0d0
        ipos=(imicro-1)*nroots_reduced+j
c j related to column
       call tce_zero(d_ex1,size_ex1)
       call c1_d1(x1(i),xp1(ipos),d_ex1,k_x1_offset,k_x1_offset,
     1            k_ex1_offset)
       call reconcilefile(d_ex1,1)
       call get_block(d_ex1,over1,1,0)
       call tce_zero(d_ex1,size_ex1)
       call c2_d2(x2(i),xp2(ipos),d_ex1,k_x2_offset,k_x2_offset,
     1            k_ex1_offset)
       call reconcilefile(d_ex1,1)
       call get_block(d_ex1,over2,1,0)
       over=over1+over2
c fix this and change header
       h_eff(i,ipos)=over
c normalization of x(imicro)
       over=0.0d0
       over1=0.0d0
       over2=0.0d0
c imicro - column
       call tce_zero(d_ex1,size_ex1)
       call c1_d1(x1(ipos),xp1(i),d_ex1,k_x1_offset,k_x1_offset,
     1            k_ex1_offset)
       call reconcilefile(d_ex1,1)
       call get_block(d_ex1,over1,1,0)
       call tce_zero(d_ex1,size_ex1)
       call c2_d2(x2(ipos),xp2(i),d_ex1,k_x2_offset,k_x2_offset,
     1            k_ex1_offset)
       call reconcilefile(d_ex1,1)
       call get_block(d_ex1,over2,1,0)
       over=over1+over2
c fix this and change header
       h_eff(ipos,i)=over
c normalization of x(imicro)

       enddo !-------------------- j-loop
       enddo !-------------------- i-loop
c h_eff(imicro,imicro)
       do i=1,nroots_reduced
       do j=1,nroots_reduced
       ipos1=(imicro-1)*nroots_reduced+i
       ipos2=(imicro-1)*nroots_reduced+j
       over=0.0d0
       over1=0.0d0
       over2=0.0d0
c nxrtials - column
       call tce_zero(d_ex1,size_ex1)
       call c1_d1(x1(ipos1),xp1(ipos2),d_ex1,
     1            k_x1_offset,k_x1_offset,
     1            k_ex1_offset)
       call reconcilefile(d_ex1,1)
       call get_block(d_ex1,over1,1,0)
       call tce_zero(d_ex1,size_ex1)
       call c2_d2(x2(ipos1),xp2(ipos2),d_ex1,
     1            k_x2_offset,k_x2_offset,
     1            k_ex1_offset)
       call reconcilefile(d_ex1,1)
       call get_block(d_ex1,over2,1,0)
       over=over1+over2
c fix this and change header
       h_eff(ipos1,ipos2)=over
c
      enddo !j ---------
      enddo !i ---------
c
      call deletefile(d_ex1)
      if (.not.ma_pop_stack(l_ex1_offset))
     1    call errquit("tce_energy: deall3",24,MA_ERR)
c
       do i=1,xdiis*nroots_reduced
       do j=1,xdiis*nroots_reduced
        h_copy(i,j)=h_eff(i,j)
       enddo
       enddo
c big number on diagonal beyond current space
       do i=imicro*nroots_reduced+1,xdiis*nroots_reduced
         h_eff(i,i)=10000000000.0d0
       enddo
c
c
        call util_dgeev(jobvl,jobvr,xdiis*nroots_reduced,h_eff,
     &             xdiis*nroots_reduced,wr,wi,
     &             vl,xdiis*nroots_reduced,vr,
     &             xdiis*nroots_reduced,work,lwork,info_x)
c
c selection
c
       do j=1,nroots_reduced ! --- i
       xmin=10000000001.0d0
       imin=0
       do i=1,xdiis*nroots_reduced
        if((wr(i).lt.xmin).and.(dabs(wi(i)).lt.0.000010d0)) then
         xmin=wr(i)
         imin=i
        end if
       enddo
       if(imin.eq.0) then
        call errquit('Eigen-solver - imaginary roots',xdiis,CALC_ERR)
       end if       
c
       omega(j)=wr(imin)
       wr(imin)=10000000001.0d0

       xmax=0.0d0
       do i=1,xdiis*nroots_reduced
        c_eff(i,j)=vr(i,imin)
        xmax=xmax+vr(i,imin)*vr(i,imin)
       enddo
       do i=1,xdiis*nroots_reduced
        c_eff(i,j)=(1.0d0/dsqrt(xmax))*c_eff(i,j)
       enddo
c
       enddo ! j ---
c
c
       do i=1,xdiis*nroots_reduced
       do j=1,xdiis*nroots_reduced
        h_eff(i,j)=h_copy(i,j)
       enddo
       enddo
c 
      return 
      end 
c
c
c
      subroutine form_correction_mr(k_x1_offset,k_x2_offset,
     &                     size_x1,size_x2,
     &                     c_eff,omega,imicro,xdiis,residual,iroot)
c
c $Id$
c
c  delta_x = (omega - H-bar) xc
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
      integer imicro,xdiis
      integer i,j,iroot
      integer d_ex1,size_ex1
      integer l_ex1_offset,k_ex1_offset
      integer k_x1_offset,k_x2_offset
      integer size_x1,size_x2
      double precision over,over1,over2
      double precision c_eff(xdiis*nroots_reduced),omega,residual
      character*255 filename
      logical nodezero         ! True if node 0
c
       nodezero=(ga_nodeid().eq.0)
c
      irrep_c=irrep_x
      irrep_d=irrep_x
c
c
      do i=1,(imicro-1)*nroots_reduced
       call daxpyfile(1,(-1.0d0)*omega*c_eff(i),x1(i),
     1                x1((imicro-1)*nroots_reduced+iroot),size_x1)
       call daxpyfile(1,c_eff(i),xp1(i),
     1                x1((imicro-1)*nroots_reduced+iroot),size_x1) 
       call daxpyfile(1,(-1.0d0)*omega*c_eff(i),x2(i),
     1                x2((imicro-1)*nroots_reduced+iroot),size_x2)
       call daxpyfile(1,c_eff(i),xp2(i),
     1                x2((imicro-1)*nroots_reduced+iroot),size_x2)
c  
c
      enddo
c 
c calculating residual
c
      over1=0.0d0
      over2=0.0d0
c
      call tce_e_offset(l_ex1_offset,k_ex1_offset,size_ex1)
      call createfile(filename,d_ex1,size_ex1)
      call tce_zero(d_ex1,size_ex1)
      call tce_filename('ex1',filename)
c
      call c1_c1(x1((imicro-1)*nroots_reduced+iroot),
     &     d_ex1,k_x1_offset,k_ex1_offset)
      call reconcilefile(d_ex1,1)
      call get_block(d_ex1,over1,1,0)
c
      call tce_zero(d_ex1,size_ex1)
      call c2_c2(x2((imicro-1)*nroots_reduced+iroot),
     &     d_ex1,k_x2_offset,k_ex1_offset)
      call reconcilefile(d_ex1,1)
      call get_block(d_ex1,over2,1,0)
      residual=dsqrt(over1+over2)
c
c 
      call deletefile(d_ex1)
      if (.not.ma_pop_stack(l_ex1_offset))
     1    call errquit("tce_energy: deall2",24,MA_ERR)
c
c JACOBI
c
      call tce_jacobi_x1(x1((imicro-1)*nroots_reduced+iroot),
     &     k_x1_offset,omega)
      call tce_jacobi_x2(x2((imicro-1)*nroots_reduced+iroot),
     &     k_x2_offset,omega)
c
c
c
      return 
      end 
c
c
c
      subroutine form_best_guess_mr(k_x1_offset,k_x2_offset,
     &                     size_x1,size_x2,
     &                     c_eff,imicro_x,xdiis,iroot)
c
c $Id$
c
c  xc(iroot) = sum(i=1,...,imicro) x(i)
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
      integer imicro_x,xdiis,residual
      integer i,j
      integer d_ex1,size_ex1
      integer l_ex1_offset,k_ex1_offset
      integer k_x1_offset,k_x2_offset
      integer size_x1,size_x2
      integer iroot
      double precision over,over1,over2
      double precision c_eff(xdiis*nroots_reduced)
      logical nodezero         ! True if node 0
c
       nodezero=(ga_nodeid().eq.0)
c

c
      irrep_c=irrep_x
      irrep_d=irrep_x
c
c
      call tce_zero(xc1(iroot),size_x1)
      call tce_zero(xc2(iroot),size_x2)
c
      do i=1,imicro_x*nroots_reduced
       call daxpyfile(1,c_eff(i),x1(i),
     1                xc1(iroot),size_x1)
       call daxpyfile(1,c_eff(i),x2(i),
     1                xc2(iroot),size_x2)
      enddo
c
c
      return 
      end 
c
c
c
      subroutine deltax_ort_norm_mr(k_x1_offset,k_x2_offset,
     1                           size_x1,size_x2,
     1                           imicro,xdiis)
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
      integer imicro,xdiis
      integer  i,j,ipos
      integer d_ex1,size_ex1
      integer l_ex1_offset,k_ex1_offset
      integer k_x1_offset,k_x2_offset
      integer size_x1,size_x2
      double precision over,over1,over2
      character*255 filename
c
      irrep_c=irrep_x
      irrep_d=irrep_x
c
      call tce_e_offset(l_ex1_offset,k_ex1_offset,size_ex1)
      call createfile(filename,d_ex1,size_ex1)
      call tce_zero(d_ex1,size_ex1)
      call tce_filename('ex1',filename)
c
      do i=1,nroots_reduced
      ipos=(imicro-1)*nroots_reduced+i
      do j=1,ipos-1
       over=0.0d0
       over1=0.0d0
       over2=0.0d0
       call tce_zero(d_ex1,size_ex1)
       call c1_d1(x1(ipos),x1(j),d_ex1,k_x1_offset,k_x1_offset,
     1            k_ex1_offset)
       call reconcilefile(d_ex1,1)
       call get_block(d_ex1,over1,1,0)
       call tce_zero(d_ex1,size_ex1)
       call c2_d2(x2(ipos),x2(j),d_ex1,k_x2_offset,k_x2_offset,
     1            k_ex1_offset)
       call reconcilefile(d_ex1,1)
       call get_block(d_ex1,over2,1,0)
       over=over1+over2
c fix this and change header
       call daxpyfile(1,(-1.0d0)*over,x1(j),
     1                x1(ipos),size_x1)
       call daxpyfile(1,(-1.0d0)*over,x2(j),
     1                x2(ipos),size_x2)
      enddo !---j
c normalization of x(imicro)
       over=0.0d0
       over1=0.0d0
       over2=0.0d0
       call tce_zero(d_ex1,size_ex1)
       call c1_c1(x1(ipos),d_ex1,k_x1_offset,k_ex1_offset)
       call reconcilefile(d_ex1,1)
       call get_block(d_ex1,over1,1,0)
       call tce_zero(d_ex1,size_ex1)
       call c2_c2(x2(ipos),d_ex1,k_x2_offset,k_ex1_offset)
       call reconcilefile(d_ex1,1)
       call get_block(d_ex1,over2,1,0)
       over=over1+over2
       over=dsqrt(over)
       call dscalfile(1.0d0/over,x1(ipos),size_x1)
       call dscalfile(1.0d0/over,x2(ipos),size_x2)
c
      enddo ! --- i 
c
      call deletefile(d_ex1)
      if (.not.ma_pop_stack(l_ex1_offset))
     1    call errquit("tce_energy: deall1",24,MA_ERR)
c
      return 
      end 
c
c 
c
c
c
c
c
c
c
c
      subroutine norm_x_vector_act_mr(k_x1_offset,k_x2_offset,
     1                           size_x1,size_x2,totr)
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
      integer iroot
      integer i,j,totr
      integer d_ex1,size_ex1
      integer l_ex1_offset,k_ex1_offset
      integer k_x1_offset,k_x2_offset
      integer size_x1,size_x2
      double precision over,over1,over2
      character*255 filename
c
      irrep_c=irrep_x
      irrep_d=irrep_x
c
      call tce_e_offset(l_ex1_offset,k_ex1_offset,size_ex1)
      call createfile(filename,d_ex1,size_ex1)
      call tce_zero(d_ex1,size_ex1)
      call tce_filename('ex1',filename)
c normalization of initial  space
      do i=1,totr
c
       do j=1,i-1
        over=0.0d0
        over1=0.0d0
        over2=0.0d0
        call tce_zero(d_ex1,size_ex1)
        call c1_d1(x1(i),x1(j),d_ex1,k_x1_offset,k_x1_offset,
     1             k_ex1_offset)
        call reconcilefile(d_ex1,1)
        call get_block(d_ex1,over1,1,0)
        call tce_zero(d_ex1,size_ex1)
        call c2_d2_act(x2(i),x2(j),d_ex1,k_x2_offset,k_x2_offset,
     1             k_ex1_offset)
        call reconcilefile(d_ex1,1)
        call get_block(d_ex1,over2,1,0)
        over=over1+over2
c fix this and change header
        call daxpyfile(1,(-1.0d0)*over,x1(j),
     1                 x1(i),size_x1)
        call daxpyfile(1,(-1.0d0)*over,x2(j),
     1                 x2(i),size_x2)
       enddo ! ---- over j
c
c normalization of x(i)
       over=0.0d0
       over1=0.0d0
       over2=0.0d0
       call tce_zero(d_ex1,size_ex1)
       call c1_c1(x1(i),d_ex1,k_x1_offset,k_ex1_offset)
       call reconcilefile(d_ex1,1)
       call get_block(d_ex1,over1,1,0)
       call tce_zero(d_ex1,size_ex1)
       call c2_c2_act(x2(i),d_ex1,k_x2_offset,k_ex1_offset)
       call reconcilefile(d_ex1,1)
       call get_block(d_ex1,over2,1,0)
       over=over1+over2
       over=dsqrt(over)
       call dscalfile(1.0d0/over,x1(i),size_x1)
       call dscalfile(1.0d0/over,x2(i),size_x2)
c
      enddo !over i ------
c
      call deletefile(d_ex1)
      if (.not.ma_pop_stack(l_ex1_offset))
     1    call errquit("tce_energy: deall1",24,MA_ERR)
c
      return 
      end 
c
c
c
c
c
      subroutine form_heff_act_mr(k_x1_offset,k_x2_offset,
     &                    size_x1,size_x2,
     &                    c_eff,h_eff,omega,imicro,xdiis,iter_eomcc)
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
      integer imicro,xdiis,residual
      integer  i,j,iter_eomcc,ipos,ipos1,ipos2,imin
      integer d_ex1,size_ex1
      integer l_ex1_offset,k_ex1_offset
      integer k_x1_offset,k_x2_offset
      integer size_x1,size_x2
      double precision over,over1,over2
      double precision h_eff(xdiis*nroots_reduced,xdiis*nroots_reduced)
      double precision c_eff(xdiis*nroots_reduced,nroots_reduced) 
      double precision h_copy(xdiis*nroots_reduced,xdiis*nroots_reduced)
      double precision omega(nroots_reduced),xmax,xmin
c
      character*1 jobvl,jobvr
      integer lwork,info_x
      double precision wr(xdiis*nroots_reduced),wi(xdiis*nroots_reduced)
      double precision vl(xdiis*nroots_reduced,xdiis*nroots_reduced)
      double precision vr(xdiis*nroots_reduced,xdiis*nroots_reduced)
      double precision work(4*xdiis*nroots_reduced)
c 
      character*255 filename
c
      logical nodezero         ! True if node 0
c
       nodezero=(ga_nodeid().eq.0)
c
      jobvl='N'
      jobvr='V'
      lwork=4*xdiis*nroots_reduced
c
      irrep_c=irrep_x
      irrep_d=irrep_x
c 
      call tce_e_offset(l_ex1_offset,k_ex1_offset,size_ex1)
      call createfile(filename,d_ex1,size_ex1)
      call tce_zero(d_ex1,size_ex1)
      call tce_filename('ex1',filename)
c
      do i=1,(imicro-1)*nroots_reduced !-------- i-loop
c
       do j=1,nroots_reduced
        over=0.0d0
        over1=0.0d0
        over2=0.0d0
        ipos=(imicro-1)*nroots_reduced+j
c j related to column
       call tce_zero(d_ex1,size_ex1)
       call c1_d1(x1(i),xp1(ipos),d_ex1,k_x1_offset,k_x1_offset,
     1            k_ex1_offset)
       call reconcilefile(d_ex1,1)
       call get_block(d_ex1,over1,1,0)
       call tce_zero(d_ex1,size_ex1)
       call c2_d2_act(x2(i),xp2(ipos),d_ex1,k_x2_offset,k_x2_offset,
     1            k_ex1_offset)
       call reconcilefile(d_ex1,1)
       call get_block(d_ex1,over2,1,0)
       over=over1+over2
c fix this and change header
       h_eff(i,ipos)=over
c normalization of x(imicro)
       over=0.0d0
       over1=0.0d0
       over2=0.0d0
c imicro - column
       call tce_zero(d_ex1,size_ex1)
       call c1_d1(x1(ipos),xp1(i),d_ex1,k_x1_offset,k_x1_offset,
     1            k_ex1_offset)
       call reconcilefile(d_ex1,1)
       call get_block(d_ex1,over1,1,0)
       call tce_zero(d_ex1,size_ex1)
       call c2_d2_act(x2(ipos),xp2(i),d_ex1,k_x2_offset,k_x2_offset,
     1            k_ex1_offset)
       call reconcilefile(d_ex1,1)
       call get_block(d_ex1,over2,1,0)
       over=over1+over2
c fix this and change header
       h_eff(ipos,i)=over
c normalization of x(imicro)

       enddo !-------------------- j-loop
       enddo !-------------------- i-loop
c h_eff(imicro,imicro)
       do i=1,nroots_reduced
       do j=1,nroots_reduced
       ipos1=(imicro-1)*nroots_reduced+i
       ipos2=(imicro-1)*nroots_reduced+j
       over=0.0d0
       over1=0.0d0
       over2=0.0d0
c nxrtials - column
       call tce_zero(d_ex1,size_ex1)
       call c1_d1(x1(ipos1),xp1(ipos2),d_ex1,
     1            k_x1_offset,k_x1_offset,
     1            k_ex1_offset)
       call reconcilefile(d_ex1,1)
       call get_block(d_ex1,over1,1,0)
       call tce_zero(d_ex1,size_ex1)
       call c2_d2_act(x2(ipos1),xp2(ipos2),d_ex1,
     1            k_x2_offset,k_x2_offset,
     1            k_ex1_offset)
       call reconcilefile(d_ex1,1)
       call get_block(d_ex1,over2,1,0)
       over=over1+over2
c fix this and change header
       h_eff(ipos1,ipos2)=over
c
      enddo !j ---------
      enddo !i ---------
c
      call deletefile(d_ex1)
      if (.not.ma_pop_stack(l_ex1_offset))
     1    call errquit("tce_energy: deall3",24,MA_ERR)
c
       do i=1,xdiis*nroots_reduced
       do j=1,xdiis*nroots_reduced
        h_copy(i,j)=h_eff(i,j)
       enddo
       enddo
c big number on diagonal beyond current space
       do i=imicro*nroots_reduced+1,xdiis*nroots_reduced
         h_eff(i,i)=10000000000.0d0
       enddo
c
c
        call util_dgeev(jobvl,jobvr,xdiis*nroots_reduced,h_eff,
     &             xdiis*nroots_reduced,wr,wi,
     &             vl,xdiis*nroots_reduced,vr,
     &             xdiis*nroots_reduced,work,lwork,info_x)
c
c selection
c
       do j=1,nroots_reduced ! --- i
       xmin=10000000001.0d0
       imin=0
       do i=1,xdiis*nroots_reduced
        if((wr(i).lt.xmin).and.(dabs(wi(i)).lt.0.000010d0)) then
         xmin=wr(i)
         imin=i
        end if
       enddo
       if(imin.eq.0) then
        call errquit('Eigen-solver - imaginary roots',xdiis,CALC_ERR)
       end if       
c
       omega(j)=wr(imin)
       wr(imin)=10000000001.0d0

       xmax=0.0d0
       do i=1,xdiis*nroots_reduced
        c_eff(i,j)=vr(i,imin)
        xmax=xmax+vr(i,imin)*vr(i,imin)
       enddo
       do i=1,xdiis*nroots_reduced
        c_eff(i,j)=(1.0d0/dsqrt(xmax))*c_eff(i,j)
       enddo
c
       enddo ! j ---
c
c
       do i=1,xdiis*nroots_reduced
       do j=1,xdiis*nroots_reduced
        h_eff(i,j)=h_copy(i,j)
       enddo
       enddo
c 
      return 
      end 
c
c
c
      subroutine form_correction_act_mr(k_x1_offset,k_x2_offset,
     &                     size_x1,size_x2,
     &                     c_eff,omega,imicro,xdiis,residual,iroot)
c
c $Id$
c
c  delta_x = (omega - H-bar) xc
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
      integer imicro,xdiis
      integer i,j,iroot
      integer d_ex1,size_ex1
      integer l_ex1_offset,k_ex1_offset
      integer k_x1_offset,k_x2_offset
      integer size_x1,size_x2
      double precision over,over1,over2
      double precision c_eff(xdiis*nroots_reduced),omega,residual
      character*255 filename
      logical nodezero         ! True if node 0
c
       nodezero=(ga_nodeid().eq.0)
c
      irrep_c=irrep_x
      irrep_d=irrep_x
c
c
      do i=1,(imicro-1)*nroots_reduced
       call daxpyfile(1,(-1.0d0)*omega*c_eff(i),x1(i),
     1                x1((imicro-1)*nroots_reduced+iroot),size_x1)
       call daxpyfile(1,c_eff(i),xp1(i),
     1                x1((imicro-1)*nroots_reduced+iroot),size_x1) 
       call daxpyfile(1,(-1.0d0)*omega*c_eff(i),x2(i),
     1                x2((imicro-1)*nroots_reduced+iroot),size_x2)
       call daxpyfile(1,c_eff(i),xp2(i),
     1                x2((imicro-1)*nroots_reduced+iroot),size_x2)
c  
c
      enddo
c 
c calculating residual
c
      over1=0.0d0
      over2=0.0d0
c
      call tce_e_offset(l_ex1_offset,k_ex1_offset,size_ex1)
      call createfile(filename,d_ex1,size_ex1)
      call tce_zero(d_ex1,size_ex1)
      call tce_filename('ex1',filename)
c
      call c1_c1(x1((imicro-1)*nroots_reduced+iroot),
     &     d_ex1,k_x1_offset,k_ex1_offset)
      call reconcilefile(d_ex1,1)
      call get_block(d_ex1,over1,1,0)
c
      call tce_zero(d_ex1,size_ex1)
      call c2_c2_act(x2((imicro-1)*nroots_reduced+iroot),
     &     d_ex1,k_x2_offset,k_ex1_offset)
      call reconcilefile(d_ex1,1)
      call get_block(d_ex1,over2,1,0)
      residual=dsqrt(over1+over2)
c
c 
      call deletefile(d_ex1)
      if (.not.ma_pop_stack(l_ex1_offset))
     1    call errquit("tce_energy: deall2",24,MA_ERR)
c
c JACOBI
c
      call tce_jacobi_x1(x1((imicro-1)*nroots_reduced+iroot),
     &     k_x1_offset,omega)
      call tce_jacobi_x2a(x2((imicro-1)*nroots_reduced+iroot),
     &     k_x2_offset,omega)
c
c
c
      return 
      end 
c
c
c
      subroutine form_best_guess_act_mr(k_x1_offset,k_x2_offset,
     &                     size_x1,size_x2,
     &                     c_eff,imicro_x,xdiis,iroot)
c
c $Id$
c
c  xc(iroot) = sum(i=1,...,imicro) x(i)
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
      integer imicro_x,xdiis,residual
      integer i,j
      integer d_ex1,size_ex1
      integer l_ex1_offset,k_ex1_offset
      integer k_x1_offset,k_x2_offset
      integer size_x1,size_x2
      integer iroot
      double precision over,over1,over2
      double precision c_eff(xdiis*nroots_reduced)
      logical nodezero         ! True if node 0
c
       nodezero=(ga_nodeid().eq.0)
c

c
      irrep_c=irrep_x
      irrep_d=irrep_x
c
c
      call tce_zero(xc1(iroot),size_x1)
      call tce_zero(xc2(iroot),size_x2)
c
      do i=1,imicro_x*nroots_reduced
       call daxpyfile(1,c_eff(i),x1(i),
     1                xc1(iroot),size_x1)
       call daxpyfile(1,c_eff(i),x2(i),
     1                xc2(iroot),size_x2)
      enddo
c
c
      return 
      end 
c
c
c
      subroutine deltax_ort_norm_act_mr(k_x1_offset,k_x2_offset,
     1                           size_x1,size_x2,
     1                           imicro,xdiis)
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
      integer imicro,xdiis
      integer  i,j,ipos
      integer d_ex1,size_ex1
      integer l_ex1_offset,k_ex1_offset
      integer k_x1_offset,k_x2_offset
      integer size_x1,size_x2
      double precision over,over1,over2
      character*255 filename
c
      irrep_c=irrep_x
      irrep_d=irrep_x
c
      call tce_e_offset(l_ex1_offset,k_ex1_offset,size_ex1)
      call createfile(filename,d_ex1,size_ex1)
      call tce_zero(d_ex1,size_ex1)
      call tce_filename('ex1',filename)
c
      do i=1,nroots_reduced
      ipos=(imicro-1)*nroots_reduced+i
      do j=1,ipos-1
       over=0.0d0
       over1=0.0d0
       over2=0.0d0
       call tce_zero(d_ex1,size_ex1)
       call c1_d1(x1(ipos),x1(j),d_ex1,k_x1_offset,k_x1_offset,
     1            k_ex1_offset)
       call reconcilefile(d_ex1,1)
       call get_block(d_ex1,over1,1,0)
       call tce_zero(d_ex1,size_ex1)
       call c2_d2_act(x2(ipos),x2(j),d_ex1,k_x2_offset,k_x2_offset,
     1            k_ex1_offset)
       call reconcilefile(d_ex1,1)
       call get_block(d_ex1,over2,1,0)
       over=over1+over2
c fix this and change header
       call daxpyfile(1,(-1.0d0)*over,x1(j),
     1                x1(ipos),size_x1)
       call daxpyfile(1,(-1.0d0)*over,x2(j),
     1                x2(ipos),size_x2)
      enddo !---j
c normalization of x(imicro)
       over=0.0d0
       over1=0.0d0
       over2=0.0d0
       call tce_zero(d_ex1,size_ex1)
       call c1_c1(x1(ipos),d_ex1,k_x1_offset,k_ex1_offset)
       call reconcilefile(d_ex1,1)
       call get_block(d_ex1,over1,1,0)
       call tce_zero(d_ex1,size_ex1)
       call c2_c2_act(x2(ipos),d_ex1,k_x2_offset,k_ex1_offset)
       call reconcilefile(d_ex1,1)
       call get_block(d_ex1,over2,1,0)
       over=over1+over2
       over=dsqrt(over)
       call dscalfile(1.0d0/over,x1(ipos),size_x1)
       call dscalfile(1.0d0/over,x2(ipos),size_x2)
c
      enddo ! --- i 
c
      call deletefile(d_ex1)
      if (.not.ma_pop_stack(l_ex1_offset))
     1    call errquit("tce_energy: deall1",24,MA_ERR)
c
      return 
      end 
c
c 

