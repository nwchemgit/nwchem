c ==> Added by by FA-03-29-12:
c Note.- To be used in shell_fock_build2() located in
c        shell_fock_build.F for cphf3() in aoresponse
c        calculations that include damping.
      subroutine shell_fock_build_nosymm(geom,basis, nder, ndens,
     $      jfac, kfac, tol2e, 
     &      g_dens, g_fock)
c Purpose : Avoid symmetrization of g_fock-FA
c
c Took oactive out for the time being, but will probably need to
c bring it back - tlw
c
c osym=.false. now means that the density matrix is anti-symmetric
c and so, at the end of this subroutine, the output fock matrix
c will be anti-symmetrized - sh, 2/15/02
c
c --> Experimental (not published yet)

      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
c
      integer basis, nder, ndens, nfock
      integer geom
      integer g_dens, g_fock
      double precision jfac(*), kfac(*), tol2e
c      logical oactive(*)        ! [input] T if atom is active
c
      integer nsh, nbf, max_at_bf, max_sh_bf, nshblocks
      integer blen
      integer lend, leri, lscr, max2e
      integer l_shmap, l_shglo, l_shghi, l_bfmap, l_rbfmap,
     $     l_bfglo, l_bfghi, l_labels, l_list, l_q4, l_bftoat
      integer k_shmap, k_shglo, k_shghi, k_bfmap, k_rbfmap,
     $     k_bfglo, k_bfghi, k_labels, k_list, k_q4, k_bftoat
      integer l_shbflo, l_shbfhi, l_eri, l_scr
      integer k_shbflo, k_shbfhi, k_eri, k_scr
      integer i, j, iat
      integer ityp, ndim, dims(3)
c
      integer 
     $     k_d_ij,  k_d_kl,  k_d_ik,  k_d_jl,  k_d_il,  k_d_jk, 
     $     k_f_ij,  k_f_kl,  k_f_ik,  k_f_jl,  k_f_il,  k_f_jk,
     $     l_d_ij,  l_d_kl,  l_d_ik,  l_d_jl,  l_d_il,  l_d_jk, 
     $     l_f_ij,  l_f_kl,  l_f_ik,  l_f_jl,  l_f_il,  l_f_jk
c
      integer maxblen           ! maximum value for blen = max no. bf in group
      integer maxsh             ! maximum no. of shells in a group
      integer maxq              ! max quartets in a request
      logical osym
      parameter (maxblen=36, maxsh=10)
      parameter (maxq=maxsh**4)
      integer ga_create_atom_blocked
      logical xc_gotxc
      external ga_create_atom_blocked,xc_gotxc
c
c Quietly modify the values of jfac and kfac so that the original
c values coming in can be consistent with the original fock_2e
c
      do i= 1, ndens
         jfac(i) = jfac(i)*4.0d00
         kfac(i) = kfac(i)*2.0d00
      enddo
c
c Get the dimensions of g_fock and eventually g_dens.
c
      call nga_inquire(g_fock, ityp, ndim, dims)
      nfock = dims(1)
c
      if (.not. bas_numbf(basis,nbf))
     $     call errquit('new_fock: could not get nbf',0, BASIS_ERR)
      if (.not. bas_numcont(basis,nsh))
     $     call errquit('new_fock: could not get nsh',0, BASIS_ERR)
      if (.not. bas_nbf_ce_max(basis,max_at_bf))
     $     call errquit('new_fock: could not get max_at_bf',0,
     &       BASIS_ERR)
      if (.not. bas_nbf_cn_max(basis,max_sh_bf))
     $     call errquit('new_fock: could not get max_sh_bf',0,
     &       BASIS_ERR)
c
      blen = min(nbf,maxblen,6*max_sh_bf) ! d(6)*6=36, 6**4=1296 quartets
c
      lend = blen*blen*ndens
      if (.not. ma_push_get(mt_dbl, lend, 'd_ij', l_d_ij, k_d_ij)) 
     $     call errquit('new_fock: d_ij', lend, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lend, 'd_kl', l_d_kl, k_d_kl)) 
     $     call errquit('new_fock: d_kl', lend, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lend, 'd_ik', l_d_ik, k_d_ik)) 
     $     call errquit('new_fock: d_ik', lend, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lend, 'd_il', l_d_il, k_d_il)) 
     $     call errquit('new_fock: d_il', lend, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lend, 'd_jk', l_d_jk, k_d_jk)) 
     $     call errquit('new_fock: d_jk', lend, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lend, 'd_jl', l_d_jl, k_d_jl)) 
     $     call errquit('new_fock: d_jl', lend, MA_ERR)
c
      lend = blen*blen*nfock
      if (.not. ma_push_get(mt_dbl, lend, 'f_ij', l_f_ij, k_f_ij)) 
     $     call errquit('new_fock: f_ij', lend, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lend, 'f_kl', l_f_kl, k_f_kl)) 
     $     call errquit('new_fock: f_kl', lend, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lend, 'f_ik', l_f_ik, k_f_ik)) 
     $     call errquit('new_fock: f_ik', lend, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lend, 'f_il', l_f_il, k_f_il)) 
     $     call errquit('new_fock: f_il', lend, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lend, 'f_jk', l_f_jk, k_f_jk)) 
     $     call errquit('new_fock: f_jk', lend, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lend, 'f_jl', l_f_jl, k_f_jl)) 
     $     call errquit('new_fock: f_jl', lend, MA_ERR)
c
      if (.not. ma_push_get(mt_int, nsh, 'shmap', l_shmap, k_shmap))
     $  call errquit('shell_fock_build:could not allocate shmap',nsh,
     &       MA_ERR)
      if (.not. ma_push_get(mt_int, nsh, 'shglo', l_shglo, k_shglo))
     $  call errquit('shell_fock_build:could not allocate blo',nsh,
     &       MA_ERR)
      if (.not. ma_push_get(mt_int, nsh, 'shghi', l_shghi, k_shghi))
     $  call errquit('shell_fock_build:could not allocate bhi',nsh,
     &       MA_ERR)
      if (.not. ma_push_get(mt_int, nsh, 'shbflo', l_shbflo, k_shbflo))
     $  call errquit('shell_fock_build:could not allocate bflo',nsh,
     &       MA_ERR)
      if (.not. ma_push_get(mt_int, nsh, 'shbfhi', l_shbfhi, k_shbfhi))
     $  call errquit('shell_fock_build:could not allocate bfhi',nsh,
     &       MA_ERR)
      if (.not. ma_push_get(mt_int, nbf, 'bfglo', l_bfglo, k_bfglo))
     $  call errquit('shell_fock_build:could not allocate blo',nbf,
     &       MA_ERR)
      if (.not. ma_push_get(mt_int, nbf, 'bfghi', l_bfghi, k_bfghi))
     $  call errquit('shell_fock_build:could not allocate bhi',nbf,
     &       MA_ERR)
      if (.not. ma_push_get(mt_int, nbf, 'bfmap', l_bfmap, k_bfmap))
     $  call errquit('shell_fock_build:could not allocate bfmap',nbf,
     &       MA_ERR)
      if (.not. ma_push_get(mt_int, nbf, 'rbfmap', l_rbfmap, k_rbfmap))
     $  call errquit('shell_fock_build:could not allocate rbfmap',nbf,
     &       MA_ERR)
      if (.not. ma_push_get(mt_int, nbf, 'bftoat', l_bftoat, k_bftoat))
     $  call errquit('shell_fock_build:could not allocate bftoat',nbf,
     &       MA_ERR)
c 
      call grad_shorder(basis, nsh, nbf, maxsh, blen,
     $     nshblocks, int_mb(k_shglo), int_mb(k_shghi),
     &     int_mb(k_shmap),
     $     int_mb(k_bfmap), int_mb(k_rbfmap), int_mb(k_bfglo),
     $     int_mb(k_bfghi), int_mb(k_shbflo), int_mb(k_shbfhi))

      call intb_mem_2e4c(max2e, lscr) ! blocking algorithm
      if (nder.eq.0) then
        leri = max(max2e,1296*100) ! 100 D quartets
      else
        leri = max(max2e*12,1296*100*12) ! 100 D quartets
      endif
c 
      if (.not. ma_push_get(mt_dbl,leri,'eri',l_eri,k_eri)) 
     $     call errquit('new_fock:could not allocate buffer',leri,
     &       MA_ERR)
      if (.not. ma_push_get(mt_dbl,lscr,'scr',l_scr,k_scr))
     $     call errquit('new_fock: scratch alloc failed', lscr,
     &       MA_ERR)
      if (.not. ma_push_get(mt_int,4*leri,'labels',l_labels,k_labels))
     $     call errquit('new_fock: could not allocate labels',4*leri,
     &       MA_ERR)
c     
      do i = 1, nbf
         if (.not. bas_bf2ce(basis, i, iat)) call errquit('bf2ce',i,
     &       MA_ERR)
         int_mb(k_bftoat+i-1) = iat
      end do

      if (.not. ma_push_get(mt_int, 4*maxq, 'list', l_list, k_list))
     $     call errquit('new_fock:could not allocate list',4*maxq,
     &       MA_ERR)
      if (.not. ma_push_get(mt_dbl, maxq, 'q4', l_q4, k_q4))
     $     call errquit('new_fock:could not allocate q4',maxq,
     &       MA_ERR)

      call nga_reorder(g_dens, .true., int_mb(k_rbfmap),
     $        .true., int_mb(k_rbfmap))
      call nga_reorder(g_fock, .true., int_mb(k_rbfmap),
     $        .true., int_mb(k_rbfmap))

      call newfock(
     $     basis, nder, ndens, nfock,
     $     g_dens, g_fock,
     $     dbl_mb(k_d_ij), dbl_mb(k_d_kl), dbl_mb(k_d_ik),  
     $     dbl_mb(k_d_jl), dbl_mb(k_d_il), dbl_mb(k_d_jk), 
     $     dbl_mb(k_f_ij), dbl_mb(k_f_kl), dbl_mb(k_f_ik),  
     $     dbl_mb(k_f_jl), dbl_mb(k_f_il), dbl_mb(k_f_jk), 
     $     lscr, dbl_mb(k_scr), leri, dbl_mb(k_eri), int_mb(k_labels), 
     $     maxq, int_mb(k_list), dbl_mb(k_q4), tol2e,
     $     nsh,  nbf, blen,
     $     nshblocks, int_mb(k_shmap), int_mb(k_shglo), int_mb(k_shghi),
     $     int_mb(k_bfglo), int_mb(k_bfghi), int_mb(k_bfmap), 
     $     int_mb(k_rbfmap),
     $     int_mb(k_bftoat), int_mb(k_shbflo), int_mb(k_shbfhi), 
c     $     oactive,
     $     jfac, kfac)

      call nga_reorder(g_dens, .true., int_mb(k_bfmap),
     $     .true., int_mb(k_bfmap))
      call nga_reorder(g_fock, .true., int_mb(k_bfmap),
     $     .true., int_mb(k_bfmap))
c
c  Clean up memory allocated in this routine
c
      if (.not. ma_chop_stack(l_d_ij))
     $  call errquit('shell_fock_build:failed chopping MA stack',555,
     &       MA_ERR)
c
c Quietly modify the values of jfac and kfac back to the original
c
      do i= 1, ndens
         jfac(i) = jfac(i)/4.0d00
         kfac(i) = kfac(i)/2.0d00
      enddo
c
      if(xc_gotxc().and.nder.eq.0)then
         call xc_newfock(geom,basis,jfac,kfac,nbf,
     %        nder, ndens, nfock, g_fock,g_dens)
      endif

c      if (osym) then
c         call ga_symmetrize(g_fock)
c      else
c         call ga_antisymmetrize(g_fock)
c      endif
c
      return
      end

      subroutine shell_fock_build2(g_fock,     ! out: Fock    matrices
     &                             g_dens,     ! in : density matrices
     &                             geom,       ! in : geom  handle
     &                             basis,      ! in : basis handle
     &                             nbf,        ! in : nr. basis functions
     &                             nvec,       ! in : nr. vecs (x,y,z)
     &                             npol,       ! in : nr. polarizations =1 RDFT =2 UDFT
     &                             ncomp,      ! in : nr. components = 1 or 2
     &                             nblock,     ! in : nr. of g_dens,g_fock blocks
     &                             noskew_uhf, ! in : = .true. if symm dens
     &                             tol2e,      ! in :
     &                             debug)      ! in : =.true. -> debugging printouts
c Purpose: Single call of shell_fock_build when having
c          two sets of densities g_dens(ipm) ipm=1,nblock  nblock=2
c          It works for unrestricted (npol=2) and 
c                         restricted (npol=1) calculations.
c Author : Fredy W. Aquino, Northwestern University (Oct 2012)
c Date   : 03-28-12
c Note1.- It works for combinations: (ncomp,nblock)=(1,[1or2]),(2,[2or4])
c        Examples: In uhf_hessv_2e2_opt() is used for (ncomp,nblock)=(1,2) RE,IM (g_x,g_ax)
c                  In uhf_hessv_2e3       is used for (ncomp,nblock)=(2,4) RE1 RE2 IM1 IM2 (g_x,g_ax)
c                  In rohf_hessv_2e3_opt  is used for (ncomp,nblock)=(2,4) RE1 RE2 IM1 IM2 (g_x,g_ax)
c Note2.- uhf_XX  are routines for unrestricted DFT calculation.
c         rohf_XX are routines for   restricted DFT calculation (even rohf means something else).           
c --> Experimental (not published yet)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
#include "bas.fh"
#include "bgj.fh"
#include "stdio.fh"
#include "util.fh"
#include "case.fh"
#include "msgids.fh"
c     case I: symmetric density matrix
c     case II: antisymmetric density matrix. note that jfac for this
c     case is zero, i.e. there is no Coulomb term because there is
c     no density. However, there are still exchange contributions from
c     HF exchange. For pure DFT we just skip the call to save CPU time.
      integer geom,basis,
     &        nbf,
     &        npol,imul,nmul,ncomp,nblock,nvec,nder,
     &        ndensity2,nfock2,
     &        shift1,shift2,shift3,
     &        nn,iset,ipm,xyz,
     &        alo(3),ahi(3),
     &        blo(3),bhi(3),
     &        dims(3),chunk(3),
     &        l_jfac,k_jfac,
     &        l_kfac,k_kfac
      logical noskew_uhf,debug
      integer g_dens(*),g_fock(*),
     &        g_dens2,g_fock2
      double precision tol2e,
     &                 jfac(4),kfac(4),
     &                 xc_xfac
      data nder /0/
      external uhf_fock_setup,
     &         shell_fock_build_nosymm,
     &         shell_fock_build_cam_nosymm
      if (.not.((ncomp.eq.1 .and. nblock.eq.1) .or. 
     &          (ncomp.eq.1 .and. nblock.eq.2) .or.
     &          (ncomp.eq.2 .and. nblock.eq.2) .or. 
     &          (ncomp.eq.2 .and. nblock.eq.4))) then
       call errquit(
     &  'sh-fockbld2: (ncomp,nblock) ne (1,[1or2]) or (2,[2or4])',
     &               nblock, MA_ERR)
       stop
      endif
      nmul=1
      if (npol.eq. 2) nmul=2
      nn=nmul*nblock*npol*nvec
c --- FA-12-02-13 -------------- START
c      if (ga_nodeid().eq.0) then
c       write(*,10) npol,nmul,nblock,nvec
c 10    format('(npol,nmul,nblock,nvec)=(',
c     &        i2,',',i2,',',i2,',',i2,')')
c      endif
c --- FA-12-02-13 -------------- END
      if (.not.MA_Push_Get(mt_dbl,nn,'hessv jfacs',l_jfac,k_jfac))
     &     call errquit('sh-fockbld2: cannot allocate jfac',
     &                  nn, MA_ERR)
      if (.not.MA_Push_Get(mt_dbl,nn,'hessv kfacs',l_kfac,k_kfac))
     &     call errquit('sh-fockbld2: cannot allocate kfac',
     &                  nn, MA_ERR)
c     initialize jfacs and kfacs (will go into the uhf_fock_setup)
c     If DFT get part of the exact exchange defined
c      xc_xfac = 1.0d0
c      if (use_theory.eq.'dft') xc_xfac = bgj_kfac() 
      xc_xfac  = bgj_kfac()       
      dims(1)  = nblock*npol*nvec
      dims(2)  = nbf
      dims(3)  = nbf
      chunk(1) = dims(1)
      chunk(2) = -1
      chunk(3) = -1
      if (.not. nga_create (MT_DBL,3,dims,'CPKS dens',chunk,
     &     g_dens2)) call errquit
     &    ('sh-fockbld2: could not allocate g_dens2',555, GA_ERR)
      call ga_zero(g_dens2)   
      if (.not. nga_create (MT_DBL, 3, dims, 'Fockv',chunk,
     &     g_fock2)) call errquit
     &    ('sh-fockbld2: could not allocate g_fock2',555, GA_ERR)
      call ga_zero(g_fock2)
      ndensity2 = nblock*npol*nvec
      nfock2    = ndensity2 
      blo(1) = 0
      bhi(1) = 0
      blo(2) = 1
      bhi(2) = nbf
      blo(3) = 1
      bhi(3) = nbf
      alo(2) = 1
      ahi(2) = nbf
      alo(3) = 1
      ahi(3) = nbf
      if (noskew_uhf) then
        jfac(1) =  1.0d0
        jfac(2) =  0.0d0
        jfac(3) =  1.0d0
        jfac(4) =  0.0d0
      else
        jfac(1) =  0.0d0
        jfac(2) =  0.0d0
        jfac(3) =  0.0d0
        jfac(4) =  0.0d0
      endif
           if (npol.eq.1) then
        kfac(1)=-0.5d0*xc_xfac
        kfac(2)=-0.5d0*xc_xfac
        kfac(3)=-0.5d0*xc_xfac
        kfac(4)=-0.5d0*xc_xfac
      else if (npol.eq.2) then
        kfac(1)= 1.0d0*xc_xfac
        kfac(2)= 1.0d0*xc_xfac
        kfac(3)= 1.0d0*xc_xfac
        kfac(4)= 1.0d0*xc_xfac
      endif ! end-if-npol
      do iset=1,npol  ! nr. of polarizations (alpha,beta)
       do ipm=1,nblock ! nr. of components
        do xyz=1,nvec ! nr. of directions (x,y,z)
          alo(1) = nvec*(iset-1)+xyz
          ahi(1) = alo(1)
          blo(1) = blo(1) + 1
          bhi(1) = blo(1)
          if (debug) then
           if (ga_nodeid().eq.0) then
            write(*,1) iset,ipm,
     &                alo(1),ahi(1),
     &                alo(2),ahi(2),
     &                alo(3),ahi(3),
     &                blo(1),bhi(1),
     &                blo(2),bhi(2),
     &                blo(3),bhi(3)
 1          format('(iset,ipm)=(',i3,',',i3,')',
     &           ' (alo-ahi)=(',i3,',',i3,',',
     &                          i3,',',i3,',',
     &                          i3,',',i3,')',
     &           ' (blo-bhi)=(',i3,',',i3,',',
     &                          i3,',',i3,',',
     &                          i3,',',i3,')')
           endif
          endif ! end-if-debug
        dbl_mb(k_jfac+blo(1)-1)=jfac(ipm)
        dbl_mb(k_kfac+blo(1)-1)=kfac(ipm)
c
c --- FA-12-02-13 -------------- START
c        if (ga_nodeid().eq.0) then
c         write(*,11) blo(1)-1,iset,ipm,xyz,
c     &               kfac(ipm),jfac(ipm)
c 11      format('(blo,iset,ipm,xyz,kfac,jfac)=(',
c     &          i2,',',i2,',',i2,',',i2,',',
c     &          f14.12,',',f8.6,')')
c        endif
c --- FA-12-02-13 -------------- END

        call nga_copy_patch('n',g_dens(ipm),alo,ahi, ! g_dens--> g_dens2
     &                          g_dens2    ,blo,bhi)
        enddo ! end-loop-xyz
       enddo ! end-loop-ipm
      enddo ! end-loop-iset

c        stop  ! FA-12-02-13

      if (debug) then
       if (ga_nodeid().eq.0)
     &    write(*,*) '---- g_dens2-nw-BEF-uhf_fock_setup --- START'
       call ga_print(g_dens2)
       if (ga_nodeid().eq.0)
     &    write(*,*) '---- g_dens2-nw-BEF-uhf_fock_setup --- END'
      endif ! end-if-debug
      if (npol.eq.2) then
       if (debug) then
        if (ga_nodeid().eq.0)
     &   write(*,*) 'BEF uhf_fock_setup4 ndensity2=',ndensity2
       endif ! end-if-debug
       call uhf_fock_setup(g_dens2,          ! duplicated size in output
     &                     g_fock2,          ! duplicated size in output
     &                     dbl_mb(k_jfac),   ! duplicated size in output
     &                     dbl_mb(k_kfac),   ! duplicated size in output
     &                     ndensity2,nfock2, ! duplicated size in output
     &                     nbf,.false.)
       if (debug) then
        if (ga_nodeid().eq.0)
     &   write(*,*) 'AFT uhf_fock_setup4 ndensity2=',ndensity2
       endif ! end-if-debug
      endif ! end-if-npol

      if (debug) then
       if (ga_nodeid().eq.0)
     &    write(*,*) '--- g_dens2-nw-AFT-uhf_fock_setup---- START'
       call ga_print(g_dens2)
       if (ga_nodeid().eq.0)
     &    write(*,*) '--- g_dens2-nw-AFT-uhf_fock_setup---- END'
       do iset = 0, nmul*npol*nblock*nvec-1
        if (ga_nodeid().eq.0) then
         write(*,3) iset,
     &              dbl_mb(k_jfac+iset),
     &              dbl_mb(k_kfac+iset)
 3       format('AFT-uhf_fock_setup-2(iset,j,k)=(',
     &          i3,',',f12.6,',',f12.6,')')
        endif
       enddo ! end-loop-iset
      endif ! end-if-debug

       if (.not.cam_exch) then
        call shell_fock_build_nosymm(
     &                           geom,basis,nder,nfock2,
     $                           dbl_mb(k_jfac),
     &                           dbl_mb(k_kfac),
     $                           tol2e,
     &                           g_dens2,  ! in: electron density
     &                           g_fock2)  ! ou: fock matrix ipm=1
       else
        call shell_fock_build_cam_nosymm(
     &                        geom,basis,nder,nfock2,
     $                        dbl_mb(k_jfac), 
     &                        dbl_mb(k_kfac),
     $                        tol2e, 
     &                        g_dens2, ! in: electron density
     &                        g_fock2) ! ou: fock matrix ipm=1
       end if ! cam_exch

       if (debug) then
        if (ga_nodeid().eq.0)
     &   write(*,*) '---------g_fock2-out-shellfockbld----- START'
        call ga_print(g_fock2)
        if (ga_nodeid().eq.0)
     &   write(*,*) '---------g_fock2-out-shellfockbld----- END'
       endif ! end-if-debug

c ===== Move g_fock2 --> g_fock(ipm)  ipm=1,nblock == START
      blo(1) = 0
      bhi(1) = 0
      blo(2) = 1
      bhi(2) = nbf
      blo(3) = 1
      bhi(3) = nbf
      alo(2) = 1
      ahi(2) = nbf
      alo(3) = 1
      ahi(3) = nbf
      do iset=1,npol  ! nr. of polarizations (alpha,beta)
       shift1=nmul*nblock*nvec*(iset-1)
       do ipm=1,nblock ! nr. of components nblock=2 for RE or 4 for RE-IM
        shift3=nvec*(ipm-1)
        do imul=1,nmul
         shift2=nblock*nvec*(imul-1)
         do xyz=1,nvec ! nr. of directions J(x,y,z) K(x,y,z) npol=2
                      !                   JK(x,y,z)          npol=1
          alo(1) = nmul*nvec*(iset-1)+nvec*(imul-1)+xyz
          ahi(1) = alo(1)
          blo(1) = shift1+shift2+shift3+xyz
          bhi(1) = blo(1)

          if (debug) then
          if (ga_nodeid().eq.0) then
           write(*,2) iset,ipm,xyz,
     &                alo(1),ahi(1),
     &                alo(2),ahi(2),
     &                alo(3),ahi(3),
     &                blo(1),bhi(1),
     &                blo(2),bhi(2),
     &                blo(3),bhi(3)
 2        format('((iset,ipm,xyz)=(',i3,',',i3,',',i3,')',
     &           ' (alo-ahi)=(',i3,',',i3,',',
     &                          i3,',',i3,',',
     &                          i3,',',i3,')',
     &           ' (blo-bhi)=(',i3,',',i3,',',
     &                          i3,',',i3,',',
     &                          i3,',',i3,')')
          endif
          endif ! end-if-debug

          call nga_copy_patch('n',g_fock2    ,blo,bhi,
     &                            g_fock(ipm),alo,ahi)
         enddo ! end-loop-xyz
        enddo ! end-loop-imul
       enddo ! end-loop-ipm
      enddo ! end-loop-iset  

c ===== Move g_fock2 --> g_fock(ipm)  ipm=1,nblock == END

      if (debug) then
       do ipm=1,nblock
        if (ga_nodeid().eq.0)
     &    write(*,*) '----- g_fock-nw(',ipm,') -------- START'
        call ga_print(g_fock(ipm))
        if (ga_nodeid().eq.0)
     &    write(*,*) '----- g_fock-nw(',ipm,') -------- END'
       enddo ! end-loop-ipm
      endif ! end-if-debug

      if      (ncomp.eq.2) then
        if      (nblock.eq.2) then
         call ga_symmetrize(g_fock(1))
         call ga_antisymmetrize(g_fock(2))
        else if (nblock.eq.4) then  
         call ga_symmetrize(g_fock(1))
         call ga_antisymmetrize(g_fock(2))
         call ga_symmetrize(g_fock(3))
         call ga_antisymmetrize(g_fock(4)) 
        else
         call errquit('sh-fockbld2: nblock ne 2 or 4',555,MA_ERR)
        endif
      else if (ncomp.eq.1) then
        if      (nblock.eq.1) then
         call ga_antisymmetrize(g_fock(1)) ! assumption could be other
                                           ! then I have to fix it
        else if (nblock.eq.2) then  
         call ga_antisymmetrize(g_fock(1)) ! it works for uhf_hessv_2e2_opt
         call ga_antisymmetrize(g_fock(2)) ! that will mimic uhf_hessv_2e2
        else
         call errquit('sh-fockbld2: nblock ne 1 or 2',555,MA_ERR)
        endif
      endif ! end-if-ncomp

      if (debug) then
       do ipm=1,nblock
       if (ga_nodeid().eq.0)
     &    write(*,*) '----- g_fock-nw-symm(',ipm,') -------- START'
       call ga_print(g_fock(ipm))
       if (ga_nodeid().eq.0)
     &    write(*,*) '----- g_fock-nw-symm(',ipm,') -------- END'
       enddo ! end-loop-ipm
      endif ! end-if-debug

      if (.not.ma_pop_stack(l_kfac))
     $  call errquit('sh-fockbld2: pop problem with l_jfac',555,
     &       MA_ERR)
      if (.not.ma_pop_stack(l_jfac))
     $  call errquit('sh-fockbld2: pop problem with l_kfac',555,
     &       MA_ERR)
        if (.not.ga_destroy(g_dens2)) call 
     &   errquit('sh-fockbld2: ga_dstry failed g_dens2',0,GA_ERR)
        if (.not.ga_destroy(g_fock2)) call 
     &   errquit('sh-fockbld2: ga_dstry failed g_dens2',0,GA_ERR)
      return
      end
c $Id$
