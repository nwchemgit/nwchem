cedo#define NBACC 1
cedo#define NBGET 1
      subroutine new_giao_2e(geom, !  in: geometry handle
     &                      basis, !  in: basis handle
     &                        nbf, !  in: nr. of basis functions
     &                      tol2e,
     &                     g_dens, !  in: e-density
     &                     g_fock, ! out: fock-matrix
     &                       kfac, !  in: exchange factor =1 HF-xch =0.2 B3LYP 
     &                       npol) !  in: nr. of polarizations
c $Id$
c                       2011-06-17 20:00     FA
      implicit none
#include "nwc_const.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "hnd_giao.fh"
#include "hnd_pointers.fh"
#include "global.fh"
#include "case.fh"
c
      integer npol,ntot
      integer geom, basis, g_fock, nbf
      integer g_dens(3)  ! g_dens(npol)
      integer g_fock1,disp
      integer i,g_dens1(4) ! for udft calc
      integer plo(3),phi(3),
     &        qlo(3),qhi(3)
      double precision tol2e,kfac,Exc(2)
      integer nExc
      integer nshell, maxang, blen, lend, leri, ii
      integer l_eri,k_eri,l_scr,k_scr,alo(4),ahi(4),g_xc
      integer       k_d_kl,  k_d_ik,  k_d_jl,  k_d_il,  k_d_jk, 
     $     k_f_ij,  k_f_ji,  k_f_ik,  k_f_jl,  k_f_il,  k_f_jk,
     $              l_d_kl,  l_d_ik,  l_d_jl,  l_d_il,  l_d_jk, 
     $     l_f_ij,  l_f_ji,  l_f_ik,  l_f_jl,  l_f_il,  l_f_jk
      integer g_loc
      logical xc_gotxc,xc_hassrhf
      external xc_gotxc,xc_hassrhf
      double precision jloc,kloc,factor_xch
      double precision factor_cam,cam_alpha1,cam_beta1
      integer acc,ind,ispin,nfock(2),coeff(2)
      integer debug_giaofock
      data coeff /1,3/
      data nfock/3,12/ ! For Unrest/Rest calc.

      debug_giaofock=0  ! =1 for debugging giaofock

      if(.not.bas_high_angular(basis,maxang))
     &         call errquit('new_giao: angm error',maxang, BASIS_ERR)
      if (.not. bas_numcont(basis,nshell))
     $     call errquit('new_giao: could not get nsh',0, BASIS_ERR)
      blen = (maxang+1)*(maxang+2)/2
      lend = blen*blen
      if (.not. ma_push_get(mt_dbl, lend, 'd_kl', l_d_kl, k_d_kl)) 
     $     call errquit('new_giao: d_kl', lend, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lend, 'd_ik', l_d_ik, k_d_ik)) 
     $     call errquit('new_giao: d_ik', lend, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lend, 'd_il', l_d_il, k_d_il)) 
     $     call errquit('new_giao: d_il', lend, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lend, 'd_jk', l_d_jk, k_d_jk)) 
     $     call errquit('new_giao: d_jk', lend, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lend, 'd_jl', l_d_jl, k_d_jl)) 
     $     call errquit('new_giao: d_jl', lend, MA_ERR)
      lend = blen*blen*3
      if (.not. ma_push_get(mt_dbl, lend, 'f_ij', l_f_ij, k_f_ij)) 
     $     call errquit('new_giao: f_ij', lend, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lend, 'f_ji', l_f_ji, k_f_ji)) 
     $     call errquit('new_giao: f_ji', lend, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lend, 'f_ik', l_f_ik, k_f_ik)) 
     $     call errquit('new_giao: f_ik', lend, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lend, 'f_il', l_f_il, k_f_il)) 
     $     call errquit('new_giao: f_il', lend, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lend, 'f_jk', l_f_jk, k_f_jk)) 
     $     call errquit('new_giao: f_jk', lend, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lend, 'f_jl', l_f_jl, k_f_jl)) 
     $     call errquit('new_giao: f_jl', lend, MA_ERR)
      leri = ngint    ! times 6 blocks of integrals for gint(i*j*k*l) length
      if (.not. ma_push_get(mt_dbl,leri*6,'eri',l_eri,k_eri)) 
     $     call errquit('new_giao:could not allocate buffer',leri,
     &       MA_ERR)
      if (.not. ma_push_get(mt_dbl,iscrln,'scr',l_scr,k_scr))
     $     call errquit('new_giao: scratch alloc failed', iscrln,
     &       MA_ERR)

c ------ create g_fock1 ---------- START
      plo(1) = -1 
      plo(2) = -1
      plo(3) =  1
      phi(1) = nbf
      phi(2) = nbf
      phi(3) =  3
      if (.not.nga_create(MT_DBL,3,phi,'Fock matrix',plo,g_fock1)) call 
     &    errquit('hnd_giaox: nga_create failed g_fock',0,GA_ERR)
      call ga_zero(g_fock1)
c     == local array ==
      if (.not. ga_duplicate (g_fock1, g_loc,'g_loc'))
     &  call errquit('new_giao_2e: duplicate failed g_loc',555, GA_ERR)
      call ga_zero(g_loc)
c ------ create g_fock1 ---------- END    
c
c     == take care of attenuated or non-attenuated 2e ==
      call ga_zero(g_fock)
      if      (npol.eq.1) then
       ntot=1 
       factor_cam =1.0d0 
       factor_xch =1.0d0  
      else if (npol.eq.2) then
       ntot=3
       factor_cam =2.0d0 
       factor_xch =2.0d0  
      else
       write(*,*) 'Error in new_giao_2e:: npol=1 or 2 npol=',npol
      endif
c
      do ispin=1,npol
        call ga_zero(g_fock1) ! reset g_fock1
c
      if (.not.cam_exch) then ! normal calculations
        call ga_zero(g_loc)
        call giaofock(basis,geom,g_dens(ntot),g_loc,
     $     dbl_mb(k_d_kl), dbl_mb(k_d_ik), dbl_mb(k_d_jl), 
     $     dbl_mb(k_d_il), dbl_mb(k_d_jk), dbl_mb(k_f_ij), 
     $     dbl_mb(k_f_ji), dbl_mb(k_f_ik), dbl_mb(k_f_jl), 
     $     dbl_mb(k_f_il), dbl_mb(k_f_jk), iscrln, dbl_mb(k_scr), 
     $     ngint, dbl_mb(k_eri), tol2e, nshell, blen, 1.d0, 0.d0) ! get full coulomb
        call ga_add(1.d0, g_loc, 1.0d0, g_fock1, g_fock1)
        call ga_zero(g_loc)
        call giaofock(basis,geom,g_dens(ispin),g_loc,
     $     dbl_mb(k_d_kl), dbl_mb(k_d_ik), dbl_mb(k_d_jl), 
     $     dbl_mb(k_d_il), dbl_mb(k_d_jk), dbl_mb(k_f_ij), 
     $     dbl_mb(k_f_ji), dbl_mb(k_f_ik), dbl_mb(k_f_jl), 
     $     dbl_mb(k_f_il), dbl_mb(k_f_jk), iscrln, dbl_mb(k_scr), 
     $     ngint, dbl_mb(k_eri), tol2e, nshell, blen, 0.d0, kfac) ! get full exchange           
        call ga_add(factor_xch, g_loc, 1.0d0, g_fock1, g_fock1)
      else
        call ga_zero(g_loc)
        call case_setflags(.false.) ! turn off attenuation
        call giaofock(basis,geom,g_dens(ntot),g_loc,
     $     dbl_mb(k_d_kl), dbl_mb(k_d_ik), dbl_mb(k_d_jl), 
     $     dbl_mb(k_d_il), dbl_mb(k_d_jk), dbl_mb(k_f_ij), 
     $     dbl_mb(k_f_ji), dbl_mb(k_f_ik), dbl_mb(k_f_jl), 
     $     dbl_mb(k_f_il), dbl_mb(k_f_jk), iscrln, dbl_mb(k_scr), 
     $     ngint, dbl_mb(k_eri), tol2e, nshell, blen, 1.d0, 0.d0)  ! get full coulomb
        call ga_add(1.d0, g_loc, 1.0d0, g_fock1, g_fock1)  ! full Coulomb
c
        call ga_zero(g_loc)
        call case_setflags(.false.) ! turn off attenuation
        cam_alpha1=factor_cam*cam_alpha
        cam_beta1 =factor_cam*cam_beta
        call giaofock(basis,geom,g_dens(ispin),g_loc,
     $     dbl_mb(k_d_kl), dbl_mb(k_d_ik), dbl_mb(k_d_jl), 
     $     dbl_mb(k_d_il), dbl_mb(k_d_jk), dbl_mb(k_f_ij), 
     $     dbl_mb(k_f_ji), dbl_mb(k_f_ik), dbl_mb(k_f_jl), 
     $     dbl_mb(k_f_il), dbl_mb(k_f_jk), iscrln, dbl_mb(k_scr), 
     $     ngint, dbl_mb(k_eri), tol2e, nshell, blen, 0.d0, kfac) ! get full exchange
c
        if (.not.xc_hassrhf()) then
           call ga_add(cam_alpha1, g_loc, 1.0d0, g_fock1, g_fock1)   ! cam_alpha*Exchange
        else
           call ga_add(1.d0, g_loc, 1.0d0, g_fock1, g_fock1)   ! 1.d0*Exchange
        end if
c
        call ga_zero(g_loc)
        call case_setflags(.true.) ! set attenuation
        call giaofock(basis,geom,g_dens(ispin),g_loc,
     $     dbl_mb(k_d_kl), dbl_mb(k_d_ik), dbl_mb(k_d_jl), 
     $     dbl_mb(k_d_il), dbl_mb(k_d_jk), dbl_mb(k_f_ij), 
     $     dbl_mb(k_f_ji), dbl_mb(k_f_ik), dbl_mb(k_f_jl), 
     $     dbl_mb(k_f_il), dbl_mb(k_f_jk), iscrln, dbl_mb(k_scr), 
     $     ngint, dbl_mb(k_eri), tol2e, nshell, blen, 0.d0, kfac) ! get attenuated exchange
c
        if (.not.xc_hassrhf()) then
           call ga_add(cam_beta1, g_loc, 1.0d0, g_fock1, g_fock1)    ! cam_beta*att_Exchange
        else
           call ga_add(-1.d0, g_loc, 1.0d0, g_fock1, g_fock1)    ! -1.d0*att_Exchange
        end if
c
        call case_setflags(.false.) ! turn off attenuation
      end if  ! cam_exch

        if (debug_giaofock.eq.1) then
         if (ga_nodeid().eq.0) then
          write(*,31) ispin
 31       format('---- giao_fock: g_fock1(',i3,') -------- START') 
         endif
         call ga_print(g_fock1)
         if (ga_nodeid().eq.0) then
          write(*,32) ispin
 32       format('---- giao_fock: g_fock1(',i3,') -------- END') 
         endif
        endif
c ------ copying g_fock1 --> g_fock ----- START
       plo(1) = 1
       phi(1) = nbf
       plo(2) = 1
       phi(2) = nbf
       plo(3) = 1
       phi(3) = 3
       disp=3*(ispin-1) 
       qlo(1) = 1
       qhi(1) = nbf
       qlo(2) = 1
       qhi(2) = nbf
       qlo(3) = disp+1
       qhi(3) = disp+3   
       call nga_copy_patch('n',g_fock1,plo,phi,
     &                          g_fock,qlo,qhi) 
c ------ copying g_fock1 --> g_fock ----- END
      enddo ! end-loop-ispin

       if (.not.ga_destroy(g_fock1)) call 
     &     errquit('new_giao_2e: ga_destroy failed g_fock1',
     &              0,GA_ERR)
c
c     == Clean up memory allocated in this routine ==
      if (.not. ma_chop_stack(l_d_kl))
     $  call errquit('new_giao_2e:failed chopping MA stack',555,
     &       MA_ERR)
c
c     Add DFT XC if running DFT calculation

      if (xc_gotxc()) then
c         if (ga_nodeid().eq.0)
c     &     write(*,*) 'Entering for xc, kfac=',kfac
         alo(1) = 1
         ahi(1) = nfock(npol) ! =3 npol=1, 12 npol=2
         alo(2) = 1
         ahi(2) = nbf
         alo(3) = 1
         ahi(3) = nbf
         if (.not.nga_create(MT_DBL,3,ahi,'xc matrix',alo,g_xc)) 
     &      call errquit('hnd_giaox: nga_create failed g_xc',0,GA_ERR)

c        if (ga_nodeid().eq.0)
c    &     write(*,*) 'check npol=',npol    
     
         if      (npol.eq.1) then
          do i=1,npol
           if (.not. ga_create(mt_dbl,nbf,nbf,
     &        'gCSSR: g_dia1',0,0,g_dens1(i)))
     $     call errquit('gCSSR: g_dia1',0,GA_ERR)
           call ga_zero(g_dens1(i))
           call ga_copy(g_dens(i),g_dens1(i)) 
          enddo
         else if (npol.eq.2) then
          do i=1,4
           if (.not. ga_create(mt_dbl,nbf,nbf,
     &        'gCSSR: g_dia1',0,0,g_dens1(i)))
     $     call errquit('gCSSR: g_dia1',0,GA_ERR)
          enddo
          ind=1
          do i=1,npol
           call ga_zero(g_dens1(ind))
           call ga_copy(g_dens(i),g_dens1(ind))            
           call ga_zero(g_dens1(ind+1))
           call ga_copy(g_dens(i),g_dens1(ind+1)) 
           ind=ind+2
          enddo
c -------- check input g_dens1 -------- START        
        if (debug_giaofock.eq.1) then
          do i=1,4
           if (ga_nodeid().eq.0)
     &      write(*,*) '---- gdens1(',i,') -------- START'
            call ga_print(g_dens1(i))
           if (ga_nodeid().eq.0)
     &      write(*,*) '---- gdens1(',i,') -------- END'
          enddo
        endif
c -------- check input g_dens1 -------- END
         endif
         call ga_zero(g_xc)
         Exc(1) = 0.0d0
         Exc(2) = 0.0d0
         nExc = 1
         call fock_xc(geom,nbf,basis,nfock(npol),
     &                  g_dens1,g_xc,
     &                  Exc,nExc,.false.) ! out: g_xc
c
        if (debug_giaofock.eq.1) then
          if (ga_nodeid().eq.0)
     &    write(*,*) '---- giao_fock: g_xc -------- START'
          call ga_print(g_xc)
          if (ga_nodeid().eq.0)
     &    write(*,*) '---- giao_fock: g_xc -------- END'
        endif

          if      (npol.eq.1) then
            do ii = 1, 3
             alo(1) = ii
             ahi(1) = ii
             alo(4) = ii
             ahi(4) = ii
             call nga_add_patch(-2.0d0,g_xc,alo,ahi,
     &                            1.d0,g_fock,alo(2),ahi(2),
     &                                 g_fock,alo(2),ahi(2))
            enddo ! end-loop-ii
          else if (npol.eq.2) then
           acc=1
           do ispin=1,npol
            ind=coeff(ispin)*3
            do ii = 1,3
             alo(1) = ind+ii
             ahi(1) = ind+ii
             alo(4) = acc
             ahi(4) = acc
             call nga_add_patch(+2.0d0,g_xc   ,alo,ahi,
     &                           1.0d0,g_fock,alo(2),ahi(2),
     &                                 g_fock,alo(2),ahi(2))
             acc=acc+1
            enddo ! end-loop-ii
           enddo ! end-loop-ispin
          endif

        if (debug_giaofock.eq.1) then
          if (ga_nodeid().eq.0)
     &    write(*,*) '---- giao_fock: g_fock-updated ----- START'
          call ga_print(g_fock)
          if (ga_nodeid().eq.0)
     &    write(*,*) '---- giao_fock: g_fock-updated -------- END'
        endif

         if (.not.ga_destroy(g_xc)) call
     &      errquit('new_giao_2e: ga_destroy failed g_xc',0,GA_ERR)
         if (npol.eq.1) then
          if (.not.ga_destroy(g_dens1(1))) call
     &       errquit('giao_aotomo: ga_destroy failed g_dens1',0,GA_ERR)
         else
          do i=1,4
          if (.not.ga_destroy(g_dens1(i))) call
     &       errquit('giao_aotomo: ga_destroy failed g_dens1',0,GA_ERR)
          enddo
         endif
      endif  ! xc_gotxc()
c
c     == clean up local array ==
      if (.not.ga_destroy(g_loc))
     $  call errquit('new_giao_2e: destroy failed g_loc',555, GA_ERR)
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c To produce output: g_fock_Coul,g_fock_Exch =========== START
      subroutine new_giao_2e_JK(
     &                       geom, !  in : geometry handle
     &                      basis, !  in : basis handle
     &                        nbf, !  in : nr. of basis functions
     &                      tol2e,
     &                     g_dens, !  in : e-density
     &                     g_fock, ! out : fock-matrix
     &                g_fock_Coul, ! out : fock-matrix with Coulomb  ONLY
     &                g_fock_Exch, ! out : fock-matrix with Exchange ONLY
     &                       kfac, !  in : exchange factor =1 HF-xch =0.2 B3LYP 
     &                       npol) !  in : nr. of polarizations
c $Id$
      implicit none
#include "nwc_const.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "hnd_giao.fh"
#include "hnd_pointers.fh"
#include "global.fh"
#include "case.fh"

      integer npol,ntot
      integer geom, basis, g_fock, nbf
      integer g_fock_Coul,g_fock_Exch
      integer g_dens(3)  ! g_dens(npol)
      integer g_fock1,disp
      integer i,g_dens1(4) ! for udft calc
      integer plo(3),phi(3),
     &        qlo(3),qhi(3)
      double precision tol2e,kfac,Exc(2)
      integer nExc
      integer nshell, maxang, blen, lend, leri, ii
      integer l_eri,k_eri,l_scr,k_scr,alo(4),ahi(4),g_xc
      integer       k_d_kl,  k_d_ik,  k_d_jl,  k_d_il,  k_d_jk, 
     $     k_f_ij,  k_f_ji,  k_f_ik,  k_f_jl,  k_f_il,  k_f_jk,
     $              l_d_kl,  l_d_ik,  l_d_jl,  l_d_il,  l_d_jk, 
     $     l_f_ij,  l_f_ji,  l_f_ik,  l_f_jl,  l_f_il,  l_f_jk
      integer g_loc
c
      logical xc_gotxc, xc_hassrhf
      external xc_gotxc, xc_hassrhf
c
      double precision jloc,kloc,factor_xch,
     &       cam_alpha1,cam_beta1,factor_cam
      integer acc,ind,ispin,nfock(2),coeff(2)
      integer debug_giaofock
      data coeff /1,3/
      data nfock/3,12/ ! For Unrest/Rest calc.

      debug_giaofock=0 ! =1 for debugging giaofock

      if(.not.bas_high_angular(basis,maxang))
     &         call errquit('new_giao: angm error',maxang, BASIS_ERR)
      if (.not. bas_numcont(basis,nshell))
     $     call errquit('new_giao: could not get nsh',0, BASIS_ERR)

      blen = (maxang+1)*(maxang+2)/2
      lend = blen*blen

      if (.not. ma_push_get(mt_dbl, lend, 'd_kl', l_d_kl, k_d_kl)) 
     $     call errquit('new_giao: d_kl', lend, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lend, 'd_ik', l_d_ik, k_d_ik)) 
     $     call errquit('new_giao: d_ik', lend, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lend, 'd_il', l_d_il, k_d_il)) 
     $     call errquit('new_giao: d_il', lend, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lend, 'd_jk', l_d_jk, k_d_jk)) 
     $     call errquit('new_giao: d_jk', lend, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lend, 'd_jl', l_d_jl, k_d_jl)) 
     $     call errquit('new_giao: d_jl', lend, MA_ERR)

      lend = blen*blen*3
      if (.not. ma_push_get(mt_dbl, lend, 'f_ij', l_f_ij, k_f_ij)) 
     $     call errquit('new_giao: f_ij', lend, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lend, 'f_ji', l_f_ji, k_f_ji)) 
     $     call errquit('new_giao: f_ji', lend, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lend, 'f_ik', l_f_ik, k_f_ik)) 
     $     call errquit('new_giao: f_ik', lend, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lend, 'f_il', l_f_il, k_f_il)) 
     $     call errquit('new_giao: f_il', lend, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lend, 'f_jk', l_f_jk, k_f_jk)) 
     $     call errquit('new_giao: f_jk', lend, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lend, 'f_jl', l_f_jl, k_f_jl)) 
     $     call errquit('new_giao: f_jl', lend, MA_ERR)

      leri = ngint    ! times 6 blocks of integrals for gint(i*j*k*l) length
 
      if (.not. ma_push_get(mt_dbl,leri*6,'eri',l_eri,k_eri)) 
     $     call errquit('new_giao:could not allocate buffer',leri,
     &       MA_ERR)
      if (.not. ma_push_get(mt_dbl,iscrln,'scr',l_scr,k_scr))
     $     call errquit('new_giao: scratch alloc failed', iscrln,
     &       MA_ERR)

c ------ create g_fock1 ---------- START
      plo(1) = -1 
      plo(2) = -1
      plo(3) =  1
      phi(1) = nbf
      phi(2) = nbf
      phi(3) =  3
      if (.not.nga_create(MT_DBL,3,phi,'Fock matrix',plo,g_fock1)) call 
     &    errquit('hnd_giaox: nga_create failed g_fock',0,GA_ERR)
      call ga_zero(g_fock1)
c     == local array ==
      if (.not. ga_duplicate (g_fock1, g_loc,'g_loc'))
     &  call errquit('new_giao_2e: duplicate failed g_loc',555, GA_ERR)
      call ga_zero(g_loc)
c ------ create g_fock1 ---------- END    

      call ga_zero(g_fock)
      if      (npol.eq.1) then
       ntot=1 
       factor_cam =1.0d0 
       factor_xch =1.0d0  
      else if (npol.eq.2) then
       ntot=3
       factor_cam =2.0d0 
       factor_xch =2.0d0  
      else
       write(*,*) 'Error in new_giao_2e:: npol=1 or 2 npol=',npol
      endif
c
c     == take care of attenuated or non-attenuated 2e ==
      call ga_zero(g_fock)
      do ispin=1,npol
       call ga_zero(g_fock1) ! reset g_fock1
      if (.not.cam_exch) then ! normal calculations
        call ga_zero(g_loc)
        call giaofock(basis,geom,g_dens(ntot),g_loc,
     $     dbl_mb(k_d_kl), dbl_mb(k_d_ik), dbl_mb(k_d_jl), 
     $     dbl_mb(k_d_il), dbl_mb(k_d_jk), dbl_mb(k_f_ij), 
     $     dbl_mb(k_f_ji), dbl_mb(k_f_ik), dbl_mb(k_f_jl), 
     $     dbl_mb(k_f_il), dbl_mb(k_f_jk), iscrln, dbl_mb(k_scr), 
     $     ngint, dbl_mb(k_eri), tol2e, nshell, blen, 1.d0, 0.d0) ! get full coulomb
        call ga_add(1.d0, g_loc, 1.0d0, g_fock1, g_fock1)
c
        call ga_zero(g_loc)
        call giaofock(basis,geom,g_dens(ispin),g_loc,
     $     dbl_mb(k_d_kl), dbl_mb(k_d_ik), dbl_mb(k_d_jl), 
     $     dbl_mb(k_d_il), dbl_mb(k_d_jk), dbl_mb(k_f_ij), 
     $     dbl_mb(k_f_ji), dbl_mb(k_f_ik), dbl_mb(k_f_jl), 
     $     dbl_mb(k_f_il), dbl_mb(k_f_jk), iscrln, dbl_mb(k_scr), 
     $     ngint, dbl_mb(k_eri), tol2e, nshell, blen, 0.d0, kfac) ! get full exchange
        call ga_add(factor_xch, g_loc, 1.0d0, g_fock1, g_fock1)
      else
        call ga_zero(g_loc)
        call case_setflags(.false.) ! turn off attenuation
        call giaofock(basis,geom,g_dens(ntot),g_loc,
     $     dbl_mb(k_d_kl), dbl_mb(k_d_ik), dbl_mb(k_d_jl), 
     $     dbl_mb(k_d_il), dbl_mb(k_d_jk), dbl_mb(k_f_ij), 
     $     dbl_mb(k_f_ji), dbl_mb(k_f_ik), dbl_mb(k_f_jl), 
     $     dbl_mb(k_f_il), dbl_mb(k_f_jk), iscrln, dbl_mb(k_scr), 
     $     ngint, dbl_mb(k_eri), tol2e, nshell, blen, 1.d0, 0.d0)  ! get full coulomb
        call ga_add(1.d0, g_loc, 1.0d0, g_fock1, g_fock1)  ! full Coulomb
c
        call ga_zero(g_loc)
        call case_setflags(.false.) ! turn off attenuation
        cam_alpha1=factor_cam*cam_alpha
        cam_beta1 =factor_cam*cam_beta
        call giaofock(basis,geom,g_dens(ispin),g_loc,
     $     dbl_mb(k_d_kl), dbl_mb(k_d_ik), dbl_mb(k_d_jl), 
     $     dbl_mb(k_d_il), dbl_mb(k_d_jk), dbl_mb(k_f_ij), 
     $     dbl_mb(k_f_ji), dbl_mb(k_f_ik), dbl_mb(k_f_jl), 
     $     dbl_mb(k_f_il), dbl_mb(k_f_jk), iscrln, dbl_mb(k_scr), 
     $     ngint, dbl_mb(k_eri), tol2e, nshell, blen, 0.d0, kfac) ! get full exchange
c
        if (.not.xc_hassrhf()) then
          call ga_add(cam_alpha1, g_loc, 1.0d0, g_fock1, g_fock1)   ! cam_alpha*Exchange
        else
          call ga_add(1.d0, g_loc, 1.0d0, g_fock1, g_fock1)   ! 1.d0*Exchange
        end if
c
        call ga_zero(g_loc)
        call case_setflags(.true.) ! set attenuation
        call giaofock(basis,geom,g_dens(ispin),g_loc,
     $     dbl_mb(k_d_kl), dbl_mb(k_d_ik), dbl_mb(k_d_jl), 
     $     dbl_mb(k_d_il), dbl_mb(k_d_jk), dbl_mb(k_f_ij), 
     $     dbl_mb(k_f_ji), dbl_mb(k_f_ik), dbl_mb(k_f_jl), 
     $     dbl_mb(k_f_il), dbl_mb(k_f_jk), iscrln, dbl_mb(k_scr), 
     $     ngint, dbl_mb(k_eri), tol2e, nshell, blen, 0.d0, kfac) ! get attenuated exchange
c
        if (.not.xc_hassrhf()) then
          call ga_add(cam_beta1, g_loc, 1.0d0, g_fock1, g_fock1)  ! cam_beta*att_Exchange
        else
          call ga_add(-1.d0, g_loc, 1.0d0, g_fock1, g_fock1)  ! -1.d0*att_Exchange
        end if
c
        call case_setflags(.false.) ! turn off attenuation
      end if  ! cam_exch

        if (debug_giaofock.eq.1) then
         if (ga_nodeid().eq.0) then
          write(*,31) ispin
 31       format('---- giao_fock: g_fock1(',i3,') -------- START') 
         endif
         call ga_print(g_fock1)
         if (ga_nodeid().eq.0) then
          write(*,32) ispin
 32       format('---- giao_fock: g_fock1(',i3,') -------- END') 
         endif
        endif
c ------ copying g_fock1 --> g_fock ----- START
       plo(1) = 1
       phi(1) = nbf
       plo(2) = 1
       phi(2) = nbf
       plo(3) = 1
       phi(3) = 3
       disp=3*(ispin-1) 
       qlo(1) = 1
       qhi(1) = nbf
       qlo(2) = 1
       qhi(2) = nbf
       qlo(3) = disp+1
       qhi(3) = disp+3   
       call nga_copy_patch('n',g_fock1,plo,phi,
     &                          g_fock,qlo,qhi) 
c ------ copying g_fock1 --> g_fock ----- END
      enddo ! end-loop-ispin

c ------- Copying Coulomb contrib --------- START
       call ga_zero(g_fock_Coul)
       call ga_copy(g_fock,g_fock_Coul) 
c ------- Copying Coulomb contrib --------- END

       if (debug_giaofock.eq.1) then
        if (ga_nodeid().eq.0) then
         write(*,33) ispin
 33      format('---- giao_fock: g_fock-scaled(',i3,') -------- START') 
        endif     
          call ga_print(g_fock)
        if (ga_nodeid().eq.0) then
         write(*,34) ispin
 34      format('---- giao_fock: g_fock-scaled(',i3,') -------- END') 
        endif
       endif
       if (.not.ga_destroy(g_fock1)) call 
     &     errquit('new_giao_2e: ga_destroy failed g_fock1',
     &              0,GA_ERR)
c
c     == Clean up memory allocated in this routine ==

      if (.not. ma_chop_stack(l_d_kl))
     $  call errquit('new_giao_2e:failed chopping MA stack',555,
     &       MA_ERR)
c
c     Add DFT XC if running DFT calculation
c
      if (xc_gotxc()) then
         if (ga_nodeid().eq.0)
     &     write(*,*) 'Entering for xc, kfac=',kfac
         alo(1) = 1
         ahi(1) = nfock(npol) ! =3 npol=1, 12 npol=2
         alo(2) = 1
         ahi(2) = nbf
         alo(3) = 1
         ahi(3) = nbf
         if (.not.nga_create(MT_DBL,3,ahi,'xc matrix',alo,g_xc)) 
     &      call errquit('hnd_giaox: nga_create failed g_xc',0,GA_ERR)

c        if (ga_nodeid().eq.0)
c    &     write(*,*) 'check npol=',npol    
     
         if      (npol.eq.1) then
          do i=1,npol
           if (.not. ga_create(mt_dbl,nbf,nbf,
     &        'gCSSR: g_dia1',0,0,g_dens1(i)))
     $     call errquit('gCSSR: g_dia1',0,GA_ERR)
           call ga_zero(g_dens1(i))
           call ga_copy(g_dens(i),g_dens1(i)) 
          enddo
         else if (npol.eq.2) then
          do i=1,4
           if (.not. ga_create(mt_dbl,nbf,nbf,
     &        'gCSSR: g_dia1',0,0,g_dens1(i)))
     $     call errquit('gCSSR: g_dia1',0,GA_ERR)
          enddo
          ind=1
          do i=1,npol
           call ga_zero(g_dens1(ind))
           call ga_copy(g_dens(i),g_dens1(ind))            
           call ga_zero(g_dens1(ind+1))
           call ga_copy(g_dens(i),g_dens1(ind+1)) 
           ind=ind+2
          enddo
c -------- check input g_dens1 -------- START        
        if (debug_giaofock.eq.1) then
          do i=1,4
           if (ga_nodeid().eq.0)
     &      write(*,*) '---- gdens1(',i,') -------- START'
            call ga_print(g_dens1(i))
           if (ga_nodeid().eq.0)
     &      write(*,*) '---- gdens1(',i,') -------- END'
          enddo
        endif
c -------- check input g_dens1 -------- END
         endif
         call ga_zero(g_xc)
         Exc(1) = 0.0d0
         Exc(2) = 0.0d0
         nExc = 2 
         call fock_xc(geom,nbf,basis,nfock(npol),
     &                g_dens1,g_xc,
     &                Exc,nExc,.false.) ! out: g_xc
        if (debug_giaofock.eq.1) then
          if (ga_nodeid().eq.0)
     &    write(*,*) '---- giao_fock: g_xc -------- START'
          call ga_print(g_xc)
          if (ga_nodeid().eq.0)
     &    write(*,*) '---- giao_fock: g_xc -------- END'
        endif
     
        call ga_zero(g_fock_Exch)

          if      (npol.eq.1) then
            do ii = 1, 3
             alo(1) = ii
             ahi(1) = ii
             alo(4) = ii
             ahi(4) = ii
             call nga_add_patch(-2.0d0,g_xc,alo,ahi,
     &                            1.d0,g_fock,alo(2),ahi(2),
     &                                 g_fock,alo(2),ahi(2))
c ------- Copying Exchange contrib --------- START
            call nga_add_patch(-2.0d0,g_xc   ,alo,ahi,
     &                          1.0d0,g_fock_Exch,alo(2),ahi(2),
     &                                g_fock_Exch,alo(2),ahi(2))
c ------- Copying Exchange contrib --------- END
            enddo ! end-loop-ii
          else if (npol.eq.2) then
           acc=1
           do ispin=1,npol
            ind=coeff(ispin)*3
            do ii = 1,3
             alo(1) = ind+ii
             ahi(1) = ind+ii
             alo(4) = acc
             ahi(4) = acc
            call nga_add_patch(+2.0d0,g_xc   ,alo,ahi,
     &                          1.0d0,g_fock,alo(2),ahi(2),
     &                                g_fock,alo(2),ahi(2))
c ------- Copying Exchange contrib --------- START
            call nga_add_patch(+2.0d0,g_xc   ,alo,ahi,
     &                          1.0d0,g_fock_Exch,alo(2),ahi(2),
     &                                g_fock_Exch,alo(2),ahi(2))
c ------- Copying Exchange contrib --------- END
 108   continue
             acc=acc+1
            enddo ! end-loop-ii
           enddo ! end-loop-ispin
          endif

        if (debug_giaofock.eq.1) then
          if (ga_nodeid().eq.0)
     &    write(*,*) '---- giao_fock: g_fock-updated ----- START'
          call ga_print(g_fock)
          if (ga_nodeid().eq.0)
     &    write(*,*) '---- giao_fock: g_fock-updated -------- END'
        endif

         if (.not.ga_destroy(g_xc)) call
     &      errquit('new_giao_2e: ga_destroy failed g_xc',0,GA_ERR)
         if (npol.eq.1) then
          if (.not.ga_destroy(g_dens1(1))) call
     &       errquit('giao_aotomo: ga_destroy failed g_dens1',0,GA_ERR)
         else
          do i=1,4
          if (.not.ga_destroy(g_dens1(i))) call
     &       errquit('giao_aotomo: ga_destroy failed g_dens1',0,GA_ERR)
          enddo
         endif
      endif  ! xc_gotxc()

c
c     == clean up local array ==
      if (.not.ga_destroy(g_loc))
     $  call errquit('new_giao_2e: destroy failed g_loc',555, GA_ERR)

      end
c To produce output: g_fock_Coul,g_fock_Exch =========== END
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      subroutine giaofock(basis, geom, g_dens, g_fock,
     $                    d_kl,  d_ik,  d_jl,  d_il,  d_jk, 
     $                    f_ij,  f_ji,  f_ik,  f_jl,  f_il,  f_jk, 
     $                    lscr, scr, leri, eri, tol2e, nshell, blen,
     &                    jfac, kfac)
c
c  This routine was essentially stolen from gradients/grad2.F.  I will work
c  to make these routines one after I get things working.
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "schwarz.fh"
#include "util.fh"
#include "bas.fh"
      common/testdata/timing(20),irepeat
      double precision timing
      integer irepeat
      integer basis, geom         ! [input] familiar handles
      integer g_dens, g_fock      ! [input]/[output] density/fock ga handle
      integer blen                ! [input] max dimension of density block
      double precision            ! [scratch] 5 blocks per possible density
     $                       d_kl(blen,blen),  
     $     d_ik(blen,blen),  d_jl(blen,blen),  
     $     d_il(blen,blen),  d_jk(blen,blen)
      double precision            ! [scratch] 6 blocks per possible fock
     $     f_ij(blen,blen,3),  f_ji(blen,blen,3),  
     $     f_ik(blen,blen,3),  f_jl(blen,blen,3),  
     $     f_il(blen,blen,3),  f_jk(blen,blen,3)
c     
      integer lscr                ! [input] Size of scratch for integrals
      double precision scr(lscr ) ! [scratch] Scratch space for integrals
      integer leri                ! [input] Size of eri buffer
      double precision eri(leri,6)! [scratch] Derivative integral buffer
      double precision tol2e      ! [input] Integral screening threshold
      double precision jfac
      double precision kfac       ! [input] kfactor depends on DFT or SCF
      integer nshell              ! [input] Number of shells
      integer  next, ijblock
      integer nxtask, task_size
      external nxtask
      integer ish, jsh, ksh, lsh, idim, jdim, kdim, ldim
      integer ibflo, ibfhi, jbflo, jbfhi, kbflo, kbfhi, lbflo, lbfhi
      logical ieqj, keql, ijeqkl, keqi
      logical ijfirst,jifirst,ikfirst,jlfirst,ilfirst,jkfirst
      integer ijhand,jihand,ikhand,jlhand,ilhand,jkhand
      integer dklhand,dikhand,djlhand,dilhand,djkhand
      integer g_densrep,g_dens0,g_fockrep,g_fock0,dorep_glob
      integer nproc,shells_tot,shells_done
      logical dorepl,dorepd,dorepf,dorepon
      double precision sij,sijkl
      logical util_mirrmat
      external util_mirrmat
c
c     smax  = schwarz_max()
c     call int_acc_std()
c
      nproc = ga_nnodes()
      task_size = (nshell*(nshell+1)/2)/8
      task_size = max(1,task_size/nproc)
      ijfirst=.true.
      jifirst=.true.
      ikfirst=.true.
      jlfirst=.true.
      ilfirst=.true.
      jkfirst=.true.
c
c     replicated DM and Fock
c
c     turned off for now: needs to be debugged

      dorepon=.true.
      dorepl=.false.
      if(ga_cluster_nnodes().gt.1.and.dorepon) then
         dorepd=util_mirrmat(1,g_dens,g_densrep,.true.,.false.)
         dorepf=util_mirrmat(1,g_fock,g_fockrep,.true.,.false.)
         dorepl=dorepd.and.dorepf
         dorep_glob=0
         if(dorepl) dorep_glob=1
         call ga_igop(375,dorep_glob,1, '+')
         dorepl=dorep_glob.eq.ga_nnodes()
         if(dorepl) then
            g_dens0=g_dens
            g_dens=g_densrep
            g_fock0=g_fock
            g_fock=g_fockrep
         else
            if(dorepd) call util_mirrstop(g_dens)
         endif
      endif
c     
c     Parallel loop over shells
c     
      ijblock = 0
      next = nxtask(nproc,task_size)
      shells_tot=0
      shells_done=0
      do ish = 1, nshell
         if (.not. bas_cn2bfr(basis, ish, ibflo, ibfhi))
     $        call errquit('giaofock: bas_cn2bfr ?', basis,
     $        BASIS_ERR)
         idim = ibfhi - ibflo + 1
         do jsh = 1, ish
            if (.not. bas_cn2bfr(basis, jsh, jbflo, jbfhi))
     $           call errquit('giaofock: bas_cn2bfr ?', basis,
     $           BASIS_ERR)
            jdim = jbfhi - jbflo + 1
            call dcopy(3*idim*jdim, 0d0, 0, f_ij, 1)
            call dcopy(3*jdim*idim, 0d0, 0, f_ji, 1)
            if (next .eq. ijblock) then
               sij = schwarz_shell(ish,jsh)
               do ksh = 1, nshell
                  if (.not. bas_cn2bfr(basis, ksh, kbflo, kbfhi))
     $                 call errquit('giaofock: bas_cn2bfr ?', basis,
     $                 BASIS_ERR)
                  call new_get_giaodens_block(g_dens, idim, 
     $                 ibflo, ibfhi, kbflo, kbfhi, d_ik,dikhand)
                  call new_get_giaodens_block(g_dens, jdim, 
     $                 jbflo, jbfhi, kbflo, kbfhi, d_jk,djkhand)
                  kdim = kbfhi - kbflo + 1
                  call dcopy(3*idim*kdim, 0d0, 0, f_ik, 1)
                  call dcopy(3*jdim*kdim, 0d0, 0, f_jk, 1)
#ifdef NBGET
                     call ga_nbwait(dikhand)
                     call ga_nbwait(djkhand)
#endif
                  do lsh = 1, ksh 
                     shells_tot=shells_tot+1
                     sijkl = sij*schwarz_shell(ksh,lsh)
                     if (sijkl .gt. min(0.1*tol2e,1d-12)) then
                        shells_done=shells_done+1
                     if (.not. bas_cn2bfr(basis, lsh, lbflo, lbfhi))
     $                   call errquit('giaofock: bas_cn2bfr ?', basis,
     $                                 BASIS_ERR)
c     
c     Get blocks of the one-particle densities
c     
                     call new_get_giaodens_block(g_dens, kdim, 
     $                    kbflo, kbfhi, lbflo, lbfhi, d_kl,dklhand)
                     call new_get_giaodens_block(g_dens, jdim, 
     $                    jbflo, jbfhi, lbflo, lbfhi, d_jl,djlhand)
                     call new_get_giaodens_block(g_dens, idim, 
     $                    ibflo, ibfhi, lbflo, lbfhi, d_il,dilhand)
c
                     ldim = lbfhi - lbflo + 1
c
                     call dcopy(3*jdim*ldim, 0d0, 0, f_jl, 1)
                     call dcopy(3*idim*ldim, 0d0, 0, f_il, 1)
c     
                     call int_giao_2e(
     $                    basis,ish,jsh,basis,ksh,lsh,
     $                    lscr,scr,leri,eri)
c     
#ifdef NBGET
                     call ga_nbwait(dklhand)
                     call ga_nbwait(djlhand)
                     call ga_nbwait(dilhand)
#endif
                     call new_giao_doit(leri,eri,tol2e,
     $                          d_kl, d_ik, d_jl, d_il, d_jk, 
     $                    f_ij, f_ji, f_ik, f_jl, f_il, f_jk, 
     $                    blen, idim, jdim, kdim, ldim, ish,
     $                    jsh, ksh, lsh, jfac, kfac)

                     call new_acc_giaofock_block(g_fock, blen, 
     $                    jbflo, jbfhi, lbflo, lbfhi, f_jl,
     F                    jlfirst,jlhand)
                     call new_acc_giaofock_block(g_fock, blen, 
     $                    ibflo, ibfhi, lbflo, lbfhi, f_il,
     F                    ilfirst,ilhand)
                  endif
                  end do
                  call new_acc_giaofock_block(g_fock, blen, 
     $                 ibflo, ibfhi, kbflo, kbfhi, f_ik,
     F                 ikfirst,ikhand)
                  call new_acc_giaofock_block(g_fock, blen, 
     $                 jbflo, jbfhi, kbflo, kbfhi, f_jk,
     F                 jkfirst,jkhand)
               end do
               call new_acc_giaofock_block(g_fock, blen, 
     $              ibflo, ibfhi, jbflo, jbfhi, f_ij,
     F              ijfirst,ijhand)
               call new_acc_giaofock_block(g_fock, blen, 
     $              jbflo, jbfhi, ibflo, ibfhi, f_ji,
     F              jifirst,jihand)
               next = nxtask(nproc,task_size)
            end if
            ijblock = ijblock + 1
         end do
      end do
c
c     call int_acc_std()

      next = nxtask(-nproc,task_size)
      if(ga_nodeid().eq.0)
     c     write(6,'(i5,a,f10.2,a)') ga_nodeid(),
     G     ' giaofock schwarz done = ',
     C     (100d0*shells_done)/shells_tot,'%'
      call ga_sync()
      if(dorepl) then
         call util_mirrstop(g_densrep)
         g_dens=g_dens0
         call util_mirrmerge(g_fockrep,g_fock0)
         call util_mirrstop(g_fockrep)
         g_fock=g_fock0
      endif
c     
      return
      end

      subroutine new_get_giaodens_block(
     $     g_dens, blen, ibflo, ibfhi, jbflo, jbfhi, buf, hand)
      implicit none
      integer  g_dens, blen, ibflo, ibfhi, jbflo, jbfhi, hand
      double precision buf(*)
c
c     call dcopy(blen*blen,0d0,0,buf,1)
#ifdef NBGET
      call ga_nbget(g_dens,ibflo,ibfhi,jbflo,jbfhi,buf,blen,hand)
#else
cold      call ga_get(g_dens,ibflo,ibfhi,jbflo,jbfhi,buf,blen)
      call ga_get(g_dens,ibflo,ibfhi,jbflo,jbfhi,
     X     buf,ibfhi-ibflo+1)
#endif
c
      end
      subroutine new_acc_giaofock_block(
     $     g_fock, blen, ibflo, ibfhi, jbflo, jbfhi, buf,
     F     first,hand)
      implicit none
      integer  g_fock, blen, ibflo, ibfhi, jbflo, jbfhi
      integer alo(3), ahi(3), bln(2)
      double precision buf(*)
      logical first
      integer hand, idim, jdim
      double precision eps
      parameter(eps=1d-12)
      double precision dabsmax
      external dabsmax
c
      idim=ibfhi-ibflo+1
      jdim=jbfhi-jbflo+1
      if(dabsmax(idim*jdim*3,buf).lt.eps)
     R     return
      alo(1) = ibflo
      ahi(1) = ibfhi
      alo(2) = jbflo
      ahi(2) = jbfhi
      alo(3) = 1
      ahi(3) = 3
      bln(1) = idim
      bln(2) = jdim
#ifdef NBACC
      if(.not.first) then
         call ga_nbwait(hand)
      endif
      first=.false.
      call nga_nbacc(g_fock,alo,ahi,buf,bln,1d0,hand)
#else
      call nga_acc(g_fock,alo,ahi,buf,bln,1d0)
#endif
c
      end
c
      subroutine new_giao_doit(leri,eri,tol2e, d_kl, d_ik, d_jl, d_il, 
     $                         d_jk, f_ij, f_ji, f_ik, f_jl, f_il, f_jk,
     $                         blen,idim,jdim,kdim,ldim,ish,jsh,ksh,lsh,
     $                         jfac,kfac)
      implicit none
      integer blen
      integer leri
      integer idim, jdim, kdim, ldim, ish, jsh, ksh, lsh
      double precision eri(leri,6)
      double precision tol2e,kfac,jfac
      double precision
     $                      d_kl(kdim,ldim), 
     $     d_ik(idim,kdim), d_jl(jdim,ldim),   
     $     d_il(idim,ldim), d_jk(jdim,kdim)
      double precision
     $     f_ij(idim,jdim,3), f_ji(jdim,idim,3), 
     $     f_ik(idim,kdim,3), f_jl(jdim,ldim,3), 
     $     f_il(idim,ldim,3), f_jk(jdim,kdim,3)
      double precision fac,er2fac,errsum,errdif,
     V     erif,erif3,facij
      logical ieqj, keql, ijeqkl, keqi, doit
      integer integ, i, j, k, l, v,m

      ieqj   = ish.eq.jsh
      ijeqkl = ish.eq.ksh.and.jsh.eq.lsh
      keql   = ksh.eq.lsh
      do v = 1, 3
         integ = 0
         do i = 1, idim
            do j = 1, jdim
               facij = 1.0d0
               doit = .true.
               if (ieqj.and.j.gt.i) doit = .false.
               if(.not.doit) then
                  integ=integ+ldim*kdim
               else
               if (ieqj.and.i.eq.j) facij = 0.5d0
               do k = 1, kdim
                  do l= 1, ldim
                     integ = integ + 1
c
c  Needs the following if statements, otherwise integrals are included double
                     doit = .true.
                     if (keql.and.l.gt.k) doit = .false.
                     if (doit) then
                        fac = facij
                        if (keql.and.k.eq.l) fac = fac*0.5d0
                        erif=eri(integ,v)*fac
                        er2fac=erif*d_kl(k,l)*2.0d0*jfac
                        f_ij(i,j,v) = f_ij(i,j,v) + er2fac
                        f_ji(j,i,v) = f_ji(j,i,v) - er2fac
                        erif3=eri(integ,v+3)*fac
                        errsum=(erif+erif3)*0.5d0*kfac
                        errdif=(erif-erif3)*0.5d0*kfac
                        f_ik(i,k,v) = f_ik(i,k,v) - d_jl(j,l)*errdif
                        f_il(i,l,v) = f_il(i,l,v) - d_jk(j,k)*errsum
                        f_jk(j,k,v) = f_jk(j,k,v) + d_il(i,l)*errsum
                        f_jl(j,l,v) = f_jl(j,l,v) + d_ik(i,k)*errdif
                     endif 
                  end do ! l
               end do ! k
               endif
            end do
         end do
      end do
      end
