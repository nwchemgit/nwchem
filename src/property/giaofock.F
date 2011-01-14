cedo#define NBACC 1
cedo#define NBGET 1
      subroutine new_giao_2e(geom,basis,nbf,tol2e,g_dens,g_fock,kfac)
c $Id$
      implicit none
#include "nwc_const.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "hnd_giao.fh"
#include "hnd_pointers.fh"
#include "global.fh"
c
      integer geom, basis, g_dens, g_fock, nbf
      double precision tol2e,kfac,Exc(2)
      integer nExc
c
      integer nshell, maxang, blen, lend, leri, ii
      integer l_eri,k_eri,l_scr,k_scr,alo(4),ahi(4),g_xc
c
      integer       k_d_kl,  k_d_ik,  k_d_jl,  k_d_il,  k_d_jk, 
     $     k_f_ij,  k_f_ji,  k_f_ik,  k_f_jl,  k_f_il,  k_f_jk,
     $              l_d_kl,  l_d_ik,  l_d_jl,  l_d_il,  l_d_jk, 
     $     l_f_ij,  l_f_ji,  l_f_ik,  l_f_jl,  l_f_il,  l_f_jk
c
      if(.not.bas_high_angular(basis,maxang))
     &         call errquit('new_giao: angm error',maxang, BASIS_ERR)
      if (.not. bas_numcont(basis,nshell))
     $     call errquit('new_giao: could not get nsh',0, BASIS_ERR)
c
      blen = (maxang+1)*(maxang+2)/2
c
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
c
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
c
      leri = ngint    ! times 6 blocks of integrals for gint(i*j*k*l) length
c 
      if (.not. ma_push_get(mt_dbl,leri*6,'eri',l_eri,k_eri)) 
     $     call errquit('new_giao:could not allocate buffer',leri,
     &       MA_ERR)
      if (.not. ma_push_get(mt_dbl,iscrln,'scr',l_scr,k_scr))
     $     call errquit('new_giao: scratch alloc failed', iscrln,
     &       MA_ERR)
c
      call giaofock(basis,geom,g_dens,g_fock,
     $     dbl_mb(k_d_kl), dbl_mb(k_d_ik), dbl_mb(k_d_jl), 
     $     dbl_mb(k_d_il), dbl_mb(k_d_jk), dbl_mb(k_f_ij), 
     $     dbl_mb(k_f_ji), dbl_mb(k_f_ik), dbl_mb(k_f_jl), 
     $     dbl_mb(k_f_il), dbl_mb(k_f_jk), iscrln, dbl_mb(k_scr), 
     $     ngint, dbl_mb(k_eri), tol2e, nshell, blen, kfac)
c
c  Clean up memory allocated in this routine
c
      if (.not. ma_chop_stack(l_d_kl))
     $  call errquit('new_giao_2e:failed chopping MA stack',555,
     &       MA_ERR)
c
c     Add DFT XC if running DFT calculation
c
      if (kfac.lt.1.0d0) then
         alo(1) = 1
         ahi(1) = 3
         alo(2) = 1
         ahi(2) = nbf
         alo(3) = 1
         ahi(3) = nbf
         if (.not.nga_create(MT_DBL,3,ahi,'xc matrix',alo,g_xc)) 
     &      call errquit('hnd_giaox: nga_create failed g_xc',0,GA_ERR)
         call ga_zero(g_xc)
         Exc(1) = 0.0d0
         Exc(2) = 0.0d0
         nExc = 1
         call fock_xc(geom,nbf,basis,3,g_dens,g_xc,Exc,nExc,.false.)
         do ii = 1, 3
            alo(1) = ii
            ahi(1) = ii
            alo(4) = ii
            ahi(4) = ii
            call nga_add_patch(-2.0d0,g_xc,alo,ahi,1.d0,g_fock,
     &                         alo(2),ahi(2),g_fock,alo(2),ahi(2))
         enddo
         if (.not.ga_destroy(g_xc)) call
     &      errquit('giao_aotomo: ga_destroy failed g_xc',0,GA_ERR)
      endif
c
      end
      subroutine giaofock(basis, geom, g_dens, g_fock,
     $                    d_kl,  d_ik,  d_jl,  d_il,  d_jk, 
     $                    f_ij,  f_ji,  f_ik,  f_jl,  f_il,  f_jk, 
     $                    lscr, scr, leri, eri, tol2e, nshell,blen,kfac)
c
c  This routine was essentially stolen from gradients/grad2.F.  I will work
c  to make these routines one after I get things working.
c
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
c
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
      double precision kfac       ! [input] kfactor depends on DFT or SCF
      integer nshell              ! [input] Number of shells
c
      integer  next, ijblock
c
      integer nxtask, task_size
      external nxtask
c
      integer ish, jsh, ksh, lsh, idim, jdim, kdim, ldim
      integer ibflo, ibfhi, jbflo, jbfhi, kbflo, kbfhi, lbflo, lbfhi
      logical ieqj, keql, ijeqkl, keqi
c
      logical ijfirst,jifirst,ikfirst,jlfirst,ilfirst,jkfirst
      integer ijhand,jihand,ikhand,jlhand,ilhand,jkhand
      integer dklhand,dikhand,djlhand,dilhand,djkhand
      integer g_densrep,g_dens0,g_fockrep,g_fock0,dorep_glob
c
      integer nproc
      logical dorepl,dorepd,dorepf,dorepon
      logical util_mirrmat
      external util_mirrmat
c
c     smax  = schwarz_max()
c     call int_acc_std()
c
      nproc = ga_nnodes()
      task_size = nshell*(nshell+1)/2
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
c
      dorepon=.false.
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
            if(dorepd) call util_mirrstop(1,g_dens)
         endif
      endif

c     
c     Parallel loop over shells
c     
      ijblock = 0
      next = nxtask(nproc,task_size)
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
     $                    jsh, ksh, lsh, kfac)
c
                     call new_acc_giaofock_block(g_fock, blen, 
     $                    jbflo, jbfhi, lbflo, lbfhi, f_jl,
     F                    jlfirst,jlhand)
                     call new_acc_giaofock_block(g_fock, blen, 
     $                    ibflo, ibfhi, lbflo, lbfhi, f_il,
     F                    ilfirst,ilhand)
c
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
c
      next = nxtask(-nproc,task_size)
      call ga_sync()
      if(dorepl) then
         call util_mirrstop(1,g_densrep)
         g_dens=g_dens0
         call util_mirrmerge(g_fockrep,g_fock0)
         call util_mirrstop(1,g_fockrep)
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
      call ga_get(g_dens,ibflo,ibfhi,jbflo,jbfhi,buf,blen)
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
     $                         kfac)
      implicit none
c
      integer blen
      integer leri
      integer idim, jdim, kdim, ldim, ish, jsh, ksh, lsh
      double precision eri(leri,6)
      double precision tol2e, kfac
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
c
      integer integ, i, j, k, l, v,m
c
      ieqj   = ish.eq.jsh
      ijeqkl = ish.eq.ksh.and.jsh.eq.lsh
      keql   = ksh.eq.lsh
      do v = 1, 3
         integ = 0
         do i = 1, idim
            do j = 1, jdim
               facij = 1.0d0
               if (ieqj.and.i.eq.j) facij = 0.5d0
               do k = 1, kdim
                  do l= 1, ldim
                     integ = integ + 1
c
c  Needs the following if statements, otherwise integrals are included double
c
                     doit = .true.
                     if (ieqj.and.j.gt.i)               doit = .false.
                     if (keql.and.l.gt.k)               doit = .false.
                     if (doit) then
                        fac = facij
                        if (keql.and.k.eq.l) fac = fac*0.5d0
                        erif=eri(integ,v)*fac
                        er2fac=erif*d_kl(k,l)*2.0d0
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
                  end do
               end do
            end do
         end do
      end do
c
      end
