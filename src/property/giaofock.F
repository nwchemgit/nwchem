      subroutine new_giao_2e(geom, basis, tol2e, g_dens, g_fock)
      implicit none
#include "nwc_const.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "hnd_giao.fh"
#include "hnd_pointers.fh"
c
      integer geom, basis, g_dens, g_fock
      double precision tol2e
c
      integer nshell, maxang, blen, lend, leri, i
      integer l_eri,k_eri,l_scr,k_scr
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
     $     ngint, dbl_mb(k_eri), tol2e, nshell, blen)
c
c  Clean up memory allocated in this routine
c
      if (.not. ma_chop_stack(l_d_kl))
     $  call errquit('new_giao_2e:failed chopping MA stack',555,
     &       MA_ERR)
c
      end
      subroutine giaofock(basis, geom, g_dens, g_fock,
     $                    d_kl,  d_ik,  d_jl,  d_il,  d_jk, 
     $                    f_ij,  f_ji,  f_ik,  f_jl,  f_il,  f_jk, 
     $                    lscr, scr, leri, eri, tol2e, nshell, blen)
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
      integer nshell              ! [input] Number of shells
c
      integer  next, nint, ijklblock
      double precision scale
c
      integer nxtask, task_size
      external nxtask
c
      integer ish, jsh, ksh, lsh, idim, jdim, kdim, ldim
      integer ibflo, ibfhi, jbflo, jbfhi, kbflo, kbfhi, lbflo, lbfhi
      double precision smax
      logical ieqj, keql
c
      integer nproc,i,j
c     
c     smax  = schwarz_max()
c     call int_acc_std()
c
      nproc = ga_nnodes()
      task_size = nshell*(nshell+1)/2
      task_size = task_size*(task_size+1)/2
      task_size = max(1,task_size/(20*nproc))
c     
c     Parallel loop over shells
c     
      ijklblock = 0
      next = nxtask(nproc,task_size)
      do ish = 1, nshell
         do jsh = 1, ish
            do ksh = 1, nshell
               do lsh = 1, ksh 
                  if (next .eq. ijklblock) then
                     if (.not. bas_cn2bfr(basis, ish, ibflo, ibfhi))
     $                   call errquit('giaofock: bas_cn2bfr ?', basis,
     $                                 BASIS_ERR)
                     idim = ibfhi - ibflo + 1
                     if (.not. bas_cn2bfr(basis, jsh, jbflo, jbfhi))
     $                   call errquit('giaofock: bas_cn2bfr ?', basis,
     $                                 BASIS_ERR)
                     jdim = jbfhi - jbflo + 1
                     if (.not. bas_cn2bfr(basis, ksh, kbflo, kbfhi))
     $                   call errquit('giaofock: bas_cn2bfr ?', basis,
     $                                 BASIS_ERR)
                     kdim = kbfhi - kbflo + 1
                     if (.not. bas_cn2bfr(basis, lsh, lbflo, lbfhi))
     $                   call errquit('giaofock: bas_cn2bfr ?', basis,
     $                                 BASIS_ERR)
                     ldim = lbfhi - lbflo + 1
c
                     ieqj = ish.eq.jsh
                     keql = ksh.eq.lsh
c     
c     Get blocks of the one-particle densities
c     
                     call new_get_giaodens_block(g_dens, blen, 
     $                    kbflo, kbfhi, lbflo, lbfhi, d_kl)
                     call new_get_giaodens_block(g_dens, blen, 
     $                    ibflo, ibfhi, kbflo, kbfhi, d_ik)
                     call new_get_giaodens_block(g_dens, blen, 
     $                    jbflo, jbfhi, lbflo, lbfhi, d_jl)
                     call new_get_giaodens_block(g_dens, blen, 
     $                    ibflo, ibfhi, lbflo, lbfhi, d_il)
                     call new_get_giaodens_block(g_dens, blen, 
     $                    jbflo, jbfhi, kbflo, kbfhi, d_jk)
c
                     call dfill(3*blen*blen, 0d0, f_ij, 1)
                     call dfill(3*blen*blen, 0d0, f_ji, 1)
                     call dfill(3*blen*blen, 0d0, f_ik, 1)
                     call dfill(3*blen*blen, 0d0, f_jl, 1)
                     call dfill(3*blen*blen, 0d0, f_il, 1)
                     call dfill(3*blen*blen, 0d0, f_jk, 1)
c     
                     call int_giao_2e(
     $                    basis,ish,jsh,basis,ksh,lsh,
     $                    lscr,scr,leri,eri)
c     
                     call new_giao_doit(leri,eri,tol2e,
     $                          d_kl, d_ik, d_jl, d_il, d_jk, 
     $                    f_ij, f_ji, f_ik, f_jl, f_il, f_jk, 
     $                    blen, idim, jdim, kdim, ldim, ieqj, keql)
c
                     call new_acc_giaofock_block(g_fock, blen, 
     $                    ibflo, ibfhi, jbflo, jbfhi, f_ij)
                     call new_acc_giaofock_block(g_fock, blen, 
     $                    jbflo, jbfhi, ibflo, ibfhi, f_ji)
                     call new_acc_giaofock_block(g_fock, blen, 
     $                    ibflo, ibfhi, kbflo, kbfhi, f_ik)
                     call new_acc_giaofock_block(g_fock, blen, 
     $                    jbflo, jbfhi, lbflo, lbfhi, f_jl)
                     call new_acc_giaofock_block(g_fock, blen, 
     $                    ibflo, ibfhi, lbflo, lbfhi, f_il)
                     call new_acc_giaofock_block(g_fock, blen, 
     $                    jbflo, jbfhi, kbflo, kbfhi, f_jk)
c                                         
                     next = nxtask(nproc,task_size)
                  end if
c     
                  ijklblock = ijklblock + 1
c     
               end do
            end do
         end do
      end do
c
c     call int_acc_std()
c
      next = nxtask(-nproc,task_size)
      call ga_sync()
c     
      return
      end

      subroutine new_get_giaodens_block(
     $     g_dens, blen, ibflo, ibfhi, jbflo, jbfhi, buf)
      implicit none
      integer  g_dens, blen, ibflo, ibfhi, jbflo, jbfhi
      double precision buf(*)
c
      call dfill(blen*blen,0d0,buf,1)
      call ga_get(g_dens,ibflo,ibfhi,jbflo,jbfhi,buf,blen)
c
      end
      subroutine new_acc_giaofock_block(
     $     g_fock, blen, ibflo, ibfhi, jbflo, jbfhi, buf)
      implicit none
      integer  g_fock, blen, ibflo, ibfhi, jbflo, jbfhi
      integer alo(3), ahi(3), bln(2)
      double precision buf(blen,blen,3)
c
      alo(1) = ibflo
      ahi(1) = ibfhi
      alo(2) = jbflo
      ahi(2) = jbfhi
      alo(3) = 1
      ahi(3) = 3
      bln(1) = blen
      bln(2) = blen
      call nga_acc(g_fock,alo,ahi,buf,bln,1d0)
c
      end
c
      subroutine new_giao_doit(leri,eri,tol2e, d_kl, d_ik, d_jl, d_il, 
     $                         d_jk, f_ij, f_ji, f_ik, f_jl, f_il, f_jk,
     $                         blen, idim, jdim, kdim, ldim, ieqj, keql)
      implicit none
c
      integer blen
      integer leri
      double precision eri(leri,6)
      double precision tol2e
      double precision
     $                      d_kl(blen,blen), 
     $     d_ik(blen,blen), d_jl(blen,blen), 
     $     d_il(blen,blen), d_jk(blen,blen)
      double precision
     $     f_ij(blen,blen,3), f_ji(blen,blen,3), 
     $     f_ik(blen,blen,3), f_jl(blen,blen,3), 
     $     f_il(blen,blen,3), f_jk(blen,blen,3)
      double precision fac
      integer idim, jdim, kdim, ldim
      logical ieqj, keql
c
      integer integ, i, j, k, l, v,m
c
      do v = 1, 3
         integ = 0
         do i = 1, idim
            do j = 1, jdim
               do k = 1, kdim
                  do l= 1, ldim
                     integ = integ + 1
c
c  Needs the following if statements, otherwise integrals are included double
c
                     fac = 1.0d0
                     if (ieqj.and.i.ge.j) fac = fac/2.0d0
                     if (keql.and.k.ge.l) fac = fac/2.0d0
            f_ij(i,j,v) = f_ij(i,j,v) + d_kl(k,l)*eri(integ,v)*2.0d0*fac
            f_ji(j,i,v) = f_ji(j,i,v) - d_kl(k,l)*eri(integ,v)*2.0d0*fac
            f_ik(i,k,v) = f_ik(i,k,v) - d_jl(j,l)*
     &                        (eri(integ,v)-eri(integ,v+3))*0.5d0*fac
            f_il(i,l,v) = f_il(i,l,v) - d_jk(j,k)*
     &                        (eri(integ,v)+eri(integ,v+3))*0.5d0*fac
            f_jk(j,k,v) = f_jk(j,k,v) - d_il(i,l)*
     &                        (-eri(integ,v)-eri(integ,v+3))*0.5d0*fac
            f_jl(j,l,v) = f_jl(j,l,v) - d_ik(i,k)*
     &                        (-eri(integ,v)+eri(integ,v+3))*0.5d0*fac
                  end do
               end do
            end do
         end do
      end do
c
      end
