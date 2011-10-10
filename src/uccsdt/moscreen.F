      subroutine moints_screen(basis, occlo, occhi, g_mos, g_mosc)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "bas.fh"
#include "mafdecls.fh"
#include "schwarz.fh"
#include "tcgmsg.fh"
#include "nwc_const.fh"
      integer basis
      integer occlo, occhi      ! [input] Only to count 3x integrals
      integer g_mos             ! [input] GA of MO vectors
      integer g_mosc            ! [output] allocated and returned
c     
c     Create |(ij|ij)|^1/2, the Schwarz screening array over
c     MO integrals.
c     
c     Assumes integrals and AO schwarz already initialized
c     
      integer nbf, nmo, nsh
      integer nshblock, shlist(2,nw_max_shells)
      integer k_mos, k_tmp1, k_tmp2, k_tmp3, k_ints
      integer l_mos, l_tmp1, l_tmp2, l_tmp3, l_ints
      integer ilo, ihi, jlo, jhi, klo, khi, llo, lhi
      integer ish, jsh, kshlo, kshhi, lshlo, lshhi
      integer idim, jdim, kdim, ldim, maxshblock, intdim
      integer junk, kshb, lshb
      integer nproc, count, next
      integer ga_type
      double precision tol
      integer nxtask
      external nxtask
      parameter (tol=1e-16)     ! Needs to be the square of required acc.
c     
c     Need high precision integrals since will take sqrt() of the result
c     
      call int_acc_high()
c     
c     Get dimension info
c     
      call ga_inquire(g_mos, ga_type, nbf, nmo)
      if (.not. bas_numcont(basis,nsh))
     $     call errquit('moints_screen: basis?', basis, BASIS_ERR)
c     
c     Allocate the result array
c     
      if (.not. ga_create(mt_dbl, nmo, nmo, 'mosc', 0, 0, g_mosc))
     $     call errquit('moscreen: failed to allocate result array',0,
     &       GA_ERR)
      call ga_zero(g_mosc)
c     
c     Get the MOs
c     
      if (.not. ma_push_get(mt_dbl, nbf*nmo, 'mos',
     $     l_mos, k_mos)) call errquit
     $     ('moints_screen: allocating mos?',nbf*nmo, MA_ERR)
      if (ga_nodeid() .eq. 0) call ga_get
     $     (g_mos, 1, nbf, 1, nmo, dbl_mb(k_mos), nbf)
      call ga_brdcst(2, dbl_mb(k_mos), 8*nmo*nbf, 0)
c     
c     Allocate temporary for accumulating the result
c     
      maxshblock = min(256,nbf)
      if (.not. ma_push_get(mt_dbl, maxshblock*nmo, 'tmp1',
     $     l_tmp1, k_tmp1)) call errquit
     $     ('moints_screen: allocating tmp?',maxshblock*nmo, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nmo*nmo, 'tmp2',
     $     l_tmp2, k_tmp2)) call errquit
     $     ('moints_screen: allocating tmp?',nmo*nmo, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nmo*nmo, 'tmp2',
     $     l_tmp3, k_tmp3)) call errquit
     $     ('moints_screen: allocating tmp?',nmo*nmo, MA_ERR)
      call dfill(nmo*nmo, 0d0, dbl_mb(k_tmp3), 1)
c     
c     Block the shells for the inner loops
c     
      call moints_shell_block(basis, maxshblock, shlist, nshblock)
c     
      nproc = ga_nnodes()
      next = nxtask(nproc, 1)
      count = 0
      do ish = 1, nsh
         if (.not. bas_cn2bfr(basis, ish, ilo, ihi))
     $        call errquit('moints_screen: basis?', 0, BASIS_ERR)
         idim = ihi-ilo+1
         do jsh = 1, ish        ! Must have ish >= jsh
            if (.not. bas_cn2bfr(basis, jsh, jlo, jhi))
     $           call errquit('moints_screen: basis?', 0, BASIS_ERR)
            jdim = jhi-jlo+1
            if (count.eq.next) then 
               if (schwarz_shell(jsh,ish)*schwarz_max() .gt. tol) then
                  do kshb = 1, nshblock
                     kshlo = shlist(1,kshb)
                     kshhi = shlist(2,kshb)
                     if (.not. bas_cn2bfr(basis, kshlo, klo, junk))
     $                    call errquit('moints_screen: basis?', 0,
     &       BASIS_ERR)
                     if (.not. bas_cn2bfr(basis, kshhi, junk, khi))
     $                    call errquit('moints_screen: basis?', 0,
     &       BASIS_ERR)
                     kdim = khi-klo+1
                     do lshb = 1, kshb ! Must have ksh >= lsh
                        lshlo = shlist(1,lshb)
                        lshhi = shlist(2,lshb)
                        if (.not. bas_cn2bfr(basis, lshlo, llo, junk))
     $                       call errquit('moints_screen: basis?', 0,
     &       BASIS_ERR)
                        if (.not. bas_cn2bfr(basis, lshhi, junk, lhi))
     $                       call errquit('moints_screen: basis?', 0,
     &       BASIS_ERR)
                        ldim = lhi-llo+1
c     
c     Allocate space to hold the integrals (ij|kl) for ij in shells
c     and kl in shell blocks.  Then generate them stored as (l,k,j,i)
c     
                        intdim = idim*jdim*kdim*ldim
                        if (.not. ma_push_get(mt_dbl, intdim, 'ints',
     $                       l_ints, k_ints)) call errquit
     $                       ('moints_screen: allocating ints', 
     $                       intdim, MA_ERR)
c     
                        call moints_aoblock(
     $                       basis, .false.,
     $                       lshlo, lshhi, llo, lhi,
     $                       kshlo, kshhi, klo, khi,
     $                       jsh,   jsh,   jlo, jhi,
     $                       ish,   ish,   ilo, ihi,
     $                       tol,
     $                       dbl_mb(k_ints))
c     
c     Now transform them to the MO basis basis one shell pair at
c     a time.
c     
                        call moints_scr2(ilo,ihi,jlo,jhi,klo,khi,
     $                       llo,lhi,
     $                       dbl_mb(k_ints), dbl_mb(k_mos), nbf, nmo,
     $                       dbl_mb(k_tmp1), dbl_mb(k_tmp2),
     $                       dbl_mb(k_tmp3))
                     end do
                  end do
               endif
               next = nxtask(nproc, 1)
            end if
            count = count + 1
         end do
      end do
c     
      next = nxtask(-nproc, 1)
      call ga_sync()
c     
      call ga_dgop(1, dbl_mb(k_tmp3), nmo*nmo, '+')
      call moints_scr3(nmo, dbl_mb(k_tmp3), occlo, occhi)
      if (ga_nodeid() .eq. 0) then
         call ga_put(g_mosc, 1, nmo, 1, nmo, dbl_mb(k_tmp3), nmo)
      endif
      call ga_sync()
      if (.not. ma_chop_stack(l_mos)) call errquit
     $     ('moints_screen: failed chopping stack', 0, MA_ERR)
c     
***   if (ga_nodeid() .eq. 0) write(6,*) ' MO SCREENING ARRAY'
***   call ga_print(g_mosc)
c     
      call int_acc_std()
c     
      end
      subroutine moints_scr3(nmo, tmp, occlo, occhi)
      implicit none
#include "global.fh"
      integer nmo, occlo, occhi
      double precision tmp(nmo, nmo)
c
c     On input tmp contains (ij|ij) which needs to be symmetrized
c     On output it contains sqrt(|(ij|ij)|).
c     Also bin the values for easy interpretation of the sparsity.
c
      integer i, j, a, b, c, bb, me, nproc, count
      double precision nabic, tci, tabci
      double precision bins(-15:3)
c
      do i = -15,3
         bins(i) = 0.0d0
      enddo
c
      do i = 1, nmo
         do j = 1, i
            tmp(j,i) = sqrt(abs(tmp(j,i) + tmp(i,j)))
            tmp(i,j) = tmp(j,i)
            if (tmp(i,j) .eq. 0.0d0) then
               bins(-15) = bins(-15) + 1
            else
               b = min(max(-15,nint(log10(tmp(j,i)))),3)
               bins(b) = bins(b) + 1
            end if
         end do
      end do
c
      if (ga_nodeid() .eq. 0) then
         write(6,11) 
 11      format(//,'log10 bins for MO schwarz screening pairs ',//)
         write(6,1) (b, bins(b), b=-15,3)
 1       format(i4,2x,1pd9.1)
      endif 
c
c     Now count how many 3-x integrals we think there are
c     according to Schwarz screening.  (ab|ci)
c
      me = ga_nodeid()
      nproc = ga_nnodes()
      count = 0
c
      do i = -15,3
         bins(i) = 0.0d0
      enddo
c
      nabic = 0.0d0
      do i = occlo, occhi
         do c = occhi+1,nmo
            if (mod(count,nproc) .eq. me) then
               tci = tmp(c,i)
               do a = occhi+1,nmo
                  do b = occhi+1,nmo
                     tabci = tmp(b,a)*tci
                     if (tabci .le. 1d-15) then
                        bins(-15) = bins(-15) + 1
                     else
                        bb = min(max(-15,nint(log10(tabci))),3)
                        bins(bb) = bins(bb) + 1
                     end if
                     nabic = nabic + 1.0d0
                  end do
               end do
            endif
            count = count + 1
         end do
      end do
      call ga_dgop(55,nabic,1,'+')
      call ga_dgop(5,bins,19,'+')
      if (ga_nodeid() .eq. 0) then
         write(6,12) 
 12      format(//,'log10 bins for 3-x integrals',//)
         write(6,1) (b, bins(b), b=-15,3)
c
         do i = 2,-15,-1
            bins(i) = bins(i) + bins(i+1)
         end do
         write(6,13) 
 13      format(//,'cumulative log10 bins for 3-x integrals',//)
         write(6,1) (b, bins(b), b=-15,3)
c     
         do i = 3,-15,-1
            bins(i) = bins(i)/nabic
         end do
         write(6,14) 
 14      format(//,'scaled cumulative log10 bins for 3-x integrals',//)
         write(6,1) (b, bins(b), b=-15,3)
      end if
      call ga_sync()
c
      end
      subroutine moints_scr2(ilo,ihi,jlo,jhi,klo,khi,llo,lhi,
     $     ints, mos, nbf, nmo, tmp1, tmp2, tmp3)
      implicit none
#include "global.fh"
      integer ilo, ihi, jlo, jhi, klo, khi, llo, lhi, nbf, nmo
      double precision ints(llo:lhi,klo:khi,jlo:jhi,ilo:ihi)
      double precision mos(nbf,nmo)
      double precision tmp1(llo:lhi,nmo), tmp2(nmo,nmo), tmp3(nmo,nmo)
      integer kdim, ldim, i, j, p, q
c
      kdim = khi-klo+1
      ldim = lhi-llo+1
c
      do i = ilo,ihi
         do j = jlo,jhi
c
c     (lk|ji) -> (lq|ji)
c
c$$$            if (ga_nodeid() .eq. 1) then
c$$$               write(6,*) ' i, j, ', i, j
c$$$               write(6,*) ' ints '
c$$$               call output(ints(llo,klo,j,i),1,ldim,1,kdim,ldim,kdim,1)
c$$$               write(6,*) ' mos '
c$$$               call output(mos,1,nbf,1,nmo,nbf,nmo,1)
c$$$            endif
            call dgemm('n', 'n', ldim, nmo, kdim, 1.0d0,
     $           ints(llo,klo,j,i), ldim, mos(klo,1), nbf, 0.0d0,
     $           tmp1(llo,1), ldim)
c
c     (lq|ji) -> (pq|ji)  (will be symmetrized at very end)
c
c$$$            if (ga_nodeid() .eq. 1) then
c$$$               write(6,*) ' tmp1 '
c$$$               call output(tmp1,1,ldim,1,nmo,ldim,nmo,1)
c$$$               write(6,*) ' tmp2 '
c$$$               call output(tmp2,1,nmo,1,nmo,nmo,nmo,1)
c$$$            endif
            call dgemm('t', 'n', nmo, nmo, ldim, 1.0d0,
     $           mos(llo,1), nbf, tmp1(llo,1), ldim,  0.0d0, 
     $           tmp2(1,1), nmo)
c
c     Transform ij also
c
            do q = 1, nmo
               do p = 1, nmo
                  tmp3(p,q) = tmp3(p,q) + tmp2(p,q)*
     $                 (mos(j,p)*mos(i,q)+mos(i,p)*mos(j,q))
               end do
            end do
c
         end do
      end do
c
      end
c $Id$
