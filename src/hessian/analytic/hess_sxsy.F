C
C $Id: hess_sxsy.F,v 1.7 2003-10-17 22:55:13 carlfahl Exp $
C
      subroutine hess_sxsy(hess, e, g_tx, oactive, ncent, nocc, 
     &                      noffset, urestr)
c
c     Accumulates the overlap derivative contribution to the hessian
c        FAC * (ei + ej) * Tij(x) * Tij(y),
c     where T(x) is -1/2 S(x) and the e's are the one-electron 
c     eigenvalues.  Note that this is called independantly for 
c     both alpha and beta with the noffset giving the offset into
c     the g_tx array.
c
      implicit none
#include "errquit.fh"
c
#include "global.fh"
#include "mafdecls.fh"
c
      integer ncent             ! [input] number of nuclear centers
      integer nocc              ! [input] number of occupied orbitals
      integer noffset           ! [input] offset into g_tx (UHF)
      double precision hess(3*ncent,3*ncent) ! [updated] hessian matrix
      double precision e(nocc)  ! [input] eigenvalues
      integer g_tx              ! [input] GA handles for overlap matrix
c                               !         first derivatives (in Tij(x) form)
      logical oactive(ncent)    ! [input] key to active atoms
      logical urestr            ! [input] calculation is unrestricted
c
      integer icent, idir, ipert, jcent, jdir, jpert, iocc, jocc
      integer ioff, joff, l_scr, k_scr
      integer lo(3), hi(3), ld(2)
      double precision deltah, fac
      integer nproc, icount, itask, mytask
      integer nxtask
      external nxtask
c
c     Setup some arrays
c
      lo(2) = 1
      hi(2) = nocc
      ld(1) = 1
      ld(2) = nocc
c
c     Prefactor of hessian contribution
c
      if (urestr) then
         fac = 4.0d0
      else
         fac = 8.0d0
      endif
c
c     Allocate local scratch space for overlap integral derivatives
c
      if (.not. ma_push_get(MT_DBL, 2*nocc,'hess_sxsy: scr',
     &     l_scr, k_scr))
     &     call errquit('hess_sxsy: cannot allocate scr',0, MA_ERR)
      ioff = k_scr - 1
      joff = ioff + nocc
c
c     Parallelize over hessian elements
c
      nproc = ga_nnodes()
      icount = 1
      itask = 0
      mytask = nxtask(nproc,icount)
c
      do icent = 1, ncent
         if (.not. oactive(icent)) goto 10
         ipert = 3*(icent-1)
         do idir = 1, 3
            ipert = ipert + 1
c
            do jcent = 1, icent
               if (.not. oactive(jcent)) goto 20
               jpert = 3*(jcent-1)
               do jdir = 1, 3
                  jpert = jpert + 1
                  if (jpert .gt. ipert) goto 30
                  if (itask .ne. mytask) goto 25
c
c     Compute increment to hess(ipert,jpert)
c
                  deltah = 0d0
c
c     Loop over occupied MO's and compute contribution a column at a time
c
                  do jocc = 1, nocc
                     lo(1) = ipert + noffset
                     hi(1) = ipert + noffset
                     lo(3) = jocc
                     hi(3) = jocc
                     call nga_get(g_tx, lo, hi,
     &                    dbl_mb(ioff+1), ld)
                     lo(1) = jpert + noffset
                     hi(1) = jpert + noffset
                     call nga_get(g_tx, lo, hi,
     &                    dbl_mb(joff+1), ld)
                     do iocc = 1, nocc
                        deltah = deltah
     &                       + (e(iocc)+e(jocc))
     &                       * dbl_mb(ioff+iocc) * dbl_mb(joff+iocc)
                     enddo
                  enddo
c
                  deltah = deltah * fac
c
c     Update hessian
c
                  hess(ipert,jpert) = hess(ipert,jpert) + deltah
                  if (ipert .ne. jpert) then
                     hess(jpert,ipert) = hess(jpert,ipert) + deltah
                  endif
c
                  mytask = nxtask(nproc,icount)
 25               continue
                  itask = itask + 1
c
 30               continue
               enddo
c
            enddo
 20         continue
         enddo
c
 10      continue
      enddo
c
      mytask = nxtask(-nproc,icount)
c
      if (.not. ma_pop_stack(l_scr))
     &     call errquit('hess_sxsy: failed to pop l_scr',0, MA_ERR)
c
      return
      end
