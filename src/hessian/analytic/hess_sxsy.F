C
C $Id: hess_sxsy.F,v 1.3 2000-08-22 17:08:34 windus Exp $
C
      subroutine hess_sxsy(hess, e, g_tx, oactive, ncent, nocc, restr)
c
c     Accumulates the overlap derivative contribution to the hessian
c        8 * (ei + ej) * Tij(x) * Tij(y),
c     where T(x) is -1/2 S(x) and the e's are the one-electron eigenvalues.
c
c     BGJ (5/00)
c
      implicit none
c
#include "global.fh"
#include "mafdecls.fh"
c
      integer ncent             ! [input] number of nuclear centers
      integer nocc              ! [input] number of occupied orbitals
      double precision hess(3*ncent,3*ncent) ! [updated] hessian matrix
      double precision e(nocc)  ! [input] eigenvalues
      integer g_tx(3*ncent)     ! [input] GA handles for overlap matrix
c                               !         first derivatives (in Tij(x) form)
      logical oactive(ncent)    ! [input] key to active atoms
      logical restr             ! [input] calculation is restricted
c
      integer icent, idir, ipert, jcent, jdir, jpert, iocc, jocc
      integer ioff, joff, l_scr, k_scr
      double precision deltah, fac
      integer nproc, icount, itask, mytask
      integer nxtask
      external nxtask
c
c     Prefactor of hessian contribution
c
      fac = 8d0
c     if (restr) fac = fac * 2d0
c
c     Allocate local scratch space for overlap integral derivatives
c
      if (.not. ma_push_get(MT_DBL, 2*nocc,'hess_sxsy: scr',
     &     l_scr, k_scr))
     &     call errquit('hess_sxsy: cannot allocate scr',0)
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
                     call ga_get(g_tx(ipert), 1, nocc, jocc, jocc,
     &                    dbl_mb(ioff+1), nocc)
                     call ga_get(g_tx(jpert), 1, nocc, jocc, jocc,
     &                    dbl_mb(joff+1), nocc)
                     do iocc = 1, nocc
                        deltah = deltah + (e(iocc)+e(jocc))
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
     &     call errquit('hess_sxsy: failed to pop l_scr',0)
c
      return
      end
