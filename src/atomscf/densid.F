      subroutine densid(dt,dold,dos,nsym,nosh,ncsh,nccup,c,damp,nconv,
     x                 nbas,nitscf,tlarge)
c.......................................................................
c
c     driver routine for density matrix processing
c.......................................................................
C$Id$
      implicit none
c
      real *8 dt(*), dold(*), dos(*), c(*)
      real *8 damp, tlarge
      integer ncsh(*), nosh(*), nccup(*),  nbas(*)
      integer nsym, nconv, nitscf
c
      real *8 vamp1, vamp2, occucl, occuop, ddiff
      integer k, nstep1, nstep2, nbas1, icount
      integer i, m, n, j
c
      nstep1 = 1
      nstep2 = 1
      k = 0
      vamp1 = 1.0d0
      vamp2 = 0.0d0
      if (nitscf.gt.1 .and. nconv.eq.0) vamp1 = 1.0d0 - damp
      if (nitscf.gt.1 .and. nconv.eq.0) vamp2 = damp
      do i = 1 , nsym
         occucl = 4*i - 2
         occuop = nccup(i)
         nbas1 = nbas(i)
         do m = 1 , nbas1
            do n = 1 , m
               k = k + 1
               dt(k) = 0.0d0
               dos(k) = 0.0d0
            enddo
         enddo
         if (ncsh(i).ne.0) call denmad(dt(nstep1),c(nstep2),ncsh(i),
     +                                 nbas1,occucl,1)
         if (nosh(i).ne.0) call denmad(dos(nstep1),c(nstep2),nosh(i),
     +                                 nbas1,occuop,ncsh(i)+1)
         nstep1 = nstep1 + nbas1*(nbas1+1)/2
         nstep2 = nstep2 + nbas1**2
      enddo
      tlarge = 0.0d0
      icount = 0
      do i = 1 , nsym
         do j = 1 , nbas(i)
            do k = 1 , j
               icount = icount + 1
               dt(icount) = (dt(icount)+dos(icount))
     +                      *vamp1 + dold(icount)*vamp2
               ddiff = dabs(dt(icount)-dold(icount))
               dold(icount) = dt(icount)
               if (ddiff.gt.tlarge) then
                  tlarge = ddiff
               end if
            enddo
         enddo
      enddo
      return
      end
