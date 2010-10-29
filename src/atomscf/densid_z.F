      subroutine densid_z(dt,dold,dos,nsym,nosh,ncsh,nccup,c,damp,nconv,
     x                 nbas,nitscf,tlarge, tag, eps)
c.......................................................................
c
c     driver routine for density matrix processing
c.......................................................................
C$Id$
      implicit none
#include "errquit.fh"
#include "atomz.fh"
#include "inp.fh"
c
      real *8 dt(*), dold(*), dos(*), c(*), eps(*)
      real *8 damp, tlarge
      integer ncsh(*), nosh(*), nccup(*),  nbas(*)
      integer nsym, nconv, nitscf
c
      character*(*) tag ! Tag being computed
      real *8 delta_n, occ_remove, occ_vec(50)
      integer ind, num_occupied, indx
      real *8 vamp1, vamp2, occucl, occuop, ddiff
      integer k, nstep1, nstep2, nbas1, icount, kk
      integer i, m, n, j
c
      delta_n = 0.0d00
      if (num_z .gt. 0) then
        if (inp_match(num_z,.true.,tag,tags_z,ind)) then
          delta_n = -z(ind)
        endif
      endif
c
      call dfill(50, 0.0d0, occ_vec, 1)
c
      k=1
      do i = 1 , nsym
        occucl = 4*i - 2
        occuop = nccup(i)
        if (delta_n.gt.0.01d00.and.nosh(i).ne.0) then
          occuop = occuop + delta_n
          delta_n = 0.0d00
        else if ((delta_n.lt.-0.01d00).and.(nosh(i).ne.0)) then
          occ_remove = max(-occuop,delta_n)
          occuop     = occuop  + occ_remove
          delta_n    = delta_n - occ_remove
        endif
        do j = 1,ncsh(i)
          occ_vec(k) = occucl
          k = k + 1
        enddo
        if (nosh(i).ne.0) then
          occ_vec(k) = occuop
          k = k + 1
        endif
      enddo
      if (delta_n.gt.0.01d00) call errquit
     &('atomscf:densid_z: cannot add electrons to closed shell atom ',
     &    911, UNKNOWN_ERR)
c
      num_occupied = k - 1
c
00001   continue
      if (delta_n.lt.-0.01d00) then
*  assume first orbital will be a 1s and will have lowest energy
        indx = 1
        do k = 2,num_occupied
          if(occ_vec(k).gt.0.0d00) then 
            if (eps(k).gt.eps(indx))indx = k
          endif
        enddo
        occ_remove    = max(-occ_vec(indx),delta_n)
        occ_vec(indx) = occ_vec(indx)  + occ_remove
        delta_n       = delta_n - occ_remove
        if (occ_remove.eq.0.0d00) call errquit
     &      ('atomscf:densid_z: you cannot have positrons in NWChem',
     &      911, INPUT_ERR)
        goto 00001
      endif
c
*      call output(occ_vec, 1, num_occupied, 1, 1, num_occupied, 1, 1)
c
      nstep1 = 1
      nstep2 = 1
      k = 0
      kk = 0
      vamp1 = 1.0d0
      vamp2 = 0.0d0
      if (nitscf.gt.1 .and. nconv.eq.0) vamp1 = 1.0d0 - damp
      if (nitscf.gt.1 .and. nconv.eq.0) vamp2 = damp
      do i = 1 , nsym
         nbas1 = nbas(i)
         do m = 1 , nbas1
            do n = 1 , m
               kk = kk + 1
               dt(kk) = 0.0d0
               dos(kk) = 0.0d0
            enddo
         enddo
         do j=1,ncsh(i)
           k = k + 1
           call denmad(dt(nstep1),c(nstep2),1,nbas1,occ_vec(k),j)
         enddo
         if (nosh(i).ne.0) then
           k = k + 1
           call denmad(dos(nstep1),c(nstep2),1,nbas1,
     &         occ_vec(k),(ncsh(i)+1))
         endif
         nstep1 = nstep1 + nbas1*(nbas1+1)/2
         nstep2 = nstep2 + nbas1**2
      enddo
c
      if (delta_n .ne. 0.0d00)
     &    call errquit
     &    ('atomscf:densid_z: adding charge failed somewhere', 911,
     &       UNKNOWN_ERR)
c
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
