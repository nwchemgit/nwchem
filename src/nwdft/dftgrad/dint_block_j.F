      subroutine dint_block_jfrc(d_ija, d_kla, ldim, frc, dint, scale,
     &     idatom, if1, il1, if2, il2, if3, il3, if4, il4)
c
c     Makes contribution to J force from block of derivative integrals
c
c     BGJ - 9/98
c
c     $Id$
c
      implicit none
c
      integer ldim
      integer if1, il1, if2, il2, if3, il3, if4, il4
      double precision scale
      double precision frc(3, *)
      double precision d_ija(ldim,ldim), d_kla(ldim,ldim)
      double precision dint(if4:il4, if3:il3, if2:il2, if1:il1, 3, 4)
      integer idatom(4)
c
      integer i, icart
      integer ip1, ip2, ip3, ip4
      double precision dE, den2
c
      do i = 1, 4
         if (idatom(i).ge.1) then
            do icart = 1, 3
               dE = 0.D0
               do ip1 = if1, il1
                  do ip2 = if2, il2
                     do ip3 = if3, il3
                        do ip4 = if4, il4
                           den2 = 0.5D0 *
     &                            (d_ija(ip1,ip2) * d_kla(ip3,ip4))
                           dE = dE
     &                        + den2 * dint(ip4,ip3,ip2,ip1,icart,i)
                        enddo
                     enddo
                  enddo
               enddo
c
               frc(icart,idatom(i)) = frc(icart,idatom(i)) + 
     &                                dE * scale
c
            enddo
         endif
      enddo 
c
      return
      end
c
c
c
      subroutine dint_block_jdmat(d_ija, d_kla, ldim, g_dj, dfij, dfkl,
     &     dint, q4, diagb, diagk, ijij, idatom, if1, il1, if2, il2,
     &     if3, il3, if4, il4)
c
c     Makes contribution to J derivative matrices from block of
c     derivative integrals
c
c     BGJ - 9/98
c
c     $Id$
c
      implicit none
cpgi$g opt=1
c
      integer ldim
      integer if1, il1, if2, il2, if3, il3, if4, il4
      logical diagb, diagk, ijij
      integer g_dj(3, *)
      double precision dfij(if2:il2, if1:il1)
      double precision dfkl(if4:il4, if3:il3)
      double precision d_ija(ldim,ldim), d_kla(ldim,ldim)
      double precision dint(if4:il4, if3:il3, if2:il2, if1:il1, 3, 4)
      double precision q4
      integer idatom(4)
c
      integer i, icart
      integer ip1, ip2, ip3, ip4
      double precision scaleb, scalek, g
c
c     Set up bra and ket scale factors; 2's for permutational symmetry,
c     and 1/2's to fix diagonal J matrix blocks because we fill out
c     the matrices at the end by adding the transpose
c
      scaleb = q4
      scalek = q4
      if (diagb) then
         scalek = scalek * 0.5d0
      else
         scaleb = scaleb * 2.0d0
      endif
      if (diagk) then
         scaleb = scaleb * 0.5d0
      else
         scalek = scalek * 2.0d0
      endif
      if (ijij) then
c        Scale (ij|ij) type integrals by 0.5 to undo double-counting
         scaleb = scaleb * 0.5d0
         scalek = scalek * 0.5d0
      endif
c
      do i = 1, 4
         if (idatom(i).ge.1) then
            do icart = 1, 3
               call dfill((il2-if2+1)*(il1-if1+1), 0.0d0, dfij, 1)
               call dfill((il4-if4+1)*(il3-if3+1), 0.0d0, dfkl, 1)
               do ip1 = if1, il1
                  do ip2 = if2, il2
                     do ip3 = if3, il3
                        do ip4 = if4, il4
                           g = dint(ip4,ip3,ip2,ip1,icart,i)
                           dfij(ip2,ip1) = dfij(ip2,ip1)
     &                         + g * d_kla(ip3-if3+1,ip4-if4+1) * scalek
                           dfkl(ip4,ip3) = dfkl(ip4,ip3)
     &                         + g * d_ija(ip1-if1+1,ip2-if2+1) * scaleb
                        enddo
                     enddo
                  enddo
               enddo
c
c     Update derivative J matrix blocks in GA's
c
c     The routine fock_upd_blk takes a temporary array for space
c     to rearrange blocks before updating.  This array is not used
c     in the case of a single matrix.  The call below relies upon
c     this and does not pass in a temporary array.  If the call
c     below is ever changed to update more than one matrix at once,
c     a temporary array will have to be provided.
c
               call fock_upd_blk(1, g_dj(icart,idatom(i)), if2, il2,
     &              if1, il1, 1.d0, dfij, 0.d0)
               call fock_upd_blk(1, g_dj(icart,idatom(i)), if4, il4,
     &              if3, il3, 1.d0, dfkl, 0.d0)
            enddo
         endif
      enddo
c
      return
      end
c
c
c
      subroutine dint_2e3c_block_jdmat(CD_coef, d_kla, nbf_cd, ldim,
     &     g_dj, db, dfkl, dint, diagkl, idatom, nf1,
     &     if3, il3, if4, il4)
c
c     Makes contribution to fitted J derivative matrices from block of
c     derivative 2e3c integrals
c
c     BGJ - 3/99
c
c     $Id$
c
      implicit none
c
      integer nbf_cd, ldim
      integer nf1, if3, il3, if4, il4
      logical diagkl
      integer g_dj(3, *)
      double precision CD_coef(*), d_kla(ldim,ldim)
      double precision db(nbf_cd, 3, *), dfkl(if4:il4, if3:il3)
      double precision dint(if4:il4, if3:il3, nf1, 3, 4)
      integer idatom(4)
c
      integer i, icart
      integer ip1, ip3, ip4
      double precision scalekl, scalecd, g
c
c     Set up scale factors: 2's for permutational symmetry,
c     and 1/2's to fix diagonal J matrix blocks because we fill out
c     the matrices at the end by adding the transpose
c
      scalekl = 1d0
      scalecd = 1d0
      if (diagkl) then
         scalecd = scalecd * 0.5d0
      else
         scalekl = scalekl * 2.0d0
      endif
c
      do i = 1, 4
         if (idatom(i).ge.1) then
            do icart = 1, 3
               call dfill((il4-if4+1)*(il3-if3+1), 0.0d0, dfkl, 1)
               do ip1 = 1, nf1
                  do ip3 = if3, il3
                     do ip4 = if4, il4
                        g = dint(ip4,ip3,ip1,icart,i)
                        dfkl(ip4,ip3) = dfkl(ip4,ip3)
     &                       + g * CD_coef(ip1) * scalecd
                        db(ip1,icart,idatom(i))
     &                       = db(ip1,icart,idatom(i))
     &                       + g * d_kla(ip3-if3+1,ip4-if4+1) * scalekl
                     enddo
                  enddo
               enddo
c
c     Update derivative J matrix block in GA
c
c     The routine fock_upd_blk takes a temporary array for space
c     to rearrange blocks before updating.  This array is not used
c     in the case of a single matrix.  The call below relies upon
c     this and does not pass in a temporary array.  If the call
c     below is ever changed to update more than one matrix at once,
c     a temporary array will have to be provided.
c
               call fock_upd_blk(1, g_dj(icart,idatom(i)), if4, il4,
     &              if3, il3, 1.d0, dfkl, 0.d0)
            enddo
         endif
      enddo
c
      return
      end
c
c
c
      subroutine dint_2e2c_block_jdmat(CD_coef, db, nbf_cd, dint,
     &     idatom, if1, il1, if2, il2)
c
c     Makes contribution to fitted J derivative matrices from block of
c     derivative 2e2c integrals
c
c     BGJ - 3/99
c
c     $Id$
c
      implicit none
c
      integer nbf_cd
      integer if1, il1, if2, il2
      double precision CD_coef(nbf_cd)
      double precision db(nbf_cd, 3, *)
      double precision dint(if2:il2, if1:il1, 3, 2)
      integer idatom(2)
c
      integer i, icart
      integer ip1, ip2
      double precision g
c
c     Don't need to do separate case for diagonal blocks since these
c     give a zero contribution
c
      do i = 1, 2
         if (idatom(i).ge.1) then
            do icart = 1, 3
               do ip1 = if1, il1
                  do ip2 = if2, il2
                     g = dint(ip2,ip1,icart,i)
#if 0
                     db(ip1,icart,idatom(i)) = db(ip1,icart,idatom(i))
     &                    - g * CD_coef(ip2)
                     db(ip2,icart,idatom(i)) = db(ip2,icart,idatom(i))
     &                    - g * CD_coef(ip1)
#else
                     db(ip1,icart,idatom(i)) = db(ip1,icart,idatom(i))
     &                    - g * CD_coef(ip2) * 0.5d0
                     db(ip2,icart,idatom(i)) = db(ip2,icart,idatom(i))
     &                    - g * CD_coef(ip1) * 0.5d0
#endif
                  enddo
               enddo
            enddo
         endif
      enddo
c
      return
      end
c
c
c
      subroutine d2int_block_jhess(d_ija, d_kla, ldim, hess, nat, d2int,
     &     scale, idatom, if1, il1, if2, il2, if3, il3, if4, il4)
c
c     Makes contribution to J hessian from block of second derivative
c     integrals
c
c     BGJ - 9/98
c
c     $Id$
c
      implicit none
c
      integer ldim, nat
      integer if1, il1, if2, il2, if3, il3, if4, il4
      double precision scale
      double precision hess(3, nat, 3, nat)
      double precision d_ija(ldim,ldim), d_kla(ldim,ldim)
      double precision d2int(78, if4:il4, if3:il3, if2:il2, if1:il1)
      integer idatom(4)
c
      integer i, icart, j, jcart, indi, indj, indab, indij
      integer ip1, ip2, ip3, ip4
      double precision d2E, den2
      double precision h(3,4,3,4)
c !!! Start BGJ temp test !!!
      double precision ti(3,4,3,4)
      character abcd(4), xyz(3)
      abcd(1) = 'A'
      abcd(2) = 'B'
      abcd(3) = 'C'
      abcd(4) = 'D'
      xyz(1) = 'x'
      xyz(2) = 'y'
      xyz(3) = 'z'
c !!! End BGJ temp test !!!
c
c      write(*,*)'*** d2int_block_jhess called: idatom',idatom
      call dfill(144, 0.0d0, h, 1)
      indab = 0
      do i = 1, 4
         do j = i, 4
            if (idatom(j).ge.1 .and. idatom(i).ge.1) then
               indij = indab
               do icart = 1, 3
                  indi = 3*(i-1) + icart
                  do 100 jcart = 1, 3
                     indj = 3*(j-1) + jcart
                     if (indj .lt. indi) goto 100
                     indij = indij + 1
                     d2E = 0.D0
                     do ip1 = if1, il1
                        do ip2 = if2, il2
                           do ip3 = if3, il3
                              do ip4 = if4, il4
                                 den2 = 0.5D0 * (d_ija(ip1,ip2)
     &                                        *  d_kla(ip3,ip4))
                                 d2E = d2E + den2 *
     &                                 d2int(indij,ip4,ip3,ip2,ip1)
c       write(*,*) abcd(i), xyz(icart), abcd(j), xyz(jcart), indij,
c     &                                d2int(indij,ip4,ip3,ip2,ip1)
                              enddo
                           enddo
                        enddo
                     enddo
                     h(icart,i,jcart,j) = h(icart,i,jcart,j)
     &                    + d2E * scale
 100              continue
               enddo
            endif
            if (i.ne.j) then
               indab = indab + 9
            else
               indab = indab + 6
            endif
         enddo
      enddo
#if 0
c     Check TI of integral 2nd derivatives
      do ip1 = if1, il1
        do ip2 = if2, il2
          do ip3 = if3, il3
            do ip4 = if4, il4
              call dfill(144, 0.0d0, ti, 1)
              indab = 0
              do i = 1, 4
                do j = i, 4
                  if (idatom(j).ge.1 .and. idatom(i).ge.1) then
                    indij = indab
                    do icart = 1, 3
                      indi = 3*(i-1) + icart
                      do 200 jcart = 1, 3
                        indj = 3*(j-1) + jcart
                        if (indj .lt. indi) goto 200
                        indij = indij + 1
c       write(*,*) abcd(i), xyz(icart), abcd(j), xyz(jcart), indij,
c     &                                d2int(indij,ip4,ip3,ip2,ip1)
                        ti(icart,i,jcart,j)=d2int(indij,ip4,ip3,ip2,ip1)
                        ti(jcart,j,icart,i)=d2int(indij,ip4,ip3,ip2,ip1)
  200                 continue
                    enddo
                  endif
                  if (i.ne.j) then
                    indab = indab + 9
                  else
                    indab = indab + 6
                  endif
                enddo
              enddo
              write(*,*)'*** 2nd deriv TI check for',ip1,ip2,ip3,ip4
              call hessti(ti,4)
            enddo
          enddo
        enddo
      enddo
#endif
c
c     Scatter the hessian contribution just formed to a full square
c
      do i = 1, 4
         do j = i, 4
            if (idatom(j).ge.1 .and. idatom(i).ge.1) then
               if (i.eq.j) then
                  do icart = 1, 2
                     do jcart = icart+1, 3
                        h(jcart,j,icart,i) = h(icart,i,jcart,j)
                     enddo
                  enddo
               else
                  do icart = 1, 3
                     do jcart = 1, 3
                        h(jcart,j,icart,i) = h(icart,i,jcart,j)
                     enddo
                  enddo
               endif
            endif
         enddo
      enddo
c
c     Increment the hessian
c
      do i = 1, 4
         do j = 1, 4
            if (idatom(j).ge.1 .and. idatom(i).ge.1) then
               do icart = 1, 3
                  do jcart = 1, 3
                     hess(icart,idatom(i),jcart,idatom(j)) =
     &               hess(icart,idatom(i),jcart,idatom(j)) +
     &                   h(icart,i,jcart,j)
                  enddo
               enddo
            endif
         enddo
      enddo
c      write(*,*)'*** TI check for hessian block'
c      call hessti(h,4)
c
      return
      end
c
c
c
      subroutine d2int_2e3c_block_jhess(coef1, d_kla, ldim, hess, nat,
     &     d2int, scale, idatom, nf1, nf3, nf4)
c
c     Makes contribution to J hessian from block of second derivative
c     3-center integrals
c
c     BGJ - 9/98
c
c     $Id$
c
      implicit none
c
      integer ldim, nat
      integer nf1, nf3, nf4
      double precision scale
      double precision hess(3, nat, 3, nat)
      double precision coef1(nf1), d_kla(ldim,ldim)
c!!! Is 45 right since dummy atom is included?
      double precision d2int(78, nf4, nf3, nf1)
      integer idatom(4)
c
      integer i, icart, j, jcart, indi, indj, indab, indij
      integer ip1, ip3, ip4
      double precision d2E, den2
      double precision h(3,4,3,4)
c !!! Start BGJ temp test !!!
      double precision ti(3,4,3,4)
      character abcd(4), xyz(3)
      abcd(1) = 'A'
      abcd(2) = 'B'
      abcd(3) = 'C'
      abcd(4) = 'D'
      xyz(1) = 'x'
      xyz(2) = 'y'
      xyz(3) = 'z'
c !!! End BGJ temp test !!!
c
c      write(*,*)'*** d2int_2e3c_block_jhess called: idatom',idatom
      call dfill(144, 0.0d0, h, 1)
      indab = 0
      do i = 1, 4
         do j = i, 4
            if (idatom(j).ge.1 .and. idatom(i).ge.1) then
               indij = indab
               do icart = 1, 3
                  indi = 3*(i-1) + icart
                  do 100 jcart = 1, 3
                     indj = 3*(j-1) + jcart
                     if (indj .lt. indi) goto 100
                     indij = indij + 1
                     d2E = 0.D0
                     do ip1 = 1, nf1
                        do ip3 = 1, nf3
                           do ip4 = 1, nf4
c !!! Is factor of 1/2 right?
                              den2 = 0.5D0 * coef1(ip1) * d_kla(ip3,ip4)
                              d2E = d2E
     &                            + den2 * d2int(indij,ip4,ip3,ip1)
c       write(*,*) abcd(i), xyz(icart), abcd(j), xyz(jcart), indij,
c     &                                d2int(indij,ip4,ip3,ip2,ip1)
                           enddo
                        enddo
                     enddo
                     h(icart,i,jcart,j) = h(icart,i,jcart,j)
     &                    + d2E * scale
 100              continue
               enddo
            endif
            if (i.ne.j) then
               indab = indab + 9
            else
               indab = indab + 6
            endif
         enddo
      enddo
#if 0
c     Check TI of integral 2nd derivatives
      do ip1 = 1, nf1
          do ip3 = 1, nf3
            do ip4 = 1, nf4
              call dfill(144, 0.0d0, ti, 1)
              indab = 0
              do i = 1, 4
                do j = i, 4
                  if (idatom(j).ge.1 .and. idatom(i).ge.1) then
                    indij = indab
                    do icart = 1, 3
                      indi = 3*(i-1) + icart
                      do 200 jcart = 1, 3
                        indj = 3*(j-1) + jcart
                        if (indj .lt. indi) goto 200
                        indij = indij + 1
c       write(*,*) abcd(i), xyz(icart), abcd(j), xyz(jcart), indij,
c     &                                d2int(indij,ip4,ip3,ip1)
                        ti(icart,i,jcart,j)=d2int(indij,ip4,ip3,ip1)
                        ti(jcart,j,icart,i)=d2int(indij,ip4,ip3,ip1)
  200                 continue
                    enddo
                  endif
                  if (i.ne.j) then
                    indab = indab + 9
                  else
                    indab = indab + 6
                  endif
                enddo
              enddo
              write(*,*)'*** 2nd deriv TI check for',ip1,ip3,ip4
              call hessti(ti,4)
            enddo
          enddo
      enddo
#endif
c
c     Scatter the hessian contribution just formed to a full square
c
      do i = 1, 4
         do j = i, 4
            if (idatom(j).ge.1 .and. idatom(i).ge.1) then
               if (i.eq.j) then
                  do icart = 1, 2
                     do jcart = icart+1, 3
                        h(jcart,j,icart,i) = h(icart,i,jcart,j)
                     enddo
                  enddo
               else
                  do icart = 1, 3
                     do jcart = 1, 3
                        h(jcart,j,icart,i) = h(icart,i,jcart,j)
                     enddo
                  enddo
               endif
            endif
         enddo
      enddo
c
c     Increment the hessian
c
      do i = 1, 4
         do j = 1, 4
            if (idatom(j).ge.1 .and. idatom(i).ge.1) then
               do icart = 1, 3
                  do jcart = 1, 3
                     hess(icart,idatom(i),jcart,idatom(j)) =
     &               hess(icart,idatom(i),jcart,idatom(j)) +
     &                   h(icart,i,jcart,j)
                  enddo
               enddo
            endif
         enddo
      enddo
c      write(*,*)'*** TI check for hessian block'
c      call hessti(h,4)
c
      return
      end
c
c
c
      subroutine d2int_2e2c_block_jhess(coef1, coef2, hess, nat,
     &     d2int, scale, idatom, nf1, nf2)
c
c     Makes contribution to J hessian from block of second derivative
c     2-center integrals
c
c     BGJ - 9/98
c
c     $Id$
c
      implicit none
c
      integer nat
      integer nf1, nf2
      double precision scale
      double precision hess(3, nat, 3, nat)
      double precision coef1(nf1), coef2(nf2)
      double precision d2int(78, nf2, nf1)
      integer idatom(4)
c
      integer i, icart, j, jcart, indi, indj, indab, indij
      integer ip1, ip2
      double precision d2E, den2
      double precision h(3,4,3,4)
c !!! Start BGJ temp test !!!
      character abcd(4), xyz(3)
      abcd(1) = 'A'
      abcd(2) = 'B'
      abcd(3) = 'C'
      abcd(4) = 'D'
      xyz(1) = 'x'
      xyz(2) = 'y'
      xyz(3) = 'z'
c !!! End BGJ temp test !!!
c
c !!! Questions:
c     Should idatom be 4 long?
c     Will the "real" centers always be 1 and 3?
c      write(*,*)'*** d2int_2e2c_block_jhess called: idatom',idatom
      call dfill(144, 0.0d0, h, 1)
      indab = 0
      do i = 1, 4
         do j = i, 4
            if (idatom(j).ge.1 .and. idatom(i).ge.1) then
               indij = indab
               do icart = 1, 3
                  indi = 3*(i-1) + icart
                  do 100 jcart = 1, 3
                     indj = 3*(j-1) + jcart
                     if (indj .lt. indi) goto 100
                     indij = indij + 1
                     d2E = 0.D0

                     do ip1 = 1, nf1
                        do ip2 = 1, nf2
c !!! Is factor of 1/2 right?
                           den2 = 0.5D0 * coef1(ip1) * coef2(ip2)
                           d2E = d2E + den2 * d2int(indij,ip2,ip1)
c       write(*,*) abcd(i), xyz(icart), abcd(j), xyz(jcart), indij,
c     &                                d2int(indij,ip2,ip1)
                        enddo
                     enddo
                     h(icart,i,jcart,j) = h(icart,i,jcart,j)
     &                    + d2E * scale
 100              continue
               enddo
            endif
            if (i.ne.j) then
               indab = indab + 9
            else
               indab = indab + 6
            endif
         enddo
      enddo

#if 0
      call dfill(36, 0.0d0, h, 1)
      indab = 0
      do i = 1, 2
         do j = i, 2
            if (idatom(j).ge.1 .and. idatom(i).ge.1) then
               indij = indab
               do icart = 1, 3
                  indi = 3*(i-1) + icart
                  do 100 jcart = 1, 3
                     indj = 3*(j-1) + jcart
                     if (indj .lt. indi) goto 100
                     indij = indij + 1
                     d2E = 0.D0
                     do ip1 = 1, nf1
                        do ip2 = 1, nf2
                           den2 = 0.5D0 * coef1(ip1) * coef2(ip2)
                           d2E = d2E + den2 * d2int(indij,ip2,ip1)
c       write(*,*) abcd(i), xyz(icart), abcd(j), xyz(jcart), indij,
c     &                                d2int(indij,ip2,ip1)
                        enddo
                     enddo
                     h(icart,i,jcart,j) = h(icart,i,jcart,j)
     &                    + d2E * scale
 100              continue
               enddo
            endif
            if (i.ne.j) then
               indab = indab + 9
            else
               indab = indab + 6
            endif
         enddo
      enddo
#endif
c
c     Scatter the hessian contribution just formed to a full square
c
      do i = 1, 4
         do j = i, 4
            if (idatom(j).ge.1 .and. idatom(i).ge.1) then
               if (i.eq.j) then
                  do icart = 1, 2
                     do jcart = icart+1, 3
                        h(jcart,j,icart,i) = h(icart,i,jcart,j)
                     enddo
                  enddo
               else
                  do icart = 1, 3
                     do jcart = 1, 3
                        h(jcart,j,icart,i) = h(icart,i,jcart,j)
                     enddo
                  enddo
               endif
            endif
         enddo
      enddo
c
c     Increment the hessian
c
      do i = 1, 4
         do j = 1, 4
            if (idatom(j).ge.1 .and. idatom(i).ge.1) then
               do icart = 1, 3
                  do jcart = 1, 3
                     hess(icart,idatom(i),jcart,idatom(j)) =
     &               hess(icart,idatom(i),jcart,idatom(j)) +
     &                   h(icart,i,jcart,j)
                  enddo
               enddo
            endif
         enddo
      enddo
c      write(*,*)'*** TI check for hessian block'
c      call hessti(h,4)
c
      return
      end
