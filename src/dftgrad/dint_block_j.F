      subroutine dint_block_jfrc(d_ija, d_kla, ldim, frc, dint, scale,
     &     idatom, if1, il1, if2, il2, if3, il3, if4, il4)
c
c     Makes contribution to J force from block of derivative integrals
c
c     BGJ - 9/98
c
c     $Id: dint_block_j.F,v 1.1 1998-09-23 21:49:10 bjohnson Exp $
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
c      write(*,*)'*** dint_block_jfrc called'
      do i = 1, 4
         if (idatom(i).ge.1) then
            do icart = 1, 3
               dE = 0.D0
               do ip1 = if1, il1
                  do ip2 = if2, il2
                     do ip3 = if3, il3
                        do ip4 = if4, il4
c
                           den2 = 0.5D0 *
     &                            (d_ija(ip1,ip2) * d_kla(ip3,ip4))
c
                           dE = dE
     &                        + den2 * dint(ip4,ip3,ip2,ip1,icart,i)
c
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
     &     dint, scale, idatom, if1, il1, if2, il2, if3, il3, if4, il4)
c
c     Makes contribution to J derivative matrices from block of
c     derivative integrals
c
c     BGJ - 9/98
c
c     $Id: dint_block_j.F,v 1.1 1998-09-23 21:49:10 bjohnson Exp $
c
      implicit none
c
      integer ldim
      integer if1, il1, if2, il2, if3, il3, if4, il4
      double precision scale
      integer g_dj(3, *)
c !!! Is dimensioning right? !!!
      double precision dfij(if2:il2, if1:il1)
      double precision dfkl(if4:il4, if3:il3)
c !!! Is dimensioning right? !!!
      double precision d_ija(ldim,ldim), d_kla(ldim,ldim)
      double precision dint(if4:il4, if3:il3, if2:il2, if1:il1, 3, 4)
      integer idatom(4)
c
      integer i, icart
      integer ip1, ip2, ip3, ip4
      double precision g
c
c      write(*,*)'*** dint_block_jdmat called'
      do i = 1, 4
         if (idatom(i).ge.1) then
            do icart = 1, 3
               call dfill((il2-if2+1)*(il1-if1+1), 0.0d0, dfij, 1)
               call dfill((il4-if4+1)*(il3-if3+1), 0.0d0, dfkl, 1)
               do ip1 = if1, il1
                  do ip2 = if2, il2
                     do ip3 = if3, il3
                        do ip4 = if4, il4
                           g = dint(ip4,ip3,ip2,ip1,icart,i) * scale
                           dfij(ip2,ip1) = dfij(ip2,ip1)
     &                              + g * d_kla(ip3,ip4)
                           dfkl(ip4,ip3) = dfkl(ip4,ip3)
     &                              + g * d_ija(ip1,ip2)
                        enddo
                     enddo
                  enddo
               enddo
c
c     Update derivative J matrix blocks in GA's
c
c     !!! what about symmetrizing?
c     !!! what about across-diagonal?
c
c     The routine fock_upd_blk takes a temporary array for space
c     to rearrange blocks before updating.  This array is not used
c     in the case of a single matrix.  The call below relies upon
c     this and does not pass in a temporary array.  If the call
c     below is ever changed to update more than one matrix at once,
c     a temporary array will have to be provided.
c
               call fock_upd_blk(1, g_dj(icart,idatom(i)), if2, il2,
     $              if1, il1, 1.d0, dfij, 0.d0)
               call fock_upd_blk(1, g_dj(icart,idatom(i)), if4, il4,
     $              if3, il3, 1.d0, dfkl, 0.d0)
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
c     $Id: dint_block_j.F,v 1.1 1998-09-23 21:49:10 bjohnson Exp $
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
      integer i, icart, j, jcart, indi, indj, indij
      integer ip1, ip2, ip3, ip4
      double precision d2E, den2
c
      write(*,*)'*** d2int_block_jhess called'
      do i = 1, 4
         if (idatom(i).ge.1) then
            do icart = 1, 3
               indi = 3*(i-1) + icart
               do j = 1, i
                  if (idatom(j).ge.1) then
                     do 100 jcart = 1, 3
                        indj = 3*(j-1) + jcart
                        if (indj .gt. indi) goto 100
                        indij = indi*(indi-1)/2 + indj
                        d2E = 0.D0
                        do ip1 = if1, il1
                           do ip2 = if2, il2
                              do ip3 = if3, il3
                                 do ip4 = if4, il4
c
                                    den2 = 0.5D0 * (d_ija(ip1,ip2)
     &                                           *  d_kla(ip3,ip4))
c
                                    d2E = d2E + den2 *
     &                                   d2int(indij,ip4,ip3,ip2,ip1)
c
                                 enddo
                              enddo
                           enddo
                        enddo
c
                        hess(icart,idatom(i),jcart,idatom(j)) =
     &                  hess(icart,idatom(i),jcart,idatom(j)) +
     &                       + d2E * scale
                        if (idatom(i) .ne. idatom(j)) then
                           hess(jcart,idatom(j),icart,idatom(i)) =
     &                     hess(jcart,idatom(j),icart,idatom(i)) +
     &                          + d2E * scale
                        endif
c
 100                 continue
                  endif
               enddo 
            enddo
         endif
      enddo 
c
      return
      end
