      subroutine fock_2e_c( nfock, tol2e, q4_shell, eri,
     $     islo, ishi, jslo, jshi,
     $     kslo, kshi, lslo, lshi,
     $     ialo, iahi, jalo, jahi,
     $     kalo, kahi, lalo, lahi,
     $     dij, dik, dli, djk, dlj, dlk, 
     $     fij, fik, fli, fjk, flj, flk )
c     
c     $Id: fock_2e_c.F,v 1.5 1995-10-21 19:39:18 d3g681 Exp $

      implicit none
c     
c     Add integrals without labels into atomic fock matrix blocks ... 
c     symmetrization will be needed later
c
c     !! THIS ROUTINE IS OBSOLETED BY ROUTINES THAT USE LABELS
c     
      integer nfock
      double precision tol2e
      double precision q4_shell
      integer islo, ishi, jslo, jshi, kslo, kshi, lslo, lshi
      integer ialo, iahi, jalo, jahi, kalo, kahi, lalo, lahi
      double precision eri(lslo:lshi,kslo:kshi,jslo:jshi,islo:ishi)
      double precision dij(nfock,ialo:iahi,jalo:jahi)
      double precision dik(nfock,ialo:iahi,kalo:kahi)
      double precision dli(nfock,lalo:lahi,ialo:iahi)
      double precision djk(nfock,jalo:jahi,kalo:kahi)
      double precision dlj(nfock,lalo:lahi,jalo:jahi)
      double precision dlk(nfock,lalo:lahi,kalo:kahi)
      double precision fij(nfock,ialo:iahi,jalo:jahi)
      double precision fik(nfock,ialo:iahi,kalo:kahi)
      double precision fli(nfock,lalo:lahi,ialo:iahi)
      double precision fjk(nfock,jalo:jahi,kalo:kahi)
      double precision flj(nfock,lalo:lahi,jalo:jahi)
      double precision flk(nfock,lalo:lahi,kalo:kahi)
c     
      integer i, j, k, l, jtop, ktop, ltop
      logical oij, okl, oikjl
      integer v
      double precision d_ij, d_jk, d_ik
      double precision f_ij, f_jk, f_ik
      double precision g
c     
c     Assume identical index ranges ==> identical shells
c     
      oij = islo.eq.jslo
      okl = kslo.eq.lslo
      oikjl = (islo.eq.kslo).and.(jslo.eq.lslo)
c     
c     Add into the fock matrix
c     
      if (nfock.gt.1) then
         do i = islo, ishi
            jtop = jshi
            if (oij) jtop = i
            do j = jslo, jtop
               ktop = kshi
               if (oikjl) ktop = i
               do k = kslo, ktop
                  ltop = lshi
                  if (okl) ltop = k
                  if (oikjl .and. k.eq.i) ltop = j
                  do l = lslo, ltop
                     g = eri(l,k,j,i)*q4_shell
                     if (abs(g) .ge. tol2e) then
                        do v=1,nfock
                           fij(v,i,j) = fij(v,i,j) + g*dlk(v,l,k)
                           flk(v,l,k) = flk(v,l,k) + g*dij(v,i,j)
                           fik(v,i,k) = fik(v,i,k) + g*dlj(v,l,j)
                           fli(v,l,i) = fli(v,l,i) + g*djk(v,j,k)
                           flj(v,l,j) = flj(v,l,j) + g*dik(v,i,k)
                           fjk(v,j,k) = fjk(v,j,k) + g*dli(v,l,i)
                        end do
                     end if
                  end do
               end do
            end do
         end do
      else
         do i = islo, ishi
            jtop = jshi
            if (oij) jtop = i
            do j = jslo, jtop
               d_ij = dij(1,i,j)
               f_ij = 0.d0
               ktop = kshi
               if (oikjl) ktop = i
               do k = kslo, ktop
                  d_jk = djk(1,j,k)
                  f_jk = 0.d0
                  d_ik = dik(1,i,k)
                  f_ik = 0.d0
                  ltop = lshi
                  if (okl) ltop = k
                  if (oikjl .and. k.eq.i) ltop = j
                  do l = lslo, ltop
                     g = eri(l,k,j,i)*q4_shell
                     if (abs(g) .ge. tol2e) then
                        f_ij       = f_ij       + g * dlk(1,l,k)
                        flk(1,l,k) = flk(1,l,k) + g * d_ij
                        f_ik       = f_ik       + g * dlj(1,l,j)
                        fli(1,l,i) = fli(1,l,i) + g * d_jk
                        flj(1,l,j) = flj(1,l,j) + g * d_ik
                        f_jk       = f_jk       + g * dli(1,l,i)
                     end if
                  end do
                  fjk(1,j,k) = fjk(1,j,k) + f_jk
                  fik(1,i,k) = fik(1,i,k) + f_ik
               end do
               fij(1,i,j) = fij(1,i,j) + f_ij
            end do
         end do
      end if

      end






