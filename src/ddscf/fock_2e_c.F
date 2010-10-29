      subroutine fock_2e_c( nfock, tol2e, q4_shell, eri,
     $     islo, ishi, jslo, jshi,
     $     kslo, kshi, lslo, lshi,
     $     ialo, iahi, jalo, jahi,
     $     kalo, kahi, lalo, lahi,
     $     dij, dik, dli, djk, dlj, dlk, 
     $     fij, fik, fli, fjk, flj, flk )
c     
c     $Id$

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
      integer i, j, k, l, jtop, ktop, ltop, v
      logical oij, okl, oikjl
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
              if (abs(eri(l,k,j,i)) .ge. tol2e) then
                g = eri(l,k,j,i)*q4_shell
                do v=1,nfock
                  fij(v,i,j) = fij(v,i,j) + g*dlk(v,l,k)
                  flk(v,l,k) = flk(v,l,k) + g*dij(v,i,j)
                  fik(v,i,k) = fik(v,i,k) + g*dlj(v,l,j)
                  fli(v,l,i) = fli(v,l,i) + g*djk(v,j,k)
                  flj(v,l,j) = flj(v,l,j) + g*dik(v,i,k)
                  fjk(v,j,k) = fjk(v,j,k) + g*dli(v,l,i)
                enddo
              end if
            enddo
          enddo
        enddo
      enddo
      
      end



c
c  Special code for nfock=1
c
      subroutine fock_2e_c1( tol2e, q4_shell, eri,
     $                       islo, ishi, jslo, jshi,
     $                       kslo, kshi, lslo, lshi,
     $                       ialo, iahi, jalo, jahi,
     $                       kalo, kahi, lalo, lahi,
     $                       dij, dik, dli, djk, dlj, dlk, 
     $                       fij, fik, fli, fjk, flj, flk )
c     
c     $Id$

      implicit none
c     
c     Add integrals without labels into atomic fock matrix blocks ... 
c     symmetrization will be needed later
c
c     !! THIS ROUTINE IS OBSOLETED BY ROUTINES THAT USE LABELS
c     
      double precision tol2e
      double precision q4_shell
      integer islo, ishi, jslo, jshi, kslo, kshi, lslo, lshi
      integer ialo, iahi, jalo, jahi, kalo, kahi, lalo, lahi
      double precision eri(lslo:lshi,kslo:kshi,jslo:jshi,islo:ishi)
      double precision dij(ialo:iahi,jalo:jahi)
      double precision dik(ialo:iahi,kalo:kahi)
      double precision dli(lalo:lahi,ialo:iahi)
      double precision djk(jalo:jahi,kalo:kahi)
      double precision dlj(lalo:lahi,jalo:jahi)
      double precision dlk(lalo:lahi,kalo:kahi)
      double precision fij(ialo:iahi,jalo:jahi)
      double precision fik(ialo:iahi,kalo:kahi)
      double precision fli(lalo:lahi,ialo:iahi)
      double precision fjk(jalo:jahi,kalo:kahi)
      double precision flj(lalo:lahi,jalo:jahi)
      double precision flk(lalo:lahi,kalo:kahi)
c     
      integer i, j, k, l, jtop, ktop, ltop
      logical oij, okl, oikjl
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
      do i = islo, ishi
        jtop = jshi
        if (oij) jtop = i
        do j = jslo, jtop
          d_ij = dij(i,j)
          f_ij = 0.d0
          ktop = kshi
          if (oikjl) ktop = i
          do k = kslo, ktop
            d_jk = djk(j,k)
            d_ik = dik(i,k)
            f_jk = 0.d0
            f_ik = 0.d0
            ltop = lshi
            if (okl) ltop = k
            if (oikjl .and. k.eq.i) ltop = j
            do l = lslo, ltop
              if (abs(eri(l,k,j,i)) .ge. tol2e) then
                g = eri(l,k,j,i)*q4_shell
                flk(l,k)   = flk(l,k)   + g * d_ij
                fli(l,i)   = fli(l,i)   + g * d_jk
                flj(l,j)   = flj(l,j)   + g * d_ik
                f_ij       = f_ij       + g * dlk(l,k)
                f_ik       = f_ik       + g * dlj(l,j)
                f_jk       = f_jk       + g * dli(l,i)
              end if
            enddo
            fjk(j,k) = fjk(j,k) + f_jk
            fik(i,k) = fik(i,k) + f_ik
          enddo
          fij(i,j) = fij(i,j) + f_ij
        enddo
      enddo
      
      end










