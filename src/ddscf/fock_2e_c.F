c     
c  Add integrals into atomic fock matrix blocks ... 
c  symmetrization will be needed later
c
      subroutine fock_2e_c( nfock, tol2e, q4_shell, eri,
     $                      islo, ishi, jslo, jshi,
     $                      kslo, kshi, lslo, lshi,
     $                      ialo, iahi, jalo, jahi,
     $                      kalo, kahi, lalo, lahi,
     $                      dij, dik, dli, djk, dlj, dlk, 
     $                      fij, fik, fli, fjk, flj, flk )
c
c$Id: fock_2e_c.F,v 1.1 1995-03-14 01:52:05 rg240 Exp $
c
      implicit none
c     
      integer nfock
      double precision tol2e
      double precision q4_shell
      double precision eri(lslo:lshi,kslo:kshi,jslo:jshi,islo:ishi)
      integer islo, ishi, jslo, jshi, kslo, kshi, lslo, lshi
      integer ialo, iahi, jalo, jahi, kalo, kahi, lalo, lahi
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
      double precision g
c
c  Assume identical index ranges ==> identical shells
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
              g = eri(l,k,j,i)*q4_shell
              if (abs(g) .ge. tol2e) then
                do v=1,nfock
                  fij(v,i,j) = fij(v,i,j) + g*dlk(v,l,k)
                  flk(v,l,k) = flk(v,l,k) + g*dij(v,i,j)
                  fik(v,i,k) = fik(v,i,k) + g*dlj(v,l,j)
                  fli(v,l,i) = fli(v,l,i) + g*djk(v,j,k)
                  flj(v,l,j) = flj(v,l,j) + g*dik(v,i,k)
                  fjk(v,j,k) = fjk(v,j,k) + g*dli(v,l,i)
                enddo
              endif
            enddo
          enddo
        enddo
      enddo

      end






