      Subroutine xc_tabcd(tol_rho, Tmat, TTmat, Amat, Bmat, Cmat, Dmat, 
     &                    chi, delchi, ipol, nq, nbf, max_at_bf, 
     &                    max_at_bf2, GRAD, basis, natoms, iniz, ifin, 
     &                    g_vxc, ibf, rchi_atom, rdelchi_atom, 
     &                    rdens_atom, cetobfr)
c
C$Id: xc_tabcd.F,v 1.8 1997-12-08 23:48:54 d3h449 Exp $
c
      implicit none
c
      Logical GRAD
      integer basis
      integer max_at_bf ! [input]
      integer max_at_bf2 ! [input]
      integer ipol  ! [input]
      integer nq    ! [input]
      integer nbf    ! [input]
      integer natoms ! [input]
      double precision tol_rho
      integer g_vxc(*)   ! GA handle [input]
      double precision Tmat(max_at_bf2), TTmat(max_at_bf2)
      double precision rchi_atom(natoms)
      double precision rdelchi_atom(natoms)
      double precision rdens_atom(natoms,natoms,ipol)
      integer cetobfr(2,natoms)
c
c     Sampling Matrices for the XC Potential & Energy
c
      double precision Amat(nq,ipol), Cmat(nq,3,ipol)
c
c     Sampling Matrices for [Products of] Basis Functions & Gradients
c
      double precision Bmat(nq,max_at_bf)
      double precision Dmat(nq,3,max_at_bf)
      integer iniz(natoms), ifin(natoms)
c
c     Basis Functions & Gradients
c
      double precision chi(nq,nbf), delchi(nq,3,nbf)

      integer ibf(nbf)
      double precision A_MAX, C_MAX, AC_MAX, FUNC_MAXI,
     &                 B_MAX, D_MAX, BD_MAX, FUNC_MAXJ
      integer iat, inizia, ifinia, nbfia, nnia, ifirst, ilast
      integer jat, inizja, ifinja, nbfja, nnja, jfirst, jlast
      integer ii, mu, mu1
      integer n
      double precision chi1
      double precision dabsmax
      external dabsmax
c
c     Compute the matrix product for the XC potential and energy:
c
c              T = transpose(A*B) + transpose(C*D)
c
      A_MAX = dabsmax(nq*ipol,Amat)
      C_MAX = dabsmax(nq*3*ipol,Cmat)
      AC_MAX = max(A_MAX,C_MAX)
      do 430 iat = 1, natoms
         inizia = iniz(iat)
         if (inizia.eq.0)goto 430
         FUNC_MAXI = max(rchi_atom(iat),rdelchi_atom(iat))
         if (((AC_MAX*FUNC_MAXI).lt.tol_rho)) goto 430
         ifinia = ifin(iat)
         ifirst = cetobfr(1,iat)
         ilast = cetobfr(2,iat)
         nnia = ifinia - inizia + 1
         nbfia = ilast - ifirst + 1
         do 405 ii = 1, ipol
            do mu = 1, nnia
               mu1 = mu+inizia-1
               if (GRAD) then
                  do n = 1, nq
                     chi1 = chi(n,mu1)
                     Bmat(n,mu) = Amat(n,ii)*chi1 +
     &                            delchi(n,1,mu1)*Cmat(n,1,ii) +
     &                            delchi(n,2,mu1)*Cmat(n,2,ii) +
     &                            delchi(n,3,mu1)*Cmat(n,3,ii)
                     Dmat(n,1,mu) = Cmat(n,1,ii)*chi1
                     Dmat(n,2,mu) = Cmat(n,2,ii)*chi1
                     Dmat(n,3,mu) = Cmat(n,3,ii)*chi1
                  enddo
               else
                  do n = 1, nq
                     Bmat(n,mu) = chi(n,mu1)*Amat(n,ii)
                  enddo
               endif
            enddo
c
            B_MAX = dabsmax(nnia*nq,Bmat)
            D_MAX = dabsmax(nnia*nq*3,Dmat)
            BD_MAX = max(B_MAX,D_MAX)
c
            do 168 jat = 1, iat
               inizja = iniz(jat)
               if (inizja.eq.0)goto 168
               FUNC_MAXJ = max(rchi_atom(jat),rdelchi_atom(jat))
               if (BD_MAX*FUNC_MAXJ.lt.tol_rho) goto 168
               ifinja = ifin(jat)
               jfirst = cetobfr(1,jat)
               jlast = cetobfr(2,jat)
               nbfja = jlast - jfirst + 1
               nnja = ifinja - inizja + 1
               call dgemm('t', 'n', nnia, nnja, nq, 1.d0, Bmat, nq,
     &                    chi(1,inizja), nq, 0.d0, Tmat, nnia)
               if (GRAD) call dgemm('t', 'n', nnia, nnja, 3*nq, 1.d0,
     &                              Dmat, 3*nq, delchi(1,1,inizja), 
     &                              3*nq, 1.d0, Tmat, nnia)
               call scat_mat(TTmat, Tmat, nbfia, nbfja, nnia, nnja,
     &                       ifirst, jfirst, ibf(inizia), ibf(inizja))
               call upd_atom_block(g_vxc(ii), basis, iat, jat, TTmat)
c
c              check to see if can skip and use ga_symmetrize
c
               if (iat.ne.jat) then
                  call transp_mat(TTmat, Tmat, nbfia, nbfja)
                  call upd_atom_block(g_vxc(ii), basis, jat, iat, Tmat)
     &                                
               endif
  168       continue
  405    continue
  430 continue
      return
      end


