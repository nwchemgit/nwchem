c
c     Computes explicit nuclear 2nd derivatives of the XC energy
c
c     BGJ - 8/98
c
c     $Id$
c 
C> \ingroup hess
C> @{
C>
C> \file xc_d2expl.F
C> Explicit 2nd derivatives wrt nuclear coordinates of the DFT 
C> functional
C>
C> \brief Computes explicit nuclear 2nd derivatives of the XC energy
C>
C> The explicit 2nd derivatives of the XC energy with respect to the
C> nuclear coordinates are defined as
C> \f{eqnarray*}{
C>   \frac{\mathrm{d}^2E_{xc}}{\mathrm{d}x\mathrm{d}y} &=& \sum_{ij}\int
C>   \frac{\mathrm{d}^2f_{xc}}{\mathrm{d}p_i\mathrm{d}p_j}
C>   \frac{\mathrm{d}p_i}{\mathrm{d}x}\frac{\mathrm{d}p_j}{\mathrm{d}y}
C>   \mathrm{d}r
C> + \sum_i\int\frac{\mathrm{d}f_{xc}}{\mathrm{d}p_i}
C>   \frac{\mathrm{d}^2p_i}{\mathrm{d}x\mathrm{d}y}
C>   \mathrm{d}r
C> \f}
C> where \f$p \in \{\rho_\alpha,\rho_\beta,\gamma_{\alpha\alpha},
C> \gamma_{\alpha\beta},\gamma_{\beta\beta}\}\f$ and the indeces 
C> \f$i\f$ and \f$j\f$ run over this set.
c
c     d2Exc             /    d2fxc    dp(i) dp(j)          / dfxc  d2p(i)
c     ----- =  sum sum  | ----------- ----- -----  +  sum  | ----- ------
c     dx dy     i   j  /  dp(i) dp(j)  dx    dy        i  /  dp(i) dx dy
c
c     where "p" represents a density parameter in the set
c     { ra, rb, gaa, gab, gbb }, and the i, j indices run over this set
c
c     First step: separable part of 2nd derivative contribution
c
      Subroutine xc_d2expl(tol_rho, scr, 
     &     Amat, Amat2, Acof2, Cmat, Cmat2, Ccof2, Mmat, Mmat2, Mcof2, 
     &     F, Pmat, ff, s, 
     &     chi, delchi, heschi, d3chi,
     &     curatoms, ncuratoms, ipol, nq, nbf, max_at_bf, GRAD, basis,
     &     natoms, iniz, ifin, drho, ddelrho, dttau, 
     &     delrho, g_dens, hess, ibf,
     &     rchi_atom, rdelchi_atom, rdens_atom, cetobfr, kske)
c
      implicit none
#include "dft2drv.fh"
#include "stdio.fh"
#include "util.fh"
#include "global.fh"
c
      logical kske !< [Input] .True. if the functional is a meta-GGA
      integer basis !< [Input] The basis set handle
      integer max_at_bf !< [Input] The maximum number of basis functions
                        !< on any atom
      integer ipol  !< [Input] The number of spin channels
      integer nq    !< [Input] The number of grid points
      integer nbf    !< [Input] The number of basis functions
      integer natoms !< [Input] The number of atoms
      integer ncuratoms !< [Input] The number of currently "active"
                        !< atoms
      integer curatoms(natoms) !< [Input] The mapping array for 
                               !< currently active atoms
      double precision tol_rho !< [Input] The tolerance on the density
      logical GRAD !< [Input] .True. if the functional is a GGA
c
c     Explicit first derivatives of density wrt current nuclei [input]
c
      double precision drho(nq,ipol,3,ncuratoms) !< [Input] Derivative
      !< of rho wrt nuclear coordinates 
      double precision ddelrho(nq,3,ipol,3,ncuratoms) !< [Input]
      !< Derivative of the density gradient wrt nuclear coordinates
      double precision dttau(nq,ipol,3,ncuratoms) !< [Input] Derivative 
      !< of the kinetic energy density wrt nuclear coordinates
c
c     Space for separable coefficients of first derivatives of density
c
      double precision Acof2(nq,ipol,3) !< [Scratch] Derivative rho
      double precision Ccof2(nq,3,ipol,3) !< [Scratch] Derivative grad
      double precision Mcof2(nq,ipol,3) !< [Scratch] Derivative tau
c
c     Spin density gradients
c
      double precision delrho(nq,3,ipol) !< [Input] Density gradient
c
      integer g_dens(ipol) !< [Input] GA handle for density matrices
c
c     Hessian matrix (updated)
c
      double precision hess(3,natoms,3,natoms) !< [In/Output] Hessian
c
      double precision scr(nq,15) !< [Scratch] matrix
      double precision rchi_atom(natoms) !< [Input] Screening parameters
      double precision rdelchi_atom(natoms) !< [Input] Screening
                                            !< parameters
      double precision rdens_atom(natoms,natoms,ipol) !< [Input]
      !< Screening parameters
      integer cetobfr(2,natoms) !< [Input] First and last basis 
      !< function of the atoms
c
      double precision Pmat(max_at_bf*max_at_bf) !< [Scratch] vector
      double precision F(max_at_bf*max_at_bf)    !< [Scratch] vector
      double precision ff(nq,3,*) !< [Scratch] arrays
      double precision s(nq,max_at_bf) !< [Scratch] arrays
c
c     Sampling Matrices for the XC Functional Derivatives
c
      double precision Amat(nq,ipol) !< [Input] Derivative of functional
      !< wrt rho
      double precision Cmat(nq,3,ipol) !< [Input] Derivative of the
      !< functional wrt the norm of the electron density gradient.
      !< The call to `transform_Cmat` changes the contents to 
      !< \f$\gamma_{\alpha\alpha}\cdot\frac{\partial \rho_\alpha}{\partial x}
      !< \ldots \gamma_{\beta\beta}\cdot\frac{\partial \rho_\beta}{\partial z}\f$
      double precision Mmat(nq,ipol) !< [Input] Derivative of functional
      !< wrt kinetic energy density
c
      double precision Amat2(nq,NCOL_AMAT2) !< [Input] 2nd derivative
      !< of functional wrt rho
      double precision Cmat2(nq,NCOL_CMAT2) !< [Input] 2nd derivative
      !< of functional wrt gamma
      double precision Mmat2(nq,NCOL_MMAT2) !< [Input] 2nd derivative
      !< of functional wrt kinetic energy density
c
c     Sampling Matrices for [Products of] Basis Functions & Gradients
c
      integer iniz(natoms) !< [Input] Start something
      integer ifin(natoms) !< [Input] End something
c
c     Basis Functions & Derivatives
c
      double precision chi(nq,nbf) !< [Input] Basis function values
      double precision delchi(nq,3,nbf) !< [Input] Basis function 
                                        !< derivative values
      double precision heschi(nq,6,nbf) !< [Input] Basis function 
                                        !< 2nd derivative values
      double precision d3chi(nq,10,nbf) !< [Input] Basis function 
                                        !< 3rd derivative values
c
      integer ibf(nbf) !< [Input] Some mapping table
c
c     local declarations
c
      logical oprint
      double precision A_MAX, C_MAX, AC_MAX, FUNC_MAXI, FUNC_MAXJ
      double precision FUNC_MAX, DELFUNC_MAX, tol_rho_tmp
      integer iatcur, jatcur
      integer iat, inizia, ifinia, nbfia, nnia, ifirst, ilast, idim
      integer jat, inizja, ifinja, nbfja, nnja, jfirst, jlast, jdim
      integer ii, mu, nu, icount
      integer n
      double precision aaa, fdchix, fdchiy, fdchiz,
     &                 ccc1, ccc2, ccc3
      double precision T(3,3),t6(6)
      integer idir, jdir
c
c     The following parameter definitions must be consistent with
c     those in routine xc_eval_basis, or this routine will not work
c
      integer iixx,iixy,iixz,
     &             iiyy,iiyz,
     &                  iizz
c
      parameter ( iixx=1,iixy=2,iixz=3,
     &                   iiyy=4,iiyz=5,
     &                          iizz=6 )
c
      integer iixxx,iixxy,iixxz,
     &              iixyy,iixyz,
     &                    iixzz,
     &                          iiyyy,iiyyz,
     &                                iiyzz,
     &                                      iizzz
c
      parameter ( iixxx=1,iixxy=2,iixxz=3,
     &                    iixyy=4,iixyz=5,
     &                            iixzz=6,
     &                                    iiyyy=7,iiyyz=8,
     &                                            iiyzz=9,
     &                                                    iizzz=10 )
c
      double precision duefac
      double precision dabsmax
      external dabsmax
c
c     d2Exc             /    d2fxc    dp(i) dp(j)          / dfxc  d2p(i)
c     ----- =  sum sum  | ----------- ----- -----  +  sum  | ----- ------
c     dx dy     i   j  /  dp(i) dp(j)  dx    dy        i  /  dp(i) dx dy
c
c     where "p" represents a density parameter in the set
c     { ra, rb, gaa, gab, gbb }, and the i, j indices run over this set
c
c     First step: separable part of 2nd derivative contribution
c
      oprint= util_print('xc_hessian',print_debug)
      do 10 iat = 1, natoms
         iatcur = curatoms(iat)
         if (iatcur.eq.0) goto 10
c
c     Form Acof2, Ccof2 for iatcur
c
         call dcopy(nq*ipol*3,drho(1,1,1,iatcur),1,Acof2,1)
         if (grad) then
            call dcopy(nq*ipol*9,ddelrho(1,1,1,1,iatcur),1,Ccof2,1)
         endif
         if (kske) then
            call dcopy(nq*ipol*3,dttau(1,1,1,iatcur),1,Mcof2,1)
         end if
c
         call xc_cpks_coeff(Acof2, Ccof2, Mcof2,
     &        Amat2, Cmat2, Cmat, Mmat2,
     &        delrho,3, ipol, nq, grad, kske, .false.)
c
         do 20 jat = 1, iat
            jatcur = curatoms(jat)
            if (jatcur.eq.0) goto 20
c
            call dfill(9,0.d0,T,1)
            if (ipol.eq.1) then
               if (.not.GRAD) then
                  do jdir = 1, 3
                     do idir = 1, 3
                     T(idir,jdir) = T(idir,jdir) +
     &                    ddot(nq,Acof2(1,1,idir),1,
     .                    drho(1,1,jdir,jatcur),1)
                     enddo
                  enddo
               else
                  do jdir = 1, 3
                     do idir = 1, 3
                        do n = 1, nq
                           T(idir,jdir) = T(idir,jdir)
     &                  + Acof2(n,1,idir)*drho(n,1,jdir,jatcur)
     &                  + Ccof2(n,1,1,idir)*ddelrho(n,1,1,jdir,jatcur)
     &                  + Ccof2(n,2,1,idir)*ddelrho(n,2,1,jdir,jatcur)
     &                  + Ccof2(n,3,1,idir)*ddelrho(n,3,1,jdir,jatcur)
                        enddo
                     enddo
                  enddo
               endif
            else
               do jdir = 1, 3
                  do idir = 1, 3
                     T(idir,jdir) = T(idir,jdir) +
     &                    ddot(nq,Acof2(1,1,idir),1,
     .                    drho(1,1,jdir,jatcur),1) +
     &                    ddot(nq,Acof2(1,2,idir),1,
     .                    drho(1,2,jdir,jatcur),1)
                  enddo
               enddo
               if (GRAD) then
                  do jdir = 1, 3
                     do idir = 1, 3
                        do n = 1, nq
                           T(idir,jdir) = T(idir,jdir)
     &                  + Ccof2(n,1,1,idir)*ddelrho(n,1,1,jdir,jatcur)
     &                  + Ccof2(n,2,1,idir)*ddelrho(n,2,1,jdir,jatcur)
     &                  + Ccof2(n,3,1,idir)*ddelrho(n,3,1,jdir,jatcur)
     &                  + Ccof2(n,1,2,idir)*ddelrho(n,1,2,jdir,jatcur)
     &                  + Ccof2(n,2,2,idir)*ddelrho(n,2,2,jdir,jatcur)
     &                  + Ccof2(n,3,2,idir)*ddelrho(n,3,2,jdir,jatcur)
                        enddo
                     enddo
                  enddo
               endif
            endif
c
c           Update Hessian block(s)
c
            do jdir = 1,3
               do idir = 1,3
                  hess(idir,iat,jdir,jat) = hess(idir,iat,jdir,jat)
     &                                    + T(idir,jdir)
                  if (iat.ne.jat) then
                     hess(jdir,jat,idir,iat) = hess(jdir,jat,idir,iat)
     &                                       + T(idir,jdir)
                  endif
               enddo
            enddo
c
 20      continue
 10   continue
c
c     Second step: remaining terms involving functional first derivatives
c                  and density parameter second derivatives
c
c     We now need Cmat in the delrho form
c
      if (GRAD) call transform_Cmat(delrho, Cmat, ipol, nq)
c
      A_MAX = dabsmax(nq*ipol,Amat)
      C_MAX = dabsmax(nq*3*ipol,Cmat)
      AC_MAX = max(A_MAX,C_MAX)
c
#if 0
      write(6,*) ' xc_d2expl: AMAT '
      call output(amat, 1, nq, 1, ipol, nq, ipol, 1)
      if (GRAD) then
         write(6,*) ' xc_d2expl: CMAT '
         call output(cmat, 1, nq, 1, 3*ipol, nq, 3*ipol, 1)
      endif
      write(6,*) ' xc_d2expl: chi '
      call output(chi, 1, nq, 1, nbf, nq, nbf, 1)
      if (GRAD) then
         write(6,*) ' xc_d2expl: delchi '
         call output(delchi, 1, nq, 1, 3*nbf, nq, 3*nbf, 1)
      endif
#endif
c
c     Screening is accomplished by:  p(r) <= |Xi(r)|*|Xj(r)|*|Dij|
c     Xi(r) is screened on desired accuracy/max(|Xj(r)|)*max(|Dij|)
c     Dij is screened on desired accuracy/max(|Xi(r)|)*max(|Xj(r)|)
c
      FUNC_MAX = dabsmax(natoms,rchi_atom)
      DELFUNC_MAX = dabsmax(natoms,rdelchi_atom)
c
      do 230 iat = 1, natoms
         inizia = iniz(iat)
         if (inizia.eq.0)goto 230
         iatcur = curatoms(iat)
         ifinia = ifin(iat)
         ifirst = cetobfr(1,iat)
         ilast = cetobfr(2,iat)
         nbfia = ilast-ifirst+1
         nnia = ifinia-inizia+1
c
c        screening parameters
c
         FUNC_MAXI = max(rchi_atom(iat),rdelchi_atom(iat))
         FUNC_MAXJ = max(FUNC_MAX,DELFUNC_MAX)
#if 0
         if (ipol.gt.1)then
            P_MAXJ_A = dabsmax(natoms,rdens_atom(1,iat,1))
            P_MAXJ_B = dabsmax(natoms,rdens_atom(1,iat,2))
            P_MAXJ = MAX(P_MAXJ_A, P_MAXJ_B)
         else
            P_MAXJ = dabsmax(natoms,rdens_atom(1,iat,1))
         endif
         if (FUNC_MAXI*FUNC_MAXJ*P_MAXJ.lt.tol_rho) goto 225
c     !!! Cutoff temporarily commented out !!!
#endif
         do 220 jat = 1, iat
            inizja = iniz(jat)
            if (inizja.eq.0)goto 220
            jatcur = curatoms(jat)
            if (iatcur .eq. 0 .and. jatcur .eq. 0) goto 220
            ifinja = ifin(jat)
            jfirst = cetobfr(1,jat)
            jlast = cetobfr(2,jat)
            nbfja = jlast-jfirst+1
            nnja = ifinja-inizja+1
c
c           screening parameters
c
            FUNC_MAXJ = max(rchi_atom(jat),rdelchi_atom(jat))
#if 0
            if (ipol.eq.1)then
               P_MAXIJ = rdens_atom(iat,jat,1)
            else
               P_MAXIJ = max(rdens_atom(iat,jat,1),
     &                       rdens_atom(iat,jat,2))
            endif
            if (FUNC_MAXI*FUNC_MAXJ*P_MAXIJ.lt.tol_rho) goto 215
c     !!! Cutoff temporarily commented out !!!
#endif
            tol_rho_tmp = tol_rho/(FUNC_MAXI*FUNC_MAXJ)
c
            do 210 ii = 1, ipol
c
c              screening parameters
c
#if 0
               P_MAXIJ = rdens_atom(iat,jat,ii)
               if (FUNC_MAXI*FUNC_MAXJ*P_MAXIJ.lt.tol_rho)goto 210
c     !!! Cutoff temporarily commented out !!!
#endif
c
               call get_atom_block(g_dens(ii), basis,
     &                             iat, jat, Pmat, idim, jdim)
c
               call gat_mat(F, Pmat, nbfia, nbfja, nnia, nnja,ifirst, 
     &                      jfirst, ibf(inizia), ibf(inizja))
c
c     Three terms to compute
c
c     First term:    Xiat(r)*hessXjat(r)*Diat,jat -> hess(jat,jat)
c             GC: delXiat(r)*hessXjat(r)*Diat,jat -> hess(jat,jat)
c             GC:    Xiat(r)*  d3Xjat(r)*Diat,jat -> hess(jat,jat)
c
               if (jatcur .ne. 0) then
                  call dfill(9,0.d0,T,1)
                  call dgemm('n','n',nq,nnja,nnia,1d0,
     A                 chi(1,inizia),nq,F,nnia,0d0,s,nq)
                  if (grad)
     G                 call dgemm('n','n',nq*3,nnja,nnia,1d0,
     A                 delchi(1,1,inizia),nq*3,F,nnia,0d0,ff,nq*3)
                  t6(1)=0d0
                  t6(2)=0d0
                  t6(3)=0d0
                  t6(4)=0d0
                  t6(5)=0d0
                  t6(6)=0d0
                  do mu=inizja,ifinja
                     if (GRAD) then
                        do n = 1, nq
                           ff(n,1,1) = 
     S                          Amat(n,ii)*s(n,mu-inizja+1)
     &                          + Cmat(n,1,ii)*ff(n,1,mu-inizja+1)
     &                          + Cmat(n,2,ii)*ff(n,2,mu-inizja+1)
     &                          + Cmat(n,3,ii)*ff(n,3,mu-inizja+1)
                        enddo
                     else
                        do n = 1, nq
                           ff(n,1,1) = Amat(n,ii)*
     *                          s(n,mu-inizja+1)
                        enddo
                     endif
                     call dgemv('t',nq, 6, 1d0,heschi(1,1,mu),nq,
     Y                    ff(1,1,1),1,1d0,t6,1)
                  enddo
                  T(1,1) = T(1,1)+t6(1)
                  T(1,2) = T(1,2)+t6(2)
                  T(1,3) = T(1,3)+t6(3)
                  T(2,2) = T(2,2)+t6(4)
                  T(2,3) = T(2,3)+t6(5)
                  T(3,3) = T(3,3)+t6(6)
                  if (GRAD) then
                     do mu = inizja, ifinja
                        do n = 1, nq
                           ccc1 = Cmat(n,1,ii)*s(n,mu-inizja+1)
                           ccc2 = Cmat(n,2,ii)*s(n,mu-inizja+1)
                           ccc3 = Cmat(n,3,ii)*s(n,mu-inizja+1)
                           T(1,1) = T(1,1) + d3chi(n,iixxx,mu)*ccc1
     &                                     + d3chi(n,iixxy,mu)*ccc2
     &                                     + d3chi(n,iixxz,mu)*ccc3
                           T(1,2) = T(1,2) + d3chi(n,iixxy,mu)*ccc1
     &                                     + d3chi(n,iixyy,mu)*ccc2
     &                                     + d3chi(n,iixyz,mu)*ccc3
                           T(1,3) = T(1,3) + d3chi(n,iixxz,mu)*ccc1
     &                                     + d3chi(n,iixyz,mu)*ccc2
     &                                     + d3chi(n,iixzz,mu)*ccc3
                           T(2,2) = T(2,2) + d3chi(n,iixyy,mu)*ccc1
     &                                     + d3chi(n,iiyyy,mu)*ccc2
     &                                     + d3chi(n,iiyyz,mu)*ccc3
                           T(2,3) = T(2,3) + d3chi(n,iixyz,mu)*ccc1
     &                                     + d3chi(n,iiyyz,mu)*ccc2
     &                                     + d3chi(n,iiyzz,mu)*ccc3
                           T(3,3) = T(3,3) + d3chi(n,iixzz,mu)*ccc1
     &                                     + d3chi(n,iiyzz,mu)*ccc2
     &                                     + d3chi(n,iizzz,mu)*ccc3
                        enddo
                     enddo
                  endif
c     
                  duefac=1d0
                  if (iat.ne.jat) duefac=2d0
c
                  T(2,1) = T(1,2)
                  T(3,1) = T(1,3)
                  T(3,2) = T(2,3)
                  do jdir = 1,3
                     do idir = 1,3
                        hess(idir,jat,jdir,jat) =
     &                  hess(idir,jat,jdir,jat) + T(idir,jdir)*duefac
                     enddo
                  enddo
               endif
c
c     Second term: hessXiat(r)*   Xjat(r)*Diat,jat -> hess(iat,iat)
c              GC: hessXiat(r)*delXjat(r)*Diat,jat -> hess(iat,iat)
c              GC:   d3Xiat(r)*   Xjat(r)*Diat,jat -> hess(iat,iat)
c
               if (iatcur .ne. 0) then
                  call dfill(9,0.d0,T,1)
                  call dgemm('n','t',nq,nnia,nnja,1d0,
     A                 chi(1,inizja),nq,F,nnia,0d0,s,nq)
                  if (grad)
     G                 call dgemm('n','t',nq*3,nnia,nnja,1d0,
     A                 delchi(1,1,inizja),nq*3,F,nnia,0d0,ff,nq*3)
                  t6(1)=0d0
                  t6(2)=0d0
                  t6(3)=0d0
                  t6(4)=0d0
                  t6(5)=0d0
                  t6(6)=0d0
                  do nu=inizia,ifinia
                     if (GRAD) then
                        do n = 1, nq
                           ff(n,1,1) = 
     S                          Amat(n,ii)*s(n,nu-inizia+1)
     &                             + Cmat(n,1,ii)*ff(n,1,nu-inizia+1)
     &                             + Cmat(n,2,ii)*ff(n,2,nu-inizia+1)
     &                             + Cmat(n,3,ii)*ff(n,3,nu-inizia+1)
                        enddo
                     else
                        do n = 1, nq
                           ff(n,1,1) = 
     A                          Amat(n,ii)*s(n,nu-inizia+1)
                        enddo
                     endif
                     call dgemv('t',nq, 6, 1d0,heschi(1,1,nu),nq,
     Y                    ff(1,1,1),1,1d0,t6,1)
                  enddo
                  T(1,1) = T(1,1)+t6(1)
                  T(1,2) = T(1,2)+t6(2)
                  T(1,3) = T(1,3)+t6(3)
                  T(2,2) = T(2,2)+t6(4)
                  T(2,3) = T(2,3)+t6(5)
                  T(3,3) = T(3,3)+t6(6)
                  if (GRAD) then
                     do nu = inizia, ifinia
c     
                        do n = 1, nq
                           ccc1 = Cmat(n,1,ii)*s(n,nu-inizia+1)
                           ccc2 = Cmat(n,2,ii)*s(n,nu-inizia+1)
                           ccc3 = Cmat(n,3,ii)*s(n,nu-inizia+1)
                           T(1,1) = T(1,1) + d3chi(n,iixxx,nu)*ccc1
     &                                     + d3chi(n,iixxy,nu)*ccc2
     &                                     + d3chi(n,iixxz,nu)*ccc3
                           T(1,2) = T(1,2) + d3chi(n,iixxy,nu)*ccc1
     &                                     + d3chi(n,iixyy,nu)*ccc2
     &                                     + d3chi(n,iixyz,nu)*ccc3
                           T(1,3) = T(1,3) + d3chi(n,iixxz,nu)*ccc1
     &                                     + d3chi(n,iixyz,nu)*ccc2
     &                                     + d3chi(n,iixzz,nu)*ccc3
                           T(2,2) = T(2,2) + d3chi(n,iixyy,nu)*ccc1
     &                                     + d3chi(n,iiyyy,nu)*ccc2
     &                                     + d3chi(n,iiyyz,nu)*ccc3
                           T(2,3) = T(2,3) + d3chi(n,iixyz,nu)*ccc1
     &                                     + d3chi(n,iiyyz,nu)*ccc2
     &                                     + d3chi(n,iiyzz,nu)*ccc3
                           T(3,3) = T(3,3) + d3chi(n,iixzz,nu)*ccc1
     &                                     + d3chi(n,iiyzz,nu)*ccc2
     &                                     + d3chi(n,iizzz,nu)*ccc3
                        enddo
                     enddo
                  endif
c
                  duefac=1d0
                  if (iat.ne.jat) duefac=2d0
c
                  T(2,1) = T(1,2)
                  T(3,1) = T(1,3)
                  T(3,2) = T(2,3)
                  do jdir = 1,3
                     do idir = 1,3
                        hess(idir,iat,jdir,iat) =
     &                  hess(idir,iat,jdir,iat) + T(idir,jdir)*duefac
                     enddo
                  enddo
               endif
c
c     Third term: delXiat(r)*del(T)Xjat(r)*Diat,jat -> hess(iat,jat)
c
               if (jatcur .ne. 0 .and. iatcur .ne. 0) then
                  call dfill(9,0.d0,T,1)
                  icount = 0
                  call dgemm('n','n',nq*3,nnja,nnia,1d0,
     A                 delchi(1,1,inizia),nq*3,F,nnia,0d0,ff,nq*3)
                  do mu = inizja, ifinja
                     do n = 1, nq
                        fdchix = Amat(n,ii)*ff(n,1,mu-inizja+1)
                        fdchiy = Amat(n,ii)*ff(n,2,mu-inizja+1)
                        fdchiz = Amat(n,ii)*ff(n,3,mu-inizja+1)
                        T(1,1) = T(1,1) + fdchix*delchi(n,1,mu)
                        T(1,2) = T(1,2) + fdchix*delchi(n,2,mu)
                        T(1,3) = T(1,3) + fdchix*delchi(n,3,mu)
                        T(2,1) = T(2,1) + fdchiy*delchi(n,1,mu)
                        T(2,2) = T(2,2) + fdchiy*delchi(n,2,mu)
                        T(2,3) = T(2,3) + fdchiy*delchi(n,3,mu)
                        T(3,1) = T(3,1) + fdchiz*delchi(n,1,mu)
                        T(3,2) = T(3,2) + fdchiz*delchi(n,2,mu)
                        T(3,3) = T(3,3) + fdchiz*delchi(n,3,mu)
                     enddo
                  enddo
                  if (GRAD) then
                     do mu = inizja, ifinja
                        call dfill(nq*6,0.d0,scr,1)
                        do nu = inizia, ifinia
                           icount = icount+1
                           aaa = F(icount)
                           if (abs(aaa).gt.tol_rho_tmp)then
                              call daxpy(nq*6,aaa,heschi(1,1,nu),1,
     &                             scr,1)
                           endif
                        enddo
                        do n = 1,nq
                           s(n,1) = Cmat(n,1,ii)*scr(n,iixx)
     &                            + Cmat(n,2,ii)*scr(n,iixy)
     &                            + Cmat(n,3,ii)*scr(n,iixz)
                           s(n,2) = Cmat(n,1,ii)*scr(n,iixy)
     &                            + Cmat(n,2,ii)*scr(n,iiyy)
     &                            + Cmat(n,3,ii)*scr(n,iiyz)
                           s(n,3) = Cmat(n,1,ii)*scr(n,iixz)
     &                            + Cmat(n,2,ii)*scr(n,iiyz)
     &                            + Cmat(n,3,ii)*scr(n,iizz)
                        enddo
                        call dgemm('t','n',3,3,nq,1d0,
     A                       s,nq,delchi(1,1,mu),nq,1d0,T,3)
                        do n = 1, nq
                           s(n,1) = Cmat(n,1,ii)*heschi(n,iixx,mu)
     &                            + Cmat(n,2,ii)*heschi(n,iixy,mu)
     &                            + Cmat(n,3,ii)*heschi(n,iixz,mu)
                           s(n,2) = Cmat(n,1,ii)*heschi(n,iixy,mu)
     &                            + Cmat(n,2,ii)*heschi(n,iiyy,mu)
     &                            + Cmat(n,3,ii)*heschi(n,iiyz,mu)
                           s(n,3) = Cmat(n,1,ii)*heschi(n,iixz,mu)
     &                            + Cmat(n,2,ii)*heschi(n,iiyz,mu)
     &                            + Cmat(n,3,ii)*heschi(n,iizz,mu)
                        enddo
                        call dgemm('t','n',3,3,nq,1d0,
     A                       ff(1,1,mu-inizja+1),nq,s,nq,1d0,T,3)
                     enddo
                  endif
c
c     This term always comes with a factor of 2 in front
c
                  duefac=2d0
c
                  do jdir = 1,3
                     do idir = 1,3
                        hess(idir,iat,jdir,jat) =
     &                  hess(idir,iat,jdir,jat) + T(idir,jdir)*duefac
                        if (iat.ne.jat) then
                           hess(jdir,jat,idir,iat) =
     &                     hess(jdir,jat,idir,iat) + T(idir,jdir)*duefac
                        endif
                     enddo
                  enddo
               endif
  210       continue
  220    continue
  230 continue
      if(oprint.and.ga_nodeid().eq.0) then
         write(luout,*) ' xc_d2expl: hess '
         call output(hess,1,3*natoms,1,3*natoms,3*natoms,3*natoms,1)
      endif
c
      return
      end
C> @}
