c
c     The following macros are conveniently used to turn various
c     2nd derivative terms on and off to help with debugging.  When
c     the code is stable they can be removed.
c
#if 1
#define TERM_0A
#define TERM_1A
#define TERM_1B
#define TERM_2B
#define TERM_3B
#else
#endif
c
c     Computes explicit nuclear 2nd derivatives of the XC energy
c
c     BGJ - 8/98
c
c     $Id: xc_d2expl.F,v 1.5 1999-05-08 13:13:14 bjohnson Exp $
c
      Subroutine xc_d2expl(tol_rho, scr, Amat, Amat2, Bmat, Cmat, Cmat2,
     &     Dmat, F, Pmat, ff, s, chi, delchi, heschi, d3chi, curatoms,
     &     ncuratoms, ipol, nq, nbf, max_at_bf, GRAD0, basis, natoms,
     &     iniz, ifin, drho, ddelrho, delrho, g_dens, hess, ibf,
     &     rchi_atom, rdelchi_atom, rdens_atom, cetobfr)
c
#include "dft2drv.fh"
c
      implicit none
c
c!!! BGJ test
      logical GRAD, grad0
c!!! BGJ test
      integer basis
      integer max_at_bf ! [input]
      integer ipol  ! [input]
      integer nq    ! [input]
      integer nbf    ! [input]
      integer natoms ! [input]
      integer ncuratoms ! Number of current "active" atoms [input]
      integer curatoms(natoms) ! Mapping array for current atoms [input]
      double precision tol_rho
c
c     Explicit first derivatives of density wrt current nuclei [input]
c
      double precision drho(nq,ipol,3,ncuratoms),
     &                 ddelrho(nq,3,ipol,3,ncuratoms) ! destroyed
c
c     Spin density gradients
c
      double precision delrho(nq,3,ipol) ! [input]
c
      integer g_dens(ipol) ! GA handle for density matrices [input]
c
c     Hessian matrix (updated)
c
      double precision hess(3,natoms,3,natoms) ! [input/output]
c
      double precision scr(nq,5,3) ! Scratch matrix [input]
      double precision rchi_atom(natoms)
      double precision rdelchi_atom(natoms)
      double precision rdens_atom(natoms,natoms,ipol)
      integer cetobfr(2,natoms)
c
      double precision Pmat(max_at_bf*max_at_bf) ! scratch vector
      double precision F(max_at_bf*max_at_bf)    ! scratch vector
      double precision ff(nq,3), s(nq)   ! scratch arrays
c
c     Sampling Matrices for the XC Functional Derivatives [destroyed]
c
      double precision Amat(nq,ipol), Cmat(nq,3,ipol),
     &                 Amat2(nq,NCOL_AMAT2), Cmat2(nq,NCOL_CMAT2)
c
c     Sampling Matrices for [Products of] Basis Functions & Gradients
c
      double precision Bmat(nq,max_at_bf)
      double precision Dmat(nq,3,max_at_bf)
      integer iniz(natoms), ifin(natoms)
c
c     Basis Functions & Derivatives
c
      double precision chi(nq,nbf), delchi(nq,3,nbf), heschi(nq,6,nbf),
     &                 d3chi(nq,10,nbf)
c
      integer ibf(nbf)
c
c     local declarations
c
      double precision A_MAX, C_MAX, AC_MAX, FUNC_MAXI,
     &                 B_MAX, D_MAX, BD_MAX, FUNC_MAXJ
      double precision FUNC_MAX, DELFUNC_MAX, tol_rho_tmp
      integer iatcur, jatcur
      integer iat, inizia, ifinia, nbfia, nnia, ifirst, ilast, idim
      integer jat, inizja, ifinja, nbfja, nnja, jfirst, jlast, jdim
      integer ii, mu, nu, icount
      integer n
      double precision a11, a12, a22, c13, c14, c15, c23, c24, c25,
     &                 c33, c34, c35, c44, c45, c55, dgaa, dgab, dgbb,
     &                 aaa, fdchix, fdchiy, fdchiz,
     &                 ccc1, ccc2, ccc3
      double precision T(3,3)
      integer idir, jdir
c
c     The following parameter definitions must be consistent with
c     those in routine new_eval_gbsets_d3, or this routine will
c     not work
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
c!!! BGJ test
#if 1
      GRAD = grad0
#else
      GRAD = .false.
#endif
c!!! BGJ test
c     For restricted, combine correction factors for total densities
c     with the functional 2nd derivatives
c
      if (ipol.eq.1) then
         do n = 1, nq
            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) * 0.25d0
            Amat2(n,D2_RA_RB) = Amat2(n,D2_RA_RB) * 0.25d0
         enddo
         if (GRAD) then
            do n = 1, nq
               Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA) * 0.125d0
               Cmat2(n,D2_RA_GAB) = Cmat2(n,D2_RA_GAB) * 0.125d0
               Cmat2(n,D2_RA_GBB) = Cmat2(n,D2_RA_GBB) * 0.125d0
               Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA) * 0.0625d0
               Cmat2(n,D2_GAA_GAB) = Cmat2(n,D2_GAA_GAB) * 0.0625d0
               Cmat2(n,D2_GAA_GBB) = Cmat2(n,D2_GAA_GBB) * 0.0625d0
               Cmat2(n,D2_GAB_GAB) = Cmat2(n,D2_GAB_GAB) * 0.0625d0
            enddo
         endif
      endif
c
c     First step: separable part of gamma 2nd derivative contribution,
c                 done here before ddelrho is overwritten and Cmat is
c                 transformed
c     !!! This could be made more efficient with some extra scratch space
#ifdef TERM_0A
      if (.not.GRAD) goto 99
c
      do 10 iat = 1, natoms
         iatcur = curatoms(iat)
c     !!! are quadratic loops necessary?
         if (iatcur.eq.0) goto 10
c
         do 20 jat = 1, iat
            jatcur = curatoms(jat)
c     !!! are quadratic loops necessary?
            if (jatcur.eq.0) goto 20
c
            call dfill(9,0.d0,T,1)
            if (ipol.eq.1) then
               do jdir = 1, 3
                  do idir = 1, 3
                     do n = 1, nq
c
c     Scale fnl derivs by 1/2 and leave out factor of 2 in dot
c     to account for total densities
c
                        c33 = Cmat(n,D1_GAA,1) + Cmat(n,D1_GAB,1)*0.5d0
                        T(idir,jdir) = T(idir,jdir)
     &     + c33*(ddelrho(n,1,1,idir,iatcur)*ddelrho(n,1,1,jdir,jatcur)
     &          + ddelrho(n,2,1,idir,iatcur)*ddelrho(n,2,1,jdir,jatcur)
     &          + ddelrho(n,3,1,idir,iatcur)*ddelrho(n,3,1,jdir,jatcur))
                     enddo
                  enddo
               enddo
            else
               do jdir = 1, 3
                  do idir = 1, 3
                     do n = 1, nq
                        c33 = Cmat(n,D1_GAA,1)*2d0
                        c44 = Cmat(n,D1_GAB,1)
                        c55 = Cmat(n,D1_GBB,1)*2d0
                        T(idir,jdir) = T(idir,jdir)
     &     + c33*(ddelrho(n,1,1,idir,iatcur)*ddelrho(n,1,1,jdir,jatcur)
     &          + ddelrho(n,2,1,idir,iatcur)*ddelrho(n,2,1,jdir,jatcur)
     &          + ddelrho(n,3,1,idir,iatcur)*ddelrho(n,3,1,jdir,jatcur))
     &     + c44*(ddelrho(n,1,1,idir,iatcur)*ddelrho(n,1,2,jdir,jatcur)
     &          + ddelrho(n,2,1,idir,iatcur)*ddelrho(n,2,2,jdir,jatcur)
     &          + ddelrho(n,3,1,idir,iatcur)*ddelrho(n,3,2,jdir,jatcur)
     &          + ddelrho(n,1,2,idir,iatcur)*ddelrho(n,1,1,jdir,jatcur)
     &          + ddelrho(n,2,2,idir,iatcur)*ddelrho(n,2,1,jdir,jatcur)
     &          + ddelrho(n,3,2,idir,iatcur)*ddelrho(n,3,1,jdir,jatcur))
     &     + c55*(ddelrho(n,1,2,idir,iatcur)*ddelrho(n,1,2,jdir,jatcur)
     &          + ddelrho(n,2,2,idir,iatcur)*ddelrho(n,2,2,jdir,jatcur)
     &          + ddelrho(n,3,2,idir,iatcur)*ddelrho(n,3,2,jdir,jatcur))
                     enddo
                  enddo
               enddo
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
 99   continue
#endif
c
c     Second step: terms involving functional second derivatives
c                  and density parameter first derivatives
c
#ifdef TERM_1A
#if 0
      write(6,*) ' xc_d2expl: AMAT2'
      call output(amat2, 1, nq, 1, NCOL_AMAT2, nq, NCOL_AMAT2, 1)
      if (GRAD) then
         write(6,*) ' xc_d2expl: CMAT2'
         call output(cmat2, 1, nq, 1, NCOL_CMAT2, nq, NCOL_CMAT2, 1)
      endif
#endif
      do 100 iat = 1, natoms
         iatcur = curatoms(iat)
c     !!! are quadratic loops necessary?
         if (iatcur.eq.0) goto 100
c
         if (GRAD) then
c
c           Transform delrho derivatives into gamma derivatives,
c           storing the result back in ddelrho
c
            if (ipol.eq.1) then
               do idir = 1, 3
                  do n = 1, nq
                    dgaa =(delrho(n,1,1)*ddelrho(n,1,1,idir,iatcur)
     &                   + delrho(n,2,1)*ddelrho(n,2,1,idir,iatcur)
     &                   + delrho(n,3,1)*ddelrho(n,3,1,idir,iatcur))*2d0
                    ddelrho(n,1,1,idir,iatcur) = dgaa
                  enddo
               enddo
            else
               do idir = 1, 3
                  do n = 1, nq
                    dgaa =(delrho(n,1,1)*ddelrho(n,1,1,idir,iatcur)
     &                   + delrho(n,2,1)*ddelrho(n,2,1,idir,iatcur)
     &                   + delrho(n,3,1)*ddelrho(n,3,1,idir,iatcur))*2d0
                    dgab = delrho(n,1,1)*ddelrho(n,1,2,idir,iatcur)
     &                   + delrho(n,2,1)*ddelrho(n,2,2,idir,iatcur)
     &                   + delrho(n,3,1)*ddelrho(n,3,2,idir,iatcur)
     &                   + delrho(n,1,2)*ddelrho(n,1,1,idir,iatcur)
     &                   + delrho(n,2,2)*ddelrho(n,2,1,idir,iatcur)
     &                   + delrho(n,3,2)*ddelrho(n,3,1,idir,iatcur)
                    dgbb =(delrho(n,1,2)*ddelrho(n,1,2,idir,iatcur)
     &                   + delrho(n,2,2)*ddelrho(n,2,2,idir,iatcur)
     &                   + delrho(n,3,2)*ddelrho(n,3,2,idir,iatcur))*2d0
                    ddelrho(n,1,1,idir,iatcur) = dgaa
                    ddelrho(n,2,1,idir,iatcur) = dgab
                    ddelrho(n,3,1,idir,iatcur) = dgbb
                  enddo
               enddo
            endif
         endif
c
c        Pre-combine derivatives of density parameters wrt iat
c        with derivatives of functional
c
         if (ipol.eq.1) then
            do n = 1, nq
               a11 = (Amat2(n,D2_RA_RA) + Amat2(n,D2_RA_RB))*2d0
               scr(n,1,1) = a11 * drho(n,1,1,iatcur)
               scr(n,1,2) = a11 * drho(n,1,2,iatcur)
               scr(n,1,3) = a11 * drho(n,1,3,iatcur)
            enddo
            if (GRAD) then
               do n = 1, nq
                  c13 = (Cmat2(n,D2_RA_GAA) + Cmat2(n,D2_RA_GAB)
     &                 + Cmat2(n,D2_RA_GBB))*2d0
                  c33 = (Cmat2(n,D2_GAA_GAA) + Cmat2(n,D2_GAA_GBB))*2d0
     &                 + Cmat2(n,D2_GAA_GAB)*4d0 + Cmat2(n,D2_GAB_GAB)
c
c     The local coefficients get some extra contributions
c
                  scr(n,1,1) = scr(n,1,1)
     &                       + c13 * ddelrho(n,1,1,1,iatcur)
                  scr(n,1,2) = scr(n,1,2)
     &                       + c13 * ddelrho(n,1,1,2,iatcur)
                  scr(n,1,3) = scr(n,1,3)
     &                       + c13 * ddelrho(n,1,1,3,iatcur)
c
c     Now for the GC coefficients
c
                  scr(n,3,1) = c13 * drho(n,1,1,iatcur)
     &                       + c33 * ddelrho(n,1,1,1,iatcur)
                  scr(n,3,2) = c13 * drho(n,1,2,iatcur)
     &                       + c33 * ddelrho(n,1,1,2,iatcur)
                  scr(n,3,3) = c13 * drho(n,1,3,iatcur)
     &                       + c33 * ddelrho(n,1,1,3,iatcur)
               enddo
            endif
         else
            do n = 1, nq
               a11 = Amat2(n,D2_RA_RA)
               a12 = Amat2(n,D2_RA_RB)
               a22 = Amat2(n,D2_RB_RB)
               scr(n,1,1) = a11 * drho(n,1,1,iatcur)
     &                    + a12 * drho(n,2,1,iatcur)
               scr(n,2,1) = a12 * drho(n,1,1,iatcur)
     &                    + a22 * drho(n,2,1,iatcur)
               scr(n,1,2) = a11 * drho(n,1,2,iatcur)
     &                    + a12 * drho(n,2,2,iatcur)
               scr(n,2,2) = a12 * drho(n,1,2,iatcur)
     &                    + a22 * drho(n,2,2,iatcur)
               scr(n,1,3) = a11 * drho(n,1,3,iatcur)
     &                    + a12 * drho(n,2,3,iatcur)
               scr(n,2,3) = a12 * drho(n,1,3,iatcur)
     &                    + a22 * drho(n,2,3,iatcur)
            enddo
            if (GRAD) then
               do n = 1, nq
                  c13 = Cmat2(n,D2_RA_GAA)
                  c14 = Cmat2(n,D2_RA_GAB)
                  c15 = Cmat2(n,D2_RA_GBB)
                  c23 = Cmat2(n,D2_RB_GAA)
                  c24 = Cmat2(n,D2_RB_GAB)
                  c25 = Cmat2(n,D2_RB_GBB)
                  c33 = Cmat2(n,D2_GAA_GAA)
                  c34 = Cmat2(n,D2_GAA_GAB)
                  c35 = Cmat2(n,D2_GAA_GBB)
                  c44 = Cmat2(n,D2_GAB_GAB)
                  c45 = Cmat2(n,D2_GAB_GBB)
                  c55 = Cmat2(n,D2_GBB_GBB)
c
c     The local coefficients get some extra contributions
c
                  scr(n,1,1) = scr(n,1,1)
     &                       + c13 * ddelrho(n,1,1,1,iatcur)
     &                       + c14 * ddelrho(n,2,1,1,iatcur)
     &                       + c15 * ddelrho(n,3,1,1,iatcur)
                  scr(n,2,1) = scr(n,2,1)
     &                       + c23 * ddelrho(n,1,1,1,iatcur)
     &                       + c24 * ddelrho(n,2,1,1,iatcur)
     &                       + c25 * ddelrho(n,3,1,1,iatcur)
                  scr(n,1,2) = scr(n,1,2)
     &                       + c13 * ddelrho(n,1,1,2,iatcur)
     &                       + c14 * ddelrho(n,2,1,2,iatcur)
     &                       + c15 * ddelrho(n,3,1,2,iatcur)
                  scr(n,2,2) = scr(n,2,2)
     &                       + c23 * ddelrho(n,1,1,2,iatcur)
     &                       + c24 * ddelrho(n,2,1,2,iatcur)
     &                       + c25 * ddelrho(n,3,1,2,iatcur)
                  scr(n,1,3) = scr(n,1,3)
     &                       + c13 * ddelrho(n,1,1,3,iatcur)
     &                       + c14 * ddelrho(n,2,1,3,iatcur)
     &                       + c15 * ddelrho(n,3,1,3,iatcur)
                  scr(n,2,3) = scr(n,2,3)
     &                       + c23 * ddelrho(n,1,1,3,iatcur)
     &                       + c24 * ddelrho(n,2,1,3,iatcur)
     &                       + c25 * ddelrho(n,3,1,3,iatcur)
c
c     Now for the GC coefficients
c
                  scr(n,3,1) = c13 * drho(n,1,1,iatcur)
     &                       + c23 * drho(n,2,1,iatcur)
     &                       + c33 * ddelrho(n,1,1,1,iatcur)
     &                       + c34 * ddelrho(n,2,1,1,iatcur)
     &                       + c35 * ddelrho(n,3,1,1,iatcur)
                  scr(n,4,1) = c14 * drho(n,1,1,iatcur)
     &                       + c24 * drho(n,2,1,iatcur)
     &                       + c34 * ddelrho(n,1,1,1,iatcur)
     &                       + c44 * ddelrho(n,2,1,1,iatcur)
     &                       + c45 * ddelrho(n,3,1,1,iatcur)
                  scr(n,5,1) = c15 * drho(n,1,1,iatcur)
     &                       + c25 * drho(n,2,1,iatcur)
     &                       + c35 * ddelrho(n,1,1,1,iatcur)
     &                       + c45 * ddelrho(n,2,1,1,iatcur)
     &                       + c55 * ddelrho(n,3,1,1,iatcur)
                  scr(n,3,2) = c13 * drho(n,1,2,iatcur)
     &                       + c23 * drho(n,2,2,iatcur)
     &                       + c33 * ddelrho(n,1,1,2,iatcur)
     &                       + c34 * ddelrho(n,2,1,2,iatcur)
     &                       + c35 * ddelrho(n,3,1,2,iatcur)
                  scr(n,4,2) = c14 * drho(n,1,2,iatcur)
     &                       + c24 * drho(n,2,2,iatcur)
     &                       + c34 * ddelrho(n,1,1,2,iatcur)
     &                       + c44 * ddelrho(n,2,1,2,iatcur)
     &                       + c45 * ddelrho(n,3,1,2,iatcur)
                  scr(n,5,2) = c15 * drho(n,1,2,iatcur)
     &                       + c25 * drho(n,2,2,iatcur)
     &                       + c35 * ddelrho(n,1,1,2,iatcur)
     &                       + c45 * ddelrho(n,2,1,2,iatcur)
     &                       + c55 * ddelrho(n,3,1,2,iatcur)
                  scr(n,3,3) = c13 * drho(n,1,3,iatcur)
     &                       + c23 * drho(n,2,3,iatcur)
     &                       + c33 * ddelrho(n,1,1,3,iatcur)
     &                       + c34 * ddelrho(n,2,1,3,iatcur)
     &                       + c35 * ddelrho(n,3,1,3,iatcur)
                  scr(n,4,3) = c14 * drho(n,1,3,iatcur)
     &                       + c24 * drho(n,2,3,iatcur)
     &                       + c34 * ddelrho(n,1,1,3,iatcur)
     &                       + c44 * ddelrho(n,2,1,3,iatcur)
     &                       + c45 * ddelrho(n,3,1,3,iatcur)
                  scr(n,5,3) = c15 * drho(n,1,3,iatcur)
     &                       + c25 * drho(n,2,3,iatcur)
     &                       + c35 * ddelrho(n,1,1,3,iatcur)
     &                       + c45 * ddelrho(n,2,1,3,iatcur)
     &                       + c55 * ddelrho(n,3,1,3,iatcur)
               enddo
            endif
         endif
c
         do 200 jat = 1, iat
            jatcur = curatoms(jat)
c     !!! are quadratic loops necessary?
            if (jatcur.eq.0) goto 200
c
            call dfill(9,0.d0,T,1)
            if (ipol.eq.1) then
               do n = 1, nq
                  T(1,1) = T(1,1) + scr(n,1,1) * drho(n,1,1,jatcur)
                  T(2,1) = T(2,1) + scr(n,1,2) * drho(n,1,1,jatcur)
                  T(3,1) = T(3,1) + scr(n,1,3) * drho(n,1,1,jatcur)
                  T(1,2) = T(1,2) + scr(n,1,1) * drho(n,1,2,jatcur)
                  T(2,2) = T(2,2) + scr(n,1,2) * drho(n,1,2,jatcur)
                  T(3,2) = T(3,2) + scr(n,1,3) * drho(n,1,2,jatcur)
                  T(1,3) = T(1,3) + scr(n,1,1) * drho(n,1,3,jatcur)
                  T(2,3) = T(2,3) + scr(n,1,2) * drho(n,1,3,jatcur)
                  T(3,3) = T(3,3) + scr(n,1,3) * drho(n,1,3,jatcur)
               enddo
               if (GRAD) then
                  do n = 1, nq
                     T(1,1) = T(1,1)
     &                      + scr(n,3,1) * ddelrho(n,1,1,1,jatcur)
                     T(2,1) = T(2,1)
     &                      + scr(n,3,2) * ddelrho(n,1,1,1,jatcur)
                     T(3,1) = T(3,1)
     &                      + scr(n,3,3) * ddelrho(n,1,1,1,jatcur)
                     T(1,2) = T(1,2)
     &                      + scr(n,3,1) * ddelrho(n,1,1,2,jatcur)
                     T(2,2) = T(2,2)
     &                      + scr(n,3,2) * ddelrho(n,1,1,2,jatcur)
                     T(3,2) = T(3,2)
     &                      + scr(n,3,3) * ddelrho(n,1,1,2,jatcur)
                     T(1,3) = T(1,3)
     &                      + scr(n,3,1) * ddelrho(n,1,1,3,jatcur)
                     T(2,3) = T(2,3)
     &                      + scr(n,3,2) * ddelrho(n,1,1,3,jatcur)
                     T(3,3) = T(3,3)
     &                      + scr(n,3,3) * ddelrho(n,1,1,3,jatcur)
                  enddo
               endif
            else
               do n = 1, nq
                  T(1,1) = T(1,1) + scr(n,1,1) * drho(n,1,1,jatcur)
     &                            + scr(n,2,1) * drho(n,2,1,jatcur)
                  T(2,1) = T(2,1) + scr(n,1,2) * drho(n,1,1,jatcur)
     &                            + scr(n,2,2) * drho(n,2,1,jatcur)
                  T(3,1) = T(3,1) + scr(n,1,3) * drho(n,1,1,jatcur)
     &                            + scr(n,2,3) * drho(n,2,1,jatcur)
                  T(1,2) = T(1,2) + scr(n,1,1) * drho(n,1,2,jatcur)
     &                            + scr(n,2,1) * drho(n,2,2,jatcur)
                  T(2,2) = T(2,2) + scr(n,1,2) * drho(n,1,2,jatcur)
     &                            + scr(n,2,2) * drho(n,2,2,jatcur)
                  T(3,2) = T(3,2) + scr(n,1,3) * drho(n,1,2,jatcur)
     &                            + scr(n,2,3) * drho(n,2,2,jatcur)
                  T(1,3) = T(1,3) + scr(n,1,1) * drho(n,1,3,jatcur)
     &                            + scr(n,2,1) * drho(n,2,3,jatcur)
                  T(2,3) = T(2,3) + scr(n,1,2) * drho(n,1,3,jatcur)
     &                            + scr(n,2,2) * drho(n,2,3,jatcur)
                  T(3,3) = T(3,3) + scr(n,1,3) * drho(n,1,3,jatcur)
     &                            + scr(n,2,3) * drho(n,2,3,jatcur)
               enddo
               if (GRAD) then
                  do n = 1, nq
                     T(1,1) = T(1,1)
     &                      + scr(n,3,1) * ddelrho(n,1,1,1,jatcur)
     &                      + scr(n,4,1) * ddelrho(n,2,1,1,jatcur)
     &                      + scr(n,5,1) * ddelrho(n,3,1,1,jatcur)
                     T(2,1) = T(2,1)
     &                      + scr(n,3,2) * ddelrho(n,1,1,1,jatcur)
     &                      + scr(n,4,2) * ddelrho(n,2,1,1,jatcur)
     &                      + scr(n,5,2) * ddelrho(n,3,1,1,jatcur)
                     T(3,1) = T(3,1)
     &                      + scr(n,3,3) * ddelrho(n,1,1,1,jatcur)
     &                      + scr(n,4,3) * ddelrho(n,2,1,1,jatcur)
     &                      + scr(n,5,3) * ddelrho(n,3,1,1,jatcur)
                     T(1,2) = T(1,2)
     &                      + scr(n,3,1) * ddelrho(n,1,1,2,jatcur)
     &                      + scr(n,4,1) * ddelrho(n,2,1,2,jatcur)
     &                      + scr(n,5,1) * ddelrho(n,3,1,2,jatcur)
                     T(2,2) = T(2,2)
     &                      + scr(n,3,2) * ddelrho(n,1,1,2,jatcur)
     &                      + scr(n,4,2) * ddelrho(n,2,1,2,jatcur)
     &                      + scr(n,5,2) * ddelrho(n,3,1,2,jatcur)
                     T(3,2) = T(3,2)
     &                      + scr(n,3,3) * ddelrho(n,1,1,2,jatcur)
     &                      + scr(n,4,3) * ddelrho(n,2,1,2,jatcur)
     &                      + scr(n,5,3) * ddelrho(n,3,1,2,jatcur)
                     T(1,3) = T(1,3)
     &                      + scr(n,3,1) * ddelrho(n,1,1,3,jatcur)
     &                      + scr(n,4,1) * ddelrho(n,2,1,3,jatcur)
     &                      + scr(n,5,1) * ddelrho(n,3,1,3,jatcur)
                     T(2,3) = T(2,3)
     &                      + scr(n,3,2) * ddelrho(n,1,1,3,jatcur)
     &                      + scr(n,4,2) * ddelrho(n,2,1,3,jatcur)
     &                      + scr(n,5,2) * ddelrho(n,3,1,3,jatcur)
                     T(3,3) = T(3,3)
     &                      + scr(n,3,3) * ddelrho(n,1,1,3,jatcur)
     &                      + scr(n,4,3) * ddelrho(n,2,1,3,jatcur)
     &                      + scr(n,5,3) * ddelrho(n,3,1,3,jatcur)
                  enddo
               endif
            endif
c
c     For restricted we must correct for the densities being *total*
c     densities, calling for a factor of 1/2 * 1/2 - a factor of 2 was
c     left out in combining the functional derivatives above, leaving
c     a factor of 1/2 remaining.
c
c!!!            if (ipol.eq.1) call dscal(9,0.5d0,T,1)
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
 200     continue
 100  continue
#endif
c
c     Third step: remaining terms involving functional first derivatives
c                 and density parameter second derivatives
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
               call gat_mat(F, Pmat, nbfia, nbfja, nnia, nnja, ifirst, 
     &                      jfirst, ibf(inizia), ibf(inizja))
c
c     Three terms to compute
c
c     First term:    Xiat(r)*hessXjat(r)*Diat,jat -> hess(jat,jat)
c             GC: delXiat(r)*hessXjat(r)*Diat,jat -> hess(jat,jat)
c             GC:    Xiat(r)*  d3Xjat(r)*Diat,jat -> hess(jat,jat)
c
#ifdef TERM_1B
               if (jatcur .ne. 0) then
                  call dfill(9,0.d0,T,1)
                  icount = 0
                  do mu = inizja, ifinja
                     call dfill(nq,0.d0,s,1)
                     if (GRAD) then
                        call dfill(nq*3,0.d0,ff,1)
                     endif
                     do nu = inizia, ifinia
                        icount = icount+1
                        aaa = F(icount)
                        if (abs(aaa).gt.tol_rho_tmp)then
                           call daxpy(nq,aaa,chi(1,nu),1,s,1)
                           if (GRAD) then
                              call daxpy(nq*3,aaa,delchi(1,1,nu),1,ff,1)
                           endif
                        endif
                     enddo
c
c     Prepare coefficient of hessXjat
c
                     if (GRAD) then
                        do n = 1, nq
                           ff(n,1) = Amat(n,ii)*s(n)
     &                             + Cmat(n,1,ii)*ff(n,1)
     &                             + Cmat(n,2,ii)*ff(n,2)
     &                             + Cmat(n,3,ii)*ff(n,3)
                        enddo
                     else
                        do n = 1, nq
                           ff(n,1) = Amat(n,ii)*s(n)
                        enddo
                     endif
c
                     do n = 1, nq
                        aaa = ff(n,1)
                        T(1,1) = T(1,1) + heschi(n,1,mu)*aaa
                        T(1,2) = T(1,2) + heschi(n,2,mu)*aaa
                        T(1,3) = T(1,3) + heschi(n,3,mu)*aaa
                        T(2,2) = T(2,2) + heschi(n,4,mu)*aaa
                        T(2,3) = T(2,3) + heschi(n,5,mu)*aaa
                        T(3,3) = T(3,3) + heschi(n,6,mu)*aaa
                     enddo
                     if (GRAD) then
                        do n = 1, nq
                           ccc1 = Cmat(n,1,ii)*s(n)
                           ccc2 = Cmat(n,2,ii)*s(n)
                           ccc3 = Cmat(n,3,ii)*s(n)
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
                     endif
                  enddo
c
                  if (iat.ne.jat) call dscal(9,2.d0,T,1)
c
                  T(2,1) = T(1,2)
                  T(3,1) = T(1,3)
                  T(3,2) = T(2,3)
                  do jdir = 1,3
                     do idir = 1,3
                        hess(idir,jat,jdir,jat) =
     &                  hess(idir,jat,jdir,jat) + T(idir,jdir)
                     enddo
                  enddo
               endif
#endif
c
c     Second term: hessXiat(r)*   Xjat(r)*Diat,jat -> hess(iat,iat)
c              GC: hessXiat(r)*delXjat(r)*Diat,jat -> hess(iat,iat)
c              GC:   d3Xiat(r)*   Xjat(r)*Diat,jat -> hess(iat,iat)
c
#ifdef TERM_2B
               if (iatcur .ne. 0) then
                  call dfill(9,0.d0,T,1)
                  do nu = inizia, ifinia
                     call dfill(nq,0.d0,s,1)
                     if (GRAD) then
                        call dfill(nq*3,0.d0,ff,1)
                     endif
                     icount = nu-inizia+1
                     do mu = inizja, ifinja
                        aaa = F(icount)
                        icount = icount+nnia
                        if (abs(aaa).gt.tol_rho_tmp)then
                           call daxpy(nq,aaa,chi(1,mu),1,s,1)
                           if (GRAD) then
                              call daxpy(nq*3,aaa,delchi(1,1,mu),1,ff,1)
                           endif
                        endif
                     enddo
c
c     Prepare coefficient of hessXiat
c
                     if (GRAD) then
                        do n = 1, nq
                           ff(n,1) = Amat(n,ii)*s(n)
     &                             + Cmat(n,1,ii)*ff(n,1)
     &                             + Cmat(n,2,ii)*ff(n,2)
     &                             + Cmat(n,3,ii)*ff(n,3)
                        enddo
                     else
                        do n = 1, nq
                           ff(n,1) = Amat(n,ii)*s(n)
                        enddo
                     endif
c
                     do n = 1, nq
                        aaa = ff(n,1)
                        T(1,1) = T(1,1) + heschi(n,1,nu)*aaa
                        T(1,2) = T(1,2) + heschi(n,2,nu)*aaa
                        T(1,3) = T(1,3) + heschi(n,3,nu)*aaa
                        T(2,2) = T(2,2) + heschi(n,4,nu)*aaa
                        T(2,3) = T(2,3) + heschi(n,5,nu)*aaa
                        T(3,3) = T(3,3) + heschi(n,6,nu)*aaa
                     enddo
                     if (GRAD) then
                        do n = 1, nq
                           ccc1 = Cmat(n,1,ii)*s(n)
                           ccc2 = Cmat(n,2,ii)*s(n)
                           ccc3 = Cmat(n,3,ii)*s(n)
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
                     endif
                  enddo
c
                  if (iat.ne.jat) call dscal(9,2.d0,T,1)
c
                  T(2,1) = T(1,2)
                  T(3,1) = T(1,3)
                  T(3,2) = T(2,3)
                  do jdir = 1,3
                     do idir = 1,3
                        hess(idir,iat,jdir,iat) =
     &                  hess(idir,iat,jdir,iat) + T(idir,jdir)
                     enddo
                  enddo
               endif
#endif
c
c     Third term: delXiat(r)*del(T)Xjat(r)*Diat,jat -> hess(iat,jat)
c
#ifdef TERM_3B
               if (jatcur .ne. 0 .and. iatcur .ne. 0) then
                  call dfill(9,0.d0,T,1)
                  icount = 0
                  do mu = inizja, ifinja
                     call dfill(nq*3,0.d0,ff,1)
                     do nu = inizia, ifinia
                        icount = icount+1
                        aaa = F(icount)
                        if (abs(aaa).gt.tol_rho_tmp)then
                           call daxpy(nq*3,aaa,delchi(1,1,nu),1,ff,1)
                        endif
                     enddo
                     do n = 1, nq
                        fdchix = Amat(n,ii)*ff(n,1)
                        fdchiy = Amat(n,ii)*ff(n,2)
                        fdchiz = Amat(n,ii)*ff(n,3)
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
c
c     This term always comes with a factor of 2 in front
c
                  call dscal(9,2.d0,T,1)
c
                  do jdir = 1,3
                     do idir = 1,3
                        hess(idir,iat,jdir,jat) =
     &                  hess(idir,iat,jdir,jat) + T(idir,jdir)
                        if (iat.ne.jat) then
                           hess(jdir,jat,idir,iat) =
     &                     hess(jdir,jat,idir,iat) + T(idir,jdir)
                        endif
                     enddo
                  enddo
               endif
#endif
  210       continue
  220    continue
  230 continue
c      write(6,*) ' xc_d2expl: hess '
c      call output(hess,1,3*natoms,1,3*natoms,3*natoms,3*natoms,1)
      return
      end
