c
c     Computes explicit nuclear 2nd derivatives of the XC energy
c
c     BGJ - 8/97
c
c     $Id: xc_d2expl.F,v 1.1 1998-09-09 13:21:13 bjohnson Exp $
c
      Subroutine xc_d2expl(tol_rho, scr, Amat, Amat2, Bmat, Cmat, Cmat2,
     &     Dmat, chi, delchi, heschi, curatoms, ncuratoms, ipol, nq,
     &     nbf, max_at_bf, GRAD, basis, natoms, iniz, ifin, 
     &     rho, delrho, g_dens, hess, ibf, rchi_atom, rdelchi_atom,
     &     rdens_atom, cetobfr)
c
#include "dft2drv.fh"
c
      implicit none
c
      Logical GRAD
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
      double precision rho(nq,ipol,3,ncuratoms),
     &                 delrho(nq,3,ipol,3,ncuratoms)
c
      integer g_dens(ipol) ! GA handle for density matrices [input]
c
c     Hessian matrix (updated)
c
      double precision hess(3,natoms,3,natoms) ! [input/output]
c
      double precision scr(2,3,nq,ipol) ! Scratch matrix [input]
      double precision rchi_atom(natoms)
      double precision rdelchi_atom(natoms)
      double precision rdens_atom(natoms,natoms,ipol)
      integer cetobfr(2,natoms)
c
c     Sampling Matrices for the XC Functional Derivatives
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
c     Basis Functions & Gradients
c
      double precision chi(nq,nbf), delchi(nq,3,nbf), heschi(nq,6,nbf)
c
      integer ibf(nbf)
c
c     local declarations
c
      double precision A_MAX, C_MAX, AC_MAX, FUNC_MAXI,
     &                 B_MAX, D_MAX, BD_MAX, FUNC_MAXJ
      integer iatcur, jatcur
      integer iat, inizia, ifinia, nbfia, nnia, ifirst, ilast
      integer jat, inizja, ifinja, nbfja, nnja, jfirst, jlast
      integer ii, mu, mu1
      integer ind2(2,2)
      integer n
      double precision a21, a22
      double precision T(3,3)
      integer idir, jdir
      double precision chi1
      double precision dabsmax
      external dabsmax
c
c     d2Exc                /      d2fxc      drho(i) drho(j)
c     ----- =  sum   sum   | --------------- ------- -------
c     dx dy   i=a,b j=a,b /  drho(i) drho(j)   dx      dy
c
c                          /  dfxc   d2rho(i)
c                 +  sum   | ------- --------
c                   i=a,b /  drho(i)  dx dy
c
      if (GRAD) then
         call errquit('xc_d2expl: 2nds not yet implemented for GC',0)
      endif
c
c     First step: terms involving functional second derivatives
c                 and density first derivatives
c
c     Indices for 2nd derivs of functionals
c
      ind2(1,1) = D2_RA_RA
      ind2(2,1) = D2_RA_RB
      ind2(1,2) = D2_RA_RB
      ind2(2,2) = D2_RB_RB
c
      do 100 iat = 1, natoms
         iatcur = curatoms(iat)
c     !!! are quadratic loops necessary?
         if (iatcur.eq.0) goto 100
c
c        Pre-combine derivatives of density wrt iat with derivatives
c        of functional
c
         do ii = 1, ipol
            do n = 1, nq
               a21 = Amat2(n,ind2(1,ii))
               a22 = Amat2(n,ind2(2,ii))
               scr(1,1,n,ii) = a21 * rho(n,ii,1,iatcur)
               scr(2,1,n,ii) = a22 * rho(n,ii,1,iatcur)
               scr(1,2,n,ii) = a21 * rho(n,ii,2,iatcur)
               scr(2,2,n,ii) = a22 * rho(n,ii,2,iatcur)
               scr(1,3,n,ii) = a21 * rho(n,ii,3,iatcur)
               scr(2,3,n,ii) = a22 * rho(n,ii,3,iatcur)
            enddo
         enddo
c
         do 200 jat = 1, iat
            jatcur = curatoms(jat)
c     !!! are quadratic loops necessary?
            if (jatcur.eq.0) goto 200
c
            call dfill(9,0.d0,T,1)
            do ii = 1, ipol
c
c     !!! Should work for restricted, but could be made more efficient
c     !!! by saving a multiply when ipol = 1 - probably not too important
c
               do n = 1, nq
                  T(1,1) = T(1,1) + scr(1,1,n,ii) * rho(n,   1,1,jatcur)
     &                            + scr(2,1,n,ii) * rho(n,ipol,1,jatcur)
                  T(2,1) = T(2,1) + scr(1,2,n,ii) * rho(n,   1,1,jatcur)
     &                            + scr(2,2,n,ii) * rho(n,ipol,1,jatcur)
                  T(3,1) = T(3,1) + scr(1,3,n,ii) * rho(n,   1,1,jatcur)
     &                            + scr(2,3,n,ii) * rho(n,ipol,1,jatcur)
                  T(1,2) = T(1,2) + scr(1,1,n,ii) * rho(n,   1,2,jatcur)
     &                            + scr(2,1,n,ii) * rho(n,ipol,2,jatcur)
                  T(2,2) = T(2,2) + scr(1,2,n,ii) * rho(n,   1,2,jatcur)
     &                            + scr(2,2,n,ii) * rho(n,ipol,2,jatcur)
                  T(3,2) = T(3,2) + scr(1,3,n,ii) * rho(n,   1,2,jatcur)
     &                            + scr(2,3,n,ii) * rho(n,ipol,2,jatcur)
                  T(1,3) = T(1,3) + scr(1,1,n,ii) * rho(n,   1,3,jatcur)
     &                            + scr(2,1,n,ii) * rho(n,ipol,3,jatcur)
                  T(2,3) = T(2,3) + scr(1,2,n,ii) * rho(n,   1,3,jatcur)
     &                            + scr(2,2,n,ii) * rho(n,ipol,3,jatcur)
                  T(3,3) = T(3,3) + scr(1,3,n,ii) * rho(n,   1,3,jatcur)
     &                            + scr(2,3,n,ii) * rho(n,ipol,3,jatcur)
               enddo
            enddo
c
c     Restricted calculations need a factor of 2 to account for the
c     beta contribution, but since for restricted the density evaluator
c     computes *total* densities, a factor of 1/2 * 1/2 is needed to
c     correct for this, giving a net factor of 1/2.
c
            if (ipol.eq.1) call dscal(9,0.5d0,T,1)
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
#if 0
      comment this out until done !!!
c
c     Second step: terms involving functional first derivatives
c                  and density second derivatives
c
      A_MAX = dabsmax(nq*ipol,Amat)
      C_MAX = dabsmax(nq*3*ipol,Cmat)
      AC_MAX = max(A_MAX,C_MAX)
c
c      write(6,*) ' xc_d2expl: AMAT '
c      call output(amat, 1, nq, 1, ipol, nq, ipol, 1)
c      write(6,*) ' xc_d2expl: CMAT '
c      call output(cmat, 1, 3*nq, 1, ipol, 3*nq, ipol, 1)
c      write(6,*) ' xc_d2expl: chi '
c      call output(chi, 1, nq, 1, nbf, nq, nbf, 1)
c      write(6,*) ' xc_d2expl: delchi '
c      call output(delchi, 1, 3*nq, 1, nbf, 3*nq, nbf, 1)
c
      do 430 iat = 1, natoms
         if (curatoms(iat).eq.0) goto 430
c !!! Is this redundant with the above line?
         inizia = iniz(iat)
         if (inizia.eq.0)goto 430
c
         FUNC_MAXI = max(rchi_atom(iat),rdelchi_atom(iat))
#if 0
         if (((AC_MAX*FUNC_MAXI).lt.tol_rho)) goto 430
c     !!! Should something be done with cutoffs here? !!!
#endif
         ifinia = ifin(iat)
         ifirst = cetobfr(1,iat)
         ilast = cetobfr(2,iat)
         nnia = ifinia - inizia + 1
         nbfia = ilast - ifirst + 1
c !!! Can move up bottom of ipol loop
         do 405 ii = 1, ipol
            do mu = 1, nnia
               mu1 = mu+inizia-1
               if (GRAD) then
                  call errquit('xc_d2expl: GC not ready yet',0)
                  do n = 1, nq
#if 0
                     chi1 = chi(n,mu1)
                     Bmat(n,mu) = Amat(n,ii)*chi1 +
     &                            delchi(n,1,mu1)*Cmat(n,1,ii) +
     &                            delchi(n,2,mu1)*Cmat(n,2,ii) +
     &                            delchi(n,3,mu1)*Cmat(n,3,ii)
                     Dmat(n,1,mu) = Cmat(n,1,ii)*chi1
                     Dmat(n,2,mu) = Cmat(n,2,ii)*chi1
                     Dmat(n,3,mu) = Cmat(n,3,ii)*chi1
#endif
                  enddo
               else
                  do n = 1, nq
                     ???
                     asdfdasfsdasd
                     Bmat(n,mu) = chi(n,mu1)*Amat(n,ii)
                  enddo
               endif
            enddo
c
c            write(6,*) ' xc_d2expl: BMAT '
c            call output(bmat, 1, nq, 1, nnia, nq, max_at_bf, 1)
c            write(6,*) ' xc_d2expl: DMAT '
c            call output(dmat, 1, 3*nq, 1, nnia, 3*nq, max_at_bf, 1)
c
#if 0
                                !!! cutoffs !!!
            B_MAX = dabsmax(nnia*nq,Bmat)
            D_MAX = dabsmax(nnia*nq*3,Dmat)
            BD_MAX = max(B_MAX,D_MAX)
#else
            B_MAX = 0D0
            D_MAX = 0D0
            BD_MAX = 0D0
#endif
c
            do 168 jat = 1, natoms
               if (curatoms(jat).eq.0) goto 168
c !!! Is this redundant with the followint lines?
               inizja = iniz(jat)
               if (inizja.eq.0)goto 168
               FUNC_MAXJ = max(rchi_atom(jat),rdelchi_atom(jat))
#if 0
               if (BD_MAX*FUNC_MAXJ.lt.tol_rho) goto 168
c     !!! Should do something with cutoffs here? !!!
#endif
               ifinja = ifin(jat)
               jfirst = cetobfr(1,jat)
               jlast = cetobfr(2,jat)
               nbfja = jlast - jfirst + 1
               nnja = ifinja - inizja + 1
c
c !!! What to do here?
               if (GRAD) call errquit('xc_d2expl: GC not ready',0)
               ???

c--------------------------------------------------------------------------


c
c     Screening is accomplished by:  p(r) <= |Xi(r)|*|Xj(r)|*|Dij|
c     Xi(r) is screened on desired accuracy/max(|Xj(r)|)*max(|Dij|)
c     Dij is screened on desired accuracy/max(|Xi(r)|)*max(|Xj(r)|)
c
      FUNC_MAX = dabsmax(natoms,rchi_atom)
      DELFUNC_MAX = dabsmax(natoms,rdelchi_atom)
c
      do 230 iat = 1, natoms
         if (curatoms(iat).eq.0) goto 230
c !!! Reconcile with following lines
         inizia = iniz(iat)
         if (inizia.eq.0)goto 230
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
            if (curatoms(jat).eq.0) goto 220
c !!! Reconcile with following lines
            inizja = iniz(jat)
            if (inizja.eq.0)goto 220
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
c              Compute Xiat(r)*Xjat(r)*Diat,jat
c
               icount = 0
               do mu = inizja, ifinja
                  call dfill(nq,0.d0,ff,1)
                  do nu = inizia, ifinia
                     icount = icount+1
                     aaa = F(icount)
                     if (abs(aaa).gt.tol_rho_tmp)then
                        call daxpy(nq,aaa,chi(1,nu),1,ff,1)
                     endif
                  enddo
                  if (iat.ne.jat) call dscal(nq,2.d0,ff,1)
                  do n = 1, nq
                     scr(n) = chi(n,mu)*ff(n)
                  enddo
c
                  ??? what to do here ???
c !!!                  call daxpy(nq,1.d0,scr,1,rho(1,i0+ii),1)
c
c                 Compute delXjat(r)*Xiat(r)*Diat,jat
c
                  if (grad)then
                     do n = 1, nq
                        delrho(n,1,ii) = delrho(n,1,ii) + 
     &                                   delchi(n,1,mu)*ff(n)
                        delrho(n,2,ii) = delrho(n,2,ii) + 
     &                                   delchi(n,2,mu)*ff(n)
                        delrho(n,3,ii) = delrho(n,3,ii) + 
     &                                   delchi(n,3,mu)*ff(n)
                     enddo
                  endif
               enddo
               if (grad)then
c
c                 Compute delXiat(r)*Xjat(r)*Diat,jat
c                  
                  do nu = inizia, ifinia
                     call dfill(nq,0.d0,ff,1)
                     icount = nu-inizia+1
                     do mu = inizja, ifinja
                        aaa = F(icount)
                        icount = icount+nnia
                        if (abs(aaa).gt.tol_rho_tmp)then
                           call daxpy(nq,aaa,chi(1,mu),1,ff,1)
                        endif
                     enddo
                     if (iat.ne.jat) call dscal(nq,2.d0,ff,1)
                     do n = 1, nq
                        delrho(n,1,ii) = delrho(n,1,ii) + 
     &                                   delchi(n,1,nu)*ff(n)
                        delrho(n,2,ii) = delrho(n,2,ii) + 
     &                                   delchi(n,2,nu)*ff(n)
                        delrho(n,3,ii) = delrho(n,3,ii) + 
     &                                   delchi(n,3,nu)*ff(n)
                     enddo
                  enddo
               endif
  210       continue
  220    continue
  225    continue
  230 continue


c--------------------------------------------------------------------------

c
c!!!               if (ipol.eq.1) call dscal(9,2.d0,T,1)
               if (ipol.eq.1) call dscal(9,0.5d0,T,1)
c
c              Update Hessian block(s)
c
               do jdir = 1,3
                  do idir = 1,3
                     hess(idir,iat,jdir,jat) = hess(idir,iat,jdir,jat)
     &                                       + T(idir,jdir)
                     if (iat.ne.jat) then
                        hess(jdir,jat,idir,iat) =
     &                  hess(jdir,jat,idir,iat) + T(idir,jdir)
                     endif
                  enddo
               enddo
c
c     !!! Can move pol loop up !!!
  168       continue
  405    continue
  430 continue
#endif
c      write(6,*) ' xc_d2expl: hess '
c      call output(hess,1,3*natoms,1,3*natoms,3*natoms,3*natoms,1)
      return
      end
