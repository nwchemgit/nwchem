c
C$Id: xc_tabcd.F,v 1.13 2005-09-27 20:15:02 edo Exp $
c
c     Modified to handle CPKS while reusing existing code
c
c     BGJ - 8/98
c
      Subroutine xc_tabcd(what,l3d,
     ,     tol_rho, Tmat, TTmat, Amat, Bmat, Cmat, Dmat, 
     &                    chi, delchi, heschi, 
     N     curatoms,ncuratoms,nmat,
     I     ipol, nq, nbf, max_at_bf, max_at_bf2, 
     G     GRAD, basis, natoms, iniz, ifin, 
     &                    g_vxc, ibf, rchi_atom, rdelchi_atom, 
     &                    rdens_atom, cetobfr,kske,Mmat)
c
c     We're using xc_tabcd for CPKS purposes
c
      implicit none
#include "dft2drv.fh"
c
      Logical GRAD
      integer what ! [in]
cwhat=0
cwhat=1 CPKS_LHS
cwhat=2 CPKS_RHS
      integer basis
      integer max_at_bf ! [input]
      integer max_at_bf2 ! [input]
      integer nmat ! Number of XC matrices (alpha + beta sets) to make [input]
c                  ! e.g. number of perturbations for CPSCF
      integer imat ! XC matrix loop index
      integer ipol  ! [input]
      integer nq    ! [input]
      integer nbf    ! [input]
      integer natoms ! [input]
      integer ncuratoms ! Number of current "active" atoms [input]
      integer curatoms(*) ! Mapping array for current atoms [input]
      integer jatcur,  nu, nu1,  indT
      double precision tol_rho
      integer g_vxc(*)   ! GA handle [input]
      double precision Tmat(*), TTmat(*)
      double precision rchi_atom(natoms)
      double precision rdelchi_atom(natoms)
      double precision rdens_atom(natoms,natoms,ipol)
      integer cetobfr(2,natoms)
c
c     Sampling Matrices for the XC Potential & Energy
c
      double precision Amat(nq,ipol,*), Cmat(nq,3,ipol,*)
      double precision Mmat(nq,ipol) 
      logical kske

c#elif defined(TABCD_CPKS_LHS)
c
c     Note: Meaning of dimensioning of Amat and Cmat changes for
c           second derivatives, simulating "overloading" of
c           Amat and Cmat
c
c     Sampling Matrices for the XC part of integrand when making
c     multiple matrices, e.g. XC part of perturbations for CPSCF
c
c      double precision Amat(nq,ipol,nmat), Cmat(nq,3,ipol,nmat)
c#elif defined(TABCD_CPKS_RHS)
c
c     For explicit nuclear derivatives of XC matrix, the same functional
c     derivative values are combined with different basis fn derivatives
c
c      double precision Amat(nq,ipol), Cmat(nq,3,ipol)
 
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
      double precision heschi(nq,6,*)
      integer ibf(nbf)
      double precision A_MAX, C_MAX, AC_MAX, FUNC_MAXI,
     &                 B_MAX, D_MAX, BD_MAX, FUNC_MAXJ
      integer iat, inizia, ifinia, nbfia, nnia, ifirst, ilast
      integer jat, inizja, ifinja, nbfja, nnja, jfirst, jlast
      integer ii, mu, mu1
      integer n,lastjat
      double precision chi1
      double precision dabsmax
      external dabsmax
      logical l3d
      integer jrsh,jrsh2,n3d,idir,jdir
      integer g_update(2)
#include "mafdecls.fh"
#include "dftpara.fh"
c
c     Indexing array for basis function hessian columns as if
c     it were a 3x3 matrix
c
      integer indh(3,3)
      logical w01,w02
      double precision ddot
      external ddot
      data indh / 1, 2, 3,
     &            2, 4, 5,
     &            3, 5, 6 /
c
c         0: l3d=.f.    & n3d=1
ccc	rhs: l3d=.true. & n3d=3
ccc	lhs: l3d=.true. & n3d=1
c
      if(what.eq.0) then
        n3d=1
      elseif(what.eq.1) then
         n3d=1
      elseif(what.eq.2) then
         n3d=3
      else
         call errquit(' wrong what value for xctabcd ',0,0)
      endif
      w01=what.eq.0.or.what.eq.1
      w02=what.eq.0.or.what.eq.2
c
c
c     Beginning of loop over multiple XC matrices
c
      do 500 imat = 1,nmat
c
c     Compute the matrix product for the XC potential and energy:
c
c              T = transpose(A*B) + transpose(C*D)
c
         A_MAX = dabsmax(nq*ipol,Amat(1,1,imat))
         if (GRAD) then
            C_MAX = dabsmax(nq*3*ipol,Cmat(1,1,1,imat))
         else
            C_MAX = 0d0
         endif
         AC_MAX = max(A_MAX,C_MAX)
c
c     repl stuff
c
      if(xcreplicated.and.dorepxc) then
         g_update(1)=k_repxc(1)
         g_update(2)=k_repxc(2)
      else
         g_update(1)=g_vxc(1)
         g_update(2)=g_vxc(2)
      endif
      do 430 iat = 1, natoms
         inizia = iniz(iat)
         if (inizia.eq.0)goto 430
         FUNC_MAXI = rchi_atom(iat)
         if(GRAD) FUNC_MAXI = max(FUNC_MAXI,rdelchi_atom(iat))
         if(what.eq.0.or.what.eq.1) then
            if (((AC_MAX*FUNC_MAXI).lt.tol_rho)) goto 430
         else
c
c     !!! Should something be done with cutoffs for other cases? !!!
c               
         endif
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


                     Bmat(n,mu) = Amat(n,ii,imat)*chi1 +
     &                    delchi(n,1,mu1)*Cmat(n,1,ii,imat) +
     &                    delchi(n,2,mu1)*Cmat(n,2,ii,imat) +
     &                    delchi(n,3,mu1)*Cmat(n,3,ii,imat)
                     Dmat(n,1,mu) = Cmat(n,1,ii,imat)*chi1
                     Dmat(n,2,mu) = Cmat(n,2,ii,imat)*chi1
                     Dmat(n,3,mu) = Cmat(n,3,ii,imat)*chi1
                     
                     if(kske) then 
                        Dmat(n,1,mu) = Dmat(n,1,mu) +
     &                       Mmat(n,ii)*delchi(n,1,mu1)
                        Dmat(n,2,mu) = Dmat(n,2,mu) +
     &                       Mmat(n,ii)*delchi(n,2,mu1)
                        Dmat(n,3,mu) = Dmat(n,3,mu) +
     &                       Mmat(n,ii)*delchi(n,3,mu1)
                     endif
                     
                  enddo
               else
                  do n = 1, nq
                     Bmat(n,mu) = chi(n,mu1)*Amat(n,ii,imat)
                  enddo
               endif
            enddo
c
c            write(6,*) ' xc_tabcd: BMAT '
c            call output(bmat, 1, nq, 1, nnia, nq, max_at_bf, 1)
c            write(6,*) ' xc_tabcd: DMAT '
c            call output(dmat, 1, 3*nq, 1, nnia, 3*nq, max_at_bf, 1)
c
            B_MAX = dabsmax(nnia*nq,Bmat)
            if (GRAD) then
               D_MAX = dabsmax(nnia*nq*3,Dmat)
            else
               D_MAX = 0d0
            endif
            BD_MAX = max(B_MAX,D_MAX)
c
            lastjat=iat
            if(what.eq.2) lastjat=natoms
            do 168 jat = 1, lastjat
               if(what.eq.2) then
c
c     To fit better into existing structure, loop over full square
c     of atom pairs and only compute nuclear derivative contribution
c     from jat.  Also, this way we only need check jatcur once and
c     for all, and don't have to check iatcur at all.
c
               jatcur = curatoms(jat)
               if (jatcur.eq.0) goto 168
            endif
               inizja = iniz(jat)
               if (inizja.eq.0) goto 168
               FUNC_MAXJ = rchi_atom(jat)
               if(grad) FUNC_MAXJ = max(rchi_atom(jat),FUNC_MAXJ)
               if(w01) then
                  if (BD_MAX*FUNC_MAXJ.lt.tol_rho) goto 168
               else
c
c     !!! Should something be done with cutoffs for other cases? !!!
c               
               endif
               ifinja = ifin(jat)
               jfirst = cetobfr(1,jat)
               jlast = cetobfr(2,jat)
               nbfja = jlast - jfirst + 1
               nnja = ifinja - inizja + 1
               if(what.eq.2) then
c
c              Loop over x, y, z directions for derivative XC mats
c
                  call dfill(n3d*max_at_bf2, 0.d0, TTmat,1) 
               endif
               do 150 jdir = 1,n3d
c
                  if(w01) then
                     call dgemm('t', 'n', nnia, nnja, nq, 1.d0, Bmat,
     &                    nq, chi(1,inizja), nq, 0.d0, Tmat, nnia)
                  else
c     Note the sign change for a nuclear derivative, and also that the
c     leading dimension of delchi must be set correctly
                     call dgemm('t', 'n', nnia, nnja, nq, -1.d0, Bmat,
     &                    nq, delchi(1,jdir,inizja), nq*3, 0.d0, Tmat,
     &                    nnia)
                  endif
                  if (GRAD) then
                     if(w01) then
                        call dgemm('t', 'n', nnia, nnja, 3*nq,
     &                       1.d0, Dmat, 3*nq, delchi(1,1,inizja), 
     &                       3*nq, 1.d0, Tmat, nnia)
                     else
                        indT = 0
                        do nu = 1, nnja
                           nu1 = nu+inizja-1
                           do mu = 1, nnia
                              indT = indT + 1
                              do idir = 1, 3
                                 Tmat(indT) = Tmat(indT)-
     *                                ddot(nq,Dmat(1,idir,mu),1,
     &                                heschi(1,indh(idir,jdir),nu1),1)
                              enddo
                           enddo
                        enddo
                     endif
                  endif
                  if(n3d.eq.1) then
                     call dfill(max_at_bf2, 0.d0, TTmat, 1)
                     call scat_mat(TTmat, Tmat, nbfia, nbfja, nnia,
     &                    nnja,ifirst,jfirst,ibf(inizia),ibf(inizja))
                  else
                     call scat_mat3(n3d,jdir,
     &                    TTmat, Tmat, nbfia, nbfja, nnia,
     &                    nnja,ifirst,jfirst,ibf(inizia),ibf(inizja))
                  endif
                     
  150          continue
               if(what.eq.0) then
                  jrsh=ii	
               elseif(what.eq.1) then
                  jrsh=imat+(ii-1)*nmat
               else
                  if (dabsmax(nbfia*nbfja*3,ttmat).lt.tol_rho)goto 168
                  jrsh=1+(jat-1)*3+(ii-1)*3*natoms
               endif
                  jrsh2=jrsh+n3d-1
                  if(l3d) then
                     call dft_3dacc(g_vxc, ttmat,
     &                    jrsh,jrsh2,
     %                 ifirst, ilast, jfirst, jlast, nbfia)
                  else
                     call upd_atom_block(g_update(ii),
     .                    basis,iat,jat,TTmat)
                  endif
c
c                 check to see if can skip and use ga_symmetrize
c
                  if ((w01.and.iat.ne.jat).or.what.eq.2) then
c     For CPKS RHS, we update with transpose even for iat = jat,
c     since that is necessary to get both contributions
c     mu * del(nu) and del(mu) * nu
                     
                     if(n3d.eq.1) then
                        call transp_mat(TTmat, Tmat,
     ,                       nbfia, nbfja)
                     else
                        call transp_mat3(n3d,TTmat, Tmat,
     ,                       nbfia, nbfja)
                     endif
                     
                     if(l3d) then
                     call dft_3dacc(g_vxc, tmat,
     &                       jrsh,jrsh2,
     %                       jfirst, jlast, ifirst, ilast, nbfja)
                  else
                     call upd_atom_block(g_update(ii), basis, jat,iat,
     &                    tmat)
                  endif
               endif
  168       continue
  405    continue
  430 continue
  500 continue
c      write(6,*)' g_vxc(1) from xc_tabcd '
c      call ga_print(g_vxc(1))
c      write(6,*)' g_vxc(2) from xc_tabcd '
c      call ga_print(g_vxc(2))
      return
      end
