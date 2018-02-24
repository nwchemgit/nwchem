c
C$Id$
c
c     Modified to handle CPKS while reusing existing code
c
c     BGJ - 8/98
c
C> \ingroup nwdft_xc
C> @{
C>
C> \brief Evaluate various matrix elements over the functionals
C>
C> This routine evaluates a variety of matrix elements that involve the
C> density functional.
C>
C> Note that there are some inconsistencies with respect to the kinetic
C> energy density. NWChem defines the kinetic energy density as
C> \f{eqnarray*}{
C>    \tau_1 &=& \sum_i \frac{1}{2}\langle\phi_i|\nabla^2|\phi_i\rangle
C> \f}
C> In most density functionals the kinetic energy density is defined as
C> \f{eqnarray*}{
C>    \tau_2 &=& \sum_i \langle\phi_i|\nabla^2|\phi_i\rangle
C> \f}
C> the reason being that the factor 1/2 can easily be absorpted into the
C> density functional parameters. As a result the density functional routines
C> evaluate 
C> \f{eqnarray*}{
C>    e &=& f(\rho,\gamma,\tau_2) \\\\
C>    Amat = &=& \frac{\partial f(\rho,\gamma,\tau_2)}{\partial\rho} \\\\
C>    Cmat = &=& \frac{\partial f(\rho,\gamma,\tau_2)}{\partial\gamma} \\\\
C>    Mmat = &=& \frac{\partial f(\rho,\gamma,\tau_2)}{\partial\tau_2} \\\\
C> \f}
C> The NWDFT module in NWChem handles this situation by conceptually wrapping
C> the functional routines used (the actual implementation just scales
C> \f$\tau\f$ by a factor 2)
C> \f{eqnarray*}{
C>   e &=& g(\rho,\gamma,\tau_1) \\\\
C>   g(\rho,\gamma,\tau_1) &=& f(\rho,\gamma,2\tau_1) \\\\
C> \f}
C> As a result, however, the density functional routines now return
C> \f{eqnarray*}{
C>    e &=& g(\rho,\gamma,\tau_1) \\\\
C>    Amat = &=& \frac{\partial g(\rho,\gamma,\tau_1)}{\partial\rho} \\\\
C>    Cmat = &=& \frac{\partial g(\rho,\gamma,\tau_1)}{\partial\gamma} \\\\
C>    Mmat = &=& \frac{\partial g(\rho,\gamma,\tau_1)}{\partial\tau_2} \\\\
C> \f}
C> where there is an inconsistency in the derivative wrt. to \f$\tau\f$ as
C> \f$g\f$ is a function of \f$\tau_1\f$ but differentiated wrt. \f$\tau_2\f$.
C>
C> In the NWXC module the differentation techniques deployed do not allow 
C> such inconsistencies and \f$g\f$ is strictly differentiated wrt.
C> \f$\tau_1\f$. This results in a factor \f$2\f$ difference that we need to
C> compensate for. 
C>
      Subroutine xc_tabcd(what,l3d_dum,
     ,     tol_rho, Tmat, TTmat, Amat, Bmat, Cmat, Dmat, 
     N     Emat,Fmat,qxyz,xyz,
     &     chi, delchi, heschi, 
     N     curatoms,ncuratoms,nmat,
     I     ipol, nq, nbf, max_at_bf, max_at_bf2, 
     G     GRAD, basis, natoms, iniz, ifin, 
     &     g_vxc, ibf, rchi_atom, rdelchi_atom, 
     &     rdens_atom, cetobfr,kske,Mmat,kslap,Gmat,Lmat,scr)
c
c     We are using xc_tabcd for CPKS purposes
c
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "dftpara.fh"
#include "dft2drv.fh"
#include "dist.fh"
#include "dft_fdist.fh"
#include "util.fh"
c
      Logical GRAD !< [Input] .True. if functional depends on density
                   !< gradient
      integer what !< [Input] What should be calculated:
                   !< - what=0: Do everything (?)
                   !< - what=1: CPKS_LHS
                   !< - what=2: CPKS_RHS
                   !< - what=3: NMR_RHS
      integer basis !< [Input] The basis set handle
      integer max_at_bf !< [Input] The maximum number of basis functions
                        !< on an atom
      integer max_at_bf2 !< [Input] The maximum number of basis
                         !< functions on an atom
      integer nmat !< [Input] Number of XC matrices (alpha + beta sets)
                   !< to compute, e.g. the number of perturbations for
                   !< CPSCF.
      integer ipol  !< [Input] The number of spin channels
      integer nq    !< [Input] The number of grid points
      integer nbf    !< [Input] The number of basis functions
      integer natoms !< [Input] The number of atoms
      integer ncuratoms !< [Input] The number of current "active" atoms
      integer curatoms(*) !< [Input] Mapping array for current atoms
      integer jatcur,  nu, nu1,  indT
      double precision tol_rho !< [Input] The electron density threshold
      integer g_vxc(*)   !< [Input] The GA handle for Fock matrix
                         !< derivatives
      double precision Tmat(*) !< [Scratch]
      double precision TTmat(*) !< [Scratch]
      double precision rchi_atom(natoms) !< [Input] The maximum basis
                                         !< function radius for each
                                         !< atom in this call
      double precision rdelchi_atom(natoms) !< [Input] The maximum 
      !< radius of the basis function gradient for each atom in 
      !< this call
      double precision rdens_atom(natoms,natoms,ipol) !< Not used
      integer cetobfr(2,natoms) !< [Input] Mapping from center to 
      !< basis functions (how different from `iniz` and `ifin`?):
      !< - cetobfr(1,*): first basis function
      !< - cetobfr(2,*): last basis function
      logical l3d_dum !< [Junk]
c
c     Sampling Matrices for the XC Potential & Energy
c
      double precision Amat(nq,ipol,*) !< [Input] The derivative wrt rho
      double precision Cmat(nq,3,ipol,*) !< [Input] The derivative wrt 
                                         !< rgamma
      double precision Mmat(nq,ipol,*) !< [Input] The derivative wrt tau
      double precision Lmat(nq,ipol,*) !< [Input] The derivative wrt lap
      double precision Gmat(nq,max_at_bf) !< [Scratch]
c
c     nmr
      double precision Emat(nq,max_at_bf) !< [Scratch] NMR only
      double precision Fmat(nq,3,max_at_bf) !< [Scratch] NMR only
      double precision qxyz(3,*) !< [Input] Grid point coordinates
      double precision xyz(3,*)  !< [Input] Nuclear coordinates
      logical kske !< [Input] .True. if functional depends on kinetic
                   !< energy density
      logical kslap !< [Input] .True. if functional depends on laplacian
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
      integer imat ! XC matrix loop index
      double precision Bmat(nq,max_at_bf) !< [Scratch]
      double precision Dmat(nq,3,max_at_bf) !< [Scratch]
      integer iniz(natoms) !< [Input] The first basis function for each
                           !< atom
      integer ifin(natoms) !< [Input] The last basis function for each
                           !< atom
c
c     Basis Functions & Gradients
c
      double precision chi(nq,nbf) !< [Input] The value of the basis
                                   !< functions at the grid points
      double precision delchi(nq,3,nbf) !< [Input] The value of the 
                                        !< gradient of the basis
                                        !< functions at the grid points
      double precision heschi(nq,6,*) !< [Input] The value of the 
                                      !< Hessian of the basis
                                      !< functions at the grid points
      integer ibf(nbf) !< [Input] The rank of the basis function for
                       !< every basis function (why do we need this?)
      double precision A_MAX, C_MAX, AC_MAX, FUNC_MAXI,
     &                 B_MAX, D_MAX, BD_MAX, FUNC_MAXJ
      integer iat, inizia, ifinia, nbfia, nnia, ifirst, ilast
      integer jat, inizja, ifinja, nbfja, nnja, jfirst, jlast
      integer ii, mu, mu1
      integer n,lastjat
      double precision chi1
      double precision scr(nq) !< [Scratch] Temporary stuff in
                               !< quadrature
      double precision dabsmax
      double precision tolrho15
      external dabsmax
      logical l3d !< [Input] .True. if XC-matrices stored in a 3D GA
      integer jrsh,jrsh2,n3d,idir,jdir
      integer g_update(2)
cnmr
      double precision Rij(3) ! vector R(jat) - R(iat)
      integer inia, iq, ix, ix1, ix2
c
c     Indexing array for basis function hessian columns as if
c     it were a 3x3 matrix
c
      integer indh(3,3)
      logical w01,w02,w013,dofull
      double precision raa,rbb
      integer nbhandl1,nbhandl2
      logical nbfirst1,nbfirst2,doitt
      integer sizeblk
#include "nwc_const.fh"
      integer nonzero(nw_max_atom),natleft,
     A     iat0,jat0
      data indh / 1, 2, 3,
     &            2, 4, 5,
     &            3, 5, 6 /
      data nbhandl1 /0./
      data nbhandl2 /0./
      save nbhandl1
      save nbhandl2
c
c         0: l3d=.f.    & n3d=1
ccc     rhs: l3d=.true. & n3d=3
ccc     lhs: l3d=.true. & n3d=1
c
      call starttimer(monitor_tabcd)
      call starttimer(monitor_screen0)
c     lingering nbacc from previous calls
      if(nbhandl1.ne.0) call ga_nbwait(nbhandl1)
      if(nbhandl2.ne.0) call ga_nbwait(nbhandl2)
      l3d = ga_ndim(g_vxc).eq.3
      natleft=0
      do  iat = 1, natoms
        if (iniz(iat).ne.0) then
          natleft=natleft+1
          nonzero(natleft)=iat
        endif
      enddo
      tolrho15=tol_rho**1.25d0
      if (what.eq.0) then
        n3d=1
      elseif (what.eq.1) then
        n3d=1
      elseif (what.eq.2) then
        n3d=3
      elseif (what.eq.3) then
        n3d=3
      else
        call errquit(' wrong what value for xctabcd ',0,0)
      endif
      w01=what.eq.0.or.what.eq.1
      w013=w01.or.what.eq.3
      w02=what.eq.0.or.what.eq.2
      dofull=what.ne.0.or.l3d
      nbfirst1=.true.
      nbfirst2=.true.
      call endtimer(monitor_screen0)
c
c     Beginning of loop over multiple XC matrices
c
      do 500 imat = 1,nmat
c
c     Compute the matrix product for the XC potential and energy:
c
c              T = transpose(A*B) + transpose(C*D)
c

        call starttimer(monitor_screen1)

        A_MAX = dabsmax(nq*ipol,Amat(1,1,imat))
        if (GRAD) then
          C_MAX = dabsmax(nq*3*ipol,Cmat(1,1,1,imat))
        else
          C_MAX = 0d0
        endif
        AC_MAX = max(A_MAX,C_MAX)

        call endtimer(monitor_screen1)

c
c     repl stuff
c
        if (xcreplicated.and.dorepxc) then
          g_update(1)=k_repxc(1)
          g_update(2)=k_repxc(2)
        else
          g_update(1)=g_vxc(1)
          g_update(2)=g_vxc(2)
        endif
        do 430 iat0=1,natleft
          call starttimer(monitor_screen2)
          iat=nonzero(iat0)
          inizia = iniz(iat)
          FUNC_MAXI = rchi_atom(iat)
          if(GRAD) FUNC_MAXI = max(FUNC_MAXI,rdelchi_atom(iat))
          doitt=(AC_MAX*FUNC_MAXI).ge.tol_rho
          call endtimer(monitor_screen2)
          if (what.eq.2.or.doitt) then
            ifinia = ifin(iat)
            ifirst = cetobfr(1,iat)
            ilast = cetobfr(2,iat)
            nnia = ifinia - inizia + 1
            nbfia = ilast - ifirst + 1
            do ii = 1, ipol
              do mu = 1, nnia

                call starttimer(monitor_mult1)

                mu1 = mu+inizia-1
                if (GRAD) then
                  do n = 1, nq
                    chi1 = chi(n,mu1)


                    Bmat(n,mu) = Amat(n,ii,imat)*chi1 +
     &                   delchi(n,1,mu1)*Cmat(n,1,ii,imat) +
     &                   delchi(n,2,mu1)*Cmat(n,2,ii,imat) +
     &                   delchi(n,3,mu1)*Cmat(n,3,ii,imat)
                    Dmat(n,1,mu) = Cmat(n,1,ii,imat)*chi1
                    Dmat(n,2,mu) = Cmat(n,2,ii,imat)*chi1
                    Dmat(n,3,mu) = Cmat(n,3,ii,imat)*chi1
                    
                    if(kske) then 
                      Dmat(n,1,mu) = Dmat(n,1,mu) +
     &                     Mmat(n,ii,imat)*delchi(n,1,mu1)
                      Dmat(n,2,mu) = Dmat(n,2,mu) +
     &                     Mmat(n,ii,imat)*delchi(n,2,mu1)
                      Dmat(n,3,mu) = Dmat(n,3,mu) +
     &                     Mmat(n,ii,imat)*delchi(n,3,mu1)
                    endif

                    if(kslap) then
                      Dmat(n,1,mu) = Dmat(n,1,mu) +
     &                     2d0*Lmat(n,ii,imat)*delchi(n,1,mu1)
                      Dmat(n,2,mu) = Dmat(n,2,mu) +
     &                     2d0*Lmat(n,ii,imat)*delchi(n,2,mu1)
                      Dmat(n,3,mu) = Dmat(n,3,mu) +
     &                     2d0*Lmat(n,ii,imat)*delchi(n,3,mu1)
                      Bmat(n,mu) = Bmat(n,mu) + Lmat(n,ii,imat)*
     &                (heschi(n,1,mu1)+heschi(n,4,mu1)+heschi(n,6,mu1))
                      Gmat(n,mu) = Lmat(n,ii,imat)*chi1
                    endif

                    
                  enddo
                else
                  do n = 1, nq
                    Bmat(n,mu) = chi(n,mu1)*Amat(n,ii,imat)
                  enddo
                endif ! GRAD
              enddo ! mu
c     Monitoring 

              call endtimer(monitor_mult1)

c
              call starttimer(monitor_screen3)
              B_MAX = dabsmax(nnia*nq,Bmat)
              if (GRAD) then
                D_MAX = dabsmax(nnia*nq*3,Dmat)
              else
                D_MAX = 0d0
              endif
              BD_MAX = max(B_MAX,D_MAX)
c
              lastjat=iat0
              if (what.eq.2) lastjat=natleft
              if (what.eq.3) lastjat=iat0-1
              call endtimer(monitor_screen3)
              do 168 jat0=1,lastjat
                jat=nonzero(jat0)
                if(what.eq.2) then
c
c     To fit better into existing structure, loop over full square
c     of atom pairs and only compute nuclear derivative contribution
c     from jat.  Also, this way we only need check jatcur once and
c     for all, and do not have to check iatcur at all.
c
                  jatcur = curatoms(jat)
                  if (jatcur.eq.0) goto 168
                endif
                call starttimer(monitor_screen4)
                inizja = iniz(jat)
                FUNC_MAXJ = rchi_atom(jat)
                if(grad) FUNC_MAXJ = max(rdelchi_atom(jat),FUNC_MAXJ)
                doitt=(BD_MAX*FUNC_MAXJ).ge.tol_rho
                call endtimer(monitor_screen4)
                if (what.eq.2.or.doitt) then
 
c     Monitoring 

                  call starttimer(monitor_mult2)

                  if (what.eq.3) then
                    Rij(1) = 0.5d0*(xyz(1,jat)-xyz(1,iat))
                    Rij(2) = 0.5d0*(xyz(2,jat)-xyz(2,iat))
                    Rij(3) = 0.5d0*(xyz(3,jat)-xyz(3,iat))
                  endif
                  ifinja = ifin(jat)
                  jfirst = cetobfr(1,jat)
                  jlast = cetobfr(2,jat)
                  nbfja = jlast - jfirst + 1
                  nnja = ifinja - inizja + 1
                  sizeblk=n3d*nbfia*nbfja
                  if (what.eq.2.or.what.eq.3) 
     Y              call dcopy(sizeblk, 0d0,0, TTmat,1) 
c
c              Loop over x, y, z directions for derivative XC mats
c
                  do jdir = 1,n3d
c
                    if (what.eq.3) then
                      ix1 = mod(jdir,3)+1
                      ix2 = mod(jdir+1,3)+1
                      raa=rij(ix1)
                      rbb=rij(ix2)
                      do iq = 1, nq
                        scr(iq) = raa*qxyz(ix2,iq) - rbb*qxyz(ix1,iq)
                      enddo
                      do inia = 1, nnia
                        do iq = 1, nq
                          Emat(iq,inia) = scr(iq)*Bmat(iq,inia)
                        enddo
                        if (GRAD) then
                          do iq = 1, nq
                            Emat(iq,inia) = Emat(iq,inia)+
     &                           (raa*Dmat(iq,ix2,inia)
     &                           -  rbb*Dmat(iq,ix1,inia))
                          enddo
                        endif
                      enddo
                      if (GRAD) then
                        do inia = 1, nnia
                          do iq = 1, nq
                             Fmat(iq,1,inia) = scr(iq)
     &                            * Dmat(iq,1,inia)
                             Fmat(iq,2,inia) = scr(iq)
     &                            * Dmat(iq,2,inia)
                             Fmat(iq,3,inia) = scr(iq)
     &                            * Dmat(iq,3,inia)
                          enddo
                        enddo
                      endif
                    endif ! what.eq.3
c
c Daniel (2-7-13): Here, we build matrix elements of the XC-kernel when
c what = 1 for CPKS RHS with fxc*drhonuc*chi_mu in Bmat.
                    if (w01) then
                      call dgemm('T', 'N', nnia, nnja, nq, 1.d0, Bmat,
     &                     nq, chi(1,inizja), nq, 0.d0, Tmat, nnia)
                    elseif (what.eq.3) then
                      call dgemm('T', 'N', nnia, nnja, nq, 1.0d0, Emat,
     &                     nq, chi(1,inizja), nq, 0.0d0, Tmat, nnia)
                    else
c     Note the sign change for a nuclear derivative, and also that the
c     leading dimension of delchi must be set correctly
                      call dgemm('T', 'N', nnia, nnja, nq, -1.d0, Bmat,
     &                     nq, delchi(1,jdir,inizja), nq*3, 0.d0, Tmat,
     &                     nnia)
                    endif
                    if (GRAD) then
                      if (w01) then
                        call dgemm('T', 'N', nnia, nnja, 3*nq,
     &                       1.d0, Dmat, 3*nq, delchi(1,1,inizja), 
     &                       3*nq, 1.d0, Tmat, nnia)
                      elseif (what.eq.3) then
                        call dgemm('T', 'N', nnia, nnja, 3*nq,
     &                       1.0d0, Fmat, 3*nq, delchi(1,1,inizja), 
     &                       3*nq, 1.0d0, Tmat, nnia)
                      else
                        indT = 0
                        do nu = 1, nnja
                          nu1 = nu+inizja-1
                          do mu = 1, nnia
                            indT = indT + 1
                            Tmat(indT) = Tmat(indT)-
     *                           ddot(nq,Dmat(1,1,mu),1,
     &                           heschi(1,indh(1,jdir),nu1),1) -
     *                           ddot(nq,Dmat(1,2,mu),1,
     &                           heschi(1,indh(2,jdir),nu1),1) -
     *                           ddot(nq,Dmat(1,3,mu),1,
     &                           heschi(1,indh(3,jdir),nu1),1)
                          enddo
                        enddo
                      endif
                    endif ! GRAD
                    if(kslap) then
                       if (w01) then
                          indT=0
                          do nu = 1, nnja
                             nu1 = nu+inizja-1
                             do mu = 1, nnia
                               indT = indT + 1
                               Tmat(indT) = Tmat(indT) +
     &                         ddot(nq,Gmat(1,mu),1,heschi(1,1,nu1),1) +
     &                         ddot(nq,Gmat(1,mu),1,heschi(1,4,nu1),1) +
     &                         ddot(nq,Gmat(1,mu),1,heschi(1,6,nu1),1)
                             enddo
                          enddo
                       endif
                    endif
c Daniel (2-7-13): For the CPKS RHS stuff, the first call has what=1,
c so n3d=1.  Also, dofull is true in that case.
                    if (n3d.eq.1) then
                      call dfill(max_at_bf2, 0.d0, TTmat, 1)
                      if (dofull) then
                        call scat_mat(TTmat, Tmat, nbfia, nbfja, nnia,
     &                       nnja,ifirst,jfirst,ibf(inizia),ibf(inizja))
                      else
                        call scat_matup(TTmat, Tmat, nbfia, nbfja, nnia,
     &                       nnja,ifirst,jfirst,ibf(inizia),ibf(inizja))
                      endif
c Daniel (2-7-13): For the CPKS RHS stuff, the second call has what=2,
c so n3d=3.  Also, dofull is true in that case.
                    else
                      call scat_mat3(n3d,jdir,
     &                     TTmat, Tmat, nbfia, nbfja, nnia,
     &                     nnja,ifirst,jfirst,ibf(inizia),ibf(inizja))
                    endif
                       
                  enddo            ! jdir (loop over x, y, z directions for nmr)
c     Monitoring 

                  call endtimer(monitor_mult2)
              
                  doitt=.true.
                  call starttimer(monitor_screen5)
                  if (what.eq.0) then
                    doitt=dabsmax(sizeblk,ttmat).ge.tolrho15
                    jrsh=ii
                  elseif (what.eq.1) then
                    doitt=.true.
                    jrsh=imat+(ii-1)*nmat
                  elseif (what.eq.3) then
                    jrsh=(ii-1)*n3d+1
                  else
                    doitt=dabsmax(sizeblk,ttmat).ge.tol_rho
                    jrsh=1+(jat-1)*3+(ii-1)*3*natoms
                  endif
                  call endtimer(monitor_screen5)
                  if (doitt) then

                    jrsh2=jrsh+n3d-1

c     Monitoring

                    call updist(monitor_size_ga_acc1, sizeblk)
                    call starttimer( monitor_comm_ga_acc1)

c Daniel (2-7-13): l3d is always true for the CPKS RHS stuff.
                    if (l3d) then
                      call dft_3dacc(g_vxc, ttmat,
     &                     jrsh,jrsh2,
     %                     ifirst, ilast, jfirst, jlast, nbfia)
                    else
                      if (dftnbacc) then
                        if (.not.nbfirst1) then
                          call starttimer( monitor_wait1)
                          call ga_nbwait(nbhandl1)
                          call endtimer( monitor_wait1)
                        endif
                        nbfirst1=.false.
                        call upd_atombl_nb(g_update(ii),
     .                       basis,iat,jat,ttmat,nbhandl1)
                      else
                        if (truerepxc) then
                           call xc_atom_blockd(dbl_mb(k_repxc(ii)),
     N                          nbf_ld,basis,iat,jat,ttmat)
                        else
                          if (truerepxc) then
                            call xc_atom_block(dbl_mb(k_repxc(ii)),
     N                           nbf_ld,basis,jat,iat,tmat)
                          else
                            call upd_atom_block(g_update(ii),
     .                           basis,iat,jat,ttmat)
                          endif
                        endif
                      endif
                    endif ! l3d
c     Monitoring

                    call endtimer( monitor_comm_ga_acc1)

c Daniel (2-7-13): l3d is always true for the CPKS RHS stuff.
                    if (what.ne.0.or.l3d) then
c
c                 check to see if can skip and use ga_symmetrize
c
c Daniel (2-7-13): This part is always done for the second call to
c xc_tabcd in CPKS RHS stuff.  For the first call, it only happens
c for off-diagonal terms.
                      if ((w013.and.iat.ne.jat).or.what.eq.2) then
c     For CPKS RHS, we update with transpose even for iat = jat,
c     since that is necessary to get both contributions
c     mu * del(nu) and del(mu) * nu
                     

                        call starttimer(monitor_comp_transp)

c Daniel (2-7-13): This happens for the first call for CPKS RHS stuff
c (n3d=1).
                        if (n3d.eq.1) then
                          call transp_mat(TTmat, Tmat,
     ,                         nbfia, nbfja)
c Daniel (2-7-13): This happens for the second call for CPKS RHS stuff
c (n3d=3).
                        else
                          if (what.eq.3) then
                            call dscal(n3d*nbfia*nbfja,-1.0d0,TTmat,1)
                          endif
                          call transp_mat3(n3d,TTmat, Tmat,
     ,                         nbfia, nbfja)
                        endif

                        call endtimer(monitor_comp_transp)

                     
c     Monitoring 

                        call starttimer(monitor_comm_ga_acc2)


c Daniel (2-7-13): For CPKS RHS stuff, l3d=true always.                       
                        if (l3d) then
                          call dft_3dacc(g_vxc, tmat,
     &                         jrsh,jrsh2,
     %                         jfirst, jlast, ifirst, ilast, nbfja)
                        else
                          if (dftnbacc) then
                            if (.not.nbfirst2) then
                              call ga_nbwait(nbhandl2)
                            endif
                            nbfirst2 = .false.
                            call upd_atombl_nb(g_update(ii),
     .                           basis,jat,iat,tmat,nbhandl2)
                          else
                            call upd_atom_block(g_update(ii),basis,
     J                           jat,iat,tmat)
                          endif
                        endif
c     Monitoring 

                        call endtimer(monitor_comm_ga_acc2)

                       
                      endif ! (w013.and.iat.ne.jat).or.what.eq.2
                    endif ! what.ne.3.or.l3d
                  endif ! doitt
                endif ! what.eq.2.or.doitt
  168         continue ! jat0 loop
            enddo ! ipol loop
          endif ! what.eq.2.or.doitt
  430   continue ! iat0 loop
  500 continue ! imat loop
      call endtimer(monitor_tabcd)

      return
      end
      subroutine upd_atombl_nb(g_array, basis, iat, jat, buf,
     $     nbhandle)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "bas.fh"
c
      integer g_array, basis, iat, jat
      integer nbhandle
      double precision buf(*)
      logical status
c
      integer ilo, ihi, jlo, jhi, idim, jdim
c
c     add atom block buf info of the matrix g_array (over basis functions)
c
      status= bas_ce2bfr(basis, iat, ilo, ihi)
      status=status.and. bas_ce2bfr(basis, jat, jlo, jhi)
      if (.not. status)
     $     call errquit('upd_atom_block: ce2bfr failed', 0, BASIS_ERR)
c
      idim = ihi - ilo + 1
      jdim = jhi - jlo + 1
c     
c     clearing lingering nbacc
c
c      if(nbhandle.ne.0) call ga_nbwait(nbhandle)
      nbhandle=0
      if (idim.gt.0 .and. jdim.gt.0) 
     $     call ga_nbacc(g_array, ilo, ihi, jlo, jhi, buf, idim,
     $        1.0d0,nbhandle)
c
      end
C>
C> @}
