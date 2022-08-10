C>
      Subroutine xc_tabcd_adft(
     ,     tol_rho, Amat, Cmat, Lmat,
     &     chi, delchi, heschi, 
     N     nmat,
     I     ipol, nq, nbf, max_at_bf, 
     G     GRAD, basis, natoms, iniz, ifin, 
     &     ibf, rchi_atom, rdelchi_atom, 
     &     cetobfr,kslap,xc_coef)
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
      integer basis !< [Input] The basis set handle
      integer max_at_bf !< [Input] The maximum number of basis functions
                        !< on an atom
      integer nmat !< [Input] Number of XC matrices (alpha + beta sets)
      integer ipol  !< [Input] The number of spin channels
      integer nq    !< [Input] The number of grid points
      integer nbf    !< [Input] The number of basis functions
      integer natoms !< [Input] The number of atoms
      double precision xc_coef(nbf,ipol)
      double precision tol_rho !< [Input] The electron density threshold
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
c
c     Sampling Matrices for the XC Potential & Energy
c
      double precision Amat(nq,ipol,*) !< [Input] The derivative wrt rho
      double precision Cmat(nq,3,ipol,*) !< [Input] The derivative wrt 
      double precision Lmat(nq,ipol,*) !< [Input] The derivative wrt lap
c
      logical kslap !< [Input] .True. if functional depends on laplacian
c
      integer imat ! XC matrix loop index
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
      double precision dabsmax
      external dabsmax
      logical l3d !< [Input] .True. if XC-matrices stored in a 3D GA
      integer jrsh,jrsh2,n3d,idir,jdir
cnmr
      integer inia, iq, ix, ix1, ix2
c
      integer nbhandl1,nbhandl2
      logical nbfirst1,nbfirst2,doitt
      integer sizeblk
#include "nwc_const.fh"
      integer nonzero(nw_max_atom),natleft,
     A     iat0,jat0
      double precision coeff(max_at_bf)
      double precision ddot
      external ddot
c
c         0: l3d=.f.    & n3d=1
ccc     rhs: l3d=.true. & n3d=3
ccc     lhs: l3d=.true. & n3d=1
c
      call starttimer(monitor_tabcd)
c     lingering nbacc from previous calls
      natleft=0
      do  iat = 1, natoms
        if (iniz(iat).ne.0) then
          natleft=natleft+1
          nonzero(natleft)=iat
        endif
      enddo
c
c     Beginning of loop over multiple XC matrices
c
      do 500 imat = 1,nmat
        call starttimer(monitor_screen1)
        A_MAX = dabsmax(nq*ipol,Amat(1,1,imat))
        if (GRAD) then
          C_MAX = dabsmax(nq*3*ipol,Cmat(1,1,1,imat))
        else
          C_MAX = 0d0
        endif
        AC_MAX = max(A_MAX,C_MAX)
        call endtimer(monitor_screen1)

        do 430 iat0=1,natleft
          call starttimer(monitor_screen2)
          iat=nonzero(iat0)
          inizia = iniz(iat)
          ifinia = ifin(iat)
          ifirst = cetobfr(1,iat)
          ilast = cetobfr(2,iat)
          nnia = ifinia - inizia + 1
          nbfia = ilast - ifirst + 1
          FUNC_MAXI = rchi_atom(iat)
          if(GRAD) FUNC_MAXI = max(FUNC_MAXI,rdelchi_atom(iat))
          doitt=(AC_MAX*FUNC_MAXI).ge.tol_rho
          call endtimer(monitor_screen2)

          if (.not.doitt) goto 430

          call starttimer(monitor_mult1)
          do ii = 1, ipol
            call dgemv('T',nq,nnia,1d0,chi(1,inizia),nq,amat(1,ii,imat),
     $                  1, 0d0, coeff,1)
            if (grad) then
              call dgemv('T',3*nq,nnia,1d0,delchi(1,1,inizia),3*nq,
     $                    cmat(1,1,ii,imat),1,1d0,coeff,1)
            endif
            if (kslap) then
              call dgemv('T',nq,nnia,1d0,heschi(1,1,inizia),6*nq,
     $                    lmat(1,ii,imat),1,1d0,coeff,1)         
              call dgemv('T',nq,nnia,1d0,heschi(1,4,inizia),6*nq,
     $                    lmat(1,ii,imat),1,1d0,coeff,1)         
              call dgemv('T',nq,nnia,1d0,heschi(1,6,inizia),6*nq,
     $                    lmat(1,ii,imat),1,1d0,coeff,1)         
            endif

            do mu=inizia,ifinia
              xc_coef(ibf(mu),ii) = xc_coef(ibf(mu),ii) + 
     $                                coeff(mu-inizia+1)
            enddo
          enddo
          call endtimer(monitor_mult1)

  430   continue ! iat0 loop
  500 continue ! imat loop
      call endtimer(monitor_tabcd)

      return
      end
