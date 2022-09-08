#define STAG 1
C> \brief Evaluate the density and their derivatives from the fitting
C>        basis on a grid
C>
C> 
      Subroutine xc_rhogen_adft(
     T     tol_rho, natoms, ipol, nq, nbf, mbf, max_at_bf, GRAD, ipol2,
     C     chi, delchi, heschi, ibf, iniz, ifin, rho, delrho, lap,
     &     rchi_atom, rdelchi_atom, cetobfr, wmax, dolap, cd_coef)
      implicit none
#include "global.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "dftpara.fh"
#include "dist.fh"
#include "dft_fdist.fh"

      Logical GRAD !< [Input] .true. when using gradient corrected
                   !< functional
      logical dolap       !< [Input] .true. if the Laplacian of the 
                           !< density is required
      integer ipol         !< [Input] no. of spin channels
      integer ipol2        !< [Input] no. of spin channels in density
                           !< - 1 if closed shell
                           !< - 3 if open shell
      integer max_at_bf
      integer nbf
      integer mbf          !< [Input] "restricted" no. of basis functions
      integer nq           !< [Input] no. of quadrature points
      integer natoms       !< [Input] no. of atoms
      double precision wmax !< [Input] max weight
      integer ibf(mbf)     !< [Input] mapping of nbf_ao -> mbf
      integer iniz(natoms) !< [Input] mapping of nbf_ao -> mbf
      integer ifin(natoms) !< [Input] mapping of nbf_ao -> mbf
      double precision tol_rho !< [Input] accuracy for rho evaluation
      double precision chi(nq,mbf)     !< [Input] function values
      double precision delchi(nq,3,mbf)!< [Input] function gradients
      double precision heschi(nq,6,mbf)!< [Input] function hessians
      double precision lap(nq,ipol2)
      double precision delrho(nq,3,ipol) !< [Output] Derivative of density
      double precision rho(nq,ipol2) !< [Output] The density
      double precision rchi_atom(natoms) !< [Input] Screening parameters
      double precision rdelchi_atom(natoms) !< [Input] Screening parameters
      double precision cd_coef(nbf,ipol)
      double precision coeff(max_at_bf)
      integer cetobfr(2,natoms) !< [Input] Centers to basis functions
c
c     local declarations
c
      integer i0, ii, mu, n, npol
      integer iat, inizia, ifirst, ilast, nbfia, nnia, iat_in
      integer ifinia, ifinja
      integer jat, inizja, jfirst, jlast, nbfja, nnja
      double precision FUNC_MAX, DELFUNC_MAX, FUNC_MAXI, P_MAXI
      double precision dabsmax
      external dabsmax
      integer g_keepd(2)
      integer nbhandl
      integer jj
      integer nonzero,nonz0
      logical zapnegatives
#ifdef DEBUG
      integer ga_nodeid
      external ga_nodeid
#endif
      call starttimer(monitor_xcrho)
      zapnegatives=.true.
c
c
c     Evaluate the charge density and its gradient at each of the
c     sampling points 
c
      npol = (ipol*(ipol+1))/2
      call dcopy(nq*npol,0d0,0,rho,1)
      if (grad) call dcopy(3*nq*ipol,0d0,0,delrho,1)
      if (dolap) call dcopy(nq*ipol2,0d0,0,lap,1)  ! total
c
      i0=ipol-1
c
      FUNC_MAX = dabsmax(natoms,rchi_atom)
      DELFUNC_MAX=0d0
      if(grad) DELFUNC_MAX = dabsmax(natoms,rdelchi_atom)
c
#ifdef STAG
      do 230 iat_in = 1+ga_nodeid(), natoms+ga_nodeid()
         iat=mod(iat_in,natoms)
         if(iat.eq.0) iat=natoms
#else
      do 230 iat = 1, natoms
#endif
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
         FUNC_MAXI = rchi_atom(iat)
         if(grad) FUNC_MAXI = max(FUNC_MAXI,rdelchi_atom(iat))
         P_MAXI = dabsmax(nbfia,cd_coef(ifirst,1))
         if (ipol.gt.1) 
     $      P_MAXI = max(P_MAXI,dabsmax(nbfia,cd_coef(ifirst,2)))
         if (FUNC_MAXI*P_MAXI.lt.tol_rho) goto 230

         do ii=1,ipol
           do mu=inizia,ifinia
             coeff(mu-inizia+1) = cd_coef(ibf(mu),ii)
           enddo

           call dgemv('N',nq,nnia,1d0,chi(1,inizia),nq,
     $                 coeff,1,1d0,rho(1,i0+ii),1)

           if (grad) then
             call dgemv('N',3*nq,nnia,1d0,delchi(1,1,inizia),
     $                   3*nq,coeff,1,1d0,delrho(1,1,ii),1)
           endif

           if (dolap) then
             call dgemv('N',nq,nnia,1d0,heschi(1,1,inizia),6*nq,
     $                   coeff,1,1d0,lap(1,ii),1)
             call dgemv('N',nq,nnia,1d0,heschi(1,4,inizia),6*nq,
     $                   coeff,1,1d0,lap(1,ii),1)
             call dgemv('N',nq,nnia,1d0,heschi(1,6,inizia),6*nq,
     $                   coeff,1,1d0,lap(1,ii),1)
           endif

         enddo
 230  continue

      if(zapnegatives) then
!$omp   parallel do collapse(2)
        do ii = 1, ipol
          do n = 1, nq
            if (rho(n,ii+i0).le.0.0d0) then
              rho(n,ii+i0) = 0.0d0
              if (grad) then
                delrho(n,1,ii) = 0.0d0
                delrho(n,2,ii) = 0.0d0
                delrho(n,3,ii) = 0.0d0
              endif
            endif
          enddo
        enddo
!$omp   end parallel do
      endif
c
      call starttimer(monitor_rhocomp2)
      if (ipol.eq.2)then
         call dcopy(nq, rho(1,2), 1, rho(1,1), 1)
         call daxpy(nq, 1.d0, rho(1,3), 1, rho(1,1), 1)
      endif
      call endtimer(monitor_rhocomp2)
c     
      call endtimer(monitor_xcrho)
      return
      end
