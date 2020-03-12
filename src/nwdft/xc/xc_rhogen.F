#define STAG 1
* $Id$
c     Modified for construction of multiple sets of perturbed densities
c     and explicit derivative densities wrt nuclei
c
c     BGJ - 8/98
c
C> \ingroup nwdft_xc
C> @{
C>
C> \file xc_rhogen.F
C> Evaluate density and derivatives
C>
C> @}
C>
C> \ingroup nwdft_xc
C> @{
C>
C> \brief Evaluate the density and their derivatives on a grid
C>
C> 
      Subroutine xc_rhogen(what_in,
     T     tol_rho, basis, g_dens, max_at_bf, 
     N     natoms, curatoms, ncuratoms, npert,       
     I     ipol, nq, nbf, mbf, GRAD, ipol2,
     &     F, Pmat, ff, ffd,
     C     chi, delchi, heschi,
     I     ibf, iniz, ifin, 
     &     rho, delrho, lap, rchi_atom,
     &     rdelchi_atom, rdens_atom, cetobfr, wmax, 
     &     ttau, kske, dolap)
      implicit none
#include "global.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "dftpara.fh"
#include "dist.fh"
#include "dft_fdist.fh"

      Logical GRAD !< [Input] .true. when using gradient corrected
                   !< functional
      Logical kske !< [Input] .true. when using kinetic energy density
                   !< functional
      integer what_in !< [Input] What to calculate
                   !< - what=0: calculate density
                   !< - what=1: calculate perturbed density (derivative
                   !<   with respect to electron position)
                   !< - what=2: calculate derivative with respect to
                   !<   nuclear coordinates
       logical dolap       !< [Input] .true. if the Laplacian of the 
                           !< density is required
c what=-1 dplot aka density contours
c what=0 dens
c what=1 pert
c what=2 nucder
      integer basis        !< [Input] basis set handle
      integer ipol         !< [Input] no. of spin channels
      integer ipol2        !< [Input] no. of spin channels in density
                           !< - 1 if closed shell
                           !< - 3 if open shell
      integer nbf          !< [Input] no. of basis functions
      integer mbf          !< [Input] "restricted" no. of basis functions
      integer max_at_bf    !< [Input] max no. bf per atom
      integer nq           !< [Input] no. of quadrature points
      integer natoms       !< [Input] no. of atoms
      double precision wmax !< [Input] max weight
      integer npert        !< [Input] number of perturbed densities
      integer curatoms(*) !< [Input] indexing array for "active" atoms
      integer ncuratoms        !< [Input] number of currently active atoms
      integer g_dens(*) !< [Input] GA handle for DM
      integer ibf(mbf)     !< [Input] mapping of nbf_ao -> mbf
      integer iniz(natoms) !< [Input] mapping of nbf_ao -> mbf
      integer ifin(natoms) !< [Input] mapping of nbf_ao -> mbf
      double precision tol_rho !< [Input] accuracy for rho evaluation
      double precision chi(nq,mbf)     !< [Input] function values
      double precision delchi(nq,3,mbf)!< [Input] function gradients
      double precision heschi(nq,6,mbf)!< [Input] function hessians
c
      double precision ttau(nq,ipol,*) !< [Output] Total Kohn-Sham K.E.density
c
       double precision lap(nq,ipol2,*)
c
      double precision delrho(nq,3,ipol,*) !< [Output] Derivative of density
      double precision Pmat(*) !< [Scratch] scratch vector
      double precision F(max_at_bf*max_at_bf) !< [Scratch] scratch vector
      double precision ff(nq,*)   !< [Scratch] scratch array
      double precision ffd(nq,3,*) !< [Scratch] scratch array
      double precision rho(nq,ipol2,*) !< [Output] The density
      double precision rchi_atom(natoms) !< [Input] Screening parameters
      double precision rdelchi_atom(natoms) !< [Input] Screening parameters
      double precision rdens_atom(natoms,natoms,ipol) !< [Input] Screening parameters
      integer cetobfr(2,natoms) !< [Input] Centers to basis functions
c
c     local declarations
c
      integer i0, ii, mu, n, npol
      integer ipert        ! perturbation loop index
      integer iat, inizia, ifirst, ilast, nbfia, nnia, iat_in
      integer ifinia, ifinja
      integer jat, inizja, jfirst, jlast, nbfja, nnja
      integer iatcur, jatcur
      double precision FUNC_MAX, DELFUNC_MAX, FUNC_MAXI, FUNC_MAXJ
      double precision P_MAXJ, P_MAXJ_A, P_MAXJ_B, P_MAXIJ
      double precision dabsmax
      integer g_keepd(2)
      integer nbhandl
      integer jj
      logical doffd,doitt
      external dabsmax
      external xc_rhoscreen
      integer xc_rhoscreen
      integer nonzero,nonz0
      integer i_nz,l_nz
      integer sizeblk, gindx
      integer what
      logical zapnegatives
#ifdef DEBUG
      integer ga_nodeid
      external ga_nodeid
#endif
      g_keepd(1)=0
      g_keepd(2)=0
      call starttimer(monitor_xcrho)
      if(what_in.ge.0) then
         what=what_in
         zapnegatives=.true.
      else
         what=-what_in
         zapnegatives=.false.
      endif
c
c     Evaluate the charge density and its gradient at each of the
c     sampling points 
c
      doffd=(what.eq.2.and.grad).or.(what.eq.0.and.kske).or.
     &      (what.eq.1.and.kske).or.(what.eq.0.and.dolap)
      
      npol = (ipol*(ipol+1))/2
c     to keep compilers quiet
      iatcur=1
      jatcur=1
c
      if(what.eq.0) then
         call dcopy(nq*npol,0d0,0,rho,1)
         if (grad) call dcopy(3*nq*ipol,0d0,0,delrho,1)
         if (kske) call dcopy(nq*ipol,0d0,0,ttau,1)  ! total
         if (dolap) call dcopy(nq*ipol2,0d0,0,lap,1)  ! total
c
c     repl stuff
c
         g_keepd(1)=g_dens(1)
         if(ipol.eq.2) g_keepd(2)=g_dens(2)
         if(xcreplicated.and.dorepdm) then
            g_dens(1)=g_repdm(1)
            if(ipol.eq.2) g_dens(2)=g_repdm(2)
         endif
      elseif(what.eq.1) then
         call dcopy(nq*ipol*npert,0d0,0,rho,1)
         if (grad) call dcopy(3*nq*ipol*npert,0d0,0,delrho,1)
         if (kske) call dcopy(nq*ipol*npert,0d0,0,ttau,1)   ! total
      elseif(what.eq.2) then
         call dcopy(nq*ipol*3*ncuratoms,0d0,0,rho,1)
         if (grad)call dcopy(3*nq*ipol*3*ncuratoms,0d0,0,delrho,1)
      else
         call errquit('wrong what value',0,0)
      endif
c
c     Screening is accomplished by:  p(r) <= |Xi(r)|*|Xj(r)|*|Dij|
c     Xi(r) is screened on desired accuracy/max(|Xj(r)|)*max(|Dij|)
c     Dij is screened on desired accuracy/max(|Xi(r)|)*max(|Xj(r)|)
c
      i0=ipol-1
c
      FUNC_MAX = dabsmax(natoms,rchi_atom)
      DELFUNC_MAX=0d0
      if(grad) DELFUNC_MAX = dabsmax(natoms,rdelchi_atom)
c
      nonzero=0
      if(dftnbget) then
         if (.not.ma_push_get
     &        (mt_int,(ipol*npert*natoms*(natoms+1))/2,
     N        'nzmap map',l_nz,i_nz))
     &        call errquit('xcrho:push_get failed', 13, MA_ERR)
         nonzero=xc_rhoscreen(grad,ipol,natoms,npert,
     I        iniz,
     W        tol_rho,wmax,
     O        int_mb(i_nz),
     R        rchi_atom,rdelchi_atom,rdens_atom)
         
         if(nonzero.eq.0) goto 1688
c
c     prefetch first DM block
c
         nonz0=1
         call xc_getdmblock(int_mb(i_nz),nonz0,natoms,cetobfr,
     G        g_dens(1),
     A        Pmat,nbhandl)
      endif

#ifdef STAG
      do 230 iat_in = 1+ga_nodeid(), natoms+ga_nodeid()
         iat=mod(iat_in,natoms)
         if(iat.eq.0) iat=natoms
#else
      do 230 iat = 1, natoms
#endif
         inizia = iniz(iat)
         if (inizia.eq.0)goto 230
         if(what.eq.2) then
            iatcur = curatoms(iat)
         endif
         ifinia = ifin(iat)
         ifirst = cetobfr(1,iat)
         ilast = cetobfr(2,iat)
         nbfia = ilast-ifirst+1
         nnia = ifinia-inizia+1
c
c        screening parameters
c
         FUNC_MAXI = rchi_atom(iat)
         if(grad)
     .     FUNC_MAXI = max(FUNC_MAXI,rdelchi_atom(iat))
         FUNC_MAXJ = max(FUNC_MAX,DELFUNC_MAX)
         if(what.lt.2) then
            if (ipol.gt.1)then
              P_MAXJ_A = dabsmax(natoms,rdens_atom(1,iat,1))
              P_MAXJ_B = dabsmax(natoms,rdens_atom(1,iat,2))
              P_MAXJ = MAX(P_MAXJ_A, P_MAXJ_B)
            else
              P_MAXJ=0d0
              do jat=1,iat
                 if(iniz(jat).ne.0) 
     .                P_MAXJ=max(P_MAXJ,rdens_atom(jat,iat,1))
              enddo
            endif
            if (wmax*FUNC_MAXI*FUNC_MAXJ*P_MAXJ.lt.tol_rho) goto 225
         endif
         do 220 jat = 1, iat
            inizja = iniz(jat)
            if (inizja.eq.0)goto 220
            if(what.eq.2) then
              jatcur = curatoms(jat)
              if (jatcur.eq.0.and.iatcur.eq.0) goto 220
            endif
            call starttimer(monitor_rscreen0)
            ifinja = ifin(jat)
            jfirst = cetobfr(1,jat)
            jlast = cetobfr(2,jat)
            nbfja = jlast-jfirst+1
            nnja = ifinja-inizja+1
c
c           screening parameters
c
            FUNC_MAXJ=rchi_atom(jat)
            if(grad)
     .        FUNC_MAXJ = max(FUNC_MAXJ,rdelchi_atom(jat))
            doitt=.true.
            if(what.lt.2) then
              P_MAXIJ = rdens_atom(iat,jat,1)
              if(ipol.eq.2) P_MAXIJ = max(P_MAXIJ,
     &             rdens_atom(iat,jat,2))
              doitt=(wmax*FUNC_MAXI*FUNC_MAXJ*P_MAXIJ.ge.tol_rho)
            endif
            call endtimer(monitor_rscreen0)
            if (doitt) then
              do ii = 1, ipol
c
c               screening parameters
c
                if((what.gt.1).or.(rdens_atom(iat,jat,ii)*
     R               wmax*FUNC_MAXI*FUNC_MAXJ.ge.tol_rho)) then
c
c                  Loop over perturbations
c
                   do 215 ipert = 1,npert
                     sizeblk=nbfia*nbfja
                     call updist(monitor_size_ga_get, sizeblk)
c
                     if(dftnbget) then
                       call starttimer(monitor_wait3)
                       call ga_nbwait(nbhandl)
                       call endtimer(monitor_wait3)
                       call gat_mat(F, Pmat, nbfia, nbfja, nnia, nnja, 
     I                      ifirst, jfirst, ibf(inizia), ibf(inizja))
                       nonz0=nonz0+1
                       if((npert*ipol).gt.1) then
                         gindx=ipert+(ii-1)*npert+1
                      
                         if(gindx.gt.npert*ipol)
     +                      gindx=mod(gindx,npert*ipol)
                       else
                         gindx=1
                       endif
                       if(nonz0.le.nonzero)
     C                   call xc_getdmblock(int_mb(i_nz),nonz0,natoms,
     C                   cetobfr,g_dens(gindx),
     A                      Pmat,nbhandl)
                       if (wmax*FUNC_MAXI*FUNC_MAXJ*
     .                     dabsmax(sizeblk,Pmat).lt.tol_rho)goto 215
                     else
                       call starttimer(monitor_gaget)
                       if(truerepdm) then
                         call xc_dmget(
     +                      dbl_mb(k_repdm(ipert+(ii-1)*npert)),nbf_ld,
     %                      ifirst, ilast, jfirst, jlast, Pmat,nbfia)
                       else
                         call ga_get(g_dens(ipert+(ii-1)*npert), 
     %                        ifirst, ilast, jfirst, jlast, Pmat,nbfia)
                       endif
                       call endtimer(monitor_gaget)
                       if (wmax*FUNC_MAXI*FUNC_MAXJ*
     .                     dabsmax(sizeblk,Pmat).lt.tol_rho)goto 215
                       call gat_mat(F, Pmat, nbfia, nbfja, nnia, nnja,
     +                      ifirst, jfirst, ibf(inizia), ibf(inizja))
                     endif
c
                     call starttimer(monitor_rhocomp)
                     if(iat.ne.jat) call dscal(nnia*nnja,2d0,F,1)
c
c                    Compute Xiat(r)*Xjat(r)*Diat,jat
c
                     call dgemm('n','n',nq,nnja,nnia,1d0,
     A                    chi(1,inizia),nq,F,nnia,0d0,ff,nq)
                     if(what.lt.2) then
                        jj=i0+ii
                        if(what.eq.1) jj=ii
                        do mu=inizja,ifinja
                           do n=1,nq
                              rho(n,jj,ipert) = rho(n,jj,ipert) + 
     P                             chi(n,mu)*ff(n,mu-inizja+1)
                           enddo
                        enddo
                     endif
                     if(doffd) then
                       call dgemm('n','n',nq*3,nnja,nnia,1d0,
     A                      delchi(1,1,inizia),nq*3,F,nnia,0d0,ffd,nq*3)
c
c                      build tau for meta GGA
c
                       if(kske) then
                         do mu=inizja,ifinja
                           do n=1,nq
                             ttau(n,ii,ipert) = ttau(n,ii,ipert)+0.5d0*(
     &                         delchi(n,1,mu)*ffd(n,1,mu-inizja+1)+
     &                         delchi(n,2,mu)*ffd(n,2,mu-inizja+1)+
     &                         delchi(n,3,mu)*ffd(n,3,mu-inizja+1))
                           enddo
                         enddo
                       endif

                       if (dolap) then
                         do mu=inizja,ifinja
                           do n=1,nq
                              ! total
                              lap(n,ii,ipert) = lap(n,ii,ipert) + 2d0*
     &                            (delchi(n,1,mu)*ffd(n,1,mu-inizja+1) +
     &                             delchi(n,2,mu)*ffd(n,2,mu-inizja+1) +
     &                             delchi(n,3,mu)*ffd(n,3,mu-inizja+1))
cold     &                             delchi(n,3,mu)*ffd(n,3,mu-inizja+1) +
cold     &                             heschi(n,1,mu)*ff(n,mu-inizja+1) +
cold     &                             heschi(n,4,mu)*ff(n,mu-inizja+1) +
cold     &                             heschi(n,6,mu)*ff(n,mu-inizja+1))
                           enddo
                         enddo
                       endif
                     endif
c
                     if((what.eq.2.and.jatcur.ne.0).or.
c
c                      We need the "grad" code at zero order in the
c                      nuclear deriv case, but we can skip this part if
c                      iat is not active
c     
     O                  (what.ne.2.and.grad)) then
                       if(what.eq.0.or.what.eq.1) then
                         call xc_dchiff(nq,inizja,ifinja,
     P                        delrho(1,1,ii,ipert),delchi,
     F                        ff)
cc DMR/Begin
                         if (dolap) then
                            do mu=inizja,ifinja
                               do n=1,nq
                                  lap(n,ii,ipert) = lap(n,ii,ipert) +
     &                            ff(n,mu-inizja+1)*(heschi(n,1,mu) +
     &                            heschi(n,4,mu) + heschi(n,6,mu))
                               enddo
                            enddo
                         endif
cc DMR/End
                       elseif(what.eq.2) then
                         call xc_dchiffp(nq,ipol2,inizja,ifinja,
     P                        rho,delchi,ff,
     M                        ii,jat)
                         if (grad) then
c     
c                          Compute nuclear gradient of delrho
c     
                           call xc_drhonuc(nq,ipol,inizja,ifinja,
     D                          delrho,heschi,delchi,ff,ffd,
     I                          ii,jat)
                         endif
                       endif
                     endif
c
                     if((what.eq.2.and.iatcur.ne.0).or.
     O                  (what.ne.2.and.grad)) then
c
c                      Compute delXiat(r)*Xjat(r)*Diat,jat
c                  
                       call dgemm('n','t',nq,nnia,nnja,1d0,
     A                      chi(1,inizja),nq,F,nnia,0d0,ff,nq)
                       if(what.lt.2) then
                         call xc_dchiff(nq,inizia,ifinia,
     P                        delrho(1,1,ii,ipert),delchi,
     F                        ff)
cc DMR/Begin
                         if (dolap) then
                            do mu=inizia,ifinia
                               do n=1,nq
                                  lap(n,ii,ipert) = lap(n,ii,ipert) +
     &                              ff(n,mu-inizia+1)*(heschi(n,1,mu) +
     &                              heschi(n,4,mu) + heschi(n,6,mu))
                               enddo
                            enddo
                         endif
cc DMR/End
                       elseif(what.eq.2) then
                         call xc_dchiffp(nq,ipol2,inizia,ifinia,
     P                        rho,delchi,ff,
     M                        ii,iat)
                         if(grad) then
                           call dgemm('n','t',nq*3,nnia,nnja,1d0,
     A                          delchi(1,1,inizja),nq*3,F,nnia,0d0,ffd,
     +                          nq*3)
c     
c                          Compute nuclear gradient of delrho
c     
                           call xc_drhonuc(nq,ipol,inizia,ifinia,
     D                          delrho,heschi,delchi,ff,ffd,
     I                          ii,iat)
                        
                         endif
                       endif
                     endif
                     call endtimer(monitor_rhocomp)
  215              continue
                endif
              enddo
            endif
  220    continue
  225    continue
  230 continue
      if(zapnegatives) then
c
c     Enforce results that are compatible with the laws of physics.
c     I.e. Rho >= 0, Tau >= 0, Rho=0 ==> Grad(Rho)=0, and
c     Tau=0 ==> Grad(Rho)=0. The reason for doing this is that density
c     functionals have been designed to be valid for physically 
c     sensible data points. They may fail spectacularly for data points
c     outside the physically sensible range. Alternatively to screening
c     the data here it might be done in the functionals themselves. 
c     However, this requires a great deal of care to make sure that the
c     right limits are taken.
c
c     Rho should be non-negative. Numerical errors in the construction
c     may leave small negative values that are formally incorrect.
c     Hence filter negative values out.
c
c     Also note that rho.eq.0 implies delrho=0, so this is enforced as
c     well. Proof: Rho is a non-negative quantity. Hence if rho.eq.0
c     then rho is at a minimum. A minimum being an extremum it follows
c     that its gradient must be zero. QED.
c
      if (what.eq.0) then
        do ii = 1, ipol
          jj = i0 + ii
          do n = 1, nq
            if (rho(n,jj,1).le.0.0d0) then
              rho(n,jj,1) = 0.0d0
              if (grad) then
                delrho(n,1,ii,1) = 0.0d0
                delrho(n,2,ii,1) = 0.0d0
                delrho(n,3,ii,1) = 0.0d0
              endif
            endif
          enddo
        enddo
      else if (what.eq.1) then
        do ii = 1, ipol
          jj = ii
          do n = 1, nq
            if (rho(n,jj,npert).le.0.0d0) then
              rho(n,jj,npert) = 0.0d0
              if (grad) then
                delrho(n,1,ii,npert) = 0.0d0
                delrho(n,2,ii,npert) = 0.0d0
                delrho(n,3,ii,npert) = 0.0d0
              endif
            endif
          enddo
        enddo
      endif
c
c     Tau should be non-negative. Numerical errors in the construction
c     may leave small negative values that can cause trouble in the
c     density functionals. Hence filter negative values out.
c
c     Also note that tau.eq.0 implies delrho=0, so this is enforced as
c     well. Proof: tau =1/2 \sum_i (\nabla\phi_i)\cdot(\nabla\phi_i)
c     hence tau.eq.0 iff (\nabla\phi_i)=0 for all i. Delrho is 
c     defined as delrho = \sum_i\nabla\phi_i^2 hence if the gradient
c     of the orbital is zero for all orbitals i then delrho must be 
c     zero as well. QED.
c
      if (kske) then
        if (what.eq.0) then
          do ii = 1, ipol
            do n = 1, nq
              if (ttau(n,ii,1).le.0.0d0) then
                ttau(n,ii,1) = 0.0d0
                if (grad) then
                  delrho(n,1,ii,1) = 0.0d0
                  delrho(n,2,ii,1) = 0.0d0
                  delrho(n,3,ii,1) = 0.0d0
                endif
              endif
            enddo
          enddo
        else if (what.eq.1) then
          do ii = 1, ipol
            do n = 1, nq
              if (ttau(n,ii,npert+1).le.0.0d0) then
                ttau(n,ii,npert+1) = 0.0d0
                if (grad) then
                  delrho(n,1,ii,npert+1) = 0.0d0
                  delrho(n,2,ii,npert+1) = 0.0d0
                  delrho(n,3,ii,npert+1) = 0.0d0
                endif
              endif
            enddo
          enddo
        endif
      endif
      endif
c
      call starttimer(monitor_rhocomp2)
      if(what.eq.0) then
c
c       Only construct total densities for regular case
c
         if (ipol.eq.2)then
            call dcopy(nq, rho(1,2,1), 1, rho(1,1,1), 1)
            call daxpy(nq, 1.d0, rho(1,3,1), 1, rho(1,1,1), 1)
         endif
      endif
      if(what.eq.0) then
         if(xcreplicated.and.dorepdm) then
            g_dens(1)=g_keepd(1)
            if(ipol.eq.2)  g_dens(2)=g_keepd(2)
         endif
      endif
      call endtimer(monitor_rhocomp2)
c     
 1688 continue
      if(dftnbget) then
         if (.not.ma_pop_stack(l_nz))
     &        call errquit('xcrho:pop_stack failed', 13, MA_ERR)
      endif

      call endtimer(monitor_xcrho)
      return
      end
c
      subroutine xc_drhonuc(nq,ipol,n0,n1,
     D     delrho,heschi,delchi,ff,ffd,
     I     ii,iat)
      implicit none
      integer nq,ipol,ii,n0,n1,iat
      double precision delrho(nq,3,ipol,3,*)
      double precision heschi(nq,6,*)
      double precision delchi(nq,3,*)
      double precision ff(nq,*)
      double precision ffd(nq,3,*)
c
      integer n,mu,mu1
c
      do mu=n0,n1
         mu1=mu-n0+1
         do n = 1, nq
            delrho(n,1,ii,1,iat) = delrho(n,1,ii,1,iat) -
     &           heschi(n,1,mu)*ff(n,mu1) - 
     -           delchi(n,1,mu)*ffd(n,1,mu1)
            delrho(n,2,ii,1,iat) = delrho(n,2,ii,1,iat) -
     &           heschi(n,2,mu)*ff(n,mu1) - 
     -           delchi(n,1,mu)*ffd(n,2,mu1)
            delrho(n,3,ii,1,iat) = delrho(n,3,ii,1,iat) -
     &           heschi(n,3,mu)*ff(n,mu1) - 
     -           delchi(n,1,mu)*ffd(n,3,mu1)
            delrho(n,1,ii,2,iat) = delrho(n,1,ii,2,iat) -
     &           heschi(n,2,mu)*ff(n,mu1) - 
     -           delchi(n,2,mu)*ffd(n,1,mu1)
            delrho(n,2,ii,2,iat) = delrho(n,2,ii,2,iat) -
     &           heschi(n,4,mu)*ff(n,mu1) - 
     -           delchi(n,2,mu)*ffd(n,2,mu1)
            delrho(n,3,ii,2,iat) = delrho(n,3,ii,2,iat) -
     &           heschi(n,5,mu)*ff(n,mu1) - 
     -           delchi(n,2,mu)*ffd(n,3,mu1)
            delrho(n,1,ii,3,iat) = delrho(n,1,ii,3,iat) -
     &           heschi(n,3,mu)*ff(n,mu1) - 
     -           delchi(n,3,mu)*ffd(n,1,mu1)
            delrho(n,2,ii,3,iat) = delrho(n,2,ii,3,iat) -
     &           heschi(n,5,mu)*ff(n,mu1) - 
     -           delchi(n,3,mu)*ffd(n,2,mu1)
            delrho(n,3,ii,3,iat) = delrho(n,3,ii,3,iat) -
     &           heschi(n,6,mu)*ff(n,mu1) - 
     -           delchi(n,3,mu)*ffd(n,3,mu1)
         enddo
      enddo
      return
      end
      subroutine xc_dchiff(nq,n0,n1,
     P     delrho,delchi,ff)
      implicit none
      integer nq,n0,n1
      double precision delrho(nq,3)
      double precision delchi(nq,3,*)
      double precision ff(nq,*)
c
      integer n,mu,mu1
c
      do mu=n0,n1
         mu1=mu-n0+1
         do n = 1, nq
            delrho(n,1) = delrho(n,1) + delchi(n,1,mu)*ff(n,mu1)
            delrho(n,2) = delrho(n,2) + delchi(n,2,mu)*ff(n,mu1)
            delrho(n,3) = delrho(n,3) + delchi(n,3,mu)*ff(n,mu1)
         enddo
      enddo
      return
      end
      subroutine xc_dchiffp(nq,ipol2,n0,n1,
     P     rho,delchi,ff,
     M     ii,iat)
      implicit none
      integer nq,ipol2,n0,n1
      double precision rho(nq,ipol2,3,*)
      double precision delchi(nq,3,*)
      double precision ff(nq,*)
      integer ii,iat
c
      integer n,mu,mu1
c
      do mu=n0,n1
         mu1=mu-n0+1
         do n = 1, nq
            rho(n,ii,1,iat) = rho(n,ii,1,iat)-delchi(n,1,mu)*ff(n,mu1)
            rho(n,ii,2,iat) = rho(n,ii,2,iat)-delchi(n,2,mu)*ff(n,mu1)
            rho(n,ii,3,iat) = rho(n,ii,3,iat)-delchi(n,3,mu)*ff(n,mu1)
         enddo
      enddo
      return
      end

c@@@@@@@@@@@@@@@@@
      integer function xc_rhoscreen(grad,ipol,natoms,npert,
     I     iniz,
     W     tol_rho,wmax,
     O     nz,
     R     rchi_atom,rdelchi_atom,rdens_atom)
      implicit none
#include "global.fh"
      logical grad
      integer ipol,natoms,npert
      integer iniz(natoms) ! mapping of nbf_ao -> mbf
      double precision tol_rho,wmax
      double precision rchi_atom(*)
      double precision rdelchi_atom(*)
      double precision rdens_atom(natoms,natoms,ipol)
c
c      output
c
      integer nz(*)
c
      integer iat,jat,inizia,inizja,iat_in
      integer ipert
      double precision FUNC_MAXI,FUNC_MAXJ,DELFUNC_MAX,FUNC_MAX,
     ,     P_MAXJ,P_MAXJ_A,P_MAXJ_B,P_MAXIJ
      integer nonzero,ii
      double precision dabsmax
      external dabsmax
c
      FUNC_MAX = dabsmax(natoms,rchi_atom)
      DELFUNC_MAX=0d0
      if(grad) DELFUNC_MAX = dabsmax(natoms,rdelchi_atom)
      
      nonzero=0
cscat
#ifdef STAG
      do iat_in = 1+ga_nodeid(), natoms+ga_nodeid()
         iat=mod(iat_in,natoms)
         if(iat.eq.0) iat=natoms
#else
      do iat = 1, natoms
#endif
         inizia = iniz(iat)
         if (inizia.ne.0) then
c     
c     screening parameters
c     
            FUNC_MAXI = rchi_atom(iat)
            if(grad)
     .           FUNC_MAXI = max(FUNC_MAXI,rdelchi_atom(iat))
            FUNC_MAXJ = max(FUNC_MAX,DELFUNC_MAX)
            if (ipol.gt.1)then
               P_MAXJ_A = dabsmax(natoms,rdens_atom(1,iat,1))
               P_MAXJ_B = dabsmax(natoms,rdens_atom(1,iat,2))
               P_MAXJ = MAX(P_MAXJ_A, P_MAXJ_B)
            else
               P_MAXJ=0d0
               do jat=1,iat
                  if(iniz(jat).ne.0) 
     .                 P_MAXJ=max(P_MAXJ,rdens_atom(jat,iat,1))
               enddo
            endif
            if (wmax*FUNC_MAXI*FUNC_MAXJ*P_MAXJ.ge.tol_rho) then
               do jat = 1, iat
                  inizja = iniz(jat)
                  if (inizja.ne.0) then
c     
c     screening parameters
c     
                     FUNC_MAXJ=rchi_atom(jat)
                     if(grad)
     .                    FUNC_MAXJ = max(FUNC_MAXJ,rdelchi_atom(jat))
                     P_MAXIJ = rdens_atom(iat,jat,1)
                     if(ipol.eq.2) P_MAXIJ = max(P_MAXIJ,
     &                    rdens_atom(iat,jat,2))
                     if(wmax*FUNC_MAXI*FUNC_MAXJ*P_MAXIJ.ge.tol_rho)then
c     
                        do ipert = 1,npert
                           do ii = 1, ipol
c     
c     screening parameters
c     
                              P_MAXIJ = rdens_atom(iat,jat,ii)
                              if (wmax*FUNC_MAXI*FUNC_MAXJ*P_MAXIJ.
     P                             ge.tol_rho) then
                                 nonzero=nonzero+1
                                 nz(nonzero)=(iat-1)*natoms+jat
                              endif
                           enddo
                        enddo
                     endif
                  endif
               enddo            ! loop jat
            endif
         endif
      enddo                     ! loop iat
      xc_rhoscreen=nonzero
      return
      end
      subroutine xc_getdmblock(nz,nonz0,natoms,cetobfr,g_dens,
     A     Pmat,nbhandl)
      implicit none
#include "dist.fh"
#include "dft_fdist.fh"
      integer nz(*)
      integer natoms
      integer nbhandl
      integer g_dens
      integer cetobfr(2,*)
      integer nonz0
      double precision pmat(*)
c
      integer ij0
      integer iat0,jat0,ifirst,ilast,jfirst,jlast,nbfia
c
      ij0=nz(nonz0)
      iat0=(ij0-1)/natoms+1
      jat0=ij0-(iat0-1)*natoms
      ifirst = cetobfr(1,iat0)
      ilast = cetobfr(2,iat0)
      jfirst = cetobfr(1,jat0)
      jlast = cetobfr(2,jat0)
      nbfia= ilast-ifirst+1
      call starttimer(monitor_ganbget)
      call ga_nbget(g_dens, 
     %     ifirst, ilast, jfirst, jlast, Pmat,nbfia,nbhandl)
      call endtimer(monitor_ganbget)
      return
      end
      subroutine xc_dmget(repdm, nbf_ld,
     %                 ilo, ihi, jlo, jhi, Pmat,nbfia)
      implicit none
      integer  ilo, ihi, jlo, jhi, nbfia,nbf_ld
      double precision pmat(nbfia,*)
      double precision repdm(*)
      integer ij,nnn
      integer i,j
c
      nnn=nbfia
      if(ilo.ne.jlo) then
         do j=jlo,jhi
            ij=((j-1)*(2*(nbf_ld+1)-j)+1)/2+ilo-j+1
            call dcopy(nnn,repdm(ij),1,pmat(1,j-jlo+1),1)
         enddo
      else
c diag block: copy only lower tr
         do j=jlo,jhi
            ij=((j-1)*(2*(nbf_ld+1)-j)+1)/2+1
            call dcopy(nnn,repdm(ij),1,pmat(j-jlo+1,j-jlo+1),1)
            nnn=nnn-1
         enddo
c     copy offdiag terms (aka transp)
         do j=1,nbfia
            do i=j+1,nbfia
               pmat(j,i)=pmat(i,j)
            enddo
         enddo
      endif
      return
      end
c
c    Essentially a clone of the old xc_rhogen.F without the new negative density screening
c    We need this for transition densities which can be -ve
c
      subroutine td_rhogen(what,
     T     tol_rho, basis, g_dens, max_at_bf, 
     N     natoms, curatoms, ncuratoms, npert,       
     I     ipol, nq, nbf, mbf, GRAD, ipol2,
     &     F, Pmat, ff, ffd,
     C     chi, delchi, heschi,
     I     ibf, iniz, ifin, 
     &     rho, delrho, lap, rchi_atom,
     &     rdelchi_atom, rdens_atom, cetobfr, wmax, 
     &     ttau, kske, dolap)
c
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "dftpara.fh"
#include "dist.fh"
#include "dft_fdist.fh"

      Logical GRAD !< [Input] .true. when using gradient corrected
                   !< functional
      Logical kske !< [Input] .true. when using kinetic energy density
                   !< functional
      integer what !< [Input] What to calculate
                   !< - what=0: calculate density
                   !< - what=1: calculate perturbed density (derivative
                   !<   with respect to electron position)
                   !< - what=2: calculate derivative with respect to
                   !<   nuclear coordinates
       logical dolap        ! true if lap is required [in]

c what=0 dens
c what=1 pert
c what=2 nucder
      integer basis        !< [Input] basis set handle
      integer ipol         !< [Input] no. of spin channels
      integer ipol2        !< [Input] no. of spin channels in density
                           !< - 1 if closed shell
                           !< - 3 if open shell
      integer nbf          !< [Input] no. of basis functions
      integer mbf          !< [Input] "restricted" no. of basis functions
      integer max_at_bf    !< [Input] max no. bf per atom
      integer nq           !< [Input] no. of quadrature points
      integer natoms       !< [Input] no. of atoms
      double precision wmax !< [Input] max weight
      integer npert        !< [Input] number of perturbed densities
      integer curatoms(*) !< [Input] indexing array for "active" atoms
      integer ncuratoms        !< [Input] number of currently active atoms
      integer g_dens(*) !< [Input] GA handle for DM
      integer ibf(mbf)     !< [Input] mapping of nbf_ao -> mbf
      integer iniz(natoms) !< [Input] mapping of nbf_ao -> mbf
      integer ifin(natoms) !< [Input] mapping of nbf_ao -> mbf
      double precision tol_rho !< [Input] accuracy for rho evaluation
      double precision chi(nq,mbf)     !< [Input] function values
      double precision delchi(nq,3,mbf)!< [Input] function gradients
      double precision heschi(nq,6,mbf)!< [Input] function hessians
c
      double precision ttau(nq,ipol,*) !< [Output] Total Kohn-Sham K.E.density
c
       double precision lap(nq,ipol2,*)
c
      double precision delrho(nq,3,ipol,*) !< [Output] Derivative of density
      double precision Pmat(*) !< [Scratch] scratch vector
      double precision F(max_at_bf*max_at_bf) !< [Scratch] scratch vector
      double precision ff(nq,*)   !< [Scratch] scratch array
      double precision ffd(nq,3,*) !< [Scratch] scratch array
      double precision rho(nq,ipol2,*) !< [Output] The density
      double precision rchi_atom(natoms) !< [Input] Screening parameters
      double precision rdelchi_atom(natoms) !< [Input] Screening parameters
      double precision rdens_atom(natoms,natoms,ipol) !< [Input] Screening parameters
      integer cetobfr(2,natoms) !< [Input] Centers to basis functions
c
c     local declarations
c
      integer i0, ii, mu, n, npol
      integer ipert        ! perturbation loop index
      integer iat, inizia, ifirst, ilast, nbfia, nnia
      integer ifinia, ifinja
      integer jat, inizja, jfirst, jlast, nbfja, nnja
      integer iatcur, jatcur
      double precision FUNC_MAX, DELFUNC_MAX, FUNC_MAXI, FUNC_MAXJ
      double precision P_MAXJ, P_MAXJ_A, P_MAXJ_B, P_MAXIJ
      double precision dabsmax
      integer g_keepd(2)
      integer nbhandl
      integer jj
      logical doffd,doitt
      external dabsmax
      external xc_rhoscreen
      integer xc_rhoscreen
      integer nonzero,nonz0
      integer i_nz,l_nz
      integer sizeblk, gindx
#ifdef DEBUG
      integer ga_nodeid
      external ga_nodeid
#endif
      call starttimer(monitor_xcrho)
      g_keepd(1)=0
      g_keepd(2)=0
c
c     Evaluate the charge density and its gradient at each of the
c     sampling points 
c
      doffd=(what.eq.2.and.grad).or.(what.eq.0.and.kske).or.
     & (what.eq.1.and.kske).or.(what.eq.0.and.dolap)
      
      npol = (ipol*(ipol+1))/2
c     to keep compilers quiet
      iatcur=1
      jatcur=1
c
      if(what.eq.0) then
         call dcopy(nq*npol,0.D0,0,rho,1)
         if (grad) call dcopy(3*nq*ipol,0.D0,0,delrho,1)
         if (kske) call dcopy(nq*ipol,0.D0,0,ttau,1)  ! total
         if (dolap) call dfill(nq*ipol2,0.D0,lap,1)  ! total
c
c     repl stuff
c
         g_keepd(1)=g_dens(1)
         if(ipol.eq.2) g_keepd(2)=g_dens(2)
         if(xcreplicated.and.dorepdm) then
            g_dens(1)=g_repdm(1)
            if(ipol.eq.2) g_dens(2)=g_repdm(2)
         endif
      elseif(what.eq.1) then
         call dfill(nq*ipol*npert,0.D0,rho,1)
         if (grad) call dfill(3*nq*ipol*npert,0.D0,delrho,1)
         if (kske) call dfill(nq*ipol*npert,0.D0,ttau,1)   ! total
      elseif(what.eq.2) then
         call dfill(nq*ipol*3*ncuratoms,0.D0,rho,1)
         if (grad)call dfill(3*nq*ipol*3*ncuratoms,0.D0,delrho,1)
      else
         call errquit('wrong what value',0,0)
      endif
c
c     Screening is accomplished by:  p(r) <= |Xi(r)|*|Xj(r)|*|Dij|
c     Xi(r) is screened on desired accuracy/max(|Xj(r)|)*max(|Dij|)
c     Dij is screened on desired accuracy/max(|Xi(r)|)*max(|Xj(r)|)
c
      i0=ipol-1
c
      FUNC_MAX = dabsmax(natoms,rchi_atom)
      DELFUNC_MAX=0d0
      if(grad) DELFUNC_MAX = dabsmax(natoms,rdelchi_atom)
c
      nonzero=0
      if(dftnbget) then
         if (.not.ma_push_get
     &        (mt_int,(ipol*npert*natoms*(natoms+1))/2,
     N        'nzmap map',l_nz,i_nz))
     &        call errquit('xcrho:push_get failed', 13, MA_ERR)
         nonzero=xc_rhoscreen(grad,ipol,natoms,npert,
     I        iniz,
     W        tol_rho,wmax,
     O        int_mb(i_nz),
     R        rchi_atom,rdelchi_atom,rdens_atom)
         
         if(nonzero.eq.0) goto 1688
c
c     prefetch first DM block
c
         nonz0=1
         call xc_getdmblock(int_mb(i_nz),nonz0,natoms,cetobfr,
     G        g_dens(1),
     A        Pmat,nbhandl)
         endif

      do 230 iat = 1, natoms
         inizia = iniz(iat)
         if (inizia.eq.0)goto 230
         if(what.eq.2) then
            iatcur = curatoms(iat)
         endif
         ifinia = ifin(iat)
         ifirst = cetobfr(1,iat)
         ilast = cetobfr(2,iat)
         nbfia = ilast-ifirst+1
         nnia = ifinia-inizia+1
c
c        screening parameters
c
         FUNC_MAXI = rchi_atom(iat)
         if(grad)
     .   FUNC_MAXI = max(FUNC_MAXI,rdelchi_atom(iat))
         FUNC_MAXJ = max(FUNC_MAX,DELFUNC_MAX)
         if(what.lt.2) then
            if (ipol.gt.1)then
            P_MAXJ_A = dabsmax(natoms,rdens_atom(1,iat,1))
            P_MAXJ_B = dabsmax(natoms,rdens_atom(1,iat,2))
            P_MAXJ = MAX(P_MAXJ_A, P_MAXJ_B)
         else
            P_MAXJ=0d0
            do jat=1,iat
               if(iniz(jat).ne.0) 
     .              P_MAXJ=max(P_MAXJ,rdens_atom(jat,iat,1))
            enddo
         endif
         if (wmax*FUNC_MAXI*FUNC_MAXJ*P_MAXJ.lt.tol_rho) goto 225
      endif
         do 220 jat = 1, iat
            inizja = iniz(jat)
            if (inizja.eq.0)goto 220
            if(what.eq.2) then
               jatcur = curatoms(jat)
               if (jatcur.eq.0.and.iatcur.eq.0) goto 220
            endif
            call starttimer(monitor_rscreen0)
            ifinja = ifin(jat)
            jfirst = cetobfr(1,jat)
            jlast = cetobfr(2,jat)
            nbfja = jlast-jfirst+1
            nnja = ifinja-inizja+1
c
c           screening parameters
c
            FUNC_MAXJ=rchi_atom(jat)
            if(grad)
     .      FUNC_MAXJ = max(FUNC_MAXJ,rdelchi_atom(jat))
            doitt=.true.
            if(what.lt.2) then
               P_MAXIJ = rdens_atom(iat,jat,1)
               if(ipol.eq.2) P_MAXIJ = max(P_MAXIJ,
     &              rdens_atom(iat,jat,2))
               doitt=(wmax*FUNC_MAXI*FUNC_MAXJ*P_MAXIJ.ge.tol_rho)
            endif
            call endtimer(monitor_rscreen0)
            if (doitt) then
            do ii = 1, ipol
c
c              screening parameters
c
            if((what.gt.1).or.(rdens_atom(iat,jat,ii)*
     R              wmax*FUNC_MAXI*FUNC_MAXJ.ge.tol_rho)) then
c
c           Loop over perturbations
c
            do 215 ipert = 1,npert
               sizeblk=nbfia*nbfja
               call updist(monitor_size_ga_get, sizeblk)
c
               if(dftnbget) then
                  call starttimer(monitor_wait3)
                  call ga_nbwait(nbhandl)
                  call endtimer(monitor_wait3)
                  call gat_mat(F, Pmat, nbfia, nbfja, nnia, nnja, 
     I                 ifirst, jfirst, ibf(inizia), ibf(inizja))
                  nonz0=nonz0+1
                  if((npert*ipol).gt.1) then
                     gindx=ipert+(ii-1)*npert+1
                    
                     if(gindx.gt.npert*ipol)gindx=mod(gindx,npert*ipol)
                  else
                     gindx=1
                  endif
                  if(nonz0.le.nonzero)
     C                 call xc_getdmblock(int_mb(i_nz),nonz0,natoms,
     C                 cetobfr,g_dens(gindx),
     A                    Pmat,nbhandl)
               if (wmax*FUNC_MAXI*FUNC_MAXJ*
     .              dabsmax(sizeblk,Pmat).lt.tol_rho)goto 215
               else
                  call starttimer(monitor_gaget)
                  if(truerepdm) then
                  call xc_dmget(dbl_mb(k_repdm(ipert+(ii-1)*npert)), 
     &                    nbf_ld,
     %                 ifirst, ilast, jfirst, jlast, Pmat,nbfia)
                  else
                  call ga_get(g_dens(ipert+(ii-1)*npert), 
     %                 ifirst, ilast, jfirst, jlast, Pmat,nbfia)
                  endif
                  call endtimer(monitor_gaget)
               if (wmax*FUNC_MAXI*FUNC_MAXJ*
     .              dabsmax(sizeblk,Pmat).lt.tol_rho)goto 215
               call gat_mat(F, Pmat, nbfia, nbfja, nnia, nnja, ifirst, 
     &                      jfirst, ibf(inizia), ibf(inizja))
               endif

c
               call starttimer(monitor_rhocomp)
               if(iat.ne.jat) call dscal(nnia*nnja,2d0,F,1)
c
c              Compute Xiat(r)*Xjat(r)*Diat,jat
c
               call dgemm('n','n',nq,nnja,nnia,1d0,
     A              chi(1,inizia),nq,F,nnia,0d0,ff,nq)
               if(what.lt.2) then
                  jj=i0+ii
                  if(what.eq.1) jj=ii
                  do mu=inizja,ifinja
                     do n=1,nq
                        rho(n,jj,ipert) = rho(n,jj,ipert) + 
     P                       chi(n,mu)*ff(n,mu-inizja+1)
                     enddo
                  enddo
               endif
               if(doffd) then
                  call dgemm('n','n',nq*3,nnja,nnia,1d0,
     A                 delchi(1,1,inizia),nq*3,F,nnia,0d0,ffd,nq*3)
c
c              build tau for meta GGA
c
                  if(kske) then
                     do mu=inizja,ifinja
                        do n=1,nq
                           ttau(n,ii,ipert) = ttau(n,ii,ipert)+0.5d0*(
     &                       delchi(n,1,mu)*ffd(n,1,mu-inizja+1)+
     &                       delchi(n,2,mu)*ffd(n,2,mu-inizja+1)+
     &                       delchi(n,3,mu)*ffd(n,3,mu-inizja+1))
                        enddo
                     enddo
                  endif

                  if (dolap) then
                     do mu=inizja,ifinja
                        do n=1,nq
                           ! total
                           lap(n,ii,ipert) = lap(n,ii,ipert) + 2d0*
     &                          (delchi(n,1,mu)*ffd(n,1,mu-inizja+1) +
     &                          delchi(n,2,mu)*ffd(n,2,mu-inizja+1) +
     &                          delchi(n,3,mu)*ffd(n,3,mu-inizja+1) +
     &                          heschi(n,1,mu)*ff(n,mu-inizja+1) +
     &                          heschi(n,4,mu)*ff(n,mu-inizja+1) +
     &                          heschi(n,6,mu)*ff(n,mu-inizja+1))
                        enddo
                     enddo
                  endif
               endif
c
               if((what.eq.2.and.jatcur.ne.0).or.
c
c     We need the "grad" code at zero order in the nuclear deriv case,
c     but we can skip this part if iat is not active
c     
     O              (what.ne.2.and.grad)) then
                  if(what.eq.0.or.what.eq.1) then
                     call xc_dchiff(nq,inizja,ifinja,
     P                    delrho(1,1,ii,ipert),delchi,
     F                    ff)
                  elseif(what.eq.2) then
                     call xc_dchiffp(nq,ipol2,inizja,ifinja,
     P                    rho,delchi,ff,
     M                    ii,jat)
                     if (grad) then
c     
c     Compute nuclear gradient of delrho
c     
                        call xc_drhonuc(nq,ipol,inizja,ifinja,
     D                       delrho,heschi,delchi,ff,ffd,
     I                       ii,jat)
                     endif
                  endif
               endif
c
               if((what.eq.2.and.iatcur.ne.0).or.
     O              (what.ne.2.and.grad)) then
c
c                 Compute delXiat(r)*Xjat(r)*Diat,jat
c                  
               call dgemm('n','t',nq,nnia,nnja,1d0,
     A              chi(1,inizja),nq,F,nnia,0d0,ff,nq)
                  if(what.lt.2) then
                     call xc_dchiff(nq,inizia,ifinia,
     P                    delrho(1,1,ii,ipert),delchi,
     F                    ff)
                  elseif(what.eq.2) then
                     call xc_dchiffp(nq,ipol2,inizia,ifinia,
     P                    rho,delchi,ff,
     M                    ii,iat)
                     if(grad) then
                        call dgemm('n','t',nq*3,nnia,nnja,1d0,
     A                    delchi(1,1,inizja),nq*3,F,nnia,0d0,ffd,nq*3)
c     
c     Compute nuclear gradient of delrho
c     
                        call xc_drhonuc(nq,ipol,inizia,ifinia,
     D                       delrho,heschi,delchi,ff,ffd,
     I                       ii,iat)
                        
                     endif
                  endif
            endif
               call endtimer(monitor_rhocomp)
  215       continue
         endif
         enddo
         endif
  220    continue
  225    continue
  230 continue

      call starttimer(monitor_rhocomp2)
      if(what.eq.0) then
c
c     Only construct total densities for regular case
c
         if (ipol.eq.2)then
            call dcopy(nq, rho(1,2,1), 1, rho(1,1,1), 1)
            call daxpy(nq, 1.d0, rho(1,3,1), 1, rho(1,1,1), 1)
         endif
      endif
      if(what.eq.0) then
         if(xcreplicated.and.dorepdm) then
            g_dens(1)=g_keepd(1)
            if(ipol.eq.2)  g_dens(2)=g_keepd(2)
         endif
      endif
      call endtimer(monitor_rhocomp2)
c     
 1688 continue
      if(dftnbget) then
         if (.not.ma_pop_stack(l_nz))
     &        call errquit('xcrho:pop_stack failed', 13, MA_ERR)
      endif

      call endtimer(monitor_xcrho)

      return
      end

      subroutine xc_delrhodotr(nq,ipol,qxyz,rho,delrho,delrhoq)
      implicit none
      integer nq,ipol                 ![in]
      double precision rho(nq,*)      ! [in]
      double precision delrho(nq,3,*) ! [in]
      double precision qxyz(3,*)      ! [in]
      double precision delrhoq(nq,*)  ! [out]
c
      integer ii,n

      if(ipol.eq.1) then
         do n = 1, nq
            delrhoq(n,1) = 3d0*rho(n,1)+
     D           delrho(n,1,1)*qxyz(1,n) +
     D           delrho(n,2,1)*qxyz(2,n) +
     D           delrho(n,3,1)*qxyz(3,n) 
         enddo
      else
      do ii=2,3
         do n = 1, nq
            delrhoq(n,ii-1) = 3d0*rho(n,ii)+
     D           delrho(n,1,ii-1)*qxyz(1,n) +
     D           delrho(n,2,ii-1)*qxyz(2,n) +
     D           delrho(n,3,ii-1)*qxyz(3,n) 
         enddo
      enddo
      endif
      return
      end
C>
C> @}




