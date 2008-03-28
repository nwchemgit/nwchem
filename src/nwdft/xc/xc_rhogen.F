* $Id: xc_rhogen.F,v 1.39 2008-03-28 18:47:17 edo Exp $
c     Modified for construction of multiple sets of perturbed densities
c     and explicit derivative densities wrt nuclei
c
c     BGJ - 8/98
      Subroutine xc_rhogen(what,
     T     tol_rho, basis, g_dens, max_at_bf, 
     N     natoms, curatoms, ncuratoms, npert,       
     I     ipol, nq, nbf, mbf, GRAD, ipol2,
     &     F, Pmat, ff, ffd,
     C     chi, delchi, heschi,
     I     ibf, iniz, ifin, 
     &     rho, delrho, rchi_atom,
     &     rdelchi_atom, rdens_atom, cetobfr, wmax, 
     &     tau, ttau, kske)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "dftpara.fh"
#include "dist.fh"
#include "dft_fdist.fh"

      Logical GRAD         ! true when use grad corrected funct.[input]
      Logical kske         ! true if tau(kohn-sham KE) used in functional [in]
      integer what  ![in]
c what=0 dens
c what=1 pert
c what=2 nucder

 
      integer basis        ! basis handle           [input]
      integer ipol         ! no. of spin states     [input]
      integer ipol2
      integer nbf          ! no. of basis functions [input]
      integer mbf          ! "restricted" no. of basis functions [input]
      integer max_at_bf    ! max no. bf per atom [input]
      integer nq           ! no. of quadrature points [input]
      integer natoms       ! no. of atoms [input]
      double precision wmax ! max weight [in]
      integer npert        ! number of perturbed densities [input]
      integer curatoms(*) ! indexing array for "active" atoms [input]
      integer ncuratoms        ! number of currently active atoms [input]
      integer g_dens(*) ! GA handle for DM [input]
      integer ibf(mbf)     ! mapping of nbf_ao -> mbf
      integer iniz(natoms) ! mapping of nbf_ao -> mbf
      integer ifin(natoms) ! mapping of nbf_ao -> mbf
      double precision tol_rho ! accuracy for rho evaluation [input]
      double precision chi(nq,mbf)     ! function values [input]
      double precision delchi(nq,3,mbf)! function gradients [input]
      double precision heschi(nq,6,mbf)! function hessians [input]
c
      double precision tau(nq,3,ipol,*)   ! Kohn-Sham K.E.density (components) [output]
      double precision ttau(nq,ipol,*)   !  Total Kohn-Sham K.E.density [output]
c
      double precision delrho(nq,3,ipol,*)
      double precision Pmat(*) ! scratch vector
      double precision F(max_at_bf*max_at_bf)    ! scratch vector
      double precision ff(nq,*)   ! scratch array
      double precision ffd(nq,3,*) ! scratch array
      double precision rho(nq,ipol2,*)
      double precision rchi_atom(natoms)
      double precision rdelchi_atom(natoms)
      double precision rdens_atom(natoms,natoms,ipol)
      integer cetobfr(2,natoms)
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
      integer sizeblk
#ifdef DEBUG
      integer ga_nodeid
      external ga_nodeid
#endif
      call starttimer(monitor_xcrho)
c
c     Evaluate the charge density and its gradient at each of the sampling 
c     points 
      doffd=(what.eq.2.and.grad).or.(what.eq.0.and.kske).or.
     & (what.eq.1.and.kske)
      
      npol = (ipol*(ipol+1))/2
c     to keep compilers quiet
      iatcur=1
      jatcur=1
c
      if(what.eq.0) then
         call dcopy(nq*npol,0.D0,0,rho,1)
         if (grad) call dcopy(3*nq*ipol,0.D0,0,delrho,1)
         if (kske) call dfill(3*nq*ipol,0.D0,tau,1) ! components
         if (kske) call dfill(nq*ipol,0.D0,ttau,1)  ! total
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
         if (kske) call dfill(3*nq*ipol*npert,0.D0,tau,1)  ! components
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
      if(dftnbget) then
         if(npert.gt.1.or.ipol.gt.1) call errquit(' crash ',0,0)
         if (.not.ma_push_get
     &        (mt_int,(natoms*(natoms+1))/2,'nzmap map',l_nz,i_nz))
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
c
c           Loop over perturbations
c
            do 215 ipert = 1,npert
            do 210 ii = 1, ipol
c
c              screening parameters
c
               if(what.lt.2) then
                  P_MAXIJ = rdens_atom(iat,jat,ii)
                  if (wmax*FUNC_MAXI*FUNC_MAXJ*P_MAXIJ.lt.tol_rho)
     G                 goto 210
               endif
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
                  if(nonz0.le.nonzero)
     C                 call xc_getdmblock(int_mb(i_nz),nonz0,natoms,
     C                 cetobfr,g_dens(ipert+(ii-1)*npert),
     A                    Pmat,nbhandl)
               if (wmax*FUNC_MAXI*FUNC_MAXJ*
     .              dabsmax(sizeblk,Pmat).lt.tol_rho)goto 210
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
     .              dabsmax(sizeblk,Pmat).lt.tol_rho)goto 210
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
                           ! components
                           tau(n,1,ii,ipert) = tau(n,1,ii,ipert)+
     &                       0.5d0*delchi(n,1,mu)*ffd(n,1,mu-inizja+1)
                           tau(n,2,ii,ipert) = tau(n,2,ii,ipert)+
     &                       0.5d0*delchi(n,2,mu)*ffd(n,2,mu-inizja+1)
                           tau(n,3,ii,ipert) = tau(n,3,ii,ipert)+
     &                       0.5d0*delchi(n,3,mu)*ffd(n,3,mu-inizja+1)

                           ! total
                           ttau(n,ii,ipert) = tau(n,1,ii,ipert) +
     &                                        tau(n,2,ii,ipert) +
     &                                        tau(n,3,ii,ipert)
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
  210       continue
  215       continue
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
      integer iat,jat,inizia,inizja
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

      do iat = 1, natoms
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
