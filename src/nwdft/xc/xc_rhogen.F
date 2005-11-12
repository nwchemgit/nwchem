* $Id: xc_rhogen.F,v 1.29 2005-11-12 00:32:47 edo Exp $
c     Modified for construction of multiple sets of perturbed densities
c     and explicit derivative densities wrt nuclei
c
c     BGJ - 8/98
      Subroutine xc_rhogen(what,
     T     tol_rho, basis, g_dens, max_at_bf, 
     N     natoms, curatoms, ncuratoms, npert,       
     I     ipol, nq, nbf, mbf, GRAD, ipol2,
     &                     F, Pmat, ff, ffd,
     C     chi, delchi, heschi,
     I     ibf, iniz, ifin, 
     &                     rho, delrho,  rchi_atom,
     &                     rdelchi_atom, rdens_atom, cetobfr,wmax, 
     &                     tau, kske)
      implicit none
#include "errquit.fh"
#include "dftpara.fh"

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
      double precision tau(nq,3,ipol)   ! KohnSham K.E.density [output]
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
      integer i0, ii, mu, nu, n, npol
      integer ipert        ! perturbation loop index
      integer iat, inizia, ifirst, ilast, nbfia, nnia, idim
      integer ifinia, ifinja, icount
      integer jat, inizja, jfirst, jlast, nbfja, nnja, jdim
      integer iatcur, jatcur
      double precision FUNC_MAX, DELFUNC_MAX, FUNC_MAXI, FUNC_MAXJ
      double precision P_MAXJ, P_MAXJ_A, P_MAXJ_B, P_MAXIJ
      double precision tol_rho_tmp, aaa
      double precision dabsmax
      integer nstride,nql,g_update(2),jjj,iii
      integer nbhandl
      logical nbfirst
      integer jj
      logical lodd,doffd,w2grad
      external dabsmax


c
c     Evaluate the charge density and its gradient at each of the sampling 
c     points 
      w2grad=what.eq.2.and.grad
      doffd=w2grad.or.
     W     (what.eq.0.and.kske)
      nbfirst=.true.
      lodd=mod(nq,2).ne.0
      
      nql=nq
      if(lodd) nql=nq-1
      npol = (ipol*(ipol+1))/2
      nstride=2
      g_update(1)=g_dens(1)	
      g_update(2)=g_dens(2)	
c     to keep compilers quiet
      iatcur=1
      jatcur=1

      if(what.eq.0) then
         call dcopy(nq*npol,0.D0,0,rho,1)
         if (grad)call dcopy(3*nq*ipol,0.D0,0,delrho,1)
         if (kske)call dfill(3*nq*ipol,0.D0,tau,1)
c
c     repl stuff
c
         if(xcreplicated.and.dorepdm) then
            g_update(1)=g_repdm(1)
            g_update(2)=g_repdm(2)
         endif
      elseif(what.eq.1) then
         call dfill(nq*ipol*npert,0.D0,rho,1)
         if (grad)call dfill(3*nq*ipol*npert,0.D0,delrho,1)
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

      do 230 iat = 1, natoms
         inizia = iniz(iat)
         if (inizia.eq.0)goto 230
         if(what.eq.2) then
            iatcur = curatoms(iat)
            iii=(iatcur-1)*3
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
         if(what.eq.0.or.what.eq.1) then
         if (ipol.gt.1)then
            P_MAXJ_A = dabsmax(natoms,rdens_atom(1,iat,1))
            P_MAXJ_B = dabsmax(natoms,rdens_atom(1,iat,2))
            P_MAXJ = MAX(P_MAXJ_A, P_MAXJ_B)
         else
            P_MAXJ=0d0
            do jat=1,natoms
            if(iniz(jat).ne.0) 
     .              P_MAXJ=max(P_MAXJ,rdens_atom(jat,iat,1))
            enddo
         endif
         if (wmax*FUNC_MAXI*FUNC_MAXJ*P_MAXJ.lt.tol_rho) goto 225
      else
c
c     !!! Cutoff temporarily commented out for other cases !!!
c
      endif
         do 220 jat = 1, iat
            inizja = iniz(jat)
            if (inizja.eq.0)goto 220
            if(what.eq.2) then
               jatcur = curatoms(jat)
               jjj=(jatcur-1)*3
               if (jatcur.eq.0.and.iatcur.eq.0) goto 220
            endif
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
            if(what.eq.0.or.what.eq.1) then
               if (ipol.eq.1)then
                  P_MAXIJ = rdens_atom(iat,jat,1)
               else
                  P_MAXIJ = max(rdens_atom(iat,jat,1),
     &                 rdens_atom(iat,jat,2))
               endif
            if (wmax*FUNC_MAXI*FUNC_MAXJ*P_MAXIJ.lt.tol_rho) goto 216
         else
c
c     !!! Cutoff temporarily commented out for other cases !!!
c
         endif
            tol_rho_tmp = tol_rho/(FUNC_MAXI*FUNC_MAXJ)
c
c
c           Loop over perturbations
c
            do 215 ipert = 1,npert
            do 210 ii = 1, ipol
c
c              screening parameters
c
               if(what.eq.0.or.what.eq.1) then
                  P_MAXIJ = rdens_atom(iat,jat,ii)
                  if (wmax*FUNC_MAXI*FUNC_MAXJ*P_MAXIJ.lt.tol_rho)
     G                 goto 210
               else
c
c     !!! Cutoff temporarily commented out for other cases !!!
c
               endif

c
               if(what.ne.1) then
                  if(dftnbget) then
                     call get_atom_block_nb(g_update(ii), basis,
     &                    iat, jat, Pmat, idim, jdim, nbhandl)
                     call ga_nbwait(nbhandl)
                  else
                     call get_atom_block(g_update(ii), basis,
     &                    iat, jat, Pmat, idim, jdim)
                  endif
               else
                  call ga_get(g_dens(ipert+(ii-1)*npert), 
     %                 ifirst, ilast, jfirst, jlast, Pmat,nbfia)
               endif
c
               if (wmax*FUNC_MAXI*FUNC_MAXJ*
     .              dabsmax(nbfia*nbfja,Pmat).lt.tol_rho)goto 210
               call gat_mat(F, Pmat, nbfia, nbfja, nnia, nnja, ifirst, 
     &                      jfirst, ibf(inizia), ibf(inizja))
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
                  if(kske) then
c     build tau for meta GGA
                     do mu=inizja,ifinja
                        do n=1,nq
                           tau(n,1,ii) = tau(n,1,ii)+
     P                       0.5d0*delchi(n,1,mu)*ffd(n,1,mu-inizja+1)
                           tau(n,2,ii) = tau(n,2,ii)+
     P                       0.5d0*delchi(n,2,mu)*ffd(n,2,mu-inizja+1)
                           tau(n,3,ii) = tau(n,3,ii)+
     P                       0.5d0*delchi(n,3,mu)*ffd(n,3,mu-inizja+1)
                        enddo
                     enddo
                  endif
               endif
               if((what.eq.2.and.jatcur.ne.0).or.
c
c     We need the "grad" code at zero order in the nuclear deriv case,
c     but we can skip this part if iat is not active
c     
     O              (what.ne.2.and.grad)) then
               do mu=inizja,ifinja
                  if(what.eq.0.or.what.eq.1) then
                     call xc_dchiff(nq,
     P                    delrho(1,1,ii,ipert),delchi(1,1,mu),
     F                    ff(1,mu-inizja+1))
                  elseif(what.eq.2) then
                     call xc_dchiffp(nq,ipol2,
     P                    rho(1,ii,1+jjj),
     D                    delchi(1,1,mu),ff(1,mu-inizja+1),
     M                    ii)
                     if (grad) then
c     
c     Compute nuclear gradient of delrho
c     
                        call xc_drhonuc(nq,ipol,ii,
     D                       delrho(1,1,1,1+jjj),
     H                       heschi(1,1,mu),delchi(1,1,mu),
     F                       ff(1,mu-inizja+1),ffd(1,1,mu-inizja+1))
                     endif
                  endif
               enddo ! close of mu loop
            endif
               if((what.eq.2.and.iatcur.ne.0).or.
     O              (what.ne.2.and.grad)) then
c
c                 Compute delXiat(r)*Xjat(r)*Diat,jat
c                  
               call dgemm('n','t',nq,nnia,nnja,1d0,
     A              chi(1,inizja),nq,F,nnia,0d0,ff,nq)
                  if(w2grad) 
     C              call dgemm('n','t',nq*3,nnia,nnja,1d0,
     A                    delchi(1,1,inizja),nq*3,F,nnia,0d0,ffd,nq*3)
               do nu=inizia,ifinia
                  if(what.lt.2) then
                     call xc_dchiff(nq,
     P                    delrho(1,1,ii,ipert),delchi(1,1,nu),
     F                    ff(1,nu-inizia+1))
                  elseif(what.eq.2) then
                     call xc_dchiffp(nq,ipol2,
     P                    rho(1,ii,1+iii),
     D                    delchi(1,1,nu),ff(1,nu-inizia+1),
     M                    ii)
                     if(grad) then
c     
c     Compute nuclear gradient of delrho
c     
                        call xc_drhonuc(nq,ipol,ii,
     D                       delrho(1,1,1,1+iii),
     H                       heschi(1,1,nu),delchi(1,1,nu),
     F                       ff(1,nu-inizia+1),ffd(1,1,nu-inizia+1))
                        
                     endif
                  endif
               enddo            ! end of nu loop
            endif
  210       continue
  215       continue
 216        continue
  220    continue
  225    continue
  230 continue

      if(what.eq.0) then
c
c     Only construct total densities for regular case
c
         if (ipol.eq.2)then
            call dcopy(nq, rho(1,2,1), 1, rho(1,1,1), 1)
            call daxpy(nq, 1.d0, rho(1,3,1), 1, rho(1,1,1), 1)
         endif
      endif
c
      return
      end
c
      subroutine get_atom_block_nb(g_array, basis, iat, jat, buf,
     $     idim, jdim,nbhandle)
C$Id: xc_rhogen.F,v 1.29 2005-11-12 00:32:47 edo Exp $
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "bas.fh"
c
      integer g_array, basis, iat, jat, idim, jdim
      integer nbhandle
      double precision buf(*)
c
      integer ilo, ihi, jlo, jhi
c
c     get the block of the matrix g_array (over basis functions)
c     corresponding to the atoms iat and jat, putting the
c     result into buf(1:idim, 1:jdim) as a dense array. 
c     idim, jdim are returned
c
      if (.not. bas_ce2bfr(basis, iat, ilo, ihi))
     $     call errquit('get_atom_block: ce2bfr failed', 0, BASIS_ERR)
      if (.not. bas_ce2bfr(basis, jat, jlo, jhi))
     $     call errquit('get_atom_block: ce2bfr failed', 0, BASIS_ERR)
      idim = ihi - ilo + 1
      jdim = jhi - jlo + 1
c
      if (idim.gt.0 .and. jdim.gt.0) then
         call ga_nbget(g_array, ilo, ihi, jlo, jhi, buf, idim,
     E        nbhandle)
      endif
*      write(6,*) ' gab: ', iat, jat
*      call output(buf, 1, idim, 1, jdim, idim, jdim, 1)
c
      end
      subroutine xc_drhonuc(nq,ipol,ii,
     D     delrho,heschi,delchi,ff,ffd)
      implicit none
      integer nq,ipol,ii
      double precision delrho(nq,3,ipol,*)
      double precision heschi(nq,6)
      double precision delchi(nq,3)
      double precision ff(*)
      double precision ffd(nq,*)
c
      integer n
c
      do n = 1, nq
         delrho(n,1,ii,1) = delrho(n,1,ii,1) -
     &        heschi(n,1)*ff(n) - delchi(n,1)*ffd(n,1)
         delrho(n,2,ii,1) = delrho(n,2,ii,1) -
     &        heschi(n,2)*ff(n) - delchi(n,1)*ffd(n,2)
         delrho(n,3,ii,1) = delrho(n,3,ii,1) -
     &        heschi(n,3)*ff(n) - delchi(n,1)*ffd(n,3)
         delrho(n,1,ii,2) = delrho(n,1,ii,2) -
     &        heschi(n,2)*ff(n) - delchi(n,2)*ffd(n,1)
         delrho(n,2,ii,2) = delrho(n,2,ii,2) -
     &        heschi(n,4)*ff(n) - delchi(n,2)*ffd(n,2)
         delrho(n,3,ii,2) = delrho(n,3,ii,2) -
     &        heschi(n,5)*ff(n) - delchi(n,2)*ffd(n,3)
         delrho(n,1,ii,3) = delrho(n,1,ii,3) -
     &        heschi(n,3)*ff(n) - delchi(n,3)*ffd(n,1)
         delrho(n,2,ii,3) = delrho(n,2,ii,3) -
     &        heschi(n,5)*ff(n) - delchi(n,3)*ffd(n,2)
         delrho(n,3,ii,3) = delrho(n,3,ii,3) -
     &        heschi(n,6)*ff(n) - delchi(n,3)*ffd(n,3)
      enddo
      return
      end
      subroutine xc_dchiff(nq,
     P     delrho,delchi,ff)
      implicit none
      integer nq
      double precision delrho(nq,3)
      double precision delchi(nq,3)
      double precision ff(*)
c
      integer n
      do n = 1, nq
         delrho(n,1) = delrho(n,1) + delchi(n,1)*ff(n)
         delrho(n,2) = delrho(n,2) + delchi(n,2)*ff(n)
         delrho(n,3) = delrho(n,3) + delchi(n,3)*ff(n)
      enddo
      return
      end
      subroutine xc_dchiffp(nq,ipol2,
     P     rho,delchi,ff,
     M     ii)
      implicit none
      integer nq,ipol2
      double precision rho(nq,ipol2,*)
      double precision delchi(nq,3)
      double precision ff(nq)
      integer ii
c
      integer n
      do n = 1, nq
         rho(n,ii,1) = rho(n,ii,1) - delchi(n,1)*ff(n)
         rho(n,ii,2) = rho(n,ii,2) - delchi(n,2)*ff(n)
         rho(n,ii,3) = rho(n,ii,3) - delchi(n,3)*ff(n)
      enddo
      return
      end
