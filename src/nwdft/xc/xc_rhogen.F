      Subroutine xc_rhogen(tol_rho, basis, g_dens, max_at_bf, natoms, 
     &                     ipol, nq, nbf, mbf, ldelley, GRAD, 
     &                     F, Pmat, ff, chi, delchi, ibf, iniz, ifin, 
     &                     rho, delrho, rhoat, scr, rchi_atom,
     &                     rdelchi_atom, rdens_atom, cetobfr)
c
c$Id: xc_rhogen.F,v 1.16 1998-07-14 18:42:58 d3h449 Exp $
c
      implicit none
c
      Logical GRAD         ! true when use grad corrected funct.[input]
      logical ldelley      ! true when Delley weights are used [input]
      integer basis        ! basis handle           [input]
      integer ipol         ! no. of spin states     [input]
      integer nbf          ! no. of basis functions [input]
      integer mbf          ! "restricted" no. of basis functions [input]
      integer max_at_bf    ! max no. bf per atom [input]
      integer nq           ! no. of quadrature points [input]
      integer natoms       ! no. of atoms [input]
      integer g_dens(2)    ! GA handle for DM [input]
      integer ibf(mbf)     ! mapping of nbf_ao -> mbf
      integer iniz(natoms) ! mapping of nbf_ao -> mbf
      integer ifin(natoms) ! mapping of nbf_ao -> mbf
      double precision tol_rho ! accuracy for rho evaluation [input]
      double precision chi(nq,mbf)     ! function values [input]
      double precision delchi(nq,3,mbf)! function gradients [input]
      double precision rhoat(nq,natoms) ! atomic density values [output]
      double precision delrho(nq,3,ipol) ! density gradients [output]
      double precision Pmat(max_at_bf*max_at_bf) ! scratch vector
      double precision F(max_at_bf*max_at_bf)    ! scratch vector
      double precision ff(nq)   ! scratch array
      double precision rho(nq,(ipol*(ipol+1))/2) ! density values [output]
      double precision scr(nq)   ! scratch array
      double precision rchi_atom(natoms)
      double precision rdelchi_atom(natoms)
      double precision rdens_atom(natoms,natoms,ipol)
      integer cetobfr(2,natoms)
c
c     local declarations
c
      integer i0, ii, mu, nu, n, npol
      integer iat, inizia, ifirst, ilast, nbfia, nnia, idim
      integer ifinia, ifinja, icount
      integer jat, inizja, jfirst, jlast, nbfja, nnja, jdim
      double precision FUNC_MAX, DELFUNC_MAX, FUNC_MAXI, FUNC_MAXJ
      double precision P_MAXJ, P_MAXJ_A, P_MAXJ_B, P_MAXIJ
      double precision tol_rho_tmp, aaa, reset
      double precision neg_d_max, neg_d_max_der, neg_d_max_at
      integer num_neg, num_neg_der, num_neg_at
      double precision dabsmax
      external dabsmax
c
c     Evaluate the charge density and its gradient at each of the sampling 
c     points 
c
c     Initial all incoming arrays:
c
      npol = (ipol*(ipol+1))/2
      call dfill(nq*npol,0.D0,rho,1)
      if (grad)call dfill(3*nq*ipol,0.D0,delrho,1)
      if (ldelley)call dfill(nq*natoms,0.D0,rhoat,1)
c
c     Screening is accomplished by:  p(r) <= |Xi(r)|*|Xj(r)|*|Dij|
c     Xi(r) is screened on desired accuracy/max(|Xj(r)|)*max(|Dij|)
c     Dij is screened on desired accuracy/max(|Xi(r)|)*max(|Xj(r)|)
c
      i0=ipol-1
c
      FUNC_MAX = dabsmax(natoms,rchi_atom)
      DELFUNC_MAX = dabsmax(natoms,rdelchi_atom)
c
      do 230 iat = 1, natoms
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
         if (ipol.gt.1)then
            P_MAXJ_A = dabsmax(natoms,rdens_atom(1,iat,1))
            P_MAXJ_B = dabsmax(natoms,rdens_atom(1,iat,2))
            P_MAXJ = MAX(P_MAXJ_A, P_MAXJ_B)
         else
            P_MAXJ = dabsmax(natoms,rdens_atom(1,iat,1))
         endif
         if (FUNC_MAXI*FUNC_MAXJ*P_MAXJ.lt.tol_rho) goto 225
c
         do 220 jat = 1, iat
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
            if (ipol.eq.1)then
               P_MAXIJ = rdens_atom(iat,jat,1)
            else
               P_MAXIJ = max(rdens_atom(iat,jat,1),
     &                       rdens_atom(iat,jat,2))
            endif
            if (FUNC_MAXI*FUNC_MAXJ*P_MAXIJ.lt.tol_rho) goto 215
            tol_rho_tmp = tol_rho/(FUNC_MAXI*FUNC_MAXJ)
c
            do 210 ii = 1, ipol
c
c              screening parameters
c
               P_MAXIJ = rdens_atom(iat,jat,ii)
               if (FUNC_MAXI*FUNC_MAXJ*P_MAXIJ.lt.tol_rho)goto 210
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
                  call daxpy(nq,1.d0,scr,1,rho(1,i0+ii),1)
c
c                 compute atomic densities for delley weights
c
                  if (ldelley)then
                     if (iat.ne.jat)then
                        call daxpy(nq,1.d0,scr,1,rhoat(1,iat),1)
                        call daxpy(nq,1.d0,scr,1,rhoat(1,jat),1)
                     else
                        call daxpy(nq,1.d0,scr,1,rhoat(1,iat),1)
                     endif
                  endif
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
  215       continue
  220    continue
  225    continue
  230 continue
c
c     Temp check for negative numerical densities:
c 
c     reset negative densities to 1d-30?
c
      reset = 1d-30
c
c     Check density rho
c
      call neg_check_rho(rho(1,ipol), delrho, nq, ipol, reset, 
     &                   num_neg, neg_d_max, grad)
c      if (num_neg.ne.0)
c     &   write(6,*)' negative numerical density components: ',
c     &   num_neg, neg_d_max
c
      if (ipol.eq.2)then
         call dcopy(nq, rho(1,2), 1, rho(1,1), 1)
         call daxpy(nq, 1.d0, rho(1,3), 1, rho(1,1), 1)
      endif
c
      return
      end
      subroutine neg_check_rho(rho, delrho, nq, ipol, reset, n_neg, 
     &                         max_neg, grad)
c
      implicit none
c
      integer nq, ipol, n_neg, i, n
      double precision rho(nq,ipol), delrho(nq,3,ipol), reset, max_neg
      logical grad
c
      n_neg = 0
      max_neg = 0.0d0
      do i = 1, ipol
         do n = 1, nq
            if (rho(n,i).lt.0.d0)then
               if (rho(n,i).lt.max_neg)max_neg = rho(n,i)
               rho(n,i) = reset
               if (grad)then
                  delrho(n,1,i) = reset
                  delrho(n,2,i) = reset
                  delrho(n,3,i) = reset
               endif
               n_neg = n_neg + 1
            endif
         enddo
      enddo
      return
      end
