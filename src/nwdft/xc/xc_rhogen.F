      Subroutine xc_rhogen(tol_rho, basis, g_dens, max_at_bf, natoms, 
     &                     ipol, nq, nbf, mbf, ldelley, GRAD, 
     &                     F, Pmat, ff, chi, delchi, ibf, iniz, ifin, 
     &                     rho, delrho, rhoat, scr)
c
c$Id: xc_rhogen.F,v 1.12 1997-04-23 21:19:20 d3h449 Exp $
c
      implicit none
c
#include "bas.fh"
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
      integer g_dens(*)    ! GA handle for DM [input]
      integer ibf(mbf)     ! 
      integer iniz(natoms) !
      integer ifin(natoms) !
      double precision tol_rho ! accuracy for rho evaluation [input]
      double precision chi(nq,mbf)     ! function values [input]
      double precision delchi(nq,3,mbf)! function gradients [input]
      double precision rhoat(nq,natoms) ! atomic density values [output]
      double precision delrho(nq,3,ipol) ! density gradients [output]
      double precision Pmat(max_at_bf*max_at_bf) ! scratch vector
      double precision F(max_at_bf*max_at_bf)    ! scratch vector
      double precision ff(nq)   ! scratch array
      double precision rho(nq,ipol*(ipol+1)/2) ! density values [output]
      double precision scr(nq)   ! scratch array
c
c     local declarations
c
      double precision rho_tol0, rho_tol1
      integer i0, ii, mu, nu, n
      integer iat, inizia, ifirst, ilast, nbfia, nnia, idim
      integer ifinia, icount
      integer jat, inizja, jfirst, jlast, nbfja, nnja, jdim
      double precision FUNC_MAX, FUNC_MAXI, FUNC_MAXJ, P_MAX
      double precision aaa
      double precision dabsmax
      external dabsmax
c
c     Evaluate the charge density and its gradient at each of the sampling 
c     points 
c
c     tolerance setting
c
      rho_tol0=tol_rho*.1d0
      rho_tol1=0.5D0*rho_tol0
c
      i0=ipol-1
c
c     Loop over products of AO basis functions to define the matrix bands.
c
      FUNC_MAX=dabsmax(nq*mbf,chi)

      do 230 iat = 1, natoms
         inizia = iniz(iat)
         ifinia = ifin(iat)
         if (inizia.ne.0)then
            if (.not.bas_ce2bfr(basis,iat,ifirst,ilast))
     &         call errquit('Exiting in rho_gen',13)
            nbfia = ilast-ifirst+1
            nnia = ifinia-inizia+1
c
            FUNC_MAXI = dabsmax(nnia*nq,chi(1,inizia))
            if (FUNC_MAXI*FUNC_MAX.lt.rho_tol0) goto 217
c
            do jat = 1, iat
               if (iniz(jat).ne.0)then
                  inizja = iniz(jat)
                  if (.not.bas_ce2bfr(basis,jat,jfirst,jlast))
     &               call errquit('Exiting in rho_gen',14)
                  nbfja = jlast-jfirst+1
                  nnja = ifin(jat)-inizja+1
                  FUNC_MAXJ = dabsmax(nnja*nq,chi(1,inizja))
c
                  if (FUNC_MAXI*FUNC_MAXJ.lt.rho_tol0) goto 216
                  do 215 ii = 1, ipol

                     call get_atom_block(g_dens(ii), basis,
     &                                   iat, jat, Pmat, idim, jdim)
                
                     P_MAX = dabsmax(idim*jdim,Pmat)
                     if ((P_MAX*FUNC_MAXI*FUNC_MAXJ).lt.tol_rho)goto 215

                     call gat_mat(F, Pmat, nbfia, nbfja, nnia, nnja,
     &                            ifirst, jfirst, ibf(inizia), 
     &                            ibf(inizja))
c
                     icount = 0
                     do mu = inizja, ifin(jat)
                        call dfill(nq,0.d0,ff,1)
                        do nu = inizia, ifinia
                           icount = icount+1
                           aaa = F(icount)
                           if (abs(aaa).gt.rho_tol0)then
                              call daxpy(nq,aaa,chi(1,nu),1,ff,1)
                           endif
                        enddo
                        if (iat.ne.jat) call dscal(nq,2.d0,ff,1)
c                       dir$ UNROLL 3
                        do n = 1, nq
                           scr(n) = chi(n,mu)*ff(n)
                        enddo
                        call daxpy(nq,1.d0,scr,1,rho(1,i0+ii),1)
                        if (ldelley)then
                           if (iat.ne.jat)then
                              call daxpy(nq,1.d0,scr,1,rhoat(1,iat),1)
                              call daxpy(nq,1.d0,scr,1,rhoat(1,jat),1)
 
                           else
                              call daxpy(nq,1.d0,scr,1,rhoat(1,iat),1)
                           endif
                        endif
                        if (grad)then
                           do n = 1, nq
                              delrho(n,1,ii) = delrho(n,1,ii) + 
     &                                         delchi(n,1,mu)*ff(n)
                              delrho(n,2,ii) = delrho(n,2,ii) + 
     &                                         delchi(n,2,mu)*ff(n)
                              delrho(n,3,ii) = delrho(n,3,ii) + 
     &                                         delchi(n,3,mu)*ff(n)
                           enddo
                        endif
                     enddo
                     if (grad)then
c                  
c                       with GRAD 2nd contrib [P(AB)X(B)]DX(A)
c                  
                        do nu = inizia, ifinia
                           call dfill(nq,0.d0,ff,1)
                           icount = nu-inizia+1
                           do mu = inizja, ifin(jat)
                              aaa = F(icount)
                              icount = icount+nnia
                              if (abs(aaa).gt.rho_tol0)then
                                 call daxpy(nq,aaa,chi(1,mu),1,ff,1)
                              endif
                           enddo
                           if (iat.ne.jat) call dscal(nq,2.d0,ff,1)
c                          dir$ UNROLL 3
                           do n = 1, nq
                              delrho(n,1,ii) = delrho(n,1,ii) + 
     &                                         delchi(n,1,nu)*ff(n)
                              delrho(n,2,ii) = delrho(n,2,ii) + 
     &                                         delchi(n,2,nu)*ff(n)
                              delrho(n,3,ii) = delrho(n,3,ii) + 
     &                                         delchi(n,3,nu)*ff(n)
                           enddo
                        enddo
                     endif
  215             continue
  216             continue
               endif
            enddo
  217       continue
         endif
  230 continue
      if (ipol.eq.2)then
         call cut_tol(nq, rho(1,3), rho_tol1)
         call cut_tol(nq, rho(1,2), rho_tol1)
         call dcopy(nq, rho(1,2), 1, rho(1,1), 1)
         call daxpy(nq, 1.d0, rho(1,3), 1, rho(1,1), 1)
      else
         call cut_tol(nq, rho, rho_tol0)
      endif
      if (ldelley) call cut_tol(nq*natoms, rhoat, rho_tol0)
      return
      end
      
      Subroutine cut_tol(n, vec, tol)
      implicit none
      integer n               ! [input]
      double precision vec(n) ! [input/ouput]
      double precision tol    ! [input]
c
c     local
c
      integer i
c
c     zeroes array elements under a given threshold
c
      do i = 1, n
         if (vec(i).lt.tol) vec(i) = tol
      enddo
      return
      end
