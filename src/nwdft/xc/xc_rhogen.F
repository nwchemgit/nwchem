      Subroutine xc_rhogen(tol_rho, basis, g_dens, max_at_bf, natoms, 
     &                     ipol, nq, nbf, mbf, ldelley, GRAD, 
     &                     F, Pmat, ff, chi, delchi, ibf, iniz, ifin, 
     &                     rho, delrho, rhoat)
c
c$Id: xc_rhogen.F,v 1.10 1997-03-20 21:35:11 d3h449 Exp $
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
      double precision delchi(3,nq,mbf)! function gradients [input]

      double precision rhoat(nq,natoms) ! atomic density values [output]
      double precision delrho(ipol,3,nq) ! density gradients [output]

      double precision Pmat(max_at_bf*max_at_bf) ! scratch vector
      double precision F(max_at_bf*max_at_bf)    ! scratch vector
      double precision ff(nq)   ! scratch array

      double precision rho(ipol*(ipol+1)/2,nq) ! density values [output]
C
C     local declarations
C
      double precision rho_tol0,rho_tol1
      integer i0,ii,mu,nu,n
      integer iat,inizia,ifirst,ilast,nbfia,nnia,idim
      integer ifinia,icount
      integer jat,inizja,jfirst,jlast,nbfja,nnja,jdim
      double precision FUNC_MAX,FUNC_MAXI,FUNC_MAX2,P_MAX
      double precision aaa
c      double precision ta
      double precision dabsmax
      external dabsmax

c
c Evaluate the charge density and its gradient at each of the sampling 
c points 
c
c
c     tolerance setting
c
      rho_tol0=tol_rho*.1d0
      rho_tol1=0.5D0*rho_tol0

      i0=ipol-1

c     Loop over products of AO basis functions to define the matrix bands.

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
            FUNC_MAX2 = FUNC_MAXI*FUNC_MAX
c
            do jat = 1, iat
               if (iniz(jat).ne.0)then
                  inizja = iniz(jat)
                  if (.not.bas_ce2bfr(basis,jat,jfirst,jlast))
     &               call errquit('Exiting in rho_gen',14)
                  do 215 ii = 1,ipol

                     call get_atom_block(g_dens(ii), basis,
     &                                   iat, jat, Pmat, idim, jdim)
                
                     P_MAX = dabsmax(idim*jdim,Pmat)
                     if ((P_MAX*FUNC_MAX2).lt.rho_tol0) goto 215

                     nbfja = jlast-jfirst+1
                     nnja = ifin(jat)-inizja+1
c
                     call gat_mat(F,Pmat,nbfia,nbfja,nnia,nnja,
     &                            ifirst,jfirst,ibf(inizia),ibf(inizja))
c
                     if (iat.ne.jat) call dscal(nq,2.d0,ff,1)
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
                           rho(i0+ii,n) = rho(i0+ii,n) + chi(n,mu)*ff(n)
                        enddo
                        if (ldelley)then
                           if (iat.ne.jat)then
                              do n = 1, nq
                                 aaa = chi(n,mu)*ff(n)
                                 rhoat(n,iat) = rhoat(n,iat) + aaa
                                 rhoat(n,jat) = rhoat(n,jat) + aaa
                              enddo
                           else
                              do n = 1, nq
                                 aaa = chi(n,mu)*ff(n)
                                 rhoat(n,iat) = rhoat(n,iat) + aaa
                              enddo
                           endif
                        endif
                        if (grad)then
                           do n = 1, nq
                              delrho(ii,1,n) = delrho(ii,1,n) + 
     &                                         delchi(1,n,mu)*ff(n)
                              delrho(ii,2,n) = delrho(ii,2,n) + 
     &                                         delchi(2,n,mu)*ff(n)
                              delrho(ii,3,n) = delrho(ii,3,n) + 
     &                                         delchi(3,n,mu)*ff(n)
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
                              delrho(ii,1,n) = delrho(ii,1,n) + 
     &                                         delchi(1,n,nu)*ff(n)
                              delrho(ii,2,n) = delrho(ii,2,n) + 
     &                                         delchi(2,n,nu)*ff(n)
                              delrho(ii,3,n) = delrho(ii,3,n) + 
     &                                         delchi(3,n,nu)*ff(n)
                           enddo
                        enddo
                     endif
  215             continue
               endif
            enddo
         endif
  230 continue
      if (ipol.eq.2)then
         call cut_tol(nq,rho(3,1),3,rho_tol1)
         call cut_tol(nq,rho(2,1),3,rho_tol1)
         do 240 n = 1, nq
            rho(1,n) = rho(2,n) + rho(3,n)
  240    continue
      else
         call cut_tol(nq,rho,1,rho_tol0)
      endif
      if (ldelley) call cut_tol(nq*natoms,rhoat,1,rho_tol0)
      return
      end
      
      Subroutine cut_tol(n,vec,inc,tol)
      implicit none
      integer n,inc           ! [input]
      double precision vec(n) ! [input/ouput]
      double precision tol    ! [input]
c
c     local
c
      integer ii,i
c
c     zeroes array elements under a given threshold
c
      i = 1
c     dir$ UNROLL 3
      do ii = 1, n
         if (vec(i).lt.tol) vec(i) = tol
         i = i+inc
      enddo
      return
      end
