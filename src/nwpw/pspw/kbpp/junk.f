*     *****************************************************
*     *                                                   *
*     *                setup                              *
*     *                                                   *
*     *****************************************************

      subroutine setup(nfft1,nfft2,nfft3,unita,unitg,G,IERR)
*
* $Id$
*
      implicit none
      integer nfft1,nfft2,nfft3
      double precision unita(3,3),unitg(3,3)
      double precision G(nfft1/2+1,nfft2,nfft3,3)





      call get_unitg(unita,unitg)

      nffth3 = nfft3/2
      nffth2 = nfft2/2
      nffth1 = nfft1/2
*     **** less confusing algorithm ****
      do k3 = -nffth3+1, nffth3
         do k2 = -nffth2+1, nffth2
            do k1 = 0,nffth1
               g1 = k1*unitg(1,1) + k2*unitg(1,2) +k3*unitg(1,3)
               g2 = k1*unitg(2,1) + k2*unitg(2,2) +k3*unitg(2,3)
               g3 = k1*unitg(3,1) + k2*unitg(3,2) +k3*unitg(3,3)
               i=k1
               j=k2
               k=k3
c              if (i .lt. 0) i = i + nfft1
               if (j .lt. 0) j = j + nfft2
               if (k .lt. 0) k = k + nfft3

               G(i+1,j+1,k+1,1) = g1
               G(i+1,j+1,k+1,2) = g2
               G(i+1,j+1,k+1,3) = g3

            end do  
         end do
      end do

      IERR=0
      RETURN
      END




      subroutine kbpp(nrmax,vp,wp,rho,f,cs,sn
     >                nfft1,nfft2,nfft3,vl,vnl,G)
      implicit none



      DIMENSION VP(NRMAX,0:3),WP(NRMAX,0:3),RHO(NRMAX)
      DIMENSION VL(nfft1/2+1,nfft2,nfft3)
      DIMENSION VNL(nfft1/2+1,nfft2,nfft3,9)
      DIMENSION VNLNRM(9)
      dimension unit(3)
      DIMENSION UNITA(3,3),UNITG(3,3)
      DIMENSION G(nfft1/2+1,nfft2,nfft3,3)
      DIMENSION F(NRMAX),CS(NRMAX),SN(NRMAX)

 
      common / HAMANN / vp,wp,rho,drho,nrho,lmax
      common / PSEUDO / vl,vnl,vnlnrm,zv,lmmax 
      common / LATTIC / g,unita,unitg,unit,omega,icube


      PI=4.0d0*DATAN(1.0d0)
      TWOPI=2.0d0*PI
      FORPI=4.0d0*PI

*     Total number of non-local pseudopotentials
      LMMAX=LMAX**2

      IF(LMMAX.GT.9) THEN
        IERR=1
        RETURN
      ENDIF
      IF(NRHO.GT.NRMAX) THEN
        IERR=3
        RETURN
      ENDIF
      IF((NRHO/2)*2.EQ.NRHO) THEN
        IERR=4
        RETURN
      ENDIF

      P0=DSQRT(FORPI)
      P1=DSQRT(3.0d0*FORPI)
      P2=DSQRT(15.0d0*FORPI)

*::::::::::::::::::  Define non-local pseudopotential  ::::::::::::::::
      DO 100 L=0,LMAX-1
        DO 100 I=1,NRHO
          VP(I,L)=VP(I,L)-VP(I,LMAX)
  100 CONTINUE

*:::::::::::::::::::::  Normarization constants  ::::::::::::::::::::::
      DO 130 L=0,LMAX-1
        DO 110 I=1,NRHO
          F(I)=VP(I,L)*WP(I,L)**2
  110   CONTINUE
        A=SIMP(NRHO,F,DRHO)
        DO 120 I=L**2+1,(L+1)**2
          VNLNRM(I)=A
  120   CONTINUE
  130 CONTINUE

*======================  Fourier transformation  ======================
      DO 700 k3=1,nfft3
      DO 700 k2=1,nfft2
      DO 700 k1=1,(nfft1/2+1)
        Q=DSQRT(G(k1,k2,k3,1)**2
     >         +G(k1,k2,k3,2)**2
     >         +G(k1,k2,k3,3)**2)

        if ((k1.eq.1).and.(k2.eq.1).and.(k3.eq.1)) go to 700

        GX=G(k1,k2,k3,1)/Q
        GY=G(k1,k2,k3,2)/Q
        GZ=G(k1,k2,k3,3)/Q
        DO 200 I=1,NRHO
          CS(I)=DCOS(Q*RHO(I))
          SN(I)=DSIN(Q*RHO(I))
  200   CONTINUE

        GO TO (600,500,400,300), LMAX+1

*::::::::::::::::::::::::::::::  d-wave  ::::::::::::::::::::::::::::::
  300   CONTINUE
        F(1)=0.0d0
        DO 310 I=2,NRHO
          A=3.0d0*(SN(I)/(Q*RHO(I))-CS(I))/(Q*RHO(I))-SN(I)
          F(I)=A*WP(I,2)*VP(I,2)
  310   CONTINUE
        D=P2*SIMP(NRHO,F,DRHO)/Q
        VNL(k1,k2,k3,5)=D*(3.0d0*GZ*GZ-1.0d0)/(2.0d0*dsqrt(3.0d0))
        VNL(k1,k2,k3,6)=D*GX*GY
        VNL(k1,k2,k3,7)=D*GY*GZ
        VNL(k1,k2,k3,8)=D*GZ*GX
        VNL(k1,k2,k3,9)=D*(GX*GX-GY*GY)/(2.0d0)

*::::::::::::::::::::::::::::::  p-wave  ::::::::::::::::::::::::::::::
  400   CONTINUE
         F(1)=0.0d0
         DO 410 I=2,NRHO
           F(I)=(SN(I)/(Q*RHO(I))-CS(I))*WP(I,1)*VP(I,1)
  410    CONTINUE
         P=P1*SIMP(NRHO,F,DRHO)/Q
         VNL(k1,k2,k3,2)=P*GX
         VNL(k1,k2,k3,3)=P*GY
         VNL(k1,k2,k3,4)=P*GZ

*::::::::::::::::::::::::::::::  s-wave  :::::::::::::::::::::::::::::::
  500   CONTINUE
        DO 510 I=1,NRHO
          F(I)=SN(I)*WP(I,0)*VP(I,0)
  510   CONTINUE
        VNL(k1,k2,k3,1)=P0*SIMP(NRHO,F,DRHO)/Q

*::::::::::::::::::::::::::::::  local  :::::::::::::::::::::::::::::::
  600   CONTINUE
        DO 610 I=1,NRHO
          F(I)=RHO(I)*VP(I,LMAX)*SN(I)
  610   CONTINUE
        VL(k1,k2,k3)=SIMP(NRHO,F,DRHO)*FORPI/Q-ZV*FORPI*CS(NRHO)/(Q*Q)

  700 CONTINUE

*:::::::::::::::::::::::::::::::  G=0  ::::::::::::::::::::::::::::::::      

      DO 800 I=1,NRHO
        F(I)=VP(I,LMAX)*RHO(I)**2
  800 CONTINUE
      VL(1,1,1)=FORPI*SIMP(NRHO,F,DRHO)+TWOPI*ZV*RHO(NRHO)**2

      DO 810 I=1,NRHO
        F(I)=RHO(I)*WP(I,0)*VP(I,0)
  810 CONTINUE
      VNL(1,1,1,1)=P0*SIMP(NRHO,F,DRHO)

      DO 820 L=2,LMMAX
        VNL(1,1,1,L)=0.0d0
  820 CONTINUE

      IERR=0
      RETURN
      END

      double precision function simp(n,y,h)
      implicit none
      integer n
      double precision y(n)
      double precision h,s
      integer n,ne,no

      ne=n/2
      no=ne+1
      S=2.0d0*dsum(no,y(1),2) + 4.0d0*dsum(ne,y(2),2)-y(1)-y(n)
      simp=s*h/3.0d0
      return
      end



      subroutine get_unitg(unita,unitg)
      implicit none

******************************************************************************
*                                                                            *
*     This routine computes primitive vectors                                *
*               in reciporocal space and the volume of primitive cell.       *
*                                                                            *
*     Inputs:                                                                *
*             unita  --- primitive vectors in coordination space             *
*                                                                            *
*     Outputs:                                                               *
*             unitg  --- primitive vectors in reciprocal space               *
*                                                                            *
*     Library:  dscal from BLAS                                              *
*                                                                            *
*     Last modification:  3/30/99  by Eric Bylaska                           *
*                                                                            *
******************************************************************************


*     ------------------
*     argument variables
*     ------------------
      double precision unita(3,3), unitg(3,3)

*     ---------------
*     local variables
*     ---------------
      double precision volume
      double precision twopi

      twopi = 8.0d0*datan(1.0d0)


*     -----------------------------------------
*     primitive vectors in the reciprocal space 
*     -----------------------------------------
      unitg(1,1) = unita(2,2)*unita(3,3) - unita(3,2)*unita(2,3)
      unitg(2,1) = unita(3,2)*unita(1,3) - unita(1,2)*unita(3,3)
      unitg(3,1) = unita(1,2)*unita(2,3) - unita(2,2)*unita(1,3)
      unitg(1,2) = unita(2,3)*unita(3,1) - unita(3,3)*unita(2,1)
      unitg(2,2) = unita(3,3)*unita(1,1) - unita(1,3)*unita(3,1)
      unitg(3,2) = unita(1,3)*unita(2,1) - unita(2,3)*unita(1,1)
      unitg(1,3) = unita(2,1)*unita(3,2) - unita(3,1)*unita(2,2)
      unitg(2,3) = unita(3,1)*unita(1,2) - unita(1,1)*unita(3,2)
      unitg(3,3) = unita(1,1)*unita(2,2) - unita(2,1)*unita(1,2)
      volume = unita(1,1)*unitg(1,1)
     >       + unita(2,1)*unitg(2,1)
     >       + unita(3,1)*unitg(3,1)
      call dscal(9,twopi/volume,unitg,1)

*     ---------------------
*     volume of a unit cell
*     ---------------------
      volume=dabs(volume)

      return
      end
