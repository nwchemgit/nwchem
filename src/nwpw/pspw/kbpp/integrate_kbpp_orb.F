*
* $Id$
*

      subroutine integrate_kbpp_orb(version,
     >                            nrho,drho,lmax,locp,
     >                            wp,rho,f,cs,sn,
     >                            nfft1,nfft2,nfft3,lmmax,
     >                            G,aorbs,
     >                            ierr)
      implicit none
      integer          version
      integer          nrho
      double precision drho
      integer          lmax
      integer          locp
      double precision wp(nrho,0:lmax)
      double precision rho(nrho)
      double precision f(nrho)
      double precision cs(nrho)
      double precision sn(nrho)

      integer nfft1,nfft2,nfft3,lmmax
      double precision G(nfft1/2+1,nfft2,nfft3,3)
      double precision aorbs(nfft1/2+1,nfft2,nfft3,lmmax)

      integer ierr

      integer np,taskid,MASTER
      parameter (MASTER=0)

*     *** local variables ****
      integer lcount,task_count,nfft3d
      integer k1,k2,k3,i,l
      double precision pi,twopi,forpi
      double precision p0,p1,p2,p3,p
      double precision gx,gy,gz,a,q,d


*     **** external functions ****
      double precision dsum,simp
      external         dsum,simp

      call Parallel_np(np)
      call Parallel_taskid(taskid)

      nfft3d = (nfft1/2+1)*nfft2*nfft3
      pi=4.0d0*datan(1.0d0)
      twopi=2.0d0*pi
      forpi=4.0d0*pi

      IF(LMMAX.GT.16) THEN
        IERR=1
        RETURN
      ENDIF
      IF((NRHO/2)*2.EQ.NRHO) THEN
        IERR=2
        RETURN
      ENDIF

      P0=DSQRT(FORPI)
      P1=DSQRT(3.0d0*FORPI)
      P2=DSQRT(15.0d0*FORPI)
      P3=DSQRT(105.0d0*FORPI)


*======================  Fourier transformation  ======================
      call dcopy(lmmax*nfft3d,0.0d0,0,aorbs,1)
      task_count = -1
      DO 700 k3=1,nfft3
      DO 700 k2=1,nfft2
      DO 700 k1=1,(nfft1/2+1)
        task_count = task_count + 1
        if (mod(task_count,np).ne.taskid) go to 700

        Q=DSQRT(G(k1,k2,k3,1)**2
     >         +G(k1,k2,k3,2)**2
     >         +G(k1,k2,k3,3)**2)

        if ((k1.eq.1).and.(k2.eq.1).and.(k3.eq.1)) go to 700

        
        GX=G(k1,k2,k3,1)/Q
        GY=G(k1,k2,k3,2)/Q
        GZ=G(k1,k2,k3,3)/Q
        DO I=1,NRHO
          CS(I)=DCOS(Q*RHO(I))
          SN(I)=DSIN(Q*RHO(I))
        END DO

        lcount = lmmax+1
        GO TO (500,400,300,200,100), LMAX+1


*::::::::::::::::::::::::::::::  g-wave  ::::::::::::::::::::::::::::::
  100   CONTINUE
*::::::::::::::::::::::::::::::  f-wave  ::::::::::::::::::::::::::::::
  200   CONTINUE
        if (locp.ne.3) then
           F(1)=0.0d0
           do I=2,NRHO
             A=SN(I)/(Q*RHO(I))
             A=15.0d0*(A-CS(I))/(Q*RHO(I))**2 - 6*A + CS(I)
             F(I)=A*WP(I,3)
           end do
           D=P3*SIMP(NRHO,F,DRHO)/Q

*          **** fy(3x2-y2) component ****
           lcount = lcount-1
           aorbs(k1,k2,k3,lcount)=D*GY*(3.0d0*(1.0d0-GZ*GZ)-4.0d0*GY*GY)
     >                          /dsqrt(24.0d0)

*          **** fxyz component ****
           lcount = lcount-1
           aorbs(k1,k2,k3,lcount)=D*GX*GY*GZ

*             **** fy(5z2-1) component ****
           lcount = lcount-1
           aorbs(k1,k2,k3,lcount)=D*GY*(5.0d0*GZ*GZ-1.0d0)
     >                          /dsqrt(40.0d0)

*          **** fz(5z2-3) component ****
           lcount = lcount-1
           aorbs(k1,k2,k3,lcount)=D*GZ*(5.0d0*GZ*GZ-3.0d0)
     >                          /dsqrt(60.0d0)

*          **** fx(5z2-1) component ****
           lcount = lcount-1
           aorbs(k1,k2,k3,lcount)=D*GX*(5.0d0*GZ*GZ-1.0d0)
     >                          /dsqrt(40.0d0)

*          **** fz(x2-y2) component ****
           lcount = lcount-1
           aorbs(k1,k2,k3,lcount)=D*GZ*(GX*GX - GY*GY)
     >                          /2.0d0

*          **** fx(x2-3y2) component ****
           lcount = lcount-1
           aorbs(k1,k2,k3,lcount)=D*GX*(4.0d0*GX*GX-3.0d0*(1.0d0-GZ*GZ))
     >                          /dsqrt(24.0d0)

        end if

*::::::::::::::::::::::::::::::  d-wave  ::::::::::::::::::::::::::::::
  300   CONTINUE
        if (locp.ne.2) then
          F(1)=0.0d0
          DO I=2,NRHO
            A=3.0d0*(SN(I)/(Q*RHO(I))-CS(I))/(Q*RHO(I))-SN(I)
            F(I)=A*WP(I,2)
          END DO
          D=P2*SIMP(NRHO,F,DRHO)/Q

*         **** dxy component ****
          lcount = lcount-1
          aorbs(k1,k2,k3,lcount)=D*GX*GY

*         **** dyz component ****
          lcount = lcount-1
          aorbs(k1,k2,k3,lcount)=D*GY*GZ

*         **** d3z2-1 component ****
          lcount = lcount-1
          aorbs(k1,k2,k3,lcount)=D*(3.0d0*GZ*GZ-1.0d0)
     >                          /(2.0d0*dsqrt(3.0d0))

*         **** dzx component ****
          lcount = lcount-1
          aorbs(k1,k2,k3,lcount)=D*GZ*GX

*         **** dx2-y2 component ****
          lcount = lcount-1
          aorbs(k1,k2,k3,lcount)=D*(GX*GX-GY*GY)/(2.0d0)
        end if

*::::::::::::::::::::::::::::::  p-wave  ::::::::::::::::::::::::::::::
  400   CONTINUE
        if (locp.ne.1) then
           F(1)=0.0d0
           DO I=2,NRHO
             F(I)=(SN(I)/(Q*RHO(I))-CS(I))*WP(I,1)
           END DO
           P=P1*SIMP(NRHO,F,DRHO)/Q

*          **** py component ****
           lcount = lcount-1
           aorbs(k1,k2,k3,lcount)=P*GY

*          **** pz component ****
           lcount = lcount-1
           aorbs(k1,k2,k3,lcount)=P*GZ

*          **** px component ****
           lcount = lcount-1
           aorbs(k1,k2,k3,lcount)=P*GX
        end if

*::::::::::::::::::::::::::::::  s-wave  :::::::::::::::::::::::::::::::
  500   CONTINUE
        if (locp.ne.0) then
          DO I=1,NRHO
            F(I)=SN(I)*WP(I,0)
          END DO
          lcount = lcount-1
          aorbs(k1,k2,k3,lcount)=P0*SIMP(NRHO,F,DRHO)/Q
        end if

*::::::::::::::::::::::::::::::  local  :::::::::::::::::::::::::::::::
  600   CONTINUE



  700 CONTINUE

      call Parallel_Vector_Sumall(lmmax*nfft3d,aorbs)
*:::::::::::::::::::::::::::::::  G=0  ::::::::::::::::::::::::::::::::      

      do l=1,lmmax
        aorbs(1,1,1,l)=0.0d0
      end do
*     *** only j0 is non-zero at zero ****
      DO  I=1,NRHO
         F(I)=RHO(I)*WP(I,0)
       END DO
      aorbs(1,1,1,1)=P0*SIMP(NRHO,F,DRHO)


      IERR=0
      RETURN
      END

