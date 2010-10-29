*
* $Id$
*

*     ******************************************
*     *                                        *
*     *   integrate_kbpp_band_stress_nonlocal  *
*     *                                        *
*     ******************************************

      subroutine integrate_kbpp_band_stress_nonlocal(version,kvec,
     >                            nrho,drho,lmax,locp,zv,
     >                            vp,wp,rho,f,cs,sn,
     >                            nfft1,nfft2,nfft3,lmmax,
     >                            G,dvnl,
     >                            ierr)
      implicit none
      integer          version
      double precision kvec(3)
      integer          nrho
      double precision drho
      integer          lmax
      integer          locp
      double precision zv
      double precision vp(nrho,0:lmax)
      double precision wp(nrho,0:lmax)
      double precision rho(nrho)
      double precision f(nrho)
      double precision cs(nrho)
      double precision sn(nrho)

      integer nfft1,nfft2,nfft3,lmmax
      double precision G(nfft1,nfft2,nfft3,3)
      double precision dvnl(nfft1,nfft2,nfft3,3,lmmax)

      integer ierr

      integer np,taskid,MASTER
      parameter (MASTER=0)

*     *** local variables ****
      integer lcount,task_count,nfft3d
      integer k1,k2,k3,i,l
      double precision pi,twopi,forpi
      double precision p0,p1,p2,p3,p,pp
      double precision gx,gy,gz,a,q,d,dd
      double precision duxdGx,duxdGy,duxdGz
      double precision duydGx,duydGy,duydGz
      double precision duzdGx,duzdGy,duzdGz
      double precision sumx,sumy,sumz
      double precision T,dTdux,dTduy,dTduz



*     **** external functions ****
      double precision dsum,simp
      external         dsum,simp

      call Parallel_np(np)
      call Parallel_taskid(taskid)

      nfft3d = (nfft1)*nfft2*nfft3
      pi=4.0d0*datan(1.0d0)
      twopi=2.0d0*pi
      forpi=4.0d0*pi

      if (lmmax.GT.16) then
        ierr=1
        return
      end if
      if ((nrho/2)*2.EQ.nrho) then
        ierr=2
        return
      end if

      p0=dsqrt(forpi)
      p1=dsqrt(3.0d0*forpi)
      p2=dsqrt(15.0d0*forpi)
      p3=dsqrt(105.0d0*forpi)

*======================  Fourier transformation  ======================
      call dcopy(3*lmmax*nfft3d,0.0d0,0,dvnl,1)
      task_count = -1
      do 700 k3=1,nfft3
      do 700 k2=1,nfft2
      do 700 k1=1,nfft1
        task_count = task_count + 1
        if (mod(task_count,np).ne.taskid) go to 700
        gx=G(k1,k2,k3,1)+kvec(1)
        gy=G(k1,k2,k3,2)+kvec(2)
        gz=G(k1,k2,k3,3)+kvec(3)

        q=dsqrt(gx**2 + gy**2 + gz**2)

        if (dabs(q).gt.1.0d-9) then

           gx=gx/q
           gy=gy/q
           gz=gz/q
           do i=1,nrho
             cs(i)=dcos(q*rho(i))
             sn(i)=dsin(q*rho(i))
           end do

*          **** calculate du_r/dG_s ****
           duxdGx = 1.0d0/q -gx*gx/q
           duxdGy = -gx*gy/q
           duxdGz = -gx*gz/q

           duydGx = -gy*gx/q
           duydGy = 1.0d0/q - gy*gy/q
           duydGz = -gy*gz/q

           duzdGx = -gz*gx/q
           duzdGy = -gz*gy/q
           duzdGz = 1.0d0/q - gz*gz/q

           lcount = lmmax+1
           GO TO (500,400,300,200), LMAX+1


*::::::::::::::::::::::::::::::  f-wave  ::::::::::::::::::::::::::::::
  200      CONTINUE
        if (locp.ne.3) then
           F(1)=0.0d0
           do i=2,nrho
             A=sn(i)/(q*rho(i))
             A=15.0d0*(A-cs(i))/(q*rho(i))**2 - 6*A + cs(i)
             f(i)=A*wp(i,3)*vp(i,3)
           end do
           D=P3*SIMP(nrho,F,drho)/q

           F(1)=0.0d0
           do i=2,nrho
             A= -60.0d0*sn(i)/(rho(i)**3 * q**5)
     >        +  60.0d0*cs(i)/(rho(i)**2 * q**4)
     >        +  27.0d0*sn(i)/(rho(i)    * q**3)
     >        -   7.0d0*cs(i)/(q**2)
     >        -   rho(i)*sn(i)/q
             f(i)=A*wp(i,3)*vp(i,3)
           end do
           DD=P3*SIMP(nrho,F,drho)

           lcount = lcount-1
           T = gx*(4.0d0*gx*gx - 3.0d0*(1.0d0-gz*gz))/dsqrt(24.0d0)
           dTdux = (12.0d0*gx*gx-3.0d0*(1.0d0-gz*gz))/dsqrt(24.0d0)
           dTduy = 0.0d0
           dTduz = 6.0d0*gx*gz/dsqrt(24.0d0)
           sumx = dTdux*duxdGx + dTduy*duydGx + dTduz*duzdGx
           sumy = dTdux*duxdGy + dTduy*duydGy + dTduz*duzdGy
           sumz = dTdux*duxdGz + dTduy*duydGz + dTduz*duzdGz
           dvnl(k1,k2,k3,1,lcount)=DD*T*gx + D*sumx
           dvnl(k1,k2,k3,2,lcount)=DD*T*gy + D*sumy
           dvnl(k1,k2,k3,3,lcount)=DD*T*gz + D*sumz

           lcount = lcount-1
           T = gy*(3.0d0*(1.0d0-gz*gz)-4.0d0*gy*gy)/dsqrt(24.0d0)
           dTdux = 0.0d0
           dTduy = (3.0d0*(1.0d0-gz*gz)-12.0d0*gy*gy)/dsqrt(24.0d0)
           dTduz = -6.0d0*gy*gz/dsqrt(24.0d0)
           sumx = dTdux*duxdGx + dTduy*duydGx + dTduz*duzdGx
           sumy = dTdux*duxdGy + dTduy*duydGy + dTduz*duzdGy
           sumz = dTdux*duxdGz + dTduy*duydGz + dTduz*duzdGz
           dvnl(k1,k2,k3,1,lcount)=DD*T*gx + D*sumx
           dvnl(k1,k2,k3,2,lcount)=DD*T*gy + D*sumy
           dvnl(k1,k2,k3,3,lcount)=DD*T*gz + D*sumz

           lcount = lcount-1
           T =gz*(gx*gx - gy*gy)/2.0d0
           dTdux =  gx*gz
           dTduy = -gy*gz
           dTduz = (gx*gx-gy*gy)/2.0d0
           sumx = dTdux*duxdGx + dTduy*duydGx + dTduz*duzdGx
           sumy = dTdux*duxdGy + dTduy*duydGy + dTduz*duzdGy
           sumz = dTdux*duxdGz + dTduy*duydGz + dTduz*duzdGz
           dvnl(k1,k2,k3,1,lcount)=DD*T*gx + D*sumx
           dvnl(k1,k2,k3,2,lcount)=DD*T*gy + D*sumy
           dvnl(k1,k2,k3,3,lcount)=DD*T*gz + D*sumz

           lcount = lcount-1
           T =gx*gy*gz
           dTdux = gy*gz
           dTduy = gx*gz
           dTduz = gx*gy
           sumx = dTdux*duxdGx + dTduy*duydGx + dTduz*duzdGx
           sumy = dTdux*duxdGy + dTduy*duydGy + dTduz*duzdGy
           sumz = dTdux*duxdGz + dTduy*duydGz + dTduz*duzdGz
           dvnl(k1,k2,k3,1,lcount)=DD*T*gx + D*sumx
           dvnl(k1,k2,k3,2,lcount)=DD*T*gy + D*sumy
           dvnl(k1,k2,k3,3,lcount)=DD*T*gz + D*sumz

           lcount = lcount-1
           T =  gx*(5.0d0*gz*gz-1.0d0)/dsqrt(40.0d0)
           dTdux = (5.0d0*gz*gz-1.0d0)/dsqrt(40.0d0)
           dTduy = 0.0d0
           dTduz = 10.0d0*gx*gz/dsqrt(40.0d0)
           sumx = dTdux*duxdGx + dTduy*duydGx + dTduz*duzdGx
           sumy = dTdux*duxdGy + dTduy*duydGy + dTduz*duzdGy
           sumz = dTdux*duxdGz + dTduy*duydGz + dTduz*duzdGz
           dvnl(k1,k2,k3,1,lcount)=DD*T*gx + D*sumx
           dvnl(k1,k2,k3,2,lcount)=DD*T*gy + D*sumy
           dvnl(k1,k2,k3,3,lcount)=DD*T*gz + D*sumz

           lcount = lcount-1
           T = gy*(5.0d0*gz*gz-1.0d0)/dsqrt(40.0d0)
           dTdux = 0.0d0
           dTduy =(5.0d0*gz*gz-1.0d0)/dsqrt(40.0d0)
           dTduz =10.0d0*gy*gz/dsqrt(40.0d0)
           sumx = dTdux*duxdGx + dTduy*duydGx + dTduz*duzdGx
           sumy = dTdux*duxdGy + dTduy*duydGy + dTduz*duzdGy
           sumz = dTdux*duxdGz + dTduy*duydGz + dTduz*duzdGz
           dvnl(k1,k2,k3,1,lcount)=DD*T*gx + D*sumx
           dvnl(k1,k2,k3,2,lcount)=DD*T*gy + D*sumy
           dvnl(k1,k2,k3,3,lcount)=DD*T*gz + D*sumz

           lcount = lcount-1
           T =gz*(5.0d0*gz*gz-3.0d0)/dsqrt(60.0d0)
           dTdux = 0.0d0
           dTduy = 0.0d0
           dTduz =(15.0d0*gz*gz -3.0d0)/dsqrt(60.0d0)
           sumx = dTdux*duxdGx + dTduy*duydGx + dTduz*duzdGx
           sumy = dTdux*duxdGy + dTduy*duydGy + dTduz*duzdGy
           sumz = dTdux*duxdGz + dTduy*duydGz + dTduz*duzdGz
           dvnl(k1,k2,k3,1,lcount)=DD*T*gx + D*sumx
           dvnl(k1,k2,k3,2,lcount)=DD*T*gy + D*sumy
           dvnl(k1,k2,k3,3,lcount)=DD*T*gz + D*sumz
        end if



*::::::::::::::::::::::::::::::  d-wave  ::::::::::::::::::::::::::::::
  300      CONTINUE
        if (locp.ne.2) then
          F(1)=0.0d0
          DO i=2,nrho
            A=3.0d0*(sn(i)/(q*rho(i))-cs(i))/(q*rho(i))-sn(i)
            f(i)=A*wp(i,2)*vp(i,2)
          END DO
          D=P2*SIMP(nrho,F,drho)/q

          F(1)=0.0d0
          DO i=2,nrho
            A= -9.0d0*sn(i)/(rho(i)**2 * q**4)
     >       +  9.0d0*cs(i)/(rho(i)    * q**3)
     >       +  4.0d0*sn(i)/(q**2)
     >       -  rho(i)*cs(i)/q
            f(i)=A*wp(i,2)*vp(i,2)
          END DO
          DD=P2*SIMP(nrho,F,drho)

          lcount = lcount-1
          T = (3.0d0*gz*gz-1.0d0)/(2.0d0*dsqrt(3.0d0))
          dTdux = 0.0d0
          dTduy = 0.0d0
          dTduz = 6.0d0*gz/(2.0d0*dsqrt(3.0d0))
          sumx = dTdux*duxdGx + dTduy*duydGx + dTduz*duzdGx
          sumy = dTdux*duxdGy + dTduy*duydGy + dTduz*duzdGy
          sumz = dTdux*duxdGz + dTduy*duydGz + dTduz*duzdGz
          dvnl(k1,k2,k3,1,lcount)=DD*T*gx + D*sumx
          dvnl(k1,k2,k3,2,lcount)=DD*T*gy + D*sumy
          dvnl(k1,k2,k3,3,lcount)=DD*T*gz + D*sumz

          lcount = lcount-1
          T = gx*gy
          dTdux = gy
          dTduy = gx
          dTduz = 0.0d0
          sumx = dTdux*duxdGx + dTduy*duydGx + dTduz*duzdGx
          sumy = dTdux*duxdGy + dTduy*duydGy + dTduz*duzdGy
          sumz = dTdux*duxdGz + dTduy*duydGz + dTduz*duzdGz
          dvnl(k1,k2,k3,1,lcount)=DD*T*gx + D*sumx
          dvnl(k1,k2,k3,2,lcount)=DD*T*gy + D*sumy
          dvnl(k1,k2,k3,3,lcount)=DD*T*gz + D*sumz

          lcount = lcount-1
          T = gy*gz
          dTdux = 0.0d0
          dTduy = gz
          dTduz = gy
          sumx = dTdux*duxdGx + dTduy*duydGx + dTduz*duzdGx
          sumy = dTdux*duxdGy + dTduy*duydGy + dTduz*duzdGy
          sumz = dTdux*duxdGz + dTduy*duydGz + dTduz*duzdGz
          dvnl(k1,k2,k3,1,lcount)=DD*T*gx + D*sumx
          dvnl(k1,k2,k3,2,lcount)=DD*T*gy + D*sumy
          dvnl(k1,k2,k3,3,lcount)=DD*T*gz + D*sumz

          lcount = lcount-1
          T = gz*gx
          dTdux = gz
          dTduy = 0.0d0
          dTduz = gx
          sumx = dTdux*duxdGx + dTduy*duydGx + dTduz*duzdGx
          sumy = dTdux*duxdGy + dTduy*duydGy + dTduz*duzdGy
          sumz = dTdux*duxdGz + dTduy*duydGz + dTduz*duzdGz
          dvnl(k1,k2,k3,1,lcount)=DD*T*gx + D*sumx
          dvnl(k1,k2,k3,2,lcount)=DD*T*gy + D*sumy
          dvnl(k1,k2,k3,3,lcount)=DD*T*gz + D*sumz

          lcount = lcount-1
          T = (gx*gx-gy*gy)/2.0d0
          dTdux = gx
          dTduy = -gy
          dTduz = 0.0d0
          sumx = dTdux*duxdGx + dTduy*duydGx + dTduz*duzdGx
          sumy = dTdux*duxdGy + dTduy*duydGy + dTduz*duzdGy
          sumz = dTdux*duxdGz + dTduy*duydGz + dTduz*duzdGz
          dvnl(k1,k2,k3,1,lcount)=DD*T*gx + D*sumx
          dvnl(k1,k2,k3,2,lcount)=DD*T*gy + D*sumy
          dvnl(k1,k2,k3,3,lcount)=DD*T*gz + D*sumz
        end if


*::::::::::::::::::::::::::::::  p-wave  ::::::::::::::::::::::::::::::
  400      CONTINUE
        if (locp.ne.1) then
           F(1)=0.0d0
           DO i=2,nrho
             f(i)=(sn(i)/(q*rho(i)) - cs(i)) * wp(i,1)*vp(i,1)
           END DO
           P=P1*SIMP(nrho,F,drho)/q

           F(1)=0.0d0
           DO i=2,nrho
             f(i)=wp(i,1)*vp(i,1)* ( -2.0d0*sn(i)/(rho(i) * q**3)
     >                              + 2.0d0*cs(i)/(q**2)
     >                              + rho(i)*sn(i)/q)
           END DO
           PP=P1*SIMP(nrho,F,drho)

           lcount = lcount-1
           T = gx
           dTdux = 1.0d0
           dTduy = 0.0d0
           dTduz = 0.0d0
           sumx = dTdux*duxdGx + dTduy*duydGx + dTduz*duzdGx
           sumy = dTdux*duxdGy + dTduy*duydGy + dTduz*duzdGy
           sumz = dTdux*duxdGz + dTduy*duydGz + dTduz*duzdGz
           dvnl(k1,k2,k3,1,lcount)= PP*T*gx + P*sumx
           dvnl(k1,k2,k3,2,lcount)= PP*T*gy + P*sumy
           dvnl(k1,k2,k3,3,lcount)= PP*T*gz + P*sumz


           lcount = lcount-1
           T = gy
           dTdux = 0.0d0
           dTduy = 1.0d0
           dTduz = 0.0d0
           sumx = dTdux*duxdGx + dTduy*duydGx + dTduz*duzdGx
           sumy = dTdux*duxdGy + dTduy*duydGy + dTduz*duzdGy
           sumz = dTdux*duxdGz + dTduy*duydGz + dTduz*duzdGz
           dvnl(k1,k2,k3,1,lcount)= PP*T*gx + P*sumx
           dvnl(k1,k2,k3,2,lcount)= PP*T*gy + P*sumy
           dvnl(k1,k2,k3,3,lcount)= PP*T*gz + P*sumz

           lcount = lcount-1
           T = gz
           dTdux = 0.0d0
           dTduy = 0.0d0
           dTduz = 1.0d0
           sumx = dTdux*duxdGx + dTduy*duydGx + dTduz*duzdGx
           sumy = dTdux*duxdGy + dTduy*duydGy + dTduz*duzdGy
           sumz = dTdux*duxdGz + dTduy*duydGz + dTduz*duzdGz
           dvnl(k1,k2,k3,1,lcount)= PP*T*gx + P*sumx
           dvnl(k1,k2,k3,2,lcount)= PP*T*gy + P*sumy
           dvnl(k1,k2,k3,3,lcount)= PP*T*gz + P*sumz
        end if



*::::::::::::::::::::::::::::::  s-wave  :::::::::::::::::::::::::::::::
  500      CONTINUE
        if (locp.ne.0) then
          DO i=1,nrho
            f(i)=wp(i,0)*vp(i,0) * ( -sn(i)/(q**2)
     >                              + rho(i)*cs(i)/q)
          END DO
          P = P0*SIMP(nrho,F,drho)
          lcount = lcount-1
          dvnl(k1,k2,k3,1,lcount) = P *gx
          dvnl(k1,k2,k3,2,lcount) = P *gy
          dvnl(k1,k2,k3,3,lcount) = P *gz
        end if


  600      CONTINUE

*:::::::::::::::::::::::::::::::  G+k=0  ::::::::::::::::::::::::::::::::      
      else

        do l=1,lmmax
          dvnl(1,1,1,1,l)=0.0d0
          dvnl(1,1,1,2,l)=0.0d0
          dvnl(1,1,1,3,l)=0.0d0
        end do

       end if

    
  700 CONTINUE

      call Parallel_Vector_SumAll(3*lmmax*nfft3d,dvnl)

      IERR=0
      RETURN
      END



