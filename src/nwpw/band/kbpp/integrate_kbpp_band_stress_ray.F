*
* $Id$
*

*     ********************************************
*     *                                          *
*     *        integrate_kbpp_band_stress_ray    *
*     *                                          *
*     ********************************************

      subroutine integrate_kbpp_band_stress_ray(version,
     >                            nrho,drho,lmax,locp,zv,
     >                            vp,wp,rho,f,cs,sn,
     >                            nray,G_ray,dvl_ray,dvnl_ray,
     >                            semicore,rho_sc_r,rho_sc_k_ray,
     >                            ierr)
      implicit none
      integer          version
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

      integer nray
      double precision G_ray(nray)
      double precision dvl_ray(nray)
      double precision dvnl_ray(nray,2,0:lmax)

      logical semicore
      double precision rho_sc_r(nrho,2)
      double precision rho_sc_k_ray(nray)
      integer ierr

#include "errquit.fh"

*     *** local variables ****
      integer np,taskid,MASTER
      parameter (MASTER=0)

      integer task_count
      integer k1,i,l
      double precision pi,twopi,forpi
      double precision p0,p1,p2,p3
      double precision a,q
      logical small_cell


*     **** external functions ****
      double precision dsum,simp
      external         dsum,simp
      logical          control_psp_semicore_small

      call Parallel_np(np)
      call Parallel_taskid(taskid)

      small_cell = control_psp_semicore_small()

      if(lmax.gt.3) THEN
         call errquit('integrate_kbpp_band_stress_ray - lmax > f',0,
     >       INPUT_ERR)
        ierr=1
        return
      end if
      if((nrho/2)*2.eq.nrho) then
       call errquit('integrate_kbpp_band_stress_ray - psp grid not odd',
     >       0,INPUT_ERR)
        ierr=2
        return
      end if

      pi=4.0d0*datan(1.0d0)
      twopi=2.0d0*pi
      forpi=4.0d0*pi
      p0=dsqrt(forpi)
      p1=dsqrt(3.0d0*forpi)
      p2=dsqrt(15.0d0*forpi)
      p3=dsqrt(105.0d0*forpi)

*======================  Fourier transformation  ======================
      call dcopy(nray,0.0d0,0,dvl_ray,1)
      call dcopy(2*(lmax+1)*nray,0.0d0,0,dvnl_ray,1)
      call dcopy(nray,0.0d0,0,rho_sc_k_ray,1)
      task_count = -1
      DO 700 k1=2,nray
        task_count = task_count + 1
        if (mod(task_count,np).ne.taskid) go to 700

        q=G_ray(k1)
        
        do i=1,nrho
          cs(i)=DCOS(q*rho(i))
          sn(i)=DSIN(q*rho(i))
        end do

        GO TO (500,400,300,200), lmax+1

*::::::::::::::::::::::::::::::  f-wave  ::::::::::::::::::::::::::::::
  200   CONTINUE
        if (locp.ne.3) then
           F(1)=0.0d0
           do i=2,nrho
             A=sn(i)/(q*rho(i))
             A=15.0d0*(A-cs(i))/(q*rho(i))**2 - 6*A + cs(i)
             f(i)=A*wp(i,3)*vp(i,3)
           end do
           dvnl_ray(k1,1,3)=p3*SIMP(nrho,F,drho)/q

           F(1)=0.0d0
           do i=2,nrho
             A= -60.0d0*sn(i)/(rho(i)**3 * q**5)
     >        +  60.0d0*cs(i)/(rho(i)**2 * q**4)
     >        +  27.0d0*sn(i)/(rho(i)    * q**3)
     >        -   7.0d0*cs(i)/(q**2)
     >        -   rho(i)*sn(i)/q
             f(i)=A*wp(i,3)*vp(i,3)
           end do
           dvnl_ray(k1,2,3)=p3*SIMP(nrho,F,drho)
        end if
*::::::::::::::::::::::::::::::  d-wave  ::::::::::::::::::::::::::::::
  300   CONTINUE
        if (locp.ne.2) then
          F(1)=0.0d0
          DO i=2,nrho
            A=3.0d0*(sn(i)/(q*rho(i))-cs(i))/(q*rho(i))-sn(i)
            f(i)=A*wp(i,2)*vp(i,2)
          END DO
          dvnl_ray(k1,1,2)=p2*SIMP(nrho,F,drho)/q

          F(1)=0.0d0
          DO i=2,nrho
            A= -9.0d0*sn(i)/(rho(i)**2 * q**4)
     >       +  9.0d0*cs(i)/(rho(i)    * q**3)
     >       +  4.0d0*sn(i)/(q**2)
     >       -  rho(i)*cs(i)/q
            f(i)=A*wp(i,2)*vp(i,2)
          END DO
          dvnl_ray(k1,2,2)=p2*SIMP(nrho,F,drho)
        end if
*::::::::::::::::::::::::::::::  p-wave  ::::::::::::::::::::::::::::::
  400   CONTINUE
        if (locp.ne.1) then
           F(1)=0.0d0
           DO i=2,nrho
             f(i)=(sn(i)/(q*rho(i)) - cs(i)) * wp(i,1)*vp(i,1)
           END DO
           dvnl_ray(k1,1,1)=p1*SIMP(nrho,F,drho)/q

           F(1)=0.0d0
           DO i=2,nrho
             f(i)=wp(i,1)*vp(i,1)* ( -2.0d0*sn(i)/(rho(i) * q**3)
     >                              + 2.0d0*cs(i)/(q**2)
     >                              + rho(i)*sn(i)/q)
           END DO
           dvnl_ray(k1,2,1)=p1*SIMP(nrho,F,drho)
        end if
*::::::::::::::::::::::::::::::  s-wave  :::::::::::::::::::::::::::::::
  500   CONTINUE
        if (locp.ne.0) then
          DO i=1,nrho
            f(i)=wp(i,0)*vp(i,0) * ( -sn(i)/(q**2) 
     >                              + rho(i)*cs(i)/q)
          END DO
          dvnl_ray(k1,1,0) = p0*SIMP(nrho,F,drho)
        end if
*::::::::::::::::::::::::::::::  local  :::::::::::::::::::::::::::::::
  600   CONTINUE

        do  i=1,nrho
          f(i)=rho(i)*vp(i,locp)*(rho(i)*cs(i)-sn(i)/q)
        end do
        dvl_ray(k1)= SIMP(nrho,f,drho)*forpi/q
     >   + zv*forpi/(q*q)*(2.0d0*cs(nrho)/q + rho(nrho)*sn(nrho))
*::::::::::::::::::::: semicore density :::::::::::::::::::::::::::::::
        if (semicore) then
           if (small_cell) then
              do  i=1,nrho
                f(i)=rho(i)*(rho_sc_r(i,1))*(rho(i)*cs(i)-sn(i)/q)
              end do
           else
              do  i=1,nrho
                f(i)=rho(i)*dsqrt(rho_sc_r(i,1))*(rho(i)*cs(i)-sn(i)/q)
              end do
           end if
           rho_sc_k_ray(k1)= SIMP(nrho,f,drho)*forpi/q
        end if
    
  700 CONTINUE

      call Parallel_Vector_SumAll(nray,rho_sc_k_ray)
      call Parallel_Vector_SumAll(nray,dvl_ray)
      call Parallel_Vector_Sumall(2*(lmax+1)*nray,dvnl_ray)

*:::::::::::::::::::::::::::::::  G=0  ::::::::::::::::::::::::::::::::      
      dvl_ray(1)      = 0.0d0
      rho_sc_k_ray(1) = 0.0d0
      do l=0,lmax
        dvnl_ray(1,1,l)=0.0d0
        dvnl_ray(1,2,l)=0.0d0
      end do

      ierr=0
      return
      end



