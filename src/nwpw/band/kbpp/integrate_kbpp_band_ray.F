*
* $Id$
*
      subroutine integrate_kbpp_band_ray(version,
     >                            nrho,drho,lmax,locp,zv,
     >                            vp,wp,rho,f,cs,sn,
     >                            nray,G_ray,vl_ray,vnl_ray,
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
      double precision vl_ray(nray)
      double precision vnl_ray(nray,0:lmax)

      logical semicore
      double precision rho_sc_r(nrho,2)
      double precision rho_sc_k_ray(nray,2)
      integer ierr

#include "errquit.fh"

*     *** local variables ****
      integer np,taskid,MASTER
      parameter (MASTER=0)

      integer lcount,task_count
      integer k1,i,l
      double precision pi,twopi,forpi
      double precision p0,p1,p2,p3,p
      double precision a,q
      logical small_cell


*     **** external functions ****
      double precision dsum,simp
      external         dsum,simp
      logical          control_psp_semicore_small
      external         control_psp_semicore_small

      call Parallel_np(np)
      call Parallel_taskid(taskid)

      small_cell = control_psp_semicore_small()

      if((nrho/2)*2.eq.nrho) then
        call errquit('integrate_kbpp_band_ray - psp grid not odd',0,
     >       INPUT_ERR)
        ierr = 2
        return
      end if

      pi=4.0d0*datan(1.0d0)
      twopi=2.0d0*pi
      forpi=4.0d0*pi
      p0=DSQRT(forpi)
      p1=DSQRT(3.0d0*forpi)
      p2=DSQRT(15.0d0*forpi)
      p3=DSQRT(105.0d0*forpi)


*======================  Fourier transformation  ======================
      call dcopy(nray,0.0d0,0,vl_ray,1)
      call dcopy((lmax+1)*nray,0.0d0,0,vnl_ray,1)
      call dcopy(2*nray,0.0d0,0,rho_sc_k_ray,1)
      task_count = -1
      DO 700 k1=2,nray
        task_count = task_count + 1
        if (mod(task_count,np).ne.taskid) go to 700

        Q=G_ray(k1)

        do i=1,nrho
          cs(i)=DCOS(Q*rho(i))
          sn(i)=DSIN(Q*rho(i))
        end do

        GO TO (500,400,300,200), lmax+1

*::::::::::::::::::::::::::::::  f-wave  ::::::::::::::::::::::::::::::
  200   CONTINUE
        if (locp.ne.3) then
           F(1)=0.0d0
           do I=2,NRHO
             A=SN(I)/(Q*RHO(I))
             A=15.0d0*(A-CS(I))/(Q*RHO(I))**2 - 6*A + CS(I)
             F(I)=A*WP(I,3)*VP(I,3)
           end do
          vnl_ray(k1,3)=P3*SIMP(NRHO,F,DRHO)/Q
        end if

*::::::::::::::::::::::::::::::  d-wave  ::::::::::::::::::::::::::::::
  300   CONTINUE
        if (locp.ne.2) then
          F(1)=0.0d0
          DO I=2,NRHO
            A=3.0d0*(SN(I)/(Q*RHO(I))-CS(I))/(Q*RHO(I))-SN(I)
            F(I)=A*WP(I,2)*VP(I,2)
          END DO
          vnl_ray(k1,2)=P2*SIMP(NRHO,F,DRHO)/Q
        end if

*::::::::::::::::::::::::::::::  p-wave  ::::::::::::::::::::::::::::::
  400   CONTINUE
        if (locp.ne.1) then
           F(1)=0.0d0
           DO I=2,NRHO
             F(I)=(SN(I)/(Q*RHO(I))-CS(I))*WP(I,1)*VP(I,1)
           END DO
           vnl_ray(k1,1)=P1*SIMP(NRHO,F,DRHO)/Q
        end if

*::::::::::::::::::::::::::::::  s-wave  :::::::::::::::::::::::::::::::
  500   CONTINUE
        if (locp.ne.0) then
          DO I=1,NRHO
            F(I)=SN(I)*WP(I,0)*VP(I,0)
          END DO
          lcount = lcount-1
        vnl_ray(k1,0)=P0*SIMP(NRHO,F,DRHO)/Q
        end if

*::::::::::::::::::::::::::  local potentials  ::::::::::::::::::::::::::
  600   CONTINUE
        do i=1,nrho
          f(i)=rho(i)*vp(i,locp)*sn(i)
        end do
        vl_ray(k1)=SIMP(nrho,f,drho)*FORPI/Q-ZV*FORPI*cs(nrho)/(Q*Q)

*::::::::::::::::::::: semicore density :::::::::::::::::::::::::::::::
        if (semicore) then
           if (small_cell) then
              do i=1,nrho
                 f(i) = rho(i)*(rho_sc_r(i,1))*sn(i)
              end do
           else
              do i=1,nrho
                 f(i) = rho(i)*dsqrt(rho_sc_r(i,1))*sn(i)
              end do
           end if
           rho_sc_k_ray(k1,1) = SIMP(nrho,f,drho)*forpi/Q

           do i=1,nrho
             f(i)=(sn(i)/(Q*rho(i))-cs(i))*rho_sc_r(i,2)*rho(i)
           end do
           P = SIMP(nrho,f,drho)*forpi/Q
           rho_sc_k_ray(k1,2)=P
        end if

  700 CONTINUE
      call Parallel_Vector_SumAll(2*nray,rho_sc_k_ray)
      call Parallel_Vector_SumAll(nray,vl_ray)
      call Parallel_Vector_Sumall((lmax+1)*nray,vnl_ray)

*:::::::::::::::::::::::::::::::  G=0  ::::::::::::::::::::::::::::::::
      do i=1,nrho
        f(i)=vp(i,locp)*rho(i)**2
      END DO
      vl_ray(1)=FORPI*SIMP(nrho,f,drho)+TWOPI*ZV*rho(nrho)**2

*     **** semicore density ****
      if (semicore) then
         if (small_cell) then
            do i=1,nrho
               f(i) = (rho_sc_r(i,1))*rho(i)**2
            end do
         else
            do i=1,nrho
               f(i) = dsqrt(rho_sc_r(i,1))*rho(i)**2
            end do
         end if
         rho_sc_k_ray(1,1) = forpi*SIMP(nrho,f,drho)
         rho_sc_k_ray(1,2) = 0.0d0
      end if

      do l=0,lmax
        vnl_ray(1,l)=0.0d0
      end do
*     *** only j0 is non-zero at zero ****
      if (locp.ne.0) then
         DO  I=1,NRHO
           F(I)=RHO(I)*WP(I,0)*VP(I,0)
         END DO
         vnl_ray(1,0)=P0*SIMP(NRHO,F,DRHO)
      end if


      IERR=0
      RETURN
      END



