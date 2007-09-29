*
* $Id: integrate_paw_ray.F,v 1.2 2007-09-29 00:33:34 bylaska Exp $
*

*     *************************************************
*     *                                               *
*     *           integrate_paw_ray                   *
*     *                                               *
*     *************************************************

      subroutine integrate_paw_ray(version,rlocal,
     >                            nrho,log_amesh,zv,
     >                            nbasis,n_ps,l_ps,
     >                            vlocal,prj_ps,
     >                            rho,f,cs,sn,
     >                            nray,G_ray,vl_ray,vnl_ray,
     >                            ierr)
      implicit none
      integer          version
      double precision rlocal
      integer          nrho
      double precision log_amesh
      double precision zv
      integer nbasis,n_ps(*),l_ps(*)

      double precision vlocal(nrho)
      double precision prj_ps(nrho,nbasis)
      double precision rho(nrho)
      double precision f(nrho)
      double precision cs(nrho)
      double precision sn(nrho)

      integer nray
      double precision G_ray(nray)
      double precision vl_ray(nray)
      double precision vnl_ray(nray,nbasis)

      integer ierr

      integer np,taskid,MASTER
      parameter (MASTER=0)

*     *** local variables ****
      integer task_count
      integer k1,i,l,ib,la,power_f
      double precision pi,twopi,forpi
      double precision p0,p1,p2,p3,p
      double precision a,q,d

*     **** Error function parameters ****
      real*8 yerf,xerf

*     **** external functions ****
      double precision log_integrate_def,util_erf
      external         log_integrate_def,util_erf


      call Parallel_np(np)
      call Parallel_taskid(taskid)

      pi=4.0d0*datan(1.0d0)
      twopi=2.0d0*pi
      forpi=4.0d0*pi

      p0=dsqrt(forpi)
      p1=dsqrt(3.0d0*forpi)
      p2=dsqrt(15.0d0*forpi)
      p3=dsqrt(105.0d0*forpi)

*======================  Fourier transformation  ======================
      call dcopy(nray,0.0d0,0,vl_ray,1)
      call dcopy(nbasis*nray,0.0d0,0,vnl_ray,1)
      task_count = -1
      DO 700 k1=2,nray
        task_count = task_count + 1
        if (mod(task_count,np).ne.taskid) go to 700

        Q=G_ray(k1)

        do i=1,nrho
          cs(i)=dcos(Q*rho(i))
          sn(i)=dsin(Q*rho(i))
        end do

        do ib=1,nbasis
          la=l_ps(ib)

*         ::::::::::::::::::::::::::::::  h-wave  ::::::::::::::::::::::::::::::
*         ::::::::::::::::::::::::::::::  g-wave  ::::::::::::::::::::::::::::::
*         ::::::::::::::::::::::::::::::  f-wave  ::::::::::::::::::::::::::::::
          if (la.eq.3) then
             power_f = 8
             f(1)=0.0d0
             do i=2,nrho
               a=sn(i)/(Q*rho(i))
               a=15.0d0*(a-cs(i))/(Q*rho(i))**2 - 6*a + cs(i)
               f(i)=A*prj_ps(i,ib)
             end do
             d=p3*log_integrate_def(power_f,f,0,rho,log_amesh,nrho)/Q
             vnl_ray(k1,ib)=d
          end if
*         ::::::::::::::::::::::::::::::  d-wave  ::::::::::::::::::::::::::::::
          if (la.eq.2) then
             power_f = 6
             f(1)=0.0d0
             do i=2,nrho
               a=3.0d0*(sn(i)/(Q*rho(i))-cs(i))/(Q*rho(i))-sn(i)
               f(i)=A*prj_ps(i,ib)
             end do
             d=p2*log_integrate_def(power_f,f,0,rho,log_amesh,nrho)/Q
             vnl_ray(k1,ib)=D
          end if
*         ::::::::::::::::::::::::::::::  p-wave  ::::::::::::::::::::::::::::::
          if (la.eq.1) then
             power_f = 4
             f(1)=0.0d0
             do i=2,nrho
               f(i)=(sn(i)/(Q*rho(i))-cs(i))*prj_ps(i,ib)
             end do
             p=p1*log_integrate_def(power_f,f,0,rho,log_amesh,nrho)/Q
             vnl_ray(k1,ib)=p
          end if
*         ::::::::::::::::::::::::::::::  s-wave  ::::::::::::::::::::::::::::::
          if (la.eq.0) then
            power_f = 2
            do i=1,nrho
              f(i)=sn(i)*prj_ps(i,ib)
            end do
            vnl_ray(k1,ib)=
     >       p0*log_integrate_def(power_f,f,0,rho,log_amesh,nrho)/Q
          end if

        end do

*::::::::::::::::::::::::::::::  local  :::::::::::::::::::::::::::::::
        if (version.eq.3) then
        do  i=1,nrho
          xerf = rho(i)/rlocal
          yerf = util_erf(xerf)
          f(I)=(rho(i)*vlocal(i)-zv*yerf)*sn(i)
        end do
ccccc        vl_ray(k1)=SIMP(NRHO,F,DRHO)*forpi/Q-zv*forpi*CS(NRHO)/(Q*Q)
        vl_ray(k1)=forpi*log_integrate_def(0,f,0,rho,log_amesh,nrho)/Q
     >            -zv*forpi*cs(nrho)/(Q*Q)
        end if
 
        if (version.eq.4) then
        do i=1,nrho
          !xerf=rho(i)/rlocal
          !yerf = util_erf(xerf)
          !f(i)=(rho(i)*vlocal(i)+zv*yerf)*sn(i)
          f(i)=(rho(i)*vlocal(i))*sn(i)
        end do
cccc        vl_ray(k1)=SIMP(NRHO,F,DRHO)*forpi/Q
        vl_ray(k1)=forpi*log_integrate_def(0,f,0,rho,log_amesh,nrho)/Q
        end if


  700 CONTINUE
      call Parallel_Vector_SumAll(nray,vl_ray)
      call Parallel_Vector_SumAll(nbasis*nray,vnl_ray)

*:::::::::::::::::::::::::::::::  G=0  ::::::::::::::::::::::::::::::::      
      if (version.eq.3) then
      do i=1,nrho
        xerf = rho(i)/rlocal
        yerf = util_erf(xerf)
        f(i)=(vlocal(i)*rho(i)-zv*yerf)*rho(i)
      end do
      vl_ray(1)=forpi*log_integrate_def(0,f,0,rho,log_amesh,nrho)
     >         +twopi*zv*rho(nrho)**2
      end if


      if (version.eq.4) then
      do i=1,nrho
        !xerf=rho(i)/rlocal
        !yerf = util_erf(xerf)
        !f(i)=(vlocal(i)*rho(i)+zv*yerf)*rho(i)
        f(i)=(vlocal(i)*rho(i))*rho(i)
      end do
      vl_ray(1)=forpi*log_integrate_def(0,f,0,rho,log_amesh,nrho)
      end if

*     *** only j0 is non-zero at zero ****
      power_f = 2
      do ib=1,nbasis
        la=l_ps(ib)
        vnl_ray(1,ib)=0.0d0
        if (la.eq.0) then
           do i=1,nrho
             f(i)=rho(i)*prj_ps(i,ib)
           end do
           vnl_ray(1,ib)=
     >      p0*log_integrate_def(power_f,f,0,rho,log_amesh,nrho)
        end if
      end do

      ierr=0
      return
      end
