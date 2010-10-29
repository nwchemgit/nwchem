*
* $Id$
*

*     **************************************
*     *                                    *
*     *  integrate_kbpp_band_stress_local  *
*     *                                    *
*     **************************************

      subroutine integrate_kbpp_band_stress_local(version,
     >                            nrho,drho,lmax,locp,zv,
     >                            vp,wp,rho,f,cs,sn,
     >                            nfft1,nfft2,nfft3,lmmax,
     >                            G,dvl,
     >                            semicore,rho_sc_r,rho_sc_k,
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

      integer nfft1,nfft2,nfft3,lmmax
      double precision G(nfft1,nfft2,nfft3,3)
      double precision dvl(nfft1,nfft2,nfft3)


      logical semicore
      double precision rho_sc_r(nrho,2)
      double precision rho_sc_k(nfft1,nfft2,nfft3,4)

      integer ierr

      integer np,taskid,MASTER
      parameter (MASTER=0)

*     *** local variables ****
      integer task_count,nfft3d
      integer k1,k2,k3,i,l
      double precision pi,twopi,forpi
      double precision q


*     **** external functions ****
      double precision simp
      external         simp

      call Parallel_np(np)
      call Parallel_taskid(taskid)

      nfft3d = (nfft1)*nfft2*nfft3
      pi=4.0d0*datan(1.0d0)
      twopi=2.0d0*pi
      forpi=4.0d0*pi

      IF((nrho/2)*2.EQ.nrho) THEN
        ierr=2
        RETURN
      ENDIF


*::::::::::::::::::  Define non-local pseudopotential  ::::::::::::::::
      do l=0,lmax
        if (l.ne.locp) then
          do I=1,nrho
            vp(i,l)=vp(i,l)-vp(i,locp)
          end do
        end if
      end do

*======================  Fourier transformation  ======================
      call dcopy(nfft3d,0.0d0,0,dvl,1)
      call dcopy(4*nfft3d,0.0d0,0,rho_sc_k,1)
      task_count = -1
      DO 700 k3=1,nfft3
      DO 700 k2=1,nfft2
      DO 700 k1=1,nfft1
        task_count = task_count + 1
        if (mod(task_count,np).ne.taskid) go to 700

        q=dsqrt(G(k1,k2,k3,1)**2
     >         +G(k1,k2,k3,2)**2
     >         +G(k1,k2,k3,3)**2)

        if ((k1.eq.1).and.(k2.eq.1).and.(k3.eq.1)) go to 700

        DO I=1,nrho
          cs(I)=dcos(q*rho(I))
          sn(I)=dsin(q*rho(I))
        END DO

*::::::::::::::::::::::::::::::  local  :::::::::::::::::::::::::::::::
  600   CONTINUE

        do  i=1,nrho
          f(i)=rho(i)*vp(i,locp)*(rho(i)*cs(i)-sn(i)/q)
        end do
        dvl(k1,k2,k3)= simp(nrho,f,drho)*forpi/q
     >   + zv*forpi/(q*q)*(2.0d0*cs(nrho)/q + rho(nrho)*sn(nrho))

 

*::::::::::::::::::::: semicore density :::::::::::::::::::::::::::::::
        if (semicore) then

           do  i=1,nrho
             f(i)=rho(i)*dsqrt(rho_sc_r(i,1))*(rho(i)*cs(i)-sn(i)/q)
           end do
           rho_sc_k(k1,k2,k3,1)= simp(nrho,f,drho)*forpi/q
        end if

  700 CONTINUE

      call Parallel_Vector_SumAll(4*nfft3d,rho_sc_k)
      call Parallel_Vector_SumAll(nfft3d,dvl)
*:::::::::::::::::::::::::::::::  G=0  ::::::::::::::::::::::::::::::::      

      dvl(1,1,1)= 0.0d0

      IERR=0
      RETURN
      END

