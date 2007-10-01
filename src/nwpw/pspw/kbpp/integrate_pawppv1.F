*
* $Id: integrate_pawppv1.F,v 1.3 2007-10-01 23:02:41 bylaska Exp $
*
*     *****************************************************
*     *                                                   *
*     *            integrate_pawppv1                      *
*     *                                                   *
*     *****************************************************
      subroutine integrate_pawppv1(version,rlocal,
     >                            nrho,log_amesh,nbasis,nmax,lmax,zv,zc,
     >                            n_ps,l_ps,
     >                            v_ps,prj_ps,
     >                            phi_ae,dphi_ae,phi_ps,dphi_ps,
     >                            rho,f,cs,sn,
     >                            nfft1,nfft2,nfft3,nproj,
     >                            G,vl,vnl,
     >                            n_prj,l_prj,m_prj,
     >                            Gijl,
     >                            nray,G_ray,vl_ray,vnl_ray,
     >                            tmp_ray,
     >                            filter,
     >                            ierr)
      implicit none
      integer          version
      double precision rlocal
      integer          nrho
      double precision log_amesh
      integer          nbasis,nmax,lmax
      double precision zv,zc
      integer          n_ps(*),l_ps(*)
      double precision v_ps(nrho)
      double precision prj_ps(nrho,nbasis)

      double precision phi_ae(nrho,nbasis)
      double precision dphi_ae(nrho,nbasis)
      double precision phi_ps(nrho,nbasis)
      double precision dphi_ps(nrho,nbasis)

      double precision rho(nrho)
      double precision f(nrho)
      double precision cs(nrho)
      double precision sn(nrho)

      integer nfft1,nfft2,nfft3,nproj
      double precision G(nfft1/2+1,nfft2,nfft3,3)
      double precision vl(nfft1/2+1,nfft2,nfft3)
      double precision vnl(nfft1/2+1,nfft2,nfft3,nproj)
      integer          n_prj(nproj),l_prj(nproj),m_prj(nproj)
      double precision Gijl(nmax,nmax,0:lmax,5)

      integer nray
      double precision G_ray(nray)
      double precision vl_ray(nray,2)
      double precision vnl_ray(nray,nbasis,2)
      double precision tmp_ray(nray)
      logical filter

      integer ierr

      integer np,taskid,MASTER
      parameter (MASTER=0)

*     *** local variables ****
      integer lcount,task_count,nfft3d
      integer k1,k2,k3,i,j,k,l,nx,la,na,lb,nb,power_f
      double precision pi,twopi,forpi
      double precision p0,p1,p2,p3,p
      double precision gx,gy,gz,a,q,d
      double precision ecut,wcut,dG,yp1
      double precision xerf,yerf

*     **** external functions ****
      double precision control_ecut,control_wcut,util_erf
      double precision nwpw_splint,log_integrate_def
      external         control_ecut,control_wcut,util_erf
      external         nwpw_splint,log_integrate_def

      call Parallel_np(np)
      call Parallel_taskid(taskid)

      nfft3d = (nfft1/2+1)*nfft2*nfft3
      pi=4.0d0*datan(1.0d0)
      twopi=2.0d0*pi
      forpi=4.0d0*pi

c      P0=dsqrt(forpi)
c      P1=dsqrt(3.0d0*forpi)
c      P2=dsqrt(15.0d0*forpi)
c      P3=dsqrt(105.0d0*forpi)

*     ***********  Normarization constants  *************
      call dcopy(5*nmax*nmax*(lmax+1),0.0d0,0,Gijl,1)

*     **** overlap ***
      do i=1,nbasis
         la = l_ps(i)
         na = n_ps(i) - la
         power_f = 2*la + 2
         do k=1,nrho
           f(k) = phi_ae(k,i)*phi_ae(k,i)
     >          - phi_ps(k,i)*phi_ps(k,i)
         end do
         d = log_integrate_def(power_f,f,0,rho,log_amesh,nrho)
         Gijl(na,na,la,2) = Gijl(na,na,la,2) + d

         do j=1,i-1
            lb = l_ps(j)
            nb = n_ps(j) - lb
            if (la.eq.lb) then
              do k=1,nrho
                f(k) = phi_ae(k,i)*phi_ae(k,j)
     >               - phi_ps(k,i)*phi_ps(k,j)
              end do
              d = log_integrate_def(power_f,f,0,rho,log_amesh,nrho)
              Gijl(na,nb,la,2) = Gijl(na,nb,la,2) + d
              Gijl(nb,na,la,2) = Gijl(nb,na,la,2) + d
            end if
         end do
      end do



*     **** add kinetic ***
      do i=1,nbasis
         la = l_ps(i)
         na = n_ps(i) - la
         power_f = 2*la
         do k=1,nrho
           f(k) = 0.5d0*(dphi_ae(k,i)*dphi_ae(k,i)
     >                  -dphi_ps(k,i)*dphi_ps(k,i))
     >       + 0.5d0*la*(la+1)*(phi_ae(k,i)*phi_ae(k,i) -
     >                          phi_ps(k,i)*phi_ps(k,i))/rho(k)**2
         end do
         d = log_integrate_def(power_f,f,0,rho,log_amesh,nrho)
         Gijl(na,na,la,3) = Gijl(na,na,la,3) + d

         do j=1,i-1
            lb = l_ps(j)
            nb = n_ps(j) - lb
            if (la.eq.lb) then
              do k=1,nrho
                f(k) = 0.5d0*(dphi_ae(k,i)*dphi_ae(k,j)
     >                       -dphi_ps(k,i)*dphi_ps(k,j))
     >            + 0.5d0*la*(la+1)*(phi_ae(k,i)*phi_ae(k,j) -
     >                               phi_ps(k,i)*phi_ps(k,j))/rho(k)**2
              end do
              d = log_integrate_def(power_f,f,0,rho,log_amesh,nrho)
              Gijl(na,nb,la,3) = Gijl(na,nb,la,3) + d
              Gijl(nb,na,la,3) = Gijl(nb,na,la,3) + d
            end if
         end do
      end do


*     **** add pseudo ***
      do i=1,nbasis
         la = l_ps(i)
         na = n_ps(i) - la
         power_f = 2*la + 2
         do k=1,nrho
           xerf = rho(k)/rlocal
           yerf = -zv*(1.0d0-util_erf(xerf))/rho(k) - zc/rho(k)
           f(k) = phi_ps(k,i)*phi_ps(k,i)*(yerf - v_ps(k))
         end do
         d = log_integrate_def(power_f,f,0,rho,log_amesh,nrho)
         Gijl(na,na,la,4) = Gijl(na,na,la,4) + d

         do j=1,i-1
            lb = l_ps(j)
            nb = n_ps(j) - lb
            if (la.eq.lb) then
              do k=1,nrho
                xerf = rho(k)/rlocal
                yerf = -zv*(1.0d0-util_erf(xerf))/rho(k) - zc/rho(k)
                f(k) = phi_ps(k,i)*phi_ps(k,j)*(yerf - v_ps(k))
              end do
              d = log_integrate_def(power_f,f,0,rho,log_amesh,nrho)
              Gijl(na,nb,la,4) = Gijl(na,nb,la,4) + d
              Gijl(nb,na,la,4) = Gijl(nb,na,la,4) + d
            end if
         end do
      end do

      l = nmax*nmax*(lmax+1)
      call dcopy(l,      Gijl(1,1,0,3),1,Gijl(1,1,0,1),1)  !** kinetic **
      call daxpy(l,1.0d0,Gijl(1,1,0,4),1,Gijl(1,1,0,1),1)  !** pseudo  **
c      call daxpy(l,1.0d0,Gijl(1,1,0,5),1,Gijl(1,1,0,1),1)
        

*     ************* compute ray fourier transforms *********************
      call integrate_paw_ray(version,rlocal,
     >                       nrho,log_amesh,zv,
     >                       nbasis,n_ps,l_ps,
     >                       v_ps,prj_ps,
     >                       rho,f,cs,sn,
     >                       nray,
     >                       G_ray,vl_ray,vnl_ray,
     >                       ierr)

*     **** filter the rays ****
      if (filter) then
         ecut = control_ecut()
         wcut = control_wcut()
         call kbpp_filter_ray(nray,G_ray,ecut,vl_ray)
         do i=1,nbasis
           call kbpp_filter_ray(nray,G_ray,wcut,vnl_ray(1,i,1))
         end do
      end if

*     **** setup cubic bsplines ****
      dG = G_ray(3)-G_ray(2)
      !yp1 = (vl_ray(3,1)-vl_ray(2,1))/dG
      !**** five point formula ***
      yp1 = ( -50.0d0*vl_ray(2,1)
     >       + 96.0d0*vl_ray(3,1)
     >       - 72.0d0*vl_ray(4,1)
     >       + 32.0d0*vl_ray(5,1)
     >       -  6.0d0*vl_ray(6,1))/(24.0d0*dG)
      call nwpw_spline(G_ray(2),vl_ray(2,1),nray-1,yp1,0.0d0,
     >                          vl_ray(2,2),tmp_ray)
      do i=1,nbasis
         call nwpw_spline(G_ray,vnl_ray(1,i,1),nray,0.0d0,0.0d0,
     >                          vnl_ray(1,i,2),tmp_ray)
      end do
*     ======================  Fourier transformation  ======================
      call dcopy(nfft3d,0.0d0,0,vl,1)
      call dcopy(nproj*nfft3d,0.0d0,0,vnl,1)
      task_count = -1
      DO 700 k3=1,nfft3
      DO 700 k2=1,nfft2
      DO 700 k1=1,(nfft1/2+1)
        task_count = task_count + 1
        if (mod(task_count,np).ne.taskid) go to 700

        Q=DSQRT(G(k1,k2,k3,1)**2
     >         +G(k1,k2,k3,2)**2
     >         +G(k1,k2,k3,3)**2)
        nx = (Q/dG) + 1.0d0

        if ((k1.eq.1).and.(k2.eq.1).and.(k3.eq.1)) go to 700
        
        GX=G(k1,k2,k3,1)/Q
        GY=G(k1,k2,k3,2)/Q
        GZ=G(k1,k2,k3,3)/Q

        lcount = 0
        do i=1,nbasis
          la=l_ps(i)

*         ::::::::::::::::::::::::::::::  h-wave  ::::::::::::::::::::::::::::::
*         ::::::::::::::::::::::::::::::  g-wave  ::::::::::::::::::::::::::::::
*         ::::::::::::::::::::::::::::::  f-wave  ::::::::::::::::::::::::::::::
          if (la.eq.3) then
          D = nwpw_splint(G_ray,vnl_ray(1,i,1),vnl_ray(1,i,2),nray,nx,Q)
           lcount = lcount+1
           vnl(k1,k2,k3,lcount)=D*GX*(4.0d0*GX*GX-3.0d0*(1.0d0-GZ*GZ))
     >                          /dsqrt(24.0d0)
           lcount = lcount+1
           vnl(k1,k2,k3,lcount)=D*GY*(3.0d0*(1.0d0-GZ*GZ)-4.0d0*GY*GY)
     >                          /dsqrt(24.0d0)
           lcount = lcount+1
           vnl(k1,k2,k3,lcount)=D*GZ*(GX*GX - GY*GY)
     >                          /2.0d0
           lcount = lcount+1
           vnl(k1,k2,k3,lcount)=D*GX*GY*GZ
           lcount = lcount+1
           vnl(k1,k2,k3,lcount)=D*GX*(5.0d0*GZ*GZ-1.0d0)
     >                          /dsqrt(40.0d0)
           lcount = lcount+1
           vnl(k1,k2,k3,lcount)=D*GY*(5.0d0*GZ*GZ-1.0d0)
     >                          /dsqrt(40.0d0)
           lcount = lcount+1
           vnl(k1,k2,k3,lcount)=D*GZ*(5.0d0*GZ*GZ-3.0d0)
     >                          /dsqrt(60.0d0)
          end if
*         ::::::::::::::::::::::::::::::  d-wave  ::::::::::::::::::::::::::::::
          if (la.eq.2) then
          D = nwpw_splint(G_ray,vnl_ray(1,i,1),vnl_ray(1,i,2),nray,nx,Q)
            lcount = lcount+1
            vnl(k1,k2,k3,lcount)=D*(3.0d0*GZ*GZ-1.0d0)
     >                            /(2.0d0*dsqrt(3.0d0))
            lcount = lcount+1
            vnl(k1,k2,k3,lcount)=D*GX*GY
            lcount = lcount+1
            vnl(k1,k2,k3,lcount)=D*GY*GZ
            lcount = lcount+1
            vnl(k1,k2,k3,lcount)=D*GZ*GX
            lcount = lcount+1
            vnl(k1,k2,k3,lcount)=D*(GX*GX-GY*GY)/(2.0d0)
          end if
*         ::::::::::::::::::::::::::::::  p-wave  ::::::::::::::::::::::::::::::
          if (la.eq.1) then
          P = nwpw_splint(G_ray,vnl_ray(1,i,1),vnl_ray(1,i,2),nray,nx,Q)
           lcount = lcount+1
           vnl(k1,k2,k3,lcount)=P*GX
           lcount = lcount+1
           vnl(k1,k2,k3,lcount)=P*GY
           lcount = lcount+1
           vnl(k1,k2,k3,lcount)=P*GZ
          end if
*         ::::::::::::::::::::::::::::::  s-wave  :::::::::::::::::::::::::::::::
          if (la.eq.0) then
          P = nwpw_splint(G_ray,vnl_ray(1,i,1),vnl_ray(1,i,2),nray,nx,Q)
            lcount = lcount+1
            vnl(k1,k2,k3,lcount)=P
          end if

        end do

*       ::::::::::::::::::::::::::::::  local  :::::::::::::::::::::::::::::::
  600   CONTINUE
        P = nwpw_splint(G_ray(2),vl_ray(2,1),vl_ray(2,2),nray-1,nx-1,Q)
        vl(k1,k2,k3)=P
       
  700 CONTINUE
      call Parallel_Vector_SumAll(nfft3d,vl)
      call Parallel_Vector_Sumall(nproj*nfft3d,vnl)

*     :::::::::::::::::::::::::::::::  G=0  ::::::::::::::::::::::::::::::::      
*     **** local potential ****
      vl(1,1,1)=vl_ray(1,1)

      do l=1,nproj
        vnl(1,1,1,l)=0.0d0
      end do


*     ***************************************    
*     **** define n_prj,l_prj, and m_prj ****
*     ***************************************
      lcount = 0
      do i=1,nbasis
        la=l_ps(i)
        na=n_ps(i) - la

        !::::::  f-wave  :::::::
        if (la.eq.3) then
          lcount = lcount+1
          n_prj(lcount) = na
          l_prj(lcount) = 3
          m_prj(lcount) = -3
    
          lcount = lcount+1
          n_prj(lcount) = na
          l_prj(lcount) = 3
          m_prj(lcount) = -2
           
          lcount = lcount+1
          n_prj(lcount) = na
          l_prj(lcount) = 3
          m_prj(lcount) = -1
     
          lcount = lcount+1
          n_prj(lcount) = na
          l_prj(lcount) = 3
          m_prj(lcount) = 0
           
          lcount = lcount+1
          n_prj(lcount) = na
          l_prj(lcount) = 3
          m_prj(lcount) = 1
           
          lcount = lcount+1
          n_prj(lcount) = na
          l_prj(lcount) = 3
          m_prj(lcount) = 2
     
          lcount = lcount+1
          n_prj(lcount) = na
          l_prj(lcount) = 3
          m_prj(lcount) = 3
        end if

        !::::  d-wave  ::::
        if (la.eq.2) then
          lcount = lcount+1
          n_prj(lcount) = na
          l_prj(lcount) = 2
          m_prj(lcount) = -2

          lcount = lcount+1
          n_prj(lcount) = na
          l_prj(lcount) = 2
          m_prj(lcount) = -1
          
          lcount = lcount+1
          n_prj(lcount) = na
          l_prj(lcount) = 2
          m_prj(lcount) = 0
          
          lcount = lcount+1
          n_prj(lcount) = na
          l_prj(lcount) = 2
          m_prj(lcount) = 1
          
          lcount = lcount+1
          n_prj(lcount) = na
          l_prj(lcount) = 2
          m_prj(lcount) = 2
        end if

        !::::  p-wave  ::::
        if (la.eq.1) then
          lcount = lcount+1
          n_prj(lcount) = na
          l_prj(lcount) = 1
          m_prj(lcount) = -1

          lcount = lcount+1
          n_prj(lcount) = na
          l_prj(lcount) = 1
          m_prj(lcount) = 0

          lcount = lcount+1
          n_prj(lcount) = na
          l_prj(lcount) = 1
          m_prj(lcount) = 1
        end if

        !::::  s-wave  ::::
        if (la.eq.0) then
          lcount = lcount+1
          n_prj(lcount) = na
          l_prj(lcount) = 0
          m_prj(lcount) = 0
          vnl(1,1,1,lcount)=vnl_ray(1,i,1) !*** j0!=0 at G=0***
        end if
      end do



      IERR=0
      RETURN
      END



