*
* $Id$
*

      subroutine integrate_kbpp_band_local_new(version,
     >                            nrho,drho,lmax,locp,zv,
     >                            vp,wp,rho,f,cs,sn,
     >                            nfft1,nfft2,nfft3,lmmax,
     >                            G,vl,
     >                            n_prj,l_prj,m_prj,vnlnrm,
     >                            semicore,rho_sc_r,rho_sc_k,
     >                            nray,G_ray,vl_ray,vnl_ray,
     >                            rho_sc_k_ray,tmp_ray,
     >                            filter,
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
      double precision vl(nfft1,nfft2,nfft3)
      integer          n_prj(lmmax),l_prj(lmmax),m_prj(lmmax)
      double precision vnlnrm(0:lmax)


      logical semicore
      double precision rho_sc_r(nrho,2)
      double precision rho_sc_k(nfft1,nfft2,nfft3,4)

      integer nray
      double precision G_ray(nray)
      double precision vl_ray(nray,2)
      double precision vnl_ray(nray,0:lmax,2)
      double precision rho_sc_k_ray(nray,2,2)
      double precision tmp_ray(nray)
      logical filter

      integer ierr

#include "bafdecls.fh"

*     *** local variables ****
      integer np,taskid,MASTER
      parameter (MASTER=0)

      integer lcount,task_count,nfft3d
      integer k1,k2,k3,i,l,nx
      double precision ecut,wcut
      double precision p
      double precision gx,gy,gz,a,q,dG,yp1

*     **** external functions ****
      double precision dsum,simp,control_ecut,control_wcut,nwpw_splint
      external         dsum,simp,control_ecut,control_wcut,nwpw_splint

      call Parallel_np(np)
      call Parallel_taskid(taskid)


      nfft3d = (nfft1)*nfft2*nfft3

      IF((NRHO/2)*2.EQ.NRHO) THEN
        WRITE(*,*)"local pseudopotential not computed: nrho is not odd"
        IERR=2
        RETURN
      ENDIF


*::::::::::::::::::  Define non-local pseudopotential  ::::::::::::::::
      do l=0,lmax
        if (l.ne.locp) then
          do i=1,nrho
            vp(i,l)=vp(i,l)-vp(i,locp)
          end do
        end if
      end do

*:::::::::::::::::::::  Normarization constants  ::::::::::::::::::::::
      lcount = 0
      do l=0,lmax
        if (l.ne.locp) then
          do i=1,nrho
            f(i)=vp(i,l)*wp(i,l)**2
          end do
          a=simp(nrho,f,drho)
          vnlnrm(l) = 1.0d0/a
        else
        vnlnrm(l) = 0.0d0
        end if
      end do

*     ************* compute ray fourier transforms *********************
      call integrate_kbpp_band_ray(version,
     >                            nrho,drho,lmax,locp,zv,
     >                            vp,wp,rho,f,cs,sn,
     >                            nray,
     >                            G_ray,vl_ray,vnl_ray,
     >                            semicore,rho_sc_r,rho_sc_k_ray,
     >                            ierr)

*     **** filter the rays ****
      if (filter) then
         ecut = control_ecut()
         wcut = control_wcut()
         call kbpp_band_filter_ray(nray,G_ray,ecut,vl_ray)
         do l=0,lmax
            if (l.ne.locp)
     >       call kbpp_band_filter_ray(nray,G_ray,wcut,vnl_ray(1,l,1))
         end do
         if (semicore) then
         call kbpp_band_filter_ray(nray,G_ray,ecut,rho_sc_k_ray(1,1,1))
         call kbpp_band_filter_ray(nray,G_ray,ecut,rho_sc_k_ray(1,2,1))
         end if
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
      do l=0,lmax
         if (l.ne.locp)
     >      call nwpw_spline(G_ray,vnl_ray(1,l,1),nray,0.0d0,0.0d0,
     >                             vnl_ray(1,l,2),tmp_ray)
      end do
      if (semicore) then
         call nwpw_spline(G_ray,rho_sc_k_ray(1,1,1),nray,0.0d0,0.0d0,
     >                          rho_sc_k_ray(1,1,2),tmp_ray)
         call nwpw_spline(G_ray,rho_sc_k_ray(1,2,1),nray,0.0d0,0.0d0,
     >                          rho_sc_k_ray(1,2,2),tmp_ray)
      end if


*======================  Fourier transformation  ======================
      call dcopy(nfft3d,0.0d0,0,vl,1)
      call dcopy(4*nfft3d,0.0d0,0,rho_sc_k,1)
      task_count = -1
      DO 700 k3=1,nfft3
      DO 700 k2=1,nfft2
      DO 700 k1=1,nfft1
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
        do i=1,nrho
          cs(i)=DCOS(Q*rho(i))
          sn(i)=DSIN(Q*rho(i))
        end do

*::::::::::::::::::::::::::::::  local  :::::::::::::::::::::::::::::::
  600   CONTINUE

        vl(k1,k2,k3)= nwpw_splint(G_ray(2),vl_ray(2,1),
     >                                     vl_ray(2,2),nray-1,nx-1,Q)

 
*::::::::::::::::::::: semicore density :::::::::::::::::::::::::::::::
        if (semicore) then
           p = nwpw_splint(G_ray,rho_sc_k_ray(1,1,1),
     >                           rho_sc_k_ray(1,1,2),nray,nx,Q)
           rho_sc_k(k1,k2,k3,1) = p

           p = nwpw_splint(G_ray,rho_sc_k_ray(1,2,1),
     >                           rho_sc_k_ray(1,2,2),nray,nx,Q)
           rho_sc_k(k1,k2,k3,2)=p*GX
           rho_sc_k(k1,k2,k3,3)=p*GY
           rho_sc_k(k1,k2,k3,4)=p*GZ

        end if
    
  700 CONTINUE
      call Parallel_Vector_SumAll(4*nfft3d,rho_sc_k)
      call Parallel_Vector_SumAll(nfft3d,vl)

*:::::::::::::::::::::::::::::::  G=0  ::::::::::::::::::::::::::::::::      

*     **** local potential ****
      vl(1,1,1)=vl_ray(1,1)

*     **** semicore density ****
      if (semicore) then
         rho_sc_k(1,1,1,1) = rho_sc_k_ray(1,1,1)
         rho_sc_k(1,1,1,2) = 0.0d0
         rho_sc_k(1,1,1,3) = 0.0d0
         rho_sc_k(1,1,1,4) = 0.0d0
      end if

*     ********************************
*     **** define n_prj and l_prj ****
*     ********************************
      lcount = lmmax+1
      GO TO (950,940,930,920), lmax+1

        !::::::  f-wave  :::::::
  920   CONTINUE
        if (locp.ne.3) then
          lcount = lcount-1
          n_prj(lcount) = 1
          l_prj(lcount) = 3
          m_prj(lcount) = -3

          lcount = lcount-1
          n_prj(lcount) = 1
          l_prj(lcount) = 3
          m_prj(lcount) = -2

          lcount = lcount-1
          n_prj(lcount) = 1
          l_prj(lcount) = 3
          m_prj(lcount) = -1

          lcount = lcount-1
          n_prj(lcount) = 1
          l_prj(lcount) = 3
          m_prj(lcount) = 0

          lcount = lcount-1
          n_prj(lcount) = 1
          l_prj(lcount) = 3
          m_prj(lcount) = 1

          lcount = lcount-1
          n_prj(lcount) = 1
          l_prj(lcount) = 3
          m_prj(lcount) = 2

          lcount = lcount-1
          n_prj(lcount) = 1
          l_prj(lcount) = 3
          m_prj(lcount) = 3
        end if


        !::::  d-wave  ::::
  930   CONTINUE
        if (locp.ne.2) then
          lcount = lcount-1
          n_prj(lcount) = 1
          l_prj(lcount) = 2
          m_prj(lcount) = -2

          lcount = lcount-1
          n_prj(lcount) = 1
          l_prj(lcount) = 2
          m_prj(lcount) = -1

          lcount = lcount-1
          n_prj(lcount) = 1
          l_prj(lcount) = 2
          m_prj(lcount) = 0

          lcount = lcount-1
          n_prj(lcount) = 1
          l_prj(lcount) = 2
          m_prj(lcount) = 1

          lcount = lcount-1
          n_prj(lcount) = 1
          l_prj(lcount) = 2
          m_prj(lcount) = 2
        end if


        !::::  p-wave  ::::
  940   CONTINUE
        if (locp.ne.1) then
          lcount = lcount-1
          n_prj(lcount) = 1
          l_prj(lcount) = 1
          m_prj(lcount) = -1

          lcount = lcount-1
          n_prj(lcount) = 1
          l_prj(lcount) = 1
          m_prj(lcount) = 0

          lcount = lcount-1
          n_prj(lcount) = 1
          l_prj(lcount) = 1
          m_prj(lcount) = 1
        end if


        !::::  s-wave  ::::
  950   CONTINUE
        if (locp.ne.0) then
          lcount = lcount-1
          n_prj(lcount) = 1
          l_prj(lcount) = 0
          m_prj(lcount) = 0
        end if

      IERR=0
      RETURN
      END

