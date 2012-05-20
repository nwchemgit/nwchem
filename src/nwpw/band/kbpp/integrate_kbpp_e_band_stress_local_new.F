*
* $Id$
*

*     *****************************************
*     *                                        *
*     *  integrate_kbpp_e_band_stress_local_new  *
*     *                                        *
*     ******************************************
      subroutine integrate_kbpp_e_band_stress_local_new(version,
     >                            nrho,drho,lmax,locp,nmax,
     >                            n_extra,n_expansion,zv,
     >                            vp,wp,rho,f,cs,sn,
     >                            nfft1,nfft2,nfft3,nprj,
     >                            G,dvl,
     >                            semicore,rho_sc_r,rho_sc_k,
     >                            nray,G_ray,dvl_ray,dvnl_ray,
     >                            rho_sc_k_ray,tmp_ray,
     >                            filter,
     >                            ierr)
      implicit none
      integer          version
      integer          nrho
      double precision drho
      integer          lmax
      integer          locp
      integer          nmax
      integer          n_extra,n_expansion(0:lmax)
      double precision zv
      double precision vp(nrho,0:lmax)
      double precision wp(nrho,0:(lmax+n_extra))
      double precision rho(nrho)
      double precision f(nrho)
      double precision cs(nrho)
      double precision sn(nrho)

      integer nfft1,nfft2,nfft3,nprj
      double precision G(nfft1,nfft2,nfft3,3)
      double precision dvl(nfft1,nfft2,nfft3)

      logical semicore
      double precision rho_sc_r(nrho,2)
      double precision rho_sc_k(nfft1,nfft2,nfft3,4)

      integer nray
      double precision G_ray(nray)
      double precision dvl_ray(nray,2)
      double precision dvnl_ray(nray,2,0:(lmax+n_extra),2)
      double precision rho_sc_k_ray(nray,2)
      double precision tmp_ray(nray)
      logical filter
      integer ierr

#include "errquit.fh"

*     *** local variables ****
      integer np,taskid,MASTER
      parameter (MASTER=0)

      integer task_count,nfft3d
      integer k1,k2,k3,i,l,nx,n,nb
      double precision q,dG,ecut,wcut
      integer indx(5,0:3)

*     **** external functions ****
      double precision dsum,simp,control_ecut,control_wcut,nwpw_splint
      external         dsum,simp,control_ecut,control_wcut,nwpw_splint

*     **** set up indx(n,l) --> to wp ****
      nb = lmax+1
      do l=0,lmax
         indx(1,l) = l
         do n=2,n_expansion(l)
            indx(n,l) = nb
            nb = nb+1
         end do
      end do

      call Parallel_np(np)
      call Parallel_taskid(taskid)

      if (lmax.gt.3) then
         call errquit('integrate_kbpp_band_stress_local_new - lmax > f',
     >       0,INPUT_ERR)
      end if
      if ((nrho/2)*2.eq.nrho) then
        call errquit(
     >  'integrate_kbpp_band_stress_local_new - psp grid not odd',
     >       0,INPUT_ERR)
      end if

      nfft3d = (nfft1)*nfft2*nfft3

*::::::::::::::::::  Define non-local pseudopotential  ::::::::::::::::
      do l=0,lmax
        if (l.ne.locp) then
          do I=1,nrho
            vp(i,l)=vp(i,l)-vp(i,locp)
          end do
        end if
      end do

*     ************* compute ray fourier transforms *********************
      call integrate_kbpp_e_band_stress_ray(version,
     >                            nrho,drho,lmax,locp,nmax,
     <                            n_extra,n_expansion,zv,
     >                            vp,wp,rho,f,cs,sn,
     >                            nray,
     >                            G_ray,dvl_ray,dvnl_ray,
     >                            semicore,rho_sc_r,rho_sc_k_ray,
     >                            ierr)
*     **** filter the rays ****
      if (filter) then
         ecut = control_ecut()
         wcut = control_wcut()
         call kbpp_band_filter_ray(nray,G_ray,ecut,dvl_ray)
         do l=0,lmax
            if (l.ne.locp) then
               do n=1,n_expansion(l)
                  call kbpp_band_filter_ray(nray,G_ray,wcut,
     >                                dvnl_ray(1,1,indx(n,l),1))
                  call kbpp_band_filter_ray(nray,G_ray,wcut,
     >                                dvnl_ray(1,2,indx(n,l),1))
               end do
            end if
         end do
         if (semicore) then
           call kbpp_filter_ray(nray,G_ray,ecut,rho_sc_k_ray)
         end if
      end if

*     **** setup cubic bsplines ****
      dG = G_ray(3)-G_ray(2)
      call nwpw_spline(G_ray,dvl_ray(1,1),nray,0.0d0,0.0d0,
     >                       dvl_ray(1,2),tmp_ray)
      do l=0,lmax
         if (l.ne.locp) then
            do n=1,n_expansion(l)
               call nwpw_spline(G_ray,dvnl_ray(1,1,indx(n,l),1),
     >                             nray,0.0d0,0.0d0,
     >                             dvnl_ray(1,1,indx(n,l),2),tmp_ray)
               call nwpw_spline(G_ray,dvnl_ray(1,2,indx(n,l),1),
     >                             nray,0.0d0,0.0d0,
     >                             dvnl_ray(1,2,indx(n,l),2),tmp_ray)
            end do
         end if
      end do
      if (semicore) then
         call nwpw_spline(G_ray,rho_sc_k_ray(1,1),nray,0.0d0,0.0d0,
     >                          rho_sc_k_ray(1,2),tmp_ray)
      end if

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
        nx = (q/dG) + 1.0d0

        if ((k1.eq.1).and.(k2.eq.1).and.(k3.eq.1)) go to 700

        do i=1,nrho
          cs(i)=dcos(q*rho(i))
          sn(i)=dsin(q*rho(i))
        end do

*::::::::::::::::::::::::::::::  local  :::::::::::::::::::::::::::::::
  600   CONTINUE

        dvl(k1,k2,k3)= nwpw_splint(G_ray(1),dvl_ray(1,1),
     >                                      dvl_ray(1,2),nray,nx,Q)
*::::::::::::::::::::: semicore density :::::::::::::::::::::::::::::::
        if (semicore) then
           rho_sc_k(k1,k2,k3,1) = nwpw_splint(G_ray,
     >                                rho_sc_k_ray(1,1),
     >                                rho_sc_k_ray(1,2),nray,nx,Q)
        end if

  700 CONTINUE
      call Parallel_Vector_SumAll(4*nfft3d,rho_sc_k)
      call Parallel_Vector_SumAll(nfft3d,dvl)

*:::::::::::::::::::::::::::::::  G=0  ::::::::::::::::::::::::::::::::      
      dvl(1,1,1)= 0.0d0
      rho_sc_k(1,1,1,1) = 0.0d0

      ierr=0
      return
      end

