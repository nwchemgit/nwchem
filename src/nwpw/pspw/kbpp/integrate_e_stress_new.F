*
* $Id$
*

*     **********************************
*     *                                *
*     *      integrate_e_stress_new    *
*     *                                *
*     **********************************

      subroutine integrate_e_stress_new(version,rlocal,
     >                            nrho,drho,lmax,locp,nmax,
     >                            n_extra,n_expansion,zv,
     >                            vp,wp,rho,f,cs,sn,
     >                            nfft3d,nprj,
     >                            G,dvl,dvnl,
     >                            semicore,rho_sc_r,rho_sc_k,
     >                            nray,G_ray,dvl_ray,dvnl_ray,
     >                            rho_sc_k_ray,tmp_ray,
     >                            filter,
     >                            ierr)
      implicit none
      integer          version
      double precision rlocal
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

      integer nfft3d,nprj
      double precision G(nfft3d,3)
      double precision dvl(nfft3d)
      double precision dvnl(nfft3d,3,nprj)

      logical semicore
      double precision rho_sc_r(nrho,2)
      double precision rho_sc_k(nfft3d,4)

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
      integer np_i,np_j,taskid_i,taskid_j,countj
      parameter (MASTER=0)

      integer lcount,pzero,zero,n,nb
      integer k1,i,l,nx
      double precision p,pp
      double precision gx,gy,gz,q,d,dd
      double precision duxdGx,duxdGy,duxdGz
      double precision duydGx,duydGy,duydGz
      double precision duzdGx,duzdGy,duzdGz
      double precision sumx,sumy,sumz,ecut,wcut,dG
      double precision T,dTdux,dTduy,dTduz
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

      call Parallel2d_np_i(np_i)
      call Parallel2d_np_j(np_j)
      call Parallel2d_taskid_i(taskid_i)
      call Parallel2d_taskid_j(taskid_j)

      if (version.ne.3) then
         call errquit('integrate_stress_new - unit cell is aperiodic',0,
     &       INPUT_ERR)
      end if
      if (lmax.gt.3) then
         call errquit('integrate_stress_new - lmax > f',0,
     >       INPUT_ERR)
      end if
      if ((nrho/2)*2.eq.nrho) then
        call errquit('integrate_stress_new - psp grid not odd',0,
     >       INPUT_ERR)
      end if


*::::::::::::::::::  Define non-local pseudopotential  ::::::::::::::::
      do l=0,lmax
        if (l.ne.locp) then
          do I=1,nrho
            vp(i,l)=vp(i,l)-vp(i,locp)
          end do
        end if
      end do

*     ************* compute ray fourier transforms *********************
      call integrate_e_stress_ray(version,rlocal,
     >                            nrho,drho,lmax,locp,nmax,
     >                            n_extra,n_expansion,zv,
     >                            vp,wp,rho,f,cs,sn,
     >                            nray,
     >                            G_ray,dvl_ray,dvnl_ray,
     >                            semicore,rho_sc_r,rho_sc_k_ray,
     >                            ierr)
*     **** filter the rays ****
      if (filter) then
         ecut = control_ecut()
         wcut = control_wcut()
         call kbpp_filter_ray(nray,G_ray,ecut,dvl_ray)
         do l=0,lmax
            if (l.ne.locp) then
              do n=1,n_expansion(l)
                 call kbpp_filter_ray(nray,G_ray,wcut,
     >                                dvnl_ray(1,1,indx(n,l),1))
                 call kbpp_filter_ray(nray,G_ray,wcut,
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
              call nwpw_spline(G_ray,dvnl_ray(1,1,indx(n,l),1),nray,
     >                   0.0d0,0.0d0,dvnl_ray(1,1,indx(n,l),2),tmp_ray)
              call nwpw_spline(G_ray,dvnl_ray(1,2,indx(n,l),1),nray,
     >                   0.0d0,0.0d0,dvnl_ray(1,2,indx(n,l),2),tmp_ray)
            end do
         end if
      end do
      if (semicore) then
         call nwpw_spline(G_ray,rho_sc_k_ray(1,1),nray,0.0d0,0.0d0,
     >                          rho_sc_k_ray(1,2),tmp_ray)
      end if

*======================  Fourier transformation  ======================
      call dcopy(nfft3d,0.0d0,0,dvl,1)
      call dcopy(3*nprj*nfft3d,0.0d0,0,dvnl,1)
      call dcopy(nfft3d,0.0d0,0,rho_sc_k,1)

*     ***** find the G==0 point in the lattice *****
      call D3dB_ijktoindexp(1,1,1,1,zero,pzero)
      
      countj = -1
      DO 700 k1=1,nfft3d

        countj = mod(countj+1,np_j)
        if (countj.ne.taskid_j) go to 700
        if ((pzero.eq.taskid_i).and.(k1.eq.zero)) go to 700

        q=dsqrt(G(k1,1)**2
     >         +G(k1,2)**2
     >         +G(k1,3)**2)
        nx = (q/dG) + 1.0d0

        
        gx=G(k1,1)/q
        gy=G(k1,2)/q
        gz=G(k1,3)/q
        DO i=1,nrho
          cs(i)=dcos(q*rho(i))
          sn(i)=dsin(q*rho(i))
        END DO

*       **** calculate du_r/dG_s ****
        duxdGx = 1.0d0/q -gx*gx/q
        duxdGy = -gx*gy/q
        duxdGz = -gx*gz/q

        duydGx = -gy*gx/q
        duydGy = 1.0d0/q - gy*gy/q
        duydGz = -gy*gz/q

        duzdGx = -gz*gx/q
        duzdGy = -gz*gy/q
        duzdGz = 1.0d0/q - gz*gz/q

        lcount = nprj+1
        GO TO (500,400,300,200), lmax+1

*::::::::::::::::::::::::::::::  f-wave  ::::::::::::::::::::::::::::::
  200   CONTINUE
        if (locp.ne.3) then
           do n=1,n_expansion(3)
           D  = nwpw_splint(G_ray,dvnl_ray(1,1,indx(n,3),1),
     >                            dvnl_ray(1,1,indx(n,3),2),nray,nx,Q)
           DD = nwpw_splint(G_ray,dvnl_ray(1,2,indx(n,3),1),
     >                            dvnl_ray(1,2,indx(n,3),2),nray,nx,Q)

           lcount = lcount-1
           T = gy*(3.0d0*(1.0d0-gz*gz)-4.0d0*gy*gy)/dsqrt(24.0d0)
           dTdux = 0.0d0
           dTduy = (3.0d0*(1.0d0-gz*gz)-12.0d0*gy*gy)/dsqrt(24.0d0)
           dTduz = -6.0d0*gy*gz/dsqrt(24.0d0)
           sumx = dTdux*duxdGx + dTduy*duydGx + dTduz*duzdGx
           sumy = dTdux*duxdGy + dTduy*duydGy + dTduz*duzdGy
           sumz = dTdux*duxdGz + dTduy*duydGz + dTduz*duzdGz
           dvnl(k1,1,lcount)=DD*T*gx + D*sumx
           dvnl(k1,2,lcount)=DD*T*gy + D*sumy
           dvnl(k1,3,lcount)=DD*T*gz + D*sumz

           lcount = lcount-1
           T =gx*gy*gz 
           dTdux = gy*gz
           dTduy = gx*gz
           dTduz = gx*gy
           sumx = dTdux*duxdGx + dTduy*duydGx + dTduz*duzdGx
           sumy = dTdux*duxdGy + dTduy*duydGy + dTduz*duzdGy
           sumz = dTdux*duxdGz + dTduy*duydGz + dTduz*duzdGz
           dvnl(k1,1,lcount)=DD*T*gx + D*sumx
           dvnl(k1,2,lcount)=DD*T*gy + D*sumy
           dvnl(k1,3,lcount)=DD*T*gz + D*sumz

           lcount = lcount-1
           T = gy*(5.0d0*gz*gz-1.0d0)/dsqrt(40.0d0)
           dTdux = 0.0d0
           dTduy =(5.0d0*gz*gz-1.0d0)/dsqrt(40.0d0)
           dTduz =10.0d0*gy*gz/dsqrt(40.0d0)
           sumx = dTdux*duxdGx + dTduy*duydGx + dTduz*duzdGx
           sumy = dTdux*duxdGy + dTduy*duydGy + dTduz*duzdGy
           sumz = dTdux*duxdGz + dTduy*duydGz + dTduz*duzdGz
           dvnl(k1,1,lcount)=DD*T*gx + D*sumx
           dvnl(k1,2,lcount)=DD*T*gy + D*sumy
           dvnl(k1,3,lcount)=DD*T*gz + D*sumz

           lcount = lcount-1
           T =gz*(5.0d0*gz*gz-3.0d0)/dsqrt(60.0d0)
           dTdux = 0.0d0
           dTduy = 0.0d0
           dTduz =(15.0d0*gz*gz -3.0d0)/dsqrt(60.0d0)
           sumx = dTdux*duxdGx + dTduy*duydGx + dTduz*duzdGx
           sumy = dTdux*duxdGy + dTduy*duydGy + dTduz*duzdGy
           sumz = dTdux*duxdGz + dTduy*duydGz + dTduz*duzdGz
           dvnl(k1,1,lcount)=DD*T*gx + D*sumx
           dvnl(k1,2,lcount)=DD*T*gy + D*sumy
           dvnl(k1,3,lcount)=DD*T*gz + D*sumz

           lcount = lcount-1
           T =  gx*(5.0d0*gz*gz-1.0d0)/dsqrt(40.0d0) 
           dTdux = (5.0d0*gz*gz-1.0d0)/dsqrt(40.0d0)
           dTduy = 0.0d0
           dTduz = 10.0d0*gx*gz/dsqrt(40.0d0)
           sumx = dTdux*duxdGx + dTduy*duydGx + dTduz*duzdGx
           sumy = dTdux*duxdGy + dTduy*duydGy + dTduz*duzdGy
           sumz = dTdux*duxdGz + dTduy*duydGz + dTduz*duzdGz
           dvnl(k1,1,lcount)=DD*T*gx + D*sumx
           dvnl(k1,2,lcount)=DD*T*gy + D*sumy
           dvnl(k1,3,lcount)=DD*T*gz + D*sumz

           lcount = lcount-1
           T =gz*(gx*gx - gy*gy)/2.0d0 
           dTdux =  gx*gz
           dTduy = -gy*gz
           dTduz = (gx*gx-gy*gy)/2.0d0
           sumx = dTdux*duxdGx + dTduy*duydGx + dTduz*duzdGx
           sumy = dTdux*duxdGy + dTduy*duydGy + dTduz*duzdGy
           sumz = dTdux*duxdGz + dTduy*duydGz + dTduz*duzdGz
           dvnl(k1,1,lcount)=DD*T*gx + D*sumx
           dvnl(k1,2,lcount)=DD*T*gy + D*sumy
           dvnl(k1,3,lcount)=DD*T*gz + D*sumz

           lcount = lcount-1
           T = gx*(4.0d0*gx*gx - 3.0d0*(1.0d0-gz*gz))/dsqrt(24.0d0)
           dTdux = (12.0d0*gx*gx-3.0d0*(1.0d0-gz*gz))/dsqrt(24.0d0)
           dTduy = 0.0d0
           dTduz = 6.0d0*gx*gz/dsqrt(24.0d0)
           sumx = dTdux*duxdGx + dTduy*duydGx + dTduz*duzdGx
           sumy = dTdux*duxdGy + dTduy*duydGy + dTduz*duzdGy
           sumz = dTdux*duxdGz + dTduy*duydGz + dTduz*duzdGz
           dvnl(k1,1,lcount)=DD*T*gx + D*sumx
           dvnl(k1,2,lcount)=DD*T*gy + D*sumy
           dvnl(k1,3,lcount)=DD*T*gz + D*sumz
           end do
        end if
*::::::::::::::::::::::::::::::  d-wave  ::::::::::::::::::::::::::::::
  300   CONTINUE
        if (locp.ne.2) then
          do n=1,n_expansion(2)
           D  = nwpw_splint(G_ray,dvnl_ray(1,1,indx(n,2),1),
     >                            dvnl_ray(1,1,indx(n,2),2),nray,nx,Q)
           DD = nwpw_splint(G_ray,dvnl_ray(1,2,indx(n,2),1),
     >                            dvnl_ray(1,2,indx(n,2),2),nray,nx,Q)

          lcount = lcount-1
          T = gx*gy
          dTdux = gy
          dTduy = gx
          dTduz = 0.0d0
          sumx = dTdux*duxdGx + dTduy*duydGx + dTduz*duzdGx
          sumy = dTdux*duxdGy + dTduy*duydGy + dTduz*duzdGy
          sumz = dTdux*duxdGz + dTduy*duydGz + dTduz*duzdGz
          dvnl(k1,1,lcount)=DD*T*gx + D*sumx
          dvnl(k1,2,lcount)=DD*T*gy + D*sumy
          dvnl(k1,3,lcount)=DD*T*gz + D*sumz

          lcount = lcount-1
          T = gy*gz
          dTdux = 0.0d0
          dTduy = gz
          dTduz = gy
          sumx = dTdux*duxdGx + dTduy*duydGx + dTduz*duzdGx
          sumy = dTdux*duxdGy + dTduy*duydGy + dTduz*duzdGy
          sumz = dTdux*duxdGz + dTduy*duydGz + dTduz*duzdGz
          dvnl(k1,1,lcount)=DD*T*gx + D*sumx
          dvnl(k1,2,lcount)=DD*T*gy + D*sumy
          dvnl(k1,3,lcount)=DD*T*gz + D*sumz

          lcount = lcount-1
          T = (3.0d0*gz*gz-1.0d0)/(2.0d0*dsqrt(3.0d0)) 
          dTdux = 0.0d0
          dTduy = 0.0d0
          dTduz = 6.0d0*gz/(2.0d0*dsqrt(3.0d0))
          sumx = dTdux*duxdGx + dTduy*duydGx + dTduz*duzdGx
          sumy = dTdux*duxdGy + dTduy*duydGy + dTduz*duzdGy
          sumz = dTdux*duxdGz + dTduy*duydGz + dTduz*duzdGz
          dvnl(k1,1,lcount)=DD*T*gx + D*sumx
          dvnl(k1,2,lcount)=DD*T*gy + D*sumy
          dvnl(k1,3,lcount)=DD*T*gz + D*sumz

          lcount = lcount-1
          T = gz*gx
          dTdux = gz
          dTduy = 0.0d0
          dTduz = gx
          sumx = dTdux*duxdGx + dTduy*duydGx + dTduz*duzdGx
          sumy = dTdux*duxdGy + dTduy*duydGy + dTduz*duzdGy
          sumz = dTdux*duxdGz + dTduy*duydGz + dTduz*duzdGz
          dvnl(k1,1,lcount)=DD*T*gx + D*sumx
          dvnl(k1,2,lcount)=DD*T*gy + D*sumy
          dvnl(k1,3,lcount)=DD*T*gz + D*sumz

          lcount = lcount-1
          T = (gx*gx-gy*gy)/2.0d0
          dTdux = gx
          dTduy = -gy
          dTduz = 0.0d0
          sumx = dTdux*duxdGx + dTduy*duydGx + dTduz*duzdGx
          sumy = dTdux*duxdGy + dTduy*duydGy + dTduz*duzdGy
          sumz = dTdux*duxdGz + dTduy*duydGz + dTduz*duzdGz
          dvnl(k1,1,lcount)=DD*T*gx + D*sumx
          dvnl(k1,2,lcount)=DD*T*gy + D*sumy
          dvnl(k1,3,lcount)=DD*T*gz + D*sumz
          end do
        end if
*::::::::::::::::::::::::::::::  p-wave  ::::::::::::::::::::::::::::::
  400   CONTINUE
        if (locp.ne.1) then
          do n=1,n_expansion(1)
           P  = nwpw_splint(G_ray,dvnl_ray(1,1,indx(n,1),1),
     >                            dvnl_ray(1,1,indx(n,1),2),nray,nx,Q)
           PP = nwpw_splint(G_ray,dvnl_ray(1,2,indx(n,1),1),
     >                            dvnl_ray(1,2,indx(n,1),2),nray,nx,Q)
           lcount = lcount-1
           T = gy
           dTdux = 0.0d0
           dTduy = 1.0d0
           dTduz = 0.0d0
           sumx = dTdux*duxdGx + dTduy*duydGx + dTduz*duzdGx
           sumy = dTdux*duxdGy + dTduy*duydGy + dTduz*duzdGy
           sumz = dTdux*duxdGz + dTduy*duydGz + dTduz*duzdGz
           dvnl(k1,1,lcount)= PP*T*gx + P*sumx
           dvnl(k1,2,lcount)= PP*T*gy + P*sumy
           dvnl(k1,3,lcount)= PP*T*gz + P*sumz

           lcount = lcount-1
           T = gz
           dTdux = 0.0d0
           dTduy = 0.0d0
           dTduz = 1.0d0
           sumx = dTdux*duxdGx + dTduy*duydGx + dTduz*duzdGx
           sumy = dTdux*duxdGy + dTduy*duydGy + dTduz*duzdGy
           sumz = dTdux*duxdGz + dTduy*duydGz + dTduz*duzdGz
           dvnl(k1,1,lcount)= PP*T*gx + P*sumx
           dvnl(k1,2,lcount)= PP*T*gy + P*sumy
           dvnl(k1,3,lcount)= PP*T*gz + P*sumz

           lcount = lcount-1
           T = gx
           dTdux = 1.0d0
           dTduy = 0.0d0
           dTduz = 0.0d0
           sumx = dTdux*duxdGx + dTduy*duydGx + dTduz*duzdGx
           sumy = dTdux*duxdGy + dTduy*duydGy + dTduz*duzdGy
           sumz = dTdux*duxdGz + dTduy*duydGz + dTduz*duzdGz
           dvnl(k1,1,lcount)= PP*T*gx + P*sumx
           dvnl(k1,2,lcount)= PP*T*gy + P*sumy
           dvnl(k1,3,lcount)= PP*T*gz + P*sumz
          end do
        end if
*::::::::::::::::::::::::::::::  s-wave  :::::::::::::::::::::::::::::::
  500   CONTINUE
        if (locp.ne.0) then
          do n=1,n_expansion(0)
          P  = nwpw_splint(G_ray,dvnl_ray(1,1,indx(n,0),1),
     >                           dvnl_ray(1,1,indx(n,0),2),nray,nx,Q)
          lcount = lcount-1
          dvnl(k1,1,lcount) = P *gx
          dvnl(k1,2,lcount) = P *gy
          dvnl(k1,3,lcount) = P *gz
          end do
        end if
*::::::::::::::::::::::::::::::  local  :::::::::::::::::::::::::::::::
  600   CONTINUE
        dvl(k1)= nwpw_splint(G_ray(1),dvl_ray(1,1),
     >                                      dvl_ray(1,2),nray,nx,Q)
*::::::::::::::::::::: semicore density :::::::::::::::::::::::::::::::
        if (semicore) then
           rho_sc_k(k1,1) = nwpw_splint(G_ray,
     >                                rho_sc_k_ray(1,1),
     >                                rho_sc_k_ray(1,2),nray,nx,Q)
        end if
  700 CONTINUE
      call D1dB_Vector_SumAll(nfft3d,dvl)
      call D1dB_Vector_Sumall(3*nprj*nfft3d,dvnl)
      call D1dB_Vector_SumAll(nfft3d,rho_sc_k)


*:::::::::::::::::::::::::::::::  G=0  ::::::::::::::::::::::::::::::::      

      ierr=0
      return
      end



