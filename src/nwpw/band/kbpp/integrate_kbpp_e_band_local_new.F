*
* $Id$
*
      subroutine integrate_kbpp_e_band_local_new(version,
     >                            nrho,drho,lmax,locp,nmax,
     >                            n_extra,n_expansion,zv,
     >                            vp,wp,rho,f,cs,sn,
     >                            nfft1,nfft2,nfft3,nprj,
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
      double precision vl(nfft1,nfft2,nfft3)
      integer          n_prj(nprj),l_prj(nprj),m_prj(nprj)
      double precision vnlnrm(nmax,nmax,0:lmax)


      logical semicore
      double precision rho_sc_r(nrho,2)
      double precision rho_sc_k(nfft1,nfft2,nfft3,4)

      integer nray
      double precision G_ray(nray)
      double precision vl_ray(nray,2)
      double precision vnl_ray(nray,0:(lmax+n_extra),2)
      double precision rho_sc_k_ray(nray,2,2)
      double precision tmp_ray(nray)
      logical filter

      integer ierr

#include "bafdecls.fh"

*     *** local variables ****
      integer np,taskid,MASTER
      parameter (MASTER=0)

      integer lcount,task_count,nfft3d
      integer k1,k2,k3,i,l,nx,n,nb,n1,n2
      double precision ecut,wcut
      double precision p,d
      double precision gx,gy,gz,a,q,dG,yp1
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

      nfft3d = (nfft1)*nfft2*nfft3

      if ((nrho/2)*2.eq.nrho) then
        write(*,*)"local pseudopotential not computed: nrho is not odd"
        ierr=2
        return
      end if


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
        if (l.eq.locp) then
           do n2 = 1, n_expansion(l)
           do n1 = n2,n_expansion(l)
              vnlnrm(n1,n2,l) = 0.0d0
              vnlnrm(n2,n1,l) = 0.0d0
           end do
           end do
        else
           do n2 = 1, n_expansion(l)
              do i=1,nrho
                 f(i)=vp(i,l)*wp(i,indx(n2,l))**2
              end do
              a=simp(nrho,f,drho)
              vnlnrm(n2,n2,l) = a
              do n1 = n2+1,n_expansion(l)
                 do i=1,nrho
                    f(i)=vp(i,l)*wp(i,indx(n1,l))*wp(i,indx(n2,l))
                 end do
                 a=simp(nrho,f,drho)
                 vnlnrm(n1,n2,l) = a
                 vnlnrm(n2,n1,l) = a
              end do
           end do
           if (n_expansion(l).eq.1) then
              vnlnrm(1,1,l) = 1/a
           else if (n_expansion(l).eq.2) then
              d = vnlnrm(1,1,l)*vnlnrm(2,2,l)
     >          - vnlnrm(1,2,l)*vnlnrm(2,1,l)
              q = vnlnrm(1,1,l)
              vnlnrm(1,1,l) = vnlnrm(2,2,l)/d
              vnlnrm(2,2,l) = q/d
              vnlnrm(1,2,l) = -vnlnrm(1,2,l)/d
              vnlnrm(2,1,l) = -vnlnrm(2,1,l)/d
           else
              call nwpw_matrix_invert(n_expansion(l),vnlnrm(1,1,l),nmax)
           end if
        end if
      end do

*     ************* compute ray fourier transforms *********************
      call integrate_kbpp_e_band_ray(version,
     >                            nrho,drho,lmax,locp,nmax,
     >                            n_extra,n_expansion,zv,
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
            if (l.ne.locp) then
              do n=1,n_expansion(l)
                call kbpp_band_filter_ray(nray,G_ray,wcut,
     >                                    vnl_ray(1,indx(n,l),1))
              end do
            end if
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
         if (l.ne.locp) then
            do n=1,n_expansion(l)
               call nwpw_spline(G_ray,vnl_ray(1,indx(n,l),1),nray,
     >                             0.0d0,0.0d0,
     >                             vnl_ray(1,indx(n,l),2),tmp_ray)
            end do
         end if
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
      lcount = nprj+1
      GO TO (950,940,930,920), lmax+1

        !::::::  f-wave  :::::::
  920   CONTINUE
        if (locp.ne.3) then
          do n=1,n_expansion(3)
             lcount = lcount-1
             n_prj(lcount) = n
             l_prj(lcount) = 3
             m_prj(lcount) = -3

             lcount = lcount-1
             n_prj(lcount) = n
             l_prj(lcount) = 3
             m_prj(lcount) = -2

             lcount = lcount-1
             n_prj(lcount) = n
             l_prj(lcount) = 3
             m_prj(lcount) = -1

             lcount = lcount-1
             n_prj(lcount) = n
             l_prj(lcount) = 3
             m_prj(lcount) = 0

             lcount = lcount-1
             n_prj(lcount) = n
             l_prj(lcount) = 3
             m_prj(lcount) = 1

             lcount = lcount-1
             n_prj(lcount) = n
             l_prj(lcount) = 3
             m_prj(lcount) = 2

             lcount = lcount-1
             n_prj(lcount) = n
             l_prj(lcount) = 3
             m_prj(lcount) = 3
          end do
        end if

        !::::  d-wave  ::::
  930   CONTINUE
        if (locp.ne.2) then
          do n=1,n_expansion(2)
             lcount = lcount-1
             n_prj(lcount) = n
             l_prj(lcount) = 2
             m_prj(lcount) = -2

             lcount = lcount-1
             n_prj(lcount) = n
             l_prj(lcount) = 2
             m_prj(lcount) = -1

             lcount = lcount-1
             n_prj(lcount) = n
             l_prj(lcount) = 2
             m_prj(lcount) = 0

             lcount = lcount-1
             n_prj(lcount) = n
             l_prj(lcount) = 2
             m_prj(lcount) = 1

             lcount = lcount-1
             n_prj(lcount) = n
             l_prj(lcount) = 2
             m_prj(lcount) = 2
          end do
        end if

        !::::  p-wave  ::::
  940   CONTINUE
        if (locp.ne.1) then
          do n=1,n_expansion(1)
             lcount = lcount-1
             n_prj(lcount) = n
             l_prj(lcount) = 1
             m_prj(lcount) = -1

             lcount = lcount-1
             n_prj(lcount) = n
             l_prj(lcount) = 1
             m_prj(lcount) = 0

             lcount = lcount-1
             n_prj(lcount) = n
             l_prj(lcount) = 1
             m_prj(lcount) = 1
          end do
        end if

        !::::  s-wave  ::::
  950   CONTINUE
        if (locp.ne.0) then
          do n=1,n_expansion(0)
             lcount = lcount-1
             n_prj(lcount) = n
             l_prj(lcount) = 0
             m_prj(lcount) = 0
          end do
        end if

      ierr=0
      return
      end
