*
* $Id$
*
      subroutine integrate_kbppv3e_new(version,rlocal,
     >                            nrho,drho,lmax,locp,nmax,
     >                            n_extra,n_expansion,zv,
     >                            vp,wp,rho,f,cs,sn,
     >                            nfft3d,nprj,
     >                            G,vl,vnl,
     >                            n_prj,l_prj,m_prj,b_prj,vnlnrm,
     >                            semicore,rho_sc_r,rho_sc_k,
     >                            nray,G_ray,vl_ray,vnl_ray,
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
      double precision vl(nfft3d)
      double precision vnl(nfft3d,nprj)
      integer          n_prj(nprj),l_prj(nprj),m_prj(nprj)
      integer          b_prj(nprj)
      double precision vnlnrm(nmax,nmax,0:lmax)

      logical semicore
      double precision rho_sc_r(nrho,2)
      double precision rho_sc_k(nfft3d,4)

      integer nray
      double precision G_ray(nray)
      double precision vl_ray(nray,2)
      double precision vnl_ray(nray,0:(lmax+n_extra),2)
      double precision rho_sc_k_ray(nray,2,2)
      double precision tmp_ray(nray)
      logical filter


      integer ierr

      integer np,taskid,MASTER
      integer np_i,np_j,taskid_i,taskid_j,countj
      parameter (MASTER=0)

*     *** local variables ****
      integer lcount
      integer k1,k2,k3,i,l,nx,pzero,zero,n,nb,n1,n2
      double precision pi,twopi,forpi
      double precision p0,p1,p2,p3,p
      double precision gx,gy,gz,a,q,d
      double precision ecut,wcut,dG,yp1
      integer indx(5,0:3)

*     **** external functions ****
      double precision dsum,simp,control_ecut,control_wcut
      double precision nwpw_splint
      external         dsum,simp,control_ecut,control_wcut
      external         nwpw_splint

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

      pi=4.0d0*datan(1.0d0)
      twopi=2.0d0*pi
      forpi=4.0d0*pi

      if (lmax.gt.3) then
        ierr=1
        return
      end if
      if ((nrho/2)*2.EQ.nrho) then
        ierr=2
        return
      end if

      P0=DSQRT(FORPI)
      P1=DSQRT(3.0d0*FORPI)
      P2=DSQRT(15.0d0*FORPI)
      P3=DSQRT(105.0d0*FORPI)

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
      call integrate_kbppv3e_ray(version,rlocal,
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
         call kbpp_filter_ray(nray,G_ray,ecut,vl_ray)
         do l=0,lmax
            if (l.ne.locp) then
              do n=1,n_expansion(l)
                 call kbpp_filter_ray(nray,G_ray,wcut,
     >                                vnl_ray(1,indx(n,l),1))
              end do
            end if
         end do
         if (semicore) then
           call kbpp_filter_ray(nray,G_ray,ecut,rho_sc_k_ray(1,1,1))
           call kbpp_filter_ray(nray,G_ray,ecut,rho_sc_k_ray(1,2,1))
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
               call nwpw_spline(G_ray,
     >                          vnl_ray(1,indx(n,l),1),
     >                          nray,0.0d0,0.0d0,
     >                          vnl_ray(1,indx(n,l),2),tmp_ray)
            end do
         end if
      end do
      if (semicore) then
         call nwpw_spline(G_ray,rho_sc_k_ray(1,1,1),nray,0.0d0,0.0d0,
     >                          rho_sc_k_ray(1,1,2),tmp_ray)
         call nwpw_spline(G_ray,rho_sc_k_ray(1,2,1),nray,0.0d0,0.0d0,
     >                          rho_sc_k_ray(1,2,2),tmp_ray)
      end if
*     ======================  Fourier transformation  ======================
      call dcopy(nfft3d,0.0d0,0,vl,1)
      call dcopy(nprj*nfft3d,0.0d0,0,vnl,1)
      call dcopy(4*nfft3d,0.0d0,0,rho_sc_k,1)

*     ***** find the G==0 point in the lattice *****
      call D3dB_ijktoindexp(1,1,1,1,zero,pzero)
      
      countj = -1
      DO 700 k1=1,nfft3d

        countj = mod(countj+1,np_j)
        if (countj.ne.taskid_j) go to 700
        if ((pzero.eq.taskid_i).and.(k1.eq.zero)) go to 700

        Q=DSQRT(G(k1,1)**2
     >         +G(k1,2)**2
     >         +G(k1,3)**2)
        nx = (Q/dG) + 1.0d0

        
        GX=G(k1,1)/Q
        GY=G(k1,2)/Q
        GZ=G(k1,3)/Q

        lcount = nprj+1
        GO TO (500,400,300,200), LMAX+1


*       ::::::::::::::::::::::::::::::  f-wave  ::::::::::::::::::::::::::::::
  200   CONTINUE
        if (locp.ne.3) then
          do n=1,n_expansion(3)
            D = nwpw_splint(G_ray,vnl_ray(1,indx(n,3),1),
     >                            vnl_ray(1,indx(n,3),2),nray,nx,Q)

             lcount = lcount-1
             vnl(k1,lcount)=D*GY*(3.0d0*(1.0d0-GZ*GZ)-4.0d0*GY*GY)
     >                            /dsqrt(24.0d0)
             lcount = lcount-1
             vnl(k1,lcount)=D*GX*GY*GZ
             lcount = lcount-1
             vnl(k1,lcount)=D*GY*(5.0d0*GZ*GZ-1.0d0)
     >                            /dsqrt(40.0d0)
             lcount = lcount-1
             vnl(k1,lcount)=D*GZ*(5.0d0*GZ*GZ-3.0d0)
     >                            /dsqrt(60.0d0)
             lcount = lcount-1
             vnl(k1,lcount)=D*GX*(5.0d0*GZ*GZ-1.0d0)
     >                            /dsqrt(40.0d0)
             lcount = lcount-1
             vnl(k1,lcount)=D*GZ*(GX*GX - GY*GY)
     >                            /2.0d0
             lcount = lcount-1
             vnl(k1,lcount)=D*GX*(4.0d0*GX*GX-3.0d0*(1.0d0-GZ*GZ))
     >                            /dsqrt(24.0d0)
          end do
        end if
*       ::::::::::::::::::::::::::::::  d-wave  ::::::::::::::::::::::::::::::
  300   CONTINUE
        if (locp.ne.2) then
          do n=1,n_expansion(2)
             D = nwpw_splint(G_ray,vnl_ray(1,indx(n,2),1),
     >                             vnl_ray(1,indx(n,2),2),nray,nx,Q)
             lcount = lcount-1
             vnl(k1,lcount)=D*GX*GY
             lcount = lcount-1
             vnl(k1,lcount)=D*GY*GZ
             lcount = lcount-1
             vnl(k1,lcount)=D*(3.0d0*GZ*GZ-1.0d0)
     >                             /(2.0d0*dsqrt(3.0d0))
             lcount = lcount-1
             vnl(k1,lcount)=D*GZ*GX
             lcount = lcount-1
             vnl(k1,lcount)=D*(GX*GX-GY*GY)/(2.0d0)
          end do
        end if
*       ::::::::::::::::::::::::::::::  p-wave  ::::::::::::::::::::::::::::::
  400   CONTINUE
        if (locp.ne.1) then
          do n=1,n_expansion(1)
             P = nwpw_splint(G_ray,vnl_ray(1,indx(n,1),1),
     >                             vnl_ray(1,indx(n,1),2),nray,nx,Q)
             lcount = lcount-1
             vnl(k1,lcount)=P*GY
             lcount = lcount-1
             vnl(k1,lcount)=P*GZ
             lcount = lcount-1
             vnl(k1,lcount)=P*GX
          end do
        end if
*       ::::::::::::::::::::::::::::::  s-wave  :::::::::::::::::::::::::::::::
  500   CONTINUE
        if (locp.ne.0) then
          do n=1,n_expansion(0)
             P = nwpw_splint(G_ray,vnl_ray(1,indx(n,0),1),
     >                             vnl_ray(1,indx(n,0),2),nray,nx,Q)
             lcount = lcount-1
             vnl(k1,lcount)=P
          end do
        end if
*       ::::::::::::::::::::::::::::::  local  :::::::::::::::::::::::::::::::
  600   CONTINUE
        P = nwpw_splint(G_ray(2),vl_ray(2,1),vl_ray(2,2),nray-1,nx-1,Q)
        vl(k1)=P
       
*       ::::::::::::::::::::: semicore density :::::::::::::::::::::::::::::::
        if (semicore) then
           P = nwpw_splint(G_ray,rho_sc_k_ray(1,1,1),
     >                           rho_sc_k_ray(1,1,2),nray,nx,Q)
           rho_sc_k(k1,1) = P

           P = nwpw_splint(G_ray,rho_sc_k_ray(1,2,1),
     >                           rho_sc_k_ray(1,2,2),nray,nx,Q)
           rho_sc_k(k1,2)=P*GX
           rho_sc_k(k1,3)=P*GY
           rho_sc_k(k1,4)=P*GZ

        end if
    
  700 CONTINUE
      call D1dB_Vector_SumAll(4*nfft3d,rho_sc_k)
      call D1dB_Vector_SumAll(nfft3d,vl)
      call D1dB_Vector_Sumall(nprj*nfft3d,vnl)



*     :::::::::::::::::::::::::::::::  G=0  ::::::::::::::::::::::::::::::::      
      if (pzero.eq.taskid_i) then

*        **** local potential ****
         vl(zero)=vl_ray(1,1)

*        **** semicore density ****
         if (semicore) then
            rho_sc_k(zero,1) = rho_sc_k_ray(1,1,1)
            rho_sc_k(zero,2) = 0.0d0
            rho_sc_k(zero,3) = 0.0d0
            rho_sc_k(zero,4) = 0.0d0
         end if

         do l=1,nprj
           vnl(zero,l)=0.0d0
         end do
*        *** only j0 is non-zero at zero ****
         if (locp.ne.0) then
            do n=1,n_expansion(0)
               vnl(zero,n_expansion(0)-n+1)=vnl_ray(1,indx(n,0),1)
            end do
         end if
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
             b_prj(lcount) = indx(n,3)+1

             lcount = lcount-1
             n_prj(lcount) = n
             l_prj(lcount) = 3
             m_prj(lcount) = -2
             b_prj(lcount) = indx(n,3)+1

             lcount = lcount-1
             n_prj(lcount) = n
             l_prj(lcount) = 3
             m_prj(lcount) = -1
             b_prj(lcount) = indx(n,3)+1

             lcount = lcount-1
             n_prj(lcount) = n
             l_prj(lcount) = 3
             m_prj(lcount) = 0
             b_prj(lcount) = indx(n,3)+1

             lcount = lcount-1
             n_prj(lcount) = n
             l_prj(lcount) = 3
             m_prj(lcount) = 1
             b_prj(lcount) = indx(n,3)+1

             lcount = lcount-1
             n_prj(lcount) = n
             l_prj(lcount) = 3
             m_prj(lcount) = 2
             b_prj(lcount) = indx(n,3)+1

             lcount = lcount-1
             n_prj(lcount) = n
             l_prj(lcount) = 3
             m_prj(lcount) = 3
             b_prj(lcount) = indx(n,3)+1
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
             b_prj(lcount) = indx(n,2)+1

             lcount = lcount-1
             n_prj(lcount) = n
             l_prj(lcount) = 2
             m_prj(lcount) = -1
             b_prj(lcount) = indx(n,2)+1

             lcount = lcount-1
             n_prj(lcount) = n
             l_prj(lcount) = 2
             m_prj(lcount) = 0
             b_prj(lcount) = indx(n,2)+1

             lcount = lcount-1
             n_prj(lcount) = n
             l_prj(lcount) = 2
             m_prj(lcount) = 1
             b_prj(lcount) = indx(n,2)+1

             lcount = lcount-1
             n_prj(lcount) = n
             l_prj(lcount) = 2
             m_prj(lcount) = 2
             b_prj(lcount) = indx(n,2)+1
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
             b_prj(lcount) = indx(n,1)+1

             lcount = lcount-1
             n_prj(lcount) = n
             l_prj(lcount) = 1
             m_prj(lcount) = 0
             b_prj(lcount) = indx(n,1)+1

             lcount = lcount-1
             n_prj(lcount) = n
             l_prj(lcount) = 1
             m_prj(lcount) = 1
             b_prj(lcount) = indx(n,1)+1
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
             b_prj(lcount) = indx(n,0)+1
          end do
        end if

      IERR=0
      RETURN
      END



