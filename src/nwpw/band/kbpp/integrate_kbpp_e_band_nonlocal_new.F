*
* $Id$
*
      subroutine integrate_kbpp_e_band_nonlocal_new(version,kvec,
     >                            nrho,drho,lmax,locp,nmax,
     >                            n_extra,n_expansion,zv,
     >                            vp,wp,rho,f,cs,sn,
     >                            nfft1,nfft2,nfft3,nprj,
     >                            G,vnl,
     >                            nray,G_ray,vnl_ray,
     >                            ierr)
      implicit none
      integer          version
      double precision kvec(3)
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
      double precision vnl(nfft1,nfft2,nfft3,nprj)

      integer nray
      double precision G_ray(nray)
      double precision vnl_ray(nray,0:(lmax+n_extra),2)

      integer ierr

*     *** local variables ****
      integer np,taskid,MASTER
      parameter (MASTER=0)

      integer lcount,task_count,nfft3d
      integer k1,k2,k3,i,l,nx,n,nb
      double precision p
      double precision gx,gy,gz,q,d,dG
      integer indx(5,0:3)

*     **** external functions ****
      double precision dsum,simp,nwpw_splint
      external         dsum,simp,nwpw_splint

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

      if (lmax.GT.3) then
        write(*,*)"non-local psp not generated: lmax exceeds 3"
        ierr=1
        return
      end if
      if ((nrho/2)*2.EQ.nrho) then
        write(*,*)"non-local psp not generated: nrho is not odd"
        ierr=2
        return
      end if

*======================  Fourier transformation  ======================
      dG = G_ray(3)-G_ray(2)
      call dcopy(nprj*nfft3d,0.0d0,0,vnl,1)
      task_count = -1
      DO 700 k3=1,nfft3
      DO 700 k2=1,nfft2
      DO 700 k1=1,nfft1
        task_count = task_count + 1
        if (mod(task_count,np).ne.taskid) go to 700
        gx=G(k1,k2,k3,1)+kvec(1)
        gy=G(k1,k2,k3,2)+kvec(2)
        gz=G(k1,k2,k3,3)+kvec(3)

        Q  =DSQRT(gx**2 + gy**2 + gz**2)
        nx = (Q/dG) + 1.0d0

        if (dabs(Q).gt.1.0d-9) then

           gx=gx/Q
           gy=gy/Q
           gz=gz/Q
           do i=1,nrho
             cs(i)=DCOS(Q*rho(i))
             sn(i)=DSIN(Q*rho(i))
           end do

           lcount = nprj+1
           GO TO (500,400,300,200), lmax+1

*::::::::::::::::::::::::::::::  f-wave  ::::::::::::::::::::::::::::::
  200      CONTINUE
           if (locp.ne.3) then
              do n=1,n_expansion(3)
                 D = nwpw_splint(G_ray,vnl_ray(1,indx(n,3),1),
     >                           vnl_ray(1,indx(n,3),2),nray,nx,Q)
                 lcount = lcount-1
                 vnl(k1,k2,k3,lcount)=D*GX*(4.0d0*GX*GX
     >                                     -3.0d0*(1.0d0-GZ*GZ))
     >                             /dsqrt(24.0d0)
                 lcount = lcount-1
                 vnl(k1,k2,k3,lcount)=D*GY*(3.0d0*(1.0d0-GZ*GZ)
     >                                     -4.0d0*GY*GY)
     >                             /dsqrt(24.0d0)
                 lcount = lcount-1
                 vnl(k1,k2,k3,lcount)=D*GZ*(GX*GX - GY*GY)
     >                             /2.0d0
                 lcount = lcount-1
                 vnl(k1,k2,k3,lcount)=D*GX*GY*GZ
                 lcount = lcount-1
                 vnl(k1,k2,k3,lcount)=D*GX*(5.0d0*GZ*GZ-1.0d0)
     >                             /dsqrt(40.0d0)
                 lcount = lcount-1
                 vnl(k1,k2,k3,lcount)=D*GY*(5.0d0*GZ*GZ-1.0d0)
     >                             /dsqrt(40.0d0)
                 lcount = lcount-1
                 vnl(k1,k2,k3,lcount)=D*GZ*(5.0d0*GZ*GZ-3.0d0)
     >                             /dsqrt(60.0d0)
              end do
           end if
*::::::::::::::::::::::::::::::  d-wave  ::::::::::::::::::::::::::::::
  300      CONTINUE
           if (locp.ne.2) then
              do n=1,n_expansion(2)
                 D = nwpw_splint(G_ray,vnl_ray(1,indx(n,2),1),
     >                           vnl_ray(1,indx(n,2),2),nray,nx,Q)
                 lcount = lcount-1
c                 vnl(k1,k2,k3,lcount)=D*(2.0d0*GZ*GZ-GX*GX-GY*GY)
c     >                          /(2.0d0*dsqrt(3.0d0))
                 vnl(k1,k2,k3,lcount)=D*(3.0d0*GZ*GZ-1.0d0)
     >                                 /(2.0d0*dsqrt(3.0d0))
                 lcount = lcount-1
                 vnl(k1,k2,k3,lcount)=D*GX*GY
                 lcount = lcount-1
                 vnl(k1,k2,k3,lcount)=D*GY*GZ
                 lcount = lcount-1
                 vnl(k1,k2,k3,lcount)=D*GZ*GX
                 lcount = lcount-1
                 vnl(k1,k2,k3,lcount)=D*(GX*GX-GY*GY)/(2.0d0)
              end do
           end if
*::::::::::::::::::::::::::::::  p-wave  ::::::::::::::::::::::::::::::
  400      CONTINUE
           if (locp.ne.1) then
             do n=1,n_expansion(1)
                P = nwpw_splint(G_ray,vnl_ray(1,indx(n,1),1),
     >                                vnl_ray(1,indx(n,1),2),nray,nx,Q)
                lcount = lcount-1
                vnl(k1,k2,k3,lcount)=P*GX
                lcount = lcount-1
                vnl(k1,k2,k3,lcount)=P*GY
                lcount = lcount-1
                vnl(k1,k2,k3,lcount)=P*GZ
             end do
           end if
*::::::::::::::::::::::::::::::  s-wave  :::::::::::::::::::::::::::::::
  500      CONTINUE
           if (locp.ne.0) then
             do n=1,n_expansion(0)
                lcount = lcount-1
                vnl(k1,k2,k3,lcount)= nwpw_splint(G_ray,
     >                                       vnl_ray(1,indx(n,0),1),
     >                                       vnl_ray(1,indx(n,0),2),
     >                                       nray,nx,Q)
             end do
           end if

  600      CONTINUE

*:::::::::::::::::::::::::::::::  G+k=0  ::::::::::::::::::::::::::::::::      
      else

         do l=1,nprj
           vnl(k1,k2,k3,l)=0.0d0
         end do
*        *** only j0 is non-zero at zero ****
         if (locp.ne.0) then
            do n=1,n_expansion(0)
               vnl(k1,k2,k3,n_expansion(0)-n+1)=vnl_ray(1,indx(n,0),1)
            end do
         end if

      end if

    
  700 CONTINUE

      call Parallel_Vector_SumAll(nprj*nfft3d,vnl)


      ierr=0
      return
      end



