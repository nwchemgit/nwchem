*
* $Id$
*
*     *****************************************************
*     *                                                   *
*     *            integrate_pawppv1                      *
*     *                                                   *
*     *****************************************************
      subroutine integrate_pawppv1(version,rlocal,
     >                            nrho,log_amesh,nbasis,nmax,lmax,icut,
     >                            zv,zc,sigma,
     >                            n_ps,l_ps,
     >                            v_ps,prj_ps,
     >                            phi_ae,dphi_ae,phi_ps,dphi_ps,
     >                            core_ae,core_ps,
     >                            rho,f1,f2,f3,f4,cs,sn,
     >                            nfft1,nfft2,nfft3,nproj,
     >                            G,vl,vlpaw,vnl,
     >                            n_prj,l_prj,m_prj,b_prj,
     >                            Gijl,
     >                            hartree_matrix,
     >                            comp_charge_matrix,comp_pot_matrix,
     >                            nray,G_ray,vl_ray,vlpaw_ray,vnl_ray,
     >                            tmp_ray,
     >                            filter,filter_scale,
     >                            ierr)
      implicit none
      integer          version
      double precision rlocal
      integer          nrho
      double precision log_amesh
      integer          nbasis,nmax,lmax,icut
      double precision zv,zc,sigma
      integer          n_ps(*),l_ps(*)
      double precision v_ps(nrho)
      double precision prj_ps(nrho,nbasis)

      double precision phi_ae(nrho,nbasis)
      double precision dphi_ae(nrho,nbasis)
      double precision phi_ps(nrho,nbasis)
      double precision dphi_ps(nrho,nbasis)

      double precision core_ae(nrho)
      double precision core_ps(nrho)

      double precision rho(nrho)
      double precision f1(nrho),f2(nrho),f3(nrho),f4(nrho)
      double precision cs(nrho)
      double precision sn(nrho)

      integer nfft1,nfft2,nfft3,nproj
      double precision G(nfft1/2+1,nfft2,nfft3,3)
      double precision vl(nfft1/2+1,nfft2,nfft3)
      double precision vlpaw(nfft1/2+1,nfft2,nfft3)
      double precision vnl(nfft1/2+1,nfft2,nfft3,nproj)
      integer          n_prj(nproj),l_prj(nproj),m_prj(nproj)
      integer          b_prj(nproj)
      double precision Gijl(nmax,nmax,0:lmax,5)
      real*8 hartree_matrix(nbasis,nbasis,nbasis,nbasis,0:(2*lmax))
      real*8 comp_charge_matrix(nbasis,nbasis,0:(2*lmax))
      real*8 comp_pot_matrix(nbasis,nbasis,0:(2*lmax))

      integer nray
      double precision G_ray(nray)
      double precision vl_ray(nray,2)
      double precision vlpaw_ray(nray,2)
      double precision vnl_ray(nray,nbasis,2)
      double precision tmp_ray(nray)
      logical filter
      double precision filter_scale

      integer ierr

#include "bafdecls.fh"
#include "errquit.fh"

      integer np,taskid,MASTER
      parameter (MASTER=0)

*     *** local variables ****
      integer lcount,task_count,nfft3d
      integer k1,k2,k3,i,j,i2,j2,k,l,nx,la,na,lb,nb,power_f
      integer sum_l,dif_l,sum_l2,dif_l2
      double precision pi,twopi,forpi
      double precision p0,p1,p2,p3,p
      double precision gx,gy,gz,a,q,d
      double precision ecut,wcut,dG,yp1
      double precision xerf,yerf
      double precision vgl,vgld(0:20)

*     **** external functions ****
      real*8   control_ecut,control_wcut,util_erf
      real*8   nwpw_splint,log_integrate_def,log_integrate_def0
      real*8   log_multipole_energy,nwpw_double_factorial
      external control_ecut,control_wcut,util_erf
      external nwpw_splint,log_integrate_def,log_integrate_def0
      external log_multipole_energy,nwpw_double_factorial

*      ****** DEBUG *****
      !call dcopy(nrho,0.0d0,0,v_ps,1)
      !rlocal = 2.0d0
       
*      ****** DEBUG *****

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

*     ****************************************************
*     **** comp_charge_matrix(nbasis,nbasis,0:2*lmax) ****
*     ****************************************************
      do i=1,nbasis
         do j=1,i
            do k=1,icut
               f1(k)=phi_ae(k,i)*phi_ae(k,j)-phi_ps(k,i)*phi_ps(k,j)
            end do
            do l=0,2*lmax
               d = log_integrate_def(2*l+2,f1,l,rho,log_amesh,icut)
               comp_charge_matrix(i,j,l) = d
               comp_charge_matrix(j,i,l) = d
            end do
         end do
      end do

*     *************************************************
*     **** comp_pot_matrix(nbasis,nbasis,0:2*lmax) ****
*     *************************************************
      do l=0,2*lmax
         k1 = l+2
         k2 = 2*l+2
         
         call nwpw_compcharge_gen_rgaussian(l,sigma,nrho,rho,f1)
         do k=1,icut
            f1(k) = f1(k)*rho(k)*rho(k)
         end do

         do i=1,nbasis
            do j=1,i
               do k=1,icut
                  f2(k) = phi_ps(k,i)*phi_ps(k,j)
               end do
               d = log_multipole_energy(l,icut,rho,
     >                                  k1,f1,
     >                                  k2,f2,log_amesh)
               comp_pot_matrix(i,j,l) = d
               comp_pot_matrix(j,i,l) = d
            end do
         end do
      end do

*     *********************************************************************
*     **** hartree_matrix(nbasis,nbasis,nbasis,nbasis,0:2*lmax)        ****
*     **** Note - This is the effective hartree matrix which includes  ****
*     **** comp_charge_matrix and comp_pot_matrix terms in it.         ****
*     *********************************************************************

      l = nbasis*nbasis*nbasis*nbasis*(2*lmax+1)
      call dcopy(l,0.0d0,0,hartree_matrix,1)
      do i=1,nbasis
         do j=1,i
            sum_l = l_ps(i)+l_ps(j)
            dif_l = abs(l_ps(i)-l_ps(j))
            k1 = sum_l + 2

            do k=1,icut
               f1(k) = phi_ae(k,i)*phi_ae(k,j)
               f3(k) = phi_ps(k,i)*phi_ps(k,j)
            end do

            do i2=1,nbasis
               do j2=1,i2
                  sum_l2 = l_ps(i2)+l_ps(j2)
                  dif_l2 = abs(l_ps(i2)-l_ps(j2))
                  k2 = sum_l2 + 2

                  do k=1,icut
                     f2(k) = phi_ae(k,i2)*phi_ae(k,j2)
                     f4(k) = phi_ps(k,i2)*phi_ps(k,j2)
                  end do

                  do l=0,2*lmax
                     d = dble(2*l+1)*nwpw_double_factorial(2*l+1)
     >                          *(sigma**(2*l+1))
                     vgl = 4.0d0*dsqrt(twopi)/d

                     if ((l.le.sum_l) .and.(l.ge.dif_l)) then
                     if ((l.le.sum_l2).and.(l.ge.dif_l2)) then
                        d = log_multipole_energy(l,icut,rho,
     >                                   k1,f1,
     >                                   k2,f2,log_amesh)
     >                    - log_multipole_energy(l,icut,rho,
     >                                   k1,f3,
     >                                   k2,f4,log_amesh)

                        hartree_matrix(i,j,i2,j2,l) = d
     >                         - 2.0d0*   comp_pot_matrix(i,j,l)
     >                                *comp_charge_matrix(i2,j2,l)
     >                         - vgl  *comp_charge_matrix(i,j,l)
     >                                *comp_charge_matrix(i2,j2,l)

                        hartree_matrix(j,i,i2,j2,l) = d
     >                         - 2.0d0*   comp_pot_matrix(j,i,l)
     >                                *comp_charge_matrix(i2,j2,l)
     >                         - vgl  *comp_charge_matrix(j,i,l)
     >                                *comp_charge_matrix(i2,j2,l)

                        hartree_matrix(i,j,j2,i2,l) = d
     >                         - 2.0d0*   comp_pot_matrix(i,j,l)
     >                                *comp_charge_matrix(j2,i2,l)
     >                         - vgl  *comp_charge_matrix(i,j,l)
     >                                *comp_charge_matrix(j2,i2,l)

                        hartree_matrix(j,i,j2,i2,l) = d
     >                         - 2.0d0*   comp_pot_matrix(j,i,l)
     >                                *comp_charge_matrix(j2,i2,l)
     >                         - vgl  *comp_charge_matrix(j,i,l)
     >                                *comp_charge_matrix(j2,i2,l)
                     end if
                     end if
                  end do
               end do
            end do
         end do
      end do



*     ******************************************************************************
*     ***********  1-electron psp operators - Normalization constants  *************
*     ******************************************************************************
      call dcopy(5*nmax*nmax*(lmax+1),0.0d0,0,Gijl,1)

*     *********************
*     **** 2 - overlap  ***
*     *********************
      do i=1,nbasis
         la = l_ps(i)
         na = n_ps(i) - la
         power_f = 2*la + 2
         do k=1,icut
           f1(k) = phi_ae(k,i)*phi_ae(k,i)
     >           - phi_ps(k,i)*phi_ps(k,i)
         end do
         d = log_integrate_def(power_f,f1,0,rho,log_amesh,icut)
         Gijl(na,na,la,2) = Gijl(na,na,la,2) + d

         do j=1,i-1
            lb = l_ps(j)
            nb = n_ps(j) - lb
            if (la.eq.lb) then
              do k=1,icut
                f1(k) = phi_ae(k,i)*phi_ae(k,j)
     >                - phi_ps(k,i)*phi_ps(k,j)
              end do
              d = log_integrate_def(power_f,f1,0,rho,log_amesh,icut)
              Gijl(na,nb,la,2) = Gijl(na,nb,la,2) + d
              Gijl(nb,na,la,2) = Gijl(nb,na,la,2) + d
            end if
         end do
      end do


*     *********************
*     **** 3 - kinetic  ***
*     *********************
      do i=1,nbasis
         la = l_ps(i)
         na = n_ps(i) - la
         power_f = 2*la
         do k=1,icut
           f1(k) = 0.5d0*(dphi_ae(k,i)*dphi_ae(k,i)
     >                   -dphi_ps(k,i)*dphi_ps(k,i))
     >           + 0.5d0*la*(la+1)*(phi_ae(k,i)*phi_ae(k,i) 
     >                             -phi_ps(k,i)*phi_ps(k,i))/rho(k)**2
         end do
         d = log_integrate_def(power_f,f1,0,rho,log_amesh,icut)
         Gijl(na,na,la,3) = Gijl(na,na,la,3) + d

         do j=1,i-1
            lb = l_ps(j)
            nb = n_ps(j) - lb
            if (la.eq.lb) then
              do k=1,icut
                f1(k) = 0.5d0*(dphi_ae(k,i)*dphi_ae(k,j)
     >                        -dphi_ps(k,i)*dphi_ps(k,j))
     >            + 0.5d0*la*(la+1)*(phi_ae(k,i)*phi_ae(k,j) -
     >                               phi_ps(k,i)*phi_ps(k,j))/rho(k)**2
              end do
              d = log_integrate_def(power_f,f1,0,rho,log_amesh,icut)
              Gijl(na,nb,la,3) = Gijl(na,nb,la,3) + d
              Gijl(nb,na,la,3) = Gijl(nb,na,la,3) + d
            end if
         end do
      end do

*     *********************************************************
*     **** 4 - Vlocal =  nloc*Vloc - ncmp*Vloc + nv*Vsrloc ****
*     ****    where Vloc = v_ps - zv*erf(r/rlocal)/r       ****
*     ****    and   Vsrloc = -zv/r - Vloc                  ****
*     *********************************************************

*     **** W2^(ion-electron) + W3^I(ion-electron)  -- add pseudo  - 4 ***
                
      q = zv/(forpi)          
      d = 2.0d0*zv/(dsqrt(twopi)*sigma)
      do i=1,nbasis                  
         la = l_ps(i)
         na = n_ps(i) - la
         Gijl(na,na,la,4) = Gijl(na,na,la,4) 
     >                    + q*comp_pot_matrix(i,i,0)
     >                    + d*Gijl(na,na,la,2)

         do j=1,i-1
            lb = l_ps(j)
            nb = n_ps(j) - lb
            if (la.eq.lb) then
              Gijl(na,nb,la,4) = Gijl(na,nb,la,4)
     >                         + q*comp_pot_matrix(i,j,0)
     >                         + d*Gijl(na,nb,la,2)
              Gijl(nb,na,la,4) = Gijl(nb,na,la,4)
     >                         + q*comp_pot_matrix(j,i,0)
     >                         + d*Gijl(nb,na,la,2)
            end if
         end do
      end do


*     **** U3^I = -tilde(n)*Vloc -- add pseudo  - 4 - PAW local psp ***
      do i=1,nbasis
         la = l_ps(i)
         na = n_ps(i) - la
         power_f = 2*la + 2
         do k=1,icut
           f1(k)= phi_ps(k,i)*phi_ps(k,i)*(-v_ps(k))
         end do
         d = log_integrate_def(power_f,f1,0,rho,log_amesh,icut)
         Gijl(na,na,la,4) = Gijl(na,na,la,4) + d

         do j=1,i-1
            lb = l_ps(j)
            nb = n_ps(j) - lb
            if (la.eq.lb) then

              do k=1,icut
                !xerf = rho(k)/rlocal
                !yerf = util_erf(xerf)
                !f1(k)= phi_ps(k,i)*phi_ps(k,j)*(-rho(k)*v_ps(k)+zv*yerf)
                f1(k)= phi_ps(k,i)*phi_ps(k,j)*(-v_ps(k))
              end do
              d = log_integrate_def(power_f,f1,0,rho,log_amesh,icut)
              Gijl(na,nb,la,4) = Gijl(na,nb,la,4) + d
              Gijl(nb,na,la,4) = Gijl(nb,na,la,4) + d

            end if
         end do
      end do

c*     **** U2^I  -- add pseudo  - 4 ***
c      call nwpw_compcharge_gen_rgaussian(0,sigma,nrho,rho,f2)
c      do k=1,icut
c         xerf = rho(k)/rlocal
c         yerf = util_erf(xerf)
c         f1(k)= f2(k)*yerf*zv
c      end do
c      d = log_integrate_def(1,f1,1,rho,log_amesh,icut)
c      a = log_integrate_def(1,f1,1,rho,log_amesh,icut-10)
c      write(*,*) "pseudo 4 d=",d,d-a,icut
c
c      do i=1,nbasis
c         la = l_ps(i)
c         na = n_ps(i) - la
c         !Gijl(na,na,la,4) = Gijl(na,na,la,4) + d*Gijl(na,na,la,2)
c
c         do j=1,i-1
c            lb = l_ps(j)
c            nb = n_ps(j) - lb
c            if (la.eq.lb) then
c              !Gijl(na,nb,la,4) = Gijl(na,nb,la,4) + d*Gijl(na,nb,la,2)
c              !Gijl(nb,na,la,4) = Gijl(nb,na,la,4) + d*Gijl(nb,na,la,2)
c            end if
c         end do
c      end do


*     **** add U1^I = -nv*(Zv/r) to pseudo - 4 ***
      do i=1,nbasis
         la = l_ps(i)
         na = n_ps(i) - la
         power_f = 2*la + 2
         do k=1,icut
c            xerf = rho(k)/sigma
c            yerf = util_erf(xerf)
c            f1(k) = (phi_ae(k,i)*phi_ae(k,i))*(-zv + zv*yerf)
            f1(k) = (phi_ae(k,i)*phi_ae(k,i))*(-zv)
         end do
         d = log_integrate_def0(power_f,f1,rho,log_amesh,icut)
         a = log_integrate_def0(power_f,f1,rho,log_amesh,icut-10)
         Gijl(na,na,la,4) = Gijl(na,na,la,4) + d

         do j=1,i-1
            lb = l_ps(j)
            nb = n_ps(j) - lb
            if (la.eq.lb) then
               do k=1,icut
c                  xerf = rho(k)/sigma
c                  yerf = util_erf(xerf)
c                  f1(k) = (phi_ae(k,i)*phi_ae(k,j))*(-zv + zv*yerf)
                  f1(k) = (phi_ae(k,i)*phi_ae(k,j))*(-zv)
               end do
               d = log_integrate_def0(power_f,f1,rho,log_amesh,icut)
               Gijl(na,nb,la,4) = Gijl(na,nb,la,4) + d
               Gijl(nb,na,la,4) = Gijl(nb,na,la,4) + d
            end if
         end do
      end do

*     **************************************
*     **** valence core matrix elements ****
*     **************************************
*     **** vcore - 5 ***
      if (Zc.gt.0.0d0) then
         do k=1,icut
            f3(k) = core_ae(k)*rho(k)**2
c            f4(k) = core_ps(k)*rho(k)**2   !*** core_ps is zero ****
         end do
         do i=1,nbasis
            la = l_ps(i)
            na = n_ps(i) - la
            power_f = 2*la + 2
            do k=1,icut
               f1(k) = phi_ae(k,i)*phi_ae(k,i)
c               f2(k) = phi_ps(k,i)*phi_ps(k,i)
            end do
            d = log_multipole_energy(0,icut,rho,power_f,
     >                              f1,2,f3,log_amesh)
c     >        - log_multipole_energy(0,icut,rho,power_f,
c     >                              f2,2,f4,log_amesh)
            Gijl(na,na,la,5) = Gijl(na,na,la,5) + d

            do j=1,i-1
               lb = l_ps(j)
               nb = n_ps(j) - lb
               if (la.eq.lb) then
                  do k=1,icut
                     f1(k) = phi_ae(k,i)*phi_ae(k,j)
c                     f2(k) = phi_ps(k,i)*phi_ps(k,j)
                  end do
                  d = log_multipole_energy(0,icut,rho,power_f,
     >                                    f1,2,f3,log_amesh)
c     >              - log_multipole_energy(0,icut,rho,power_f,
c     >                                    f2,2,f4,log_amesh)
                  Gijl(na,nb,la,5) = Gijl(na,nb,la,5) + d
                  Gijl(nb,na,la,5) = Gijl(nb,na,la,5) + d
               end if
            end do
         end do

*        **** add -Vzc to vcore - 5 ***
         do i=1,nbasis
            la = l_ps(i)
            na = n_ps(i) - la
            power_f = 2*la + 2
            do k=1,icut
               f1(k) = phi_ae(k,i)*phi_ae(k,i)/rho(k)
            end do
            d = Zc*log_integrate_def(power_f,f1,0,rho,log_amesh,icut)
            Gijl(na,na,la,5) = Gijl(na,na,la,5) - d

            do j=1,i-1
               lb = l_ps(j)
               nb = n_ps(j) - lb
               if (la.eq.lb) then
                  do k=1,icut
                     f1(k) = phi_ae(k,i)*phi_ae(k,j)/rho(k)
                  end do
                  d = Zc*log_integrate_def(power_f,f1,0,
     >                                     rho,log_amesh,icut)
                  Gijl(na,nb,la,5) = Gijl(na,nb,la,5) - d
                  Gijl(nb,na,la,5) = Gijl(nb,na,la,5) - d
               end if
            end do
         end do
      end if


      l = nmax*nmax*(lmax+1)
      call dcopy(l,      Gijl(1,1,0,3),1,Gijl(1,1,0,1),1)  !** kinetic **
      call daxpy(l,1.0d0,Gijl(1,1,0,4),1,Gijl(1,1,0,1),1)  !** Vsrloc - short range pseudo  **
      call daxpy(l,1.0d0,Gijl(1,1,0,5),1,Gijl(1,1,0,1),1)  !** vcore - valence-core   **


*     ************* compute ray fourier transforms *********************
      call integrate_paw_ray(version,rlocal,sigma,
     >                       nrho,log_amesh,zv,
     >                       nbasis,n_ps,l_ps,
     >                       v_ps,prj_ps,
     >                       rho,f1,cs,sn,
     >                       nray,
     >                       G_ray,vl_ray,vlpaw_ray,vnl_ray,
     >                       ierr)

*     **** filter the rays ****
      if (filter) then
         ecut = control_ecut()*filter_scale
         wcut = control_wcut()*filter_scale
         call kbpp_filter_ray(nray,G_ray,ecut,vl_ray)
         call kbpp_filter_ray(nray,G_ray,ecut,vlpaw_ray)
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

      !**** five point formula ***
      yp1 = ( -50.0d0*vlpaw_ray(2,1)
     >       + 96.0d0*vlpaw_ray(3,1)
     >       - 72.0d0*vlpaw_ray(4,1)
     >       + 32.0d0*vlpaw_ray(5,1)
     >       -  6.0d0*vlpaw_ray(6,1))/(24.0d0*dG)
      call nwpw_spline(G_ray(2),vlpaw_ray(2,1),nray-1,yp1,0.0d0,
     >                          vlpaw_ray(2,2),tmp_ray)

      do i=1,nbasis
         call nwpw_spline(G_ray,vnl_ray(1,i,1),nray,0.0d0,0.0d0,
     >                          vnl_ray(1,i,2),tmp_ray)
      end do
*     ======================  Fourier transformation  ======================
      call dcopy(nfft3d,0.0d0,0,vl,1)
      call dcopy(nfft3d,0.0d0,0,vlpaw,1)
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
           D=nwpw_splint(G_ray,vnl_ray(1,i,1),vnl_ray(1,i,2),nray,nx,Q)
           lcount = lcount+1
           vnl(k1,k2,k3,lcount)=D*GY*(3.0d0*(1.0d0-GZ*GZ)-4.0d0*GY*GY)
     >                          /dsqrt(24.0d0)
           lcount = lcount+1
           vnl(k1,k2,k3,lcount)=D*GX*GY*GZ
           lcount = lcount+1
           vnl(k1,k2,k3,lcount)=D*GY*(5.0d0*GZ*GZ-1.0d0)
     >                          /dsqrt(40.0d0)
           lcount = lcount+1
           vnl(k1,k2,k3,lcount)=D*GZ*(5.0d0*GZ*GZ-3.0d0)
     >                          /dsqrt(60.0d0)
           lcount = lcount+1
           vnl(k1,k2,k3,lcount)=D*GX*(5.0d0*GZ*GZ-1.0d0)
     >                          /dsqrt(40.0d0)
           lcount = lcount+1
           vnl(k1,k2,k3,lcount)=D*GZ*(GX*GX - GY*GY)
     >                          /2.0d0
           lcount = lcount+1
           vnl(k1,k2,k3,lcount)=D*GX*(4.0d0*GX*GX-3.0d0*(1.0d0-GZ*GZ))
     >                          /dsqrt(24.0d0)

cccccc
c           lcount = lcount+1
c           vnl(k1,k2,k3,lcount)=D*GX*(4.0d0*GX*GX-3.0d0*(1.0d0-GZ*GZ))
c     >                          /dsqrt(24.0d0)
c           lcount = lcount+1
c           vnl(k1,k2,k3,lcount)=D*GY*(3.0d0*(1.0d0-GZ*GZ)-4.0d0*GY*GY)
c     >                          /dsqrt(24.0d0)
c           lcount = lcount+1
c           vnl(k1,k2,k3,lcount)=D*GZ*(GX*GX - GY*GY)
c     >                          /2.0d0
c           lcount = lcount+1
c           vnl(k1,k2,k3,lcount)=D*GX*GY*GZ
c           lcount = lcount+1
c           vnl(k1,k2,k3,lcount)=D*GX*(5.0d0*GZ*GZ-1.0d0)
c     >                          /dsqrt(40.0d0)
c           lcount = lcount+1
c           vnl(k1,k2,k3,lcount)=D*GY*(5.0d0*GZ*GZ-1.0d0)
c     >                          /dsqrt(40.0d0)
c           lcount = lcount+1
c           vnl(k1,k2,k3,lcount)=D*GZ*(5.0d0*GZ*GZ-3.0d0)
c     >                          /dsqrt(60.0d0)
          end if
*         ::::::::::::::::::::::::::::::  d-wave  ::::::::::::::::::::::::::::::
*         ::::::::::::::::::::::: -1 phase factor added ::::::::::::::::::::::::
          if (la.eq.2) then
            D=nwpw_splint(G_ray,vnl_ray(1,i,1),vnl_ray(1,i,2),nray,nx,Q)
            lcount = lcount+1
            vnl(k1,k2,k3,lcount)=-D*GX*GY
            lcount = lcount+1
            vnl(k1,k2,k3,lcount)=-D*GY*GZ
            lcount = lcount+1
            vnl(k1,k2,k3,lcount)=-D*(3.0d0*GZ*GZ-1.0d0)
     >                            /(2.0d0*dsqrt(3.0d0))
            lcount = lcount+1
            vnl(k1,k2,k3,lcount)=-D*GZ*GX
            lcount = lcount+1
            vnl(k1,k2,k3,lcount)=-D*(GX*GX-GY*GY)/(2.0d0)

c            lcount = lcount+1
c            vnl(k1,k2,k3,lcount)=D*(3.0d0*GZ*GZ-1.0d0)
c     >                            /(2.0d0*dsqrt(3.0d0))
c            lcount = lcount+1
c            vnl(k1,k2,k3,lcount)=D*GX*GY
c            lcount = lcount+1
c            vnl(k1,k2,k3,lcount)=D*GY*GZ
c            lcount = lcount+1
c            vnl(k1,k2,k3,lcount)=D*GZ*GX
c            lcount = lcount+1
c            vnl(k1,k2,k3,lcount)=D*(GX*GX-GY*GY)/(2.0d0)
          end if
*         ::::::::::::::::::::::::::::::  p-wave  ::::::::::::::::::::::::::::::
*         ::::::::::::::::::::::: -1 phase factor added ::::::::::::::::::::::::
          if (la.eq.1) then
           P=nwpw_splint(G_ray,vnl_ray(1,i,1),vnl_ray(1,i,2),nray,nx,Q)
           lcount = lcount+1
           vnl(k1,k2,k3,lcount)=-P*GY
           lcount = lcount+1
           vnl(k1,k2,k3,lcount)=-P*GZ
           lcount = lcount+1
           vnl(k1,k2,k3,lcount)=-P*GX
          end if
*         ::::::::::::::::::::::::::::::  s-wave  :::::::::::::::::::::::::::::::
          if (la.eq.0) then
            P=nwpw_splint(G_ray,vnl_ray(1,i,1),vnl_ray(1,i,2),nray,nx,Q)
            lcount = lcount+1
            vnl(k1,k2,k3,lcount)=P
          end if

        end do

*       ::::::::::::::::::::::::::::::  local  :::::::::::::::::::::::::::::::
  600   CONTINUE
        P = nwpw_splint(G_ray(2),vl_ray(2,1),vl_ray(2,2),nray-1,nx-1,Q)
        vl(k1,k2,k3)=P
        P = nwpw_splint(G_ray(2),vlpaw_ray(2,1),vlpaw_ray(2,2),
     >                  nray-1,nx-1,Q)
        vlpaw(k1,k2,k3)=P

       
  700 CONTINUE
      call Parallel_Vector_SumAll(nfft3d,vl)
      call Parallel_Vector_SumAll(nfft3d,vlpaw)
      call Parallel_Vector_Sumall(nproj*nfft3d,vnl)

*     :::::::::::::::::::::::::::::::  G=0  ::::::::::::::::::::::::::::::::      
*     **** local potential ****
      vl(1,1,1)    = vl_ray(1,1)
      vlpaw(1,1,1) = vlpaw_ray(1,1)

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
          b_prj(lcount) = i
    
          lcount = lcount+1
          n_prj(lcount) = na
          l_prj(lcount) = 3
          m_prj(lcount) = -2
          b_prj(lcount) = i
           
          lcount = lcount+1
          n_prj(lcount) = na
          l_prj(lcount) = 3
          m_prj(lcount) = -1
          b_prj(lcount) = i
     
          lcount = lcount+1
          n_prj(lcount) = na
          l_prj(lcount) = 3
          m_prj(lcount) = 0
          b_prj(lcount) = i
           
          lcount = lcount+1
          n_prj(lcount) = na
          l_prj(lcount) = 3
          m_prj(lcount) = 1
          b_prj(lcount) = i
           
          lcount = lcount+1
          n_prj(lcount) = na
          l_prj(lcount) = 3
          m_prj(lcount) = 2
          b_prj(lcount) = i
     
          lcount = lcount+1
          n_prj(lcount) = na
          l_prj(lcount) = 3
          m_prj(lcount) = 3
          b_prj(lcount) = i
        end if

        !::::  d-wave  ::::
        if (la.eq.2) then
          lcount = lcount+1
          n_prj(lcount) = na
          l_prj(lcount) = 2
          m_prj(lcount) = -2
          b_prj(lcount) = i

          lcount = lcount+1
          n_prj(lcount) = na
          l_prj(lcount) = 2
          m_prj(lcount) = -1
          b_prj(lcount) = i
          
          lcount = lcount+1
          n_prj(lcount) = na
          l_prj(lcount) = 2
          m_prj(lcount) = 0
          b_prj(lcount) = i
          
          lcount = lcount+1
          n_prj(lcount) = na
          l_prj(lcount) = 2
          m_prj(lcount) = 1
          b_prj(lcount) = i
          
          lcount = lcount+1
          n_prj(lcount) = na
          l_prj(lcount) = 2
          m_prj(lcount) = 2
          b_prj(lcount) = i
        end if

        !::::  p-wave  ::::
        if (la.eq.1) then
          lcount = lcount+1
          n_prj(lcount) = na
          l_prj(lcount) = 1
          m_prj(lcount) = -1
          b_prj(lcount) = i

          lcount = lcount+1
          n_prj(lcount) = na
          l_prj(lcount) = 1
          m_prj(lcount) = 0
          b_prj(lcount) = i

          lcount = lcount+1
          n_prj(lcount) = na
          l_prj(lcount) = 1
          m_prj(lcount) = 1
          b_prj(lcount) = i
        end if

        !::::  s-wave  ::::
        if (la.eq.0) then
          lcount = lcount+1
          n_prj(lcount) = na
          l_prj(lcount) = 0
          m_prj(lcount) = 0
          b_prj(lcount) = i
          vnl(1,1,1,lcount)=vnl_ray(1,i,1) !*** j0!=0 at G=0***
        end if
      end do




      IERR=0
      RETURN
      END



