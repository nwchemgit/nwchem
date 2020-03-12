*
* $Id$
*

      subroutine integrate_kbppv3e_ray(version,rlocal,
     >                            nrho,drho,lmax,locp,nmax,
     >                            n_extra,n_expansion,zv,
     >                            vp,wp,rho,f,cs,sn,
     >                            nray,G_ray,vl_ray,vnl_ray,
     >                            semicore,rho_sc_r,rho_sc_k_ray,
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

      integer nray
      double precision G_ray(nray)
      double precision vl_ray(nray)
      double precision vnl_ray(nray,0:(lmax+n_extra))

      logical semicore
      double precision rho_sc_r(nrho,2)
      double precision rho_sc_k_ray(nray,2)
      integer ierr

      integer np,taskid,MASTER
      parameter (MASTER=0)

*     *** local variables ****
      logical fast_erf,small_cell
      integer task_count
      integer k1,i,l,n,nb
      double precision pi,twopi,forpi
      double precision p0,p1,p2,p3,p
      double precision a,q,d
      integer indx(5,0:3)

*     **** Error function parameters ****
      real*8 yerf,xerf
      real*8 c1,c2,c3,c4,c5,c6
      parameter (c1=0.07052307840d0,c2=0.04228201230d0)
      parameter (c3=0.00927052720d0)
      parameter (c4=0.00015201430d0,c5=0.00027656720d0)
      parameter (c6=0.00004306380d0)


*     **** external functions ****
      logical          control_fast_erf,control_psp_semicore_small
      external         control_fast_erf,control_psp_semicore_small
      double precision dsum,simp,util_erf
      external         dsum,simp,util_erf


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

      fast_erf = control_fast_erf()
      small_cell = control_psp_semicore_small()

      pi=4.0d0*datan(1.0d0)
      twopi=2.0d0*pi
      forpi=4.0d0*pi

      if ((nrho/2)*2.EQ.nrho) then
        ierr=2
        return
      end if

      P0=DSQRT(FORPI)
      P1=DSQRT(3.0d0*FORPI)
      P2=DSQRT(15.0d0*FORPI)
      P3=DSQRT(105.0d0*FORPI)

*======================  Fourier transformation  ======================
      call dcopy(nray,0.0d0,0,vl_ray,1)
      call dcopy((lmax+1+n_extra)*nray,0.0d0,0,vnl_ray,1)
      call dcopy(2*nray,0.0d0,0,rho_sc_k_ray,1)
      task_count = -1
      DO 700 k1=2,nray
        task_count = task_count + 1
        if (mod(task_count,np).ne.taskid) go to 700

        Q=G_ray(k1)

        DO I=1,NRHO
          CS(I)=DCOS(Q*RHO(I))
          SN(I)=DSIN(Q*RHO(I))
        END DO

        GO TO (500,400,300,200), lmax+1


*::::::::::::::::::::::::::::::  f-wave  ::::::::::::::::::::::::::::::
  200   CONTINUE
        if (locp.ne.3) then
           do n=1,n_expansion(3)
              F(1)=0.0d0
              do I=2,NRHO
                A=SN(I)/(Q*RHO(I))
                A=15.0d0*(A-CS(I))/(Q*RHO(I))**2 - 6*A + CS(I)
                F(I)=A*WP(I,indx(n,3))*VP(I,3)
              end do
              D=P3*SIMP(NRHO,F,DRHO)/Q
              vnl_ray(k1,indx(n,3))=D
           end do
        end if
*::::::::::::::::::::::::::::::  d-wave  ::::::::::::::::::::::::::::::
  300   CONTINUE
        if (locp.ne.2) then
          do n=1,n_expansion(2)
             F(1)=0.0d0
             DO I=2,NRHO
               A=3.0d0*(SN(I)/(Q*RHO(I))-CS(I))/(Q*RHO(I))-SN(I)
               F(I)=A*WP(I,indx(n,2))*VP(I,2)
             END DO
             D=P2*SIMP(NRHO,F,DRHO)/Q
             vnl_ray(k1,indx(n,2))=D
          end do
        end if
*::::::::::::::::::::::::::::::  p-wave  ::::::::::::::::::::::::::::::
  400   CONTINUE
        if (locp.ne.1) then
           do n=1,n_expansion(1)
              F(1)=0.0d0
              DO I=2,NRHO
                F(I)=(SN(I)/(Q*RHO(I))-CS(I))*WP(I,indx(n,1))*VP(I,1)
              END DO
              P=P1*SIMP(NRHO,F,DRHO)/Q
              vnl_ray(k1,indx(n,1))=P
           end do
        end if
*::::::::::::::::::::::::::::::  s-wave  :::::::::::::::::::::::::::::::
  500   CONTINUE
        if (locp.ne.0) then
          do n=1,n_expansion(0)
             DO I=1,NRHO
               F(I)=SN(I)*WP(I,indx(n,0))*VP(I,0)
             END DO
             vnl_ray(k1,indx(n,0))=P0*SIMP(NRHO,F,DRHO)/Q
          end do
        end if

*::::::::::::::::::::::::::::::  local  :::::::::::::::::::::::::::::::
  600   CONTINUE


        if (version.eq.3) then
        DO  I=1,NRHO
          F(I)=RHO(I)*VP(I,locp)*SN(I)
        END DO
        vl_ray(k1)=SIMP(NRHO,F,DRHO)*FORPI/Q-ZV*FORPI*CS(NRHO)/(Q*Q)
        end if
 
        if (version.eq.4) then
        if (fast_erf) then
           do I=1,NRHO
             xerf=RHO(I)/rlocal
             yerf = (1.0d0
     >            + xerf*(c1 + xerf*(c2
     >            + xerf*(c3 + xerf*(c4
     >            + xerf*(c5 + xerf*c6))))))**4
             yerf = (1.0d0 - 1.0d0/yerf**4)
             F(I)=(RHO(I)*VP(I,locp)+ZV*yerf)*SN(I)
           end do
        else
           do I=1,NRHO
             xerf=RHO(I)/rlocal
             yerf = util_erf(xerf)
             F(I)=(RHO(I)*VP(I,locp)+ZV*yerf)*SN(I)
           end do
        end if
        vl_ray(k1)=SIMP(NRHO,F,DRHO)*FORPI/Q
        end if


*::::::::::::::::::::: semicore density :::::::::::::::::::::::::::::::
        if (semicore) then
           if (small_cell) then
              do i=1,nrho
                 f(i) = rho(i)*rho_sc_r(i,1)*sn(i)
              end do
           else
              do i=1,nrho
                 f(i) = rho(i)*dsqrt(rho_sc_r(i,1))*sn(i)
              end do
           end if
           rho_sc_k_ray(k1,1) = SIMP(nrho,f,drho)*forpi/Q

           do i=1,nrho
             f(i)=(sn(i)/(Q*rho(i))-cs(i))*rho_sc_r(i,2)*rho(i)
           end do
           P = SIMP(nrho,f,drho)*forpi/Q
           rho_sc_k_ray(k1,2)=P
        end if
    
  700 CONTINUE
      call Parallel_Vector_SumAll(2*nray,rho_sc_k_ray)
      call Parallel_Vector_SumAll(nray,vl_ray)
      call Parallel_Vector_Sumall((lmax+1+n_extra)*nray,vnl_ray)

*:::::::::::::::::::::::::::::::  G=0  ::::::::::::::::::::::::::::::::      
      if (version.eq.3) then
      DO I=1,NRHO
        F(I)=VP(I,locp)*RHO(I)**2
      END DO
      vl_ray(1)=FORPI*SIMP(NRHO,F,DRHO)+TWOPI*ZV*RHO(NRHO)**2
      end if

      if (version.eq.4) then
      if (fast_erf) then
         do I=1,NRHO
           xerf=RHO(I)/rlocal
           yerf = (1.0d0
     >          + xerf*(c1 + xerf*(c2
     >          + xerf*(c3 + xerf*(c4
     >          + xerf*(c5 + xerf*c6))))))**4
           yerf = (1.0d0 - 1.0d0/yerf**4)
           F(I)=(VP(I,locp)*RHO(I)+ZV*yerf)*RHO(I)
         end do
      else
         do I=1,NRHO
           xerf=RHO(I)/rlocal
           yerf = util_erf(xerf)
           F(I)=(VP(I,locp)*RHO(I)+ZV*yerf)*RHO(I)
         end do
      end if
      vl_ray(1)=FORPI*SIMP(NRHO,F,DRHO)
      end if

*     **** semicore density ****
      if (semicore) then
         if (small_cell) then
            do i=1,nrho
               f(i) = rho_sc_r(i,1)*rho(i)**2
            end do
         else
            do i=1,nrho
               f(i) = dsqrt(rho_sc_r(i,1))*rho(i)**2
            end do
         end if
         rho_sc_k_ray(1,1) = forpi*SIMP(nrho,f,drho)
         rho_sc_k_ray(1,2) = 0.0d0
      end if

      do l=0,lmax
         do n=1,n_expansion(l)
            vnl_ray(1,indx(n,l))=0.0d0
         end do
      end do
*     *** only j0 is non-zero at zero ****
      if (locp.ne.0) then
         do n=1,n_expansion(0)
            DO  I=1,NRHO
               F(I)=RHO(I)*WP(I,indx(n,0))*VP(I,0)
            END DO
            vnl_ray(1,indx(n,0))=P0*SIMP(NRHO,F,DRHO)
         end do
      end if

      IERR=0
      RETURN
      END



