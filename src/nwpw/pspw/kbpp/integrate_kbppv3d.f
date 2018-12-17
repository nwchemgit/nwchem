*
* $Id$
*

      subroutine integrate_kbppv3d(version,rlocal,
     >                            nrho,drho,lmax,locp,zv,
     >                            vp,wp,rho,f,cs,sn,
     >                            nfft3d,lmmax,
     >                            G,vl,vnl,
     >                            n_prj,l_prj,m_prj,b_prj,vnlnrm,
     >                            semicore,rho_sc_r,rho_sc_k,
     >                            ierr)
      implicit none
      integer          version
      double precision rlocal
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

      integer nfft3d,lmmax
      double precision G(nfft3d,3)
      double precision vl(nfft3d)
      double precision vnl(nfft3d,lmmax)
      integer          n_prj(lmmax),l_prj(lmmax),m_prj(lmmax)
      integer          b_prj(lmmax)
      double precision vnlnrm(0:lmax)

      logical semicore
      double precision rho_sc_r(nrho,2)
      double precision rho_sc_k(nfft3d,4)

      integer ierr

      integer np,taskid,MASTER
      integer np_i,np_j,taskid_i,taskid_j,countj
      parameter (MASTER=0)

*     *** local variables ****
      logical fast_erf,small_cell
      integer lcount
      integer k1,k2,k3,i,l,pzero,zero
      double precision pi,twopi,forpi
      double precision p0,p1,p2,p3,p
      double precision gx,gy,gz,a,q,d
      double precision tollz
      parameter(tollz=1d-16)

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

      call Parallel2d_np_i(np_i)
      call Parallel2d_np_j(np_j)
      call Parallel2d_taskid_i(taskid_i)
      call Parallel2d_taskid_j(taskid_j)

      fast_erf   = control_fast_erf()
      small_cell = control_psp_semicore_small()

      pi=4.0d0*datan(1.0d0)
      twopi=2.0d0*pi
      forpi=4.0d0*pi

      IF(LMMAX.GT.16) THEN
        IERR=1
        RETURN
      ENDIF
      IF((NRHO/2)*2.EQ.NRHO) THEN
        IERR=2
        RETURN
      ENDIF

      P0=DSQRT(FORPI)
      P1=DSQRT(3.0d0*FORPI)
      P2=DSQRT(15.0d0*FORPI)
      P3=DSQRT(105.0d0*FORPI)

*::::::::::::::::::  Define non-local pseudopotential  ::::::::::::::::
      do l=0,lmax
        if (l.ne.locp) then
          do I=1,nrho
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
          vnlnrm(l) = (1.0d0/a)
        else
          vnlnrm(l) = 0.0d0
        end if
      end do

*======================  Fourier transformation  ======================
      call dcopy(nfft3d,0.0d0,0,vl,1)
      call dcopy(lmmax*nfft3d,0.0d0,0,vnl,1)
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
        if (abs(q).lt.tollz) go to 700

        
        GX=G(k1,1)/Q
        GY=G(k1,2)/Q
        GZ=G(k1,3)/Q
        DO I=1,NRHO
          CS(I)=DCOS(Q*RHO(I))
          SN(I)=DSIN(Q*RHO(I))
        END DO

        lcount = lmmax+1
        GO TO (500,400,300,200), LMAX+1


*::::::::::::::::::::::::::::::  f-wave  ::::::::::::::::::::::::::::::
  200   CONTINUE
        if (locp.ne.3) then
           F(1)=0.0d0
           do I=2,NRHO
             A=SN(I)/(Q*RHO(I))
             A=15.0d0*(A-CS(I))/(Q*RHO(I))**2 - 6*A + CS(I)
             F(I)=A*WP(I,3)*VP(I,3)
           end do
           D=P3*SIMP(NRHO,F,DRHO)/Q

           lcount = lcount-1
           vnl(k1,lcount)=D*GY*(3.0d0*(1.0d0-GZ*GZ)-4.0d0*GY*GY)
     >                          /dsqrt(24.0d0)
           lcount = lcount-1
           vnl(k1,lcount)=D*GX*GY*GZ
           lcount = lcount-1
           vnl(k1,lcount)=D*GY*(5.0d0*GZ*GZ-1.0d0)
     >                          /dsqrt(40.0d0)
           lcount = lcount-1
           vnl(k1,lcount)=D*GZ*(5.0d0*GZ*GZ-3.0d0)
     >                          /dsqrt(60.0d0)
           lcount = lcount-1
           vnl(k1,lcount)=D*GX*(5.0d0*GZ*GZ-1.0d0)
     >                          /dsqrt(40.0d0)
           lcount = lcount-1
           vnl(k1,lcount)=D*GZ*(GX*GX - GY*GY)
     >                          /2.0d0
           lcount = lcount-1
           vnl(k1,lcount)=D*GX*(4.0d0*GX*GX-3.0d0*(1.0d0-GZ*GZ))
     >                          /dsqrt(24.0d0)
c           lcount = lcount-1
c           vnl(k1,lcount)=D*GX*(4.0d0*GX*GX-3.0d0*(1.0d0-GZ*GZ))
c     >                          /dsqrt(24.0d0)
c           lcount = lcount-1
c           vnl(k1,lcount)=D*GY*(3.0d0*(1.0d0-GZ*GZ)-4.0d0*GY*GY)
c     >                          /dsqrt(24.0d0)
c           lcount = lcount-1
c           vnl(k1,lcount)=D*GZ*(GX*GX - GY*GY)
c     >                          /2.0d0
c           lcount = lcount-1
c           vnl(k1,lcount)=D*GX*GY*GZ
c           lcount = lcount-1
c           vnl(k1,lcount)=D*GX*(5.0d0*GZ*GZ-1.0d0)
c     >                          /dsqrt(40.0d0)
c           lcount = lcount-1
c           vnl(k1,lcount)=D*GY*(5.0d0*GZ*GZ-1.0d0)
c     >                          /dsqrt(40.0d0)
c           lcount = lcount-1
c           vnl(k1,lcount)=D*GZ*(5.0d0*GZ*GZ-3.0d0)
c     >                          /dsqrt(60.0d0)
        end if



*::::::::::::::::::::::::::::::  d-wave  ::::::::::::::::::::::::::::::
  300   CONTINUE
        if (locp.ne.2) then
          F(1)=0.0d0
          DO I=2,NRHO
            A=3.0d0*(SN(I)/(Q*RHO(I))-CS(I))/(Q*RHO(I))-SN(I)
            F(I)=A*WP(I,2)*VP(I,2)
          END DO
          D=P2*SIMP(NRHO,F,DRHO)/Q

          lcount = lcount-1
          vnl(k1,lcount)=D*GX*GY
          lcount = lcount-1
          vnl(k1,lcount)=D*GY*GZ
          lcount = lcount-1
          vnl(k1,lcount)=D*(3.0d0*GZ*GZ-1.0d0)
     >                          /(2.0d0*dsqrt(3.0d0))
          lcount = lcount-1
          vnl(k1,lcount)=D*GZ*GX
          lcount = lcount-1
          vnl(k1,lcount)=D*(GX*GX-GY*GY)/(2.0d0)

c          lcount = lcount-1
c          vnl(k1,lcount)=D*(3.0d0*GZ*GZ-1.0d0)
c     >                          /(2.0d0*dsqrt(3.0d0))
c          lcount = lcount-1
c          vnl(k1,lcount)=D*GX*GY
c          lcount = lcount-1
c          vnl(k1,lcount)=D*GY*GZ
c          lcount = lcount-1
c          vnl(k1,lcount)=D*GZ*GX
c          lcount = lcount-1
c          vnl(k1,lcount)=D*(GX*GX-GY*GY)/(2.0d0)
        end if

*::::::::::::::::::::::::::::::  p-wave  ::::::::::::::::::::::::::::::
  400   CONTINUE
        if (locp.ne.1) then
           F(1)=0.0d0
           DO I=2,NRHO
             F(I)=(SN(I)/(Q*RHO(I))-CS(I))*WP(I,1)*VP(I,1)
           END DO
           P=P1*SIMP(NRHO,F,DRHO)/Q
           lcount = lcount-1
           vnl(k1,lcount)=P*GY
           lcount = lcount-1
           vnl(k1,lcount)=P*GZ
           lcount = lcount-1
           vnl(k1,lcount)=P*GX
        end if

*::::::::::::::::::::::::::::::  s-wave  :::::::::::::::::::::::::::::::
  500   CONTINUE
        if (locp.ne.0) then
          DO I=1,NRHO
            F(I)=SN(I)*WP(I,0)*VP(I,0)
          END DO
          lcount = lcount-1
          vnl(k1,lcount)=P0*SIMP(NRHO,F,DRHO)/Q
        end if

*::::::::::::::::::::::::::::::  local  :::::::::::::::::::::::::::::::
  600   CONTINUE


        if (version.eq.3) then
        DO  I=1,NRHO
          F(I)=RHO(I)*VP(I,locp)*SN(I)
        END DO
        vl(k1)=SIMP(NRHO,F,DRHO)*FORPI/Q-ZV*FORPI*CS(NRHO)/(Q*Q)
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
        vl(k1)=SIMP(NRHO,F,DRHO)*FORPI/Q
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
           rho_sc_k(k1,1) = SIMP(nrho,f,drho)*forpi/Q

           do i=1,nrho
             f(i)=(sn(i)/(Q*rho(i))-cs(i))*rho_sc_r(i,2)*rho(i)
           end do
           P = SIMP(nrho,f,drho)*forpi/Q
           rho_sc_k(k1,2)=P*GX
           rho_sc_k(k1,3)=P*GY
           rho_sc_k(k1,4)=P*GZ

        end if
    
  700 CONTINUE
      call D1dB_Vector_SumAll(4*nfft3d,rho_sc_k)
      call D1dB_Vector_SumAll(nfft3d,vl)
      call D1dB_Vector_Sumall(lmmax*nfft3d,vnl)


*:::::::::::::::::::::::::::::::  G=0  ::::::::::::::::::::::::::::::::      
      if (pzero.eq.taskid_i) then

         if (version.eq.3) then
         DO I=1,NRHO
           F(I)=VP(I,locp)*RHO(I)**2
         END DO
         vl(zero)=FORPI*SIMP(NRHO,F,DRHO)+TWOPI*ZV*RHO(NRHO)**2
         end if

         if (version.eq.4) then
         if (fast_erf) then
            do I=1,NRHO
              xerf=RHO(I)/rlocal
              yerf = (1.0d0
     >             + xerf*(c1 + xerf*(c2
     >             + xerf*(c3 + xerf*(c4
     >             + xerf*(c5 + xerf*c6))))))**4
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
         vl(zero)=FORPI*SIMP(NRHO,F,DRHO)
         end if

*        **** semicore density ****
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
            rho_sc_k(zero,1) = forpi*SIMP(nrho,f,drho)
            rho_sc_k(zero,2) = 0.0d0
            rho_sc_k(zero,3) = 0.0d0
            rho_sc_k(zero,4) = 0.0d0
         end if

         do l=1,lmmax
           vnl(zero,l)=0.0d0
         end do
*        *** only j0 is non-zero at zero ****
         if (locp.ne.0) then
            DO  I=1,NRHO
              F(I)=RHO(I)*WP(I,0)*VP(I,0)
            END DO
            vnl(zero,1)=P0*SIMP(NRHO,F,DRHO)
         end if

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
          b_prj(lcount) = 4
    
          lcount = lcount-1
          n_prj(lcount) = 1
          l_prj(lcount) = 3
          m_prj(lcount) = -2
          b_prj(lcount) = 4
           
          lcount = lcount-1
          n_prj(lcount) = 1
          l_prj(lcount) = 3
          m_prj(lcount) = -1
          b_prj(lcount) = 4
     
          lcount = lcount-1
          n_prj(lcount) = 1
          l_prj(lcount) = 3
          m_prj(lcount) = 0
          b_prj(lcount) = 4
           
          lcount = lcount-1
          n_prj(lcount) = 1
          l_prj(lcount) = 3
          m_prj(lcount) = 1
          b_prj(lcount) = 4
           
          lcount = lcount-1
          n_prj(lcount) = 1
          l_prj(lcount) = 3
          m_prj(lcount) = 2
          b_prj(lcount) = 4
     
          lcount = lcount-1
          n_prj(lcount) = 1
          l_prj(lcount) = 3
          m_prj(lcount) = 3
          b_prj(lcount) = 4
        end if


        !::::  d-wave  ::::
  930   CONTINUE
        if (locp.ne.2) then
          lcount = lcount-1
          n_prj(lcount) = 1
          l_prj(lcount) = 2
          m_prj(lcount) = -2
          b_prj(lcount) = 3

          lcount = lcount-1
          n_prj(lcount) = 1
          l_prj(lcount) = 2
          m_prj(lcount) = -1
          b_prj(lcount) = 3
          
          lcount = lcount-1
          n_prj(lcount) = 1
          l_prj(lcount) = 2
          m_prj(lcount) = 0
          b_prj(lcount) = 3
          
          lcount = lcount-1
          n_prj(lcount) = 1
          l_prj(lcount) = 2
          m_prj(lcount) = 1
          b_prj(lcount) = 3
          
          lcount = lcount-1
          n_prj(lcount) = 1
          l_prj(lcount) = 2
          m_prj(lcount) = 2
          b_prj(lcount) = 3
        end if


        !::::  p-wave  ::::
  940   CONTINUE
        if (locp.ne.1) then
          lcount = lcount-1
          n_prj(lcount) = 1
          l_prj(lcount) = 1
          m_prj(lcount) = -1
          b_prj(lcount) = 2

          lcount = lcount-1
          n_prj(lcount) = 1
          l_prj(lcount) = 1
          m_prj(lcount) = 0
          b_prj(lcount) = 2

          lcount = lcount-1
          n_prj(lcount) = 1
          l_prj(lcount) = 1
          m_prj(lcount) = 1
          b_prj(lcount) = 2
        end if


        !::::  s-wave  ::::
  950   CONTINUE
        if (locp.ne.0) then
          lcount = lcount-1
          n_prj(lcount) = 1
          l_prj(lcount) = 0
          m_prj(lcount) = 0
          b_prj(lcount) = 1
        end if

      IERR=0
      RETURN
      END



