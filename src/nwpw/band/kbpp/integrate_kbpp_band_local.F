*
* $Id$
*

      subroutine integrate_kbpp_band_local(version,
     >                            nrho,drho,lmax,locp,zv,
     >                            vp,wp,rho,f,cs,sn,
     >                            nfft1,nfft2,nfft3,lmmax,
     >                            G,vl,
     >                            n_prj,l_prj,m_prj,vnlnrm,
     >                            semicore,rho_sc_r,rho_sc_k,
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

      integer ierr

      integer np,taskid,MASTER
      parameter (MASTER=0)

*     *** local variables ****
      integer lcount,task_count,nfft3d
      integer k1,k2,k3,i,l
      double precision pi,twopi,forpi
      double precision p0,p1,p2,p3,p
      double precision gx,gy,gz,a,q,d
      logical small_cell


*     **** external functions ****
      double precision dsum,simp
      external         dsum,simp
      logical          control_psp_semicore_small
      external         control_psp_semicore_small

#include "bafdecls.fh"
      logical value

      call Parallel_np(np)
      call Parallel_taskid(taskid)

      small_cell = control_psp_semicore_small()

      nfft3d = (nfft1)*nfft2*nfft3
      pi=4.0d0*datan(1.0d0)
      twopi=2.0d0*pi
      forpi=4.0d0*pi

      IF((NRHO/2)*2.EQ.NRHO) THEN
        WRITE(*,*)"local pseudopotential not computed: nrho is not odd"
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

*======================  Fourier transformation  ======================
      call dcopy(nfft3d,0.0d0,0,VL,1)
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

        if ((k1.eq.1).and.(k2.eq.1).and.(k3.eq.1)) go to 700

        
        GX=G(k1,k2,k3,1)/Q
        GY=G(k1,k2,k3,2)/Q
        GZ=G(k1,k2,k3,3)/Q
        DO I=1,NRHO
          CS(I)=DCOS(Q*RHO(I))
          SN(I)=DSIN(Q*RHO(I))
        END DO

*::::::::::::::::::::::::::::::  local  :::::::::::::::::::::::::::::::
  600   CONTINUE

        DO  I=1,NRHO
          F(I)=RHO(I)*VP(I,locp)*SN(I)
        END DO
        vl(k1,k2,k3)=SIMP(nrho,f,drho)*forpi/Q-ZV*forpi*cs(nrho)/(Q*Q)
 

*::::::::::::::::::::: semicore density :::::::::::::::::::::::::::::::
        if (semicore) then
           if (small_cell) then
              do i=1,nrho
                 f(i) = rho(i)*(rho_sc_r(i,1))*sn(i)
              end do
           else
              do i=1,nrho
                 f(i) = rho(i)*dsqrt(rho_sc_r(i,1))*sn(i)
              end do
           end if
           rho_sc_k(k1,k2,k3,1) = SIMP(nrho,f,drho)*forpi/Q

           do i=1,nrho
             f(i)=(sn(i)/(Q*rho(i))-cs(i))*rho_sc_r(i,2)*rho(i)
           end do
           P = SIMP(nrho,f,drho)*forpi/Q
           rho_sc_k(k1,k2,k3,2)=P*GX
           rho_sc_k(k1,k2,k3,3)=P*GY
           rho_sc_k(k1,k2,k3,4)=P*GZ

        end if
    
  700 CONTINUE

      call Parallel_Vector_SumAll(4*nfft3d,rho_sc_k)
      call Parallel_Vector_SumAll(nfft3d,VL)

*:::::::::::::::::::::::::::::::  G=0  ::::::::::::::::::::::::::::::::      

      DO I=1,NRHO
        F(I)=VP(I,locp)*RHO(I)**2
      END DO
      VL(1,1,1)=FORPI*SIMP(NRHO,F,DRHO)+TWOPI*ZV*RHO(NRHO)**2

*     **** semicore density ****
      if (semicore) then
         if (small_cell) then
            do i=1,nrho
               f(i) = (rho_sc_r(i,1))*rho(i)**2
            end do
         else
            do i=1,nrho
               f(i) = dsqrt(rho_sc_r(i,1))*rho(i)**2
            end do
         end if
         rho_sc_k(1,1,1,1) = forpi*SIMP(nrho,f,drho)
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

