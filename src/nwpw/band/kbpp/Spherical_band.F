      subroutine Spherical_band(l,m,kvec,
     >                    nfft1,nfft2,nfft3,
     >                    G,
     >                    Ylm)
      implicit none
      integer          l,m
      double precision kvec(3)

      integer nfft1,nfft2,nfft3
      double precision G(nfft1,nfft2,nfft3,3)
      complex*16 Ylm(nfft1,nfft2,nfft3)

      integer np,taskid,MASTER
      parameter (MASTER=0)

*     *** local variables ****
      integer task_count,nfft3d
      integer k1,k2,k3
      real*8 pi,twopi,forpi
      real*8 psin,pcos,gsin
      real*8 gx,gy,gz,q,tmp
      real*8 scof(7,0:3)
      complex*16 ihat,expm

      call Parallel_np(np)
      call Parallel_taskid(taskid)

      nfft3d = nfft1*nfft2*nfft3

      pi=4.0d0*datan(1.0d0)
      twopi=2.0d0*pi
      forpi=4.0d0*pi 
      ihat=dcmplx(0.0d0,1.0d0)
      scof(1,0)=dsqrt(1.0d0/forpi)
      scof(1,1)=dsqrt(0.75d0/pi)
      scof(2,1)=dsqrt(0.75/2.0d0/forpi)
      scof(3,1)= -scof(1,1)
      scof(1,2)=dsqrt(15.0d0/32.0d0/pi)
      scof(2,2)=dsqrt(15.0d0/8.0d0/pi)     
      scof(3,2)=dsqrt(5.0d0/16.0d0/pi)
      scof(4,2)=-scof(2,2)
      scof(5,2)=scof(2,1)
      scof(1,3)=dsqrt(35.0d0/64.0d0/pi)
      scof(2,3)=dsqrt(105.d0/32.0d0/pi)
      scof(3,3)=dsqrt(21.0d0/64.0d0/pi)
      scof(4,3)=dsqrt(7.0d0/16.0d0/pi)
      scof(5,3)= -scof(3,3)
      scof(6,3)=scof(3,2)
      scof(7,3)= -scof(3,1)

*======================  Fourier transformation  ======================
      call dcopy(2*nfft3d,0.0d0,0,Ylm,1)
      task_count = -1
      DO 700 k3=1,nfft3
      DO 700 k2=1,nfft2
      DO 700 k1=1,nfft1
        task_count = task_count + 1
        if (mod(task_count,np).ne.taskid) go to 700

        gx = G(k1,k2,k3,1) + kvec(1)
        gy = G(k1,k2,k3,2) + kvec(2)
        gz = G(k1,k2,k3,3) + kvec(3)
        Q=DSQRT(gx**2 + gy**2 + gz**2)

        if (dabs(Q).lt.1.0d-9) then
           gx=0.0d0
           gy=0.0d0
           gz=0.0d0
        else
           gx=gx/Q
           gy=gy/Q
           gz=gz/Q
        end if

        gsin=dsqrt(gx*gx+gy*gy)
        if (dabs(gsin).lt.1.e-9) then
           psin=0.0d0
           pcos=1.0d0
        else
           psin=gy/gsin
           pcos=gx/gsin
        end if

*::::::::::::::::::::::::::::::  f-wave  ::::::::::::::::::::::::::::::

        if ((l.eq.3).and.(m.eq.-3)) then 
          if (pcos.lt.1.e-15) then 
            tmp=pi*0.5
          else
           tmp=datan(psin/pcos)
          end if
          expm=dcmplx(dcos(3.0d0*tmp),-dsin(3.0d0*tmp))
          Ylm(k1,k2,k3)=scof(1,3)*gsin*gsin*gsin*expm
        end if
        if ((l.eq.3).and.(m.eq.-2)) then
           expm=dcmplx(pcos*pcos-psin*psin,-2.0d0*pcos*psin) 
           Ylm(k1,k2,k3)=scof(2,3)*gz*gsin*gsin*expm
        end if  
        if ((l.eq.3).and.(m.eq.-1)) then
           expm=dcmplx(pcos,-psin)
           Ylm(k1,k2,k3)=scof(3,3)*gsin*(5.0d0*gz*gz-1.0d0)*expm
        end if           

        if ((l.eq.3).and.(m.eq.0)) then
           Ylm(k1,k2,k3)=scof(4,3)*(5.0d0*gz*gz*gz-3.0d0*gz)
        end if
    
        if ((l.eq.3).and.(m.eq.1)) then 
           expm=dcmplx(pcos,psin)
           Ylm(k1,k2,k3)=scof(5,3)*gsin*(5.0d0*gz*gz-1.0d0)*expm
        end if

        if ((l.eq.3).and.(m.eq.2)) then
           expm=dcmplx(pcos*pcos-psin*psin,2.0d0*pcos*psin) 
           Ylm(k1,k2,k3)=scof(6,3)*expm*gsin*gsin*gz 
        end if

        if ((l.eq.3).and.(m.eq.3)) then 
          if (pcos.lt.1.e-15) then 
            tmp=pi*0.5
          else
           tmp=datan(psin/pcos)
          end if
          expm=dcmplx(dcos(3.0d0*tmp),dsin(3.0d0*tmp))
          Ylm(k1,k2,k3)=scof(7,3)*gsin*gsin*gsin*expm
        end if


*::::::::::::::::::::::::::::::  d-wave  ::::::::::::::::::::::::::::::


        if ((l.eq.2).and.(m.eq.-2)) then
          expm=dcmplx(pcos*pcos-psin*psin,-psin*pcos*2.0d0) 
          Ylm(k1,k2,k3)=scof(1,2)*expm*gsin*gsin
        end if

        if ((l.eq.2).and.(m.eq.-1)) then
          expm=dcmplx(pcos,-psin)
          Ylm(k1,k2,k3)=scof(2,2)*expm*gsin*gz
        end if

        if ((l.eq.2).and.(m.eq.0)) then
          Ylm(k1,k2,k3)=scof(3,2)*(3.0d0*gz*gz-1.0d0)
        end if

        if ((l.eq.2).and.(m.eq.1)) then
          expm=dcmplx(pcos,psin)
          Ylm(k1,k2,k3)=scof(4,2)*expm*gsin*gz
        end if

        if ((l.eq.2).and.(m.eq.2)) then
          expm=dcmplx(pcos*pcos-psin*psin,2.0d0*pcos*psin)
          Ylm(k1,k2,k3)=scof(5,2)*expm*gsin*gsin
        end if

*::::::::::::::::::::::::::::::  p-wave  ::::::::::::::::::::::::::::::

        if ((l.eq.1).and.(m.eq.-1)) 
     >    Ylm(k1,k2,k3)=scof(1,1)*gsin*dcmplx(pcos,-psin)

        if ((l.eq.1).and.(m.eq.0)) 
     >    Ylm(k1,k2,k3)=scof(2,1)*gz

        if ((l.eq.1).and.(m.eq.1)) 
     >    Ylm(k1,k2,k3)=scof(3,1)*gsin*dcmplx(pcos,psin)


*::::::::::::::::::::::::::::::  s-wave  :::::::::::::::::::::::::::::::

        if ((l.eq.0).and.(m.eq.0)) 
     >    Ylm(k1,k2,k3)=scof(1,0)
  700 CONTINUE
      call Parallel_Vector_Sumall(nfft3d,Ylm)
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Generate the value of the spherical harmonic for
C  the unit cartesian vector (x,y,z)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      complex*16 function sphharm2(L,M,x,y,z)
      implicit none
      integer L,M,I,AM
      real*8 x,y,z,pleg
      real*8 plm,fact1,fact2,sim,sre,PHI
      external pleg
cccccccccccccccccccccccccccccc
      AM=ABS(M)
      if (AM.GT.L) THEN
        SPHHARM2=DCMPLX(0.0D0,0.0D0)
        GOTO 777
      END IF   
      plm=Pleg(z,L,AM)
      fact1=1.0d0
      do i=2,(l-AM)
        fact1=fact1*dble(i)
      end do
      fact1=fact1*(2*l+1)
      fact2=1.0d0
      do i=2,(l+AM)
        fact2=fact2*dble(i)
      end do
      fact2=fact2*16.0d0*datan(1.0d0)
      sim=dsqrt(fact1/fact2)
      PHI=DBLE(M)*DATAN2(Y,X)
      sre=sim*dcos(phi)*plm
      sim=sim*dsin(phi)*plm
      sphharm2=dcmplx(sre,sim)
777   return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c calculates the Associated Legendre Function
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function pleg(x,L,M)
      implicit none
#include "errquit.fh"
      real*8 x
      integer L,M
ccccccc locals cccccccccc
      real*8 pmm,pmmp1,PLL,SX2,FACT
      INTEGER LL,J
      if (m.gt.L) then
        CALL ERRQUIT("m>l in call assoc. legendre function",0,0)
      end if 
      IF ((X.GT.1.0D0).OR.(X.LT.-1.D0)) THEN
        CALL ERRQUIT("bad x in call assoc. legendre function",0,0)
      END IF
      sx2=(1.0d0-x)*(1.0d0+x)
      if (sx2.lt.0.0d0) sx2=0.0d0
      sx2=dsqrt(sx2)
      pmm=1.0d0
      pll=0.0d0
      fact=1.0d0
      do j=1,m
        pmm = -pmm*fact*sx2
        fact = fact +2.0d0
      end do
      if (l.eq.m) then
        pleg=pmm
        goto 777
      else
        pmmp1=x*pmm*dble(2*m+1)
        if (l.eq.(m+1)) then
          pleg=pmmp1
          GOTO 777
        else
          do ll=m+2,L
            pll=x*dble(2*ll-1)*pmmp1
            pll=pll-dble(ll+m-1)*pmm
            pll=pll/dble(ll-m)
            pmm=pmmp1
            pmmp1=pll
          end do
          pleg=pll
        end if
      end if
777   return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccc
      complex*16 function sphharm(l,m,gx,gy,gz)
      implicit none
      integer l,m
      real*8 gx,gy,gz     
*     *** local variables ****
      real*8 pi,twopi,forpi
      real*8 psin,pcos,gsin
      real*8 q,tmp
      real*8 scof(7,0:3)
      complex*16 ihat,expm
      integer am

      am=abs(m)
      if (am.gt.l) then
        sphharm=dcmplx(0.d0,0.d0)
        return
      end if
      pi=4.0d0*datan(1.0d0)
      twopi=2.0d0*pi
      forpi=4.0d0*pi 
      ihat=dcmplx(0.0d0,1.0d0)
      scof(1,0)=dsqrt(1.0d0/forpi)
      scof(1,1)=dsqrt(0.75d0/pi)
      scof(2,1)=dsqrt(0.75/2.0d0/forpi)
      scof(3,1)= -scof(1,1)
      scof(1,2)=dsqrt(15.0d0/32.0d0/pi)
      scof(2,2)=dsqrt(15.0d0/8.0d0/pi)     
      scof(3,2)=dsqrt(5.0d0/16.0d0/pi)
      scof(4,2)=-scof(2,2)
      scof(5,2)=scof(2,1)
      scof(1,3)=dsqrt(35.0d0/64.0d0/pi)
      scof(2,3)=dsqrt(105.d0/32.0d0/pi)
      scof(3,3)=dsqrt(21.0d0/64.0d0/pi)
      scof(4,3)=dsqrt(7.0d0/16.0d0/pi)
      scof(5,3)= -scof(3,3)
      scof(6,3)=scof(3,2)
      scof(7,3)= -scof(3,1)

*======================  Fourier transformation  ======================
        gsin=dsqrt(gx*gx+gy*gy)
        if (dabs(gsin).lt.1.e-9) then
           psin=0.0d0
           pcos=1.0d0
        else
           psin=gy/gsin
           pcos=gx/gsin
        end if

*::::::::::::::::::::::::::::::  f-wave  ::::::::::::::::::::::::::::::
      if (l.eq.3) then
        if ((m.eq.-3)) then 
          if (pcos.lt.1.e-15) then 
            tmp=pi*0.5
          else
            tmp=datan(psin/pcos)
          end if
          expm=dcmplx(dcos(3.0d0*tmp),-dsin(3.0d0*tmp))
          sphharm=scof(1,3)*gsin*gsin*gsin*expm
          return
        end if
        if ((m.eq.-2)) then
           expm=dcmplx(pcos*pcos-psin*psin,-2.0d0*pcos*psin) 
           sphharm=scof(2,3)*gz*gsin*gsin*expm
           return
        end if  
        if ((m.eq.-1)) then
           expm=dcmplx(pcos,-psin)
           sphharm=scof(3,3)*gsin*(5.0d0*gz*gz-1.0d0)*expm
           return
        end if           

        if ((m.eq.0)) then
           sphharm=scof(4,3)*(5.0d0*gz*gz*gz-3.0d0*gz)
           return
        end if
    
        if ((m.eq.1)) then 
           expm=dcmplx(pcos,psin)
           sphharm=scof(5,3)*gsin*(5.0d0*gz*gz-1.0d0)*expm
           return
        end if

        if ((m.eq.2)) then
           expm=dcmplx(pcos*pcos-psin*psin,2.0d0*pcos*psin) 
           sphharm=scof(6,3)*expm*gsin*gsin*gz 
           return
        end if

        if ((m.eq.3)) then 
          if (pcos.lt.1.e-15) then 
            tmp=pi*0.5
          else
           tmp=datan(psin/pcos)
          end if
          expm=dcmplx(dcos(3.0d0*tmp),dsin(3.0d0*tmp))
          sphharm=scof(7,3)*gsin*gsin*gsin*expm
          return
        end if
      end if

*::::::::::::::::::::::::::::::  d-wave  ::::::::::::::::::::::::::::::

      if (l.eq.2) then
        if ((m.eq.-2)) then
          expm=dcmplx(pcos*pcos-psin*psin,-psin*pcos*2.0d0) 
          sphharm=scof(1,2)*expm*gsin*gsin
          return
        end if

        if ((m.eq.-1)) then
          expm=dcmplx(pcos,-psin)
          sphharm=scof(2,2)*expm*gsin*gz
          return
        end if

        if ((m.eq.0)) then
          sphharm=scof(3,2)*(3.0d0*gz*gz-1.0d0)
          return
        end if

        if ((m.eq.1)) then
          expm=dcmplx(pcos,psin)
          sphharm=scof(4,2)*expm*gsin*gz
          return
        end if

        if ((m.eq.2)) then
          expm=dcmplx(pcos*pcos-psin*psin,2.0d0*pcos*psin)
          sphharm=scof(5,2)*expm*gsin*gsin
          return
        end if
      end if
*::::::::::::::::::::::::::::::  p-wave  ::::::::::::::::::::::::::::::
      if (l.eq.1) then
        if ((m.eq.-1)) then
          sphharm=scof(1,1)*gsin*dcmplx(pcos,-psin)
          return
        end if
        if ((m.eq.0)) then
          sphharm=scof(2,1)*gz
          return
        end if
        if ((m.eq.1)) then
          sphharm=scof(3,1)*gsin*dcmplx(pcos,psin)
          return
        end if
      end if  
*::::::::::::::::::::::::::::::  s-wave  :::::::::::::::::::::::::::::::
      if ((l.eq.0).and.(m.eq.0)) then 
          sphharm=scof(1,0)
          return
      end if
      write(*,*)"Error did not find l or m in list"
      write(*,*)"l = ",l," m= ",m
      call flush(6)
      sphharm=dcmplx(0.0d0,0.0d0)
      return 
      end
c $Id$
