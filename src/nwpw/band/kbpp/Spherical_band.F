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

        if (dabs(Q).lt.1.0d-12) then
           gx=0.0d0
           gy=0.0d0
           gz=0.0d0
        else
           gx=gx/Q
           gy=gy/Q
           gz=gz/Q
        end if

        gsin=dsqrt(gx*gx+gy*gy)
        if (dabs(gsin).lt.1.e-12) then
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
      call C3dB_Vector_Sumall(nfft3d,Ylm)


      return
      end



