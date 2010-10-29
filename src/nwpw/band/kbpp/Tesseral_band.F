*
* $Id$
*
      subroutine Tesseral_band(l,m,kvec,
     >                    nfft1,nfft2,nfft3,
     >                    G,
     >                    Ylm)
      implicit none
      integer          l,m
      double precision kvec(3)

      integer nfft1,nfft2,nfft3
      double precision G(nfft1,nfft2,nfft3,3)
      double precision Ylm(nfft1,nfft2,nfft3)

      integer np,taskid,MASTER
      parameter (MASTER=0)

*     *** local variables ****
      integer task_count,nfft3d
      integer k1,k2,k3
      double precision pi,twopi,forpi
      double precision p0,p1,p2,p3
      double precision gx,gy,gz,q



      call Parallel_np(np)
      call Parallel_taskid(taskid)

      nfft3d = nfft1*nfft2*nfft3
      pi=4.0d0*datan(1.0d0)
      twopi=2.0d0*pi
      forpi=4.0d0*pi

      P0=DSQRT(1.0d0/FORPI)
      P1=DSQRT(3.0d0/FORPI)
      P2=DSQRT(15.0d0/FORPI)
      P3=DSQRT(105.0d0/FORPI)


*======================  Fourier transformation  ======================
      call dcopy(nfft3d,0.0d0,0,Ylm,1)
      task_count = -1
      DO 700 k3=1,nfft3
      DO 700 k2=1,nfft2
      DO 700 k1=1,nfft1
        task_count = task_count + 1
        if (mod(task_count,np).ne.taskid) go to 700

        gx = G(k1,k2,k3,1) + kvec(1)
        gy = G(k1,k2,k3,2) + kvec(2)
        gz = G(k1,k2,k3,3) + kvec(3)
        Q=DSQRT(gx*gx + gy*gy + gz*gz)

        if (dabs(Q).lt.1.0d-12) then
	    Ylm(k1,k2,k3)=0.0d0
            if (l.eq.0) Ylm(k1,k2,k3)=P0
            goto 700
        end if
        gx=gx/Q
        gy=gy/Q
        gz=gz/Q

*::::::::::::::::::::::::::::::  f-wave  ::::::::::::::::::::::::::::::
C           if ((l.eq.3).and.(m.eq.-3)) 
C     >       Ylm(k1,k2,k3)=P3*GX*(4.0d0*GX*GX-3.0d0*(1.0d0-GZ*GZ))
C     >                          /dsqrt(24.0d0)

           if ((l.eq.3).and.(m.eq.-3)) 
     >       Ylm(k1,k2,k3)=P3*GX*(GX*GX-GY*GY)
     >                          /dsqrt(24.0d0)

        
C           if ((l.eq.3).and.(m.eq.-2)) 
C     >     Ylm(k1,k2,k3)=P3*GY*(3.0d0*(1.0d0-GZ*GZ)-4.0d0*GY*GY)
C     >                          /dsqrt(24.0d0)
         
            if ((l.eq.3).and.(m.eq.-2)) 
     >     Ylm(k1,k2,k3)=P3*GY*(GX*GX-GY*GY)
     >                          /dsqrt(24.0d0)
         
           if ((l.eq.3).and.(m.eq.-1)) 
     >     Ylm(k1,k2,k3)=P3*GZ*(GX*GX - GY*GY)
     >                          /2.0d0
          
           if ((l.eq.3).and.(m.eq.0)) 
     >     Ylm(k1,k2,k3)=P3*GX*GY*GZ
           
c           if ((l.eq.3).and.(m.eq.1)) 
c     >     Ylm(k1,k2,k3)=P3*GX*(5.0d0*GZ*GZ-1.0d0)
c     >                          /dsqrt(40.0d0)

           if ((l.eq.3).and.(m.eq.1)) 
     >     Ylm(k1,k2,k3)=P3*GX*(4.0d0*GZ*GZ-GY*GY-GX*GX)
     >                          /dsqrt(40.0d0)


C           if ((l.eq.3).and.(m.eq.2)) 
C     >     Ylm(k1,k2,k3)=P3*GY*(5.0d0*GZ*GZ-1.0d0)
C     >                          /dsqrt(40.0d0)

            if ((l.eq.3).and.(m.eq.1)) 
     >      Ylm(k1,k2,k3)=P3*GY*(4.0d0*GZ*GZ-GY*GY-GX*GX)
     >                          /dsqrt(40.0d0)

c           if ((l.eq.3).and.(m.eq.3)) 
c     >     Ylm(k1,k2,k3)=P3*GZ*(5.0d0*GZ*GZ-3.0d0)
c     >                          /dsqrt(60.0d0)

            if ((l.eq.3).and.(m.eq.3)) 
     >      Ylm(k1,k2,k3)=P3*GZ*(2.0d0*GZ*GZ-3.0d0*GY*GY-3.0d0*GX*GX)
     >                          /dsqrt(60.0d0)



*::::::::::::::::::::::::::::::  d-wave  ::::::::::::::::::::::::::::::
c           if ((l.eq.2).and.(m.eq.-2)) 
c     >     Ylm(k1,k2,k3)=P2*(3.0d0*GZ*GZ-1.0d0)
c     >                          /(2.0d0*dsqrt(3.0d0))
           if ((l.eq.2).and.(m.eq.-2)) 
     >     Ylm(k1,k2,k3)=P2*(2.0d0*GZ*GZ-GX*GX-GY*GY)
     >                          /(2.0d0*dsqrt(3.0d0))

           if ((l.eq.2).and.(m.eq.-1)) 
     >     Ylm(k1,k2,k3)=P2*GX*GY

           if ((l.eq.2).and.(m.eq.0)) 
     >     Ylm(k1,k2,k3)=P2*GY*GZ

           if ((l.eq.2).and.(m.eq.1)) 
     >     Ylm(k1,k2,k3)=P2*GZ*GX

           if ((l.eq.2).and.(m.eq.2)) 
     >     Ylm(k1,k2,k3)=P2*(GX*GX-GY*GY)/(2.0d0)


*::::::::::::::::::::::::::::::  p-wave  ::::::::::::::::::::::::::::::
           if ((l.eq.1).and.(m.eq.-1)) 
     >     Ylm(k1,k2,k3)=P1*GX

           if ((l.eq.1).and.(m.eq.0)) 
     >     Ylm(k1,k2,k3)=P1*GY

           if ((l.eq.1).and.(m.eq.1)) 
     >     Ylm(k1,k2,k3)=P1*GZ


*::::::::::::::::::::::::::::::  s-wave  :::::::::::::::::::::::::::::::
           if ((l.eq.0).and.(m.eq.0)) 
     >     Ylm(k1,k2,k3)=P0



  700 CONTINUE
      call Parallel_Vector_Sumall(nfft3d,Ylm)


      return
      end



