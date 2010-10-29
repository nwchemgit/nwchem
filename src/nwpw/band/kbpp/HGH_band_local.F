*
* $Id$
*


      subroutine HGH_band_local(version,
     >                     zv,rloc,C1,C2,C3,C4,
     >                     nfft1,nfft2,nfft3,
     >                     G,
     >                     vl)
      implicit none
      integer          version
      double precision zv,rloc,C1,C2,C3,C4

      integer nfft1,nfft2,nfft3
      double precision G(nfft1,nfft2,nfft3,3)
      double precision vl(nfft1,nfft2,nfft3)


      integer np,taskid,MASTER
      parameter (MASTER=0)

*     *** local variables ****
      integer task_count,nfft3d
      integer k1,k2,k3
      double precision pi,twopi,forpi
      double precision b0,b1
      double precision x,x2,x4,x6,a,f,Q


*     **** Error function parameters ****


      call Parallel_np(np)
      call Parallel_taskid(taskid)
      nfft3d = nfft1*nfft2*nfft3
      pi=4.0d0*datan(1.0d0)
      twopi=2.0d0*pi
      forpi=4.0d0*pi
      b0 = -forpi*zv
      b1 = dsqrt(8.0d0*pi*pi*pi)*(rloc*rloc*rloc)
      
*======================  Fourier transformation  ======================
      call dcopy(nfft3d,0.0d0,0,vl,1)
      task_count = -1
      DO 700 k3=1,nfft3
      DO 700 k2=1,nfft2
      DO 700 k1=1,nfft1
        task_count = task_count + 1
        if (mod(task_count,np).ne.taskid) go to 700

        if ((k1.eq.1).and.(k2.eq.1).and.(k3.eq.1)) go to 700


        Q=DSQRT(G(k1,k2,k3,1)**2
     >         +G(k1,k2,k3,2)**2
     >         +G(k1,k2,k3,3)**2)

        x = Q*rloc
        x2 = x*x
        x4 = x2*x2
        x6 = x4*x2
        a = dexp(-0.5d0*x2)

        f = b0*a/(Q*Q)
     >    + b1*a
     >        *( C1 
     >         + C2*(3.0d0-x2)
     >         + C3*(15.0d0 - 10.0d0*x2 + x4)
     >         + C4*(105.0d0 - 105.0d0*x2 + 21.0d0*x4 - x6))

        vl(k1,k2,k3)= f
    
  700 continue


      call Parallel_Vector_SumAll(nfft3d,vl)


      vl(1,1,1)=b1*(C1+C2*3.0d0+C3*15.0D0+C4*105.0D0)-b0*0.5d0*rloc*rloc
      
      return
      end



