*
* $Id: HGH_band_local.F,v 1.4 2007-02-01 01:52:05 d3p708 Exp $
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
      double precision b0,b1,scal
      double precision x,x2,x3,x4,x6,a,f,Q

cccccc  external cccccccccccccccccccccccc
      real*8 lattice_omega
      external lattice_omega

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
      call dcopy(nfft3d,0.0d0,0,VL,1)
      task_count = -1
      DO 700 k3=1,nfft3
      DO 700 k2=1,nfft2
      DO 700 k1=1,nfft1
        task_count = task_count + 1
        if (mod(task_count,np).ne.taskid) go to 700

        Q=DSQRT(G(k1,k2,k3,1)**2
     >         +G(k1,k2,k3,2)**2
     >         +G(k1,k2,k3,3)**2)

        if ( (k1*k2*k3).eq.1) then
           vl(k1,k2,k3)=0.0d0
           cycle
        end if
        x = Q*rloc
        x2 = x*x
        x4 = x2*x2
        x6 = x4*x2
        a = dexp(-0.5d0*x2)

        f = b0/(Q*Q)
     >    + b1
     >        *( C1 
     >         + C2*(3.0d0-x2)
     >         + C3*(15.0d0 - 10.0d0*x2 + x2*x2)
     >         + C4*(105.0d0 - x2*(105.0d0-x2*(21.0d0 - x2))))

        vl(k1,k2,k3)= f*a
    
  700 continue
   
      call C3dB_Vector_SumAll(nfft3d,VL)
 
      return
      end



