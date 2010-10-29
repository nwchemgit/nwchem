*
* $Id$
*


      subroutine HGH_nonlocal(n,l,
     >                        r,
     >                        nfft1,nfft2,nfft3,
     >                        G,
     >                        vnl)
      implicit none
      integer          n,l
      double precision r

      integer nfft1,nfft2,nfft3
      double precision G(nfft1/2+1,nfft2,nfft3,3)
      double precision vnl(nfft1/2+1,nfft2,nfft3)

      integer np,taskid,MASTER
      parameter (MASTER=0)

*     *** local variables ****
      integer task_count,nfft3d
      integer k1,k2,k3
      double precision pi,twopi,forpi
      double precision P0,coef
      double precision r2,r3,r4,r5,r6,r7,r8,r9
      double precision g1,g2,g3,g4,a,poly


      call Parallel_np(np)
      call Parallel_taskid(taskid)

      nfft3d = (nfft1/2+1)*nfft2*nfft3
      pi=4.0d0*datan(1.0d0)
      twopi=2.0d0*pi
      forpi=4.0d0*pi

      P0 = pi**(1.25d0)

      r2 = r*r
      r3 = r2*r
      r4 = r3*r
      r5 = r4*r
      r6 = r5*r
      r7 = r6*r
      r8 = r7*r
      r9 = r8*r

      coef = 0.0d0
      if ((l.eq.0).and.(n.eq.1)) coef=4.0d0*dsqrt(2.0d0*r3)*P0
      if ((l.eq.0).and.(n.eq.2)) coef=8.0d0*dsqrt(2.0d0*r3/15.0d0)*P0
      if ((l.eq.0).and.(n.eq.3)) 
     >   coef=(16.0d0/3.0d0)*dsqrt(2.0d0*r3/105.0d0)*P0

      if ((l.eq.1).and.(n.eq.1)) coef=8.0d0*dsqrt(r5/3.0d0)*P0
      if ((l.eq.1).and.(n.eq.2)) coef=16.0d0*dsqrt(r5/105.0d0)*P0
      if ((l.eq.1).and.(n.eq.3)) 
     >   coef=(32.0d0/3.0d0)*dsqrt(r5/1155.0d0)*P0

      if ((l.eq.2).and.(n.eq.1)) coef=8.0d0*dsqrt(2*r7/15.0d0)*P0
      if ((l.eq.2).and.(n.eq.2)) 
     >    coef=(16.0d0/3.0d0)*dsqrt(2*r7/105.0d0)*P0

      if ((l.eq.3).and.(n.eq.1)) coef=16.0d0*dsqrt(r9/105.0d0)*P0


*======================  Fourier transformation  ======================
      poly = 0.0d0
      call dcopy(nfft3d,0.0d0,0,vnl,1)
      task_count = -1
      DO 700 k3=1,nfft3
      DO 700 k2=1,nfft2
      DO 700 k1=1,(nfft1/2+1)
        task_count = task_count + 1
        if (mod(task_count,np).ne.taskid) go to 700

        g2=(G(k1,k2,k3,1)**2
     >     +G(k1,k2,k3,2)**2
     >     +G(k1,k2,k3,3)**2)
        g1= dsqrt(g2)
        g3 = g2*g1
        g4 = g2*g2


        a = dexp(-0.5d0*g2*r2)

*       **** f-wave  ****
        if ((l.eq.3).and.(n.eq.1)) poly = g3
    

*       **** d-wave ****
        if ((l.eq.2).and.(n.eq.1)) poly = g2
        if ((l.eq.2).and.(n.eq.2)) poly = g2*(7.0d0-g2*r2)
        

*       **** p-wave ****
        if ((l.eq.1).and.(n.eq.1)) poly = g1
        if ((l.eq.1).and.(n.eq.2)) poly = g1*(5.0d0-g2*r2)
        if ((l.eq.1).and.(n.eq.3)) poly = g1*(35.0d0-14.0d0*g2*r2+g4*r4)


*       **** s-wave ****
        if ((l.eq.0).and.(n.eq.1)) poly = 1.0d0
        if ((l.eq.0).and.(n.eq.2)) poly = (3.0d0-g2*r2)
        if ((l.eq.0).and.(n.eq.3)) poly = (15.0d0-10.0d0*g2*r2+g4*r4)
        

        vnl(k1,k2,k3)=coef*a*poly

  700 CONTINUE
      call Parallel_Vector_SumAll(nfft3d,vnl)

      return
      end



