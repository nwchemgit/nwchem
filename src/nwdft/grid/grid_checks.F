      subroutine grid_checks(exex,ll)
      implicit none
#include "cdft.fh"
      integer ll(*) ! [out] for lindh grid
      double precision exex(*) ! [out] for lindh grid
c
      double precision pi,x
      double precision pts(MXrad)
      double precision wghts(MXrad)
      integer itype,n0,L
      double precision s_err,alpha,gg,nn
      double precision grid_sint,range_mk,grid_norm
      external grid_sint,grid_norm
c      double precision norm2
c      norm2(x)=sqrt(dacos(-1d0)/(x*x*x))*0.25d0
c
      pi=dacos(-1d0)
c
      do itype=1,ntypes
         n0=nrad(itype)
         range_mk=rad_cutoff(1,itype)
             call xc_murakn(
     .            pts, 
     &            wghts,n0,range_mk,
     .            .true.)
             alpha=exex(itype)
             l=ll(itype)
             gg=grid_sint(n0,pts,wghts,alpha,l)
c             s_err = gg - norm2(alpha)
             nn= grid_norm(l,alpha)
             s_err = gg - nn
             write(0,*) ' **********'
             write(0,1) alpha,l,'  err  ',s_err/nn
             write(0,*) ' **********'
 1           format(' exp ',d10.2,' l ',i2,A,d10.2,A,d10.2)
      enddo
      return
      end
      
      double precision function grid_sint(n0,pts,wghts,alpha,l)
      implicit none
      integer n0,l
      double precision pts(n0)
      double precision wghts(n0)
      double precision alpha
c
      integer i
c
      grid_sint=0d0
c      write(0,*) ' n0 ',n0
      do i=1,n0
         if(l.eq.0) then
         grid_sint=grid_sint+wghts(i)*exp(-alpha*pts(i)*pts(i))
         else
         grid_sint=grid_sint+wghts(i)*exp(-alpha*pts(i)*pts(i))*
     .        pts(i)**l
         endif
c         write(0,11) wghts(i),pts(i),grid_sint
c 11       format(' int ',3F20.16)
      enddo
      return
      end
      double precision function grid_norm(l,alpha)
      implicit none
      integer l
      double precision alpha
c
      logical even
      integer nff,n,n_in
      double precision sqrtpi
      integer grid_fact,grid_ffact
      external grid_fact,grid_ffact
c
      sqrtpi=sqrt(dacos(-1d0))
      even=mod(l,2).eq.0
      if(even) then
         n=(l+2)/2
c     
c     compute 2n-1 factorial ... 
c     
         n_in=2*n-1
         nff=grid_ffact(n_in)
         write(0,*) ' 2n-1',2*n-1
         write(0,*) l,' even fact ',nff
         grid_norm=sqrtpi*alpha**(-n-0.5d0)*nff*2d0**(-n-1)
      else
c     
c     compute n factorial ... 
c     
         n=(l+1)/2
         nff=grid_fact(n)
         write(0,*) l,' odd fact ',nff
         grid_norm=alpha**(-n-1)*nff*0.5d0
      endif
      return
      end
      integer function grid_fact(n)
      integer n
      integer i,pp
      pp=1
      if(n.gt.1) then
         pp=1
         do i=2,n
            pp=i*pp
         enddo
      endif
      grid_fact=pp
      write(0,*) ' fl ',grid_fact
      return
      end
      integer function grid_ffact(n)
      integer n
      integer i,pp
      pp=1
      write(0,*) ' ff: n ',n
      if(n.gt.1) then
         pp=1
         do i=3,n,2
            pp=i*pp
         enddo
      endif
      grid_ffact=pp
      write(0,*) ' ffl ',grid_ffact
      return
      end
      
