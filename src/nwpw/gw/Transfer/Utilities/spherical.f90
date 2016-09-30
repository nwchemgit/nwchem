      double precision function tripleYlmint(l1,l2,l3,m1,m2)
! evaluates \int d\Omega Y_{l1}^{m1} Y_{l2}^{m2} Y_{l3}^{m3} 
! where m3=-m1-m2
! Requires functions lnfactorial and lngamma
      implicit none
      integer l1,l2,l3,m1,m2
      double precision wig1,wig2,xterm,pi
      parameter (pi=3.1415926535897932384626433832795028841972d0)
      double precision wigner3j
      external wigner3j
      tripleYlmint=0.d0
      if (l3.lt.abs(l1-l2).or.l3.gt.l1+l2) return
      if (m1.lt.-l1.or.m1.gt.l1) return
      if (m2.lt.-l2.or.m2.gt.l2) return
      wig1=wigner3j(l1,l2,l3,0,0)
      wig2=wigner3j(l1,l2,l3,m1,m2)
      xterm=sqrt(dble((2*l1+1)*(2*l2+1)*(2*l3+1))/(4.d0*pi))
      tripleYlmint=xterm*wig1*wig2
      return
      end

!***********************************************************************

      double precision function wigner3j(j1,j2,j3,m1,m2)
! evaluates | j1 j2 j3 |
!           | m1 m2 m3 |
! where m3=-m1-m2
      implicit none
      integer j1,j2,j3,m1,m2,m3
      integer it,itmax,itmin
      double precision xterm1,xterm2,sterm
      double precision triangle
      double precision lnfactorial
      external lnfactorial
      m3=-m1-m2
      wigner3j=0.d0
      if (j3.lt.abs(j1-j2).or.j3.gt.j1+j2) return
      if (m1.lt.-j1.or.m1.gt.j1) return
      if (m2.lt.-j2.or.m2.gt.j2) return
      if (m3.lt.-j3.or.m3.gt.j3) return
      triangle=(lnfactorial(j1+j2-j3)+lnfactorial(j1-j2+j3)  &
&             +lnfactorial(-j1+j2+j3)-lnfactorial(j1+j2+j3+1))/2.d0
      xterm1=(lnfactorial(j1+m1)+lnfactorial(j1-m1)  &
&            +lnfactorial(j2+m2)+lnfactorial(j2-m2)  &
&            +lnfactorial(j3+m3)+lnfactorial(j3-m3))/2.d0
      xterm2=0.d0
      itmin=max(0,j2-j3-m1,j1-j3+m2)
      itmax=min(j1+j2-j3,j1-m1,j2+m2)
      do it=itmin,itmax
        sterm=exp(triangle+xterm1-lnfactorial(it)  &
&                -lnfactorial(j3-j2+it+m1)-lnfactorial(j3-j1+it-m2)  &
&                -lnfactorial(j1+j2-j3-it)-lnfactorial(j1-it-m1)  &
&                -lnfactorial(j2-it+m2))
        xterm2=xterm2+((-1)**it)*sterm
      enddo
      wigner3j=((-1)**(j1-j2-m3))*xterm2
      return
      end

!***********************************************************************

      function spharm(ll,mm,rr)
! Spherical harmonic Y_ll^mm for radius vector rr
! cosine of polar angle theta = rr(3)/|rr|
! azimuthal angle phi = atan2(rr(2),rr(1))
      implicit none
      double complex :: spharm
      double precision :: rr(3)
      integer :: ll,mm
      double precision :: costheta,cosphi,sinphi,rmag
      double complex :: expphi
      integer :: im,mmm
      double precision :: fact,prefact,pi,lpolyn
      parameter (pi=3.1415926535897932384626433832795028841972d0)
      external lpolyn
      
      mmm=abs(mm)
      if (ll.lt.mmm) then
        write(6,*) 'bad input parameters in function spharm'
        stop
      endif
      fact=1.d0
      do im=ll-mmm+1,ll+mmm
        fact=fact*im
      enddo
      costheta=rr(3)/sqrt(dot_product(rr,rr))
      rmag=sqrt(dot_product(rr(1:2),rr(1:2)))
      if (rmag.le.0.d0) then
        if (mmm.eq.0) then
          expphi=(1.d0,0.d0)
        else
          expphi=(0.d0,0.d0)
        endif
      else
        cosphi=rr(1)/rmag
        sinphi=rr(2)/rmag
        expphi=cosphi+(0.d0,1.d0)*sinphi
      endif
      prefact=sqrt((2.d0*ll+1.d0)/(fact*4.d0*pi))
      spharm=prefact*lpolyn(ll,mmm,costheta)*expphi**mmm
      if (mm.lt.0) spharm=conjg(spharm)*(-1)**mmm

      return
      end

!***********************************************************************

      function lpolyn(ll,mm,xx)
! Legendre Polynomial P_ll^mm (xx)
      implicit none
      double precision lpolyn,xx
      integer ll,mm
      double precision term,fact,lpstart,lpl,lplm1,lplm2
      integer im,il

      if (ll.lt.mm.or.abs(xx).gt.1.d0) then
        write(6,*) 'bad input parameters for function lpolyn'
        stop
      endif
      term=sqrt(1.d0-xx*xx)
      lpstart=1.d0
      do im=1,mm
        lpstart=-lpstart*term*(2*im-1)
      enddo
      if (ll.eq.mm) then
        lpolyn=lpstart
      else
        lplm1=lpstart
        lplm2=0.d0
        do il=mm+1,ll
          lpl=(xx*(2*il-1)*lplm1-(il+mm-1)*lplm2)/dble(il-mm)
          lplm2=lplm1
          lplm1=lpl
        enddo
        lpolyn=lpl
      endif

      return
      end
 
!***********************************************************************

      double precision function lngamma(xx)
      implicit double precision (a-h,o-z)
      double precision xx,x, cc(6), gam, sumg, work, sq2pi
      integer i
      data cc/76.18009172947146d0,-86.50532032941677d0,  &
&             24.01409824083091d0,-1.231739572450155d0,  &
&             0.1208650973866179d-2,-0.5395239384953d-5/
      x=xx-1.d0
      sq2pi=2.5066282746310005d0
      sumg=1.000000000190015d0
      do i=1,6
        sumg=sumg+cc(i)/(x+i)
      enddo
      work=x+5.5
      lngamma=(x+0.5)*log(work)-work+log(sq2pi*sumg)
      return
      end

!***********************************************************************

      double precision function lnfactorial(ii)
      implicit none
      integer ii
      double precision factstore(0:100),lngamma
      external lngamma
      save factstore
      data factstore/101*-1.d0/
      if (ii.le.100) then
        if (factstore(ii).lt.0.d0) factstore(ii)=lngamma(dble(ii+1))
        lnfactorial=factstore(ii)
      else
        lnfactorial=lngamma(dble(ii+1))
      endif
      return
      end

!***********************************************************************

      function sphbesselj(xx,nn)
      implicit none
      double precision sphbesselj,xx
      integer nn
      integer ii,np,sigfig
      double precision tiny,dj0,dj1,djnn,djn,djnp1,djnm1,scale

      if (nn.lt.0) then
        write(6,*) 'bad argument for order n of j_n(xx) in sphbesselj'
        write(6,*) 'nn = ',nn
        stop
      endif
      sigfig=12
      tiny=1.d-10
      if (xx.lt.tiny) then
        dj0=1.d0-xx**2/6.d0
      else
        dj0=sin(xx)/xx
      endif
      if (nn.eq.0) then
        sphbesselj=dj0
        return
      endif
      if (xx.lt.tiny) then
        sphbesselj=1.d0
        do ii=1,nn
          sphbesselj=sphbesselj*xx/dble(2*ii+1)
        enddo
      elseif (xx.lt.2*nn) then
        djn=tiny
        djnp1=0.d0
        np=int(sqrt(dble(nn*sigfig**2)))+1
        do ii=nn+np,1,-1
          djnm1=(2*ii+1)*djn/xx-djnp1
          if (ii.eq.nn) djnn=djn
          djnp1=djn
          djn=djnm1
        enddo
        scale=djn/dj0
        sphbesselj=djnn/scale
      else
        if (xx.lt.tiny) then
          dj1=xx/3.d0-4.d0*xx**3/30.d0
        else
          dj1=sin(xx)/(xx**2)-cos(xx)/xx
        endif
        djnm1=dj0
        djn=dj1
        do ii=1,nn-1
          djnp1=(2*ii+1)*djn/xx-djnm1
          djnm1=djn
          djn=djnp1
        enddo
        sphbesselj=djn
      endif

      return
      end
