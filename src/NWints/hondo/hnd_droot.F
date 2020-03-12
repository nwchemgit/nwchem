      subroutine hnd_droot
c $Id$
      implicit none
#include "hnd_rys.fh"
      integer ierr
c
c Wrapper routine. 
c
      goto (1,1,1,2,3) nroots
      call hnd_droot1(ierr)
      if(ierr.ne.0) call hnd_droot_new
      return
    1 call hnd_rt123
      return
    2 call hnd_root4
      return
    3 call hnd_root5
      return
 678  continue
      end
c
      subroutine hnd_droot1(ierr)
      implicit none
      integer ierr
#ifdef DBLEPREC_RYS
      double precision ff, r, ww, c, s, a, rt, pt5, zero, one, four
      double precision x, wsum, dum, root, poly
#else
      real*16 ff, r, ww, c, s, a, rt, pt5, zero, one, four
      real*16 x, wsum, dum, root, poly
#endif
#include "hnd_rys.fh"
      integer n, n1, nn, i, j, k, k1, jmax, j1, m, info
      common/hnd_rysff/ff(2*maxrys+1)
      common/hnd_rysrw/r(maxrys,maxrys),ww(maxrys,maxrys)
      dimension c(maxrys+1,maxrys+1),s(maxrys+1,maxrys+1)
      dimension a(maxrys+1),rt(maxrys+1)
#ifdef DBLEPREC_RYS
      data pt5,zero,one,four /0.5d+00,0.0d+00,1.0d+00,4.0d+00/
#else
      data pt5,zero,one,four /0.5q+00,0.0q+00,1.0q+00,4.0q+00/
#endif
c
c     this version uses christoffel formula for weights.
c     ith root of the jth rys polynomial is returned in r(i,j) with
c     the corresponding weight factor in ww(i,j).   j=1,2,...,n
c
      ierr=0
      n=nroots
      x=yy
      if(n.lt.2) n=2
      n1=n+1
      nn=n+n
      call hnd_rysfun(x,nn)
      do 10 i=1,n1
      do 10 j=1,n1
   10 s(i,j)=ff(i+j-1)
      call hnd_ryssmt(c,s,n1)
      do 20 i=1,n
      do 20 j=1,i
      ww(i,j)= zero
   20 r(i,j)= zero
      wsum=ff(1)
      ww(1,1)=wsum
      r(1,1)=ff(2)/wsum
      dum= sqrt(c(2,3)**2-four *c(1,3)*c(3,3))
      r(1,2)=   pt5 *(-c(2,3)-dum)/c(3,3)
      r(2,2)=   pt5 *(-c(2,3)+dum)/c(3,3)
      if(n.eq.2) go to 70
      do 25 i=3,n1
   25 rt(i)=  one
      rt(1)=r(1,2)
      rt(2)=r(2,2)
      do 60 k=3,n
      k1=k+1
      do 30 i=1,k1
   30 a(i)=c(i,k1)
      call hnd_rysnod(a,rt,k,ierr)
      if(ierr.ne.0) return
      do 50 i=1,k
   50 r(i,k)=rt(i)
   60 continue
   70 do 150 k=2,n
      jmax=k-1
      do 150 i=1,k
      root=r(i,k)
      dum=  one  /ff(1)
      do 110 j=1,jmax
      j1=j+1
      poly=c(j1,j1)
      do 100 m=1,j
  100 poly=poly*root+c(j1-m,j1)
  110 dum=dum+poly*poly
  150 ww(i,k)=  one  /dum
      do 160 k=1,nroots
      dum=r(k,nroots)
      u9(k)=dum/(  one  -dum)
  160 w9(k)=ww(k,nroots)
      return
      end
c
      subroutine hnd_rysfun(x,n)
      implicit none
#include "hnd_rys.fh"
#ifdef DBLEPREC_RYS
      double precision ff, pt5, one, two, x
      double precision e, xx, facmin, term 
      double precision a, sum, fac, t, s, tmax
#else
      real*16 ff, pt5, one, two, x
      real*16 e, xx, facmin, term 
      real*16 a, sum, fac, t, s, tmax
#endif
      integer n, m
      common/hnd_rysff/ff(2*maxrys+1)
#ifdef DBLEPREC_RYS
      data pt5,one,two /0.5d+00,1.0d+00,2.0d+00/
      e=0.5409855304296342219319112d-78
      xx=x+x
      facmin=xx
      if(facmin.lt.2*180.2160d+00) e= exp(-x)
      if(facmin.gt.   90.0000d+00) go to 100
#else
      data pt5,one,two /0.5q+00,1.0q+00,2.0q+00/
      e=0.5409855304296342219319112q-78
      xx=x+x
      facmin=xx
      if(facmin.lt.2*180.2160q+00) e= exp(-x)
      if(facmin.gt.   90.0000q+00) go to 100
#endif
      term=one
      sum =one
      fac=n
      fac=fac+pt5
   10 fac=fac+one
      term=term*x/fac
      sum=sum+term
      if(fac.le.facmin) go to 10
      t=term
      s=sum
      if(t.gt.s*hnd_acc) go to 10
      fac=n+n+1
      ff(n+1)=sum*e/fac
      m=n-1
      fac=m+m+1
   20 if(m.lt.0) return
      ff(m+1)=(e+xx*ff(m+2))/fac
      m=m-1
      fac=fac-two
      go to 20
c
c     use asymptotic expansion for large arguments.
c
#ifdef DBLEPREC_RYS
  100 a= sqrt(.7853981633974483096156608d+00/x)
#else
  100 a= sqrt(.7853981633974483096156608q+00/x)
#endif
      tmax=a*hnd_acc/e
      term=one/xx
      sum=term
      fac=one
  110 fac=fac-two
      term=fac*term/xx
      sum=term+sum
      t=term
      if( abs(t).gt.tmax) go to 110
      ff(1)=a-e*sum
      fac=-one
      m=0
  120 if(m.eq.n) return
      m=m+1
      fac=fac+two
      ff(m+1)=(fac*ff(m)-e)/xx
      go to 120
      end
c
      subroutine hnd_ryssmt(c,s,n)
      implicit none
#include "hnd_rys.fh"
#ifdef DBLEPREC_RYS
      double precision c, s, v, y, zero, one
      double precision fac, dot
#else
      real*16 c, s, v, y, zero, one
      real*16 fac, dot
#endif
      integer n, i, j, k, kmax,ii,jj
      dimension c(maxrys+1,maxrys+1),s(maxrys+1,maxrys+1)
      dimension v(maxrys+1),y(maxrys+1)
#ifdef DBLEPREC_RYS
      data zero,one /0.0d+00,1.0d+00/
#else
      data zero,one /0.0q+00,1.0q+00/
#endif
c
c     routine returns an n by n triangular matrix c such that
c     c(transpose)sc=i,  where i is an n by n identity matrix.
c
      do 10 i=1,n
      do 10 j=1,i
   10 c(i,j)= zero
      do 100 j=1,n
      kmax=j-1
      fac=s(j,j)
      if(kmax.eq.0) go to 60
      do 20 k=1,kmax
      v(k)= zero
   20 y(k)=s(k,j)
      do 50 k=1,kmax
      dot= zero
      do 30 i=1,k
   30 dot=c(i,k)*y(i)+dot
      do 40 i=1,k
   40 v(i)=v(i)-dot*c(i,k)
   50 fac=fac-dot*dot
   60 fac=one/ sqrt(fac)
      c(j,j)=fac
      if(kmax.eq.0) go to 100
      do 70 k=1,kmax
   70 c(k,j)=fac*v(k)
  100 continue
      return
      end
c
      subroutine hnd_rysnod(a,rt,k,ierr)
      implicit none
#include "hnd_rys.fh"
#ifdef DBLEPREC_RYS
      double precision a, rt, zero, tol, r1, p1, r2, p2
      double precision prod, r3, p3, r4, p4, r5, p5, r6, p6
      double precision r, dr, delta, ff5, ff6
#else
      real*16 a, rt, zero, tol, r1, p1, r2, p2
      real*16 prod, r3, p3, r4, p4, r5, p5, r6, p6
      real*16 r, dr, delta, ff5, ff6
#endif
#include "stdio.fh"
#include "errquit.fh"
      integer ierr
      integer k, k1, m, i
      dimension a(maxrys+1),rt(maxrys+1)
#ifdef DBLEPREC_RYS
      data zero /0.0d+00/
#else
      data zero /0.0q+00/
#endif
c
c     routine returns rt(i) the ith root of a polynomial of order
c     k whose mth coefficient is stored in a(m+1).  it is assumed that
c     the initial values in rt bracket the final values.
c
      ierr=0
#ifdef DBLEPREC_RYS
      tol=1.0d-09
#else
      tol=1.0q-11
#endif
      k1=k+1
      r2= zero
      p2=a(1)
      do 100 m=1,k
      r1=r2
      p1=p2
      r2=rt(m)
      p2=a(k1)
      do 10 i=1,k
   10 p2=p2*r2+a(k1-i)
      prod=p1*p2
      if(prod.lt. zero) go to 20
      write(luout,15) m,k
   15 format(//,' root number ',i3,' was not found for polynomial',
     1 ' of order ',i3,//)
#if 1
      ierr=-1
      return
#else
      call errquit('hnd_rysnod: root not found', 555,
     &       INT_ERR)
#endif
   20 r5=r1
      p5=p1
      r6=r2
      p6=p2
   30 r3=r5
      p3=p5
      r4=r6
      p4=p6
      r =(r3*p4-r4*p3)/(p4-p3)
      dr=r4-r3
      delta=dr
      if( abs(delta).lt.tol) go to 90
#ifdef DBLEPREC_RYS
      dr=0.0625d+00*dr
#else
      dr=0.0625q+00*dr
#endif
      r5=r-dr
      if(r5.lt.r3) r5=r3
      r6=r+dr
      if(r6.gt.r4) r6=r4
      p5=a(k1)
      p6=p5
      do i=1,k
         p5=p5*r5+a(k1-i)
         p6=p6*r6+a(k1-i)
      enddo
   45 prod=p5*p6
      if(prod.lt. zero) go to 30
      prod=p3*p5
      if(prod.gt. zero) go to 60
#ifdef DBLEPREC_RYS
      r5=0.25d+00*r3+0.75d+00*r5
#else
      r5=0.25q+00*r3+0.75q+00*r5
#endif
      p5=a(k1)
      do 50 i=1,k
   50 p5=p5*r5+a(k1-i)
      go to 45
#ifdef DBLEPREC_RYS
   60 r6=0.25d+00*r4+0.75d+00*r6
#else
   60 r6=0.25q+00*r4+0.75q+00*r6
#endif
      p6=a(k1)
      do 70 i=1,k
   70 p6=p6*r6+a(k1-i)
      go to 45
   90 rt(m)=r
  100 continue
      return
      end
c
      subroutine hnd_droot_new_set
      implicit none
#include "hnd_rys.fh"
      double precision tmp_rts(maxleg),tmp_wts(maxleg)
      double precision a(maxleg),b(maxleg),buff(maxleg)
      double precision tmp
      integer i,j,n2,ierr
c
c     Z values for which we switch to cheapest algorithm
c
c     Empirical switch from large to small based on Ishida
c
      data large_z / 37.0d0, 44.0d0, 51.0d0, 58.0d0, 64.0d0, 69.0d0,
     &               75.0d0, 80.0d0, 86.0d0, 91.0d0, 96.0d0,101.0d0,
     &              105.0d0,106.0d0,111.0d0,116.0d0,121.0d0,126.0d0,
     &              131.0d0,136.0d0 /
c
c     Unique auxiliary Legendre (Shifted) quadrature grid dimensions to 
c     approximate roots and weights for Gauss-Rys quadrature
c
      data leg_aux / 20,25,30,35,40,45,50,55,60,65,70,75 /
c
c     Different Gauss-Rys roots can have the same grid dimension
c     I.e. make appropriate mapping from 20 to 12 unique grids
c
      data rys_leg / 1,2,3,3,4,5,5,5,6,7,7,8,8,9,10,10,11,11,12,12 /
c
c     Roots and weights of 2n-point Gauss-Hermite quadrature independent of Z
c     Always an even function, so only store half of the data
c
      do i=1,maxrys
         n2=2*i
         do j=1,n2
            a(j)=0.0d0
            b(j)=0.5d0*dble(j-1)
         enddo
         b(1)=dsqrt(4.0d0*datan(1.0d0))
         call gaussh(n2,a,b,hnd_acc,tmp_rts,tmp_wts,ierr,buff)
         do j=1,i  
            rts_large(j,i)=tmp_rts(j)*tmp_rts(j)
            wts_large(j,i)=tmp_wts(j)
         enddo
      enddo
c
c     Generate roots and weights for shifted Legendre quadrature
c 
c     First set up a and b. As these don't change, do just the largest grid
c
      n2=leg_aux(maxuniqleg)
      do i=2,n2
         a(i)=0.5d0
         tmp=dble(i-1)
         b(i)=0.25d0/(4.0d0-(1.0d0/(tmp*tmp)))
      enddo
      a(1)=0.5d0
      b(1)=1.0d0
c
c     Generate roots and weights for Legendre quadrature
c
      do i=1,maxuniqleg
         call gaussh(leg_aux(i),a,b,hnd_acc,tmp_rts,tmp_wts,ierr,buff)
         do j=1,leg_aux(i)
            rts_leg(j,i)=tmp_rts(j)*tmp_rts(j)
            wts_leg(j,i)=tmp_wts(j)
         enddo
      enddo
c
      return
      end
c
      subroutine hnd_droot_new
      implicit none
#include "hnd_rys.fh"
      double precision tmp,p0,p1,p2,tmp_rts,tmp_wts,a,b,buff
      dimension tmp(maxleg),p0(maxleg),p1(maxleg),p2(maxleg)
      dimension tmp_rts(maxrys),tmp_wts(maxrys),buff(maxrys)
      dimension a(maxrys),b(maxrys)
      double precision temp,rfac,wfac
      integer ierr,i,naux,map,offset
c
      call hnd_droot_new_set ! initialize rts_large,wts_large
c
      if (yy.ge.large_z(nroots)) then
c 
c     Large Z values, use roots and weights from 2n-point Gauss-Hermite quadrature
c
         rfac=1.0d0/yy
         wfac=sqrt(rfac)
         do i=1,nroots
            temp=rfac*rts_large(i,nroots)
            u9(i)=temp/(1.0d0-temp)
            w9(i)=wfac*wts_large(i,nroots)
         enddo
      else  
c 
c     Use Golub-Welsh method
c
         offset=rys_leg(nroots)
         do i=1,leg_aux(offset)
            tmp(i)=wts_leg(i,offset)*exp(-yy*rts_leg(i,offset))
         enddo
c        call sti(nroots,leg_aux(offset),rts_leg(1,nroots),tmp,a,b,ierr,
         call sti(nroots,leg_aux(offset),rts_leg(1,offset),tmp,a,b,ierr,
     &            p0,p1,p2)
         call gaussh(nroots,a,b,hnd_acc,tmp_rts,tmp_wts,ierr,buff)
         do i=1,nroots
            temp=tmp_rts(i)
            u9(i)=temp/(1.0d0-temp)
            w9(i)=tmp_wts(i)
         enddo
      endif
c
      return
      end
c
      subroutine gaussh(n,alpha,beta,eps,zero,weigh,ierr,e)
c
c Given  n  and a measure  dlambda, this routine generates the n-point
c Gaussian quadrature formula
c
c     integral over supp(dlambda) of f(x)dlambda(x)
c
c        = sum from k=1 to k=n of w(k)f(x(k)) + R(n;f).
c
c The nodes are returned as  zero(k)=x(k) and the weights as
c weight(k)=w(k), k=1,2,...,n. The user has to supply the recursion
c coefficients  alpha(k), beta(k), k=0,1,2,...,n-1, for the measure
c dlambda. The routine computes the nodes as eigenvalues, and the
c weights in term of the first component of the respective normalized
c eigenvectors of the n-th order Jacobi matrix associated with  dlambda.
c It uses a translation and adaptation of the algol procedure  imtql2,
c Numer. Math. 12, 1968, 377-383, by Martin and Wilkinson, as modified
c by Dubrulle, Numer. Math. 15, 1970, 450. See also Handbook for
c Autom. Comput., vol. 2 - Linear Algebra, pp.241-248, and the eispack
c routine  imtql2.
c
c        Input:  n - - the number of points in the Gaussian quadrature
c                      formula; type integer
c                alpha,beta - - arrays of dimension  n  to be filled
c                      with the values of  alpha(k-1), beta(k-1), k=1,2,
c                      ...,n
c                eps - the relative accuracy desired in the nodes
c                      and weights
c
c        Output: zero- array of dimension  n  containing the Gaussian
c                      nodes (in increasing order)  zero(k)=x(k), k=1,2,
c                      ...,n
c                weight - array of dimension  n  containing the
c                      Gaussian weights  weight(k)=w(k), k=1,2,...,n
c                ierr- an error flag equal to  0  on normal return,
c                      equal to  i  if the QR algorithm does not
c                      converge within 30 iterations on evaluating the
c                      i-th eigenvalue, equal to  -1  if  n  is not in
c                      range, and equal to  -2  if one of the beta's is
c                      negative.
c
c The array  e  is needed for working space.
c
      implicit none
      integer n,ierr,k,l,m,i,j,ii,mml
      double precision alpha,beta,eps,zero,weigh,e,dp,dg,dr,
     *ds,dc,df,db
      dimension alpha(n),beta(n),zero(n),weigh(n),e(n)
c
      if(n.lt.1) then
        ierr=-1
        return
      end if
      ierr=0
      zero(1)=alpha(1)
      if(beta(1).lt.0.0d0) then
        ierr=-2
        return
      end if
      weigh(1)=beta(1)
      if (n.eq.1) return
      weigh(1)=1.0d0
      e(n)=0.0d0
      do 100 k=2,n
        zero(k)=alpha(k)
        if(beta(k).lt.0.0d0) then
          ierr=-2
          return
        end if
        e(k-1)=sqrt(beta(k))
        weigh(k)=0.0d0
  100 continue
      do 240 l=1,n
        j=0
  105   do 110 m=l,n
          if(m.eq.n) goto 120
          if(abs(e(m)).le.eps*(abs(zero(m))+abs(zero(m+1))))
     *      goto 120
  110   continue
  120   dp=zero(l)
        if(m.eq.l) goto 240
        if(j.eq.30) goto 400
        j=j+1
        dg=(zero(l+1)-dp)/(2.0d0*e(l))
        dr=sqrt(dg*dg+1.0d0) 
        dg=zero(m)-dp+e(l)/(dg+sign(dr,dg))
        ds=1.0d0
        dc=1.0d0
        dp=0.0d0 
        mml=m-l
        do 200 ii=1,mml
          i=m-ii
          df=ds*e(i)
          db=dc*e(i)
          if(abs(df).lt.abs(dg)) goto 150
          dc=dg/df
          dr=sqrt(dc*dc+1.0d0)
          e(i+1)=df*dr
          ds=1.0d0/dr
          dc=dc*ds
          goto 160
  150     ds=df/dg
          dr=sqrt(ds*ds+1.0d0)
          e(i+1)=dg*dr
          dc=1.0d0/dr
          ds=ds*dc
  160     dg=zero(i+1)-dp
          dr=(zero(i)-dg)*ds+2.0d0*dc*db
          dp=ds*dr
          zero(i+1)=dg+dp
          dg=dc*dr-db
          df=weigh(i+1)
          weigh(i+1)=ds*weigh(i)+dc*df
          weigh(i)=dc*weigh(i)-ds*df
  200   continue
        zero(l)=zero(l)-dp
        e(l)=dg
        e(m)=0.0d0
        goto 105
  240 continue
      do 300 ii=2,n
        i=ii-1
        k=i
        dp=zero(i)
        do 260 j=ii,n
          if(zero(j).ge.dp) goto 260
          k=j
          dp=zero(j)
  260   continue
        if(k.eq.i) goto 300
        zero(k)=zero(i)
        zero(i)=dp
        dp=weigh(i)
        weigh(i)=weigh(k)
        weigh(k)=dp
  300 continue
      do 310 k=1,n
        weigh(k)=beta(1)*weigh(k)*weigh(k)
  310 continue
      return
  400 ierr=l
      return
      end
c
      subroutine sti(n,ncap,x,w,alpha,beta,ierr,p0,p1,p2)
c
c This routine applies ``Stieltjes's procedure'' (cf. Section 2.1 of
c W. Gautschi,``On generating orthogonal polynomials'', SIAM J. Sci.
c Statist. Comput. 3, 1982, 289-317) to generate the recursion
c coefficients  alpha(k), beta(k) , k=0,1,...,n-1, for the discrete
c (monic) orthogonal polynomials associated with the inner product
c
c     (f,g)=sum over k from 1 to ncap of w(k)*f(x(k))*g(x(k)).
c
c The integer  n  must be between  1  and  ncap, inclusive; otherwise,
c there is an error exit with  ierr=1. The results are stored in the
c arrays  alpha, beta; the arrays  p0, p1, p2  are working arrays.
c
c If there is a threat of underflow or overflow in the calculation
c of the coefficients  alpha(k)  and  beta(k), the routine exits with
c the error flag  ierr  set equal to  -k  (in the case of underflow)
c or  +k  (in the case of overflow), where  k  is the recursion index
c for which the problem occurs. The former [latter] can often be avoided
c by multiplying all weights  w(k)  by a sufficiently large [small]
c scaling factor prior to entering the routine, and, upon exit, divide
c the coefficient  beta(0)  by the same factor.
c
c This routine should be used with caution if  n  is relatively close
c to  ncap, since there is a distinct possibility of numerical
c instability developing. (See W. Gautschi,``Is the recurrence relation
c for orthogonal polynomials always stable?'', BIT, 1993, to appear.)
c In that case, the routine  lancz  should be used.
c
c The routine uses the function subroutine  r1mach.
c
      implicit none
      integer ierr,n,ncap,nm1,m,k
      double precision x,w,alpha,beta,p0,p1,p2,dtiny,dhuge,dsum0,dsum1,
     *                 dsum2,dt
      dimension x(ncap),w(ncap),alpha(n),beta(n),p0(ncap),p1(ncap),
     *          p2(ncap)
c
      dtiny=1.0d-40
      dhuge=1.0d40
      ierr=0
      if(n.le.0 .or. n.gt.ncap) then
        ierr=1
        return
      end if
      nm1=n-1
      dsum0=0.0d0
      dsum1=0.0d0
      do 10 m=1,ncap
        dsum0=dsum0+w(m)
        dsum1=dsum1+w(m)*x(m)
   10 continue
      alpha(1)=dsum1/dsum0
      beta(1)=dsum0
      if(n.eq.1) return
      do 20 m=1,ncap
        p1(m)=0.0d0
        p2(m)=1.0d0
   20 continue
      do 40 k=1,nm1
        dsum1=0.0d0
        dsum2=0.0d0
        do 30 m=1,ncap
          if(w(m).eq.0.0d0) goto 30
          p0(m)=p1(m)
          p1(m)=p2(m)
          p2(m)=(x(m)-alpha(k))*p1(m)-beta(k)*p0(m)
          if(abs(p2(m)).gt.dhuge .or. abs(dsum2).gt.dhuge) then
            ierr=k
            return
          end if
          dt=w(m)*p2(m)*p2(m)
          dsum1=dsum1+dt
          dsum2=dsum2+dt*x(m)
   30   continue
        if(abs(dsum1).lt.dtiny) then
          ierr=-k
          return
        end if
        alpha(k+1)=dsum2/dsum1
        beta(k+1)=dsum1/dsum0
        dsum0=dsum1
   40 continue
      return
      end
