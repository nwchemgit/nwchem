      subroutine hnd_diag(a,vec,eig,ndim,order,righthand)
c
c $Id: hnd_diag.F,v 1.1 2004-07-21 17:28:27 bert Exp $
c
      implicit none
c
      double precision a(6), h(3,3), vec(3,3), eig(3), big(6)
      double precision zero, one, conv, xx, test, t
      integer i, j, ij, ndim, jj
      logical order, righthand
      data zero   /0.0d+00/
      data one    /1.0d+00/
      data conv   /1.0d-10/
c
c     ----- diagonalize the matrix
c
      do 20 i=1,3
         do 10 j=1,3
   10       vec(j,i)=zero
   20    vec(i,i)=one
      ij=0
      do 25 i=1,3
      do 25 j=1,i
         ij=ij+1
         h(i,j)=a(ij)
   25    h(j,i)=a(ij)
      call hnd_diajac(h,vec,eig,ndim,ndim,big)
c
      if (order) then
c
c     ----- put the principal components of the tensor and -----
c           the corresponding eigenvectors in increasing
c           order according to the absolute value of the
c           eigenvalue ( the principal component ) .
c
         do 50 i=1,3
            jj=i
            do 30 j=i,3
               if( abs(eig(j)).lt. abs(eig(jj))) jj=j
   30          continue
            if(jj.eq.i) go to 50
c
            xx=eig(jj)
            eig(jj)=eig(i)
            eig(i)=xx
            do 40 j=1,3
               xx=vec(j,jj)
               vec(j,jj)=vec(j,i)
               vec(j,i)=xx
   40       continue
c
   50    continue
      endif
c
      if (righthand) then
c
c     ----- check for right handedness, correct if not -----
c
         test = vec(1,3)*( vec(2,1)*vec(3,2) - vec(3,1)*vec(2,2) )
     1        + vec(2,3)*( vec(3,1)*vec(1,2) - vec(1,1)*vec(3,2) )
     2        + vec(3,3)*( vec(1,1)*vec(2,2) - vec(2,1)*vec(1,2) )
         if(test.gt.zero) return
         if( abs(eig(1)-eig(2)).gt.conv) go to 160
         t = eig(1)
         eig(1) = eig(2)
         eig(2) = t
         do 150 i=1,3
            t = vec(i,1)
            vec(i,1) = vec(i,2)
            vec(i,2) = t
  150    continue
         return
  160    if( abs(eig(2)-eig(3)).gt.conv) go to 180
         t = eig(2)
         eig(2) = eig(3)
         eig(3) = t
         do 170 i=1,3
            t = vec(i,2)
            vec(i,2) = vec(i,3)
            vec(i,3) = t
  170    continue
         return
  180    do 190 i=1,3
            vec(i,3) = - vec(i,3)
  190    continue
      endif
c
      return
      end
c
      subroutine hnd_jacobi(a,v,d,n,np,nrot,b,z)
      implicit none
      double precision a(np,np),v(np,np),d(np),b(np),z(np)
      double precision zero, one, pt2, pt5, hundrd
      double precision sm, thres, h, g, theta, t, s, c, tau
      integer n, np, nrot, ip, iq, i, j
      data zero   /0.0d+00/
      data one    /1.0d+00/
      data pt2    /0.2d+00/
      data pt5    /0.5d+00/
      data hundrd /1.0d+02/
c
      do 12 ip=1,n
         do 11 iq=1,n
            v(iq,ip)=zero
   11    continue
         v(ip,ip)=one
   12 continue
      do 13 ip=1,n
         b(ip)=a(ip,ip)
         d(ip)=b(ip)
         z(ip)=zero
   13 continue
      nrot=0
      if(n.eq.1) return
      do 24 i=1,50
         sm=zero
         do 15 ip=1,n-1
            do 14 iq=ip+1,n
               sm=sm+abs(a(ip,iq))
   14       continue
   15    continue
         if(sm.eq.zero) return
         if(i.lt.4) then
            thres=pt2*sm/n**2
         else
            thres=zero
         endif
         do 22 ip=1,n-1
            do 21 iq=ip+1,n
               g=hundrd*abs(a(ip,iq))
               if((i.gt.4).and.(abs(d(ip))+g.eq.abs(d(ip))).and.
     1                         (abs(d(iq))+g.eq.abs(d(iq)))) then
                  a(ip,iq)=zero
               else if(abs(a(ip,iq)).gt.thres) then
                  h=d(iq)-d(ip)
                  if(abs(h)+g.eq.abs(h)) then
                     t=a(ip,iq)/h
                  else
                     theta=pt5*h/a(ip,iq)
                     t=one/(abs(theta)+sqrt(one+theta**2))
                     if(theta.lt.zero) t=-t
                  endif
                  c=one/sqrt(one+t**2)
                  s=t*c
                  tau=s/(one+c)
                  h=t*a(ip,iq)
                  z(ip)=z(ip)-h
                  z(iq)=z(iq)+h
                  d(ip)=d(ip)-h
                  d(iq)=d(iq)+h
                  a(ip,iq)=zero
                  do 16 j=1,ip-1
                     g=a(j,ip)
                     h=a(j,iq)
                     a(j,ip)=g-s*(h+g*tau)
                     a(j,iq)=h+s*(g-h*tau)
   16             continue
                  do 17 j=ip+1,iq-1
                     g=a(ip,j)
                     h=a(j,iq)
                     a(ip,j)=g-s*(h+g*tau)
                     a(j,iq)=h+s*(g-h*tau)
   17             continue
                  do 18 j=iq+1,n
                     g=a(ip,j)
                     h=a(iq,j)
                     a(ip,j)=g-s*(h+g*tau)
                     a(iq,j)=h+s*(g-h*tau)
   18             continue
                  do 19 j=1,n
                     g=v(j,ip)
                     h=v(j,iq)
                     v(j,ip)=g-s*(h+g*tau)
                     v(j,iq)=h+s*(g-h*tau)
   19             continue
                  nrot=nrot+1
               endif
   21       continue
   22    continue
         do 23 ip=1,n
            b(ip)=b(ip)+z(ip)
            d(ip)=b(ip)
            z(ip)=zero
   23    continue
   24 continue
      write(6,*) '50 iterations should never happen'
      return
      end
c
      subroutine hnd_jacsrt(v,d,n,np)
      implicit none
      integer n, np
      double precision v(np,np),d(np)
      integer k, i, j
      double precision p
      if(n.eq.1) return
      do 13 i=1,n-1
         k=i
         p=d(i)
         do 11 j=i+1,n
            if(d(j).ge.p) then
               k=j
               p=d(j)
            endif
   11    continue
         if(k.ne.i) then
            d(k)=d(i)
            d(i)=p
            do 12 j=1,n
               p=v(j,i)
               v(j,i)=v(j,k)
               v(j,k)=p
   12       continue
         endif
   13 continue
      return
      end
