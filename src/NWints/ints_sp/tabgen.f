C$PRAGMA SUN OPT=2
      subroutine tabgen
c $Id$
c
c        *****  computes and tabulates f0(x) to f5(x)           *****
c        *****  in range x = -0.24 to x = 26.4                  *****
c        *****  in units of x = 0.08                            *****
c        *****  the two electron integral sp routines           *****
c        *****  the table is generated only once for each entry *****
c
      implicit none
      double precision c, ppp(350)
      common/tabint/ c(1000,6)
      double precision pt184, pt5, six, tenm7, four, two, done, pt886
      integer mm, l, i, m, j, notrms
      double precision q, qqq, a, term, ptlsum, b, temp1, temp2
      double precision approx, fimult, fiprop
      data pt184,pt5/ 0.184d0,0.50d0/
      data six,tenm7/6.0d0,1.0d-20 /
      data four,two,done/4.0d0,2.0d0,1.0d0/
      data pt886/0.8862269254527d0/
c
      q=-done
      do 30 mm=1,6
      m=mm-1
      q=q+done
      qqq = -0.24d0
      do 20 i=1,340
      qqq = qqq+0.08d0
      a=q
c        *****  change limit of approximate solution.           *****
      if(qqq-15.0d0) 1,1,10
    1 a=a+pt5
      term=done/a
      ptlsum=term
      do 2 l=2,50
      a=a+done
      term=term*qqq/a
      ptlsum=ptlsum+term
      if( dabs(term/ptlsum)-tenm7)3,2,2
    2 continue
    3 ppp(i)=pt5*ptlsum* dexp(-qqq)
      go to 20
   10 b=a+pt5
      a=a-pt5
      approx=pt886/(dsqrt(qqq)*qqq**m)
      if(m.eq.0) go to 13
      do 12 l=1,m
      b=b-done
   12 approx=approx*b
   13 fimult=pt5* dexp(-qqq)/qqq
      fiprop=fimult/approx
      term=done
      ptlsum=term
      notrms=qqq
      notrms=notrms+m
      do 14 l=2,notrms
      term=term*a/qqq
      ptlsum=ptlsum+term
      if( dabs(term*fiprop/ptlsum)-tenm7)15,15,14
   14 a=a-done
   15 ppp(i)=approx-fimult*ptlsum
   20 continue
      do 30 i=1,333
      j=i+2
      c(i,mm)=ppp(j)
      c(i+333,mm)=ppp(j+1)-ppp(j)
      temp1=-two*ppp(j)+ppp(j+1)+ppp(j-1)
      temp2=six*ppp(j)-four*ppp(j+1)-four*ppp(j-1)+ppp(j-2)+ppp(j+2)
   30 c(i+666,mm) = (temp1-pt184*temp2)/six
      return
      end
