      double precision function linesrch(n,x,f,g)
*
* $Id$
*
      implicit none
      integer n
      double precision x(n)
      double precision f(n)
      double precision g(n)
      integer ia, ib, ic
      double precision dx, tol, prec
      double precision step
      logical linebracket
      external linebracket

      write(6,*)
      write(6,*)
      do ia=max(1,(n-4)),n
        write(6,899) ia,x(ia),f(ia),g(ia)
 899    format(i5,2f24.12,e20.6)
      enddo
      if (n.eq.2) then
        dx = (g(n) - g(n-1))/(x(n)-x(n-1))
        write(6,884) dx
 884    format('Curvature:',e15.4)
        if ((dx.lt.0.d0).and.(f(n).gt.f(n-1))) then
          linesrch = -x(n)
          return
        endif
      endif
      tol = 1.d-13
      step = (g(n)*x(n-1) - g(n-1)*x(n))/(g(n) - g(n-1))
c
c  Bracketing
c
      if (n.gt.2) then
        if (linebracket(n,x,f,ia,ib,ic)) then
          write(6,900) ib,x(ib),ia,x(ia),(f(ia)-f(ib)),
     $                          ic,x(ic),(f(ic)-f(ib))
 900      format(  'Minimum:',10x,i5,f22.14
     $           /,'Left Bracket:',5x,i5,f22.14,5x,e14.6,
     $           /,'Right Bracket:',4x,i5,f22.14,5x,e14.6)
          write(6,771) step,n,n-1
 771      format(5x,'Step:',f20.10,' from ',i3,'-',i3)
c
c  If last step is not the minimum cannot use it
c  in NR extrapolation, use current brackets and 
c  sign of gradient for extrapolation.
c
          if (ib.ne.n) then
            if (g(ib).gt.0) then
              step = (g(ia)*x(ib) - g(ib)*x(ia))/(g(ia) - g(ib))
              write(6,771) step,ia,ib
            else 
              step = (g(ic)*x(ib) - g(ib)*x(ic))/(g(ic) - g(ib))
              write(6,771) step,ib,ic
            endif
          endif
c
c  Apply bracketing criterion only if brackets 
c  are greater than sqrt(relative machine-precision)*x
c  and difference in function values is greater than machine precision
c
          prec = abs(sqrt(tol)*x(ib))
          if (((step.lt.x(ia)).and.(abs(x(ib)-x(ia)).gt.prec).and.
     $         (abs(f(ia)-f(ib)).gt.(tol*f(ib)))).or.
     $         ((step.gt.x(ic)).and.(abs(x(ic)-x(ia)).gt.prec).and.
     $         (abs(f(ic)-f(ib)).gt.(tol*f(ib))))) then
            if (g(ib).gt.0) then
              step = (x(ia)+x(ib))/2.d0
            else
              step = (x(ib)+x(ic))/2.d0
            endif
            write(6,772) step
 772        format('Step outside bracket -- bisection step',f20.10)
          endif
        else
c
c  Cannot bracket, assume function is monotonic 
c  and use points with least function value
c
          write(6,774) ib,x(ib),f(ib)
 774      format('Minimum: ',i5,2f20.10)
          if ((ia.ne.0).and.(g(ib).lt.0.0d0)) then
            dx = (g(ia) - g(ib))/(x(ia)-x(ib))
            write(6,775) ia,x(ia),f(ia),dx
 775        format('Extrapolate right:',i5,2f20.10,5x,
     $             '2nd deriv',2x,e12.3)
            step = x(ib) + (x(ib) - x(1))
          else if ((ic.ne.0).and.(g(ib).gt.0.d0)) then
            dx = (g(ic) - g(ib))/(x(ic)-x(ib))
            write(6,776) ic,x(ic),f(ic),dx
 776        format('Extrapolate left:',i5,2f20.10,5x,
     $             '2nd deriv',2x,e12.3)
            step = x(ib) - (x(1) - x(ib))
          else
            write(6,773)
 773        format('Internal error....')
          endif
        endif
      else if (n.eq.1) then
        step = 1.d0
      endif
      linesrch = step
      write(6,881)
 881  format(//)
      return
      end






      logical function linebracket(n,x,f,ia,ib,ic)
      implicit none
      integer n
      double precision x(n),f(n)
      integer ia,ib,ic
      double precision xmin, ymin, del
      integer i

      ia = 0
      ic = 0
      ib = 1
      xmin = f(1)
      do i=2,n
        if (f(i).lt.xmin) then
          xmin = f(i)
          ib = i
        endif
      enddo
      xmin = -1.d0*abs(xmin)
      ymin = abs(xmin)
      do i=1,n
        del = x(i) - x(ib)
        if ((del.gt.0.d0).and.(del.lt.ymin)) then
          ymin = del
          ic = i
        else if ((del.lt.0.d0).and.(del.gt.xmin)) then
          xmin = del
          ia = i
        endif
      enddo
      linebracket = ((ia.ne.0).and.(ic.ne.0))
      return
      end










      double precision function linesrch1(n,x,f,g)
      implicit none
#include "global.fh"
      integer n
      double precision x(n)
      double precision f(n)
      double precision g(n)
      double precision a, b1, b2, c1, c2, f1, f2, x0
      double precision e1, e2, err
      save e1, e2
      data err/0.d0/
      
      if (n.eq.1) then
        linesrch1 = 1.d0
        return
      else
        if (n.gt.2) then
          err = max(abs((f(n)-e1)/f(n)),abs((f(n)-e2)/f(n)))
C          if (ga_nodeid().eq.0) write(6,901) x(n),e1,e2,f(n),err
 901      format(//,'Predictions @ ',f12.6,5x,2f20.10,
     $           /,'Obs.',f20.10,
     $           /,'Relative error:',e12.2)
        endif
        a = (g(n)-g(n-1))/(2.d0*(x(n)-x(n-1)))
        b1 = g(n-1) - 4.d0*a*x(n-1)
        b2 = g(n) - 4.d0*a*x(n)
        c1 = f(n-1) - a*x(n-1)**2 - b1*x(n-1)
        c2 = f(n) - a*x(n)**2 - b2*x(n)
        f1 = a*x(n-1)**2 + b1*x(n-1) + c1
        f2 = a*x(n)**2 + b2*x(n) + c2
C        if (ga_nodeid().eq.0)
C    $     write(6,900) a,b1,b2,c1,c2,f1,f(n-1),f2,f(n)
 900    format('A:   ',5x,f20.10,/,
     $         'B(1):',5x,f20.10,/,
     $         'B(2):',5x,f20.10,/,
     $         'C(1):',5x,f20.10,/,
     $         'C(2):',5x,f20.10,/,
     $         'F(1):',5x,f20.10,'    f(1):',f20.10,/,
     $         'F(2):',5x,f20.10,'    f(2):',f20.10)
C        if (err.lt.1.d-1) then
          x0 = (g(n)*x(n-1) - g(n-1)*x(n))/(g(n) - g(n-1))
          e1 = a*x0**2 + b1*x0 + c1
          e2 = a*x0**2 + b2*x0 + c2
          linesrch1 = x0
C        else
          
      endif
      return
      end


        
