!234567890
      double precision function ei(x,err)
      implicit none
      double precision x,err
      double precision gamma
      parameter (gamma=0.57721566490153286060606512090082402431042d0)
      double precision eisum,xfact,xterm,xtermold
      integer ii

      if (abs(x).lt.gamma*err) then
        ei=gamma+log(abs(x))
      elseif (abs(x).lt.abs(log(err))) then
        eisum=0.d0
        xfact=1.d0
        ii=0
 100    continue
          ii=ii+1
          xfact=xfact*x/dble(ii)
          xterm=xfact/dble(ii)
          eisum=eisum+xterm
          if (abs(xterm).lt.eisum*err) goto 200
          if (abs(xterm).lt.gamma*err) goto 200
        goto 100
 200    continue
        ei=eisum+gamma+log(abs(x))
      else
        eisum=1.d0
        xterm=1.d0
        ii=0
 300    continue
          ii=ii+1
          xtermold=xterm
          xterm=xterm*ii/x
          if (xterm.lt.err) goto 400
          if (xterm.lt.xtermold) then
            eisum=eisum+xterm
          else
            eisum=eisum-xtermold
            goto 400
          endif
        goto 300
 400    continue
        ei=exp(x)*eisum/x
      endif

      return
      end
