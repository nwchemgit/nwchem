      double precision function double_factorial(n)

      integer n
      integer n11(-1:16)
      data (n11(i), i=-1,16) /1,1,1,2,3,8,15,48,105,384,945,
     >                        3840,10395,46080,135135,645120,
     >                        2027025,10321920/
                                         
      if( n.ge.(-1) .and. n.le.16) then
        double_factorial = n11(n)
      else
        call errquit("too big parameter in double_factorial",1,0)
      end if

      end
c $Id$
