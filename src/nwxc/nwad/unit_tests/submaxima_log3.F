      subroutine submaxima_log3(x,y,z,f,df,df2,df3) 
      implicit none 
      double precision x,y,z,f,df(3),df2(3*4/2),df3(3*4*5/6) 
      f = log(3.333333333333333d-1*(sin(z)+sin(y)+sin(x)))
      df(1) = 1.0d+0*cos(x)/(sin(z)+sin(y)+sin(x))
      df(2) = 1.0d+0*cos(y)/(sin(z)+sin(y)+sin(x))
      df(3) = 1.0d+0*cos(z)/(sin(z)+sin(y)+sin(x))
      df2(1) = -1.0d+0*sin(x)/(sin(z)+sin(y)+sin(x))-1.0d+0*cos(x)**2/(s
     1   in(z)+sin(y)+sin(x))**2
      df2(2) = -1.0d+0*cos(x)*cos(y)/(sin(z)+sin(y)+sin(x))**2
      df2(3) = -1.0d+0*sin(y)/(sin(z)+sin(y)+sin(x))-1.0d+0*cos(y)**2/(s
     1   in(z)+sin(y)+sin(x))**2
      df2(4) = -1.0d+0*cos(x)*cos(z)/(sin(z)+sin(y)+sin(x))**2
      df2(5) = -1.0d+0*cos(y)*cos(z)/(sin(z)+sin(y)+sin(x))**2
      df2(6) = -1.0d+0*sin(z)/(sin(z)+sin(y)+sin(x))-1.0d+0*cos(z)**2/(s
     1   in(z)+sin(y)+sin(x))**2
      df3(1) = -1.0d+0*cos(x)/(sin(z)+sin(y)+sin(x))+3.0d+0*cos(x)*sin(x
     1   )/(sin(z)+sin(y)+sin(x))**2+2.0d+0*cos(x)**3/(sin(z)+sin(y)+sin
     2   (x))**3
      df3(2) = 1.0d+0*sin(x)*cos(y)/(sin(z)+sin(y)+sin(x))**2+2.0d+0*cos
     1   (x)**2*cos(y)/(sin(z)+sin(y)+sin(x))**3
      df3(3) = 1.0d+0*cos(x)*sin(y)/(sin(z)+sin(y)+sin(x))**2+2.0d+0*cos
     1   (x)*cos(y)**2/(sin(z)+sin(y)+sin(x))**3
      df3(4) = -1.0d+0*cos(y)/(sin(z)+sin(y)+sin(x))+3.0d+0*cos(y)*sin(y
     1   )/(sin(z)+sin(y)+sin(x))**2+2.0d+0*cos(y)**3/(sin(z)+sin(y)+sin
     2   (x))**3
      df3(5) = 1.0d+0*sin(x)*cos(z)/(sin(z)+sin(y)+sin(x))**2+2.0d+0*cos
     1   (x)**2*cos(z)/(sin(z)+sin(y)+sin(x))**3
      df3(6) = 2.0d+0*cos(x)*cos(y)*cos(z)/(sin(z)+sin(y)+sin(x))**3
      df3(7) = 1.0d+0*sin(y)*cos(z)/(sin(z)+sin(y)+sin(x))**2+2.0d+0*cos
     1   (y)**2*cos(z)/(sin(z)+sin(y)+sin(x))**3
      df3(8) = 1.0d+0*cos(x)*sin(z)/(sin(z)+sin(y)+sin(x))**2+2.0d+0*cos
     1   (x)*cos(z)**2/(sin(z)+sin(y)+sin(x))**3
      df3(9) = 1.0d+0*cos(y)*sin(z)/(sin(z)+sin(y)+sin(x))**2+2.0d+0*cos
     1   (y)*cos(z)**2/(sin(z)+sin(y)+sin(x))**3
      df3(10) = -1.0d+0*cos(z)/(sin(z)+sin(y)+sin(x))+3.0d+0*cos(z)*sin(
     1   z)/(sin(z)+sin(y)+sin(x))**2+2.0d+0*cos(z)**3/(sin(z)+sin(y)+si
     2   n(x))**3
      end 
