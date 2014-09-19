      subroutine submaxima_sinh3(x,y,z,f,df,df2,df3) 
      implicit none 
      double precision x,y,z,f,df(3),df2(3*4/2),df3(3*4*5/6) 
      f = sinh(sin(z)+sin(y)+sin(x)+4.0d+0)
      df(1) = cos(x)*cosh(sin(z)+sin(y)+sin(x)+4.0d+0)
      df(2) = cos(y)*cosh(sin(z)+sin(y)+sin(x)+4.0d+0)
      df(3) = cos(z)*cosh(sin(z)+sin(y)+sin(x)+4.0d+0)
      df2(1) = cos(x)**2*sinh(sin(z)+sin(y)+sin(x)+4.0d+0)-sin(x)*cosh(s
     1   in(z)+sin(y)+sin(x)+4.0d+0)
      df2(2) = cos(x)*cos(y)*sinh(sin(z)+sin(y)+sin(x)+4.0d+0)
      df2(3) = cos(y)**2*sinh(sin(z)+sin(y)+sin(x)+4.0d+0)-sin(y)*cosh(s
     1   in(z)+sin(y)+sin(x)+4.0d+0)
      df2(4) = cos(x)*cos(z)*sinh(sin(z)+sin(y)+sin(x)+4.0d+0)
      df2(5) = cos(y)*cos(z)*sinh(sin(z)+sin(y)+sin(x)+4.0d+0)
      df2(6) = cos(z)**2*sinh(sin(z)+sin(y)+sin(x)+4.0d+0)-sin(z)*cosh(s
     1   in(z)+sin(y)+sin(x)+4.0d+0)
      df3(1) = -3*cos(x)*sin(x)*sinh(sin(z)+sin(y)+sin(x)+4.0d+0)+cos(x)
     1   **3*cosh(sin(z)+sin(y)+sin(x)+4.0d+0)-cos(x)*cosh(sin(z)+sin(y)
     2   +sin(x)+4.0d+0)
      df3(2) = cos(x)**2*cos(y)*cosh(sin(z)+sin(y)+sin(x)+4.0d+0)-sin(x)
     1   *cos(y)*sinh(sin(z)+sin(y)+sin(x)+4.0d+0)
      df3(3) = cos(x)*cos(y)**2*cosh(sin(z)+sin(y)+sin(x)+4.0d+0)-cos(x)
     1   *sin(y)*sinh(sin(z)+sin(y)+sin(x)+4.0d+0)
      df3(4) = -3*cos(y)*sin(y)*sinh(sin(z)+sin(y)+sin(x)+4.0d+0)+cos(y)
     1   **3*cosh(sin(z)+sin(y)+sin(x)+4.0d+0)-cos(y)*cosh(sin(z)+sin(y)
     2   +sin(x)+4.0d+0)
      df3(5) = cos(x)**2*cos(z)*cosh(sin(z)+sin(y)+sin(x)+4.0d+0)-sin(x)
     1   *cos(z)*sinh(sin(z)+sin(y)+sin(x)+4.0d+0)
      df3(6) = cos(x)*cos(y)*cos(z)*cosh(sin(z)+sin(y)+sin(x)+4.0d+0)
      df3(7) = cos(y)**2*cos(z)*cosh(sin(z)+sin(y)+sin(x)+4.0d+0)-sin(y)
     1   *cos(z)*sinh(sin(z)+sin(y)+sin(x)+4.0d+0)
      df3(8) = cos(x)*cos(z)**2*cosh(sin(z)+sin(y)+sin(x)+4.0d+0)-cos(x)
     1   *sin(z)*sinh(sin(z)+sin(y)+sin(x)+4.0d+0)
      df3(9) = cos(y)*cos(z)**2*cosh(sin(z)+sin(y)+sin(x)+4.0d+0)-cos(y)
     1   *sin(z)*sinh(sin(z)+sin(y)+sin(x)+4.0d+0)
      df3(10) = -3*cos(z)*sin(z)*sinh(sin(z)+sin(y)+sin(x)+4.0d+0)+cos(z
     1   )**3*cosh(sin(z)+sin(y)+sin(x)+4.0d+0)-cos(z)*cosh(sin(z)+sin(y
     2   )+sin(x)+4.0d+0)
      end 
