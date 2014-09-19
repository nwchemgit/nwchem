      subroutine submaxima_cosh3(x,y,z,f,df,df2,df3) 
      implicit none 
      double precision x,y,z,f,df(3),df2(3*4/2),df3(3*4*5/6) 
      f = cosh(sin(z)+sin(y)+sin(x)+4.0d+0)
      df(1) = cos(x)*sinh(sin(z)+sin(y)+sin(x)+4.0d+0)
      df(2) = cos(y)*sinh(sin(z)+sin(y)+sin(x)+4.0d+0)
      df(3) = cos(z)*sinh(sin(z)+sin(y)+sin(x)+4.0d+0)
      df2(1) = cos(x)**2*cosh(sin(z)+sin(y)+sin(x)+4.0d+0)-sin(x)*sinh(s
     1   in(z)+sin(y)+sin(x)+4.0d+0)
      df2(2) = cos(x)*cos(y)*cosh(sin(z)+sin(y)+sin(x)+4.0d+0)
      df2(3) = cos(y)**2*cosh(sin(z)+sin(y)+sin(x)+4.0d+0)-sin(y)*sinh(s
     1   in(z)+sin(y)+sin(x)+4.0d+0)
      df2(4) = cos(x)*cos(z)*cosh(sin(z)+sin(y)+sin(x)+4.0d+0)
      df2(5) = cos(y)*cos(z)*cosh(sin(z)+sin(y)+sin(x)+4.0d+0)
      df2(6) = cos(z)**2*cosh(sin(z)+sin(y)+sin(x)+4.0d+0)-sin(z)*sinh(s
     1   in(z)+sin(y)+sin(x)+4.0d+0)
      df3(1) = cos(x)**3*sinh(sin(z)+sin(y)+sin(x)+4.0d+0)-cos(x)*sinh(s
     1   in(z)+sin(y)+sin(x)+4.0d+0)-3*cos(x)*sin(x)*cosh(sin(z)+sin(y)+
     2   sin(x)+4.0d+0)
      df3(2) = cos(x)**2*cos(y)*sinh(sin(z)+sin(y)+sin(x)+4.0d+0)-sin(x)
     1   *cos(y)*cosh(sin(z)+sin(y)+sin(x)+4.0d+0)
      df3(3) = cos(x)*cos(y)**2*sinh(sin(z)+sin(y)+sin(x)+4.0d+0)-cos(x)
     1   *sin(y)*cosh(sin(z)+sin(y)+sin(x)+4.0d+0)
      df3(4) = cos(y)**3*sinh(sin(z)+sin(y)+sin(x)+4.0d+0)-cos(y)*sinh(s
     1   in(z)+sin(y)+sin(x)+4.0d+0)-3*cos(y)*sin(y)*cosh(sin(z)+sin(y)+
     2   sin(x)+4.0d+0)
      df3(5) = cos(x)**2*cos(z)*sinh(sin(z)+sin(y)+sin(x)+4.0d+0)-sin(x)
     1   *cos(z)*cosh(sin(z)+sin(y)+sin(x)+4.0d+0)
      df3(6) = cos(x)*cos(y)*cos(z)*sinh(sin(z)+sin(y)+sin(x)+4.0d+0)
      df3(7) = cos(y)**2*cos(z)*sinh(sin(z)+sin(y)+sin(x)+4.0d+0)-sin(y)
     1   *cos(z)*cosh(sin(z)+sin(y)+sin(x)+4.0d+0)
      df3(8) = cos(x)*cos(z)**2*sinh(sin(z)+sin(y)+sin(x)+4.0d+0)-cos(x)
     1   *sin(z)*cosh(sin(z)+sin(y)+sin(x)+4.0d+0)
      df3(9) = cos(y)*cos(z)**2*sinh(sin(z)+sin(y)+sin(x)+4.0d+0)-cos(y)
     1   *sin(z)*cosh(sin(z)+sin(y)+sin(x)+4.0d+0)
      df3(10) = cos(z)**3*sinh(sin(z)+sin(y)+sin(x)+4.0d+0)-cos(z)*sinh(
     1   sin(z)+sin(y)+sin(x)+4.0d+0)-3*cos(z)*sin(z)*cosh(sin(z)+sin(y)
     2   +sin(x)+4.0d+0)
      end 
