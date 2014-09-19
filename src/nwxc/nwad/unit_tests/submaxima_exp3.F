      subroutine submaxima_exp3(x,y,z,f,df,df2,df3) 
      implicit none 
      double precision x,y,z,f,df(3),df2(3*4/2),df3(3*4*5/6) 
      f = exp(sin(z)+sin(y)+sin(x))
      df(1) = cos(x)*exp(sin(z)+sin(y)+sin(x))
      df(2) = cos(y)*exp(sin(z)+sin(y)+sin(x))
      df(3) = exp(sin(z)+sin(y)+sin(x))*cos(z)
      df2(1) = cos(x)**2*exp(sin(z)+sin(y)+sin(x))-sin(x)*exp(sin(z)+sin
     1   (y)+sin(x))
      df2(2) = cos(x)*cos(y)*exp(sin(z)+sin(y)+sin(x))
      df2(3) = cos(y)**2*exp(sin(z)+sin(y)+sin(x))-sin(y)*exp(sin(z)+sin
     1   (y)+sin(x))
      df2(4) = cos(x)*exp(sin(z)+sin(y)+sin(x))*cos(z)
      df2(5) = cos(y)*exp(sin(z)+sin(y)+sin(x))*cos(z)
      df2(6) = exp(sin(z)+sin(y)+sin(x))*cos(z)**2-exp(sin(z)+sin(y)+sin
     1   (x))*sin(z)
      df3(1) = -3*cos(x)*sin(x)*exp(sin(z)+sin(y)+sin(x))+cos(x)**3*exp(
     1   sin(z)+sin(y)+sin(x))-cos(x)*exp(sin(z)+sin(y)+sin(x))
      df3(2) = cos(x)**2*cos(y)*exp(sin(z)+sin(y)+sin(x))-sin(x)*cos(y)*
     1   exp(sin(z)+sin(y)+sin(x))
      df3(3) = cos(x)*cos(y)**2*exp(sin(z)+sin(y)+sin(x))-cos(x)*sin(y)*
     1   exp(sin(z)+sin(y)+sin(x))
      df3(4) = -3*cos(y)*sin(y)*exp(sin(z)+sin(y)+sin(x))+cos(y)**3*exp(
     1   sin(z)+sin(y)+sin(x))-cos(y)*exp(sin(z)+sin(y)+sin(x))
      df3(5) = cos(x)**2*exp(sin(z)+sin(y)+sin(x))*cos(z)-sin(x)*exp(sin
     1   (z)+sin(y)+sin(x))*cos(z)
      df3(6) = cos(x)*cos(y)*exp(sin(z)+sin(y)+sin(x))*cos(z)
      df3(7) = cos(y)**2*exp(sin(z)+sin(y)+sin(x))*cos(z)-sin(y)*exp(sin
     1   (z)+sin(y)+sin(x))*cos(z)
      df3(8) = cos(x)*exp(sin(z)+sin(y)+sin(x))*cos(z)**2-cos(x)*exp(sin
     1   (z)+sin(y)+sin(x))*sin(z)
      df3(9) = cos(y)*exp(sin(z)+sin(y)+sin(x))*cos(z)**2-cos(y)*exp(sin
     1   (z)+sin(y)+sin(x))*sin(z)
      df3(10) = -3*exp(sin(z)+sin(y)+sin(x))*cos(z)*sin(z)+exp(sin(z)+si
     1   n(y)+sin(x))*cos(z)**3-exp(sin(z)+sin(y)+sin(x))*cos(z)
      end 
