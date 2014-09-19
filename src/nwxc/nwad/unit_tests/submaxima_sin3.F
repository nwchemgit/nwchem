      subroutine submaxima_sin3(x,y,z,f,df,df2,df3) 
      implicit none 
      double precision x,y,z,f,df(3),df2(3*4/2),df3(3*4*5/6) 
      f = sin(sin(z)+sin(y)+sin(x))
      df(1) = cos(x)*cos(sin(z)+sin(y)+sin(x))
      df(2) = cos(y)*cos(sin(z)+sin(y)+sin(x))
      df(3) = cos(z)*cos(sin(z)+sin(y)+sin(x))
      df2(1) = -cos(x)**2*sin(sin(z)+sin(y)+sin(x))-sin(x)*cos(sin(z)+si
     1   n(y)+sin(x))
      df2(2) = -cos(x)*cos(y)*sin(sin(z)+sin(y)+sin(x))
      df2(3) = -cos(y)**2*sin(sin(z)+sin(y)+sin(x))-sin(y)*cos(sin(z)+si
     1   n(y)+sin(x))
      df2(4) = -cos(x)*cos(z)*sin(sin(z)+sin(y)+sin(x))
      df2(5) = -cos(y)*cos(z)*sin(sin(z)+sin(y)+sin(x))
      df2(6) = -cos(z)**2*sin(sin(z)+sin(y)+sin(x))-sin(z)*cos(sin(z)+si
     1   n(y)+sin(x))
      df3(1) = 3*cos(x)*sin(x)*sin(sin(z)+sin(y)+sin(x))-cos(x)**3*cos(s
     1   in(z)+sin(y)+sin(x))-cos(x)*cos(sin(z)+sin(y)+sin(x))
      df3(2) = sin(x)*cos(y)*sin(sin(z)+sin(y)+sin(x))-cos(x)**2*cos(y)*
     1   cos(sin(z)+sin(y)+sin(x))
      df3(3) = cos(x)*sin(y)*sin(sin(z)+sin(y)+sin(x))-cos(x)*cos(y)**2*
     1   cos(sin(z)+sin(y)+sin(x))
      df3(4) = 3*cos(y)*sin(y)*sin(sin(z)+sin(y)+sin(x))-cos(y)**3*cos(s
     1   in(z)+sin(y)+sin(x))-cos(y)*cos(sin(z)+sin(y)+sin(x))
      df3(5) = sin(x)*cos(z)*sin(sin(z)+sin(y)+sin(x))-cos(x)**2*cos(z)*
     1   cos(sin(z)+sin(y)+sin(x))
      df3(6) = -cos(x)*cos(y)*cos(z)*cos(sin(z)+sin(y)+sin(x))
      df3(7) = sin(y)*cos(z)*sin(sin(z)+sin(y)+sin(x))-cos(y)**2*cos(z)*
     1   cos(sin(z)+sin(y)+sin(x))
      df3(8) = cos(x)*sin(z)*sin(sin(z)+sin(y)+sin(x))-cos(x)*cos(z)**2*
     1   cos(sin(z)+sin(y)+sin(x))
      df3(9) = cos(y)*sin(z)*sin(sin(z)+sin(y)+sin(x))-cos(y)*cos(z)**2*
     1   cos(sin(z)+sin(y)+sin(x))
      df3(10) = 3*cos(z)*sin(z)*sin(sin(z)+sin(y)+sin(x))-cos(z)**3*cos(
     1   sin(z)+sin(y)+sin(x))-cos(z)*cos(sin(z)+sin(y)+sin(x))
      end 
