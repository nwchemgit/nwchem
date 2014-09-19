      subroutine submaxima_cos3(x,y,z,f,df,df2,df3) 
      implicit none 
      double precision x,y,z,f,df(3),df2(3*4/2),df3(3*4*5/6) 
      f = cos(cos(z)+cos(y)+cos(x))
      df(1) = sin(x)*sin(cos(z)+cos(y)+cos(x))
      df(2) = sin(y)*sin(cos(z)+cos(y)+cos(x))
      df(3) = sin(z)*sin(cos(z)+cos(y)+cos(x))
      df2(1) = cos(x)*sin(cos(z)+cos(y)+cos(x))-sin(x)**2*cos(cos(z)+cos
     1   (y)+cos(x))
      df2(2) = -sin(x)*sin(y)*cos(cos(z)+cos(y)+cos(x))
      df2(3) = cos(y)*sin(cos(z)+cos(y)+cos(x))-sin(y)**2*cos(cos(z)+cos
     1   (y)+cos(x))
      df2(4) = -sin(x)*sin(z)*cos(cos(z)+cos(y)+cos(x))
      df2(5) = -sin(y)*sin(z)*cos(cos(z)+cos(y)+cos(x))
      df2(6) = cos(z)*sin(cos(z)+cos(y)+cos(x))-sin(z)**2*cos(cos(z)+cos
     1   (y)+cos(x))
      df3(1) = -sin(x)**3*sin(cos(z)+cos(y)+cos(x))-sin(x)*sin(cos(z)+co
     1   s(y)+cos(x))-3*cos(x)*sin(x)*cos(cos(z)+cos(y)+cos(x))
      df3(2) = -sin(x)**2*sin(y)*sin(cos(z)+cos(y)+cos(x))-cos(x)*sin(y)
     1   *cos(cos(z)+cos(y)+cos(x))
      df3(3) = -sin(x)*sin(y)**2*sin(cos(z)+cos(y)+cos(x))-sin(x)*cos(y)
     1   *cos(cos(z)+cos(y)+cos(x))
      df3(4) = -sin(y)**3*sin(cos(z)+cos(y)+cos(x))-sin(y)*sin(cos(z)+co
     1   s(y)+cos(x))-3*cos(y)*sin(y)*cos(cos(z)+cos(y)+cos(x))
      df3(5) = -sin(x)**2*sin(z)*sin(cos(z)+cos(y)+cos(x))-cos(x)*sin(z)
     1   *cos(cos(z)+cos(y)+cos(x))
      df3(6) = -sin(x)*sin(y)*sin(z)*sin(cos(z)+cos(y)+cos(x))
      df3(7) = -sin(y)**2*sin(z)*sin(cos(z)+cos(y)+cos(x))-cos(y)*sin(z)
     1   *cos(cos(z)+cos(y)+cos(x))
      df3(8) = -sin(x)*sin(z)**2*sin(cos(z)+cos(y)+cos(x))-sin(x)*cos(z)
     1   *cos(cos(z)+cos(y)+cos(x))
      df3(9) = -sin(y)*sin(z)**2*sin(cos(z)+cos(y)+cos(x))-sin(y)*cos(z)
     1   *cos(cos(z)+cos(y)+cos(x))
      df3(10) = -sin(z)**3*sin(cos(z)+cos(y)+cos(x))-sin(z)*sin(cos(z)+c
     1   os(y)+cos(x))-3*cos(z)*sin(z)*cos(cos(z)+cos(y)+cos(x))
      end 
