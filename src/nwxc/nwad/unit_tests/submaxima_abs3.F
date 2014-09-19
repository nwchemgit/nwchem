      subroutine submaxima_abs3(x,y,z,f,df,df2,df3) 
      implicit none 
      double precision x,y,z,f,df(3),df2(3*4/2),df3(3*4*5/6) 
      f = 3.333333333333333d-1*abs(sin(z)+sin(y)+sin(x))
      df(1) = 3.333333333333333d-1*cos(x)*(sin(z)+sin(y)+sin(x))/abs(sin
     1   (z)+sin(y)+sin(x))
      df(2) = 3.333333333333333d-1*cos(y)*(sin(z)+sin(y)+sin(x))/abs(sin
     1   (z)+sin(y)+sin(x))
      df(3) = 3.333333333333333d-1*cos(z)*(sin(z)+sin(y)+sin(x))/abs(sin
     1   (z)+sin(y)+sin(x))
      df2(1) = -3.333333333333333d-1*sin(x)*(sin(z)+sin(y)+sin(x))/abs(s
     1   in(z)+sin(y)+sin(x))
      df2(2) = 0.0d+0
      df2(3) = -3.333333333333333d-1*sin(y)*(sin(z)+sin(y)+sin(x))/abs(s
     1   in(z)+sin(y)+sin(x))
      df2(4) = 0.0d+0
      df2(5) = 0.0d+0
      df2(6) = -3.333333333333333d-1*sin(z)*(sin(z)+sin(y)+sin(x))/abs(s
     1   in(z)+sin(y)+sin(x))
      df3(1) = -3.333333333333333d-1*cos(x)*(sin(z)+sin(y)+sin(x))/abs(s
     1   in(z)+sin(y)+sin(x))
      df3(2) = 0.0d+0
      df3(3) = 0
      df3(4) = -3.333333333333333d-1*cos(y)*(sin(z)+sin(y)+sin(x))/abs(s
     1   in(z)+sin(y)+sin(x))
      df3(5) = 0.0d+0
      df3(6) = 0
      df3(7) = 0.0d+0
      df3(8) = 0
      df3(9) = 0
      df3(10) = -3.333333333333333d-1*cos(z)*(sin(z)+sin(y)+sin(x))/abs(
     1   sin(z)+sin(y)+sin(x))
      end 
