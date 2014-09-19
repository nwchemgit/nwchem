      subroutine submaxima_addx3(x,y,z,a,f,df,df2,df3) 
      implicit none 
      double precision x,y,z,a,f,df(3),df2(3*4/2),df3(3*4*5/6) 
      f = sin(z)+sin(y)+sin(x)+a
      df(1) = cos(x)
      df(2) = cos(y)
      df(3) = cos(z)
      df2(1) = -sin(x)
      df2(2) = 0
      df2(3) = -sin(y)
      df2(4) = 0
      df2(5) = 0
      df2(6) = -sin(z)
      df3(1) = -cos(x)
      df3(2) = 0
      df3(3) = 0
      df3(4) = -cos(y)
      df3(5) = 0
      df3(6) = 0
      df3(7) = 0
      df3(8) = 0
      df3(9) = 0
      df3(10) = -cos(z)
      end 
