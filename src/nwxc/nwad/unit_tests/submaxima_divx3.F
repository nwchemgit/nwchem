      subroutine submaxima_divx3(x,y,z,a,f,df,df2,df3) 
      implicit none 
      double precision x,y,z,a,f,df(3),df2(3*4/2),df3(3*4*5/6) 
      f = sin(x)*sin(y)*sin(z)/a
      df(1) = cos(x)*sin(y)*sin(z)/a
      df(2) = sin(x)*cos(y)*sin(z)/a
      df(3) = sin(x)*sin(y)*cos(z)/a
      df2(1) = -sin(x)*sin(y)*sin(z)/a
      df2(2) = cos(x)*cos(y)*sin(z)/a
      df2(3) = -sin(x)*sin(y)*sin(z)/a
      df2(4) = cos(x)*sin(y)*cos(z)/a
      df2(5) = sin(x)*cos(y)*cos(z)/a
      df2(6) = -sin(x)*sin(y)*sin(z)/a
      df3(1) = -cos(x)*sin(y)*sin(z)/a
      df3(2) = -sin(x)*cos(y)*sin(z)/a
      df3(3) = -cos(x)*sin(y)*sin(z)/a
      df3(4) = -sin(x)*cos(y)*sin(z)/a
      df3(5) = -sin(x)*sin(y)*cos(z)/a
      df3(6) = cos(x)*cos(y)*cos(z)/a
      df3(7) = -sin(x)*sin(y)*cos(z)/a
      df3(8) = -cos(x)*sin(y)*sin(z)/a
      df3(9) = -sin(x)*cos(y)*sin(z)/a
      df3(10) = -sin(x)*sin(y)*cos(z)/a
      end 
