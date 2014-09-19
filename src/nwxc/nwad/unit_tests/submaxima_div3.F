      subroutine submaxima_div3(x,y,z,f,df,df2,df3) 
      implicit none 
      double precision x,y,z,f,df(3),df2(3*4/2),df3(3*4*5/6) 
      f = sin(z)/sin(y)+sin(z)/sin(x)+sin(y)/sin(z)+sin(x)/sin(z)+sin(y)
     1   /sin(x)+sin(x)/sin(y)
      df(1) = -cos(x)*sin(z)/sin(x)**2+cos(x)/sin(z)-cos(x)*sin(y)/sin(x
     1   )**2+cos(x)/sin(y)
      df(2) = -cos(y)*sin(z)/sin(y)**2+cos(y)/sin(z)-sin(x)*cos(y)/sin(y
     1   )**2+cos(y)/sin(x)
      df(3) = -sin(y)*cos(z)/sin(z)**2-sin(x)*cos(z)/sin(z)**2+cos(z)/si
     1   n(y)+cos(z)/sin(x)
      df2(1) = sin(z)/sin(x)+2*cos(x)**2*sin(z)/sin(x)**3-sin(x)/sin(z)+
     1   sin(y)/sin(x)+2*cos(x)**2*sin(y)/sin(x)**3-sin(x)/sin(y)
      df2(2) = -cos(x)*cos(y)/sin(y)**2-cos(x)*cos(y)/sin(x)**2
      df2(3) = sin(z)/sin(y)+2*cos(y)**2*sin(z)/sin(y)**3-sin(y)/sin(z)-
     1   sin(y)/sin(x)+sin(x)/sin(y)+2*sin(x)*cos(y)**2/sin(y)**3
      df2(4) = -cos(x)*cos(z)/sin(z)**2-cos(x)*cos(z)/sin(x)**2
      df2(5) = -cos(y)*cos(z)/sin(z)**2-cos(y)*cos(z)/sin(y)**2
      df2(6) = -sin(z)/sin(y)-sin(z)/sin(x)+sin(y)/sin(z)+sin(x)/sin(z)+
     1   2*sin(y)*cos(z)**2/sin(z)**3+2*sin(x)*cos(z)**2/sin(z)**3
      df3(1) = -5*cos(x)*sin(z)/sin(x)**2-6*cos(x)**3*sin(z)/sin(x)**4-c
     1   os(x)/sin(z)-5*cos(x)*sin(y)/sin(x)**2-6*cos(x)**3*sin(y)/sin(x
     2   )**4-cos(x)/sin(y)
      df3(2) = sin(x)*cos(y)/sin(y)**2+cos(y)/sin(x)+2*cos(x)**2*cos(y)/
     1   sin(x)**3
      df3(3) = cos(x)*sin(y)/sin(x)**2+cos(x)/sin(y)+2*cos(x)*cos(y)**2/
     1   sin(y)**3
      df3(4) = -5*cos(y)*sin(z)/sin(y)**2-6*cos(y)**3*sin(z)/sin(y)**4-c
     1   os(y)/sin(z)-5*sin(x)*cos(y)/sin(y)**2-6*sin(x)*cos(y)**3/sin(y
     2   )**4-cos(y)/sin(x)
      df3(5) = sin(x)*cos(z)/sin(z)**2+cos(z)/sin(x)+2*cos(x)**2*cos(z)/
     1   sin(x)**3
      df3(6) = 0
      df3(7) = sin(y)*cos(z)/sin(z)**2+cos(z)/sin(y)+2*cos(y)**2*cos(z)/
     1   sin(y)**3
      df3(8) = cos(x)*sin(z)/sin(x)**2+cos(x)/sin(z)+2*cos(x)*cos(z)**2/
     1   sin(z)**3
      df3(9) = cos(y)*sin(z)/sin(y)**2+cos(y)/sin(z)+2*cos(y)*cos(z)**2/
     1   sin(z)**3
      df3(10) = -5*sin(y)*cos(z)/sin(z)**2-5*sin(x)*cos(z)/sin(z)**2-6*s
     1   in(y)*cos(z)**3/sin(z)**4-6*sin(x)*cos(z)**3/sin(z)**4-cos(z)/s
     2   in(y)-cos(z)/sin(x)
      end 
