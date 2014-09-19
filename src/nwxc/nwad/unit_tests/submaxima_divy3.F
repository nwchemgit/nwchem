      subroutine submaxima_divy3(x,y,z,a,f,df,df2,df3) 
      implicit none 
      double precision x,y,z,a,f,df(3),df2(3*4/2),df3(3*4*5/6) 
      f = a/(sin(x)*sin(y)*sin(z))
      df(1) = -a*cos(x)/(sin(x)**2*sin(y)*sin(z))
      df(2) = -a*cos(y)/(sin(x)*sin(y)**2*sin(z))
      df(3) = -a*cos(z)/(sin(x)*sin(y)*sin(z)**2)
      df2(1) = a/(sin(x)*sin(y)*sin(z))+2*a*cos(x)**2/(sin(x)**3*sin(y)*
     1   sin(z))
      df2(2) = a*cos(x)*cos(y)/(sin(x)**2*sin(y)**2*sin(z))
      df2(3) = a/(sin(x)*sin(y)*sin(z))+2*a*cos(y)**2/(sin(x)*sin(y)**3*
     1   sin(z))
      df2(4) = a*cos(x)*cos(z)/(sin(x)**2*sin(y)*sin(z)**2)
      df2(5) = a*cos(y)*cos(z)/(sin(x)*sin(y)**2*sin(z)**2)
      df2(6) = a/(sin(x)*sin(y)*sin(z))+2*a*cos(z)**2/(sin(x)*sin(y)*sin
     1   (z)**3)
      df3(1) = -5*a*cos(x)/(sin(x)**2*sin(y)*sin(z))-6*a*cos(x)**3/(sin(
     1   x)**4*sin(y)*sin(z))
      df3(2) = -a*cos(y)/(sin(x)*sin(y)**2*sin(z))-2*a*cos(x)**2*cos(y)/
     1   (sin(x)**3*sin(y)**2*sin(z))
      df3(3) = -a*cos(x)/(sin(x)**2*sin(y)*sin(z))-2*a*cos(x)*cos(y)**2/
     1   (sin(x)**2*sin(y)**3*sin(z))
      df3(4) = -5*a*cos(y)/(sin(x)*sin(y)**2*sin(z))-6*a*cos(y)**3/(sin(
     1   x)*sin(y)**4*sin(z))
      df3(5) = -a*cos(z)/(sin(x)*sin(y)*sin(z)**2)-2*a*cos(x)**2*cos(z)/
     1   (sin(x)**3*sin(y)*sin(z)**2)
      df3(6) = -a*cos(x)*cos(y)*cos(z)/(sin(x)**2*sin(y)**2*sin(z)**2)
      df3(7) = -a*cos(z)/(sin(x)*sin(y)*sin(z)**2)-2*a*cos(y)**2*cos(z)/
     1   (sin(x)*sin(y)**3*sin(z)**2)
      df3(8) = -a*cos(x)/(sin(x)**2*sin(y)*sin(z))-2*a*cos(x)*cos(z)**2/
     1   (sin(x)**2*sin(y)*sin(z)**3)
      df3(9) = -a*cos(y)/(sin(x)*sin(y)**2*sin(z))-2*a*cos(y)*cos(z)**2/
     1   (sin(x)*sin(y)**2*sin(z)**3)
      df3(10) = -5*a*cos(z)/(sin(x)*sin(y)*sin(z)**2)-6*a*cos(z)**3/(sin
     1   (x)*sin(y)*sin(z)**4)
      end 
