      subroutine submaxima_powy3(x,y,z,a,f,df,df2,df3) 
      implicit none 
      double precision x,y,z,a,f,df(3),df2(3*4/2),df3(3*4*5/6) 
      f = a**(sin(z)+sin(y)+sin(x))
      df(1) = a**(sin(z)+sin(y)+sin(x))*log(a)*cos(x)
      df(2) = a**(sin(z)+sin(y)+sin(x))*log(a)*cos(y)
      df(3) = a**(sin(z)+sin(y)+sin(x))*log(a)*cos(z)
      df2(1) = a**(sin(z)+sin(y)+sin(x))*log(a)**2*cos(x)**2-a**(sin(z)+
     1   sin(y)+sin(x))*log(a)*sin(x)
      df2(2) = a**(sin(z)+sin(y)+sin(x))*log(a)**2*cos(x)*cos(y)
      df2(3) = a**(sin(z)+sin(y)+sin(x))*log(a)**2*cos(y)**2-a**(sin(z)+
     1   sin(y)+sin(x))*log(a)*sin(y)
      df2(4) = a**(sin(z)+sin(y)+sin(x))*log(a)**2*cos(x)*cos(z)
      df2(5) = a**(sin(z)+sin(y)+sin(x))*log(a)**2*cos(y)*cos(z)
      df2(6) = a**(sin(z)+sin(y)+sin(x))*log(a)**2*cos(z)**2-a**(sin(z)+
     1   sin(y)+sin(x))*log(a)*sin(z)
      df3(1) = -3*a**(sin(z)+sin(y)+sin(x))*log(a)**2*cos(x)*sin(x)+a**(
     1   sin(z)+sin(y)+sin(x))*log(a)**3*cos(x)**3-a**(sin(z)+sin(y)+sin
     2   (x))*log(a)*cos(x)
      df3(2) = a**(sin(z)+sin(y)+sin(x))*log(a)**3*cos(x)**2*cos(y)-a**(
     1   sin(z)+sin(y)+sin(x))*log(a)**2*sin(x)*cos(y)
      df3(3) = a**(sin(z)+sin(y)+sin(x))*log(a)**3*cos(x)*cos(y)**2-a**(
     1   sin(z)+sin(y)+sin(x))*log(a)**2*cos(x)*sin(y)
      df3(4) = -3*a**(sin(z)+sin(y)+sin(x))*log(a)**2*cos(y)*sin(y)+a**(
     1   sin(z)+sin(y)+sin(x))*log(a)**3*cos(y)**3-a**(sin(z)+sin(y)+sin
     2   (x))*log(a)*cos(y)
      df3(5) = a**(sin(z)+sin(y)+sin(x))*log(a)**3*cos(x)**2*cos(z)-a**(
     1   sin(z)+sin(y)+sin(x))*log(a)**2*sin(x)*cos(z)
      df3(6) = a**(sin(z)+sin(y)+sin(x))*log(a)**3*cos(x)*cos(y)*cos(z)
      df3(7) = a**(sin(z)+sin(y)+sin(x))*log(a)**3*cos(y)**2*cos(z)-a**(
     1   sin(z)+sin(y)+sin(x))*log(a)**2*sin(y)*cos(z)
      df3(8) = a**(sin(z)+sin(y)+sin(x))*log(a)**3*cos(x)*cos(z)**2-a**(
     1   sin(z)+sin(y)+sin(x))*log(a)**2*cos(x)*sin(z)
      df3(9) = a**(sin(z)+sin(y)+sin(x))*log(a)**3*cos(y)*cos(z)**2-a**(
     1   sin(z)+sin(y)+sin(x))*log(a)**2*cos(y)*sin(z)
      df3(10) = -3*a**(sin(z)+sin(y)+sin(x))*log(a)**2*cos(z)*sin(z)+a**
     1   (sin(z)+sin(y)+sin(x))*log(a)**3*cos(z)**3-a**(sin(z)+sin(y)+si
     2   n(x))*log(a)*cos(z)
      end 
