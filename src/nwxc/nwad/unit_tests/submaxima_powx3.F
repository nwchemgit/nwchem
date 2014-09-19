      subroutine submaxima_powx3(x,y,z,a,f,df,df2,df3) 
      implicit none 
      double precision x,y,z,a,f,df(3),df2(3*4/2),df3(3*4*5/6) 
      f = (sin(z)+sin(y)+sin(x)+3.01d+0)**a
      df(1) = a*cos(x)*(sin(z)+sin(y)+sin(x)+3.01d+0)**(a-1)
      df(2) = a*cos(y)*(sin(z)+sin(y)+sin(x)+3.01d+0)**(a-1)
      df(3) = a*cos(z)*(sin(z)+sin(y)+sin(x)+3.01d+0)**(a-1)
      df2(1) = (a-1)*a*cos(x)**2*(sin(z)+sin(y)+sin(x)+3.01d+0)**(a-2)-a
     1   *sin(x)*(sin(z)+sin(y)+sin(x)+3.01d+0)**(a-1)
      df2(2) = (a-1)*a*cos(x)*cos(y)*(sin(z)+sin(y)+sin(x)+3.01d+0)**(a-
     1   2)
      df2(3) = (a-1)*a*cos(y)**2*(sin(z)+sin(y)+sin(x)+3.01d+0)**(a-2)-a
     1   *sin(y)*(sin(z)+sin(y)+sin(x)+3.01d+0)**(a-1)
      df2(4) = (a-1)*a*cos(x)*cos(z)*(sin(z)+sin(y)+sin(x)+3.01d+0)**(a-
     1   2)
      df2(5) = (a-1)*a*cos(y)*cos(z)*(sin(z)+sin(y)+sin(x)+3.01d+0)**(a-
     1   2)
      df2(6) = (a-1)*a*cos(z)**2*(sin(z)+sin(y)+sin(x)+3.01d+0)**(a-2)-a
     1   *sin(z)*(sin(z)+sin(y)+sin(x)+3.01d+0)**(a-1)
      df3(1) = -a*cos(x)*(sin(z)+sin(y)+sin(x)+3.01d+0)**(a-1)-3*(a-1)*a
     1   *cos(x)*sin(x)*(sin(z)+sin(y)+sin(x)+3.01d+0)**(a-2)+(a-2)*(a-1
     2   )*a*cos(x)**3*(sin(z)+sin(y)+sin(x)+3.01d+0)**(a-3)
      df3(2) = (a-2)*(a-1)*a*cos(x)**2*cos(y)*(sin(z)+sin(y)+sin(x)+3.01
     1   d+0)**(a-3)-(a-1)*a*sin(x)*cos(y)*(sin(z)+sin(y)+sin(x)+3.01d+0
     2   )**(a-2)
      df3(3) = (a-2)*(a-1)*a*cos(x)*cos(y)**2*(sin(z)+sin(y)+sin(x)+3.01
     1   d+0)**(a-3)-(a-1)*a*cos(x)*sin(y)*(sin(z)+sin(y)+sin(x)+3.01d+0
     2   )**(a-2)
      df3(4) = -a*cos(y)*(sin(z)+sin(y)+sin(x)+3.01d+0)**(a-1)-3*(a-1)*a
     1   *cos(y)*sin(y)*(sin(z)+sin(y)+sin(x)+3.01d+0)**(a-2)+(a-2)*(a-1
     2   )*a*cos(y)**3*(sin(z)+sin(y)+sin(x)+3.01d+0)**(a-3)
      df3(5) = (a-2)*(a-1)*a*cos(x)**2*cos(z)*(sin(z)+sin(y)+sin(x)+3.01
     1   d+0)**(a-3)-(a-1)*a*sin(x)*cos(z)*(sin(z)+sin(y)+sin(x)+3.01d+0
     2   )**(a-2)
      df3(6) = (a-2)*(a-1)*a*cos(x)*cos(y)*cos(z)*(sin(z)+sin(y)+sin(x)+
     1   3.01d+0)**(a-3)
      df3(7) = (a-2)*(a-1)*a*cos(y)**2*cos(z)*(sin(z)+sin(y)+sin(x)+3.01
     1   d+0)**(a-3)-(a-1)*a*sin(y)*cos(z)*(sin(z)+sin(y)+sin(x)+3.01d+0
     2   )**(a-2)
      df3(8) = (a-2)*(a-1)*a*cos(x)*cos(z)**2*(sin(z)+sin(y)+sin(x)+3.01
     1   d+0)**(a-3)-(a-1)*a*cos(x)*sin(z)*(sin(z)+sin(y)+sin(x)+3.01d+0
     2   )**(a-2)
      df3(9) = (a-2)*(a-1)*a*cos(y)*cos(z)**2*(sin(z)+sin(y)+sin(x)+3.01
     1   d+0)**(a-3)-(a-1)*a*cos(y)*sin(z)*(sin(z)+sin(y)+sin(x)+3.01d+0
     2   )**(a-2)
      df3(10) = -a*cos(z)*(sin(z)+sin(y)+sin(x)+3.01d+0)**(a-1)-3*(a-1)*
     1   a*cos(z)*sin(z)*(sin(z)+sin(y)+sin(x)+3.01d+0)**(a-2)+(a-2)*(a-
     2   1)*a*cos(z)**3*(sin(z)+sin(y)+sin(x)+3.01d+0)**(a-3)
      end 
