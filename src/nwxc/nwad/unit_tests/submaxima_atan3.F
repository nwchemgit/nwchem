      subroutine submaxima_atan3(x,y,z,f,df,df2,df3) 
      implicit none 
      double precision x,y,z,f,df(3),df2(3*4/2),df3(3*4*5/6) 
      f = atan(2.0d+0*(sin(z)+sin(y)+sin(x)))
      df(1) = 2.0d+0*cos(x)/(4.0d+0*(sin(z)+sin(y)+sin(x))**2+1)
      df(2) = 2.0d+0*cos(y)/(4.0d+0*(sin(z)+sin(y)+sin(x))**2+1)
      df(3) = 2.0d+0*cos(z)/(4.0d+0*(sin(z)+sin(y)+sin(x))**2+1)
      df2(1) = -2.0d+0*sin(x)/(4.0d+0*(sin(z)+sin(y)+sin(x))**2+1)-1.6d+
     1   1*cos(x)**2*(sin(z)+sin(y)+sin(x))/(4.0d+0*(sin(z)+sin(y)+sin(x
     2   ))**2+1)**2
      df2(2) = -1.6d+1*cos(x)*cos(y)*(sin(z)+sin(y)+sin(x))/(4.0d+0*(sin
     1   (z)+sin(y)+sin(x))**2+1)**2
      df2(3) = -2.0d+0*sin(y)/(4.0d+0*(sin(z)+sin(y)+sin(x))**2+1)-1.6d+
     1   1*cos(y)**2*(sin(z)+sin(y)+sin(x))/(4.0d+0*(sin(z)+sin(y)+sin(x
     2   ))**2+1)**2
      df2(4) = -1.6d+1*cos(x)*cos(z)*(sin(z)+sin(y)+sin(x))/(4.0d+0*(sin
     1   (z)+sin(y)+sin(x))**2+1)**2
      df2(5) = -1.6d+1*cos(y)*cos(z)*(sin(z)+sin(y)+sin(x))/(4.0d+0*(sin
     1   (z)+sin(y)+sin(x))**2+1)**2
      df2(6) = -2.0d+0*sin(z)/(4.0d+0*(sin(z)+sin(y)+sin(x))**2+1)-1.6d+
     1   1*cos(z)**2*(sin(z)+sin(y)+sin(x))/(4.0d+0*(sin(z)+sin(y)+sin(x
     2   ))**2+1)**2
      df3(1) = -2.0d+0*cos(x)/(4.0d+0*(sin(z)+sin(y)+sin(x))**2+1)+4.8d+
     1   1*cos(x)*sin(x)*(sin(z)+sin(y)+sin(x))/(4.0d+0*(sin(z)+sin(y)+s
     2   in(x))**2+1)**2-1.6d+1*cos(x)**3/(4.0d+0*(sin(z)+sin(y)+sin(x))
     3   **2+1)**2+2.56d+2*cos(x)**3*(sin(z)+sin(y)+sin(x))**2/(4.0d+0*(
     4   sin(z)+sin(y)+sin(x))**2+1)**3
      df3(2) = 1.6d+1*sin(x)*cos(y)*(sin(z)+sin(y)+sin(x))/(4.0d+0*(sin(
     1   z)+sin(y)+sin(x))**2+1)**2-1.6d+1*cos(x)**2*cos(y)/(4.0d+0*(sin
     2   (z)+sin(y)+sin(x))**2+1)**2+2.56d+2*cos(x)**2*cos(y)*(sin(z)+si
     3   n(y)+sin(x))**2/(4.0d+0*(sin(z)+sin(y)+sin(x))**2+1)**3
      df3(3) = 1.6d+1*cos(x)*sin(y)*(sin(z)+sin(y)+sin(x))/(4.0d+0*(sin(
     1   z)+sin(y)+sin(x))**2+1)**2-1.6d+1*cos(x)*cos(y)**2/(4.0d+0*(sin
     2   (z)+sin(y)+sin(x))**2+1)**2+2.56d+2*cos(x)*cos(y)**2*(sin(z)+si
     3   n(y)+sin(x))**2/(4.0d+0*(sin(z)+sin(y)+sin(x))**2+1)**3
      df3(4) = -2.0d+0*cos(y)/(4.0d+0*(sin(z)+sin(y)+sin(x))**2+1)+4.8d+
     1   1*cos(y)*sin(y)*(sin(z)+sin(y)+sin(x))/(4.0d+0*(sin(z)+sin(y)+s
     2   in(x))**2+1)**2-1.6d+1*cos(y)**3/(4.0d+0*(sin(z)+sin(y)+sin(x))
     3   **2+1)**2+2.56d+2*cos(y)**3*(sin(z)+sin(y)+sin(x))**2/(4.0d+0*(
     4   sin(z)+sin(y)+sin(x))**2+1)**3
      df3(5) = 1.6d+1*sin(x)*cos(z)*(sin(z)+sin(y)+sin(x))/(4.0d+0*(sin(
     1   z)+sin(y)+sin(x))**2+1)**2-1.6d+1*cos(x)**2*cos(z)/(4.0d+0*(sin
     2   (z)+sin(y)+sin(x))**2+1)**2+2.56d+2*cos(x)**2*cos(z)*(sin(z)+si
     3   n(y)+sin(x))**2/(4.0d+0*(sin(z)+sin(y)+sin(x))**2+1)**3
      df3(6) = 2.56d+2*cos(x)*cos(y)*cos(z)*(sin(z)+sin(y)+sin(x))**2/(4
     1   .0d+0*(sin(z)+sin(y)+sin(x))**2+1)**3-1.6d+1*cos(x)*cos(y)*cos(
     2   z)/(4.0d+0*(sin(z)+sin(y)+sin(x))**2+1)**2
      df3(7) = 1.6d+1*sin(y)*cos(z)*(sin(z)+sin(y)+sin(x))/(4.0d+0*(sin(
     1   z)+sin(y)+sin(x))**2+1)**2-1.6d+1*cos(y)**2*cos(z)/(4.0d+0*(sin
     2   (z)+sin(y)+sin(x))**2+1)**2+2.56d+2*cos(y)**2*cos(z)*(sin(z)+si
     3   n(y)+sin(x))**2/(4.0d+0*(sin(z)+sin(y)+sin(x))**2+1)**3
      df3(8) = 1.6d+1*cos(x)*sin(z)*(sin(z)+sin(y)+sin(x))/(4.0d+0*(sin(
     1   z)+sin(y)+sin(x))**2+1)**2-1.6d+1*cos(x)*cos(z)**2/(4.0d+0*(sin
     2   (z)+sin(y)+sin(x))**2+1)**2+2.56d+2*cos(x)*cos(z)**2*(sin(z)+si
     3   n(y)+sin(x))**2/(4.0d+0*(sin(z)+sin(y)+sin(x))**2+1)**3
      df3(9) = 1.6d+1*cos(y)*sin(z)*(sin(z)+sin(y)+sin(x))/(4.0d+0*(sin(
     1   z)+sin(y)+sin(x))**2+1)**2-1.6d+1*cos(y)*cos(z)**2/(4.0d+0*(sin
     2   (z)+sin(y)+sin(x))**2+1)**2+2.56d+2*cos(y)*cos(z)**2*(sin(z)+si
     3   n(y)+sin(x))**2/(4.0d+0*(sin(z)+sin(y)+sin(x))**2+1)**3
      df3(10) = -2.0d+0*cos(z)/(4.0d+0*(sin(z)+sin(y)+sin(x))**2+1)+4.8d
     1   +1*cos(z)*sin(z)*(sin(z)+sin(y)+sin(x))/(4.0d+0*(sin(z)+sin(y)+
     2   sin(x))**2+1)**2-1.6d+1*cos(z)**3/(4.0d+0*(sin(z)+sin(y)+sin(x)
     3   )**2+1)**2+2.56d+2*cos(z)**3*(sin(z)+sin(y)+sin(x))**2/(4.0d+0*
     4   (sin(z)+sin(y)+sin(x))**2+1)**3
      end 
