      subroutine submaxima_pow3(x,y,z,f,df,df2,df3) 
      implicit none 
      double precision x,y,z,f,df(3),df2(3*4/2),df3(3*4*5/6) 
      f = (sin(y)+sin(x)+2.0099999999999998d+0)**(sin(z)+sin(y))
      df(1) = cos(x)*(sin(y)+sin(x)+2.0099999999999998d+0)**(sin(z)+sin(
     1   y)-1)*(sin(z)+sin(y))
      df(2) = (sin(y)+sin(x)+2.0099999999999998d+0)**(sin(z)+sin(y))*(co
     1   s(y)*(sin(z)+sin(y))/(sin(y)+sin(x)+2.0099999999999998d+0)+cos(
     2   y)*log(sin(y)+sin(x)+2.0099999999999998d+0))
      df(3) = (sin(y)+sin(x)+2.0099999999999998d+0)**(sin(z)+sin(y))*log
     1   (sin(y)+sin(x)+2.0099999999999998d+0)*cos(z)
      df2(1) = cos(x)**2*(sin(y)+sin(x)+2.0099999999999998d+0)**(sin(z)+
     1   sin(y)-2)*(sin(z)+sin(y)-1)*(sin(z)+sin(y))-sin(x)*(sin(y)+sin(
     2   x)+2.0099999999999998d+0)**(sin(z)+sin(y)-1)*(sin(z)+sin(y))
      df2(2) = cos(x)*(sin(y)+sin(x)+2.0099999999999998d+0)**(sin(z)+sin
     1   (y)-1)*(sin(z)+sin(y))*(cos(y)*(sin(z)+sin(y)-1)/(sin(y)+sin(x)
     2   +2.0099999999999998d+0)+cos(y)*log(sin(y)+sin(x)+2.009999999999
     3   9998d+0))+cos(x)*cos(y)*(sin(y)+sin(x)+2.0099999999999998d+0)**
     4   (sin(z)+sin(y)-1)
      df2(3) = (sin(y)+sin(x)+2.0099999999999998d+0)**(sin(z)+sin(y))*(c
     1   os(y)*(sin(z)+sin(y))/(sin(y)+sin(x)+2.0099999999999998d+0)+cos
     2   (y)*log(sin(y)+sin(x)+2.0099999999999998d+0))**2+(sin(y)+sin(x)
     3   +2.0099999999999998d+0)**(sin(z)+sin(y))*(-sin(y)*(sin(z)+sin(y
     4   ))/(sin(y)+sin(x)+2.0099999999999998d+0)-cos(y)**2*(sin(z)+sin(
     5   y))/(sin(y)+sin(x)+2.0099999999999998d+0)**2-sin(y)*log(sin(y)+
     6   sin(x)+2.0099999999999998d+0)+2*cos(y)**2/(sin(y)+sin(x)+2.0099
     7   999999999998d+0))
      df2(4) = cos(x)*(sin(y)+sin(x)+2.0099999999999998d+0)**(sin(z)+sin
     1   (y)-1)*log(sin(y)+sin(x)+2.0099999999999998d+0)*cos(z)*(sin(z)+
     2   sin(y))+cos(x)*(sin(y)+sin(x)+2.0099999999999998d+0)**(sin(z)+s
     3   in(y)-1)*cos(z)
      df2(5) = (sin(y)+sin(x)+2.0099999999999998d+0)**(sin(z)+sin(y))*lo
     1   g(sin(y)+sin(x)+2.0099999999999998d+0)*cos(z)*(cos(y)*(sin(z)+s
     2   in(y))/(sin(y)+sin(x)+2.0099999999999998d+0)+cos(y)*log(sin(y)+
     3   sin(x)+2.0099999999999998d+0))+cos(y)*(sin(y)+sin(x)+2.00999999
     4   99999998d+0)**(sin(z)+sin(y)-1)*cos(z)
      df2(6) = (sin(y)+sin(x)+2.0099999999999998d+0)**(sin(z)+sin(y))*lo
     1   g(sin(y)+sin(x)+2.0099999999999998d+0)**2*cos(z)**2-(sin(y)+sin
     2   (x)+2.0099999999999998d+0)**(sin(z)+sin(y))*log(sin(y)+sin(x)+2
     3   .0099999999999998d+0)*sin(z)
      df3(1) = cos(x)**3*(sin(y)+sin(x)+2.0099999999999998d+0)**(sin(z)+
     1   sin(y)-3)*(sin(z)+sin(y)-2)*(sin(z)+sin(y)-1)*(sin(z)+sin(y))-3
     2   *cos(x)*sin(x)*(sin(y)+sin(x)+2.0099999999999998d+0)**(sin(z)+s
     3   in(y)-2)*(sin(z)+sin(y)-1)*(sin(z)+sin(y))-cos(x)*(sin(y)+sin(x
     4   )+2.0099999999999998d+0)**(sin(z)+sin(y)-1)*(sin(z)+sin(y))
      df3(2) = -sin(x)*(sin(y)+sin(x)+2.0099999999999998d+0)**(sin(z)+si
     1   n(y)-1)*(sin(z)+sin(y))*(cos(y)*(sin(z)+sin(y)-1)/(sin(y)+sin(x
     2   )+2.0099999999999998d+0)+cos(y)*log(sin(y)+sin(x)+2.00999999999
     3   99998d+0))+cos(x)**2*(sin(y)+sin(x)+2.0099999999999998d+0)**(si
     4   n(z)+sin(y)-2)*(sin(z)+sin(y)-1)*(sin(z)+sin(y))*(cos(y)*(sin(z
     5   )+sin(y)-2)/(sin(y)+sin(x)+2.0099999999999998d+0)+cos(y)*log(si
     6   n(y)+sin(x)+2.0099999999999998d+0))+cos(x)**2*cos(y)*(sin(y)+si
     7   n(x)+2.0099999999999998d+0)**(sin(z)+sin(y)-2)*(sin(z)+sin(y))+
     8   cos(x)**2*cos(y)*(sin(y)+sin(x)+2.0099999999999998d+0)**(sin(z)
     9   +sin(y)-2)*(sin(z)+sin(y)-1)-sin(x)*cos(y)*(sin(y)+sin(x)+2.009
     :   9999999999998d+0)**(sin(z)+sin(y)-1)
      df3(3) = cos(x)*(sin(y)+sin(x)+2.0099999999999998d+0)**(sin(z)+sin
     1   (y)-1)*(sin(z)+sin(y))*(cos(y)*(sin(z)+sin(y)-1)/(sin(y)+sin(x)
     2   +2.0099999999999998d+0)+cos(y)*log(sin(y)+sin(x)+2.009999999999
     3   9998d+0))**2+cos(x)*(sin(y)+sin(x)+2.0099999999999998d+0)**(sin
     4   (z)+sin(y)-1)*(sin(z)+sin(y))*(-sin(y)*(sin(z)+sin(y)-1)/(sin(y
     5   )+sin(x)+2.0099999999999998d+0)-cos(y)**2*(sin(z)+sin(y)-1)/(si
     6   n(y)+sin(x)+2.0099999999999998d+0)**2-sin(y)*log(sin(y)+sin(x)+
     7   2.0099999999999998d+0)+2*cos(y)**2/(sin(y)+sin(x)+2.00999999999
     8   99998d+0))+2*cos(x)*cos(y)*(sin(y)+sin(x)+2.0099999999999998d+0
     9   )**(sin(z)+sin(y)-1)*(cos(y)*(sin(z)+sin(y)-1)/(sin(y)+sin(x)+2
     :   .0099999999999998d+0)+cos(y)*log(sin(y)+sin(x)+2.00999999999999
     ;   98d+0))-cos(x)*sin(y)*(sin(y)+sin(x)+2.0099999999999998d+0)**(s
     <   in(z)+sin(y)-1)
      df3(4) = (sin(y)+sin(x)+2.0099999999999998d+0)**(sin(z)+sin(y))*(c
     1   os(y)*(sin(z)+sin(y))/(sin(y)+sin(x)+2.0099999999999998d+0)+cos
     2   (y)*log(sin(y)+sin(x)+2.0099999999999998d+0))**3+3*(sin(y)+sin(
     3   x)+2.0099999999999998d+0)**(sin(z)+sin(y))*(cos(y)*(sin(z)+sin(
     4   y))/(sin(y)+sin(x)+2.0099999999999998d+0)+cos(y)*log(sin(y)+sin
     5   (x)+2.0099999999999998d+0))*(-sin(y)*(sin(z)+sin(y))/(sin(y)+si
     6   n(x)+2.0099999999999998d+0)-cos(y)**2*(sin(z)+sin(y))/(sin(y)+s
     7   in(x)+2.0099999999999998d+0)**2-sin(y)*log(sin(y)+sin(x)+2.0099
     8   999999999998d+0)+2*cos(y)**2/(sin(y)+sin(x)+2.0099999999999998d
     9   +0))+(sin(y)+sin(x)+2.0099999999999998d+0)**(sin(z)+sin(y))*(-c
     :   os(y)*(sin(z)+sin(y))/(sin(y)+sin(x)+2.0099999999999998d+0)+3*c
     ;   os(y)*sin(y)*(sin(z)+sin(y))/(sin(y)+sin(x)+2.0099999999999998d
     <   +0)**2+2*cos(y)**3*(sin(z)+sin(y))/(sin(y)+sin(x)+2.00999999999
     =   99998d+0)**3-cos(y)*log(sin(y)+sin(x)+2.0099999999999998d+0)-6*
     >   cos(y)*sin(y)/(sin(y)+sin(x)+2.0099999999999998d+0)-3*cos(y)**3
     ?   /(sin(y)+sin(x)+2.0099999999999998d+0)**2)
      df3(5) = cos(x)**2*(sin(y)+sin(x)+2.0099999999999998d+0)**(sin(z)+
     1   sin(y)-2)*log(sin(y)+sin(x)+2.0099999999999998d+0)*cos(z)*(sin(
     2   z)+sin(y)-1)*(sin(z)+sin(y))+cos(x)**2*(sin(y)+sin(x)+2.0099999
     3   999999998d+0)**(sin(z)+sin(y)-2)*cos(z)*(sin(z)+sin(y))-sin(x)*
     4   (sin(y)+sin(x)+2.0099999999999998d+0)**(sin(z)+sin(y)-1)*log(si
     5   n(y)+sin(x)+2.0099999999999998d+0)*cos(z)*(sin(z)+sin(y))+cos(x
     6   )**2*(sin(y)+sin(x)+2.0099999999999998d+0)**(sin(z)+sin(y)-2)*c
     7   os(z)*(sin(z)+sin(y)-1)-sin(x)*(sin(y)+sin(x)+2.009999999999999
     8   8d+0)**(sin(z)+sin(y)-1)*cos(z)
      df3(6) = cos(x)*(sin(y)+sin(x)+2.0099999999999998d+0)**(sin(z)+sin
     1   (y)-1)*log(sin(y)+sin(x)+2.0099999999999998d+0)*cos(z)*(sin(z)+
     2   sin(y))*(cos(y)*(sin(z)+sin(y)-1)/(sin(y)+sin(x)+2.009999999999
     3   9998d+0)+cos(y)*log(sin(y)+sin(x)+2.0099999999999998d+0))+cos(x
     4   )*(sin(y)+sin(x)+2.0099999999999998d+0)**(sin(z)+sin(y)-1)*cos(
     5   z)*(cos(y)*(sin(z)+sin(y)-1)/(sin(y)+sin(x)+2.0099999999999998d
     6   +0)+cos(y)*log(sin(y)+sin(x)+2.0099999999999998d+0))+cos(x)*cos
     7   (y)*(sin(y)+sin(x)+2.0099999999999998d+0)**(sin(z)+sin(y)-2)*co
     8   s(z)*(sin(z)+sin(y))+cos(x)*cos(y)*(sin(y)+sin(x)+2.00999999999
     9   99998d+0)**(sin(z)+sin(y)-1)*log(sin(y)+sin(x)+2.00999999999999
     :   98d+0)*cos(z)
      df3(7) = (sin(y)+sin(x)+2.0099999999999998d+0)**(sin(z)+sin(y))*lo
     1   g(sin(y)+sin(x)+2.0099999999999998d+0)*cos(z)*(cos(y)*(sin(z)+s
     2   in(y))/(sin(y)+sin(x)+2.0099999999999998d+0)+cos(y)*log(sin(y)+
     3   sin(x)+2.0099999999999998d+0))**2+(sin(y)+sin(x)+2.009999999999
     4   9998d+0)**(sin(z)+sin(y))*log(sin(y)+sin(x)+2.0099999999999998d
     5   +0)*cos(z)*(-sin(y)*(sin(z)+sin(y))/(sin(y)+sin(x)+2.0099999999
     6   999998d+0)-cos(y)**2*(sin(z)+sin(y))/(sin(y)+sin(x)+2.009999999
     7   9999998d+0)**2-sin(y)*log(sin(y)+sin(x)+2.0099999999999998d+0)+
     8   2*cos(y)**2/(sin(y)+sin(x)+2.0099999999999998d+0))+2*cos(y)*(si
     9   n(y)+sin(x)+2.0099999999999998d+0)**(sin(z)+sin(y)-1)*cos(z)*(c
     :   os(y)*(sin(z)+sin(y))/(sin(y)+sin(x)+2.0099999999999998d+0)+cos
     ;   (y)*log(sin(y)+sin(x)+2.0099999999999998d+0))+(sin(y)+sin(x)+2.
     <   0099999999999998d+0)**(sin(z)+sin(y))*(-sin(y)*cos(z)/(sin(y)+s
     =   in(x)+2.0099999999999998d+0)-cos(y)**2*cos(z)/(sin(y)+sin(x)+2.
     >   0099999999999998d+0)**2)
      df3(8) = -cos(x)*(sin(y)+sin(x)+2.0099999999999998d+0)**(sin(z)+si
     1   n(y)-1)*log(sin(y)+sin(x)+2.0099999999999998d+0)*sin(z)*(sin(z)
     2   +sin(y))+cos(x)*(sin(y)+sin(x)+2.0099999999999998d+0)**(sin(z)+
     3   sin(y)-1)*log(sin(y)+sin(x)+2.0099999999999998d+0)**2*cos(z)**2
     4   *(sin(z)+sin(y))-cos(x)*(sin(y)+sin(x)+2.0099999999999998d+0)**
     5   (sin(z)+sin(y)-1)*sin(z)+2*cos(x)*(sin(y)+sin(x)+2.009999999999
     6   9998d+0)**(sin(z)+sin(y)-1)*log(sin(y)+sin(x)+2.009999999999999
     7   8d+0)*cos(z)**2
      df3(9) = -(sin(y)+sin(x)+2.0099999999999998d+0)**(sin(z)+sin(y))*l
     1   og(sin(y)+sin(x)+2.0099999999999998d+0)*sin(z)*(cos(y)*(sin(z)+
     2   sin(y))/(sin(y)+sin(x)+2.0099999999999998d+0)+cos(y)*log(sin(y)
     3   +sin(x)+2.0099999999999998d+0))+(sin(y)+sin(x)+2.00999999999999
     4   98d+0)**(sin(z)+sin(y))*log(sin(y)+sin(x)+2.0099999999999998d+0
     5   )**2*cos(z)**2*(cos(y)*(sin(z)+sin(y))/(sin(y)+sin(x)+2.0099999
     6   999999998d+0)+cos(y)*log(sin(y)+sin(x)+2.0099999999999998d+0))-
     7   cos(y)*(sin(y)+sin(x)+2.0099999999999998d+0)**(sin(z)+sin(y)-1)
     8   *sin(z)+2*cos(y)*(sin(y)+sin(x)+2.0099999999999998d+0)**(sin(z)
     9   +sin(y)-1)*log(sin(y)+sin(x)+2.0099999999999998d+0)*cos(z)**2
      df3(10) = -3*(sin(y)+sin(x)+2.0099999999999998d+0)**(sin(z)+sin(y)
     1   )*log(sin(y)+sin(x)+2.0099999999999998d+0)**2*cos(z)*sin(z)+(si
     2   n(y)+sin(x)+2.0099999999999998d+0)**(sin(z)+sin(y))*log(sin(y)+
     3   sin(x)+2.0099999999999998d+0)**3*cos(z)**3-(sin(y)+sin(x)+2.009
     4   9999999999998d+0)**(sin(z)+sin(y))*log(sin(y)+sin(x)+2.00999999
     5   99999998d+0)*cos(z)
      end 
