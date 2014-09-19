      subroutine submaxima_erf3(x,y,z,f,df,df2,df3) 
      implicit none 
      double precision x,y,z,f,df(3),df2(3*4/2),df3(3*4*5/6) 
      f = erf(6.666666666666666d-1*(sin(z)+sin(y)+sin(x)))
      df(1) = 7.52252778063675d-1*cos(x)*exp(-4.444444444444444d-1*(sin(
     1   z)+sin(y)+sin(x))**2)
      df(2) = 7.52252778063675d-1*cos(y)*exp(-4.444444444444444d-1*(sin(
     1   z)+sin(y)+sin(x))**2)
      df(3) = 7.52252778063675d-1*exp(-4.444444444444444d-1*(sin(z)+sin(
     1   y)+sin(x))**2)*cos(z)
      df2(1) = -6.686691360565999d-1*cos(x)**2*exp(-4.444444444444444d-1
     1   *(sin(z)+sin(y)+sin(x))**2)*(sin(z)+sin(y)+sin(x))-7.5225277806
     2   3675d-1*sin(x)*exp(-4.444444444444444d-1*(sin(z)+sin(y)+sin(x))
     3   **2)
      df2(2) = -6.686691360565999d-1*cos(x)*cos(y)*exp(-4.44444444444444
     1   4d-1*(sin(z)+sin(y)+sin(x))**2)*(sin(z)+sin(y)+sin(x))
      df2(3) = -6.686691360565999d-1*cos(y)**2*exp(-4.444444444444444d-1
     1   *(sin(z)+sin(y)+sin(x))**2)*(sin(z)+sin(y)+sin(x))-7.5225277806
     2   3675d-1*sin(y)*exp(-4.444444444444444d-1*(sin(z)+sin(y)+sin(x))
     3   **2)
      df2(4) = -6.686691360565999d-1*cos(x)*exp(-4.444444444444444d-1*(s
     1   in(z)+sin(y)+sin(x))**2)*cos(z)*(sin(z)+sin(y)+sin(x))
      df2(5) = -6.686691360565999d-1*cos(y)*exp(-4.444444444444444d-1*(s
     1   in(z)+sin(y)+sin(x))**2)*cos(z)*(sin(z)+sin(y)+sin(x))
      df2(6) = -6.686691360565999d-1*exp(-4.444444444444444d-1*(sin(z)+s
     1   in(y)+sin(x))**2)*cos(z)**2*(sin(z)+sin(y)+sin(x))-7.5225277806
     2   3675d-1*exp(-4.444444444444444d-1*(sin(z)+sin(y)+sin(x))**2)*si
     3   n(z)
      df3(1) = 5.943725653836444d-1*cos(x)**3*exp(-4.444444444444444d-1*
     1   (sin(z)+sin(y)+sin(x))**2)*(sin(z)+sin(y)+sin(x))**2+2.00600740
     2   81698d+0*cos(x)*sin(x)*exp(-4.444444444444444d-1*(sin(z)+sin(y)
     3   +sin(x))**2)*(sin(z)+sin(y)+sin(x))-6.686691360565999d-1*cos(x)
     4   **3*exp(-4.444444444444444d-1*(sin(z)+sin(y)+sin(x))**2)-7.5225
     5   2778063675d-1*cos(x)*exp(-4.444444444444444d-1*(sin(z)+sin(y)+s
     6   in(x))**2)
      df3(2) = 5.943725653836444d-1*cos(x)**2*cos(y)*exp(-4.444444444444
     1   444d-1*(sin(z)+sin(y)+sin(x))**2)*(sin(z)+sin(y)+sin(x))**2+6.6
     2   86691360565999d-1*sin(x)*cos(y)*exp(-4.444444444444444d-1*(sin(
     3   z)+sin(y)+sin(x))**2)*(sin(z)+sin(y)+sin(x))-6.686691360565999d
     4   -1*cos(x)**2*cos(y)*exp(-4.444444444444444d-1*(sin(z)+sin(y)+si
     5   n(x))**2)
      df3(3) = 5.943725653836444d-1*cos(x)*cos(y)**2*exp(-4.444444444444
     1   444d-1*(sin(z)+sin(y)+sin(x))**2)*(sin(z)+sin(y)+sin(x))**2+6.6
     2   86691360565999d-1*cos(x)*sin(y)*exp(-4.444444444444444d-1*(sin(
     3   z)+sin(y)+sin(x))**2)*(sin(z)+sin(y)+sin(x))-6.686691360565999d
     4   -1*cos(x)*cos(y)**2*exp(-4.444444444444444d-1*(sin(z)+sin(y)+si
     5   n(x))**2)
      df3(4) = 5.943725653836444d-1*cos(y)**3*exp(-4.444444444444444d-1*
     1   (sin(z)+sin(y)+sin(x))**2)*(sin(z)+sin(y)+sin(x))**2+2.00600740
     2   81698d+0*cos(y)*sin(y)*exp(-4.444444444444444d-1*(sin(z)+sin(y)
     3   +sin(x))**2)*(sin(z)+sin(y)+sin(x))-6.686691360565999d-1*cos(y)
     4   **3*exp(-4.444444444444444d-1*(sin(z)+sin(y)+sin(x))**2)-7.5225
     5   2778063675d-1*cos(y)*exp(-4.444444444444444d-1*(sin(z)+sin(y)+s
     6   in(x))**2)
      df3(5) = 5.943725653836444d-1*cos(x)**2*exp(-4.444444444444444d-1*
     1   (sin(z)+sin(y)+sin(x))**2)*cos(z)*(sin(z)+sin(y)+sin(x))**2+6.6
     2   86691360565999d-1*sin(x)*exp(-4.444444444444444d-1*(sin(z)+sin(
     3   y)+sin(x))**2)*cos(z)*(sin(z)+sin(y)+sin(x))-6.686691360565999d
     4   -1*cos(x)**2*exp(-4.444444444444444d-1*(sin(z)+sin(y)+sin(x))**
     5   2)*cos(z)
      df3(6) = 5.943725653836444d-1*cos(x)*cos(y)*exp(-4.444444444444444
     1   d-1*(sin(z)+sin(y)+sin(x))**2)*cos(z)*(sin(z)+sin(y)+sin(x))**2
     2   -6.686691360565999d-1*cos(x)*cos(y)*exp(-4.444444444444444d-1*(
     3   sin(z)+sin(y)+sin(x))**2)*cos(z)
      df3(7) = 5.943725653836444d-1*cos(y)**2*exp(-4.444444444444444d-1*
     1   (sin(z)+sin(y)+sin(x))**2)*cos(z)*(sin(z)+sin(y)+sin(x))**2+6.6
     2   86691360565999d-1*sin(y)*exp(-4.444444444444444d-1*(sin(z)+sin(
     3   y)+sin(x))**2)*cos(z)*(sin(z)+sin(y)+sin(x))-6.686691360565999d
     4   -1*cos(y)**2*exp(-4.444444444444444d-1*(sin(z)+sin(y)+sin(x))**
     5   2)*cos(z)
      df3(8) = 5.943725653836444d-1*cos(x)*exp(-4.444444444444444d-1*(si
     1   n(z)+sin(y)+sin(x))**2)*cos(z)**2*(sin(z)+sin(y)+sin(x))**2+6.6
     2   86691360565999d-1*cos(x)*exp(-4.444444444444444d-1*(sin(z)+sin(
     3   y)+sin(x))**2)*sin(z)*(sin(z)+sin(y)+sin(x))-6.686691360565999d
     4   -1*cos(x)*exp(-4.444444444444444d-1*(sin(z)+sin(y)+sin(x))**2)*
     5   cos(z)**2
      df3(9) = 5.943725653836444d-1*cos(y)*exp(-4.444444444444444d-1*(si
     1   n(z)+sin(y)+sin(x))**2)*cos(z)**2*(sin(z)+sin(y)+sin(x))**2+6.6
     2   86691360565999d-1*cos(y)*exp(-4.444444444444444d-1*(sin(z)+sin(
     3   y)+sin(x))**2)*sin(z)*(sin(z)+sin(y)+sin(x))-6.686691360565999d
     4   -1*cos(y)*exp(-4.444444444444444d-1*(sin(z)+sin(y)+sin(x))**2)*
     5   cos(z)**2
      df3(10) = 5.943725653836444d-1*exp(-4.444444444444444d-1*(sin(z)+s
     1   in(y)+sin(x))**2)*cos(z)**3*(sin(z)+sin(y)+sin(x))**2+2.0060074
     2   081698d+0*exp(-4.444444444444444d-1*(sin(z)+sin(y)+sin(x))**2)*
     3   cos(z)*sin(z)*(sin(z)+sin(y)+sin(x))-6.686691360565999d-1*exp(-
     4   4.444444444444444d-1*(sin(z)+sin(y)+sin(x))**2)*cos(z)**3-7.522
     5   52778063675d-1*exp(-4.444444444444444d-1*(sin(z)+sin(y)+sin(x))
     6   **2)*cos(z)
      end 
