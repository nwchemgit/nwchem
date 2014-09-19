      subroutine submaxima_tan3(x,y,z,f,df,df2,df3) 
      implicit none 
      double precision x,y,z,f,df(3),df2(3*4/2),df3(3*4*5/6) 
      double precision sec 
      external         sec 
      f = tan(1.57d+0*sin(x)*sin(y)*sin(z))
      df(1) = 1.57d+0*cos(x)*sin(y)*sin(z)*sec(1.57d+0*sin(x)*sin(y)*sin
     1   (z))**2
      df(2) = 1.57d+0*sin(x)*cos(y)*sin(z)*sec(1.57d+0*sin(x)*sin(y)*sin
     1   (z))**2
      df(3) = 1.57d+0*sin(x)*sin(y)*cos(z)*sec(1.57d+0*sin(x)*sin(y)*sin
     1   (z))**2
      df2(1) = 4.9298000000000003d+0*cos(x)**2*sin(y)**2*sin(z)**2*sec(1
     1   .57d+0*sin(x)*sin(y)*sin(z))**2*tan(1.57d+0*sin(x)*sin(y)*sin(z
     2   ))-1.57d+0*sin(x)*sin(y)*sin(z)*sec(1.57d+0*sin(x)*sin(y)*sin(z
     3   ))**2
      df2(2) = 4.9298000000000003d+0*cos(x)*sin(x)*cos(y)*sin(y)*sin(z)*
     1   *2*sec(1.57d+0*sin(x)*sin(y)*sin(z))**2*tan(1.57d+0*sin(x)*sin(
     2   y)*sin(z))+1.57d+0*cos(x)*cos(y)*sin(z)*sec(1.57d+0*sin(x)*sin(
     3   y)*sin(z))**2
      df2(3) = 4.9298000000000003d+0*sin(x)**2*cos(y)**2*sin(z)**2*sec(1
     1   .57d+0*sin(x)*sin(y)*sin(z))**2*tan(1.57d+0*sin(x)*sin(y)*sin(z
     2   ))-1.57d+0*sin(x)*sin(y)*sin(z)*sec(1.57d+0*sin(x)*sin(y)*sin(z
     3   ))**2
      df2(4) = 4.9298000000000003d+0*cos(x)*sin(x)*sin(y)**2*cos(z)*sin(
     1   z)*sec(1.57d+0*sin(x)*sin(y)*sin(z))**2*tan(1.57d+0*sin(x)*sin(
     2   y)*sin(z))+1.57d+0*cos(x)*sin(y)*cos(z)*sec(1.57d+0*sin(x)*sin(
     3   y)*sin(z))**2
      df2(5) = 4.9298000000000003d+0*sin(x)**2*cos(y)*sin(y)*cos(z)*sin(
     1   z)*sec(1.57d+0*sin(x)*sin(y)*sin(z))**2*tan(1.57d+0*sin(x)*sin(
     2   y)*sin(z))+1.57d+0*sin(x)*cos(y)*cos(z)*sec(1.57d+0*sin(x)*sin(
     3   y)*sin(z))**2
      df2(6) = 4.9298000000000003d+0*sin(x)**2*sin(y)**2*cos(z)**2*sec(1
     1   .57d+0*sin(x)*sin(y)*sin(z))**2*tan(1.57d+0*sin(x)*sin(y)*sin(z
     2   ))-1.57d+0*sin(x)*sin(y)*sin(z)*sec(1.57d+0*sin(x)*sin(y)*sin(z
     3   ))**2
      df3(1) = 1.5479572d+1*cos(x)**3*sin(y)**3*sin(z)**3*sec(1.57d+0*si
     1   n(x)*sin(y)*sin(z))**2*tan(1.57d+0*sin(x)*sin(y)*sin(z))**2-1.4
     2   7894d+1*cos(x)*sin(x)*sin(y)**2*sin(z)**2*sec(1.57d+0*sin(x)*si
     3   n(y)*sin(z))**2*tan(1.57d+0*sin(x)*sin(y)*sin(z))+7.73978600000
     4   0001d+0*cos(x)**3*sin(y)**3*sin(z)**3*sec(1.57d+0*sin(x)*sin(y)
     5   *sin(z))**4-1.57d+0*cos(x)*sin(y)*sin(z)*sec(1.57d+0*sin(x)*sin
     6   (y)*sin(z))**2
      df3(2) = 1.5479572d+1*cos(x)**2*sin(x)*cos(y)*sin(y)**2*sin(z)**3*
     1   sec(1.57d+0*sin(x)*sin(y)*sin(z))**2*tan(1.57d+0*sin(x)*sin(y)*
     2   sin(z))**2-4.9298000000000003d+0*sin(x)**2*cos(y)*sin(y)*sin(z)
     3   **2*sec(1.57d+0*sin(x)*sin(y)*sin(z))**2*tan(1.57d+0*sin(x)*sin
     4   (y)*sin(z))+9.859600000000001d+0*cos(x)**2*cos(y)*sin(y)*sin(z)
     5   **2*sec(1.57d+0*sin(x)*sin(y)*sin(z))**2*tan(1.57d+0*sin(x)*sin
     6   (y)*sin(z))+7.739786000000001d+0*cos(x)**2*sin(x)*cos(y)*sin(y)
     7   **2*sin(z)**3*sec(1.57d+0*sin(x)*sin(y)*sin(z))**4-1.57d+0*sin(
     8   x)*cos(y)*sin(z)*sec(1.57d+0*sin(x)*sin(y)*sin(z))**2
      df3(3) = 1.5479572d+1*cos(x)*sin(x)**2*cos(y)**2*sin(y)*sin(z)**3*
     1   sec(1.57d+0*sin(x)*sin(y)*sin(z))**2*tan(1.57d+0*sin(x)*sin(y)*
     2   sin(z))**2-4.9298000000000003d+0*cos(x)*sin(x)*sin(y)**2*sin(z)
     3   **2*sec(1.57d+0*sin(x)*sin(y)*sin(z))**2*tan(1.57d+0*sin(x)*sin
     4   (y)*sin(z))+9.859600000000001d+0*cos(x)*sin(x)*cos(y)**2*sin(z)
     5   **2*sec(1.57d+0*sin(x)*sin(y)*sin(z))**2*tan(1.57d+0*sin(x)*sin
     6   (y)*sin(z))+7.739786000000001d+0*cos(x)*sin(x)**2*cos(y)**2*sin
     7   (y)*sin(z)**3*sec(1.57d+0*sin(x)*sin(y)*sin(z))**4-1.57d+0*cos(
     8   x)*sin(y)*sin(z)*sec(1.57d+0*sin(x)*sin(y)*sin(z))**2
      df3(4) = 1.5479572d+1*sin(x)**3*cos(y)**3*sin(z)**3*sec(1.57d+0*si
     1   n(x)*sin(y)*sin(z))**2*tan(1.57d+0*sin(x)*sin(y)*sin(z))**2-1.4
     2   7894d+1*sin(x)**2*cos(y)*sin(y)*sin(z)**2*sec(1.57d+0*sin(x)*si
     3   n(y)*sin(z))**2*tan(1.57d+0*sin(x)*sin(y)*sin(z))+7.73978600000
     4   0001d+0*sin(x)**3*cos(y)**3*sin(z)**3*sec(1.57d+0*sin(x)*sin(y)
     5   *sin(z))**4-1.57d+0*sin(x)*cos(y)*sin(z)*sec(1.57d+0*sin(x)*sin
     6   (y)*sin(z))**2
      df3(5) = 1.5479572d+1*cos(x)**2*sin(x)*sin(y)**3*cos(z)*sin(z)**2*
     1   sec(1.57d+0*sin(x)*sin(y)*sin(z))**2*tan(1.57d+0*sin(x)*sin(y)*
     2   sin(z))**2-4.9298000000000003d+0*sin(x)**2*sin(y)**2*cos(z)*sin
     3   (z)*sec(1.57d+0*sin(x)*sin(y)*sin(z))**2*tan(1.57d+0*sin(x)*sin
     4   (y)*sin(z))+9.859600000000001d+0*cos(x)**2*sin(y)**2*cos(z)*sin
     5   (z)*sec(1.57d+0*sin(x)*sin(y)*sin(z))**2*tan(1.57d+0*sin(x)*sin
     6   (y)*sin(z))+7.739786000000001d+0*cos(x)**2*sin(x)*sin(y)**3*cos
     7   (z)*sin(z)**2*sec(1.57d+0*sin(x)*sin(y)*sin(z))**4-1.57d+0*sin(
     8   x)*sin(y)*cos(z)*sec(1.57d+0*sin(x)*sin(y)*sin(z))**2
      df3(6) = 1.5479572d+1*cos(x)*sin(x)**2*cos(y)*sin(y)**2*cos(z)*sin
     1   (z)**2*sec(1.57d+0*sin(x)*sin(y)*sin(z))**2*tan(1.57d+0*sin(x)*
     2   sin(y)*sin(z))**2+1.47894d+1*cos(x)*sin(x)*cos(y)*sin(y)*cos(z)
     3   *sin(z)*sec(1.57d+0*sin(x)*sin(y)*sin(z))**2*tan(1.57d+0*sin(x)
     4   *sin(y)*sin(z))+7.739786000000001d+0*cos(x)*sin(x)**2*cos(y)*si
     5   n(y)**2*cos(z)*sin(z)**2*sec(1.57d+0*sin(x)*sin(y)*sin(z))**4+1
     6   .57d+0*cos(x)*cos(y)*cos(z)*sec(1.57d+0*sin(x)*sin(y)*sin(z))**
     7   2
      df3(7) = 1.5479572d+1*sin(x)**3*cos(y)**2*sin(y)*cos(z)*sin(z)**2*
     1   sec(1.57d+0*sin(x)*sin(y)*sin(z))**2*tan(1.57d+0*sin(x)*sin(y)*
     2   sin(z))**2-4.9298000000000003d+0*sin(x)**2*sin(y)**2*cos(z)*sin
     3   (z)*sec(1.57d+0*sin(x)*sin(y)*sin(z))**2*tan(1.57d+0*sin(x)*sin
     4   (y)*sin(z))+9.859600000000001d+0*sin(x)**2*cos(y)**2*cos(z)*sin
     5   (z)*sec(1.57d+0*sin(x)*sin(y)*sin(z))**2*tan(1.57d+0*sin(x)*sin
     6   (y)*sin(z))+7.739786000000001d+0*sin(x)**3*cos(y)**2*sin(y)*cos
     7   (z)*sin(z)**2*sec(1.57d+0*sin(x)*sin(y)*sin(z))**4-1.57d+0*sin(
     8   x)*sin(y)*cos(z)*sec(1.57d+0*sin(x)*sin(y)*sin(z))**2
      df3(8) = 1.5479572d+1*cos(x)*sin(x)**2*sin(y)**3*cos(z)**2*sin(z)*
     1   sec(1.57d+0*sin(x)*sin(y)*sin(z))**2*tan(1.57d+0*sin(x)*sin(y)*
     2   sin(z))**2-4.9298000000000003d+0*cos(x)*sin(x)*sin(y)**2*sin(z)
     3   **2*sec(1.57d+0*sin(x)*sin(y)*sin(z))**2*tan(1.57d+0*sin(x)*sin
     4   (y)*sin(z))+9.859600000000001d+0*cos(x)*sin(x)*sin(y)**2*cos(z)
     5   **2*sec(1.57d+0*sin(x)*sin(y)*sin(z))**2*tan(1.57d+0*sin(x)*sin
     6   (y)*sin(z))+7.739786000000001d+0*cos(x)*sin(x)**2*sin(y)**3*cos
     7   (z)**2*sin(z)*sec(1.57d+0*sin(x)*sin(y)*sin(z))**4-1.57d+0*cos(
     8   x)*sin(y)*sin(z)*sec(1.57d+0*sin(x)*sin(y)*sin(z))**2
      df3(9) = 1.5479572d+1*sin(x)**3*cos(y)*sin(y)**2*cos(z)**2*sin(z)*
     1   sec(1.57d+0*sin(x)*sin(y)*sin(z))**2*tan(1.57d+0*sin(x)*sin(y)*
     2   sin(z))**2-4.9298000000000003d+0*sin(x)**2*cos(y)*sin(y)*sin(z)
     3   **2*sec(1.57d+0*sin(x)*sin(y)*sin(z))**2*tan(1.57d+0*sin(x)*sin
     4   (y)*sin(z))+9.859600000000001d+0*sin(x)**2*cos(y)*sin(y)*cos(z)
     5   **2*sec(1.57d+0*sin(x)*sin(y)*sin(z))**2*tan(1.57d+0*sin(x)*sin
     6   (y)*sin(z))+7.739786000000001d+0*sin(x)**3*cos(y)*sin(y)**2*cos
     7   (z)**2*sin(z)*sec(1.57d+0*sin(x)*sin(y)*sin(z))**4-1.57d+0*sin(
     8   x)*cos(y)*sin(z)*sec(1.57d+0*sin(x)*sin(y)*sin(z))**2
      df3(10) = 1.5479572d+1*sin(x)**3*sin(y)**3*cos(z)**3*sec(1.57d+0*s
     1   in(x)*sin(y)*sin(z))**2*tan(1.57d+0*sin(x)*sin(y)*sin(z))**2-1.
     2   47894d+1*sin(x)**2*sin(y)**2*cos(z)*sin(z)*sec(1.57d+0*sin(x)*s
     3   in(y)*sin(z))**2*tan(1.57d+0*sin(x)*sin(y)*sin(z))+7.7397860000
     4   00001d+0*sin(x)**3*sin(y)**3*cos(z)**3*sec(1.57d+0*sin(x)*sin(y
     5   )*sin(z))**4-1.57d+0*sin(x)*sin(y)*cos(z)*sec(1.57d+0*sin(x)*si
     6   n(y)*sin(z))**2
      end 
      double precision function sec(x) 
      implicit double precision (a-z) 
      sec=1.d0/cos(x) 
      return 
      end 
