      subroutine submaxima_tanh3(x,y,z,f,df,df2,df3) 
      implicit none 
      double precision x,y,z,f,df(3),df2(3*4/2),df3(3*4*5/6) 
      double precision sech 
      external         sech 
      f = tanh(6.666666666666666d-1*(sin(z)+sin(y)+sin(x)))
      df(1) = 6.666666666666666d-1*cos(x)*sech(6.666666666666666d-1*(sin
     1   (z)+sin(y)+sin(x)))**2
      df(2) = 6.666666666666666d-1*cos(y)*sech(6.666666666666666d-1*(sin
     1   (z)+sin(y)+sin(x)))**2
      df(3) = 6.666666666666666d-1*cos(z)*sech(6.666666666666666d-1*(sin
     1   (z)+sin(y)+sin(x)))**2
      df2(1) = -8.888888888888888d-1*cos(x)**2*sech(6.666666666666666d-1
     1   *(sin(z)+sin(y)+sin(x)))**2*tanh(6.666666666666666d-1*(sin(z)+s
     2   in(y)+sin(x)))-6.666666666666666d-1*sin(x)*sech(6.6666666666666
     3   66d-1*(sin(z)+sin(y)+sin(x)))**2
      df2(2) = -8.888888888888888d-1*cos(x)*cos(y)*sech(6.66666666666666
     1   6d-1*(sin(z)+sin(y)+sin(x)))**2*tanh(6.666666666666666d-1*(sin(
     2   z)+sin(y)+sin(x)))
      df2(3) = -8.888888888888888d-1*cos(y)**2*sech(6.666666666666666d-1
     1   *(sin(z)+sin(y)+sin(x)))**2*tanh(6.666666666666666d-1*(sin(z)+s
     2   in(y)+sin(x)))-6.666666666666666d-1*sin(y)*sech(6.6666666666666
     3   66d-1*(sin(z)+sin(y)+sin(x)))**2
      df2(4) = -8.888888888888888d-1*cos(x)*cos(z)*sech(6.66666666666666
     1   6d-1*(sin(z)+sin(y)+sin(x)))**2*tanh(6.666666666666666d-1*(sin(
     2   z)+sin(y)+sin(x)))
      df2(5) = -8.888888888888888d-1*cos(y)*cos(z)*sech(6.66666666666666
     1   6d-1*(sin(z)+sin(y)+sin(x)))**2*tanh(6.666666666666666d-1*(sin(
     2   z)+sin(y)+sin(x)))
      df2(6) = -8.888888888888888d-1*cos(z)**2*sech(6.666666666666666d-1
     1   *(sin(z)+sin(y)+sin(x)))**2*tanh(6.666666666666666d-1*(sin(z)+s
     2   in(y)+sin(x)))-6.666666666666666d-1*sin(z)*sech(6.6666666666666
     3   66d-1*(sin(z)+sin(y)+sin(x)))**2
      df3(1) = 1.1851851851851851d+0*cos(x)**3*sech(6.666666666666666d-1
     1   *(sin(z)+sin(y)+sin(x)))**2*tanh(6.666666666666666d-1*(sin(z)+s
     2   in(y)+sin(x)))**2+2.6666666666666666d+0*cos(x)*sin(x)*sech(6.66
     3   6666666666666d-1*(sin(z)+sin(y)+sin(x)))**2*tanh(6.666666666666
     4   666d-1*(sin(z)+sin(y)+sin(x)))-5.925925925925926d-1*cos(x)**3*s
     5   ech(6.666666666666666d-1*(sin(z)+sin(y)+sin(x)))**4-6.666666666
     6   666666d-1*cos(x)*sech(6.666666666666666d-1*(sin(z)+sin(y)+sin(x
     7   )))**2
      df3(2) = 1.1851851851851851d+0*cos(x)**2*cos(y)*sech(6.66666666666
     1   6666d-1*(sin(z)+sin(y)+sin(x)))**2*tanh(6.666666666666666d-1*(s
     2   in(z)+sin(y)+sin(x)))**2+8.888888888888888d-1*sin(x)*cos(y)*sec
     3   h(6.666666666666666d-1*(sin(z)+sin(y)+sin(x)))**2*tanh(6.666666
     4   666666666d-1*(sin(z)+sin(y)+sin(x)))-5.925925925925926d-1*cos(x
     5   )**2*cos(y)*sech(6.666666666666666d-1*(sin(z)+sin(y)+sin(x)))**
     6   4
      df3(3) = 1.1851851851851851d+0*cos(x)*cos(y)**2*sech(6.66666666666
     1   6666d-1*(sin(z)+sin(y)+sin(x)))**2*tanh(6.666666666666666d-1*(s
     2   in(z)+sin(y)+sin(x)))**2+8.888888888888888d-1*cos(x)*sin(y)*sec
     3   h(6.666666666666666d-1*(sin(z)+sin(y)+sin(x)))**2*tanh(6.666666
     4   666666666d-1*(sin(z)+sin(y)+sin(x)))-5.925925925925926d-1*cos(x
     5   )*cos(y)**2*sech(6.666666666666666d-1*(sin(z)+sin(y)+sin(x)))**
     6   4
      df3(4) = 1.1851851851851851d+0*cos(y)**3*sech(6.666666666666666d-1
     1   *(sin(z)+sin(y)+sin(x)))**2*tanh(6.666666666666666d-1*(sin(z)+s
     2   in(y)+sin(x)))**2+2.6666666666666666d+0*cos(y)*sin(y)*sech(6.66
     3   6666666666666d-1*(sin(z)+sin(y)+sin(x)))**2*tanh(6.666666666666
     4   666d-1*(sin(z)+sin(y)+sin(x)))-5.925925925925926d-1*cos(y)**3*s
     5   ech(6.666666666666666d-1*(sin(z)+sin(y)+sin(x)))**4-6.666666666
     6   666666d-1*cos(y)*sech(6.666666666666666d-1*(sin(z)+sin(y)+sin(x
     7   )))**2
      df3(5) = 1.1851851851851851d+0*cos(x)**2*cos(z)*sech(6.66666666666
     1   6666d-1*(sin(z)+sin(y)+sin(x)))**2*tanh(6.666666666666666d-1*(s
     2   in(z)+sin(y)+sin(x)))**2+8.888888888888888d-1*sin(x)*cos(z)*sec
     3   h(6.666666666666666d-1*(sin(z)+sin(y)+sin(x)))**2*tanh(6.666666
     4   666666666d-1*(sin(z)+sin(y)+sin(x)))-5.925925925925926d-1*cos(x
     5   )**2*cos(z)*sech(6.666666666666666d-1*(sin(z)+sin(y)+sin(x)))**
     6   4
      df3(6) = 1.1851851851851851d+0*cos(x)*cos(y)*cos(z)*sech(6.6666666
     1   66666666d-1*(sin(z)+sin(y)+sin(x)))**2*tanh(6.666666666666666d-
     2   1*(sin(z)+sin(y)+sin(x)))**2-5.925925925925926d-1*cos(x)*cos(y)
     3   *cos(z)*sech(6.666666666666666d-1*(sin(z)+sin(y)+sin(x)))**4
      df3(7) = 1.1851851851851851d+0*cos(y)**2*cos(z)*sech(6.66666666666
     1   6666d-1*(sin(z)+sin(y)+sin(x)))**2*tanh(6.666666666666666d-1*(s
     2   in(z)+sin(y)+sin(x)))**2+8.888888888888888d-1*sin(y)*cos(z)*sec
     3   h(6.666666666666666d-1*(sin(z)+sin(y)+sin(x)))**2*tanh(6.666666
     4   666666666d-1*(sin(z)+sin(y)+sin(x)))-5.925925925925926d-1*cos(y
     5   )**2*cos(z)*sech(6.666666666666666d-1*(sin(z)+sin(y)+sin(x)))**
     6   4
      df3(8) = 1.1851851851851851d+0*cos(x)*cos(z)**2*sech(6.66666666666
     1   6666d-1*(sin(z)+sin(y)+sin(x)))**2*tanh(6.666666666666666d-1*(s
     2   in(z)+sin(y)+sin(x)))**2+8.888888888888888d-1*cos(x)*sin(z)*sec
     3   h(6.666666666666666d-1*(sin(z)+sin(y)+sin(x)))**2*tanh(6.666666
     4   666666666d-1*(sin(z)+sin(y)+sin(x)))-5.925925925925926d-1*cos(x
     5   )*cos(z)**2*sech(6.666666666666666d-1*(sin(z)+sin(y)+sin(x)))**
     6   4
      df3(9) = 1.1851851851851851d+0*cos(y)*cos(z)**2*sech(6.66666666666
     1   6666d-1*(sin(z)+sin(y)+sin(x)))**2*tanh(6.666666666666666d-1*(s
     2   in(z)+sin(y)+sin(x)))**2+8.888888888888888d-1*cos(y)*sin(z)*sec
     3   h(6.666666666666666d-1*(sin(z)+sin(y)+sin(x)))**2*tanh(6.666666
     4   666666666d-1*(sin(z)+sin(y)+sin(x)))-5.925925925925926d-1*cos(y
     5   )*cos(z)**2*sech(6.666666666666666d-1*(sin(z)+sin(y)+sin(x)))**
     6   4
      df3(10) = 1.1851851851851851d+0*cos(z)**3*sech(6.666666666666666d-
     1   1*(sin(z)+sin(y)+sin(x)))**2*tanh(6.666666666666666d-1*(sin(z)+
     2   sin(y)+sin(x)))**2+2.6666666666666666d+0*cos(z)*sin(z)*sech(6.6
     3   66666666666666d-1*(sin(z)+sin(y)+sin(x)))**2*tanh(6.66666666666
     4   6666d-1*(sin(z)+sin(y)+sin(x)))-5.925925925925926d-1*cos(z)**3*
     5   sech(6.666666666666666d-1*(sin(z)+sin(y)+sin(x)))**4-6.66666666
     6   6666666d-1*cos(z)*sech(6.666666666666666d-1*(sin(z)+sin(y)+sin(
     7   x)))**2
      end 
      double precision function sech(x) 
      implicit double precision (a-z) 
      sech=2.d0/(exp(x)+exp(-x)) 
      return 
      end 
