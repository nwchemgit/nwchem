      subroutine submaxima_sqrt3(x,y,z,f,df,df2,df3) 
      implicit none 
      double precision x,y,z,f,df(3),df2(3*4/2),df3(3*4*5/6) 
      f = 5.773502691896257d-1*(sin(z)+sin(y)+sin(x))**5.0d-1
      df(1) = 2.8867513459481287d-1*cos(x)/(sin(z)+sin(y)+sin(x))**5.0d-
     1   1
      df(2) = 2.8867513459481287d-1*cos(y)/(sin(z)+sin(y)+sin(x))**5.0d-
     1   1
      df(3) = 2.8867513459481287d-1*cos(z)/(sin(z)+sin(y)+sin(x))**5.0d-
     1   1
      df2(1) = -2.8867513459481287d-1*sin(x)/(sin(z)+sin(y)+sin(x))**5.0
     1   d-1-1.4433756729740643d-1*cos(x)**2/(sin(z)+sin(y)+sin(x))**1.5
     2   d+0
      df2(2) = -1.4433756729740643d-1*cos(x)*cos(y)/(sin(z)+sin(y)+sin(x
     1   ))**1.5d+0
      df2(3) = -2.8867513459481287d-1*sin(y)/(sin(z)+sin(y)+sin(x))**5.0
     1   d-1-1.4433756729740643d-1*cos(y)**2/(sin(z)+sin(y)+sin(x))**1.5
     2   d+0
      df2(4) = -1.4433756729740643d-1*cos(x)*cos(z)/(sin(z)+sin(y)+sin(x
     1   ))**1.5d+0
      df2(5) = -1.4433756729740643d-1*cos(y)*cos(z)/(sin(z)+sin(y)+sin(x
     1   ))**1.5d+0
      df2(6) = -2.8867513459481287d-1*sin(z)/(sin(z)+sin(y)+sin(x))**5.0
     1   d-1-1.4433756729740643d-1*cos(z)**2/(sin(z)+sin(y)+sin(x))**1.5
     2   d+0
      df3(1) = -2.8867513459481287d-1*cos(x)/(sin(z)+sin(y)+sin(x))**5.0
     1   d-1+4.330127018922193d-1*cos(x)*sin(x)/(sin(z)+sin(y)+sin(x))**
     2   1.5d+0+2.1650635094610965d-1*cos(x)**3/(sin(z)+sin(y)+sin(x))**
     3   2.5d+0
      df3(2) = 1.4433756729740643d-1*sin(x)*cos(y)/(sin(z)+sin(y)+sin(x)
     1   )**1.5d+0+2.1650635094610965d-1*cos(x)**2*cos(y)/(sin(z)+sin(y)
     2   +sin(x))**2.5d+0
      df3(3) = 1.4433756729740643d-1*cos(x)*sin(y)/(sin(z)+sin(y)+sin(x)
     1   )**1.5d+0+2.1650635094610965d-1*cos(x)*cos(y)**2/(sin(z)+sin(y)
     2   +sin(x))**2.5d+0
      df3(4) = -2.8867513459481287d-1*cos(y)/(sin(z)+sin(y)+sin(x))**5.0
     1   d-1+4.330127018922193d-1*cos(y)*sin(y)/(sin(z)+sin(y)+sin(x))**
     2   1.5d+0+2.1650635094610965d-1*cos(y)**3/(sin(z)+sin(y)+sin(x))**
     3   2.5d+0
      df3(5) = 1.4433756729740643d-1*sin(x)*cos(z)/(sin(z)+sin(y)+sin(x)
     1   )**1.5d+0+2.1650635094610965d-1*cos(x)**2*cos(z)/(sin(z)+sin(y)
     2   +sin(x))**2.5d+0
      df3(6) = 2.1650635094610965d-1*cos(x)*cos(y)*cos(z)/(sin(z)+sin(y)
     1   +sin(x))**2.5d+0
      df3(7) = 1.4433756729740643d-1*sin(y)*cos(z)/(sin(z)+sin(y)+sin(x)
     1   )**1.5d+0+2.1650635094610965d-1*cos(y)**2*cos(z)/(sin(z)+sin(y)
     2   +sin(x))**2.5d+0
      df3(8) = 1.4433756729740643d-1*cos(x)*sin(z)/(sin(z)+sin(y)+sin(x)
     1   )**1.5d+0+2.1650635094610965d-1*cos(x)*cos(z)**2/(sin(z)+sin(y)
     2   +sin(x))**2.5d+0
      df3(9) = 1.4433756729740643d-1*cos(y)*sin(z)/(sin(z)+sin(y)+sin(x)
     1   )**1.5d+0+2.1650635094610965d-1*cos(y)*cos(z)**2/(sin(z)+sin(y)
     2   +sin(x))**2.5d+0
      df3(10) = -2.8867513459481287d-1*cos(z)/(sin(z)+sin(y)+sin(x))**5.
     1   0d-1+4.330127018922193d-1*cos(z)*sin(z)/(sin(z)+sin(y)+sin(x))*
     2   *1.5d+0+2.1650635094610965d-1*cos(z)**3/(sin(z)+sin(y)+sin(x))*
     3   *2.5d+0
      end 
