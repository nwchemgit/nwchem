      function chelp_atom_rad(i)
      implicit none
#include "errquit.fh"
      integer i
c
      double precision radius(107)
      double precision chelp_atom_rad
      character*30 pname

      pname = "chelp_atom_rad"
c
      chelp_atom_rad = 0.0d0
c
      radius(  1)=0.030d0
      radius(  2)=0.122d0
      radius(  3)=0.123d0
      radius(  4)=0.089d0
      radius(  5)=0.088d0
      radius(  6)=0.077d0
      radius(  7)=0.070d0
      radius(  8)=0.066d0
      radius(  9)=0.058d0
      radius( 10)=0.160d0
      radius( 11)=0.140d0
      radius( 12)=0.136d0
      radius( 13)=0.125d0
      radius( 14)=0.117d0
      radius( 15)=0.110d0
      radius( 16)=0.104d0
      radius( 17)=0.099d0
      radius( 18)=0.191d0
      radius( 19)=0.203d0
      radius( 20)=0.174d0
      radius( 21)=0.144d0
      radius( 22)=0.132d0
      radius( 23)=0.122d0
      radius( 24)=0.119d0
      radius( 25)=0.117d0
      radius( 26)=0.117d0
      radius( 27)=0.116d0
      radius( 28)=0.115d0
      radius( 29)=0.117d0
      radius( 30)=0.125d0
      radius( 31)=0.125d0
      radius( 32)=0.122d0
      radius( 33)=0.121d0
      radius( 34)=0.117d0
      radius( 35)=0.114d0
      radius( 36)=0.198d0
      radius( 37)=0.222d0
      radius( 38)=0.192d0
      radius( 39)=0.162d0
      radius( 40)=0.145d0
      radius( 41)=0.134d0
      radius( 42)=0.129d0
      radius( 43)=0.127d0
      radius( 44)=0.124d0
      radius( 45)=0.125d0
      radius( 46)=0.128d0
      radius( 47)=0.134d0
      radius( 48)=0.141d0
      radius( 49)=0.150d0
      radius( 50)=0.140d0
      radius( 51)=0.141d0
      radius( 52)=0.137d0
      radius( 53)=0.133d0
      radius( 54)=0.209d0
      radius( 55)=0.235d0
      radius( 56)=0.198d0
      radius( 57)=0.169d0
      radius( 58)=0.165d0
      radius( 59)=0.165d0
      radius( 60)=0.164d0
      radius( 61)=0.165d0
      radius( 62)=0.166d0
      radius( 63)=0.165d0
      radius( 64)=0.161d0
      radius( 65)=0.159d0
      radius( 66)=0.159d0
      radius( 67)=0.158d0
      radius( 68)=0.157d0
      radius( 69)=0.156d0
      radius( 70)=0.156d0
      radius( 71)=0.156d0
      radius( 72)=0.144d0
      radius( 73)=0.134d0
      radius( 74)=0.130d0
      radius( 75)=0.128d0
      radius( 76)=0.126d0
      radius( 77)=0.126d0
      radius( 78)=0.129d0
      radius( 79)=0.134d0
      radius( 80)=0.144d0
      radius( 81)=0.155d0
      radius( 82)=0.154d0
      radius( 83)=0.152d0
      radius( 84)=0.153d0
      radius( 85)=0.150d0
      radius( 86)=0.220d0
      radius( 87)=0.324d0
      radius( 88)=0.268d0
      radius( 89)=0.225d0
      radius( 90)=0.216d0
      radius( 91)=0.193d0
      radius( 92)=0.166d0
      radius( 93)=0.157d0
      radius( 94)=0.181d0
      radius( 95)=0.221d0
      radius( 96)=0.143d0
      radius( 97)=0.142d0
      radius( 98)=0.140d0
      radius( 99)=0.139d0
      radius(100)=0.138d0
      radius(101)=0.137d0
      radius(102)=0.136d0
      radius(103)=0.134d0
      radius(104)=0.130d0
      radius(105)=0.130d0
      radius(106)=0.130d0
      radius(107)=0.130d0
 
      if(i.lt.107) then
        chelp_atom_rad = radius(i)   
      else
        call errquit("out of bounds in "//pname,0,0)
      end if

      end function
c $Id$
