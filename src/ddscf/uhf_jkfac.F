      subroutine uhf_jkfac(jfac,kfac)
c
      implicit none
      double precision jfac(4),kfac(4)
c
      jfac(1)=1d0
      jfac(2)=0d0
      jfac(3)=1d0
      jfac(4)=0d0
      kfac(1)=0d0
      kfac(2)=1d0
      kfac(3)=0d0
      kfac(4)=1d0
c
      return
      end
c $Id$
