      subroutine rohf_jkfac(jfac,kfac)
*
* $Id: rohf_jkfac.F,v 1.1 2003-01-27 23:17:10 edo Exp $
*
      implicit none
      double precision jfac(3),kfac(3)
c
      jfac(1)=1d0
      jfac(2)=1d0
      jfac(3)=0d0
      kfac(1)=-.5d0
      kfac(2)=0d0
      kfac(3)=-.5d0
      return
      end
      subroutine uhf_jkfac(jfac,kfac)
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
      return
      end
