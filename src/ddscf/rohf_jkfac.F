      subroutine rohf_jkfac(jfac,kfac)
*
* $Id$
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
c
      return
      end
