C$Id$
************************************************************************
*                                                                      *
      subroutine ecp_so_ang (l,angint,soint,n_ang)
*                                                                      *
*   Perform sum over m for spin-orbit integrals. This includes the 2/l *
*   factor in the s.l formalism - see the last piece of code. To use   *
*   sigma instead of s, change from two to one.                        *
*   Note that the sqrt((1+delta_m,0)(1+delta_m',0)) term is included.  *
*                                                                      *
*   Argument (status) - description                                    *
*                                                                      *
*   l (inp) - angular momentum of projector                            *
*   angint (inp) - angular integrals for LH projector                  *
*   soint (out) - angular integrals with spin-orbit term added         *
*   n_ang (inp) - number of angular integrals                          *
*                                                                      *
*   Written by K. G. Dyall                                             *
*                                                                      *
************************************************************************
      implicit none
#include "ecp_consts.fh"
      integer l,ll,m,n_ang
      double precision angint,soint,wa
*
      dimension angint(-l:l,n_ang),soint(-l:l,n_ang,3)
*
      ll = 2*l+1
      call dfill (ll*n_ang*3,zero,soint,1)
      if (l .eq. 0) return
      do m = 1,l
*
*     z spin-orbit component
*
        wa = m
        call daxpy (n_ang,wa,angint(-m,1),ll,soint(m,1,1),ll)
        call daxpy (n_ang,-wa,angint(m,1),ll,soint(-m,1,1),ll)
*
        wa = (l+m)*(l-m+1)
        wa = sqrt(wa)
        if (m .ne. 1) wa = wa*half
*
*     x,y spin-orbit components
*
        call daxpy (n_ang,-wa,angint(m-1,1),ll,soint(m,1,2),ll)
        call daxpy (n_ang,wa,angint(m,1),ll,soint(m-1,1,2),ll)
        call daxpy (n_ang,-wa,angint(-m,1),ll,soint(m-1,1,3),ll)
        call daxpy (n_ang,wa,angint(m-1,1),ll,soint(-m,1,3),ll)
*
        if (m .gt. 1) then
          call daxpy (n_ang,wa,angint(-m,1),ll,soint(1-m,1,2),ll)
          call daxpy (n_ang,-wa,angint(1-m,1),ll,soint(-m,1,2),ll)
          call daxpy (n_ang,-wa,angint(1-m,1),ll,soint(m,1,3),ll)
          call daxpy (n_ang,wa,angint(m,1),ll,soint(1-m,1,3),ll)
        end if
      end do
*
      return
      end
