C$Id: ecp_so_ang.F,v 1.4 1998-04-03 01:40:39 mg201 Exp $
************************************************************************
*                                                                      *
      subroutine ecp_so_ang (l,angint,soint,n_ang)
*                                                                      *
*   Perform sum over m for spin-orbit integrals. This includes the 2/l *
*   factor in the s.l formalism - see the last piece of code. To use   *
*   sigma instead of s, change from two to one.                        *
*                                                                      *
*   Argument (status) - description                                    *
*                                                                      *
*   l (inp) - angular momentum of projector                            *
*   angint - angular integrals for LH projector                        *
*   soint - angular integrals with spin-orbit term added               *
*   n_ang - number of angular integrals                                *
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
      wa = l
*      wa = two/wa
      wa = one/wa
      call dscal (ll*n_ang*3,wa,soint,1)
*
      return
      end
