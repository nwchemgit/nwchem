C$Id
************************************************************************
*                                                                      *
      subroutine ecp_so_ang (l,angint,soint,n_ang)
*                                                                      *
*   Perform sum over m for spin-orbit integrals                        *
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
        if (m .eq. 1) wa = wa+wa
*
*     x spin-orbit component
*
        call daxpy (n_ang,wa,angint(-m,1),ll,soint(m-1,1,3),ll)
        call daxpy (n_ang,-wa,angint(m-1,1),ll,soint(-m,1,3),ll)
*
*     y spin-orbit component
*
        call daxpy (n_ang,wa,angint(m-1,1),ll,soint(m,1,2),ll)
        call daxpy (n_ang,-wa,angint(m,1),ll,soint(m-1,1,2),ll)
        call daxpy (n_ang,wa,angint(-m,1),ll,soint(1-m,1,2),ll)
        call daxpy (n_ang,-wa,angint(1-m,1),ll,soint(-m,1,2),ll)
      end do
*
      return
      end
