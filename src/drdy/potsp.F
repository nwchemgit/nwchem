C $Id$
      Subroutine drdy_potsp (x,vsp,dx,amass,natom,ndim)
C
      Implicit None
C
#include "nwc_drdyP.fh"
C
      Integer natom,ndim
      Double Precision x(ndim),vsp,dx(ndim),amass(ndim)
      Double Precision dummy
C
C Transform to cartesian coordinates
      Call drdy_trans (2,ndim,amass,x,dx)
C do single point calculation (is fourth dummy OK?)
      Call nwc_gen (x,vsp,dx,dummy,dummy,dummy,natom,ndim,1,1,
     &   DRDY_CODE_SPENERGY)
C  Transform coordinates and gradient to mass-weight cartesians
      Call drdy_trans (1,ndim,amass,x,dx)
      return
      end
