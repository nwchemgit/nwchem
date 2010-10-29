C $Id$
      Subroutine drdy_pot2 (x,v,dx,fmat,hess,xs,amass,
     * natom,ndim,maxdim,maxbar)
C
      Implicit none
C
#include "nwc_drdyP.fh"
C
      Integer natom,ndim,maxdim,maxbar
      Double Precision x(ndim),v,dx(ndim),fmat(maxdim,ndim),
     * hess(maxbar),xs(ndim),amass(ndim)
C
C Transform to cartesian coordinates
      Call drdy_trans (2,ndim,amass,x,dx)
C potential routine specific for NWChem
      Call nwc_gen (x,v,dx,fmat,hess,xs,natom,ndim,maxdim,maxbar,
     &   DRDY_CODE_HESSIAN)
C  Transform coordinates and gradient to mass-weight cartesians
      Call drdy_trans (1,ndim,amass,x,dx)
C  Transform hessian matrix to mass-weighted coordinates
      Call drdy_tranf (1,maxdim,ndim,amass,fmat)
C
      return
      end
