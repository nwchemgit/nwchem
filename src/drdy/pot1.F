C $Id: pot1.F,v 1.3 2001-01-19 17:43:13 windus Exp $
      Subroutine drdy_pot1 (x,v,dx,xs,amass,natom,ndim)
      Implicit none
      Integer natom,ndim
      Double Precision x(ndim),v,dx(ndim),xs(ndim),amass(ndim)
C
C Transform to cartesian coordinates
      Call drdy_trans (2,ndim,amass,x,dx)
C potential routine specific for NWChem
      Call nwc_s1 (x,v,dx,xs,amass,natom,ndim)
C  Transform coordinates and gradient to mass-weight cartesians
      Call drdy_trans (1,ndim,amass,x,dx)
C
      return
      end
