C $Id: potsp.F,v 1.1 2001-01-19 17:43:13 windus Exp $
      Subroutine drdy_potsp (x,vsp,dx,amass,natom,ndim)
C
      Implicit None
      Integer natom,ndim
      Double Precision x(ndim),vsp,dx(ndim),amass(ndim)
C
C Transform to cartesian coordinates
      Call drdy_trans (2,ndim,amass,x,dx)
C potential routine specific for Gaussian
C for the moment don't do single point calculation
C     Call gsssp (x,vsp,natom,ndim)
C  Transform coordinates and gradient to mass-weight cartesians
      Call drdy_trans (1,ndim,amass,x,dx)
      return
      end
