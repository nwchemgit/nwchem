C $Id: pot2.F,v 1.2 2001-01-08 15:52:00 windus Exp $
      Subroutine drdy_pot2
     &      (x,v,dx,fmat,hess,xs,amass,natom,n3,n3tm,ndbar)
      Implicit Double Precision (A-H,O-Z)
      Dimension x(n3),dx(n3),fmat(n3tm,n3),hess(ndbar),xs(n3),amass(n3)
C
C potential routine specific for NWChem
      Call nwc_s2 (x,v,dx,fmat,hess,xs,amass,natom,n3,n3tm,ndbar)
C
      return
      end
