      Subroutine drdy_pot1 (x,v,dx,xs,amass,natom,n3)
      Implicit Double Precision (A-H,O-Z)
      Dimension x(n3),dx(n3),xs(n3),amass(n3)
C
C potential routine specific for Gaussian
      Call nwc_s1 (x,v,dx,xs,amass,natom,n3)
C
      return
      end
