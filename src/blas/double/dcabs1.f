      double precision function dcabs1(z)
*
* $Id: dcabs1.f,v 1.2 1997-03-17 21:20:46 d3e129 Exp $
*
      double complex z,zz
      double precision t(2)
      equivalence (zz,t(1))
      zz = z
      dcabs1 = dabs(t(1)) + dabs(t(2))
      return
      end
