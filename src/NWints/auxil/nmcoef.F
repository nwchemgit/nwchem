      Subroutine nmcoef(z,cc,l,np,norm)
c $Id: nmcoef.F,v 1.1 2000-02-04 20:20:11 mg201 Exp $

      Implicit none
c passed
      Integer l,np,norm
      Double precision z(np),cc(np)
c local
      Integer l_lim,n,n1,n2
      Double precision z2,z4  !  sums of exponents
      Double precision f      !  double factorial
      Double precision S,sum  !  overlap
      Double precision PI     !  pi
      Parameter (PI=3.1415926535898D0)
#include "basdeclsP.fh"
c
c Generally, the normalization (norm=1) of the Gaussian primitives and 
c contractions is based on the value of the integral
c
c   /   2i           2                      / PI \ 1/2  /  1 \ i
c   |  x   exp( - p x  ) dx  =  (2i - 1)!! | ---- |    | ---- |    .
c   /                                       \ p  /      \ 2p /
c
c Note, however, that the double factorial is not included in the 
c normalization. Thus, the normalization depends only on the total
c angular  momentum (l = i+j+k) and not on the individual angular 
c momentum indices (i,j,k).
c
c
c The factorial is added in for only the total angular momentum not the 
c individual components.
c
c******************************************************************************

      if( norm.eq.BasNorm_UN )then
 
c Do not normalize the GTFs.

        return

      else if( norm.eq.BasNorm_STD )then

c Normalize the GTFs so that the integral over all space of its
c square is unity as detailed above.

        if (np .eq. 1) cc(1) = 1.0d0 ! RJH
        l_lim = 2*l - 1
        f = 1.0d00
        do n = l_lim,2,-2
          f = f*dble(n)
        end do
        do n = 1,np
          cc(n) = cc(n)*((2.D0*z(n)/PI)**0.75D0)*(4.D0*z(n))**(0.5D0*l)
          cc(n) = cc(n) / sqrt(f)
        end do

c Normalize the contractions.

        if (np.gt.1) then
          sum = 0.D0
          do 20 n1 = 1,np
            do 21 n2 = 1,np
              S = ( (PI/(z(n1) + z(n2)))**1.5D0 )/
     &            ( (2.D0*(z(n1) + z(n2)))**l )
              sum = sum + cc(n1)*cc(n2)*S
   21       continue
   20     continue
c...
c add in (2l-1)!!
          f = 1.0d00
          l_lim = 2*l - 1
          do 22 n  = l_lim,2,-2
            f = f * dble(n)
   22     continue
          sum = sum * f
c..
          f = 1/sqrt(sum)
          do 30 n = 1,np
            cc(n) = cc(n)*f
   30     continue
        endif

      else if( norm.eq.BasNorm_2c )then

c Normalize the GTFs by requiring (a|a) = 1. That is, its 2-ctr ERI must 
c be equal to unity. This is deMon's convention for the CD basis set.

       if( np.ne.1 .or. ( l.lt.0 .or. l.gt.2 ) )then

        write(*,*)
        write(*,*) 'NMCOEF:  For (a|a) = 1 type of normalization,'
        write(*,*) '         this routine has not been generalized '
        write(*,*) '         to handle contracted GTFs or those with'
        write(*,*) '         angular momentum greater than D.'
        stop

c This can be generalized by computing 2-ctr ERIs. However, care must 
c be used to select a consistent value in the case when the values are 
c not all equal for each member of the given shell.

       end if

       if( l.eq.0 )then

        f = sqrt(2.D0)*(PI/z(1))**2.5D0

        cc(1) = 1/sqrt(f)

       else if( l.eq.1 )then

        f = (sqrt(2.D0)*(PI/z(1))**2.5D0)/(12.D0*z(1))

        cc(1) = 1/sqrt(f)

       else if( l.eq.2 )then

        f = (sqrt(2.D0)*(PI/z(1))**2.5D0)/(80.D0*z(1)**2)

        cc(1) = 1/sqrt(f)

       end if

      else if( norm.eq.BasNorm_rel )then

c Normalize the GTFs so that the integral over all space of its
c square is unity as detailed above, but do not normalize the 
c contractions.

        l_lim = 2*l - 1
        f = 1.0d00
        do n = l_lim,2,-2
          f = f*dble(n)
        end do
        do n = 1,np
          z2 = z(n)+z(n)
          z4 = z2+z2
          cc(n) = cc(n)*sqrt((sqrt(z2/PI)**3)*(z4**l)/f)
        end do

      end if

      end
