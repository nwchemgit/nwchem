c $Id$
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     hondo : 2e- giao memory routine
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine hnd_giao_2emem(maxang)
      implicit none
#include "apiP.fh"
#include "hnd_pointers.fh"
c
      integer maxnum, modtyp, maxang

      maxnum = ((maxang+1)*(maxang+2))/2
      modtyp = maxang+2
      maxxyz = 255 
c     maxxyz = 255 / 3
c
c     ----- size of gijkl =
c                         1   if s  shells
c                        81   if p  shells
c                      1296   if d  shells
c                     10000   if f  shells
c                     50625   if g  shells
c
      ngint  = (maxnum**4)
c
c     x(igijkl) = ( 2-d , 4 centers ) integrals ( 3 sets for -giao- )
c     x(ignkl ) = ( 2-d , 3 centers ) integrals
c     x(ignm  ) = ( 2-d , 2 centers ) integrals
c     x(idij  ) = contraction density for -ij- charge distribution
c     x(idkl  ) = contraction density for -kl- charge distribution
c     x(ib00  ) = -b00-
c     x(ib01  ) = -b01-
c     x(ib10  ) = -b10-
c     x(ic00  ) = -c00-
c     x(id00  ) = -d00-
c     x(if00  ) = -f00-
c     x(iabv  ) = -ab- vector for primitive integrals
c     x(icv   ) = -cv- vector for primitive integrals
c     x(irw   ) = -rw- vector for rys roots and weights
c
      igijkl=1
      ignkl =igijkl+( modtyp**2       * modtyp**2       )*maxxyz*3*3
      ignm  =ignkl +( modtyp**2       *(modtyp+modtyp-1))*maxxyz*3
      ib00  =ignm  +((modtyp+modtyp-1)*(modtyp+modtyp-1))*maxxyz*3
      ib01  =ib00  +((modtyp+modtyp-1)                  )*maxxyz
      ib10  =ib01  +((modtyp+modtyp-1)                  )*maxxyz
      ic00  =ib10  +((modtyp+modtyp-1)                  )*maxxyz
      id00  =ic00  +(  1                                )*maxxyz*3
      if00  =id00  +(  1                                )*maxxyz*3
      iabv  =if00  +(  1                                )*maxxyz
      icv   =iabv  +(  5                                )*maxxyz
      irw   =icv   +( 12                                )*maxxyz
      iscrln=irw   +(  2                                )*maxxyz
c
      return
      end
