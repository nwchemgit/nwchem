      Subroutine grid_setspac_params
c
c$Id: grid_setspac_params.F,v 1.1 1999-04-23 16:30:36 d3h449 Exp $
c
      implicit none
c      
#include "cdft.fh"
#include "util.fh"
c
      double precision eps
      parameter (eps = 1.d-20)
      logical oprint
c
      integer iti, itj, chij
      double precision uij, aij
c
c     Set print options.
c
      oprint = util_print('quadrature', print_high)
c
c     Compute Becke s parameter for center pairs.
c      
      if (.not.ldelley) then 
         call dfill(ntypes*ntypes,0.d0,aratio,1)
         do 310 iti = 1, ntypes
            if (bsrad_atom_type(iti).le.EPS) then
               aratio(iti,1) = -1234.d0
               goto 310
            endif
            do 300 itj = 1, iti-1
               if (bsrad_atom_type(itj).le.EPS) then
                  aratio(iti,itj) = -1234.d0
                  goto 300
               endif
c            
c               chij = sqrt(bsrad_atom_type(iti)/bsrad_atom_type(itj))
               chij = bsrad_atom_type(iti)/bsrad_atom_type(itj)
c
c              uij = (bsrad_atom_type(iti - bsrad_atom_type(itj))/(bsrad_atom_type(iti) + bsrad_atom_type(itj))
c
               uij = (chij-1.d0)/(chij+1.d0)
               aij = uij/(uij*uij - 1.D0)
c               write(*,*)' bsrad_atom_type(iti), bsrad_atom_type(itj), chij, uij, aij = ',
c     &                     bsrad_atom_type(iti), bsrad_atom_type(itj), chij, uij, aij
               if (aij.lt.-0.5d0)then
                  aij = -0.5d0
               elseif (aij.gt.0.5d0)then
                  aij = 0.5d0
               endif
c
               aratio(iti,itj) = aij
               aratio(itj,iti) = -aij
c
  300       continue
  310    continue
c
c        temp writes
c
         write(6,*)' aratio '
         call output(aratio, 1, ntypes, 1, ntypes, 
     &               dft_ntags_bsmx, dft_ntags_bsmx, 1)
         write(6,*)' iatype_pt_chg ', 
     &             (iatype_pt_chg(iti),iti=1,ncenters)
         write(6,*)' iatype ', (iatype(iti),iti=1,ncenters)
      endif
      return
      end
