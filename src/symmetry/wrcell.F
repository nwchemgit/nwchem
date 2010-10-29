c**********************************************************************
c
c   reads input deck for the cell parameters
c
c**********************************************************************
      subroutine wrcell(geom,cdist,cang,itype)
C$Id$
      implicit double precision(a-h,o-z)
#include "errquit.fh"
      integer geom
      dimension cdist(3),cang(3)
      logical geom_get_user_units
      external geom_get_user_units
      character*12 units
c
      if (.not. geom_get_user_units(geom, units))
     $     call errquit('wrcell: user units?',0, GEOM_ERR)
c     

      if(itype.eq.3) then
         write(*,13) units
         write(*,14) (cdist(i), i=1,3)
         write(*,15) (cang(i), i=1,3)
      elseif(itype.eq.2) then
         write(*,13) units
         write(*,16) (cdist(i), i=1,2)
         write(*,17) cang(3)
      endif
c
 13   format(12x,'LATTICE PARAMETERS IN ',A12,' AND DEGREES:',/)
 14   format(12x,'AXIS:',7x,'A= ',f8.4,6x,'B= ',f8.4,7x,'C= ',f8.4)
 15   format(12x,'ANGLES: ALPHA= ',f8.4,3x,'BETA= ',f8.4,3x,'gamma= '
     $     ,f8.4)
 16   format(15x,'AXIS:',7x,'A= ',f8.4,6x,'B= ',f8.4)
 17   format(15x,'ANGLES: GAMMA= ',f8.4)
      return
      end
