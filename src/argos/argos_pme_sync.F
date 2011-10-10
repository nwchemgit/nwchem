      subroutine argos_pme_sync(number,nodes)
c
      implicit none
c
      integer number
      integer nodes(number)
c
      external mxcombv1,argos_pme_synop
      integer mxcombv1,argos_pme_synop
c
      integer ires
      real*8 datum,buffer
c
      datum=0.0d0
      ires=mxcombv1(datum,argos_pme_synop,8,1,number,nodes,123,buffer)
c
      return
      end
c $Id$
