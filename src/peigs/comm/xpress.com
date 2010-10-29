*
* $Id$
*
C	This is a common block provided by Express for certain message passing 
C	communication configurations. 
c
      integer nocare
      integer norder
      integer nonode
      integer ihost
      integer ialnod
      integer ialprc
      common/xpress/ nocare,norder,nonode,ihost,ialnod,ialprc
