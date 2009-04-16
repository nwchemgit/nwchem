*
* $Id: xpress.com,v 1.2 1999-07-28 00:39:06 d3e129 Exp $
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
