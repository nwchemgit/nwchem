c
c $Id: cm_cons.inc,v 1.1 2008-10-01 22:32:26 marat Exp $
c

      integer conlist

      real*8  consdist 

      dimension conlist(mxcons2,3),consdist(mxcons)

      common/cm_cons/conlist,consdist
