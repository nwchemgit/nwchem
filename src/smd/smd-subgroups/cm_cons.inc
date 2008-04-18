c
c $Id: cm_cons.inc,v 1.1 2008-04-18 17:48:12 marat Exp $
c

      integer conlist

      real*8  consdist 

      dimension conlist(mxcons2,3),consdist(mxcons)

      common/cm_cons/conlist,consdist
