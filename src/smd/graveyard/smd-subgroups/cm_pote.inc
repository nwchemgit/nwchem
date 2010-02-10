c
c $Id: cm_pote.inc,v 1.1 2008-04-18 17:48:12 marat Exp $
c

      integer potkey

      real*8  potpar 

      dimension potkey(mxpote),potpar(mxpote,3)

      common/cm_pote/potkey,potpar

