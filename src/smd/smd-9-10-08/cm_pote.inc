c
c $Id: cm_pote.inc,v 1.1 2008-10-01 22:32:27 marat Exp $
c

      integer potkey

      real*8  potpar 

      dimension potkey(mxpote),potpar(mxpote,3)

      common/cm_pote/potkey,potpar

