c
c $Id: cm_bond.inc,v 1.1 2008-10-01 22:32:26 marat Exp $
c

      integer bonlist

      real*8  bonddist,bondfrce

      dimension bonlist(mxbond2,3),bonddist(mxbond),bondfrce(mxbond)

      common/cm_bond/bonlist,bonddist,bondfrce
