      integer bonlist

      real*8  bonddist,bondfrce

      dimension bonlist(mxbond2,3),bonddist(mxbond),bondfrce(mxbond)

      common/cm_bond/bonlist,bonddist,bondfrce
