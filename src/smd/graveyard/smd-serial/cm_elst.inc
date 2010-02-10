c
c $Id: cm_elst.inc,v 1.1 2008-04-18 17:40:30 marat Exp $
c

      integer epoint,elist 

      dimension epoint(mxatms),elist(mxatms*mxelist)

      common/cm_elst/epoint,elist
