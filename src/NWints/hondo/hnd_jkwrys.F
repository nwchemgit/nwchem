      subroutine hnd_jkwrys(rwv,abv,numg)
c
c $Id: hnd_jkwrys.F,v 1.1 2004-07-21 16:24:36 bert Exp $
c
      implicit double precision (a-h,o-z)
#include "hnd_rys.fh"
      dimension rwv(2,numg,nroots),abv(5,numg)
c
      do 120 ng=1,numg
      yy=abv(5,ng)
      call hnd_droot
      do 110 nr=1,nroots
      rwv(1,ng,nr)=u9(nr)
      rwv(2,ng,nr)=w9(nr)
  110 continue
  120 continue
      return
      end
