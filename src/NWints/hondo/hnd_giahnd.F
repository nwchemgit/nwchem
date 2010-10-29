      subroutine hnd_giahnd(axyz,aprim,acoef,a_prim,a_gen,la,
     &                      bxyz,bprim,bcoef,b_prim,b_gen,lb,
     &                      cxyz,cprim,ccoef,c_prim,c_gen,lc,
     &                      dxyz,dprim,dcoef,d_prim,d_gen,ld,
     &                      ieqj,keql,giao_eri,leri,scr,lscr)
c
c $Id$
c
      implicit none
#include "nwc_const.fh"
#include "hnd_tol.fh"
#include "hnd_giao.fh"
#include "hnd_pointers.fh"
#include "hnd_rys.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
      common/testdata/timing(20),irepeat
      double precision timing
      integer irepeat
c
c     Input parameters
c
      integer a_prim, a_gen, la
      integer b_prim, b_gen, lb
      integer c_prim, c_gen, lc
      integer d_prim, d_gen, ld
      integer leri, lscr
      double precision axyz(3), aprim(a_prim), acoef(a_prim,a_gen)
      double precision bxyz(3), bprim(a_prim), bcoef(b_prim,b_gen)
      double precision cxyz(3), cprim(a_prim), ccoef(c_prim,c_gen)
      double precision dxyz(3), dprim(a_prim), dcoef(d_prim,d_gen)
      double precision giao_eri(*)
      double precision scr(lscr)
      logical ieqj,keql
c
c     ----- Zero out giao_eri integral block
c
c     timing(1)=timing(1)-util_wallsec()
      call dcopy(leri*6,0.0d0,0,giao_eri,1)
c     timing(1)=timing(1)+util_wallsec()
c
c     ----- get one-electron charge distributions -----
c
      iieqjj = ieqj
      kkeqll = keql
c     timing(2)=timing(2)-util_wallsec()
      call hnd_oechrg(axyz,aprim,acoef,a_prim,a_gen,la,
     &                bxyz,bprim,bcoef,b_prim,b_gen,lb,
     &                cxyz,cprim,ccoef,c_prim,c_gen,lc,
     &                dxyz,dprim,dcoef,d_prim,d_gen,ld)
c     timing(2)=timing(2)+util_wallsec()
c
c     ----- two-electron integrals -----
c
      nroots=(lit+ljt+lkt+lmt-1)/2
c
c     ----- -giao- (ii,jj//kk,ll) -----
c
c     timing(11)=timing(11)-util_wallsec()
      call hnd_giaspd(giao_eri,scr(igijkl),scr(ignkl),scr(ignm),
     2                scr(ib00),scr(ib01),scr(ib10),scr(ic00),scr(id00),
     3                scr(if00),scr(iabv),scr(icv),
     4                scr(irw))
c     timing(11)=timing(11)+util_wallsec()
c
      return
      end
