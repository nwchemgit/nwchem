      subroutine esp_esp(ene,occ,dns,pot,x,q,grid,val)
c
c $Id: esp_esp.F,v 1.7 1998-08-18 16:53:53 d3j191 Exp $
c
      implicit none
c
#include "esp_params.fh"
#include "esp_common.fh"
#include "global.fh"
#include "bas.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "util.fh"
c
      logical int_normalize
      external int_normalize
c
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
c
      real*8 ene(nbf),occ(nbf),dns(nbf,nbf),pot(nbf,nbf)
      real*8 x(3,natoms),q(natoms)
      real*8 grid(3,mxgrid),val(mxgrid)
c
      integer l_bt,ibtptr,l_st,istptr
c
      integer i,j,ish,jsh,ix,jx,ifr,max1e,mem1,nshell
      integer ilo,ihi,jlo,jhi
      real*8 dist,vt
      integer imin
      real*8 dmin
c
c     calculate density matrix, occupation and energies
c
      call esp_denmat(dns,occ,ene)
c
c     initialize integrals
c
      if(.not.int_normalize(ibasis))
     + call errquit('esp: int_normalize failed',me)
      call int_init(irtdb,1,ibasis)
c
c     get number of shells
c
      if(.not.bas_numcont(ibasis,nshell))
     + call errquit('esp: bas_numcont failed',me)
c
c     allocate memory for the integrals
c
      call int_mem_1e(max1e,mem1)
c
      if(.not.ma_push_get(mt_dbl,max1e,'bt',l_bt,ibtptr))
     + call errquit('esp: ma_push_get b failed',me)
      if(.not.ma_push_get(mt_dbl,mem1 ,'st',l_st,istptr))
     + call errquit('esp: ma_push_get s failed',me)
c
c     determine gridpoints handled by current node
c
      if(np.gt.1) then
      ngfr=me*(ngrid/np)+1
      ngto=(me+1)*(ngrid/np)
      if(me.eq.np) ngto=ngrid
      else
      ngfr=1
      ngto=ngrid
      endif
c
c     get electrostatic potential on the grid points
c
      do 2 i=ngfr,ngto
      do 3 jsh=1,nshell
      do 4 ish=1,jsh
      if(.not.bas_cn2bfr(ibasis,ish,ilo,ihi))
     + call errquit('esp: bas_cn2bfr failed',me)
      if(.not.bas_cn2bfr(ibasis,jsh,jlo,jhi))
     + call errquit('esp: bas_cn2bfr failed',me)
      call int_1er(grid(1,i),1,one,ibasis,jsh,ibasis,ish,mem1,
     + dbl_mb(istptr),max1e,dbl_mb(ibtptr))
      ifr=ibtptr
      do 5 jx=jlo,jhi
      do 6 ix=ilo,ihi
      pot(jx,ix)=dbl_mb(ifr)
      pot(ix,jx)=pot(jx,ix)
      ifr=ifr+1
    6 continue
    5 continue
    4 continue
    3 continue
      val(i)=0.0d0
      do 333 ix=1,nbf
      do 332 jx=1,nbf
      val(i)=val(i)+dns(ix,jx)*pot(ix,jx)
  332 continue
  333 continue
c
      vt=val(i)
      imin=0
      do 9 ix=1,natoms
      dist=sqrt((grid(1,i)-x(1,ix))*(grid(1,i)-x(1,ix))+
     + (grid(2,i)-x(2,ix))*(grid(2,i)-x(2,ix))+
     + (grid(3,i)-x(3,ix))*(grid(3,i)-x(3,ix)))
      if(imin.eq.0) then
      dmin=dist
      imin=ix
      endif
      if(dmin.gt.dist) then
      dmin=dist
      imin=ix
      endif
      val(i)=val(i)+q(ix)/dist
    9 continue
      if(me.eq.0.and.util_print('where',print_debug)) then
      if(i.eq.1) write(lfnout,1003)
 1003 format(//,' Electrostatic potential on grid',/)
      write(lfnout,1004) imin,i,(grid(j,i),j=1,3),val(i)-vt,vt,val(i)
 1004 format(2i5,3f8.4,3f12.6)
      endif
    2 continue
c
c     global sum of value array
c
      if(np.gt.1) then
      if(ngfr.gt.1) then
      do 10 i=1,ngfr-1
      val(i)=zero
   10 continue
      endif
      if(ngto.lt.ngrid) then
      do 11 i=ngto+1,ngrid
      val(i)=zero
   11 continue
      endif
      call ga_dgop(mre_003,val,ngrid,'+')
      endif
c
c     deallocate memory for the integrals
c
      if(.not.ma_pop_stack(l_st))
     + call errquit('esp: ma_pop_stack st failed',me)
      if(.not.ma_pop_stack(l_bt))
     + call errquit('esp: ma_pop_stack bt failed',me)
c
c     integral termination
c
      call int_terminate()
c
      return
      end
