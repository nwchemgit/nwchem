      subroutine esp_esp(ene,occ,dns,pot,x,q,grid,val)
c
c $Id: esp_esp.F,v 1.5 1998-01-13 19:38:09 d3j191 Exp $
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
c
      logical movecs_read,movecs_read_header
      real*8 ddot
      external movecs_read,movecs_read_header,ddot
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
c     read orbital occupations and energies
c
      call util_file_name('movecs',.false.,.false.,movfil)
      if(.not.movecs_read_header(movfil,title,basnam,scftyp,nbf,nset,
     + nmo,2)) call errquit('esp: movecs_read_header failed',me)
      if(.not.movecs_read(movfil,1,occ,ene,lg_v))
     + call errquit('eps: movecs_read failed',me)
c
      nocc=0
      charge=chgnuc
      do 1 i=1,nbf
      if(occ(i).gt.small) charge=charge-occ(i)
    1 continue
c
c     get AO density matrix
c
      call esp_denmat(igeom,ibasis,lg_v,lg_d,occ)
c
c     get density matrix
c
      call ga_get(lg_d,1,nbf,1,nbf,dns,nbf)
c
c     initialize integrals
c
      call int_init(irtdb,1,ibasis)
c
c     number of shells
c
      if(.not.bas_numcont(ibasis,nshell))
     + call errquit('esp: bas_numcont failed',me)
c
c     allocate memory for the integrals
c
      call int_mem_1e(max1e,mem1)
      if(.not.ma_push_get(mt_dbl,max1e,'bt',l_bt,ibtptr))
     + call errquit('esp: ma_push_get b failed',me)
      if(.not.ma_push_get(mt_dbl,mem1 ,'st',l_st,istptr))
     + call errquit('esp: ma_push_get s failed',me)
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
      val(i)=ddot(nbf*nbf,dns,1,pot,1)
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
