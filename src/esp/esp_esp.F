      subroutine esp_esp(ene,occ,dns,x,q,grid,val)
c
c $Id: esp_esp.F,v 1.9 1999-01-27 18:18:00 d3j191 Exp $
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
#include "schwarz.fh"
c
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
c
      real*8 ene(nbf),occ(nbf),dns(nbf,mbf)
      real*8 x(3,natoms),q(natoms)
      real*8 grid(3,mxgrid),val(mxgrid)
c
      integer l_bt,ibtptr,l_st,istptr
c
      integer i,j,ish,jsh,ix,jx,ifr
      integer ilo,ihi,jlo,jhi
      real*8 dist,vt
      integer imin,nval,jshi
      real*8 dmin,fact
      character*10 today,now
c
c     calculate density matrix, occupation and energies
c
      call esp_denmat(occ,ene)
c
c     allocate memory for the integrals
c
      if(.not.ma_push_get(mt_dbl,max1e,'bt',l_bt,ibtptr))
     + call errquit('esp: ma_push_get b failed',me)
      if(.not.ma_push_get(mt_dbl,mem1 ,'st',l_st,istptr))
     + call errquit('esp: ma_push_get s failed',me)
c
      if(lscrn.ne.0) call schwarz_init(igeom,ibasis)
c
c     get electrostatic potential on the grid points
c
      do 2 i=1,ngrid
      val(i)=0.0d0
    2 continue
c
      jshi=0
      if(me.eq.0) then
      call util_file_name('er',.false.,.false.,erfil)
      open(unit=lfner,file=erfil,form='unformatted',status='old',
     + err=12)
      read(lfner) jshi,nval
      if(nval.ne.ngrid) call errquit('Incompatible er file',0)
      read(lfner) (val(i),i=1,ngrid)
      close(unit=lfner)
   12 continue
      write(lfnout,1000) jshi
 1000 format(/,' Recovering from shell ',i6)
      endif
      call ga_brdcst(mre_005,jshi,ma_sizeof(mt_int,1,mt_byte),0)
c
      do 3 jsh=jshi+1,nshell
      if(.not.bas_cn2bfr(ibasis,jsh,jlo,jhi))
     + call errquit('esp: bas_cn2bfr failed',me)
      call esp_getden(dns,1,jhi,jlo,jhi)
      do 4 ish=1,jsh
      if(lscrn.ne.0) then
      if(schwarz_shell(ish,jsh).lt.1.0d-5) goto 4
      endif
      fact=1.0d0
      if(ish.ne.jsh) fact=2.0d0
      if(.not.bas_cn2bfr(ibasis,ish,ilo,ihi))
     + call errquit('esp: bas_cn2bfr failed',me)     
      do 5 i=ngfr,ngto
      call int_1er(grid(1,i),1,one,ibasis,jsh,ibasis,ish,mem1,
     + dbl_mb(istptr),max1e,dbl_mb(ibtptr))
      ifr=ibtptr
      do 6 jx=1,jhi-jlo+1
      do 7 ix=ilo,ihi
      val(i)=val(i)+fact*dns(ix,jx)*dbl_mb(ifr)
      ifr=ifr+1
    7 continue
    6 continue
    5 continue
    4 continue
      call ga_dgop(mre_003,val,ngrid,'+')
      if(me.eq.0) then
      open(unit=lfner,file=erfil,form='unformatted',status='unknown')
      rewind(unit=lfner)
      write(lfner) jsh,ngrid
      write(lfner) (val(i),i=1,ngrid)
      close(unit=lfner)
      if(util_print('output',print_medium)) then
      call arg_swatch(today,now)
      write(lfnout,1001) today,now,jsh,nshell
 1001 format(/,2a10,' er shell ',i6,' of ',i6)
      call util_flush(lfnout)
      endif
      else
      do 13 i=1,ngrid
      val(i)=zero
   13 continue
      endif
    3 continue
c
      if(me.eq.0) then
      open(unit=lfner,file=erfil,form='unformatted',status='unknown')
      close(unit=lfner,status='delete')
      write(lfnout,1002)
 1002 format(/,' Recovery file deleted',/)
      endif
c
      do 8 i=ngfr,ngto
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
    8 continue
c
      if(lscrn.ne.0) call schwarz_tidy()
c
c     global sum of value array
c
      if(np.gt.1) then
      call ga_dgop(mre_006,val,ngrid,'+')
      endif
c
c     deallocate memory for the integrals
c
      if(.not.ga_destroy(lg_d))
     + call errquit('esp_denmat: ga_destroy lg_d failed',me)
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
