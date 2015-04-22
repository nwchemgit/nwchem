      subroutine esp_fastesp(ene,occ,dns,x,q,grid,val)
c
c $Id: esp_esp.F,v 1.22 2003-10-20 21:50:48 d3j191 Exp $
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
      real*8 ene(2*nbf),occ(2*nbf),dns(mbf,mbf)
      real*8 x(3,natoms),q(natoms)
      real*8 grid(3,mxgrid),val(mxgrid)
c
      integer l_bt,ibtptr,l_st,istptr
c
      integer i,j,ish,jsh,ix
      integer ilo,ihi,jlo,jhi
      real*8 dist,vt
      integer imin,nval,jshi
      real*8 dmin,fact
      character*10 today,now
      logical ldummy(1)
c
c      integer lmax2, lencart, lenprod
c      parameter (lmax2 = 2*3) ! Twice max angular momentum
c      parameter (lencart = ((lmax2+1)*(lmax2+2))/2)
c      parameter (lenprod = ((lmax2+1)*(lmax2+2)*(lmax2+3))/6)
c      real*8 dwk(lencart,-lmax2:lmax2,0:lmax2),
c     $     diwk(lenprod, -lmax2:lmax2,0:lmax2)
c
      real*8 dwk(2548),diwk(7644)
      integer npts
      if(bas_is_spherical(ibasis))
     + call md_abort('fastesp does not work with spherical basis ',me)
c
      call xlm_init()
      call anl_fit_init()
      call xlm_coeff_inv(6,dwk,diwk)
c
c     calculate density matrix, occupation and energies
c
      call esp_denmat(occ,ene)
c
c     allocate memory for the integrals
c
      if(.not.ma_push_get(mt_dbl,max1e,'bt',l_bt,ibtptr))
     + call md_abort('esp: ma_push_get b failed',me)
      if(.not.ma_push_get(mt_dbl,mem1 ,'st',l_st,istptr))
     + call md_abort('esp: ma_push_get s failed',me)
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
      if(nval.ne.ngrid) call md_abort('Incompatible er file',0)
      read(lfner) (val(i),i=1,ngrid)
      close(unit=lfner)
   12 continue
      write(lfnout,1000) jshi
 1000 format(/,' Recovering from shell ',i6,/)
      endif
      call ga_brdcst(mre_005,jshi,ma_sizeof(mt_int,1,mt_byte),0)
c
      do 3 jsh=jshi+1,nshell
      if(.not.bas_cn2bfr(ibasis,jsh,jlo,jhi))
     + call md_abort('esp: bas_cn2bfr failed',me)
      do 4 ish=1,jsh
      if(lscrn.ne.0) then
      if(schwarz_shell(ish,jsh).lt.scrtol) goto 4
      endif
      fact=1.0d0
      if(ish.ne.jsh) fact=2.0d0
      if(.not.bas_cn2bfr(ibasis,ish,ilo,ihi))
     + call md_abort('esp: bas_cn2bfr failed',me)
      npts=ngto-ngfr+1
      call esp_getden(dns,jlo,jhi,ilo,ihi)
      call potential_from_shell_pair(ibasis,jsh,ish,dns,mbf,
     + scrtol,dwk,diwk,6,npts,grid(1,ngfr),val(ngfr),fact,
     + .false.,ldummy,1)
    4 continue
      call ga_dgop(mre_003,val,ngrid,'+')
      if(me.eq.0) then
      open(unit=lfner,file=erfil,form='unformatted',status='unknown')
      rewind(unit=lfner)
      write(lfner) jsh,ngrid
      write(lfner) (val(i),i=1,ngrid)
      close(unit=lfner)
      if(util_print('output',print_medium)) then
      call swatch(today,now)
      write(lfnout,1001) today,now,jsh,nshell
 1001 format(2a10,' er shell ',i6,' of ',i6)
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
     + call md_abort('esp_denmat: ga_destroy lg_d failed',me)
      if(.not.ma_pop_stack(l_st))
     + call md_abort('esp: ma_pop_stack st failed',me)
      if(.not.ma_pop_stack(l_bt))
     + call md_abort('esp: ma_pop_stack bt failed',me)
c
c     integral termination
c
      call int_terminate()
c
      if(me.eq.0) then
      call util_file_name('grid',.false.,.false.,grdfil)
      open(unit=lfngrd,file=grdfil,form='formatted',status='unknown')
      rewind(lfngrd)
      write(lfngrd,'(i10,f20.10)') ngrid,charge
      do 10 i=1,ngrid
      write(lfngrd,'(4f20.10)') (grid(j,i),j=1,3),val(i)
   10 continue
      close(unit=lfngrd)
      endif
c
      return
      end
