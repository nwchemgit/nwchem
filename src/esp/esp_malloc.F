      subroutine esp_malloc
c
c $Id$
c
      implicit none
c
#include "esp_params.fh"
#include "esp_common.fh"
#include "bas.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      logical int_normalize
      external int_normalize
c
      integer m,n,ish,ilo,ihi
c
c     allocate memory:
c
c     l_t itptr 't'    : 16*natoms : atomic tags
c     l_x ixptr 'x'    :  3*natoms : atomic coordinates
c     l_q iqptr 'q'    :    natoms : atomic charges
c     l_qf iqfptr 'qf' :  4*natoms : fitted atomic charges
c
      ndim=2*natoms+nconst
      if(.not.ma_push_get(mt_byte,natoms*16,'t',l_t,itptr))
     + call md_abort('esp: ma_push_get t failed',me)
      if(.not.ma_push_get(mt_dbl,natoms*3,'x',l_x,ixptr))
     + call md_abort('esp: ma_push_get x failed',me)
      if(.not.ma_push_get(mt_dbl,ndim,'q',l_q,iqptr))
     + call md_abort('esp: ma_push_get q failed',me)
      if(.not.ma_push_get(mt_dbl,ndim,'qk',l_qk,iqkptr))
     + call md_abort('esp: ma_push_get qk failed',me)
      if(.not.ma_push_get(mt_dbl,6*ndim,'qf',l_qf,iqfptr))
     + call md_abort('esp: ma_push_get qf failed',me)
      if(.not.ma_push_get(mt_int,natoms,'n',l_n,inptr))
     + call md_abort('esp: ma_push_get n failed',me)
      if(.not.ma_push_get(mt_int,natoms,'l',l_l,ilptr))
     + call md_abort('esp: ma_push_get l failed',me)
c
c     get atomic coordinates
c
      call esp_getx(byte_mb(itptr),dbl_mb(ixptr),dbl_mb(iqptr),
     + int_mb(inptr))
c
      if(me.eq.0.and.newgrd.eq.0) then
        newgrd=1
        call util_file_name('grid',.false.,.false.,grdfil)
        open(unit=lfngrd,file=grdfil,form='formatted',
     +       status='old',err=9)
        read(lfngrd,'(i10)',err=9) ngrid
        mxgrid=ngrid
        newgrd=0
        close(unit=lfngrd)
    9   continue
      endif
      call ga_brdcst(mre_011,newgrd,ma_sizeof(mt_int,1,mt_byte),0)
      call ga_brdcst(mre_012,mxgrid,ma_sizeof(mt_int,1,mt_byte),0)
      ngrid=mxgrid
c
      if(newgrd.ne.0.or.denmat.ne.0) then
c
c     initialize integrals
c
        if(.not.int_normalize(irtdb,ibasis))
     +   call md_abort('esp: int_normalize failed',me)
        call int_init(irtdb,1,ibasis)
c
c       get number of shells
c
        if(.not.bas_numcont(ibasis,nshell))
     +   call md_abort('esp: bas_numcont failed',me)
        call int_mem_1e(max1e,mem1)
c
        mbf=0
        do 1 ish=1,nshell
          if(.not.bas_cn2bfr(ibasis,ish,ilo,ihi))
     +     call md_abort('esp: bas_cn2bfr failed',me)
          mbf=max(mbf,ihi-ilo+1)
    1   continue
      else
        nbf=1
        mbf=1
      endif
c
c     allocate memory:
c
c     l_d idptr 'd' : nbf*nbf : density matrix
c     l_e ieptr 'e' :   2*nbf : orbital energies
c     l_o ioptr 'o' :   2*nbf : orbital occupation
c
      if(.not.ma_push_get(mt_dbl,mbf*mbf,'d',l_d,idptr))
     + call md_abort('esp: ma_push_get d failed',me)
      if(.not.ma_push_get(mt_dbl,2*nbf,'e',l_e,ieptr))
     + call md_abort('esp: ma_push_get e failed',me)
      if(.not.ma_push_get(mt_dbl,2*nbf,'o',l_o,ioptr))
     + call md_abort('esp: ma_push_get o failed',me)
c
c     allocate memory:
c
c     l_a iaptr 'A' : (natoms+nconst)^2 : linear coefficient matrix
c     l_b ibptr 'b' : (natoms+nconst)   : solution vector
c     l_c icptr 'c' : (natoms+nconst)   : result vector
c
      if(.not.ma_push_get(mt_dbl,ndim*ndim,'A',l_a,iaptr))
     + call md_abort('esp: ma_push_get A failed',me)
      if(.not.ma_push_get(mt_dbl,ndim*ndim,'Ak',l_ak,iakptr))
     + call md_abort('esp: ma_push_get Ak failed',me)
      if(.not.ma_push_get(mt_dbl,ndim,'b',l_b,ibptr))
     + call md_abort('esp: ma_push_get b failed',me)
      if(.not.ma_push_get(mt_dbl,ndim,'c',l_c,icptr))
     + call md_abort('esp: ma_push_get c failed',me)
c
c     check how much memory is available
c
      m=ma_inquire_avail(mt_byte)/2
c
c     reserve memory for several workarrays
c
      m=m-ndim*(4*ma_sizeof(mt_int,1,mt_byte)+
     + ma_sizeof(mt_dbl,1,mt_byte))
      m=m-(max1e+mem1)*ma_sizeof(mt_dbl,1,mt_byte)
      n=max(2*nbf,256)
      m=m-(3*n+1)*n*ma_sizeof(mt_dbl,1,mt_byte)
      if(newgrd.ne.0) call esp_grdcnt(dbl_mb(ixptr),dbl_mb(iqptr))
c
c     allocate memory for the grid and potential
c
c     l_g igptr 'g' : 3*mxgrid : grid point coordinates
c     l_v ivptr 'v' :   mxgrid : grid point values
c
      if(.not.ma_push_get(mt_dbl,mxgrid*3,'g',l_g,igptr))
     + call md_abort('esp: ma_push_get g failed',me)
      if(.not.ma_push_get(mt_dbl,mxgrid,'v',l_v,ivptr))
     + call md_abort('esp: ma_push_get v failed',me)
c
      if(me.eq.0.and.newgrd.eq.0)
     + call esp_getg(dbl_mb(igptr),dbl_mb(ivptr))
c
      return
      end
