      subroutine esp_malloc
c
c $Id: esp_malloc.F,v 1.3 1998-08-18 16:53:54 d3j191 Exp $
c
      implicit none
c
#include "esp_params.fh"
#include "esp_common.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "mafdecls.fh"
c
      integer i,j
c
c     allocate memory:
c
c     l_t itptr 't'    : 16*natoms : atomic tags
c     l_x ixptr 'x'    :  3*natoms : atomic coordinates
c     l_q iqptr 'q'    :    natoms : atomic charges
c     l_qf iqfptr 'qf' :  4*natoms : fitted atomic charges
c
      ndim=natoms+nconst
      if(.not.ma_push_get(mt_byte,natoms*16,'t',l_t,itptr))
     + call errquit('esp: ma_push_get t failed',me)
      if(.not.ma_push_get(mt_dbl,natoms*3,'x',l_x,ixptr))
     + call errquit('esp: ma_push_get x failed',me)
      if(.not.ma_push_get(mt_dbl,ndim,'q',l_q,iqptr))
     + call errquit('esp: ma_push_get q failed',me)
      if(.not.ma_push_get(mt_dbl,ndim,'qk',l_qk,iqkptr))
     + call errquit('esp: ma_push_get qk failed',me)
      if(.not.ma_push_get(mt_dbl,4*ndim,'qf',l_qf,iqfptr))
     + call errquit('esp: ma_push_get qf failed',me)
      if(.not.ma_push_get(mt_int,natoms,'n',l_n,inptr))
     + call errquit('esp: ma_push_get n failed',me)
c
c     get atomic coordinates
c
      call esp_getx(byte_mb(itptr),dbl_mb(ixptr),dbl_mb(iqptr),
     + int_mb(inptr))
c
c     allocate memory for the grid and potential
c
c     l_g igptr 'g' : 3*mxgrid : grid point coordinates
c     l_v ivptr 'v' :   mxgrid : grid point values
c
      if(.not.ma_push_get(mt_dbl,mxgrid*3,'g',l_g,igptr))
     + call errquit('esp: ma_push_get g failed',me)
      if(.not.ma_push_get(mt_dbl,mxgrid,'v',l_v,ivptr))
     + call errquit('esp: ma_push_get v failed',me)
c
c     allocate memory:
c
c     l_d idptr 'd' : nbf*nbf : density matrix
c     l_p ipptr 'p' : nbf*nbf : potential matrix
c     l_e ieptr 'e' :     nbf : orbital energies
c     l_o ioptr 'o' :     nbf : orbital occupation
c
      if(.not.ma_push_get(mt_dbl,nbf*nbf,'d',l_d,idptr))
     + call errquit('esp: ma_push_get d failed',me)
      if(.not.ma_push_get(mt_dbl,nbf*nbf,'p',l_p,ipptr))
     + call errquit('esp: ma_push_get p failed',me)
      if(.not.ma_push_get(mt_dbl,nbf,'e',l_e,ieptr))
     + call errquit('esp: ma_push_get e failed',me)
      if(.not.ma_push_get(mt_dbl,nbf,'o',l_o,ioptr))
     + call errquit('esp: ma_push_get o failed',me)
c
c     allocate memory:
c
c     l_a iaptr 'A' : (natoms+nconst)^2 : linear coefficient matrix
c     l_b ibptr 'b' : (natoms+nconst)   : solution vector
c     l_c icptr 'c' : (natoms+nconst)   : result vector
c
      if(.not.ma_push_get(mt_dbl,ndim*ndim,'A',l_a,iaptr))
     + call errquit('esp: ma_push_get A failed',me)
      if(.not.ma_push_get(mt_dbl,ndim*ndim,'Ak',l_ak,iakptr))
     + call errquit('esp: ma_push_get Ak failed',me)
      if(.not.ma_push_get(mt_dbl,ndim,'b',l_b,ibptr))
     + call errquit('esp: ma_push_get b failed',me)
      if(.not.ma_push_get(mt_dbl,ndim,'c',l_c,icptr))
     + call errquit('esp: ma_push_get c failed',me)
c
      return
      end

