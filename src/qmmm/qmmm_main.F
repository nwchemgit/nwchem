      subroutine qmmm_forces(inrtdb,mwm,nwm,mwa,nwa,iwz,xw,fw,
     + msa,nsa,isat,isdt,isq,isz,xs,fs,energy,lnkatm,nbq,lesp)
c
c $Id: qmmm_main.F,v 1.17 2004-09-30 06:46:30 marat Exp $
c
      implicit none
#include "errquit.fh"
c
#include "rtdb.fh"
#include "geom.fh"
#include "util.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "qmmm_common.fh"
c
      external task_gradient,esp
      logical task_gradient,esp
c
      integer cf_element
      external cf_element
c
      integer inrtdb
      integer mwm,nwm,mwa,nwa,msa,nsa,lnkatm,nbq
      integer iwz(mwm),isat(msa),isdt(msa),isq(msa),isz(msa)
      real*8 energy,xw(mwm,3,mwa),fw(mwm,3,mwa),xs(msa,3),fs(msa,3)
c
      integer i,np,ioffq,ioffz,igeom
      integer i_ioff,l_ioff,i_c,l_c,i_q,l_q,i_t,l_t
      integer i_i,l_i,i_j,l_j,i_m,l_m
      integer nzatom,nqatom,nzatomt,nqatomt
      logical lesp
      character*32 optimize
      logical qmmm_init
c
      if (.not. rtdb_get(irtdb,'qmmm:init',mt_log,1,qmmm_init))
     $     qmmm_init=.false.

      if(.not.qmmm_init) return

      if (.not.rtdb_cget(inrtdb,'qmmm:optimize',1,optimize))
     + call errquit('qmmm_input: failed to get qmmm:optimize',
     &       0,RTDB_ERR)

      if(optimize.eq.'qm') then
        write(*,*) "optimize qm"
        energy = 0
        return
      end if
      
      call qmmm_get_qm_coord_mm(inrtdb)
      return


      return
      end

