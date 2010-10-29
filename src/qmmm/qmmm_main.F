      subroutine qmmm_forces(inrtdb,mwm,nwm,mwa,nwa,iwz,xw,fw,
     + msa,nsa,isat,isdt,isq,isz,xs,fs,energy)
c
c $Id$
c
      implicit none
#include "errquit.fh"
c
#include "rtdb.fh"
#include "geom.fh"
#include "util.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "qmmm.fh"
#include "qmmm_params.fh"
c
c
      integer inrtdb
      integer mwm,nwm,mwa,nwa,msa,nsa
      integer iwz(mwm),isat(msa),isdt(msa),isq(msa),isz(msa)
      real*8 energy,xw(mwm,3,mwa),fw(mwm,3,mwa),xs(msa,3),fs(msa,3)
c
      character*32 interface
      character*32 optimization
      character*255 theory
      character*30 operation
      logical qmmm_init
      character*32 pname
      double precision eatoms
      logical status
      logical do_grad

      logical qmmm_energy_gradient
      external qmmm_energy_gradient

      pname = 'qmmm_main'
      if(qmmm_print_debug()) 
     >  write(*,*) "in "//pname
c
c     make sure qmmm is initialized
c     -----------------------------
      if (.not. rtdb_get(inrtdb,'qmmm:init',mt_log,1,qmmm_init))
     $     qmmm_init=.false.
      if(.not.qmmm_init) return
c
c     nothing should be calculated if
c     qm module (e.g. driver is in charge)
c     -----------------------------------
      interface = qmmm_get_interface()
      if(interface.eq.'qm') then
      if(qmmm_print_debug())
     > write(*,*) "exiting out of qmmm_forces since interface is set qm"
        energy = 0
        return
      end if

      if (.not. rtdb_cget(inrtdb, 'qmmm:operation', 1, operation))
     $     operation = ' '

      if(operation.eq.'energy') then
        do_grad=.false.
      else
        do_grad=.true.
      end if

      call  qmmm_bq_data_reload()
      if(.not.qmmm_energy_gradient(inrtdb,do_grad)) 
     $   call errquit(pname//'failed qmmm_energy_gradient',0,0)

      if(do_grad) call qmmm_cons_fixed()

      if (.not.rtdb_get(inrtdb,'qmmm:uqmatm',mt_dbl,1,eatoms))
     $     call errquit('qmmm: failed getting  ref energy',0,RTDB_ERR)

      if (.not. rtdb_get(inrtdb,'qmmm:qm_energy',mt_dbl,1,energy))
     $     call errquit('qmmm: failed getting qm energy',0,RTDB_ERR)

      energy = energy-eatoms
      energy=energy*cau2kj
 
      return
      end

