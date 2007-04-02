*
* $Id: force_fd.F,v 1.2 2007-04-02 15:38:40 bylaska Exp $
*
      logical function force_fd(rtdb)
      implicit none
      integer rtdb
      logical value

#include "mafdecls.fh"
#include "rtdb.fh"
#include "inp.fh"
#include "util.fh"
#include "geom.fh"
#include "errquit.fh"

      logical oprint
      integer taskid,MASTER
      parameter (MASTER=0)

      integer i
      integer r1(2),geom,nion
      integer r2(2),few(2),fvl(2),fvnl(2),fxc(2),fLJ(2),fQ(2)


      double precision alpha,E

      integer date_len
      character*40 date_string

*     *** external functions ****
      logical  pspw_noit_energy
      external pspw_noit_energy

      call Parallel_taskid(taskid)
      oprint = (taskid.eq.MASTER)

*     **** initialize and read in geometry object ****
      value=.true.
      value     = value.and.geom_create(geom,'geometry')
      value     = value.and.geom_rtdb_load(rtdb,geom,'geometry')
      value     = value.and.geom_ncent(geom,nion)
      value     = value.and.
     >            MA_alloc_get(mt_dbl,(3*nion),'r1',r1(2),r1(1))
      value     = value.and.
     >            MA_alloc_get(mt_dbl,(3*nion),'r2',r2(2),r2(1))

      value     = value.and.
     >            MA_alloc_get(mt_dbl,(3*nion),'fvl',fvl(2),fvl(1))
      value     = value.and.
     >            MA_alloc_get(mt_dbl,(3*nion),'fvnl',fvnl(2),fvnl(1))
      value     = value.and.
     >            MA_alloc_get(mt_dbl,(3*nion),'few',few(2),few(1))
      value     = value.and.
     >            MA_alloc_get(mt_dbl,(3*nion),'fxc',fxc(2),fxc(1))
      value     = value.and.
     >            MA_alloc_get(mt_dbl,(3*nion),'fQ',fLJ(2),fLJ(1))
      value     = value.and.
     >            MA_alloc_get(mt_dbl,(3*nion),'fQ',fQ(2),fQ(1))
      if(.not.value) call errquit
     E   (' stressfd: geom malloc failure ',0,0)

      if (.not. geom_cart_coords_get(geom,dbl_mb(r1(1))))
     >   call errquit('force_fd: geom cart?',0,GEOM_ERR)


*     **** get finite difference step and initialize ****
      if (.not.
     >    rtdb_get(rtdb,'pspw:force_alpha',mt_dbl,1,alpha))
     >  alpha = 0.00100d0
      call dcopy(3*nion,0.0d0,0,dbl_mb(fvl(1)),1)
      call dcopy(3*nion,0.0d0,0,dbl_mb(fvnl(1)),1)
      call dcopy(3*nion,0.0d0,0,dbl_mb(few(1)),1)
      call dcopy(3*nion,0.0d0,0,dbl_mb(fxc(1)),1)
      call dcopy(3*nion,0.0d0,0,dbl_mb(fLJ(1)),1)
      call dcopy(3*nion,0.0d0,0,dbl_mb(fQ(1)),1)


      do i=1,(3*nion)


*     *****************************
*     **** + finite difference ****
*     *****************************
      call dcopy(3*nion,dbl_mb(r1(1)),1,dbl_mb(r2(1)),1)
      dbl_mb(r2(1)+i-1) = dbl_mb(r2(1)+i-1) + alpha

*     **** update geometry object ****
      call util_date(date_string)
      date_len = inp_strlen(date_string)
      if (oprint) write(6,*)
      if (oprint) write(6,
     &   '(1x,a12,i4,1x,f8.3,1x,a4,1x,a10,f10.1,6x,a5,a27)')
     &   'force i:',i,dbl_mb(r2(1)+i-1),'(+)',
     &   'wall time:',util_wallsec(),
     &   'date:',date_string(1:date_len)
      if (oprint) call util_flush(6)
      if (.not. geom_cart_coords_set(geom,dbl_mb(r2(1))))
     >   call errquit('force_fd: geom cart?',0,GEOM_ERR)
      value = value.and.geom_rtdb_store(rtdb,geom,'geometry')
      value = value.and.pspw_noit_energy(rtdb)


      value = value.and.
     >        rtdb_get(rtdb,'pspw:E_ewald',mt_dbl,1,E)
      dbl_mb(few(1)+i-1) = E

      value = value.and.
     >        rtdb_get(rtdb,'pspw:E_local',mt_dbl,1,E)
      dbl_mb(fvl(1)+i-1) = E

      value = value.and.
     >        rtdb_get(rtdb,'pspw:E_nonlocal',mt_dbl,1,E)
      dbl_mb(fvnl(1)+i-1) = E

      value = value.and.
     >        rtdb_get(rtdb,'pspw:E_xc',mt_dbl,1,E)
      dbl_mb(fxc(1)+i-1) = E

      value = value.and.
     >        rtdb_get(rtdb,'pspw:E_LJ',mt_dbl,1,E)
      dbl_mb(fLJ(1)+i-1) = E

      value = value.and.
     >        rtdb_get(rtdb,'pspw:E_Q',mt_dbl,1,E)
      dbl_mb(fQ(1)+i-1) = E


*     *****************************
*     **** - finite difference ****
*     *****************************
      call dcopy(3*nion,dbl_mb(r1(1)),1,dbl_mb(r2(1)),1)
      dbl_mb(r2(1)+i-1) = dbl_mb(r2(1)+i-1) - alpha

*     **** update geometry object ****
      call util_date(date_string)
      date_len = inp_strlen(date_string)
      if (oprint) write(6,*)
      if (oprint) write(6,
     &   '(1x,a12,i4,1x,f8.3,1x,a4,1x,a10,f10.1,6x,a5,a27)')
     &   'force i:',i,dbl_mb(r2(1)+i-1),'(-)',
     &   'wall time:',util_wallsec(),
     &   'date:',date_string(1:date_len)
      if (oprint) call util_flush(6)
      if (.not. geom_cart_coords_set(geom,dbl_mb(r2(1))))
     >   call errquit('force_fd: geom cart?',0,GEOM_ERR)
      value = value.and.geom_rtdb_store(rtdb,geom,'geometry')
    
      value = pspw_noit_energy(rtdb)


      value = value.and.rtdb_get(rtdb,'pspw:E_ewald',mt_dbl,1,E)
      dbl_mb(few(1)+i-1) = dbl_mb(few(1)+i-1) - E

      value = value.and.rtdb_get(rtdb,'pspw:E_local',mt_dbl,1,E)
      dbl_mb(fvl(1)+i-1) = dbl_mb(fvl(1)+i-1) - E

      value = value.and.rtdb_get(rtdb,'pspw:E_nonlocal',mt_dbl,1,E)
      dbl_mb(fvnl(1)+i-1) = dbl_mb(fvnl(1)+i-1) - E

      value = value.and.rtdb_get(rtdb,'pspw:E_xc',mt_dbl,1,E)
      dbl_mb(fxc(1)+i-1) = dbl_mb(fxc(1)+i-1) - E

      value = value.and.rtdb_get(rtdb,'pspw:E_LJ',mt_dbl,1,E)
      dbl_mb(fLJ(1)+i-1) = dbl_mb(fLJ(1)+i-1) - E

      value = value.and.rtdb_get(rtdb,'pspw:E_Q',mt_dbl,1,E)
      dbl_mb(fQ(1)+i-1) = dbl_mb(fQ(1)+i-1) - E

      end do
      call dscal(3*nion,(-0.5d0/alpha),dbl_mb(few(1)),1)
      call dscal(3*nion,(-0.5d0/alpha),dbl_mb(fvl(1)),1)
      call dscal(3*nion,(-0.5d0/alpha),dbl_mb(fvnl(1)),1)
      call dscal(3*nion,(-0.5d0/alpha),dbl_mb(fxc(1)),1)
      call dscal(3*nion,(-0.5d0/alpha),dbl_mb(fLJ(1)),1)
      call dscal(3*nion,(-0.5d0/alpha),dbl_mb(fQ(1)),1)

      call write_force(nion,dbl_mb(few(1)),"ewald FD")
      call write_force(nion,dbl_mb(fvl(1)),"local FD")
      call write_force(nion,dbl_mb(fvnl(1)),"nonlocal FD")
      call write_force(nion,dbl_mb(fxc(1)),"xc FD")
      call write_force(nion,dbl_mb(fLJ(1)),"LJ FD")
      call write_force(nion,dbl_mb(fQ(1)), "Q FD")

      call daxpy(3*nion,1.0d0,dbl_mb(fvl(1)),1,dbl_mb(few(1)),1)
      call daxpy(3*nion,1.0d0,dbl_mb(fvnl(1)),1,dbl_mb(few(1)),1)
      call daxpy(3*nion,1.0d0,dbl_mb(fxc(1)),1,dbl_mb(few(1)),1)
      call daxpy(3*nion,1.0d0,dbl_mb(fLJ(1)),1,dbl_mb(few(1)),1)
      call daxpy(3*nion,1.0d0,dbl_mb(fQ(1)),1,dbl_mb(few(1)),1)
      call write_force(nion,dbl_mb(few(1)),"total FD")



      if (.not. geom_cart_coords_set(geom,dbl_mb(r1(1))))
     >   call errquit('force_fd: geom cart?',0,GEOM_ERR)
      !value = value.and.geom_rtdb_delete(rtdb,'geometry')
      value = value.and.geom_rtdb_store(rtdb,geom,'geometry')
      value = value.and.geom_destroy(geom)
      value = value.and.MA_free_heap(r1(2))
      value = value.and.MA_free_heap(r2(2))
      value = value.and.MA_free_heap(few(2))
      value = value.and.MA_free_heap(fvl(2))
      value = value.and.MA_free_heap(fvnl(2))
      value = value.and.MA_free_heap(fxc(2))
      value = value.and.MA_free_heap(fLJ(2))
      value = value.and.MA_free_heap(fQ(2))
 
      force_fd = value
      return 
      end

      subroutine write_force(nion,fion,name)
      implicit none
      integer nion
      real*8  fion(3,*)
      character*(*) name

      integer i,k
      character*4 ion_aname
      external    ion_aname

      write(*,1100) name
      write(*,1190) (i,(fion(k,i),k=1,3),i=1,nion)
      write(*,*)

      return
 1100 FORMAT(//A,' Ion Forces:')
 1190 FORMAT(5X, I4, ' (',3F11.5,' )')
      end
