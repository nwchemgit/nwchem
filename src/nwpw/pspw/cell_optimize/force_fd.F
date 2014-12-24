*
* $Id$
*
      logical function force_fd(rtdb)
      implicit none
      integer rtdb
      logical value

#include "bafdecls.fh"
#include "btdb.fh"
#include "inp.fh"
#include "util.fh"
#include "geom.fh"
#include "errquit.fh"

      logical oprint
      integer taskid,MASTER
      parameter (MASTER=0)

      integer i
      integer r1(2),geom,nion
      integer r2(2),few(2),fvl(2),fvnl(2),fxc(2),fLJ(2),fQ(2),fspring(2)
      integer fdftu(2),fdisp(2)
      integer fkloc(2),fkcore(2),fvalcore(2),fioncore(2),fcmpvl(2)
      integer fvlocatom(2)
      integer fhatom(2),fcmp(2),fcmppw(2),fxcatom(2),fpxcatom(2),fall(2)


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
     >            BA_alloc_get(mt_dbl,(3*nion),'r1',r1(2),r1(1))
      value     = value.and.
     >            BA_alloc_get(mt_dbl,(3*nion),'r2',r2(2),r2(1))

      value     = value.and.
     >            BA_alloc_get(mt_dbl,(3*nion),'fvl',fvl(2),fvl(1))
      value     = value.and.
     >            BA_alloc_get(mt_dbl,(3*nion),'fvnl',fvnl(2),fvnl(1))
      value     = value.and.
     >            BA_alloc_get(mt_dbl,(3*nion),'few',few(2),few(1))
      value     = value.and.
     >            BA_alloc_get(mt_dbl,(3*nion),'fxc',fxc(2),fxc(1))
      value     = value.and.
     >            BA_alloc_get(mt_dbl,(3*nion),'fJ',fLJ(2),fLJ(1))
      value     = value.and.
     >            BA_alloc_get(mt_dbl,(3*nion),'fQ',fQ(2),fQ(1))
      value     = value.and.
     >     BA_alloc_get(mt_dbl,(3*nion),'fspring',fspring(2),fspring(1))
      value     = value.and.
     >     BA_alloc_get(mt_dbl,(3*nion),'fdftu',fdftu(2),fdftu(1))
      value     = value.and.
     >     BA_alloc_get(mt_dbl,(3*nion),'fdisp',fdisp(2),fdisp(1))
      value     = value.and.
     >     BA_alloc_get(mt_dbl,(3*nion),'fkloc',fkloc(2),fkloc(1))
      value     = value.and.
     >     BA_alloc_get(mt_dbl,(3*nion),'fkcore',fkcore(2),fkcore(1))
      value     = value.and.
     > BA_alloc_get(mt_dbl,(3*nion),'fvalcore',fvalcore(2),fvalcore(1))
      value     = value.and.
     > BA_alloc_get(mt_dbl,(3*nion),'fioncore',fioncore(2),fioncore(1))
      value     = value.and.
     >     BA_alloc_get(mt_dbl,(3*nion),'fcmpvl',fcmpvl(2),fcmpvl(1))
      value     = value.and.
     >     BA_alloc_get(mt_dbl,(3*nion),
     >                  'fvlocatom',fvlocatom(2),fvlocatom(1))
      value     = value.and.
     >     BA_alloc_get(mt_dbl,(3*nion),'fhatom',fhatom(2),fhatom(1))
      value     = value.and.
     >     BA_alloc_get(mt_dbl,(3*nion),'fcmp',fcmp(2),fcmp(1))
      value     = value.and.
     >     BA_alloc_get(mt_dbl,(3*nion),'fcmppw',fcmppw(2),fcmppw(1))
      value     = value.and.
     >     BA_alloc_get(mt_dbl,(3*nion),'fxcatom',fxcatom(2),fxcatom(1))
      value     = value.and.
     >  BA_alloc_get(mt_dbl,(3*nion),'fpxcatom',fpxcatom(2),fpxcatom(1))
      value     = value.and.
     >     BA_alloc_get(mt_dbl,(3*nion),'fall',fall(2),fall(1))

      if(.not.value) call errquit
     E   (' stressfd: geom malloc failure ',0,0)

      if (.not. geom_cart_coords_get(geom,dbl_mb(r1(1))))
     >   call errquit('force_fd: geom cart?',0,GEOM_ERR)


*     **** get finite difference step and initialize ****
      if (.not.
     >    btdb_get(rtdb,'pspw:force_alpha',mt_dbl,1,alpha))
     >  alpha = 0.00100d0
      call dcopy(3*nion,0.0d0,0,dbl_mb(fvl(1)),1)
      call dcopy(3*nion,0.0d0,0,dbl_mb(fvnl(1)),1)
      call dcopy(3*nion,0.0d0,0,dbl_mb(few(1)),1)
      call dcopy(3*nion,0.0d0,0,dbl_mb(fxc(1)),1)
      call dcopy(3*nion,0.0d0,0,dbl_mb(fLJ(1)),1)
      call dcopy(3*nion,0.0d0,0,dbl_mb(fQ(1)),1)
      call dcopy(3*nion,0.0d0,0,dbl_mb(fspring(1)),1)
      call dcopy(3*nion,0.0d0,0,dbl_mb(fdftu(1)),1)
      call dcopy(3*nion,0.0d0,0,dbl_mb(fdisp(1)),1)
      call dcopy(3*nion,0.0d0,0,dbl_mb(fkloc(1)),1)
      call dcopy(3*nion,0.0d0,0,dbl_mb(fkcore(1)),1)
      call dcopy(3*nion,0.0d0,0,dbl_mb(fvalcore(1)),1)
      call dcopy(3*nion,0.0d0,0,dbl_mb(fioncore(1)),1)
      call dcopy(3*nion,0.0d0,0,dbl_mb(fcmpvl(1)),1)
      call dcopy(3*nion,0.0d0,0,dbl_mb(fvlocatom(1)),1)
      call dcopy(3*nion,0.0d0,0,dbl_mb(fhatom(1)),1)
      call dcopy(3*nion,0.0d0,0,dbl_mb(fcmp(1)),1)
      call dcopy(3*nion,0.0d0,0,dbl_mb(fcmppw(1)),1)
      call dcopy(3*nion,0.0d0,0,dbl_mb(fxcatom(1)),1)
      call dcopy(3*nion,0.0d0,0,dbl_mb(fpxcatom(1)),1)
      call dcopy(3*nion,0.0d0,0,dbl_mb(fall(1)),1)


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
     >        btdb_get(rtdb,'pspw:E_ewald',mt_dbl,1,E)
      dbl_mb(few(1)+i-1) = E

      value = value.and.
     >        btdb_get(rtdb,'pspw:E_local',mt_dbl,1,E)
      dbl_mb(fvl(1)+i-1) = E

      value = value.and.
     >        btdb_get(rtdb,'pspw:E_nonlocal',mt_dbl,1,E)
      dbl_mb(fvnl(1)+i-1) = E

      value = value.and.
     >        btdb_get(rtdb,'pspw:E_xc',mt_dbl,1,E)
      dbl_mb(fxc(1)+i-1) = E

      value = value.and.
     >        btdb_get(rtdb,'pspw:E_LJ',mt_dbl,1,E)
      dbl_mb(fLJ(1)+i-1) = E

      value = value.and.
     >        btdb_get(rtdb,'pspw:E_Q',mt_dbl,1,E)
      dbl_mb(fQ(1)+i-1) = E

      value = value.and.
     >        btdb_get(rtdb,'pspw:E_spring',mt_dbl,1,E)
      dbl_mb(fspring(1)+i-1) = E

      value = value.and.
     >        btdb_get(rtdb,'pspw:E_dftu',mt_dbl,1,E)
      dbl_mb(fdftu(1)+i-1) = E

      value = value.and.
     >        btdb_get(rtdb,'pspw:E_disp',mt_dbl,1,E)
      dbl_mb(fdisp(1)+i-1) = E

      value = value.and.
     >        btdb_get(rtdb,'pspw:kin_loc',mt_dbl,1,E)
      dbl_mb(fkloc(1)+i-1) = E

      value = value.and.
     >        btdb_get(rtdb,'pspw:kin_core',mt_dbl,1,E)
      dbl_mb(fkcore(1)+i-1) = E

      value = value.and.
     >        btdb_get(rtdb,'pspw:val_core',mt_dbl,1,E)
      dbl_mb(fvalcore(1)+i-1) = E

      value = value.and.
     >        btdb_get(rtdb,'pspw:ion_core',mt_dbl,1,E)
      dbl_mb(fioncore(1)+i-1) = E

      value = value.and.
     >        btdb_get(rtdb,'pspw:E_ncmp_vloc',mt_dbl,1,E)
      dbl_mb(fcmpvl(1)+i-1) = E

      value = value.and.
     >        btdb_get(rtdb,'pspw:vloc_atom',mt_dbl,1,E)
      dbl_mb(fvlocatom(1)+i-1) = E

      value = value.and.
     >        btdb_get(rtdb,'pspw:E_hartree_atom',mt_dbl,1,E)
      dbl_mb(fhatom(1)+i-1) = E

      value = value.and.
     >        btdb_get(rtdb,'pspw:E_cmp_cmp',mt_dbl,1,E)
      dbl_mb(fcmp(1)+i-1) = E

      value = value.and.
     >        btdb_get(rtdb,'pspw:E_cmp_pw',mt_dbl,1,E)
      dbl_mb(fcmppw(1)+i-1) = E

      value = value.and.
     >        btdb_get(rtdb,'pspw:E_exc_atom',mt_dbl,1,E)
      dbl_mb(fxcatom(1)+i-1) = E

      value = value.and.
     >        btdb_get(rtdb,'pspw:E_pxc_atom',mt_dbl,1,E)
      dbl_mb(fpxcatom(1)+i-1) = E

      value = value.and.
     >        btdb_get(rtdb,'pspw:E_all',mt_dbl,1,E)
      dbl_mb(fall(1)+i-1) = E

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


      value = value.and.btdb_get(rtdb,'pspw:E_ewald',mt_dbl,1,E)
      dbl_mb(few(1)+i-1) = dbl_mb(few(1)+i-1) - E

      value = value.and.btdb_get(rtdb,'pspw:E_local',mt_dbl,1,E)
      dbl_mb(fvl(1)+i-1) = dbl_mb(fvl(1)+i-1) - E

      value = value.and.btdb_get(rtdb,'pspw:E_nonlocal',mt_dbl,1,E)
      dbl_mb(fvnl(1)+i-1) = dbl_mb(fvnl(1)+i-1) - E

      value = value.and.btdb_get(rtdb,'pspw:E_xc',mt_dbl,1,E)
      dbl_mb(fxc(1)+i-1) = dbl_mb(fxc(1)+i-1) - E

      value = value.and.btdb_get(rtdb,'pspw:E_LJ',mt_dbl,1,E)
      dbl_mb(fLJ(1)+i-1) = dbl_mb(fLJ(1)+i-1) - E

      value = value.and.btdb_get(rtdb,'pspw:E_Q',mt_dbl,1,E)
      dbl_mb(fQ(1)+i-1) = dbl_mb(fQ(1)+i-1) - E

      value = value.and.btdb_get(rtdb,'pspw:E_spring',mt_dbl,1,E)
      dbl_mb(fspring(1)+i-1) = dbl_mb(fspring(1)+i-1) - E

      value = value.and.btdb_get(rtdb,'pspw:E_dftu',mt_dbl,1,E)
      dbl_mb(fdftu(1)+i-1) = dbl_mb(fdftu(1)+i-1) - E

      value = value.and.btdb_get(rtdb,'pspw:E_disp',mt_dbl,1,E)
      dbl_mb(fdisp(1)+i-1) = dbl_mb(fdisp(1)+i-1) - E

      value = value.and.btdb_get(rtdb,'pspw:kin_loc',mt_dbl,1,E)
      dbl_mb(fkloc(1)+i-1) = dbl_mb(fkloc(1)+i-1) - E

      value = value.and.btdb_get(rtdb,'pspw:kin_core',mt_dbl,1,E)
      dbl_mb(fkcore(1)+i-1) = dbl_mb(fkcore(1)+i-1) - E

      value = value.and.btdb_get(rtdb,'pspw:val_core',mt_dbl,1,E)
      dbl_mb(fvalcore(1)+i-1) = dbl_mb(fvalcore(1)+i-1) - E

      value = value.and.btdb_get(rtdb,'pspw:ion_core',mt_dbl,1,E)
      dbl_mb(fioncore(1)+i-1) = dbl_mb(fioncore(1)+i-1) - E

      value = value.and.btdb_get(rtdb,'pspw:E_ncmp_vloc',mt_dbl,1,E)
      dbl_mb(fcmpvl(1)+i-1) = dbl_mb(fcmpvl(1)+i-1) - E

      value = value.and.btdb_get(rtdb,'pspw:vloc_atom',mt_dbl,1,E)
      dbl_mb(fvlocatom(1)+i-1) = dbl_mb(fvlocatom(1)+i-1) - E

      value = value.and.btdb_get(rtdb,'pspw:E_hartree_atom',mt_dbl,1,E)
      dbl_mb(fhatom(1)+i-1) = dbl_mb(fhatom(1)+i-1) - E

      value = value.and.btdb_get(rtdb,'pspw:E_cmp_cmp',mt_dbl,1,E)
      dbl_mb(fcmp(1)+i-1) = dbl_mb(fcmp(1)+i-1) - E

      value = value.and.btdb_get(rtdb,'pspw:E_cmp_pw',mt_dbl,1,E)
      dbl_mb(fcmppw(1)+i-1) = dbl_mb(fcmppw(1)+i-1) - E

      value = value.and.btdb_get(rtdb,'pspw:E_exc_atom',mt_dbl,1,E)
      dbl_mb(fxcatom(1)+i-1) = dbl_mb(fxcatom(1)+i-1) - E

      value = value.and.btdb_get(rtdb,'pspw:E_pxc_atom',mt_dbl,1,E)
      dbl_mb(fpxcatom(1)+i-1) = dbl_mb(fpxcatom(1)+i-1) - E

      value = value.and.btdb_get(rtdb,'pspw:E_all',mt_dbl,1,E)
      dbl_mb(fall(1)+i-1) = dbl_mb(fall(1)+i-1) - E

      end do
      call dscal(3*nion,(-0.5d0/alpha),dbl_mb(few(1)),1)
      call dscal(3*nion,(-0.5d0/alpha),dbl_mb(fvl(1)),1)
      call dscal(3*nion,(-0.5d0/alpha),dbl_mb(fvnl(1)),1)
      call dscal(3*nion,(-0.5d0/alpha),dbl_mb(fxc(1)),1)
      call dscal(3*nion,(-0.5d0/alpha),dbl_mb(fLJ(1)),1)
      call dscal(3*nion,(-0.5d0/alpha),dbl_mb(fQ(1)),1)
      call dscal(3*nion,(-0.5d0/alpha),dbl_mb(fspring(1)),1)
      call dscal(3*nion,(-0.5d0/alpha),dbl_mb(fdftu(1)),1)
      call dscal(3*nion,(-0.5d0/alpha),dbl_mb(fdisp(1)),1)
      call dscal(3*nion,(-0.5d0/alpha),dbl_mb(fkloc(1)),1)
      call dscal(3*nion,(-0.5d0/alpha),dbl_mb(fkcore(1)),1)
      call dscal(3*nion,(-0.5d0/alpha),dbl_mb(fvalcore(1)),1)
      call dscal(3*nion,(-0.5d0/alpha),dbl_mb(fioncore(1)),1)
      call dscal(3*nion,(-0.5d0/alpha),dbl_mb(fcmpvl(1)),1)
      call dscal(3*nion,(-0.5d0/alpha),dbl_mb(fhatom(1)),1)
      call dscal(3*nion,(-0.5d0/alpha),dbl_mb(fvlocatom(1)),1)
      call dscal(3*nion,(-0.5d0/alpha),dbl_mb(fcmp(1)),1)
      call dscal(3*nion,(-0.5d0/alpha),dbl_mb(fcmppw(1)),1)
      call dscal(3*nion,(-0.5d0/alpha),dbl_mb(fxcatom(1)),1)
      call dscal(3*nion,(-0.5d0/alpha),dbl_mb(fpxcatom(1)),1)
      call dscal(3*nion,(-0.5d0/alpha),dbl_mb(fall(1)),1)

      call write_force(nion,dbl_mb(few(1)),"ewald FD")
      call write_force(nion,dbl_mb(fvl(1)),"local FD")
      call write_force(nion,dbl_mb(fvnl(1)),"nonlocal FD")
      call write_force(nion,dbl_mb(fxc(1)),"xc FD")
      call write_force(nion,dbl_mb(fLJ(1)),"LJ FD")
      call write_force(nion,dbl_mb(fQ(1)), "Q FD")
      call write_force(nion,dbl_mb(fspring(1)), "spring FD")
      call write_force(nion,dbl_mb(fdftu(1)), "DFTU FD")
      call write_force(nion,dbl_mb(fdisp(1)), "Dispersion FD")
      call write_force(nion,dbl_mb(fkloc(1)), "kinetic loc FD")
      call write_force(nion,dbl_mb(fkcore(1)), "kinetic core FD")
      call write_force(nion,dbl_mb(fvalcore(1)), "valence core FD")
      call write_force(nion,dbl_mb(fioncore(1)), "ion core FD")
      call write_force(nion,dbl_mb(fcmpvl(1)), "ncmpVloc FD")
      call write_force(nion,dbl_mb(fvlocatom(1)), "Vloc Atom FD")
      call write_force(nion,dbl_mb(fhatom(1)), "Hartree Atom FD")
      call write_force(nion,dbl_mb(fcmp(1)), "cmp-cmp FD")
      call write_force(nion,dbl_mb(fcmppw(1)), "cmp-pw FD")
      call write_force(nion,dbl_mb(fxcatom(1)), "exc_atom FD")
      call write_force(nion,dbl_mb(fpxcatom(1)), "pxc_atom FD")
      call write_force(nion,dbl_mb(fall(1)), "E all FD")

c      call daxpy(3*nion,1.0d0,dbl_mb(fvl(1)),1,dbl_mb(few(1)),1)
c      call daxpy(3*nion,1.0d0,dbl_mb(fvnl(1)),1,dbl_mb(few(1)),1)
c      call daxpy(3*nion,1.0d0,dbl_mb(fxc(1)),1,dbl_mb(few(1)),1)
c      call daxpy(3*nion,1.0d0,dbl_mb(fLJ(1)),1,dbl_mb(few(1)),1)
c      call daxpy(3*nion,1.0d0,dbl_mb(fQ(1)),1,dbl_mb(few(1)),1)
c      call daxpy(3*nion,1.0d0,dbl_mb(fspring(1)),1,dbl_mb(few(1)),1)
c      call daxpy(3*nion,1.0d0,dbl_mb(fdftu(1)),1,dbl_mb(few(1)),1)
c      call daxpy(3*nion,1.0d0,dbl_mb(fdisp(1)),1,dbl_mb(few(1)),1)
c
c      call daxpy(3*nion,-1.0d0,dbl_mb(fcmpvl(1)),1,dbl_mb(few(1)),1)
c      call daxpy(3*nion,-1.0d0,dbl_mb(fhatom(1)),1,dbl_mb(few(1)),1)
c      call daxpy(3*nion,-1.0d0,dbl_mb(fcmp(1)),1,dbl_mb(few(1)),1)
c      call daxpy(3*nion,-1.0d0,dbl_mb(fcmppw(1)),1,dbl_mb(few(1)),1)
c      call daxpy(3*nion,1.0d0,dbl_mb(fxcatom(1)),1,dbl_mb(few(1)),1)
c      call daxpy(3*nion,-1.0d0,dbl_mb(fpxcatom(1)),1,dbl_mb(few(1)),1)
c      call write_force(nion,dbl_mb(few(1)),"total FD")

c      call dcopy(3*nion,dbl_mb(fcmpvl(1)),      1,dbl_mb(few(1)),1)
c      call dcopy(3*nion,dbl_mb(fvlocatom(1)),   1,dbl_mb(few(1)),1)
c      call daxpy(3*nion,1.0d0,dbl_mb(fhatom(1)),1,dbl_mb(few(1)),1)
c      call daxpy(3*nion,1.0d0,dbl_mb(fcmp(1)),  1,dbl_mb(few(1)),1)
c      call daxpy(3*nion,1.0d0,dbl_mb(fcmppw(1)),1,dbl_mb(few(1)),1)
c      call daxpy(3*nion,1.0d0,dbl_mb(fxcatom(1)),1,dbl_mb(few(1)),1)

      !call dcopy(3*nion,0.0d0,0,dbl_mb(few(1)),1)
      call daxpy(3*nion,1.0d0,dbl_mb(fvl(1)),1,dbl_mb(few(1)),1)
      call daxpy(3*nion,1.0d0,dbl_mb(fkloc(1)),1,dbl_mb(few(1)),1)
      call daxpy(3*nion,1.0d0,dbl_mb(fkcore(1)),1,dbl_mb(few(1)),1)
      call daxpy(3*nion,1.0d0,dbl_mb(fvalcore(1)),1,dbl_mb(few(1)),1)
      call daxpy(3*nion,1.0d0,dbl_mb(fioncore(1)),1,dbl_mb(few(1)),1)
      call daxpy(3*nion,1.0d0,dbl_mb(fvlocatom(1)),1,dbl_mb(few(1)),1)
      call daxpy(3*nion,1.0d0,dbl_mb(fhatom(1)),1,dbl_mb(few(1)),1)
      call daxpy(3*nion,1.0d0,dbl_mb(fxcatom(1)),1,dbl_mb(few(1)),1)
      call daxpy(3*nion,1.0d0,dbl_mb(fcmp(1)),1,dbl_mb(few(1)),1)
      call daxpy(3*nion,1.0d0,dbl_mb(fcmppw(1)),1,dbl_mb(few(1)),1)
      call daxpy(3*nion,1.0d0,dbl_mb(fcmpvl(1)),1,dbl_mb(few(1)),1)

      call write_force(nion,dbl_mb(few(1)),
     > "TEST FD:kloc+vloc_at+core+h_at+xc_at+cmp+cmp_pw+cmpvl+vl+ewald")




      if (.not. geom_cart_coords_set(geom,dbl_mb(r1(1))))
     >   call errquit('force_fd: geom cart?',0,GEOM_ERR)
      !value = value.and.geom_rtdb_delete(rtdb,'geometry')
      value = value.and.geom_rtdb_store(rtdb,geom,'geometry')
      value = value.and.geom_destroy(geom)
      value = value.and.BA_free_heap(r1(2))
      value = value.and.BA_free_heap(r2(2))
      value = value.and.BA_free_heap(few(2))
      value = value.and.BA_free_heap(fvl(2))
      value = value.and.BA_free_heap(fvnl(2))
      value = value.and.BA_free_heap(fxc(2))
      value = value.and.BA_free_heap(fLJ(2))
      value = value.and.BA_free_heap(fQ(2))
      value = value.and.BA_free_heap(fspring(2))
      value = value.and.BA_free_heap(fdftu(2))
      value = value.and.BA_free_heap(fdisp(2))
      value = value.and.BA_free_heap(fkloc(2))
      value = value.and.BA_free_heap(fkcore(2))
      value = value.and.BA_free_heap(fvalcore(2))
      value = value.and.BA_free_heap(fioncore(2))
      value = value.and.BA_free_heap(fcmpvl(2))
      value = value.and.BA_free_heap(fvlocatom(2))
      value = value.and.BA_free_heap(fhatom(2))
      value = value.and.BA_free_heap(fcmp(2))
      value = value.and.BA_free_heap(fcmppw(2))
      value = value.and.BA_free_heap(fxcatom(2))
      value = value.and.BA_free_heap(fpxcatom(2))
      value = value.and.BA_free_heap(fall(2))
 
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
