*
* $Id$
*
      logical function stress_fd(rtdb)
      implicit none
      integer rtdb
      logical value

#include "bafdecls.fh"
#include "btdb.fh"
#include "inp.fh"
#include "util.fh"
#include "geom.fh"
#include "stdio.fh"

      logical oprint,lstress_numerical,noit_energy
      integer taskid,MASTER
      parameter (MASTER=0)

      integer l,i,j
      integer r1(2),geom,nion


      double precision alpha,E
      double precision unita(3,3)
      double precision unita2(3,3)
      double precision unita3(3,3)

      double precision stress(3,3)
      double precision ke_stress(3,3)
      double precision xc_stress(3,3)
      double precision ewald_stress(3,3)
      double precision hartree_stress(3,3)
      double precision local_stress(3,3)
      double precision nonlocal_stress(3,3)
      double precision alpha_stress(3,3)
      double precision RIJ_stress(3,3)
      double precision ewaldG_stress(3,3)
      double precision ehsic_stress(3,3)
      double precision excsic_stress(3,3)
      double precision hfx_stress(3,3)
      double precision LJ_stress(3,3)
      double precision Q_stress(3,3)
      double precision spring_stress(3,3)
      double precision dftu_stress(3,3)
      double precision disp_stress(3,3)

      double precision kloc_stress(3,3)
      double precision kcore_stress(3,3)
      double precision valcore_stress(3,3)
      double precision ioncore_stress(3,3)
      double precision cmpvl_stress(3,3)
      double precision vlocatom_stress(3,3)
      double precision hatom_stress(3,3)
      double precision cmp_stress(3,3)
      double precision cmppw_stress(3,3)
      double precision xcatom_stress(3,3)
      double precision pxcatom_stress(3,3)
      double precision all_stress(3,3)

      double precision lattice(6)
      double precision lattice2(6)
      double precision lattice3(6)

      double precision lstress(6)
      double precision ke_lstress(6)
      double precision xc_lstress(6)
      double precision ewald_lstress(6)
      double precision hartree_lstress(6)
      double precision local_lstress(6)
      double precision nonlocal_lstress(6)
      double precision alpha_lstress(6)
      double precision RIJ_lstress(6)
      double precision ewaldG_lstress(6)
      double precision ehsic_lstress(6)
      double precision excsic_lstress(6)
      double precision hfx_lstress(6)
      double precision LJ_lstress(6)
      double precision Q_lstress(6)
      double precision spring_lstress(6)
      double precision dftu_lstress(6)
      double precision disp_lstress(6)

      double precision kloc_lstress(6)
      double precision kcore_lstress(6)
      double precision valcore_lstress(6)
      double precision ioncore_lstress(6)
      double precision cmpvl_lstress(6)
      double precision vlocatom_lstress(6)
      double precision hatom_lstress(6)
      double precision cmp_lstress(6)
      double precision cmppw_lstress(6)
      double precision xcatom_lstress(6)
      double precision pxcatom_lstress(6)
      double precision all_lstress(6)

      integer date_len
      character*40 date_string

*     *** external functions ****
      logical  pspw_noit_energy,pspw_energy
      external pspw_noit_energy,pspw_energy

      call Parallel_taskid(taskid)
      oprint = (taskid.eq.MASTER)

*     **** initialize and read in geometry object ****
      value=.true.
      value     = value.and.geom_create(geom,'geometry')
      value     = value.and.geom_rtdb_load(rtdb,geom,'geometry')
      value     = value.and.geom_ncent(geom,nion)
      value     = value.and.
     >            BA_alloc_get(mt_dbl,(3*nion),'r1',r1(2),r1(1))
      if(.not.value) call errquit
     E   (' stressfd: geom malloc failure ',0,0)
      call load_cellgeometry(geom,dbl_mb(r1(1)),unita)

      if (.not.
     >    btdb_get(rtdb,'nwpw:lstress_numerical',
     >             mt_log,1,lstress_numerical))
     >  lstress_numerical = .false.

      if (.not.
     >    btdb_get(rtdb,'nwpw:noit_stress_numerical',
     >             mt_log,1,noit_energy))
     >  noit_energy = .true.


*     **** get finite difference step and initialize ****
      if (.not.
     >    btdb_get(rtdb,'pspw:stress_alpha',mt_dbl,1,alpha))
     >  alpha = 0.00100d0


*     **** print out finite difference stress options ****
      if (oprint) then
         write(luout,*)
         if (lstress_numerical) then
            write(luout,*) "numerical finite difference pspw lstress"
         else
            write(luout,*) "numerical finite difference pspw stress"
         end if
         if (noit_energy) then
            write(luout,*) "  - noit_energy"
         else
            write(luout,*) "  - energy mimimized"
         end if
         write(luout,*) "  - alpha=",alpha
         write(luout,*)
      end if




      if (.not.lstress_numerical) then
      call dcopy(9,0.0d0,0,ke_stress,      1)
      call dcopy(9,0.0d0,0,xc_stress,      1)
      call dcopy(9,0.0d0,0,hartree_stress, 1)
      call dcopy(9,0.0d0,0,ewald_stress,   1)
      call dcopy(9,0.0d0,0,local_stress,   1)
      call dcopy(9,0.0d0,0,nonlocal_stress,1)
      call dcopy(9,0.0d0,0,ehsic_stress,1)
      call dcopy(9,0.0d0,0,excsic_stress,1)
      call dcopy(9,0.0d0,0,hfx_stress,1)
      call dcopy(9,0.0d0,0,LJ_stress,1)
      call dcopy(9,0.0d0,0,Q_stress,1)
      call dcopy(9,0.0d0,0,spring_stress,1)
      call dcopy(9,0.0d0,0,dftu_stress,1)
      call dcopy(9,0.0d0,0,disp_stress,1)

      call dcopy(9,0.0d0,0,kloc_stress,1)
      call dcopy(9,0.0d0,0,kcore_stress,1)
      call dcopy(9,0.0d0,0,valcore_stress,1)
      call dcopy(9,0.0d0,0,ioncore_stress,1)
      call dcopy(9,0.0d0,0,cmpvl_stress,1)
      call dcopy(9,0.0d0,0,vlocatom_stress,1)
      call dcopy(9,0.0d0,0,hatom_stress,1)
      call dcopy(9,0.0d0,0,cmp_stress,1)
      call dcopy(9,0.0d0,0,cmppw_stress,1)
      call dcopy(9,0.0d0,0,xcatom_stress,1)
      call dcopy(9,0.0d0,0,pxcatom_stress,1)
      call dcopy(9,0.0d0,0,all_stress,1)


      do j=1,3
      do i=1,3


*     *****************************
*     **** + finite difference ****
*     *****************************
      call dcopy(9,unita,1,unita2,1)
      unita2(i,j) = unita2(i,j) + alpha

*     **** update geometry object ****
      call util_date(date_string)
      date_len = inp_strlen(date_string)
      if (oprint) write(luout,*)
      if (oprint) write(luout,
     &   '(1x,a12,i2,i2,1x,f8.3,1x,a4,1x,a10,f10.1,6x,a5,a27)')
     &   'amatrix i,j:',i,j,unita2(i,j),'(+)',
     &   'wall time:',util_wallsec(),
     &   'date:',date_string(1:date_len)
      if (oprint) call util_flush(6)
      call store_cellgeometry(geom,dbl_mb(r1(1)),unita2)
      value = value.and.geom_rtdb_store(rtdb,geom,'geometry')
      if (noit_energy) then
         value = value.and.pspw_noit_energy(rtdb)
      else
         value = value.and.pspw_energy(rtdb)
      end if


      value = value.and.
     >        btdb_get(rtdb,'pspw:E_ke',mt_dbl,1,E)
      ke_stress(i,j) = E

      value = value.and.
     >        btdb_get(rtdb,'pspw:E_xc',mt_dbl,1,E)
      xc_stress(i,j) = E

      value = value.and.
     >        btdb_get(rtdb,'pspw:E_hartree',mt_dbl,1,E)
      hartree_stress(i,j) = E

      value = value.and.
     >        btdb_get(rtdb,'pspw:E_ewald',mt_dbl,1,E)
      ewald_stress(i,j) = E

      value = value.and.
     >        btdb_get(rtdb,'pspw:E_local',mt_dbl,1,E)
      local_stress(i,j) = E

      value = value.and.
     >        btdb_get(rtdb,'pspw:E_nonlocal',mt_dbl,1,E)
      nonlocal_stress(i,j) = E

      value = value.and.
     >        btdb_get(rtdb,'pspw:E_hsic',mt_dbl,1,E)
      ehsic_stress(i,j) = E

      value = value.and.
     >        btdb_get(rtdb,'pspw:E_xcsic',mt_dbl,1,E)
      excsic_stress(i,j) = E

      value = value.and.
     >        btdb_get(rtdb,'pspw:E_hfx',mt_dbl,1,E)
      hfx_stress(i,j) = E

      value = value.and.
     >        btdb_get(rtdb,'pspw:E_LJ',mt_dbl,1,E)
      LJ_stress(i,j) = E

      value = value.and.
     >        btdb_get(rtdb,'pspw:E_Q',mt_dbl,1,E)
      Q_stress(i,j) = E

      value = value.and.
     >        btdb_get(rtdb,'pspw:E_spring',mt_dbl,1,E)
      spring_stress(i,j) = E

      value = value.and.
     >        btdb_get(rtdb,'pspw:E_dftu',mt_dbl,1,E)
      dftu_stress(i,j) = E

      value = value.and.
     >        btdb_get(rtdb,'pspw:E_disp',mt_dbl,1,E)
      disp_stress(i,j) = E

      value = value.and.
     >        btdb_get(rtdb,'pspw:kin_loc',mt_dbl,1,E)
      kloc_stress(i,j) = E

      value = value.and.
     >        btdb_get(rtdb,'pspw:kin_core',mt_dbl,1,E)
      kcore_stress(i,j) = E

      value = value.and.
     >        btdb_get(rtdb,'pspw:val_core',mt_dbl,1,E)
      valcore_stress(i,j) = E

      value = value.and.
     >        btdb_get(rtdb,'pspw:ion_core',mt_dbl,1,E)
      ioncore_stress(i,j) = E

      value = value.and.
     >        btdb_get(rtdb,'pspw:E_ncmp_vloc',mt_dbl,1,E)
      cmpvl_stress(i,j) = E

      value = value.and.
     >        btdb_get(rtdb,'pspw:vloc_atom',mt_dbl,1,E)
      vlocatom_stress(i,j) = E

      value = value.and.
     >        btdb_get(rtdb,'pspw:E_hartree_atom',mt_dbl,1,E)
      hatom_stress(i,j) = E

      value = value.and.
     >        btdb_get(rtdb,'pspw:E_cmp_cmp',mt_dbl,1,E)
      cmp_stress(i,j) = E

      value = value.and.
     >        btdb_get(rtdb,'pspw:E_cmp_pw',mt_dbl,1,E)
      cmppw_stress(i,j) = E

      value = value.and.
     >        btdb_get(rtdb,'pspw:E_exc_atom',mt_dbl,1,E)
      xcatom_stress(i,j) = E

      value = value.and.
     >        btdb_get(rtdb,'pspw:E_pxc_atom',mt_dbl,1,E)
      pxcatom_stress(i,j) = E

      value = value.and.
     >        btdb_get(rtdb,'pspw:E_all',mt_dbl,1,E)
      all_stress(i,j) = E




*     *****************************
*     **** - finite difference ****
*     *****************************
      call dcopy(9,unita,1,unita3,1)
      unita3(i,j) = unita3(i,j) - alpha

*     **** update geometry object ****
      call util_date(date_string)
      date_len = inp_strlen(date_string)
      if (oprint) write(luout,*)
      if (oprint) write(luout,
     &   '(1x,a12,i2,i2,1x,f8.3,1x,a4,1x,a10,f10.1,6x,a5,a27)')
     &   'amatrix i,j:',i,j,unita3(i,j),'(-)',
     &   'wall time:',util_wallsec(),
     &   'date:',date_string(1:date_len)
      if (oprint) call util_flush(6)
      call store_cellgeometry(geom,dbl_mb(r1(1)),unita3)
      value = value.and.geom_rtdb_store(rtdb,geom,'geometry')
      if (noit_energy) then
         value = value.and.pspw_noit_energy(rtdb)
      else
         value = value.and.pspw_energy(rtdb)
      end if


      value = value.and.btdb_get(rtdb,'pspw:E_ke',mt_dbl,1,E)
      ke_stress(i,j) = ke_stress(i,j) - E

      value = value.and.btdb_get(rtdb,'pspw:E_xc',mt_dbl,1,E)
      xc_stress(i,j) = xc_stress(i,j) - E

      value = value.and.btdb_get(rtdb,'pspw:E_hartree',mt_dbl,1,E)
      hartree_stress(i,j) = hartree_stress(i,j) - E

      value = value.and.btdb_get(rtdb,'pspw:E_ewald',mt_dbl,1,E)
      ewald_stress(i,j) = ewald_stress(i,j) - E

      value = value.and.btdb_get(rtdb,'pspw:E_local',mt_dbl,1,E)
      local_stress(i,j) = local_stress(i,j) - E

      value = value.and.btdb_get(rtdb,'pspw:E_nonlocal',mt_dbl,1,E)
      nonlocal_stress(i,j) = nonlocal_stress(i,j) - E

      value = value.and.btdb_get(rtdb,'pspw:E_hsic',mt_dbl,1,E)
      ehsic_stress(i,j) = ehsic_stress(i,j) - E

      value = value.and.btdb_get(rtdb,'pspw:E_xcsic',mt_dbl,1,E)
      excsic_stress(i,j) = excsic_stress(i,j) - E

      value = value.and.btdb_get(rtdb,'pspw:E_hfx',mt_dbl,1,E)
      hfx_stress(i,j) = hfx_stress(i,j) - E

      value = value.and.btdb_get(rtdb,'pspw:E_LJ',mt_dbl,1,E)
      LJ_stress(i,j) = LJ_stress(i,j) - E

      value = value.and.btdb_get(rtdb,'pspw:E_Q',mt_dbl,1,E)
      Q_stress(i,j) = Q_stress(i,j) - E

      value = value.and.btdb_get(rtdb,'pspw:E_spring',mt_dbl,1,E)
      spring_stress(i,j) = spring_stress(i,j) - E

      value = value.and.btdb_get(rtdb,'pspw:E_dftu',mt_dbl,1,E)
      dftu_stress(i,j) = dftu_stress(i,j) - E

      value = value.and.btdb_get(rtdb,'pspw:E_disp',mt_dbl,1,E)
      disp_stress(i,j) = disp_stress(i,j) - E

      value = value.and.btdb_get(rtdb,'pspw:kin_loc',mt_dbl,1,E)
      kloc_stress(i,j) = kloc_stress(i,j) - E

      value = value.and.btdb_get(rtdb,'pspw:kin_core',mt_dbl,1,E)
      kcore_stress(i,j) = kcore_stress(i,j) - E

      value = value.and.btdb_get(rtdb,'pspw:val_core',mt_dbl,1,E)
      valcore_stress(i,j) = valcore_stress(i,j) - E

      value = value.and.btdb_get(rtdb,'pspw:ion_core',mt_dbl,1,E)
      ioncore_stress(i,j) = ioncore_stress(i,j) - E

      value = value.and.btdb_get(rtdb,'pspw:E_ncmp_vloc',mt_dbl,1,E)
      cmpvl_stress(i,j) = cmpvl_stress(i,j) - E

      value = value.and.btdb_get(rtdb,'pspw:vloc_atom',mt_dbl,1,E)
      vlocatom_stress(i,j) = vlocatom_stress(i,j) - E

      value = value.and.btdb_get(rtdb,'pspw:E_hartree_atom',mt_dbl,1,E)
      hatom_stress(i,j) = hatom_stress(i,j) - E

      value = value.and.btdb_get(rtdb,'pspw:E_cmp_cmp',mt_dbl,1,E)
      cmp_stress(i,j) = cmp_stress(i,j) - E

      value = value.and.btdb_get(rtdb,'pspw:E_cmp_pw',mt_dbl,1,E)
      cmppw_stress(i,j) = cmppw_stress(i,j) - E

      value = value.and.btdb_get(rtdb,'pspw:E_exc_atom',mt_dbl,1,E)
      xcatom_stress(i,j) = xcatom_stress(i,j) - E

      value = value.and.btdb_get(rtdb,'pspw:E_pxc_atom',mt_dbl,1,E)
      pxcatom_stress(i,j) = pxcatom_stress(i,j) - E

      value = value.and.btdb_get(rtdb,'pspw:E_all',mt_dbl,1,E)
      all_stress(i,j) = all_stress(i,j) - E

      end do
      end do
      call dscal(9,(0.5d0/alpha),ke_stress,1)
      call dscal(9,(0.5d0/alpha),xc_stress,1)
      call dscal(9,(0.5d0/alpha),hartree_stress,1)
      call dscal(9,(0.5d0/alpha),ewald_stress,1)
      call dscal(9,(0.5d0/alpha),local_stress,1)
      call dscal(9,(0.5d0/alpha),nonlocal_stress,1)
      call dscal(9,(0.5d0/alpha),ehsic_stress,1)
      call dscal(9,(0.5d0/alpha),excsic_stress,1)
      call dscal(9,(0.5d0/alpha),hfx_stress,1)
      call dscal(9,(0.5d0/alpha),LJ_stress,1)
      call dscal(9,(0.5d0/alpha),Q_stress,1)
      call dscal(9,(0.5d0/alpha),spring_stress,1)
      call dscal(9,(0.5d0/alpha),dftu_stress,1)
      call dscal(9,(0.5d0/alpha),disp_stress,1)

      call dscal(9,(0.5d0/alpha),kloc_stress,1)
      call dscal(9,(0.5d0/alpha),kcore_stress,1)
      call dscal(9,(0.5d0/alpha),valcore_stress,1)
      call dscal(9,(0.5d0/alpha),ioncore_stress,1)
      call dscal(9,(0.5d0/alpha),cmpvl_stress,1)
      call dscal(9,(0.5d0/alpha),vlocatom_stress,1)
      call dscal(9,(0.5d0/alpha),hatom_stress,1)
      call dscal(9,(0.5d0/alpha),cmp_stress,1)
      call dscal(9,(0.5d0/alpha),cmppw_stress,1)
      call dscal(9,(0.5d0/alpha),xcatom_stress,1)
      call dscal(9,(0.5d0/alpha),pxcatom_stress,1)
      call dscal(9,(0.5d0/alpha),all_stress,1)

      call write_stress(ke_stress,"ke FD")
      call write_stress(hartree_stress,"coulomb FD")
      call write_stress(local_stress,"local FD")
      call write_stress(nonlocal_stress,"nonlocal FD")
      call write_stress(xc_stress,"xc FD")
      call write_stress(ewald_stress,"ewald FD")
      call write_stress(ehsic_stress,"Eh SIC FD")
      call write_stress(excsic_stress,"Exc SIC FD")
      call write_stress(hfx_stress,"HFX FD")
      call write_stress(LJ_stress,"LJ FD")
      call write_stress(Q_stress,"Q FD")
      call write_stress(spring_stress,"spring FD")
      call write_stress(dftu_stress,"DFTU FD")
      call write_stress(disp_stress,"Dispersion FD")

      call write_stress(kloc_stress,"kinetic loc FD")
      call write_stress(kcore_stress,"kinetic core FD")
      call write_stress(valcore_stress,"valence core FD")
      call write_stress(ioncore_stress, "ion core FD")
      call write_stress(cmpvl_stress, "ncmpVloc FD")
      call write_stress(vlocatom_stress,"Vloc Atom FD")
      call write_stress(hatom_stress, "Hartree atom FD")
      call write_stress(cmp_stress, "cmp-cmp FD")
      call write_stress(cmppw_stress,"cmp-pw FD")
      call write_stress(xcatom_stress,"xc_atom FD")
      call write_stress(pxcatom_stress,"pxc_atom FD")
      call write_stress(all_stress, "E all FD")

      !call dcopy(9,ke_stress,1,stress,1)
      !call daxpy(9,1.0d0,xc_stress,1,stress,1)
      !call daxpy(9,1.0d0,hartree_stress,1,stress,1)
      !call daxpy(9,1.0d0,ewald_stress,1,stress,1)
      !call daxpy(9,1.0d0,local_stress,1,stress,1)
      !call daxpy(9,1.0d0,nonlocal_stress,1,stress,1)
      !call daxpy(9,1.0d0,ehsic_stress,1,stress,1)
      !call daxpy(9,1.0d0,excsic_stress,1,stress,1)
      !call daxpy(9,1.0d0,hfx_stress,1,stress,1)
      !call daxpy(9,1.0d0,LJ_stress,1,stress,1)
      !call daxpy(9,1.0d0,Q_stress,1,stress,1)
      !call daxpy(9,1.0d0,spring_stress,1,stress,1)
      !call daxpy(9,1.0d0,dftu_stress,1,stress,1)
      !call daxpy(9,1.0d0,disp_stress,1,stress,1)
      !call write_stress(stress,"TOTAL FD")
      call dcopy(9,all_stress,1,stress,1)

*     **************************************
*     **** output stress for use driver ****
*     **************************************
      value = value.and.btdb_put(rtdb,'pspw:stress',mt_dbl,9,stress)
      
      else
         call cell_unita_abc_abg(unita,lattice)

         do i=1,6
*           *****************************
*           **** + finite difference ****
*           *****************************
            call dcopy(6,lattice,1,lattice2,1)
            lattice2(i) = lattice2(i) + alpha
            call cell_abc_abg_unita(lattice2,unita2)

*           **** update geometry object ****
            call util_date(date_string)
            date_len = inp_strlen(date_string)
            if (oprint) write(luout,*)
            if (oprint) write(luout,
     &      '(1x,a12,i2,1x,f8.3,1x,a4,1x,a10,f10.1,6x,a5,a27)')
     &      'lattice i:',i,lattice2(i),'(+)',
     &      'wall time:',util_wallsec(),
     &      'date:',date_string(1:date_len)
            if (oprint) call util_flush(6)
            call store_cellgeometry(geom,dbl_mb(r1(1)),unita2)
            value = value.and.geom_rtdb_store(rtdb,geom,'geometry')
            if (noit_energy) then
               value = value.and.pspw_noit_energy(rtdb)
            else
               value = value.and.pspw_energy(rtdb)
            end if

            value = value.and.
     >              btdb_get(rtdb,'pspw:E_ke',mt_dbl,1,E)
            ke_lstress(i) = E
            value = value.and.
     >              btdb_get(rtdb,'pspw:E_xc',mt_dbl,1,E)
            xc_lstress(i) = E
            value = value.and.
     >              btdb_get(rtdb,'pspw:E_hartree',mt_dbl,1,E)
            hartree_lstress(i) = E
            value = value.and.
     >              btdb_get(rtdb,'pspw:E_ewald',mt_dbl,1,E)
            ewald_lstress(i) = E
            value = value.and.
     >              btdb_get(rtdb,'pspw:E_local',mt_dbl,1,E)
            local_lstress(i) = E
            value = value.and.
     >              btdb_get(rtdb,'pspw:E_nonlocal',mt_dbl,1,E)
            nonlocal_lstress(i) = E
            value = value.and.
     >              btdb_get(rtdb,'pspw:E_hsic',mt_dbl,1,E)
            ehsic_lstress(i) = E
            value = value.and.
     >              btdb_get(rtdb,'pspw:E_xcsic',mt_dbl,1,E)
            excsic_lstress(i) = E
            value = value.and.
     >              btdb_get(rtdb,'pspw:E_hfx',mt_dbl,1,E)
            hfx_lstress(i) = E
            value = value.and.
     >              btdb_get(rtdb,'pspw:E_LJ',mt_dbl,1,E)
            LJ_lstress(i) = E
            value = value.and.
     >              btdb_get(rtdb,'pspw:E_Q',mt_dbl,1,E)
            Q_lstress(i) = E
            value = value.and.
     >              btdb_get(rtdb,'pspw:E_spring',mt_dbl,1,E)
            spring_lstress(i) = E
            value = value.and.
     >              btdb_get(rtdb,'pspw:E_dftu',mt_dbl,1,E)
            dftu_lstress(i) = E
            value = value.and.
     >              btdb_get(rtdb,'pspw:E_disp',mt_dbl,1,E)
            disp_lstress(i) = E

            value = value.and.btdb_get(rtdb,'pspw:kin_loc',mt_dbl,1,E)
            kloc_lstress(i) = E

            value = value.and.btdb_get(rtdb,'pspw:kin_core',mt_dbl,1,E)
            kcore_lstress(i) = E

            value = value.and.btdb_get(rtdb,'pspw:val_core',mt_dbl,1,E)
            valcore_lstress(i) = E

            value = value.and.btdb_get(rtdb,'pspw:ion_core',mt_dbl,1,E)
            ioncore_lstress(i) = E

            value = value.and.
     >              btdb_get(rtdb,'pspw:E_ncmp_vloc',mt_dbl,1,E)
            cmpvl_lstress(i) = E

            value = value.and.btdb_get(rtdb,'pspw:vloc_atom',mt_dbl,1,E)
            vlocatom_lstress(i) = E

            value = value.and.
     >              btdb_get(rtdb,'pspw:E_hartree_atom',mt_dbl,1,E)
            hatom_lstress(i) = E

            value = value.and.btdb_get(rtdb,'pspw:E_cmp_cmp',mt_dbl,1,E)
            cmp_lstress(i) = E

            value = value.and.btdb_get(rtdb,'pspw:E_cmp_pw',mt_dbl,1,E)
            cmppw_lstress(i) = E

            value = value.and.
     >              btdb_get(rtdb,'pspw:E_exc_atom',mt_dbl,1,E)
            xcatom_lstress(i) = E

            value = value.and.
     >              btdb_get(rtdb,'pspw:E_pxc_atom',mt_dbl,1,E)
            pxcatom_lstress(i) = E

            value = value.and.btdb_get(rtdb,'pspw:E_all',mt_dbl,1,E)
            all_lstress(i) = E


*           *****************************
*           **** - finite difference ****
*           *****************************
            call dcopy(6,lattice,1,lattice3,1)
            lattice3(i) = lattice3(i) - alpha
            call cell_abc_abg_unita(lattice3,unita3)

*           **** update geometry object ****
            call util_date(date_string)
            date_len = inp_strlen(date_string)
            if (oprint) write(luout,*)
            if (oprint) write(luout,
     &         '(1x,a12,i2,1x,f8.3,1x,a4,1x,a10,f10.1,6x,a5,a27)')
     &         'lattice i:',i,lattice3(i),'(-)',
     &         'wall time:',util_wallsec(),
     &         'date:',date_string(1:date_len)
            if (oprint) call util_flush(6)
            call store_cellgeometry(geom,dbl_mb(r1(1)),unita3)
            value = value.and.geom_rtdb_store(rtdb,geom,'geometry')
            if (noit_energy) then
               value = value.and.pspw_noit_energy(rtdb)
            else
               value = value.and.pspw_energy(rtdb)
            end if

            value = value.and.btdb_get(rtdb,'pspw:E_ke',mt_dbl,1,E)
            ke_lstress(i) = ke_lstress(i) - E

            value = value.and.btdb_get(rtdb,'pspw:E_xc',mt_dbl,1,E)
            xc_lstress(i) = xc_lstress(i) - E

            value = value.and.btdb_get(rtdb,'pspw:E_hartree',mt_dbl,1,E)
            hartree_lstress(i) = hartree_lstress(i) - E

            value = value.and.btdb_get(rtdb,'pspw:E_ewald',mt_dbl,1,E)
            ewald_lstress(i) = ewald_lstress(i) - E

            value = value.and.btdb_get(rtdb,'pspw:E_local',mt_dbl,1,E)
            local_lstress(i) = local_lstress(i) - E

            value = value.and.
     >              btdb_get(rtdb,'pspw:E_nonlocal',mt_dbl,1,E)
            nonlocal_lstress(i) = nonlocal_lstress(i) - E

            value = value.and.btdb_get(rtdb,'pspw:E_hsic',mt_dbl,1,E)
            ehsic_lstress(i) = ehsic_lstress(i) - E

            value = value.and.btdb_get(rtdb,'pspw:E_xcsic',mt_dbl,1,E)
            excsic_lstress(i) = excsic_lstress(i) - E

            value = value.and.btdb_get(rtdb,'pspw:E_hfx',mt_dbl,1,E)
            hfx_lstress(i) = hfx_lstress(i) - E

            value = value.and.btdb_get(rtdb,'pspw:E_LJ',mt_dbl,1,E)
            LJ_lstress(i) = LJ_lstress(i) - E

            value = value.and.btdb_get(rtdb,'pspw:E_Q',mt_dbl,1,E)
            Q_lstress(i) = Q_lstress(i) - E

            value = value.and.btdb_get(rtdb,'pspw:E_spring',mt_dbl,1,E)
            spring_lstress(i) = spring_lstress(i) - E

            value = value.and.btdb_get(rtdb,'pspw:E_dftu',mt_dbl,1,E)
            dftu_lstress(i) = dftu_lstress(i) - E

            value = value.and.btdb_get(rtdb,'pspw:E_disp',mt_dbl,1,E)
            disp_lstress(i) = disp_lstress(i) - E

            value = value.and.btdb_get(rtdb,'pspw:kin_loc',mt_dbl,1,E)
            kloc_lstress(i) = kloc_lstress(i) - E

            value = value.and.btdb_get(rtdb,'pspw:kin_core',mt_dbl,1,E)
            kcore_lstress(i) = kcore_lstress(i) - E

            value = value.and.btdb_get(rtdb,'pspw:val_core',mt_dbl,1,E)
            valcore_lstress(i) = valcore_lstress(i) - E

            value = value.and.btdb_get(rtdb,'pspw:ion_core',mt_dbl,1,E)
            ioncore_lstress(i) = ioncore_lstress(i) - E

            value = value.and.
     >              btdb_get(rtdb,'pspw:E_ncmp_vloc',mt_dbl,1,E)
            cmpvl_lstress(i) = cmpvl_lstress(i) - E

            value = value.and.btdb_get(rtdb,'pspw:vloc_atom',mt_dbl,1,E)
            vlocatom_lstress(i) = vlocatom_lstress(i) - E

            value = value.and.
     >              btdb_get(rtdb,'pspw:E_hartree_atom',mt_dbl,1,E)
            hatom_lstress(i) = hatom_lstress(i) - E

            value = value.and.btdb_get(rtdb,'pspw:E_cmp_cmp',mt_dbl,1,E)
            cmp_lstress(i) = cmp_lstress(i) - E

            value = value.and.btdb_get(rtdb,'pspw:E_cmp_pw',mt_dbl,1,E)
            cmppw_lstress(i) = cmppw_lstress(i) - E

            value = value.and.
     >              btdb_get(rtdb,'pspw:E_exc_atom',mt_dbl,1,E)
            xcatom_lstress(i) = xcatom_lstress(i) - E

            value = value.and.
     >              btdb_get(rtdb,'pspw:E_pxc_atom',mt_dbl,1,E)
            pxcatom_lstress(i) = pxcatom_lstress(i) - E

            value = value.and.btdb_get(rtdb,'pspw:E_all',mt_dbl,1,E)
            all_lstress(i) = all_lstress(i) - E

         end do

         call dscal(6,(0.5d0/alpha),ke_lstress,1)
         call dscal(6,(0.5d0/alpha),xc_lstress,1)
         call dscal(6,(0.5d0/alpha),hartree_lstress,1)
         call dscal(6,(0.5d0/alpha),ewald_lstress,1)
         call dscal(6,(0.5d0/alpha),local_lstress,1)
         call dscal(6,(0.5d0/alpha),nonlocal_lstress,1)
         call dscal(6,(0.5d0/alpha),ehsic_lstress,1)
         call dscal(6,(0.5d0/alpha),excsic_lstress,1)
         call dscal(6,(0.5d0/alpha),hfx_lstress,1)
         call dscal(6,(0.5d0/alpha),LJ_lstress,1)
         call dscal(6,(0.5d0/alpha),Q_lstress,1)
         call dscal(6,(0.5d0/alpha),spring_lstress,1)
         call dscal(6,(0.5d0/alpha),dftu_lstress,1)
         call dscal(6,(0.5d0/alpha),disp_lstress,1)

         call dscal(6,(0.5d0/alpha),kloc_lstress,1)
         call dscal(6,(0.5d0/alpha),kcore_lstress,1)
         call dscal(6,(0.5d0/alpha),valcore_lstress,1)
         call dscal(6,(0.5d0/alpha),ioncore_lstress,1)
         call dscal(6,(0.5d0/alpha),cmpvl_lstress,1)
         call dscal(6,(0.5d0/alpha),vlocatom_lstress,1)
         call dscal(6,(0.5d0/alpha),hatom_lstress,1)
         call dscal(6,(0.5d0/alpha),cmp_lstress,1)
         call dscal(6,(0.5d0/alpha),cmppw_lstress,1)
         call dscal(6,(0.5d0/alpha),xcatom_lstress,1)
         call dscal(6,(0.5d0/alpha),pxcatom_lstress,1)
         call dscal(6,(0.5d0/alpha),all_lstress,1)

         call write_lstress(ke_lstress,"ke FD")
         call write_lstress(hartree_lstress,"coulomb FD")
         call write_lstress(local_lstress,"local FD")
         call write_lstress(nonlocal_lstress,"nonlocal FD")
         call write_lstress(xc_lstress,"xc FD")
         call write_lstress(ewald_lstress,"ewald FD")
         call write_lstress(ehsic_lstress,"Eh SIC FD")
         call write_lstress(excsic_lstress,"Exc SIC FD")
         call write_lstress(hfx_lstress,"HFX FD")
         call write_lstress(LJ_lstress,"LJ FD")
         call write_lstress(Q_lstress,"Q FD")
         call write_lstress(spring_lstress,"spring FD")
         call write_lstress(dftu_lstress,"dftu FD")
         call write_lstress(disp_lstress,"Dispersion FD")

         call write_lstress(kloc_lstress,"kinetic loc FD")
         call write_lstress(kcore_lstress,"kinetic core FD")
         call write_lstress(valcore_lstress,"valence core FD")
         call write_lstress(ioncore_lstress, "ion core FD")
         call write_lstress(cmpvl_lstress, "ncmpVloc FD")
         call write_lstress(vlocatom_lstress,"Vloc Atom FD")
         call write_lstress(hatom_lstress, "Hartree atom FD")
         call write_lstress(cmp_lstress, "cmp-cmp FD")
         call write_lstress(cmppw_lstress,"cmp-pw FD")
         call write_lstress(xcatom_lstress,"xc_atom FD")
         call write_lstress(pxcatom_lstress,"pxc_atom FD")
         call write_lstress(all_lstress, "E all FD")

         !call dcopy(6,ke_lstress,1,lstress,1)
         !call daxpy(6,1.0d0,xc_lstress,1,lstress,1)
         !call daxpy(6,1.0d0,hartree_lstress,1,lstress,1)
         !call daxpy(6,1.0d0,ewald_lstress,1,lstress,1)
         !call daxpy(6,1.0d0,local_lstress,1,lstress,1)
         !call daxpy(6,1.0d0,nonlocal_lstress,1,lstress,1)
         !call daxpy(6,1.0d0,ehsic_lstress,1,lstress,1)
         !call daxpy(6,1.0d0,excsic_lstress,1,lstress,1)
         !call daxpy(6,1.0d0,hfx_lstress,1,lstress,1)
         !call daxpy(6,1.0d0,LJ_lstress,1,lstress,1)
         !call daxpy(6,1.0d0,Q_lstress,1,lstress,1)
         !call daxpy(6,1.0d0,spring_lstress,1,lstress,1)
         !call daxpy(6,1.0d0,dftu_lstress,1,lstress,1)
         !call daxpy(6,1.0d0,disp_lstress,1,lstress,1)
         !call write_lstress(lstress,"TOTAL FD")
         call dcopy(6,all_lstress,1,lstress,1)

*        **************************************
*        **** output lstress for use driver ****
*        **************************************
         value = value.and.
     >           btdb_put(rtdb,'pspw:lstress',mt_dbl,6,lstress)
      end if

      call store_cellgeometry(geom,dbl_mb(r1(1)),unita)
      !value = value.and.geom_rtdb_delete(rtdb,'geometry')
      value = value.and.geom_rtdb_store(rtdb,geom,'geometry')
      value = value.and.geom_destroy(geom)
      value = value.and.BA_free_heap(r1(2))
 
      stress_fd = value
      return 
      end


