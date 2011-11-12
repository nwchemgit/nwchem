c
c $Id: pspw_Efield_grad.F 19707 2010-10-29 17:59:36Z d3y133 $
c


********************************************************************
*
*	pspw_gen_Efield_grad module: used to generate electric gradients fields at atoms
*
*
*	Data structure
*
********************************************************************

*     ***********************************
*     *					*
*     *		pspw_gen_Efield_grad	*
*     *					*
*     ***********************************
      subroutine pspw_gen_Efield_grad(rtdb,ispin,ne,psi,dng)
      implicit none
      integer rtdb
      integer ispin,ne(2)
      complex*16 psi(*)
      complex*16 dng(*)

#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"

*     **** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)

      logical value,doefield
      integer nion
      integer efg(2),efg_ion(2),efg_smooth(2),efg_cor(2)

*     **** external functions ****
      integer  ion_nion,control_version
      external ion_nion,control_version

      if (.not.rtdb_get(rtdb,'nwpw:efield',mt_log,1,doefield)) 
     >  doefield=.true.

      if (doefield) then
       
       call Parallel_taskid(taskid)
       nion = ion_nion()
       value =           MA_push_get(mt_dbl,9*nion,'efg',
     >                               efg(2),efg(1))
       value = value.and.MA_push_get(mt_dbl,9*nion,'efg_ion',
     >                               efg_ion(2),efg_ion(1))
       value = value.and.MA_push_get(mt_dbl,9*nion,'efg_smooth',
     >                               efg_smooth(2),efg_smooth(1))
       value = value.and.MA_push_get(mt_dbl,9*nion,'efg_cor',
     >                               efg_cor(2),efg_cor(1))
       if (.not.value) 
     >    call errquit('pspw_gen_Efield_grad:out of stack memory',0,
     >                 MA_ERR)

       call dcopy(9*nion,0.0d0,0,dbl_mb(efg_ion(1)),1)
       call dcopy(9*nion,0.0d0,0,dbl_mb(efg_smooth(1)),1)
       call dcopy(9*nion,0.0d0,0,dbl_mb(efg_cor(1)),1)
       call dcopy(9*nion,0.0d0,0,dbl_mb(efg(1)),1)

       call ewald_efg(dbl_mb(efg_ion(1)))
       call coulomb_efg(dng,dbl_mb(efg_smooth(1)))
       call psp_efg_atoms(ispin,ne,psi,dbl_mb(efg_cor(1)))

       call daxpy(9*nion,1.0d0,dbl_mb(efg_ion(1)),1,dbl_mb(efg(1)),1)
       call daxpy(9*nion,1.0d0,dbl_mb(efg_smooth(1)),1,dbl_mb(efg(1)),1)
       call daxpy(9*nion,1.0d0,dbl_mb(efg_cor(1)),1,dbl_mb(efg(1)),1)


       if (taskid.eq.MASTER) then
         call pspw_print_Efield_grad(luout,
     >                               nion,
     >                               dbl_mb(efg(1)),
     >                               dbl_mb(efg_ion(1)),
     >                               dbl_mb(efg_smooth(1)),
     >                               dbl_mb(efg_cor(1)))
       end if
       value =           MA_pop_stack(efg_cor(2))
       value = value.and.MA_pop_stack(efg_smooth(2))
       value = value.and.MA_pop_stack(efg_ion(2))
       value = value.and.MA_pop_stack(efg(2))
       if (.not.value) 
     >    call errquit('pspw_gen_Efield_grad:popping of stack memory',0,
     >                 MA_ERR)

      end if

      return
      end

*     ***********************************
*     *					*
*     *	    pspw_print_Efield_grad	*
*     *					*
*     ***********************************
      subroutine pspw_print_Efield_grad(unit,nion,efg,
     >                                  efg_ion,efg_smooth,efg_cor)
      implicit none
      integer unit
      integer nion
      real*8 efg(3,3,*)
      real*8 efg_ion(3,3,*)
      real*8 efg_smooth(3,3,*)
      real*8 efg_cor(3,3,*)

*     **** local variables ****
      integer ii,j,k

*     **** external functions ****
      character*4 ion_aname
      external    ion_aname

      write(unit,1300)
      write(unit,1190)(ii,ion_aname(ii),
     >                ((efg(j,k,ii),k=1,3),j=1,3),
     >                  ii=1,nion)
      write(unit,*)
      write(unit,1191)(ii,ion_aname(ii),
     >                ((efg_ion(j,k,ii),k=1,3),j=1,3),
     >                  ii=1,nion)
      write(unit,*)
      write(unit,1192)(ii,ion_aname(ii),
     >                ((efg_smooth(j,k,ii),k=1,3),j=1,3),
     >                  ii=1,nion)

      write(unit,*)
      write(unit,1193)(ii,ion_aname(ii),
     >                ((efg_cor(j,k,ii),k=1,3),j=1,3),
     >                  ii=1,nion)

 1300 FORMAT(//'=== Electric Field gradient at Atoms ==='/)
 1190 FORMAT(3X,I4,1x,A4,' Electric Field Gradient =(',3F11.5,' )',
     >      /12X,        '                          (',3F11.5,' )',
     >      /12X,        '                          (',3F11.5,' )')
 1191 FORMAT(3X,I4,1x,A4,' Electric Field Gradient =(',3F11.5,' )',
     >      /12X,        ' (ionic)                  (',3F11.5,' )',
     >      /12X,        '                          (',3F11.5,' )')
 1192 FORMAT(3X,I4,1x,A4,' Electric Field Gradient =(',3F11.5,' )',
     >      /12X,        ' (electronic)             (',3F11.5,' )',
     >      /12X,        '                          (',3F11.5,' )')
 1193 FORMAT(3X,I4,1x,A4,' Electric Field Gradient =(',3F11.5,' )',
     >      /12X,        ' (paw)                    (',3F11.5,' )',
     >      /12X,        '                          (',3F11.5,' )')

      return
      end

