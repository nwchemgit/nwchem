c
c $Id$
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

#include "bafdecls.fh"
#include "btdb.fh"
#include "stdio.fh"
#include "errquit.fh"

*     **** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)

      logical value,doefield
      integer nion
      integer efg(2)   ! total efg
      integer efg_ion(2) ! ionic part
      integer efg_smoothr(2),efg_smoothi(2) ! real and imag parts of smooth electronic part
      integer efg_paw(2)  ! paw part

*     **** external functions ****
      integer  ion_nion,control_version
      external ion_nion,control_version

      if (.not.btdb_get(rtdb,'nwpw:efield_grad',mt_log,1,doefield)) 
     >  doefield=.false.

      if (doefield) then
       
       call Parallel_taskid(taskid)
       nion = ion_nion()
       value =           BA_push_get(mt_dbl,9*nion,'efg',
     >                               efg(2),efg(1))
       value = value.and.BA_push_get(mt_dbl,9*nion,'efg_ion',
     >                               efg_ion(2),efg_ion(1))
       value = value.and.BA_push_get(mt_dbl,9*nion,'efg_smooth_real',
     >                               efg_smoothr(2),efg_smoothr(1))
       value = value.and.BA_push_get(mt_dbl,9*nion,'efg_smooth_imag',
     >                               efg_smoothi(2),efg_smoothi(1))
       value = value.and.BA_push_get(mt_dbl,9*nion,'efg_paw',
     >                               efg_paw(2),efg_paw(1))
       if (.not.value) 
     >    call errquit('pspw_gen_Efield_grad:out of stack memory',0,
     >                 MA_ERR)

       call dcopy(9*nion,0.0d0,0,dbl_mb(efg_ion(1)),1)
       call dcopy(9*nion,0.0d0,0,dbl_mb(efg_smoothr(1)),1)
       call dcopy(9*nion,0.0d0,0,dbl_mb(efg_smoothi(1)),1)
       call dcopy(9*nion,0.0d0,0,dbl_mb(efg_paw(1)),1)
       call dcopy(9*nion,0.0d0,0,dbl_mb(efg(1)),1)

       call ewald_efg(dbl_mb(efg_ion(1)))   ! ionic part
       call coulomb_efg(dng,dbl_mb(efg_smoothr(1)),   ! smooth electronic part
     &                      dbl_mb(efg_smoothi(1)) )
       call psp_efg_atoms(ispin,ne,psi,dbl_mb(efg_paw(1)))  ! paw part

       call daxpy(9*nion,1.d0,dbl_mb(efg_ion(1)),1,dbl_mb(efg(1)),1)
!       call daxpy(9*nion,1.d0,dbl_mb(efg_smoothr(1)),1,dbl_mb(efg(1)),1)
!       call daxpy(9*nion,1.d0,dbl_mb(efg_paw(1)),1,dbl_mb(efg(1)),1)

       if (taskid.eq.MASTER) then
         call pspw_final_Efield_grad(luout,
     >                               nion,
     >                               dbl_mb(efg(1)),
     >                               dbl_mb(efg_ion(1)),
     >                               dbl_mb(efg_smoothr(1)),
     >                               dbl_mb(efg_paw(1)))
       end if



       value =           BA_pop_stack(efg_paw(2))
       value = value.and.BA_pop_stack(efg_smoothi(2))
       value = value.and.BA_pop_stack(efg_smoothr(2))
       value = value.and.BA_pop_stack(efg_ion(2))
       value = value.and.BA_pop_stack(efg(2))
       if (.not.value) 
     >    call errquit('pspw_gen_Efield_grad:popping of stack memory',0,
     >                 MA_ERR)

      end if

      return
      end

*     ***********************************
*     *					*
*     *	    pspw_final_Efield_grad	*
*     *					*
*     ***********************************
      subroutine pspw_final_Efield_grad(unit,nion,efg,
     >                                  efg_ion,efg_smooth,efg_paw)
      implicit none

      integer unit
      integer nion
      real*8 efg(3,3,*)
      real*8 efg_ion(3,3,*)
      real*8 efg_smooth(3,3,*)
      real*8 efg_paw(3,3,*)
      real*8 eta,eigval(3),vec(3,3),tmp1(3),tmp2(3)
      real*8 eigvec(3,3),cq,eval(3),evec(3,3)
      integer ind(3)

*     **** local variables ****
      integer ii,i,j,k

*     **** external functions ****
      character*4 ion_aname
      external    ion_aname

*     **** print efg tensors for atoms ****
      write(unit,1300)
      do ii=1,nion
       write(unit,1190) ii,ion_aname(ii),((efg(j,k,ii),k=1,3),j=1,3)
       write(unit,1191) ii,ion_aname(ii),((efg_ion(j,k,ii),k=1,3),j=1,3)
       write(unit,1192) ii,ion_aname(ii),
     &         ((efg_smooth(j,k,ii),k=1,3),j=1,3)
       write(unit,1193) ii,ion_aname(ii),
     &         ((efg_paw(j,k,ii),k=1,3),j=1,3)
       write(unit,*)
      end do  ! ii

*     **** compute principal components and spectroscopic parameters ****
      write(unit,1194)
      do ii = 1,nion
        do i=1,3
         ind(i) = 0
         do j=1,3
            vec(i,j) = efg(i,j,ii)
         end do
        end do
        call util_diag_jacobi(unit,3,3,vec,eval,evec,tmp1,tmp2) ! diagonalize matrix
        do i = 1,3
           eigval(i) = eval(i)
        end do
        call util_sorteps(3, dabs(eval), ind, 1.d-4)
        do i=1,3
        !  write(unit,*) i,eigval(ind(i)),ind(i)
           eigval(i) = eigval(ind(i))
        end do
        write(unit,*)

        eta = (eigval(1)-eigval(2))/eigval(3) ! asymmetry parameter
        if (abs(eigval(3)) <= 1.d-4) eta = 0.d0
        cq = eigval(3)  ! quadupolar coupling constant
        write(unit,1199) ii,ion_aname(ii)
        write(unit,1195) eigval(1),eigval(2),eigval(3),eta,cq
        write(unit,1196) ((eigvec(i,j),j=1,3),i=1,3)
        write(unit,*)
      enddo ! ii
c
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
 1194 FORMAT(1x,'Principal components (a.u.) and orientation ',
     2       10x,'eta',10x,'Cq'/,1x,86(1h-))
 1195 FORMAT(1x,3f11.5,15x,f11.5,2x,f11.5/)
 1196 FORMAT(1x,3F11.5)
 1199 FORMAT(3X,I4,1x,A4)

      return
      end
