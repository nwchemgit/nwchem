c
c $Id$
c


********************************************************************
*
*	pspw_gen_Efield module: used to generate electric fields at atoms
*
*
*	Data structure
*
********************************************************************

*     ***********************************
*     *					*
*     *		pspw_gen_Efield		*
*     *					*
*     ***********************************
      subroutine pspw_gen_Efield(rtdb,ispin,dn,dng)
      implicit none
      integer rtdb
      integer ispin
      real*8  dn(*)
      complex*16 dng(*)

#include "bafdecls.fh"
#include "btdb.fh"
#include "errquit.fh"
#include "stdio.fh"


*     **** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)

      logical value,doefield
      integer ii,i,j,k,l,n,npack0,nfft3d,n2ft3d,nion
      integer nx,ny,nxh,np1,np2,np3
      integer Gx(2),Gy(2),Gz(2)
      integer tmp1(2),G(3),E1(2),Ea(2),Ee(2),vc(2)
      integer tmpx(2),tmpy(2),tmpz(2)
      real*8  coef,x,gg,fourpi,omega,scal1,scal2
      real*8  w,wx,wy,wz,p(3)

*     **** external functions ****
      integer  ion_nion,ion_katm,G_indx
      integer  control_version
      real*8   lattice_omega,psp_zv,ion_rion
      external ion_nion, ion_katm,G_indx
      external control_version
      external lattice_omega,psp_zv,ion_rion

      if (.not.btdb_get(rtdb,'nwpw:efield',mt_log,1,doefield)) 
     >  doefield=.true.

      if (doefield) then
       
       call Parallel_taskid(taskid)
       call Pack_npack(0,npack0)
       call D3dB_nfft3d(1,nfft3d)
       n2ft3d = 2*nfft3d
       nion = ion_nion()

*      **** get memory from stack ****
       value =           BA_push_get(mt_dcpl,nfft3d,'vc',vc(2),vc(1))
       value = value.and.BA_push_get(mt_dcpl,nfft3d,
     >                               'tmpx',tmpx(2),tmpx(1))
       value = value.and.BA_push_get(mt_dcpl,nfft3d,
     >                               'tmpy',tmpy(2),tmpy(1))
       value = value.and.BA_push_get(mt_dcpl,nfft3d,
     >                               'tmpz',tmpz(2),tmpz(1))
       value = value.and.BA_push_get(mt_dcpl,nfft3d,
     >                               'tmp1',tmp1(2),tmp1(1))
       value = value.and.BA_push_get(mt_dbl,nfft3d,'Gx',Gx(2),Gx(1))
       value = value.and.BA_push_get(mt_dbl,nfft3d,'Gy',Gy(2),Gy(1))
       value = value.and.BA_push_get(mt_dbl,nfft3d,'Gz',Gz(2),Gz(1))
       value = value.and.BA_push_get(mt_dbl,(3*nion),'E1',E1(2),E1(1))
       value = value.and.BA_push_get(mt_dbl,(3*nion),'Ea',Ea(2),Ea(1))
       value = value.and.BA_push_get(mt_dbl,(3*nion),'Ee',Ee(2),Ee(1))
       if (.not. value) call errquit('out of stack memory',0, MA_ERR)


*      *** generate electronic part ***
       if (control_version().eq.3) then
          call coulomb_v(dng,dcpl_mb(vc(1)))
          scal2 = 1.0d0/dsqrt(lattice_omega())
          scal2 = 1.0d0
c          call Pack_c_SMul(0,scal2,dcpl_mb(vc(1)),dcpl_mb(vc(1)))
          call Pack_c_SMul1(0,scal2,dcpl_mb(vc(1)))
       end if

       if (control_version().eq.4)  then
          call D3dB_rr_Sum(1,dn(1),dn(1+(ispin-1)*n2ft3d),
     >                    dcpl_mb(tmp1(1)))
          call coulomb2_v(dcpl_mb(tmp1(1)),dcpl_mb(vc(1)))

*         **** fft potential  ****
          call D3dB_nx(1,np1)
          call D3dB_ny(1,np2)
          call D3dB_nz(1,np3)
          scal1=1.0d0/dble(np1*np2*np3)
c          call D3dB_r_SMul(1,scal1,dcpl_mb(vc(1)),dcpl_mb(vc(1)))
          call D3dB_r_SMul1(1,scal1,dcpl_mb(vc(1)))
          call D3dB_r_Zero_Ends(1,dcpl_mb(vc(1)))
          call D3dB_rc_fft3f(1,dcpl_mb(vc(1)))
          call Pack_c_pack(0,dcpl_mb(vc(1)))

       end if
       call dcopy(3*nion,0.0d0,0,dbl_mb(Ee(1)),1)

*      **** Efield at mesh points ****
       call phafac()
      call D3dB_t_Copy(1,dbl_mb(G_indx(1)),dbl_mb(Gx(1)))
      call D3dB_t_Copy(1,dbl_mb(G_indx(2)),dbl_mb(Gy(1)))
      call D3dB_t_Copy(1,dbl_mb(G_indx(3)),dbl_mb(Gz(1)))
      call Pack_t_pack(0,dbl_mb(Gx(1)))
      call Pack_t_pack(0,dbl_mb(Gy(1)))
      call Pack_t_pack(0,dbl_mb(Gz(1)))

      call Pack_itc_Mul(0,dbl_mb(Gx(1)),dcpl_mb(vc(1)),dcpl_mb(tmpx(1)))
      call Pack_itc_Mul(0,dbl_mb(Gy(1)),dcpl_mb(vc(1)),dcpl_mb(tmpy(1)))
      call Pack_itc_Mul(0,dbl_mb(Gz(1)),dcpl_mb(vc(1)),dcpl_mb(tmpz(1)))
       do ii=1,nion
c       p(1) = ion_rion(1,ii)
c       p(2) = ion_rion(2,ii)
c       p(3) = ion_rion(3,ii)
c        do l=1,nfft3d
c          w = dbl_mb(Gx+l-1)*p(1)
c     >      + dbl_mb(Gy+l-1)*p(2)
c     >      + dbl_mb(Gz+l-1)*p(3)
c          dcpl_mb(tmp1(1)+l-1)=dcmplx(dcos(w),-dsin(w))
c        end do
c        call Pack_c_pack(0,dcpl_mb(tmp1(1)))
        call strfac_pack(0,ii,dcpl_mb(tmp1(1)))
       
        call Pack_cc_dot(0,dcpl_mb(tmp1(1)),dcpl_mb(tmpx(1)),wx)
        call Pack_cc_dot(0,dcpl_mb(tmp1(1)),dcpl_mb(tmpy(1)),wy)
        call Pack_cc_dot(0,dcpl_mb(tmp1(1)),dcpl_mb(tmpz(1)),wz)
        dbl_mb(Ee(1)+     (ii-1)*3) = wx
        dbl_mb(Ee(1)+ 1 + (ii-1)*3) = wy
        dbl_mb(Ee(1)+ 2 + (ii-1)*3) = wz
      end do



*      *** generate ion part ***
       call dcopy(3*nion,0.0d0,0,dbl_mb(Ea(1)),1)
       if (control_version().eq.3) call ewald_f(dbl_mb(Ea(1)))
       if (control_version().eq.4) call ion_ion_f(dbl_mb(Ea(1)))
       do ii=1,nion
       do j=1,3
         dbl_mb(Ea(1)+j-1+(ii-1)*3)= dbl_mb(Ea(1)+j-1+(ii-1)*3)
     >                             /psp_zv(ion_katm(ii))
       end do
       end do

*      *** Add up electronic and ion parts ***
       do ii=1,nion
       do j=1,3
        dbl_mb(E1(1)+j-1+(ii-1)*3)= dbl_mb(Ea(1)+j-1+(ii-1)*3)
     >                            + dbl_mb(Ee(1)+j-1+(ii-1)*3)
       end do
       end do


       if (taskid.eq.MASTER) then
         call pspw_print_Efield(luout,nion,dbl_mb(E1(1)),
     >                                dbl_mb(Ea(1)),
     >                                dbl_mb(Ee(1)))
       end if

       value =           BA_pop_stack(Ee(2))
       value = value.and.BA_pop_stack(Ea(2))
       value = value.and.BA_pop_stack(E1(2))
       value = value.and.BA_pop_stack(Gz(2))
       value = value.and.BA_pop_stack(Gy(2))
       value = value.and.BA_pop_stack(Gx(2))
       value = value.and.BA_pop_stack(tmp1(2))
       value = value.and.BA_pop_stack(tmpz(2))
       value = value.and.BA_pop_stack(tmpy(2))
       value = value.and.BA_pop_stack(tmpx(2))
       value = value.and.BA_pop_stack(vc(2))
       if (.not. value) 
     >    call errquit('pspw_gen_Efield:popping of stack memory',0,
     >    MA_ERR)


      end if

      return
      end

*     ***********************************
*     *					*
*     *		pspw_print_Efield	*
*     *					*
*     ***********************************
      subroutine pspw_print_Efield(unit,nion,E1,Ea,Ee)
      implicit none
      integer unit
      integer nion
      real*8 E1(3,*)
      real*8 Ea(3,*)
      real*8 Ee(3,*)

*     **** local variables ****
      integer i,k

*     **** external functions ****
      character*4 ion_aname
      external    ion_aname

      write(unit,1300)
      write(unit,1190)(i,ion_aname(i),
     >                  (E1(k,i),k=1,3),
     >                  (Ea(k,i),k=1,3),
     >                  (Ee(k,i),k=1,3),
     >                  i=1,nion)
 1300 FORMAT(//'=== Electric Field at Atoms ==='/)
 1190 FORMAT(3X,I4,1x,A4,' Atomic Electric Field =(',3F11.5,' )',
     >      /12X,        '        (ion)          =(',3F11.5,' )',
     >      /12X,        '        (electronic)   =(',3F11.5,' )')

      return
      end

