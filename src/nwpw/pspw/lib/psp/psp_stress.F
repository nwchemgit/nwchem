*
* $Id$
*


#define TCGMSG

*     ***********************************
*     *					*
*     *	     psp_stress_init 		*
*     *					*
*     ***********************************

      subroutine psp_stress_init()      
      implicit none

#include "bafdecls.fh"
#include "psp.fh"
#include "errquit.fh"

      integer npack1,npack0
      logical value

*     **** external functions *****
      integer  ion_nkatm
      external ion_nkatm
 

      call Pack_npack(1,npack1)
      call Pack_npack(0,npack0)

      value = BA_alloc_get(mt_dbl,(npsp*npack0),
     >                    'dvl',dvl(2),dvl(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,(3*nmax_max*lmmax_max*npsp*npack1),
     >                    'dvnl',dvnl(2),dvnl(1))

      if (.not. value) call errquit('out of heap memory',0, MA_ERR)

      call dcopy(npsp*npack0,0.0d0,0,dbl_mb(dvl(1)),1)
      call dcopy(3*nmax_max*lmmax_max*npsp*npack1,
     >           0.0d0,0,
     >           dbl_mb(dvnl(1)),1)
      return
      end



*     ***********************************
*     *					*
*     *	     psp_stress_end 		*
*     *					*
*     ***********************************

      subroutine psp_stress_end()      
      implicit none

#include "errquit.fh"
#include "bafdecls.fh"
#include "psp.fh"


      logical value

      value =           BA_free_heap(dvl(2))
      value = value.and.BA_free_heap(dvnl(2))

      if (.not.value) 
     > call errquit('psp_stress_end:freeing heap',0, MA_ERR)
 
      return
      end


*     ***********************************
*     *					*
*     *	 	   v_local_euv  	*
*     *					*
*     ***********************************

      subroutine v_local_euv(dng,euv)
      implicit none
      complex*16 dng(*)
      real*8     euv(3,3)

#include "bafdecls.fh"
#include "psp.fh"
#include "errquit.fh"

*     *** local variables ***
      integer nfft3d,npack0
      integer ii,ia,u,v,s
      integer exi(2),tmp1(2),tmp2(2)
      integer G(2,3),vll(2)
      logical value
      real*8 elocal,ftmp(3)
      real*8 hm(3,3),Bus(3,3)
      real*8 ss,sum,pi,fourpi

*     **** common block used for coulomb.f ****
      integer vc_indx,vc_hndl
      common / vc_block / vc_indx,vc_hndl


*     **** external functions ****
      integer  G_indx,ion_nion,ion_katm
      real*8   lattice_omega,lattice_unitg
      external G_indx,ion_nion,ion_katm
      external lattice_omega,lattice_unitg

      call nwpw_timing_start(5)

      call D3dB_nfft3d(1,nfft3d)
      call Pack_npack(0,npack0)

      pi     = 4.0d0*datan(1.0d0)
      fourpi = 4.0d0*pi
      ss   = 1.0d0/(2.0d0*pi)

*     *** define hm ****
      do v=1,3
      do u=1,3
         hm(u,v) = ss*lattice_unitg(u,v)
      end do
      end do

*     **** average Kohn-Sham v_local energy ****
      value = BA_push_get(mt_dcpl,npack0,'vll',vll(2),vll(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)
      call v_local(dcpl_mb(vll(1)),.false.,dng,ftmp)
      call Pack_cc_dot(0,dng,dcpl_mb(vll(1)),elocal)
      value = BA_pop_stack(vll(2))
      if (.not. value) call errquit('error poping stack memory',0,
     &       MA_ERR)


      value = BA_push_get(mt_dcpl,nfft3d,'exi', exi(2), exi(1))
      value = value.and.
     >        BA_push_get(mt_dbl, npack0,'tmp1',tmp1(2),tmp1(1))
      value = value.and.
     >        BA_push_get(mt_dbl, npack0,'tmp2',tmp2(2),tmp2(1))
      value = value.and.
     >        BA_push_get(mt_dbl, nfft3d,'G1',G(2,1),G(1,1))
      value = value.and.
     >        BA_push_get(mt_dbl, nfft3d,'G2',G(2,2),G(1,2))
      value = value.and.
     >        BA_push_get(mt_dbl, nfft3d,'G3',G(2,3),G(1,3))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

*     **** define Gx,Gy and Gz in packed space ****
      call D3dB_t_Copy(1,dbl_mb(G_indx(1)),dbl_mb(G(1,1)))
      call D3dB_t_Copy(1,dbl_mb(G_indx(2)),dbl_mb(G(1,2)))
      call D3dB_t_Copy(1,dbl_mb(G_indx(3)),dbl_mb(G(1,3)))
      call Pack_t_pack(0,dbl_mb(G(1,1)))
      call Pack_t_pack(0,dbl_mb(G(1,2)))
      call Pack_t_pack(0,dbl_mb(G(1,3)))


      call dcopy(9,0.0d0,0,Bus,1)
      do ii=1,ion_nion()
        ia=ion_katm(ii)

*       **** structure factor and local pseudopotential ****
        !call strfac(ii,dcpl_mb(exi(1)))
        !call Pack_c_pack(0,dcpl_mb(exi(1)))
        call strfac_pack(0,ii,dcpl_mb(exi(1)))
        
*       **** tmp2(G) = Real[ dconjg(dng(G))*exi(G) ] ****
        call Pack_cct_conjgMul(0,dng,
     >                          dcpl_mb(exi(1)),
     >                          dbl_mb(tmp2(1)))

*       **** tmp2(G) = tmp2(G)*(dvl(G))
c        call Pack_tt_Mul(0,dbl_mb(tmp2(1)),
c     >                     dbl_mb(dvl(1)+(ia-1)*npack0),
c     >                     dbl_mb(tmp2(1)))
        call Pack_tt_Mul2(0,dbl_mb(dvl(1)+(ia-1)*npack0),
     >                      dbl_mb(tmp2(1)))

*       **** tmp2(G) = tmp2(G)/G ****
        ss     = 1.0d0/fourpi
        call Pack_t_SMul(0,ss,dbl_mb(vc_indx),dbl_mb(tmp1(1)))
c        call Pack_tt_Sqrt(0,dbl_mb(tmp1(1)),dbl_mb(tmp1(1)))
c        call Pack_tt_Mul(0,dbl_mb(tmp1(1)),
c     >                     dbl_mb(tmp2(1)),
c     >                     dbl_mb(tmp2(1)))
        call Pack_tt_Sqrt1(0,dbl_mb(tmp1(1)))
        call Pack_tt_Mul2(0,dbl_mb(tmp1(1)),
     >                      dbl_mb(tmp2(1)))
        
*       **** Bus = Bus - Sum(G) tmp2(G)*Gu*Gs ***
        do u=1,3
        do s=u,3
          call Pack_tt_Mul(0,dbl_mb(G(1,u)),
     >                       dbl_mb(G(1,s)),
     >                       dbl_mb(tmp1(1)))
          call Pack_tt_dot(0,dbl_mb(tmp1(1)),dbl_mb(tmp2(1)),sum)

          Bus(u,s) = Bus(u,s) - sum
        end do
        end do

      end do
      value =           BA_pop_stack(G(2,3))
      value = value.and.BA_pop_stack(G(2,2))
      value = value.and.BA_pop_stack(G(2,1))
      value = value.and.BA_pop_stack(tmp2(2))
      value = value.and.BA_pop_stack(tmp1(2))
      value = value.and.BA_pop_stack(exi(2))
      if (.not. value) call errquit('error poping stack memory',0,
     &       MA_ERR)

      do u=1,3
      do s=u+1,3
         Bus(s,u) = Bus(u,s)
      end do
      end do
      do v=1,3
      do u=1,3
         euv(u,v) = -elocal*hm(u,v)
         do s=1,3
            euv(u,v) = euv(u,v) + Bus(u,s)*hm(s,v)
         end do
      end do
      end do


      call nwpw_timing_end(5)
      return 
      end


*     ***********************************
*     *					*
*     *	 	   v_nonlocal_euv_2	*
*     *					*
*     ***********************************

      subroutine v_nonlocal_euv_2(ispin,ne,psi1,euv)
      implicit none
      integer    ispin,ne(2)
      complex*16 psi1(*)
      real*8 euv(3,3)

#include "bafdecls.fh"
#include "psp.fh"
#include "errquit.fh"


*     *** local variables ***
      integer nfft3d,npack1,shift,shift2,np_j
      integer nproj,l_prj
      integer i,ii,ia,k,l,nn
      integer s,u,v
      real*8  omega,Bus(3,3),hm(3,3)
      real*8  pi,scal
      integer exi(2),vtmp(2),tmp1(2),sw1(2),sw2(2),sw3(2)
      integer G(2,3),Gijl_indx
      logical value,sd_function

*     **** external functions ****
      integer  ion_nion,ion_katm,G_indx,psi_data_get_ptr
      real*8   lattice_omega,lattice_unitg
      external ion_nion,ion_katm,G_indx,psi_data_get_ptr
      external lattice_omega,lattice_unitg

      call nwpw_timing_start(6)
      call Parallel2d_np_j(np_j)
      
      
*     **** allocate local memory ****
      nn = ne(1)+ne(2)
      call D3dB_nfft3d(1,nfft3d)
      call Pack_npack(1,npack1)

      value = BA_push_get(mt_dcpl,nfft3d,'exi', exi(2), exi(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,nfft3d,'vtmp',vtmp(2),vtmp(1))
      value = value.and.
     >        BA_push_get(mt_dbl, npack1,'tmp1',tmp1(2),tmp1(1))
      value = value.and.
     >        BA_push_get(mt_dbl, nfft3d,'Gx',G(2,1),G(1,1))
      value = value.and.
     >        BA_push_get(mt_dbl, nfft3d,'Gy',G(2,2),G(1,2))
      value = value.and.
     >        BA_push_get(mt_dbl, nfft3d,'Gz',G(2,3),G(1,3))
      value = value.and.
     >        BA_push_get(mt_dbl,nn*nmax_max*lmmax_max,
     >                    'sw1',sw1(2),sw1(1))
      value = value.and.
     >        BA_push_get(mt_dbl,nn*nmax_max*lmmax_max,
     >                    'sw2',sw2(2),sw2(1))
      value = value.and.
     >        BA_push_get(mt_dbl,9*nn,
     >                    'sw3',sw3(2),sw3(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

      
*     **** define Gx,Gy and Gz in packed space ****
      call D3dB_t_Copy(1,dbl_mb(G_indx(1)),dbl_mb(G(1,1)))
      call D3dB_t_Copy(1,dbl_mb(G_indx(2)),dbl_mb(G(1,2)))
      call D3dB_t_Copy(1,dbl_mb(G_indx(3)),dbl_mb(G(1,3)))
      call Pack_t_pack(1,dbl_mb(G(1,1)))
      call Pack_t_pack(1,dbl_mb(G(1,2)))
      call Pack_t_pack(1,dbl_mb(G(1,3)))

      omega = lattice_omega()

      
*     ***********************
*     **** calculate Bus ****
*     ***********************
      call dcopy(9,0.0d0,0,Bus,1)
      do ii=1,ion_nion()
        ia=ion_katm(ii)
        nproj=int_mb(nprj(1)+ia-1)

      
        if (nproj.gt.0) then
*       **** structure factor and local pseudopotential ****
        !call strfac(ii,dcpl_mb(exi(1)))
        !call Pack_c_pack(1,dcpl_mb(exi(1)))
        call strfac_pack(1,ii,dcpl_mb(exi(1)))

  
*       *********************************************
*       **** calculate F^(lm)_I = <psi|vnl(nlm)> ****
*       *********************************************
        do l=1,nproj

           !shift = vnl(1)+(l-1)*npack1+(ia-1)*npack1*nmax_max*lmmax_max
           shift = psi_data_get_ptr(int_mb(vnl(1)+ia-1),l)
           l_prj = int_mb(l_projector(1)+(l-1)
     >                                  + (ia-1)*(nmax_max*lmmax_max))
           !sd_function = .not.and(l_prj,1)
#ifdef GCC4
           k = iand(l_prj,1)
#else
           k = and(l_prj,1)
#endif
           sd_function = (k.eq.0)

      
*          *** current function is s or d ****
           if (sd_function) then
              call Pack_tc_Mul(1,dbl_mb(shift),
     >                        dcpl_mb(exi(1)),
     >                        dcpl_mb(vtmp(1)))
*          *** current function is p or f ****
           else
              call Pack_itc_Mul(1,dbl_mb(shift),
     >                        dcpl_mb(exi(1)),
     >                        dcpl_mb(vtmp(1)))
           end if
           call Pack_cc_indot(1,nn,
     >                      psi1,
     >                      dcpl_mb(vtmp(1)),
     >                      dbl_mb(sw1(1)+(l-1)*nn))

        end do !** l**
        call D3dB_Vector_Sumall((nn*nproj),dbl_mb(sw1(1)))


      
*       **** sw2 = Gijl*sw1 ******
        Gijl_indx = psi_data_get_ptr(int_mb(Gijl(1)+ia-1),1)
        call Multiply_Gijl_sw1(nn,
     >                         nproj,
     >                         int_mb(nmax(1)+ia-1),
     >                         int_mb(lmax(1)+ia-1),
     >                         int_mb(n_projector(1)
     >                                + (ia-1)*(nmax_max*lmmax_max)),
     >                         int_mb(l_projector(1)
     >                                + (ia-1)*(nmax_max*lmmax_max)),
     >                         int_mb(m_projector(1)
     >                                + (ia-1)*(nmax_max*lmmax_max)),
     >                         dbl_mb(Gijl_indx),
     >                         dbl_mb(sw1(1)),
     >                         dbl_mb(sw2(1)))
      
        if (ispin.eq.1) call dscal(nn*nproj,2.0d0,dbl_mb(sw2(1)),1)
      
 

*       **********************************
*       **** calculate dF^(lm)_I/dhus ****
*       **********************************

        do l=1,nproj
            l_prj = int_mb(l_projector(1)
     >                    + (l-1)
     >                    + (ia-1)*(nmax_max*lmmax_max))
#ifdef GCC4
           k = iand(l_prj,1)
#else
           k = and(l_prj,1)
#endif
           sd_function = (k.eq.0)
      


           do s=1,3
           do u=1,3
              shift2 = dvnl(1) + (s-1)*npack1
     >                         + (l-1)*3*npack1
     >                         +(ia-1)*3*npack1*lmmax_max*nmax_max
              call Pack_tt_Mul(1,dbl_mb(shift2),
     >                         dbl_mb(G(1,u)),
     >                         dbl_mb(tmp1(1)))

*             *** current function is s or d ****
              if (sd_function) then
                call Pack_tc_Mul(1,dbl_mb(tmp1(1)),
     >                        dcpl_mb(exi(1)),
     >                        dcpl_mb(vtmp(1)))
*             *** current function is p or f ****
              else
                 call Pack_itc_Mul(1,dbl_mb(tmp1(1)),
     >                        dcpl_mb(exi(1)),
     >                        dcpl_mb(vtmp(1)))
              end if
              call Pack_cc_ndot(1,nn,
     >                        psi1,
     >                        dcpl_mb(vtmp(1)),
     >                        dbl_mb(sw3(1)+(u-1)*nn
     >                                     +(s-1)*nn*3))
          end do
          end do

          do i=1,nn
            do s=1,3
            do u=1,3

             Bus(u,s) = Bus(u,s) 
     >         - 2.0d0/(omega)
     >          * dbl_mb(sw2(1)+(i-1)+(l-1)*nn)
     >          * dbl_mb(sw3(1)+(i-1)
     >                         +(u-1)*nn
     >                         +(s-1)*nn*3)
            end do
            end do
          end do

        end do !** l **

        end if
      end do !** ii **
      if (np_j.gt.1) call D1dB_Vector_SumAll(9,Bus)

      value =           BA_pop_stack(sw3(2))
      value = value.and.BA_pop_stack(sw2(2))
      value = value.and.BA_pop_stack(sw1(2))
      value = value.and.BA_pop_stack(G(2,3))
      value = value.and.BA_pop_stack(G(2,2))
      value = value.and.BA_pop_stack(G(2,1))
      value = value.and.BA_pop_stack(tmp1(2))
      value = value.and.BA_pop_stack(vtmp(2))
      value = value.and.BA_pop_stack(exi(2))
      if (.not. value) 
     >  call errquit('v_nonlocal_euv_2:error poping stack',0,
     &       MA_ERR)


*     *** define hm ****
      pi   = 4.0d0*datan(1.0d0)
      scal = 1.0d0/(2.0d0*pi)
      do v=1,3
      do u=1,3
         hm(u,v) = scal*lattice_unitg(u,v)
      end do
      end do

*     *** calculate euv = Sum(s) hm(s,v)*Bus(u,s)
      call dcopy(9,0.0d0,0,euv,1)
      do u=1,3
      do v=1,3
         do s=1,3
            euv(u,v) = euv(u,v) + Bus(u,s)*hm(s,v)
         end do
      end do
      end do

      call nwpw_timing_end(6)
      return 
      end



*     ***********************************
*     *					*
*     *	  	psp_stress_read	        *
*     *					*
*     ***********************************

      subroutine psp_stress_read(fname,
     >                       version,
     >                       nfft,unita,
     >                       npack0,dvl,
     >                       npack1,nprj,dvnl,
     >                       semicore,dncore,
     >                       tmp,tmp2,
     >                       ierr)
      implicit none 
      character*50 fname
      integer version
      integer nfft(3)
      real*8  unita(3,3)
      integer npack0
      real*8 dvl(*)
      integer npack1,nprj
      real*8 dvnl(npack1,3,nprj)
      logical semicore
      real*8 dncore(*)
      complex*16 tmp(*)
      real*8     tmp2(*)
      integer ierr

#ifdef MPI
      include 'mpif.h'
      integer mpierr
#endif
#ifdef TCGMSG
#include "tcgmsg.fh"
#include "msgtypesf.h"
#endif

*    *** local variables ***
      logical pio
      integer MASTER,taskid,taskid_i,taskid_p,com_p
      parameter(MASTER=0)
      integer i,n,l
      integer msglen
      character*255 full_filename

      logical  control_parallel_io
      external control_parallel_io

      call Parallel_taskid(taskid)
      call Parallel2d_taskid_i(taskid_i)

      pio = control_parallel_io()
      if (pio) then
         taskid_p = taskid_i
         com_p = 1
      else
         taskid_p = taskid
         com_p = 0
      end if
  
*     **** open fname binary file ****
      if (taskid_p.eq.MASTER) then
         call util_file_name_noprefix(fname,.false.,
     >                             .false.,
     >                       full_filename)
         l = index(full_filename,' ') - 1
         call openfile(5,full_filename,l,'r',l)
         call iread(5,version,1)
         call iread(5,nfft,3)
         call dread(5,unita,9)
      end if


*     **** send header data to all processors ****
      msglen = 1
      call Parallela_Brdcst_ivalues(com_p,MASTER,msglen,version)
      msglen = 3
      call Parallela_Brdcst_ivalues(com_p,MASTER,msglen,nfft)
      msglen = 9
      call Parallela_Brdcst_values(com_p,MASTER,msglen,unita)

 
*     *** read in vl 3d block ***
      if (pio) then
         call D3dB_t_read_pio(1,5,tmp2,tmp,-1)
      else
         call D3dB_t_read(1,5,tmp2,tmp,-1)
      end if
      call Pack_t_pack(0,tmp2)
      call Pack_t_Copy(0,tmp2,dvl)

*     **** read in vnl 3d blocks ****
      do n=1,nprj
      do i=1,3
         if (pio) then
            call D3dB_t_read_pio(1,5,tmp2,tmp,-1)
         else
            call D3dB_t_read(1,5,tmp2,tmp,-1)
         end if
         call Pack_t_pack(1,tmp2)
         call Pack_t_Copy(1,tmp2,dvnl(1,i,n))
      end do
      end do

*     **** read in semicore density block ****
      if (semicore) then
         !write(*,*) "reading in semicore block"   !debug
         if (pio) then
            call D3dB_t_read_pio(1,5,tmp2,tmp,-1)
         else
            call D3dB_t_read(1,5,tmp2,tmp,-1)
         end if
         call Pack_t_pack(0,tmp2)
         call Pack_t_Copy(0,tmp2,dncore)
      end if

*     *** close fname binary file ***
      if (taskid_p.eq.MASTER) then
c       close(11)
         call closefile(5)
      end if

      ierr = 0
      return
      end

*     ***********************************
*     *					*
*     *	     psp_stress_readall 	*
*     *					*
*     ***********************************

      subroutine psp_stress_readall()
      implicit none

#include "bafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "errquit.fh"      

#include "psp.fh"

*     **** semicore common block ****
      integer ncore(2),rcore(2)
      integer semicore(2)
      common / ccore / ncore,rcore,semicore


*     **** local variables ****
      integer ngp(3)
      real*8  unita(3,3)
      integer version,nfft3d,npack1,npack0
      integer ia,l,nproj
      character*12 boundry
      integer tmp(2),tmp2(2),ierr
      logical value,found,correct_box
      character*5  element
      character*50 fname

*     **** parallel i/o variable ****
      integer MASTER,taskid
      parameter(MASTER=0)

*     **** external functions ****
      logical      nwpw_filefind,control_print
      integer      control_ngrid
      integer      psp_nprj
      real*8       control_unita
      character*12 control_boundry
      character*4  ion_atom
      external     nwpw_filefind,control_print
      external     control_ngrid
      external     psp_nprj
      external     control_unita
      external     control_boundry
      external     ion_atom
      

      call D3dB_nfft3d(1,nfft3d)
      call Pack_npack(1,npack1)
      call Pack_npack(0,npack0)
      call Parallel_taskid(taskid)

      value = BA_push_get(mt_dbl,(2*nfft3d),'tmp',tmp(2),tmp(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

      value = BA_push_get(mt_dbl,(nfft3d),'tmp2',tmp2(2),tmp2(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

*     **** read pseudopotentials ****
      do ia=1,npsp

*       **** define formatted psp name ****
         element = '     '
         element = ion_atom(ia)
         l = index(element,' ') - 1
         fname = element(1:l)//'.vpp2'
         nproj=int_mb(nprj(1)+ia-1)
    

*        **** not finished ****
         found = .false.
         do while (.not. found)
         if (nwpw_filefind(fname)) then
            call psp_stress_read(fname,
     >                  version,
     >                  ngp,unita,
     >                  npack0,
     >                  dbl_mb(dvl(1) + (ia-1)*npack0),
     >                  npack1,nproj,
     >                  dbl_mb(dvnl(1) 
     >                        +(ia-1)*3*npack1*nmax_max*lmmax_max),
     >                  log_mb(semicore(1)+ia),
     >                  dbl_mb(ncore(1) + npack0 + (ia-1)*5*npack0),
     >                  dbl_mb(tmp(1)),dbl_mb(tmp2(1)),
     >                  ierr)

           if (ierr.gt.0) go to 9000
           boundry = control_boundry()
           l =index(boundry,' ') - 1
           if ((boundry(1:l).ne.'periodic') .or. (version.ne.3)) then
              if ((taskid.eq.MASTER).and.
     >            (control_print(print_medium))) then
                 write(luout,*) 
     >             "Error reading formatted pseudoptental stress file:"
                 write(luout,*) "Formatted pseudopotential was not",
     >                   " generated with periodic boundary condtions"
              end if
           end if

           correct_box = .true.
           if ( (ngp(1).ne.control_ngrid(1)) .or.
     >       (ngp(2).ne.control_ngrid(2)) .or.
     >       (ngp(3).ne.control_ngrid(3)) .or.
     >       (unita(1,1).ne.control_unita(1,1)) .or.
     >       (unita(2,1).ne.control_unita(2,1)) .or.
     >       (unita(3,1).ne.control_unita(3,1)) .or.
     >       (unita(1,2).ne.control_unita(1,2)) .or.
     >       (unita(2,2).ne.control_unita(2,2)) .or.
     >       (unita(3,2).ne.control_unita(3,2)) .or.
     >       (unita(1,3).ne.control_unita(1,3)) .or.
     >       (unita(2,3).ne.control_unita(2,3)) .or.
     >       (unita(3,3).ne.control_unita(3,3)) .or.
     >       ((boundry(1:l).eq.'periodic').and.(version.ne.3)).or.
     >       ((boundry(1:l).eq.'aperiodic').and.(version.ne.4))) then
              correct_box = .false.
              if ((taskid.eq.MASTER).and.
     >            (control_print(print_medium))) then
              write(luout,*) 
     >         "pseudopotential is not correctly formatted:",fname
              end if
           end if
           if (correct_box) found = .true.
         end if

*        **** generate formatted pseudopotential atom.vpp2 *****
         if (.not.found) then
             call psp_stress_formatter_auto(ion_atom(ia))
         end if

         end do !*** do while ****
 
      end do
 9000 value = BA_pop_stack(tmp2(2))
      value = BA_pop_stack(tmp(2))

      return
      end
