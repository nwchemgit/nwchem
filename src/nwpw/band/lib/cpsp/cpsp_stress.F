*
* $Id$
*


#define TCGMSG

*     ***********************************
*     *					*
*     *	     cpsp_stress_init 		*
*     *					*
*     ***********************************

      subroutine cpsp_stress_init()      
      implicit none

#include "bafdecls.fh"
#include "cpsp_common.fh"
#include "errquit.fh"

      integer npack1,npack0,nbrillq
      logical value

*     **** external functions *****
      integer  ion_nkatm,brillioun_nbrillq
      external ion_nkatm,brillioun_nbrillq
 

      call Cram_npack(0,npack0)
      call Cram_max_npack(npack1)
      nbrillq = brillioun_nbrillq()

      value = BA_alloc_get(mt_dbl,(npsp*npack0),
     >                    'dvl',dvl(2),dvl(1))
      value = value.and.
     >        BA_alloc_get(mt_int,npsp,'dvnl',dvnl(2),dvnl(1))

      if (.not. value) 
     > call errquit('cpsp_stress_init:out of heap memory',0,MA_ERR)

      call dcopy(npsp*npack0,0.0d0,0,dbl_mb(dvl(1)),1)
      return
      end



*     ***********************************
*     *					*
*     *	     cpsp_stress_end 		*
*     *					*
*     ***********************************

      subroutine cpsp_stress_end()      
      implicit none

#include "errquit.fh"
#include "bafdecls.fh"
#include "cpsp_common.fh"


      logical value

      value =           BA_free_heap(dvl(2))
      value = value.and.BA_free_heap(dvnl(2))

      if (.not.value) 
     > call errquit('cpsp_stress_end:freeing heap',0, MA_ERR)
 
      return
      end


*     ***********************************
*     *					*
*     *	 	   cpsp_v_local_euv  	*
*     *					*
*     ***********************************

      subroutine cpsp_v_local_euv(dng,euv)
      implicit none
      complex*16 dng(*)
      real*8     euv(3,3)

#include "bafdecls.fh"
#include "cpsp_common.fh"
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
      common / c_vc_block / vc_indx,vc_hndl


*     **** external functions ****
      integer  c_G_indx,ion_nion,ion_katm
      real*8   lattice_omega,lattice_unitg
      external c_G_indx,ion_nion,ion_katm
      external lattice_omega,lattice_unitg

      call nwpw_timing_start(5)

      call C3dB_nfft3d(1,nfft3d)
      call Cram_npack(0,npack0)

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
      if (.not. value) 
     > call errquit('cpsp_v_local_euv:out of stack memory',0,MA_ERR)
      call cpsp_v_local(dcpl_mb(vll(1)),.false.,dng,ftmp)
      call Cram_cc_dot(0,dng,dcpl_mb(vll(1)),elocal)
      value = BA_pop_stack(vll(2))
      if (.not. value) 
     > call errquit('cpsp_v_local_euv:error popping stack',0,MA_ERR)


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
      if (.not. value) 
     > call errquit('cpsp_v_local_euv:out of stack memory',0,MA_ERR)

*     **** define Gx,Gy and Gz in packed space ****
      call C3dB_t_Copy(1,dbl_mb(c_G_indx(1)),dbl_mb(G(1,1)))
      call C3dB_t_Copy(1,dbl_mb(c_G_indx(2)),dbl_mb(G(1,2)))
      call C3dB_t_Copy(1,dbl_mb(c_G_indx(3)),dbl_mb(G(1,3)))
      call Cram_r_pack(0,dbl_mb(G(1,1)))
      call Cram_r_pack(0,dbl_mb(G(1,2)))
      call Cram_r_pack(0,dbl_mb(G(1,3)))


      call dcopy(9,0.0d0,0,Bus,1)
      do ii=1,ion_nion()
        ia=ion_katm(ii)

*       **** structure factor and local pseudopotential ****
        call cstrfac_pack(0,ii,dcpl_mb(exi(1)))
        
*       **** tmp2(G) = Real[ dconjg(dng(G))*exi(G) ] ****
        call Cram_ccr_conjgMul(0,dng,
     >                          dcpl_mb(exi(1)),
     >                          dbl_mb(tmp2(1)))

*       **** tmp2(G) = tmp2(G)*(dvl(G))
c        call Cram_rr_Mul(0,dbl_mb(tmp2(1)),
c     >                     dbl_mb(dvl(1)+(ia-1)*npack0),
c     >                     dbl_mb(tmp2(1)))
        call Cram_rr_Mul2(0,dbl_mb(dvl(1)+(ia-1)*npack0),
     >                      dbl_mb(tmp2(1)))

*       **** tmp2(G) = tmp2(G)/G ****
        ss     = 1.0d0/fourpi
        call Cram_r_SMul(0,ss,dbl_mb(vc_indx),dbl_mb(tmp1(1)))
c        call Cram_rr_Sqrt(0,dbl_mb(tmp1(1)),dbl_mb(tmp1(1)))
c        call Cram_rr_Mul(0,dbl_mb(tmp1(1)),
c     >                     dbl_mb(tmp2(1)),
c     >                     dbl_mb(tmp2(1)))
        call Cram_rr_Sqrt1(0,dbl_mb(tmp1(1)))
        call Cram_rr_Mul2(0,dbl_mb(tmp1(1)),dbl_mb(tmp2(1)))
        
*       **** Bus = Bus - Sum(G) tmp2(G)*Gu*Gs ***
        do u=1,3
        do s=u,3
          call Cram_rr_Mul(0,dbl_mb(G(1,u)),
     >                       dbl_mb(G(1,s)),
     >                       dbl_mb(tmp1(1)))
          call Cram_rr_dot(0,dbl_mb(tmp1(1)),dbl_mb(tmp2(1)),sum)

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
      if (.not. value) call errquit('error popping stack memory',0,
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
*     *	    cpsp_v_nonlocal_euv_2	*
*     *					*
*     ***********************************

      subroutine cpsp_v_nonlocal_euv_2(ispin,ne,psi1_tag,euv)
      implicit none
      integer    ispin,ne(2)
      integer psi1_tag
      real*8 euv(3,3)

#include "bafdecls.fh"
#include "cpsp_common.fh"
#include "errquit.fh"


*     *** local variables ***
      integer nfft3d,npack1,npack,shift,shift2,nbrillq,nb
      integer nproj,l_prj,psi1_shift,occ1_tag,occ1_shift,occ_shift
      integer i,ii,ia,k,l,n,nn
      integer s,u,v
      real*8  omega,Bus(3,3),hm(3,3),kx,ky,kz,Aus(3,3)
      complex*16 ctmp,cxr
      real*8  pi,scal,weight
      integer exi(2),vtmp(2),tmp1(2),zsw1(2),zsw2(2),zsw3(2)
      integer G(2,3)
      logical value,sd_function

*     **** external functions ****
      logical  is_sORd
      integer  ion_nion,ion_katm,c_G_indx
      integer  Pneb_nbrillq,cpsp_projector_get_ptr
      integer  cpsi_data_get_chnk,cpsi_data_get_next
      real*8   brillioun_weight,brillioun_k
      real*8   lattice_omega,lattice_unitg
      external is_sORd
      external ion_nion,ion_katm,c_G_indx
      external Pneb_nbrillq,cpsp_projector_get_ptr
      external cpsi_data_get_chnk,cpsi_data_get_next
      external brillioun_weight,brillioun_k
      external lattice_omega,lattice_unitg

      call nwpw_timing_start(6)
      
      occ1_tag = cpsi_data_get_next(psi1_tag)
      
      omega = lattice_omega()

*     **** allocate local memory ****
      nn         = ne(1)+ne(2)
      nbrillq = Pneb_nbrillq()
      call C3dB_nfft3d(1,nfft3d)
      call Cram_max_npack(npack1)

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
     >        BA_push_get(mt_dcpl,nn*nprj_max,
     >                    'zsw1',zsw1(2),zsw1(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,nn*nprj_max,
     >                    'zsw2',zsw2(2),zsw2(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,9*nn,
     >                    'zsw3',zsw3(2),zsw3(1))
      if (.not. value) 
     > call errquit('cpsp_v_nonlocal_euv_2:out of stack',0,MA_ERR)

      
      
*     ***********************
*     **** calculate Bus ****
*     ***********************
      call dcopy(9,0.0d0,0,Bus,1)
      do ii=1,ion_nion()
        ia=ion_katm(ii)
        nproj = int_mb(nprj(1)+ia-1)

        if (nproj.gt.0) then

        do nb=1,nbrillq
          call dcopy(9,0.0d0,0,Aus,1)

          psi1_shift = cpsi_data_get_chnk(psi1_tag,nb)
          if (occ1_tag.gt.0) 
     >       occ1_shift = cpsi_data_get_chnk(occ1_tag,nb)
          call Cram_npack(nb,npack)
          kx = brillioun_k(1,nb)
          ky = brillioun_k(2,nb)
          kz = brillioun_k(3,nb)
          weight = brillioun_weight(nb)
          call C3dB_t_Copy(1,dbl_mb(c_G_indx(1)),dbl_mb(G(1,1)))
          call C3dB_t_Copy(1,dbl_mb(c_G_indx(2)),dbl_mb(G(1,2)))
          call C3dB_t_Copy(1,dbl_mb(c_G_indx(3)),dbl_mb(G(1,3)))
          call Cram_r_pack(nb,dbl_mb(G(1,1)))
          call Cram_r_pack(nb,dbl_mb(G(1,2)))
          call Cram_r_pack(nb,dbl_mb(G(1,3)))
          call daxpy(npack,1.0d0,kx,0,dbl_mb(G(1,1)),1)
          call daxpy(npack,1.0d0,ky,0,dbl_mb(G(1,2)),1)
          call daxpy(npack,1.0d0,kz,0,dbl_mb(G(1,3)),1)

*         **** structure factor and local pseudopotential ****
          call cstrfac_pack(nb,ii,dcpl_mb(exi(1)))
          call cstrfac_k(ii,nb,cxr)
          call zscal(npack,cxr,dcpl_mb(exi(1)),1)

  
*       *********************************************
*       **** calculate F^(lm)_I = <psi|vnl(nlm)> ****
*       *********************************************
        do l=1,nproj
           shift = cpsp_projector_get_ptr(int_mb(vnl(1)+ia-1),nb,l)
           l_prj = int_mb(l_projector(1)+(l-1)+(ia-1)*jmmax_max)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccc change this for actinides where we might have g's ccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#ifdef GCC4
           k = iand(l_prj,1)
#else
           k = and(l_prj,1)
#endif
           sd_function = (k.eq.0)


*          *** current function is s or d ****
           if (sd_function) then
              call Cram_rc_Mul(nb,dbl_mb(shift),
     >                        dcpl_mb(exi(1)),
     >                        dcpl_mb(vtmp(1)))
*          *** current function is p or f ****
           else
              call Cram_irc_Mul(nb,dbl_mb(shift),
     >                        dcpl_mb(exi(1)),
     >                        dcpl_mb(vtmp(1)))
           end if
           call Cram_cc_inzdot(nb,nn,
     >                        dbl_mb(psi1_shift), 
     >                        dcpl_mb(vtmp(1)),
     >                        dcpl_mb(zsw1(1)+(l-1)*nn))
        end do
        call C3dB_Vector_SumAll((2*nn*nproj),dcpl_mb(zsw1(1)))

*       **** zsw2 = Gijl*zsw1 ******
        call Multiply_Gijl_zsw1(nn,
     >                         nproj,
     >                         int_mb(nmax(1)+ia-1),
     >                         int_mb(lmax(1)+ia-1),
     >                         int_mb(n_projector(1)+(ia-1)*jmmax_max),
     >                         int_mb(l_projector(1)+(ia-1)*jmmax_max),
     >                         int_mb(m_projector(1)+(ia-1)*jmmax_max),
     >                         dbl_mb(Gijl(1)+(ia-1)*gij_stride),
     >                         dcpl_mb(zsw1(1)),
     >                         dcpl_mb(zsw2(1)))

        if (ispin.eq.1) call dscal(2*nn*nproj,2.0d0,dcpl_mb(zsw2(1)),1)
        do n=1,nn*nproj
          dcpl_mb(zsw1(1)+n-1) = dconjg(dcpl_mb(zsw2(1)+n-1))
        end do


*       **********************************
*       **** calculate dF^(lm)_I/dhus ****
*       **********************************
        do l=1,nproj
           l_prj = int_mb(l_projector(1)+(l-1)+(ia-1)*jmmax_max)
#ifdef GCC4
           k = iand(l_prj,1)
#else
           k = and(l_prj,1)
#endif
           sd_function = (k.eq.0)

           do s=1,3
              shift2 = cpsp_projector_get_ptr(
     >                  int_mb(dvnl(1)+ia-1),nb,(3*(l-1)+s))
              do u=1,3
                 call Cram_rr_Mul(nb,dbl_mb(shift2),
     >                            dbl_mb(G(1,u)),
     >                            dbl_mb(tmp1(1)))

*                *** current function is s or d ****
                 if (sd_function) then
                   call Cram_rc_Mul(nb,dbl_mb(tmp1(1)),
     >                              dcpl_mb(exi(1)),
     >                              dcpl_mb(vtmp(1)))

*                *** current function is p or f ****
                 else
                    call Cram_irc_Mul(nb,dbl_mb(tmp1(1)),
     >                                dcpl_mb(exi(1)),
     >                                dcpl_mb(vtmp(1)))
                 end if
                 call Cram_cc_nzdot(nb,nn,
     >                              dbl_mb(psi1_shift),
     >                              dcpl_mb(vtmp(1)),
     >                              dcpl_mb(zsw3(1)+(u-1)*nn
     >                                             +(s-1)*nn*3))
             end do
          end do

          if (occ1_tag.gt.0) then
             occ_shift = occ1_shift
             do i=1,nn
               do s=1,3
               do u=1,3
                ctmp = dcpl_mb(zsw1(1)+(i-1)+(l-1)*nn)
     >                *dcpl_mb(zsw3(1)+(i-1)
     >                                +(u-1)*nn
     >                                +(s-1)*nn*3)

                Bus(u,s) = Bus(u,s)
     >            - dbl_mb(occ_shift)*weight*2.0d0*dble(ctmp)/(omega)
               end do
               end do
               occ_shift = occ_shift + 1
             end do
          else
             do i=1,nn
               do s=1,3
               do u=1,3
                ctmp = dcpl_mb(zsw1(1)+(i-1)+(l-1)*nn)
     >                *dcpl_mb(zsw3(1)+(i-1)
     >                                +(u-1)*nn
     >                                +(s-1)*nn*3)

                Bus(u,s) = Bus(u,s) 
     >            - weight*2.0d0*dble(ctmp)/(omega)
                Aus(u,s) = 
     >            - weight*2.0d0*dble(ctmp)/(omega)
               end do
               end do
             end do
          end if

        end do !** l **
        end do  !** nb **

        end if
      end do !** ii **
      call K1dB_Vector_SumAll(9,Bus)

      value =           BA_pop_stack(zsw3(2))
      value = value.and.BA_pop_stack(zsw2(2))
      value = value.and.BA_pop_stack(zsw1(2))
      value = value.and.BA_pop_stack(G(2,3))
      value = value.and.BA_pop_stack(G(2,2))
      value = value.and.BA_pop_stack(G(2,1))
      value = value.and.BA_pop_stack(tmp1(2))
      value = value.and.BA_pop_stack(vtmp(2))
      value = value.and.BA_pop_stack(exi(2))
      if (.not. value) 
     >  call errquit('cpsp_v_nonlocal_euv_2:error popping stack',0,
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
*     *	  	cpsp_stress_read        *
*     *					*
*     ***********************************

      subroutine cpsp_stress_read(fname,
     >                       version,
     >                       nfft,unita,
     >                       npack0,dvl,
     >                       npack1,lmmax,lmmax_max,dvnl_tag,
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
      integer npack1,lmmax,lmmax_max
      integer dvnl_tag
      logical semicore
      real*8 dncore(*)
      real*8 tmp(*)
      real*8 tmp2(*)
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
      integer MASTER,taskid,taskid_k
      parameter(MASTER=0)
      integer i,n,l
      integer msglen
      character*255 full_filename

      real*8 kv(3)
      integer nbrillioun,nb,nbq,pk

      integer  brillioun_nbrillioun,brillioun_nbrillq
      integer  cpsp_projector_alloc
      real*8   brillioun_all_k
      external brillioun_nbrillioun,brillioun_nbrillq
      external cpsp_projector_alloc
      external brillioun_all_k

      call Parallel_taskid(taskid)
      call Parallel3d_taskid_k(taskid_k)
  
*     **** open fname binary file ****
      if (taskid.eq.MASTER) then
         call util_file_name_noprefix(fname,.false.,
     >                             .false.,
     >                       full_filename)
         l = index(full_filename,' ') - 1
         call openfile(5,full_filename,l,'r',l)
         call iread(5,version,1)
         call iread(5,nfft,3)
         call dread(5,unita,9)
         call iread(5,nbrillioun,1)
         ierr = 0
         if (nbrillioun.eq.brillioun_nbrillioun()) then
            do nb=1,nbrillioun
               call dread(5,kv,3)
               if ((brillioun_all_k(1,nb).ne.kv(1)).or.
     >             (brillioun_all_k(2,nb).ne.kv(2)).or.
     >             (brillioun_all_k(3,nb).ne.kv(3))) ierr = 1
            end do
         else
            ierr = 1
            call closefile(5)
         end if
      end if

      msglen = 1
      call Parallel_Brdcst_ivalues(MASTER,msglen,ierr)
      if (ierr.ne.0) then
         return
      end if


*     **** send header data to all processors ****
      msglen = 1
      call Parallel_Brdcst_ivalues(MASTER,msglen,version)
      msglen = 3
      call Parallel_Brdcst_ivalues(MASTER,msglen,nfft)
      msglen = 9
      call Parallel_Brdcst_values(MASTER,msglen,unita)

 
*     *** read in dvl 3d block ***
      call C3dB_r_read(1,5,tmp2,tmp,-1,-1,.true.)
      call Cram_r_pack(0,tmp2)
      call Cram_r_Copy(0,tmp2,dvl)

*     **** read in dvnl 3d blocks ****
      if (lmmax.gt.0) then
      dvnl_tag = cpsp_projector_alloc(brillioun_nbrillq(),
     >                                3*lmmax,npack1)
      do nb=1,brillioun_nbrillioun()
         call K1dB_ktoqp(nb,nbq,pk)
         do n=1,lmmax
         do i=1,3
            call C3dB_r_read(1,5,tmp2,tmp,-1,pk,.true.)
            if (pk.eq.taskid_k) then
               call Cram_r_pack(nbq,tmp2)
               call cpsp_projector_add(dvnl_tag,nbq,(3*(n-1)+i),tmp2)
            end if
         end do
         end do
      end do
      end if

*     **** read in semicore density block ****
      if (semicore) then
         !write(*,*) "reading in semicore block"   !debug
         call C3dB_r_Read(1,5,tmp2,tmp,-1,-1,.true.)
         call Cram_r_pack(0,tmp2)
         call Cram_r_Copy(0,tmp2,dncore)
      end if

*     *** close fname binary file ***
      if (taskid.eq.MASTER) then
c       close(11)
         call closefile(5)
      end if

      ierr = 0
      return
      end

*     ***********************************
*     *					*
*     *	     cpsp_stress_readall 	*
*     *					*
*     ***********************************

      subroutine cpsp_stress_readall()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"      
#include "cpsp_common.fh"
#include "c_semicore_common.fh"



*     **** local variables ****
      integer ngp(3)
      real*8  unita(3,3)
      integer version,nfft3d,npack1,npack0,nbrill
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
      logical      nwpw_filefind
      integer      control_ngrid,brillioun_nbrillioun
      integer      psp_lmmax,cpsp_nprj
      real*8       control_unita
      character*12 control_boundry
      character*4  ion_atom
      external     nwpw_filefind
      external     control_ngrid,brillioun_nbrillioun
      external     psp_lmmax
      external     control_unita
      external     control_boundry
      external     ion_atom
      external     cpsp_nprj      

      call Parallel_taskid(taskid)

      call C3dB_nfft3d(1,nfft3d)
      call Cram_npack(0,npack0)
      call Cram_max_npack(npack1)
      nbrill = brillioun_nbrillioun()

ccc corrected to complex by pjn 11-7-06
      value = BA_push_get(mt_dbl,(4*nfft3d),'tmp',tmp(2),tmp(1))
      if (.not. value) 
     > call errquit('cpsp_stress_readall:out of stack memory',0,MA_ERR)

      value = BA_push_get(mt_dbl,(4*nfft3d),'tmp2',tmp2(2),tmp2(1))
      if (.not. value) 
     > call errquit('cpsp_stress_readall:out of stack memory',0,MA_ERR)

*     **** read pseudopotentials ****
      do ia=1,npsp

*       **** define formatted psp name ****
         element = '     '
         element = ion_atom(ia)
         l = index(element,' ') - 1
         fname = element(1:l)//'.cpp2'
ccccccccccccc possible error here!!!!!! CCCCCCCCCCCCCCCCCCc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c         !nproj = psp_lmmax(ia)
         nproj=cpsp_nprj(ia)
    

*        **** not finished ****
         found = .false.
         do while (.not. found)
         if (nwpw_filefind(fname)) then
            call cpsp_stress_read(fname,
     >                  version,
     >                  ngp,unita,
     >                  npack0,
     >                  dbl_mb(dvl(1) + (ia-1)*npack0),
     >                  npack1,nproj,lmmax_max,
     >                  int_mb(dvnl(1)+ia-1),
     >                  log_mb(semicore(1)+ia),
     >                  dbl_mb(ncore(1) + npack0 + (ia-1)*5*npack0),
     >                  dbl_mb(tmp(1)),dbl_mb(tmp2(1)),
     >                  ierr)


*          **** set semicore(0) ****
           if (log_mb(semicore(1)+ia)) log_mb(semicore(1)) = .true.
           if (ierr.gt.1) go to 9000

*          **************************************************************
*          ***** logic for finding out if psp is correctly formatted ****
*          **************************************************************
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
              if (taskid.eq.MASTER) then
              write(6,*) "pseudopotential is not correctly formatted:",
     >                    fname
              end if
              if ((ierr.eq.0).and.(nproj.gt.0)) then
                 call cpsp_projector_dealloc(int_mb(dvnl(1)+ia-1))
              end if
           end if
           if (correct_box) found = .true.
           if (ierr.eq.1)   then
              found = .false.
              if (taskid.eq.MASTER) then
              write(6,*)
              write(6,*) "pseudopotential is not correctly formatted-",
     >                   "bad brillioun zone:",fname
              end if
           end if

         end if

*        **** generate formatted pseudopotential atom.cpp2 *****
         if (.not.found) then
             call cpsp_stress_formatter_auto(ion_atom(ia))
         end if

         end do !*** do while ****
 
      end do
 9000 value = BA_pop_stack(tmp2(2))
      value = value.and.BA_pop_stack(tmp(2))
      if (.not. value) 
     >  call errquit('cpsp_stress_readall:error popping stack',0,
     &       MA_ERR)
      return
      end
