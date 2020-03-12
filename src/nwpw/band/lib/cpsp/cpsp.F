*
* $Id$
*

#define TCGMSG

*     ****************************************************
*     *                                                  *
*     *                cpsp_init                         *
*     *                                                  *
*     ****************************************************
      subroutine cpsp_init()      
      implicit none

#include "bafdecls.fh"
#include "cpsp_common.fh"
#include "errquit.fh"

*     ***** local variables ****
      logical value
      integer axsize
      integer npack0,npack1,nbrillq,nion

*     **** external functions *****
      logical  two_component_pseudopotential,control_pspspin
      external two_component_pseudopotential,control_pspspin
      integer  ion_nkatm, brillioun_nbrillq,ion_nion
      external ion_nkatm, brillioun_nbrillq,ion_nion
 
      npsp = ion_nkatm()
      nbrillq = brillioun_nbrillq()

      call Cram_npack(0,npack0)
      call Cram_max_npack(npack1)

*     **** allocate projector datastructure ****
      call cpsp_projector_init(5*npsp)

*     ***** allocate local potentitals *****
      value = BA_alloc_get(mt_dbl,(npsp*npack0),'vl',vl(2),vl(1))


*     ***** allocate projector tag lists *****
      value = value.and.
     >        BA_alloc_get(mt_int,npsp,'vnl',vnl(2),vnl(1))
      value = value.and.
     >        BA_alloc_get(mt_int,npsp,'vnlso',vnlso(2),vnlso(1))

*     **** nonlocal coeffiencients *****
      axsize=npsp*gij_stride
      value = value.and.
     >      BA_alloc_get(mt_dbl,axsize,'Gijl',Gijl(2),Gijl(1))
      axsize=npsp*kij_stride
      value = value.and.
     >      BA_alloc_get(mt_dbl,axsize,'Kijl',Kijl(2),Kijl(1))

      value = value.and.
     >        BA_alloc_get(mt_int,(npsp),'nprj',nprj(2),nprj(1))

      axsize=npsp*jmmax_max
      value = value.and.
     >        BA_alloc_get(mt_int,axsize,
     >                     'n_projector',n_projector(2),n_projector(1))
      value = value.and.
     >        BA_alloc_get(mt_int,axsize,
     >                     'l_projector',l_projector(2),l_projector(1))
      value = value.and.
     >        BA_alloc_get(mt_int,axsize,
     >                     'm_projector',m_projector(2),m_projector(1))


      value = value.and.
     >        BA_alloc_get(mt_dbl,(npsp),'zv',zv(2),zv(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,(npsp),'amass',amass(2),amass(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,(npsp*(lmax_max+1)),'rc',rc(2),rc(1))
      value = value.and.
     >        BA_alloc_get(mt_int,(npsp),'lmmax',lmmax(2),lmmax(1))
      value = value.and.
     >        BA_alloc_get(mt_int,(npsp),'lmax',lmax(2),lmax(1))
      value = value.and.
     >        BA_alloc_get(mt_int,(npsp),'locp',locp(2),locp(1))
      value = value.and.
     >        BA_alloc_get(mt_int,(npsp),'nmax',nmax(2),nmax(1))
      value = value.and.
     >        BA_alloc_get(mt_int,(npsp),
     >                     'psp_type',psp_type(2),psp_type(1))

*     **** setup pspspin structure - used for generating antiferromagnetic structures ****
      pspspin = control_pspspin()
      if (pspspin) then
         nion  = ion_nion()
         value = value.and.
     >   BA_alloc_get(mt_log,nion,'pspspin_upions',
     >                pspspin_upions(2),pspspin_upions(1))
         value = value.and.
     >   BA_alloc_get(mt_log,nion,'pspspin_downions',
     >                pspspin_downions(2),pspspin_downions(1))
         if (.not. value)
     >   call errquit('psp_init:out of heap memory',0, MA_ERR)

         call control_set_pspspin(pspspin_upscale,pspspin_downscale,
     >                            pspspin_upl,pspspin_downl,
     >                            nion,
     >                            log_mb(pspspin_upions(1)),
     >                            log_mb(pspspin_downions(1)))
      end if

      if (.not. value) 
     >   call errquit('cpsp_init:out of heap memory',0, MA_ERR)

      call dcopy(npsp*npack0,0.0d0,0,dbl_mb(vl(1)),1)
      call dcopy(npsp,       0.0d0,0,dbl_mb(zv(1)),1)
      call dcopy(npsp,       0.0d0,0,dbl_mb(amass(1)),1)
      call dcopy(npsp*(lmax_max+1),    0.0d0,0,dbl_mb(rc(1)),1)

*     **** allocate semicore data ****
      call c_semicore_init()
      return
      end

*     ****************************************************
*     *                                                  *
*     *                cpsp_proj_init                    *
*     *                                                  *
*     ****************************************************
      subroutine cpsp_proj_init()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"
#include "cpsp_common.fh"

*     ***** local variables ****
      logical value
      integer axsize,npack1

*     ***** external functions *****
      integer  cpsp_nprj_max
      external cpsp_nprj_max

      nprj_max = cpsp_nprj_max()
      call Cram_max_npack(npack1)

      axsize=npack1*nprj_max
      value = BA_alloc_get(mt_dcpl,axsize,'prjtmp',prjtmp(2),prjtmp(1))
      if (.not. value) 
     >   call errquit('cpsp_proj_init:out of heap memory',0,MA_ERR)

      return
      end

*     ******************************************
*     *                                        *
*     *             cpsp_end                   *
*     *                                        *
*     ******************************************
      subroutine cpsp_end()      
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"
#include "cpsp_common.fh"

      logical value

*     **** external functions ****


*     **** deallocate projector data ****
      call cpsp_projector_end()

*     **** deallocate semicore data ****
      call c_semicore_end()

      value = BA_free_heap(prjtmp(2))
      value = value.and.BA_free_heap(vl(2))
      value = value.and.BA_free_heap(vnlso(2))
      value = value.and.BA_free_heap(vnl(2))
      value = value.and.BA_free_heap(Kijl(2))
      value = value.and.BA_free_heap(Gijl(2))
      value = value.and.BA_free_heap(nprj(2))
      value = value.and.BA_free_heap(n_projector(2))
      value = value.and.BA_free_heap(l_projector(2))
      value = value.and.BA_free_heap(m_projector(2))
      value = value.and.BA_free_heap(zv(2))
      value = value.and.BA_free_heap(amass(2))
      value = value.and.BA_free_heap(rc(2))
      value = value.and.BA_free_heap(lmmax(2))
      value = value.and.BA_free_heap(lmax(2))
      value = value.and.BA_free_heap(locp(2))
      value = value.and.BA_free_heap(nmax(2))
      value = value.and.BA_free_heap(psp_type(2))
      if (pspspin) then
         value = value.and.BA_free_heap(pspspin_upions(2))
         value = value.and.BA_free_heap(pspspin_downions(2))
      end if
      if (.not. value)
     >  call errquit('cpsp_end:error freeing heap memory',0,MA_ERR)

      return
      end


*     ***********************************
*     *					*
*     *	 	   cpsp_zv		*
*     *					*
*     ***********************************
      real*8 function cpsp_zv(ia)
      implicit none
      integer ia

#include "bafdecls.fh"
#include "cpsp_common.fh"


      cpsp_zv = dbl_mb(zv(1)+ia-1)
      return
      end


*     ***********************************
*     *					*
*     *	 	   cpsp_amass		*
*     *					*
*     ***********************************
      real*8 function cpsp_amass(ia)
      implicit none
      integer ia

#include "bafdecls.fh"
#include "cpsp_common.fh"

      cpsp_amass = dbl_mb(amass(1)+ia-1)
      return
      end


*     ***********************************
*     *					*
*     *	 	   cpsp_rc		*
*     *					*
*     ***********************************
      real*8 function cpsp_rc(i,ia)
      implicit none
      integer i,ia

#include "bafdecls.fh"
#include "cpsp_common.fh"

      cpsp_rc = dbl_mb(rc(1) + i + (lmax_max+1)*(ia-1))
      return
      end

*     ***********************************
*     *					*
*     *	 	   cpsp_atom		*
*     *					*
*     ***********************************
      character*2 function cpsp_atom(ia)
      implicit none
      integer  ia

#include "cpsp_common.fh"

      cpsp_atom = atom(ia)
      return
      end


*     ***********************************
*     *                                 *
*     *            psp_comment          *
*     *                                 *
*     ***********************************
      character*(*) function cpsp_comment(ia)
      implicit none
      integer  ia

#include "cpsp_common.fh"


      cpsp_comment = comment(ia)
      return
      end




*     ***********************************
*     *					*
*     *	 	   cpsp_lmmax		*
*     *					*
*     ***********************************
      integer function cpsp_lmmax(ia)
      implicit none
      integer  ia

#include "bafdecls.fh"
#include "cpsp_common.fh"


      cpsp_lmmax = int_mb(lmmax(1)+ia-1)
      return
      end

*     ***********************************
*     *                                 *
*     *            cpsp_nprj             *
*     *                                 *
*     ***********************************
      integer function cpsp_nprj(ia)
      implicit none
      integer  ia

#include "bafdecls.fh"
#include "cpsp_common.fh"

      cpsp_nprj  = int_mb(nprj(1)+ia-1)
      return
      end

*     ***********************************
*     *                                 *
*     *            cpsp_nprj_max        *
*     *                                 *
*     ***********************************
      integer function cpsp_nprj_max()
      implicit none

#include "bafdecls.fh"
#include "cpsp_common.fh"

      integer ia,nprjmax,nprjtmp

      nprjmax = 0
      do ia=1,npsp
         nprjtmp = (int_mb(nprj(1)+ia-1))
         if (int_mb(psp_type(1)+ia-1).eq.7) nprjtmp = 2*nprjtmp
         if (nprjtmp.gt.nprjmax) nprjmax = nprjtmp
      end do

      cpsp_nprj_max  = nprjmax
      return
      end

*     ***********************************
*     *                                 *
*     *            cpsp_psp_type        *
*     *                                 *
*     ***********************************
      integer function cpsp_psp_type(ia)
      implicit none
      integer  ia

#include "bafdecls.fh"
#include "cpsp_common.fh"

      cpsp_psp_type  = int_mb(psp_type(1)+ia-1)
      return
      end


*     ***********************************
*     *					*
*     *	 	   cpsp_lmax		*
*     *					*
*     ***********************************
      integer function cpsp_lmax(ia)
      implicit none
      integer  ia

#include "bafdecls.fh"
#include "cpsp_common.fh"


      cpsp_lmax = int_mb(lmax(1)+ia-1)
      return
      end

*     ***********************************
*     *					*
*     *	 	   cpsp_locp		*
*     *					*
*     ***********************************
      integer function cpsp_locp(ia)
      implicit none
      integer  ia

#include "bafdecls.fh"
#include "cpsp_common.fh"


      cpsp_locp = int_mb(locp(1)+ia-1)
      return
      end

*     ***********************************
*     *					*
*     *	 	   cpsp_npsp		*
*     *					*
*     ***********************************
      integer function cpsp_npsp()
      implicit none

#include "cpsp_common.fh"


      cpsp_npsp = npsp
      return
      end


*     ***********************************
*     *					*
*     *	 	   cpsp_v_local         *
*     *					*
*     ***********************************
      subroutine cpsp_v_local(vl_out,move,dng,fion)
      implicit none
      complex*16 vl_out(*)
      logical    move
      complex*16 dng(*)
      real*8     fion(3,*)

#include "bafdecls.fh"
#include "cpsp_common.fh"
#include "errquit.fh"

*     *** local variables ***
      integer nfft3d,npack0
      integer i,ii,ia
      integer exi(2),vtmp(2),xtmp(2),G(3)
      integer Gx(2),Gy(2),Gz(2)
      logical value

*     **** external functions ****
      integer  c_G_indx,ion_nion,ion_katm
      external c_G_indx,ion_nion,ion_katm

      call nwpw_timing_start(5)
      call C3dB_nfft3d(1,nfft3d)
      call Cram_npack(0,npack0)

      value = BA_push_get(mt_dcpl,npack0,'exi', exi(2), exi(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,nfft3d,'vtmp',vtmp(2),vtmp(1))
      if (.not. value) call errquit('cpsp_v_local: pushing stack',0,
     &       MA_ERR)
      
      if (move) then
        value = BA_push_get(mt_dbl, nfft3d,'xtmp',xtmp(2),xtmp(1))
        value = value.and.
     >          BA_push_get(mt_dbl, nfft3d,'Gx',Gx(2),Gx(1))
        value = value.and.
     >          BA_push_get(mt_dbl, nfft3d,'Gy',Gy(2),Gy(1))
        value = value.and.
     >          BA_push_get(mt_dbl, nfft3d,'Gz',Gz(2),Gz(1))
        if (.not. value) call errquit('cpsp_v_local: pushing stack',1,
     &       MA_ERR)
        G(1)  = c_G_indx(1)
        G(2)  = c_G_indx(2)
        G(3)  = c_G_indx(3)

*       **** define Gx,Gy and Gz in packed space ****
        call C3dB_t_Copy(1,dbl_mb(G(1)),dbl_mb(Gx(1)))
        call C3dB_t_Copy(1,dbl_mb(G(2)),dbl_mb(Gy(1)))
        call C3dB_t_Copy(1,dbl_mb(G(3)),dbl_mb(Gz(1)))
        call Cram_r_pack(0,dbl_mb(Gx(1)))
        call Cram_r_pack(0,dbl_mb(Gy(1)))
        call Cram_r_pack(0,dbl_mb(Gz(1)))
      end if

      call dcopy((2*npack0),0.0d0,0,vl_out,1)
      do ii=1,ion_nion()
        ia=ion_katm(ii)

*       **** structure factor and local pseudopotential ****
        call cstrfac_pack(0,ii,dcpl_mb(exi(1)))
       
*       **** add to local psp ****
        call Cram_rc_Mul(0,dbl_mb(vl(1)+npack0*(ia-1)),
     >                   dcpl_mb(exi(1)),
     >                   dcpl_mb(vtmp(1)))
c        call Cram_cc_Sum(0,vl_out,dcpl_mb(vtmp(1)),vl_out)
        call Cram_cc_Sum2(0,dcpl_mb(vtmp(1)),vl_out)

        if (move) then
       
          do i=1,npack0
             dbl_mb(xtmp(1)+i-1) 
     >          = dimag(dng(i))* dble(dcpl_mb(vtmp(1)+i-1))
     >           - dble(dng(i))*dimag(dcpl_mb(vtmp(1)+i-1))
          end do
         call Cram_rr_dot(0,dbl_mb(Gx(1)),dbl_mb(xtmp(1)),fion(1,ii))
         call Cram_rr_dot(0,dbl_mb(Gy(1)),dbl_mb(xtmp(1)),fion(2,ii))
         call Cram_rr_dot(0,dbl_mb(Gz(1)),dbl_mb(xtmp(1)),fion(3,ii))
   
       end if
        

      end do

      if (move) then
        value = BA_pop_stack(Gz(2))
        value = value.and.BA_pop_stack(Gy(2))
        value = value.and.BA_pop_stack(Gx(2))
        value = value.and.BA_pop_stack(xtmp(2))
        if (.not. value) call errquit('cpsp_v_local: popping stack',0,
     &       MA_ERR)
      end if
      value = BA_pop_stack(vtmp(2))
      value = value.and.BA_pop_stack(exi(2))
      if (.not. value) call errquit('cpsp_v_local: popping stack',1,
     &       MA_ERR)

      call nwpw_timing_end(5)
      return 
      end


*     ***********************************
*     *					*
*     *	 	   cpsp_f_vlocal  	*
*     *					*
*     ***********************************

      subroutine cpsp_f_vlocal(dng,fion)
      implicit none
      complex*16 dng(*)
      real*8     fion(3,*)

#include "bafdecls.fh"
#include "cpsp_common.fh"
#include "errquit.fh"


*     *** local variables ***
      integer nfft3d,npack0
      integer i,ii,ia
      integer exi(2),vtmp(2),xtmp(2),G(3)
      integer Gx(2),Gy(2),Gz(2)
      logical value

*     **** external functions ****
      integer  c_G_indx,ion_nion,ion_katm
      external c_G_indx,ion_nion,ion_katm

      call nwpw_timing_start(5)
      
      call C3dB_nfft3d(1,nfft3d)
      call Cram_npack(0,npack0)
      value = BA_push_get(mt_dcpl,npack0,'exi', exi(2), exi(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,nfft3d,'vtmp',vtmp(2),vtmp(1))
      value = value.and.
     >        BA_push_get(mt_dbl, nfft3d,'xtmp',xtmp(2),xtmp(1))
      value = value.and.
     >        BA_push_get(mt_dbl, nfft3d,'Gx',Gx(2),Gx(1))
      value = value.and.
     >        BA_push_get(mt_dbl, nfft3d,'Gy',Gy(2),Gy(1))
      value = value.and.
     >        BA_push_get(mt_dbl, nfft3d,'Gz',Gz(2),Gz(1))
      if (.not. value) call errquit('cpsp_f_vlocal:pushing stack',0,
     &       MA_ERR)
      G(1)  = c_G_indx(1)
      G(2)  = c_G_indx(2)
      G(3)  = c_G_indx(3)

*     **** define Gx,Gy and Gz in packed space ****
      call C3dB_t_Copy(1,dbl_mb(G(1)),dbl_mb(Gx(1)))
      call C3dB_t_Copy(1,dbl_mb(G(2)),dbl_mb(Gy(1)))
      call C3dB_t_Copy(1,dbl_mb(G(3)),dbl_mb(Gz(1)))
      call Cram_r_pack(0,dbl_mb(Gx(1)))
      call Cram_r_pack(0,dbl_mb(Gy(1)))
      call Cram_r_pack(0,dbl_mb(Gz(1)))

      do ii=1,ion_nion()
        ia=ion_katm(ii)

*       **** structure factor and local pseudopotential ****
c        call cstrfac(ii,dcpl_mb(exi(1)))
c        call Cram_c_pack(0,dcpl_mb(exi(1)))
        call cstrfac_pack(0,ii,dcpl_mb(exi(1)))
        
*       **** add to local psp ****
        call Cram_rc_Mul(0,dbl_mb(vl(1)+npack0*(ia-1)),
     >                   dcpl_mb(exi(1)),
     >                   dcpl_mb(vtmp(1)))

        do i=1,npack0
           dbl_mb(xtmp(1)+i-1) 
     >        = dimag(dng(i))* dble(dcpl_mb(vtmp(1)+i-1))
     >         - dble(dng(i))*dimag(dcpl_mb(vtmp(1)+i-1))
        end do

        call Cram_rr_dot(0,dbl_mb(Gx(1)),dbl_mb(xtmp(1)),fion(1,ii))
        call Cram_rr_dot(0,dbl_mb(Gy(1)),dbl_mb(xtmp(1)),fion(2,ii))
        call Cram_rr_dot(0,dbl_mb(Gz(1)),dbl_mb(xtmp(1)),fion(3,ii))
 
      end do
      value = BA_pop_stack(Gz(2))
      value = value.and.BA_pop_stack(Gy(2))
      value = value.and.BA_pop_stack(Gx(2))
      value = value.and.BA_pop_stack(xtmp(2))
      value = value.and.BA_pop_stack(vtmp(2))
      value = value.and.BA_pop_stack(exi(2))
      if (.not. value) call errquit('cpsp_f_vlocal:popping stack',1,
     &       MA_ERR)

      call nwpw_timing_end(5)
      return 
      end



*     ***********************************
*     *					*
*     *	 	 cpsp_v_nonlocal  	*
*     *					*
*     ***********************************

      subroutine cpsp_v_nonlocal(ispin,ne,psi1_tag,psi2_tag,move,fion)
      implicit none
      integer    ispin,ne(2)
      integer    psi1_tag
      integer    psi2_tag
      logical move
      real*8 fion(3,*)

#include "bafdecls.fh"
#include "cpsp_common.fh"
#include "errquit.fh"

*     *** local variables ***
      complex*16 one,mone,ione
      integer nfft3d,G(3),npack1,npack
      integer ii,ia,l,n,nn,nb,nbq,neall,nbrillq,nion
      integer shift,l_prj,nproj
      integer psi1_shift,psi2_shift,psi_shift,nshift
      integer occ_tag,occ_shift,occ1_shift
      real*8  omega,weight,scal,wfrac
      complex*16 cxr
      integer exi(2),xtmp(2),zsw1(2),zsw2(2),sum(2),ftmp(2)
      integer Gx(2),Gy(2),Gz(2),shfts,shftd
      logical value,sd_function

*     **** external functions ****
      logical  is_sORd,cpsi_spin_orbit
      integer  ion_nion,ion_katm,c_G_indx,Pneb_nbrillq
      integer  cpsp_projector_get_ptr,cpsi_data_get_chnk
      integer  cpsi_data_get_next
      real*8   lattice_omega,brillioun_weight
      external is_sORd,cpsi_spin_orbit
      external ion_nion,ion_katm,c_G_indx,Pneb_nbrillq
      external cpsp_projector_get_ptr,cpsi_data_get_chnk
      external cpsi_data_get_next
      external lattice_omega,brillioun_weight

      one  = dcmplx( 1.0d0,0.0d0)
      mone = dcmplx(-1.0d0,0.0d0)
      ione = dcmplx( 0.0d0,1.0d0)
      call nwpw_timing_start(6)

      
*     **** allocate local memory ****
      nn = ne(1)+ne(2)
      neall = ne(1)+ne(2)
      nion  = ion_nion()
      nbrillq = Pneb_nbrillq()
      call C3dB_nfft3d(1,nfft3d)
      call Cram_max_npack(npack1)
      nshift = 2*npack1

      value = BA_push_get(mt_dcpl,npack1,'exi',exi(2),exi(1))
      value = value.and.BA_push_get(mt_dcpl,nn*nprj_max,
     >                              'zsw1',zsw1(2),zsw1(1))
      value = value.and.BA_push_get(mt_dcpl,nn*nprj_max,
     >                              'zsw2',zsw2(2),zsw2(1))
      if (.not.value) call errquit('cpsp_v_nonlocal:pushing stack',0,
     >       MA_ERR)

ccccccccccccccc if move we do this ccccccccccccccccccccccccc
      if (move) then
        occ_tag = cpsi_data_get_next(psi1_tag)
        value = BA_push_get(mt_dbl, nfft3d,'xtmp',xtmp(2),xtmp(1))
        value = value.and.
     >          BA_push_get(mt_dbl, nfft3d,'Gx',Gx(2),Gx(1))
        value = value.and.
     >          BA_push_get(mt_dbl, nfft3d,'Gy',Gy(2),Gy(1))
        value = value.and.
     >          BA_push_get(mt_dbl, nfft3d,'Gz',Gz(2),Gz(1))
        value = value.and.
     >        BA_push_get(mt_dbl,3*nn,'sum',sum(2),sum(1))
        value = value.and.
     >        BA_push_get(mt_dbl,3*nion,'ftmp',ftmp(2),ftmp(1))
        if (.not.value) call errquit('cpsp_v_nonlocal:pushing stack',1,
     &       MA_ERR)
        G(1)  = c_G_indx(1)
        G(2)  = c_G_indx(2)
        G(3)  = c_G_indx(3)
        call dcopy(3*nion,0.0d0,0,dbl_mb(ftmp(1)),1)
      end if
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      omega = lattice_omega()
      scal  = 1.0d0/omega

      do nbq=1,nbrillq
        psi1_shift = cpsi_data_get_chnk(psi1_tag,nbq)
        psi2_shift = cpsi_data_get_chnk(psi2_tag,nbq)
        call Cram_npack(nbq,npack)

      do 500 ii=1,nion
        ia=ion_katm(ii)

        nproj = int_mb(nprj(1)+ia-1)
        if (nproj.gt.0) then 

          if (int_mb(psp_type(1)+ia-1).eq.7) then

            call cpsp_v_nonlocal_rel(nbq,ii,ia,nproj,npack1,nbrillq,
     >      nn,ne,
     >      dcpl_mb(exi(1)),dcpl_mb(prjtmp(1)),
     >      dcpl_mb(zsw1(1)),
     >      dbl_mb(psi1_shift),dbl_mb(psi2_shift),
     >      dbl_mb(Kijl(1)+(ia-1)*kij_stride),scal)

            if (move) then

*              **** define Gx,Gy and Gz in packed space ****
               call C3dB_t_Copy(1,dbl_mb(G(1)),dbl_mb(Gx(1)))
               call C3dB_t_Copy(1,dbl_mb(G(2)),dbl_mb(Gy(1)))
               call C3dB_t_Copy(1,dbl_mb(G(3)),dbl_mb(Gz(1)))
               call Cram_r_pack(nbq,dbl_mb(Gx(1)))
               call Cram_r_pack(nbq,dbl_mb(Gy(1)))
               call Cram_r_pack(nbq,dbl_mb(Gz(1)))

              shftd=nproj*npack1
              shfts=ne(1)*npack1*2
              weight = brillioun_weight(nbq)
              if (occ_tag.gt.0) 
     >           occ1_shift = cpsi_data_get_chnk(occ_tag,nbq)
              do l=1,nproj          
                 psi_shift = psi1_shift
                 do n=1,ne(1)

*                   **** routine should be vectorized!!!****
                    call Cram_zccr_Multiply2(nbq,
     >                     dcpl_mb(zsw1(1)+(l-1)*ne(1)+n-1),
     >                     dbl_mb(psi_shift),
     >                     dcpl_mb(prjtmp(1)+(l-1)*npack1),
     >                     dbl_mb(xtmp(1)))
                    call Cram_zccr_Multiply2Add(nbq,
     >                     dcpl_mb(zsw1(1)+(l-1)*ne(1)+n-1),
     >                     dbl_mb(psi_shift+shfts),
     >                     dcpl_mb(prjtmp(1)+(l-1)*npack1+shftd),
     >                     dbl_mb(xtmp(1)))
                    psi_shift = psi_shift + nshift

                    call Cram_rr_idot(nbq,dbl_mb(Gx(1)),dbl_mb(xtmp(1)),
     >                                dbl_mb(sum(1)+3*(n-1)))
                    call Cram_rr_idot(nbq,dbl_mb(Gy(1)),dbl_mb(xtmp(1)),
     >                                dbl_mb(sum(1)+1+3*(n-1)))
                    call Cram_rr_idot(nbq,dbl_mb(Gz(1)),dbl_mb(xtmp(1)),
     >                                dbl_mb(sum(1)+2+3*(n-1)))
                 end do !**n**

                 call C3dB_Vector_Sumall((3*ne(1)),dbl_mb(sum(1)))

                 wfrac = 1.0d0
                 if (occ_tag.gt.0) occ_shift = occ1_shift
                 do n=1,ne(1)
                    if (occ_tag.gt.0) then
                       wfrac = dbl_mb(occ_shift)
                       occ_shift = occ_shift+1
                    end if
                    dbl_mb(ftmp(1)+3*(ii-1)) = 
     >                         dbl_mb(ftmp(1)+3*(ii-1))
     >                         + 2.0d0*weight*wfrac
     >                                *dbl_mb(sum(1)+3*(n-1))
                    dbl_mb(ftmp(1)+3*(ii-1)+1) = 
     >                         dbl_mb(ftmp(1)+3*(ii-1)+1)
     >                         + 2.0d0*weight*wfrac
     >                                *dbl_mb(sum(1)+1+3*(n-1))
                    dbl_mb(ftmp(1)+3*(ii-1)+2) = 
     >                         dbl_mb(ftmp(1)+3*(ii-1)+2)
     >                         + 2.0d0*weight*wfrac
     >                                *dbl_mb(sum(1)+2+3*(n-1))
                 end do !** ne(1) **

              end do !** l **
            end if !** move **
            goto 500

          end if

*       **** structure factor pseudopotential ****
            call cstrfac_pack(nbq,ii,dcpl_mb(exi(1)))
            call cstrfac_k(ii,nbq,cxr)
c             call Cram_c_ZMul(nbq,cxr,dcpl_mb(exi(1)),dcpl_mb(exi(1)))
            call zscal(npack,cxr,dcpl_mb(exi(1)),1)


*          **** generate zsw1's and projectors ****
            do l=1,nproj

              shift = cpsp_projector_get_ptr(
     >                     int_mb(vnl(1)+ia-1),nbq,l)
              l_prj = int_mb(l_projector(1)+(l-1)
     >                      +(ia-1)*jmmax_max)

              sd_function=.true.
              if (mod(l_prj,2).ne.0) then
                sd_function=.false.
              end if

*             **** phase factor does not matter therefore ****
*             **** (-i)^l is the same as (i)^l in the     ****
*             **** Rayleigh scattering formula            ****

c             *** current function is s or d ****
              if (sd_function) then
                 call Cram_rc_Mul(nbq,dbl_mb(shift),
     >                               dcpl_mb(exi(1)),
     >                               dcpl_mb(prjtmp(1)+(l-1)*npack1))
*             *** current function is p or f ****
              else
                 call Cram_irc_Mul(nbq,dbl_mb(shift),
     >                                dcpl_mb(exi(1)),
     >                                dcpl_mb(prjtmp(1)+(l-1)*npack1))
              end if

*             **** compute nnXnproj matrix zsw1 = <psi1|prj> ****
              call Cram_cc_inzdot(nbq,nn,
     >                      dbl_mb(psi1_shift),
     >                      dcpl_mb(prjtmp(1)+(l-1)*npack1),
     >                      dcpl_mb(zsw1(1)+(l-1)*nn))

*           ***** scale psp by factor - used for generating antiferromagnetic structures ****
*           **** nwchem input: pspspin up/down scale l ion_numbers                       ****
            if (pspspin) then
               if (log_mb(pspspin_upions(1)+ii-1).and.
     >            (l_prj.eq.pspspin_upl))
     >            call dscal(2*ne(1),pspspin_upscale,
     >                       dcpl_mb(zsw1(1)+(l-1)*nn),1)
               if (log_mb(pspspin_downions(1)+ii-1).and.
     >            (l_prj.eq.pspspin_downl))
     >            call dscal(2*ne(2),pspspin_downscale,
     >                       dcpl_mb(zsw1(1)+(l-1)*nn+ne(1)),1)
            end if

           end do !**l**

           call C3dB_Vector_SumAll((2*nn*nproj),dcpl_mb(zsw1(1)))


*          **** zsw2 = Gijl*zsw1 ******
           call Multiply_Gijl_zsw1(nn,
     >                         nproj,
     >                         int_mb(nmax(1)+ia-1),
     >                         int_mb(lmax(1)+ia-1),
     >                         int_mb(n_projector(1)
     >                                + (ia-1)*jmmax_max),
     >                         int_mb(l_projector(1)
     >                                + (ia-1)*jmmax_max),
     >                         int_mb(m_projector(1)
     >                                + (ia-1)*jmmax_max),
     >                         dbl_mb(Gijl(1)
     >                         + (ia-1)*gij_stride),
     >                         dcpl_mb(zsw1(1)),
     >                         dcpl_mb(zsw2(1)))


*          **** do Kleinman-Bylander Multiplication ****
           call dscal(2*nn*nproj,scal,dcpl_mb(zsw2(1)),1)
           call ZGEMM('N','C',npack,nn,nproj,
     >                mone,
     >                dcpl_mb(prjtmp(1)), npack1,
     >                dcpl_mb(zsw2(1)),   nn,
     >                one,
     >                dbl_mb(psi2_shift), npack1)


           if (move) then
              weight = brillioun_weight(nbq)
              if (occ_tag.gt.0) 
     >           occ1_shift = cpsi_data_get_chnk(occ_tag,nbq)
              if (ispin.eq.1) 
     >           call dscal(2*nn*nproj,2.0d0,dcpl_mb(zsw2(1)),1)
          
              do l=1,nproj
          
                 psi_shift = psi1_shift
                 do n=1,nn

*                   **** routine should be vectorized!!!****
                    call Cram_zccr_Multiply2(nbq,
     >                                 dcpl_mb(zsw2(1)+(l-1)*nn+n-1),
     >                                 dbl_mb(psi_shift),
     >                                 dcpl_mb(prjtmp(1)+(l-1)*npack1),
     >                                 dbl_mb(xtmp(1)))
                    psi_shift = psi_shift + nshift
c                    do i=1,npack
c                       ctmp = psi1(i+(n-1)*npack1
c     >                              +(nb-1)*neall*npack1)
c     >                      *dconjg(dcpl_mb(prjtmp(1)+(l-1)*npack1+i-1))
c     >                      *dconjg(dcpl_mb(zsw2(1)+(l-1)*nn+n-1))
c                       dbl_mb(xtmp(1)+i-1) = dimag(ctmp)
c                    end do


*                   **** define Gx,Gy and Gz in packed space ****
                    call C3dB_t_Copy(1,dbl_mb(G(1)),dbl_mb(Gx(1)))
                    call C3dB_t_Copy(1,dbl_mb(G(2)),dbl_mb(Gy(1)))
                    call C3dB_t_Copy(1,dbl_mb(G(3)),dbl_mb(Gz(1)))
                    call Cram_r_pack(nbq,dbl_mb(Gx(1)))
                    call Cram_r_pack(nbq,dbl_mb(Gy(1)))
                    call Cram_r_pack(nbq,dbl_mb(Gz(1)))

                    call Cram_rr_idot(nbq,dbl_mb(Gx(1)),dbl_mb(xtmp(1)),
     >                                dbl_mb(sum(1)+3*(n-1)))
                    call Cram_rr_idot(nbq,dbl_mb(Gy(1)),dbl_mb(xtmp(1)),
     >                                dbl_mb(sum(1)+1+3*(n-1)))
                    call Cram_rr_idot(nbq,dbl_mb(Gz(1)),dbl_mb(xtmp(1)),
     >                                dbl_mb(sum(1)+2+3*(n-1)))
                 end do !**n**

                 call C3dB_Vector_Sumall(3*(nn),dbl_mb(sum(1)))

                 wfrac = 1.0d0
                 if (occ_tag.gt.0) occ_shift = occ1_shift
                 do n=1,nn
                    if (occ_tag.gt.0) then
                       wfrac = dbl_mb(occ_shift)
                       occ_shift = occ_shift+1
                    end if
                    dbl_mb(ftmp(1)+3*(ii-1)) = 
     >                         dbl_mb(ftmp(1)+3*(ii-1))
     >                         + 2.0d0*weight*wfrac
     >                                *dbl_mb(sum(1)+3*(n-1))
                    dbl_mb(ftmp(1)+3*(ii-1)+1) = 
     >                         dbl_mb(ftmp(1)+3*(ii-1)+1)
     >                         + 2.0d0*weight*wfrac
     >                                *dbl_mb(sum(1)+1+3*(n-1))
                    dbl_mb(ftmp(1)+3*(ii-1)+2) = 
     >                         dbl_mb(ftmp(1)+3*(ii-1)+2)
     >                         + 2.0d0*weight*wfrac
     >                                *dbl_mb(sum(1)+2+3*(n-1))
                 end do !** nn **

              end do !** l **
           end if !** move **
       
        end if !** nproj>0**

 500  continue
      end do !** nbq **

      if (move) then
        call K1dB_Vector_SumAll(3*nion,dbl_mb(ftmp(1)))
        call daxpy(3*nion,1.0d0,dbl_mb(ftmp(1)),1,fion,1)
        value = BA_pop_stack(ftmp(2))
        value = value.and.BA_pop_stack(sum(2))
        value = value.and.BA_pop_stack(Gz(2))
        value = value.and.BA_pop_stack(Gy(2))
        value = value.and.BA_pop_stack(Gx(2))
        value = value.and.BA_pop_stack(xtmp(2))
        if (.not.value) call errquit('cpsp_v_nonlocal:popping stack',2,
     &       MA_ERR)
      end if

      value =           BA_pop_stack(zsw2(2))
      value = value.and.BA_pop_stack(zsw1(2))
      value = value.and.BA_pop_stack(exi(2))
      if (.not.value) call errquit('cpsp_v_nonlocal:popping stack',3,
     &       MA_ERR)

      call nwpw_timing_end(6)

      return 
      end

*     ************************************************
*     *                                              *
*     *              Multiply_Gijl_zsw1              *
*     *                                              *
*     ************************************************
cccccccccccccccccccc
c
cccccccccccccccccccc
      subroutine Multiply_Gijl_zsw1(nn,nprj,nmax,lmax,
     >                             n_prj,l_prj,m_prj,
     >                             G,
     >                             zsw1,zsw2)
      implicit none
      integer nn
      integer nprj,nmax,lmax
      integer n_prj(nprj)
      integer l_prj(nprj)
      integer m_prj(nprj)
      real*8  G(nmax,nmax,0:lmax)
      complex*16 zsw1(nn,nprj)
      complex*16 zsw2(nn,nprj)

      !**** local variables ****
      integer a,b,na,nb,la,lb,ma,mb


      call dcopy(2*nn*nprj,0.0d0,0,zsw2,1)
      do b=1,nprj
         lb = l_prj(b)
         mb = m_prj(b)

         do a=1,nprj
            la = l_prj(a)
            ma = m_prj(a)

            if ((la.eq.lb).and.(ma.eq.mb)) then
              na = n_prj(a)
              nb = n_prj(b)
              call daxpy(2*nn,G(nb,na,la),zsw1(1,a),1,zsw2(1,b),1)
            end if

         end do
      end do
      return
      end


*     ***********************************
*     *					*
*     *	 	 cpsp_v_nonlocal_orb  	*
*     *					*
*     ***********************************

      subroutine cpsp_v_nonlocal_orb(nb,orb1,orb2)
      implicit none
      integer    nb
      complex*16 orb1(*)
      complex*16 orb2(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "cpsp_common.fh"


*     *** local variables ***
      complex*16 one,mone,ione
c      parameter  (one=(1.0d0,0.0d0), mone=(-1.0d0,0.0d0))
c      parameter  (ione=(0.0d0,1.0d0))
      integer nfft3d,npack1,npack
      integer ii,ia,l,ne1
      integer shift,l_prj,nproj
      real*8  omega,scal
      complex*16 cxr
      integer exi(2),zsw1(2),zsw2(2)
      logical value,sd_function

*     **** external functions ****
      logical  is_sORd,cpsi_spin_orbit
      integer  ion_nion,ion_katm
      integer  cpsp_projector_get_ptr
      integer  cpsi_neq
      real*8   lattice_omega
      external is_sORd,cpsi_spin_orbit
      external ion_nion,ion_katm
      external cpsp_projector_get_ptr
      external lattice_omega
      external cpsi_neq
      one=dcmplx(1.0d0,0.0d0)
      mone=dcmplx(-1.0d0,0.d0)
      ione=dcmplx(0.0d0,1.d0)

      if (cpsi_spin_orbit()) then
        call cpsp_v_nonlocal_orb_2com(nb,orb1,orb2)
        return
      end if

      call nwpw_timing_start(6)
      
*     **** allocate local memory ****
      call C3dB_nfft3d(1,nfft3d)
      call Cram_max_npack(npack1)
c      nbrill = brillioun_nbrillioun()
      ne1=cpsi_neq(1)
      value = BA_push_get(mt_dcpl,npack1,'exi', exi(2), exi(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,nprj_max,
     >                    'zsw1',zsw1(2),zsw1(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,nprj_max,
     >                    'zsw2',zsw2(2),zsw2(1))
      if (.not.value) 
     > call errquit('cpsp_v_nonlocal_orb:pushing stack',0,MA_ERR)


      omega = lattice_omega()
      scal  = 1.0d0/omega

      do 500 ii=1,ion_nion()
         ia=ion_katm(ii)
         nproj = int_mb(nprj(1)+ia-1)

         if (nproj.gt.0) then

            if ( int_mb(psp_type(1)+ia-1) .eq. 7) then
              call cpsp_v_nonlocal_rel_orb(nb,orb1,orb2,
     >             dcpl_mb(zsw1(1)),
     >             dbl_mb(Kijl(1)+(ia-1)*kij_stride),
     >             dcpl_mb(exi(1)),
     >             nfft3d,ia,ii,ne1,
     >             npack1,nproj)
              goto 500
            end if


*           **** structure factor ****
            call Cram_npack(nb,npack)
            call cstrfac_pack(nb,ii,dcpl_mb(exi(1)))
            call cstrfac_k(ii,nb,cxr)
cc            call Cram_c_ZMul(nb,cxr,dcpl_mb(exi(1)),dcpl_mb(exi(1)))
            call zscal(npack,cxr,dcpl_mb(exi(1)),1)

            do l=1,nproj

c              shift = vnl(1)+(l-1) *npack1
c     >                      +(nb-1)*npack1*vnl_stride
c     >                      +(ia-1)*npack1*vnl_stride*nbrill
              shift = cpsp_projector_get_ptr(
     >                     int_mb(vnl(1)+ia-1),nb,l)
              l_prj = int_mb(l_projector(1)+(l-1)
     >                                     +(ia-1)*jmmax_max)

              sd_function=.true.
              if (mod(l_prj,2).ne.0) then
                sd_function=.false.
              end if

*             **** phase factor does not matter therefore ****
*             **** (-i)^l is the same as (i)^l in the     ****
*             **** Rayleigh scattering formula            ****
*             *** current function is s or d ****
              if (sd_function) then
                 call Cram_rc_Mul(nb,dbl_mb(shift),
     >                               dcpl_mb(exi(1)),
     >                               dcpl_mb(prjtmp(1)+(l-1)*npack1))
*             *** current function is p or f ****
              else
                 call Cram_irc_Mul(nb,dbl_mb(shift),
     >                                dcpl_mb(exi(1)),
     >                                dcpl_mb(prjtmp(1)+(l-1)*npack1))
              end if

*             **** compute 1Xnproj matrix zsw1 = <psi1|prj> ****
              call Cram_cc_izdot(nb,
     >                      orb1,
     >                      dcpl_mb(prjtmp(1)+(l-1)*npack1),
     >                      dcpl_mb(zsw1(1)+(l-1)))
            end do !**l**
            call C3dB_Vector_SumAll((2*nproj),dcpl_mb(zsw1(1)))

*           **** zsw2 = Gijl*zsw1 ******
            call Multiply_Gijl_zsw1(1,
     >                         nproj,
     >                         int_mb(nmax(1)+ia-1),
     >                         int_mb(lmax(1)+ia-1),
     >                         int_mb(n_projector(1)
     >                                + (ia-1)*jmmax_max),
     >                         int_mb(l_projector(1)
     >                                + (ia-1)*jmmax_max),
     >                         int_mb(m_projector(1)
     >                                + (ia-1)*jmmax_max),
     >                         dbl_mb(Gijl(1)
     >                         + (ia-1)*gij_stride),
     >                         dcpl_mb(zsw1(1)),
     >                         dcpl_mb(zsw2(1)))

*           **** do Kleinman-Bylander Multiplication ****
            call dscal(2*nproj,scal,dcpl_mb(zsw2(1)),1)
            call ZGEMM('N','C',npack,1,nproj,
     >                 mone,
     >                 dcpl_mb(prjtmp(1)), npack1,
     >                 dcpl_mb(zsw2(1)),   1,
     >                 one,
     >                 orb2, npack1)

         end if !** nproj>0 **
500   continue !** ii ****

      value =           BA_pop_stack(zsw2(2))
      value = value.and.BA_pop_stack(zsw1(2))
      value = value.and.BA_pop_stack(exi(2))
      if (.not.value) 
     > call errquit('cpsp_v_nonlocal_orb:popping stack',3,MA_ERR)

      call nwpw_timing_end(6)

      return 
      end



*     ***********************************
*     *					*
*     *	       cpsp_f_vnonlocal		*
*     *					*
*     ***********************************

      subroutine cpsp_f_vnonlocal(ispin,ne,psi1_tag,fion)
      implicit none
      integer    ispin,ne(2)
      integer    psi1_tag
      real*8 fion(3,*)

#include "bafdecls.fh"
#include "cpsp_common.fh"
#include "errquit.fh"


*     *** local variables ***
      integer nfft3d,G(3),npack1,npack,nbrillq,neall,nion
      integer ii,ia,l,n,nn,nbq,l_prj,nproj,shift
      integer psi1_shift,psi_shift,nshift
      real*8  omega,weight,scal,wfrac
      complex*16 cxr
      integer exi(2),xtmp(2),zsw1(2),zsw2(2),sum(2),ftmp(2)
      integer Gx(2),Gy(2),Gz(2),it
      integer occ1_tag,occ1_shift,occ_shift
      logical value,sd_function

*     **** external functions ****
      logical  is_sORd
      integer  ion_nion,ion_katm,c_G_indx,Pneb_nbrillq
      integer  cpsp_projector_get_ptr,cpsi_data_get_chnk
      integer  cpsi_data_get_next
      real*8   lattice_omega,brillioun_weight
      external is_sORd
      external ion_nion,ion_katm,c_G_indx,Pneb_nbrillq
      external cpsp_projector_get_ptr,cpsi_data_get_chnk
      external cpsi_data_get_next
      external lattice_omega,brillioun_weight

      call nwpw_timing_start(6)
      
*     **** allocate local memory ****
      nn = ne(1)+ne(2)
      neall   = ne(1)+ne(2)
      nion    = ion_nion()
      nbrillq = Pneb_nbrillq()
      call C3dB_nfft3d(1,nfft3d)
      call Cram_max_npack(npack1)
      nshift = 2*npack1
      occ1_tag = cpsi_data_get_next(psi1_tag)

      value = BA_push_get(mt_dcpl,npack1,'exi', exi(2), exi(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,nn*nprj_max,
     >                    'zsw1',zsw1(2),zsw1(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,nn*nprj_max,
     >                    'zsw2',zsw2(2),zsw2(1))
      if (.not.value) call errquit('cpsp_v_nonlocal:pushing stack',0,
     &       MA_ERR)

      value = BA_push_get(mt_dbl, nfft3d,'xtmp',xtmp(2),xtmp(1))
      value = value.and.
     >        BA_push_get(mt_dbl, nfft3d,'Gx',Gx(2),Gx(1))
      value = value.and.
     >        BA_push_get(mt_dbl, nfft3d,'Gy',Gy(2),Gy(1))
      value = value.and.
     >        BA_push_get(mt_dbl, nfft3d,'Gz',Gz(2),Gz(1))

      value = value.and.
     >      BA_push_get(mt_dbl,3*nn,'sum',sum(2),sum(1))
      value = value.and.
     >      BA_push_get(mt_dbl,3*nion,'ftmp',ftmp(2),ftmp(1))
      if (.not.value) call errquit('cpsp_v_nonlocal:pushing stack',1,
     &       MA_ERR)
      call dcopy(3*nion,0.0d0,0,dbl_mb(ftmp(1)),1)
      G(1)  = c_G_indx(1)
      G(2)  = c_G_indx(2)
      G(3)  = c_G_indx(3)

      omega = lattice_omega()
      scal  = 1.0d0/omega

      do nbq = 1,nbrillq
        call Cram_npack(nbq,npack)
        weight = brillioun_weight(nbq)
        psi1_shift = cpsi_data_get_chnk(psi1_tag,nbq)
        if (occ1_tag.gt.0) 
     >     occ1_shift = cpsi_data_get_chnk(occ1_tag,nbq)


*       **** define Gx,Gy and Gz in packed space ****
        call C3dB_t_Copy(1,dbl_mb(G(1)),dbl_mb(Gx(1)))
        call C3dB_t_Copy(1,dbl_mb(G(2)),dbl_mb(Gy(1)))
        call C3dB_t_Copy(1,dbl_mb(G(3)),dbl_mb(Gz(1)))
        call Cram_r_pack(nbq,dbl_mb(Gx(1)))
        call Cram_r_pack(nbq,dbl_mb(Gy(1)))
        call Cram_r_pack(nbq,dbl_mb(Gz(1)))

        do 100 ii=1,nion
           ia=ion_katm(ii)
           it=int_mb(psp_type(1)+ia-1)
           nproj = int_mb(nprj(1)+ia-1)
           if (nproj.gt.0) then
              if (it.eq.7) then 
                call cpsp_f_nonlocal_rel(nbq,ii,ia,nproj,npack1,ne(1),
     >          dcpl_mb(exi(1)),dcpl_mb(zsw1(1)),
     >          dbl_mb(psi1_shift),dcpl_mb(prjtmp(1)),
     >          dbl_mb(xtmp(1)),dbl_mb(sum(1)),
     >          dbl_mb(Kijl(1)+(ia-1)*kij_stride),
     >          dbl_mb(Gx(1)),dbl_mb(Gy(1)),dbl_mb(Gz(1)),
     >          dbl_mb(ftmp(1)+3*(ii-1)),
     >          dbl_mb(ftmp(1)+3*(ii-1)+1),
     >          dbl_mb(ftmp(1)+3*(ii-1)+2),
     >          weight,scal)
                goto 100 
              end if

*             **** structure factor ****
              call cstrfac_pack(nbq,ii,dcpl_mb(exi(1)))
              call cstrfac_k(ii,nbq,cxr)
c              call Cram_c_ZMul(nbq,cxr,dcpl_mb(exi(1)),dcpl_mb(exi(1)))
              call zscal(npack,cxr,dcpl_mb(exi(1)),1)

              do l=1,nproj
c                 shift = vnl(1)
c     >                      +(l-1) *npack1
c     >                      +(nb-1)*npack1*vnl_stride
c     >                      +(ia-1)*npack1*vnl_stride*nbrill
                 shift = cpsp_projector_get_ptr(
     >                     int_mb(vnl(1)+ia-1),nbq,l)
                 l_prj = int_mb(l_projector(1)+(l-1)
     >                          +(ia-1)*jmmax_max)
              
                sd_function=.true.
                if (mod(l_prj,2).ne.0) then
                  sd_function=.false.
                end if

*                **** phase factor does not matter therefore ****
*                **** (-i)^l is the same as (i)^l in the     ****
*                **** Rayleigh scattering formula            ****
*                *** current function is s or d ****
                 if (sd_function) then
                    call Cram_rc_Mul(nbq,dbl_mb(shift),
     >                                  dcpl_mb(exi(1)),
     >                                  dcpl_mb(prjtmp(1)+(l-1)*npack1))
*                *** current function is p or f ****
                 else
                    call Cram_irc_Mul(nbq,dbl_mb(shift),
     >                                   dcpl_mb(exi(1)),
     >                                  dcpl_mb(prjtmp(1)+(l-1)*npack1))
                 end if

*                **** compute nnXnproj matrix zsw1 = <psi1|prj> ****
                 call Cram_cc_inzdot(nbq,nn,
     >                         dbl_mb(psi1_shift),
     >                         dcpl_mb(prjtmp(1)+(l-1)*npack1),
     >                         dcpl_mb(zsw1(1)+(l-1)*nn))

*                ***** scale psp by factor - used for generating antiferromagnetic structures ****
*                **** nwchem input: pspspin up/down scale l ion_numbers                       ****
                 if (pspspin) then
                    if (log_mb(pspspin_upions(1)+ii-1).and.
     >                 (l_prj.eq.pspspin_upl))
     >                 call dscal(2*ne(1),pspspin_upscale,
     >                            dcpl_mb(zsw1(1)+(l-1)*nn),1)
                    if (log_mb(pspspin_downions(1)+ii-1).and.
     >                 (l_prj.eq.pspspin_downl))
     >                 call dscal(2*ne(2),pspspin_downscale,
     >                            dcpl_mb(zsw1(1)+(l-1)*nn+ne(1)),1)
                 end if
              end do
              call C3dB_Vector_SumAll((2*nn*nproj),dcpl_mb(zsw1(1)))


*             **** do kleinman-bylander multiplication ****
*             **** sw2 = Gijl*sw1 ******
              call Multiply_Gijl_zsw1(nn,
     >                          nproj,
     >                          int_mb(nmax(1)+ia-1),
     >                          int_mb(lmax(1)+ia-1),
     >                          int_mb(n_projector(1)
     >                                 + (ia-1)*jmmax_max),
     >                          int_mb(l_projector(1)
     >                                 + (ia-1)*jmmax_max),
     >                          int_mb(m_projector(1)
     >                                 + (ia-1)*jmmax_max),
     >                          dbl_mb(Gijl(1)
     >                         + (ia-1)*gij_stride),
     >                          dcpl_mb(zsw1(1)),
     >                          dcpl_mb(zsw2(1)))

              call dscal(2*nn*nproj,scal,dcpl_mb(zsw2(1)),1)
              if (ispin.eq.1) 
     >           call dscal(2*nn*nproj,2.0d0,dcpl_mb(zsw2(1)),1)

              do l=1,nproj

                 psi_shift = psi1_shift
                 do n=1,nn
                    call Cram_zccr_Multiply2(nbq,
     >                                 dcpl_mb(zsw2(1)+(l-1)*nn+n-1),
     >                                 dbl_mb(psi_shift),
     >                                 dcpl_mb(prjtmp(1)+(l-1)*npack1),
     >                                 dbl_mb(xtmp(1)))
                    psi_shift = psi_shift + nshift
c                    do i=1,npack
c                       ctmp = psi1(i+(n-1)*npack1
c     >                              +(nb-1)*neall*npack1)
c     >                      *dconjg(dcpl_mb(prjtmp(1)+(l-1)*npack1+i-1))
c     >                      *dconjg(dcpl_mb(zsw2(1)+(l-1)*nn+n-1))
c                       dbl_mb(xtmp(1)+i-1) = dimag(ctmp)
c                    end do
 
                    call Cram_rr_idot(nbq,dbl_mb(Gx(1)),dbl_mb(xtmp(1)),
     >                                   dbl_mb(sum(1)+3*(n-1)))
                    call Cram_rr_idot(nbq,dbl_mb(Gy(1)),dbl_mb(xtmp(1)),
     >                                dbl_mb(sum(1)+1+3*(n-1)))
                    call Cram_rr_idot(nbq,dbl_mb(Gz(1)),dbl_mb(xtmp(1)),
     >                                dbl_mb(sum(1)+2+3*(n-1)))
                 end do
                 call C3dB_Vector_Sumall(3*(nn),dbl_mb(sum(1)))
  
                 wfrac = 1.0d0
                 if (occ1_tag.gt.0) occ_shift = occ1_shift
                 do n=1,nn
                    if (occ1_tag.gt.0) then
                       wfrac = dbl_mb(occ_shift)
                       occ_shift = occ_shift+1
                    end if
                    dbl_mb(ftmp(1)+3*(ii-1)) = 
     >                         dbl_mb(ftmp(1)+3*(ii-1))
     >                         + 2.0d0*wfrac*weight
     >                           *dbl_mb(sum(1)+3*(n-1))
                    dbl_mb(ftmp(1)+3*(ii-1)+1) = 
     >                         dbl_mb(ftmp(1)+3*(ii-1)+1)
     >                         + 2.0d0*wfrac*weight
     >                           *dbl_mb(sum(1)+1+3*(n-1))
                    dbl_mb(ftmp(1)+3*(ii-1)+2) = 
     >                         dbl_mb(ftmp(1)+3*(ii-1)+2)
     >                         + 2.0d0*wfrac*weight
     >                          *dbl_mb(sum(1)+2+3*(n-1))
                 end do
        
              end do !** l **

           end if !** nproj>0) ***
100      continue !** ii **

      end do !** nbq ***

      call K1dB_Vector_SumAll(3*nion,dbl_mb(ftmp(1)))
      call daxpy(3*nion,1.0d0,dbl_mb(ftmp(1)),1,fion,1)
      value =           BA_pop_stack(ftmp(2))
      value = value.and.BA_pop_stack(sum(2))
      value = value.and.BA_pop_stack(Gz(2))
      value = value.and.BA_pop_stack(Gy(2))
      value = value.and.BA_pop_stack(Gx(2))
      value = value.and.BA_pop_stack(xtmp(2))
      if (.not.value) call errquit('cpsp_f_vnonlocal:popping stack',2,
     &       MA_ERR)

      value =           BA_pop_stack(zsw2(2))
      value = value.and.BA_pop_stack(zsw1(2))
      value = value.and.BA_pop_stack(exi(2))
      if (.not.value) call errquit('cpsp_f_vnonlocal:popping stack',3,
     &       MA_ERR)

      call nwpw_timing_end(6)
      

      return 
      end

*     *************************************************************
*     *                                                           *
*     *                   cpsp_read                               *
*     *                                                           *
*     *************************************************************
cccccccccccccccccccccccccccccccccccccccccccccccccc
c Read data from psp file *.cpp
cccccccccccccccccccccccccccccccccccccccccccccccccc 
      subroutine cpsp_read(fname,comment,
     >                       psp_type,
     >                       version,
     >                       nfft,unita,
     >                       atom,amass,zv,lmmax,lmax,locp,nmax,
     >                       rc,
     >                       nprj,n_projector,l_projector,m_projector,
     >                       Gijl,
     >                       Kijl,
     >                       nfft3d,npack0,npack1,lmmax_max,nmax_max,
     >                       vl,vnl_tag,vnl2_tag,
     >                       semicore,rcore,ncore,
     >                       tmp,tmp2,
     >                       ierr)
      implicit none 
      character*(*) comment
      character*50 fname
      integer psp_type
      integer version
      integer nfft(3)
      real*8  unita(3,3)
      character*2 atom
      real*8 amass,zv
      integer lmmax
      integer lmax
      integer locp
      integer nmax
      real*8 rc(*)
      integer nprj,n_projector(*),l_projector(*),m_projector(*)
      real*8 Gijl(*),Kijl(*)
      integer nfft3d,npack0,npack1,lmmax_max,nmax_max
      real*8 vl(*)
      integer vnl_tag
      integer vnl2_tag
      logical semicore
      real*8  rcore
      real*8  ncore(*)
      real*8     tmp(*)
      real*8     tmp2(*)
      integer ierr

#include "bafdecls.fh"
#include "btdb.fh"
#include "util.fh"

#ifdef TCGMSG
#include "tcgmsg.fh"
#include "msgtypesf.h"
#endif

*    *** local variables ***
      integer MASTER,taskid,taskid_k
      parameter(MASTER=0)
      integer n,l,nb,nbq,pk
      integer msglen,rsize
      integer iatom(2)
      character*255 full_filename
      real*8 kv(3)
      integer nbrillioun
      logical mprint,value

*     ***** external functions ****
      integer  brillioun_nbrillioun,cpsp_projector_alloc
      integer  brillioun_nbrillq
      real*8   brillioun_all_k
      logical  control_print
      external brillioun_nbrillioun,cpsp_projector_alloc
      external brillioun_nbrillq
      external brillioun_all_k
      external control_print

      call Parallel_taskid(taskid)
      call Parallel3d_taskid_k(taskid_k)
      mprint = (taskid.eq.MASTER).and.control_print(print_medium)

*     **** open fname binary file ****
      
      if (taskid.eq.MASTER) then
         call util_file_name_noprefix(fname,.false.,
     >                             .false.,
     >                       full_filename)
         l = index(full_filename,' ') - 1
         call openfile(5,full_filename,l,'r',l)
         call cread(5,comment,80)
         call iread(5,psp_type,1)
         call iread(5,version,1)
         call iread(5,nfft,3)
         call dread(5,unita,9)
         call cread(5,atom,2)
         call dread(5,amass,1)
         call dread(5,zv,1)
         call iread(5,lmax,1)
         call iread(5,locp,1)
         call iread(5,nmax,1)
         lmmax=(lmax+1)**2 - (2*locp+1)
         amass = amass*1822.89d0
         call dread(5,rc,(lmax+1))
         call iread(5,nprj,1)
         if (nprj.gt.0) then
           call iread(5,n_projector,nprj)
           call iread(5,l_projector,nprj)
           call iread(5,m_projector,nprj)
           if (psp_type.eq.7) then
              call dread(5,Kijl,nprj)
           else
ccc!** number of matrix elements = nmax*nmax*(lmax+1) **
             rsize=nmax*nmax*(lmax+1)
             call dread(5,Gijl,rsize)
             if (psp_type.eq.1) then
               call dread(5,Kijl,rsize)
             end if
ccc!** number of matrix elements = nmax*nmax*(lmax+1) **
           end if
         end if  
         call dread(5,rcore,1)
         call iread(5,nbrillioun,1)
         ierr = 0
         if (nbrillioun.eq.brillioun_nbrillioun()) then
           do nb=1,nbrillioun
               call dread(5,kv,3)
               if ((brillioun_all_k(1,nb).ne.kv(1)).or.  
     >             (brillioun_all_k(2,nb).ne.kv(2)).or.
     >             (brillioun_all_k(3,nb).ne.kv(3))) 
     >          ierr = 1
           end do
           if (ierr.eq.1) then
              if (mprint) then
                 write(*,*)"Brillioun Zone Vectors do not match!"
                 call flush(6)
              end if
           end if
         else
           if (mprint) then
           write(*,*)"Brillioun Zone Points do not match!"
           write(*,*)"NB = ",nbrillioun," not equal ",
     >       brillioun_nbrillioun()
           call flush(6)
           end if
           ierr = 1
         end if
      end if

      msglen = 1
      call Parallel_Brdcst_ivalues(MASTER,msglen,ierr)
      if (ierr.ne.0) then
         if (taskid.eq.MASTER) call closefile(5)
         return
      end if

*     **** send header data to all processors ****
      msglen = 1
      call Parallel_Brdcst_ivalues(MASTER,msglen,version)
      call Parallel_Brdcst_ivalues(MASTER,msglen,psp_type)
      msglen = 3
      call Parallel_Brdcst_ivalues(MASTER,msglen,nfft)
      msglen = 9
      call Parallel_Brdcst_values(MASTER,msglen,unita)

      iatom(1) = ichar(atom(1:1))
      iatom(2) = ichar(atom(2:2))
      msglen = 2
      call BRDCST(9+MSGCHR,iatom,mitob(msglen),MASTER)
      atom(1:1) = char(iatom(1))
      atom(2:2) = char(iatom(2))

      msglen = 1
      call Parallel_Brdcst_values(MASTER,msglen,amass)
      msglen = 1
      call Parallel_Brdcst_values(MASTER,msglen,zv)
      msglen = 1
      call Parallel_Brdcst_ivalues(MASTER,msglen,lmax)
      call Parallel_Brdcst_ivalues(MASTER,msglen,locp)
      call Parallel_Brdcst_ivalues(MASTER,msglen,nmax)
      call Parallel_Brdcst_ivalues(MASTER,msglen,nprj)
      lmmax=(lmax+1)**2 - (2*locp+1)

      msglen=lmax+1
      call Parallel_Brdcst_values(MASTER,msglen,rc)

      msglen=nprj
      call Parallel_Brdcst_ivalues(MASTER,msglen,n_projector)
      call Parallel_Brdcst_ivalues(MASTER,msglen,l_projector)
      call Parallel_Brdcst_ivalues(MASTER,msglen,m_projector)

      if (psp_type.eq.7) then
       msglen=nprj
       call Parallel_Brdcst_values(MASTER,msglen,Kijl)
      else
       msglen=nmax*nmax*(lmax+1)
       call Parallel_Brdcst_values(MASTER,msglen,Gijl)
       if(psp_type.eq.1) call Parallel_Brdcst_values(MASTER,msglen,Kijl)
      end if
      msglen=1
      call Parallel_Brdcst_values(MASTER,msglen,rcore)


*     **** determine semicore value ****
      if (rcore.gt.0.0d0) then
         semicore = .true.
      else
         semicore = .false.
      end if

*     *** read in vl 3d block ***
      call C3dB_r_read(1,5,tmp2,tmp,-1,-1,.true.)
      call Cram_r_pack(0,tmp2)
      call Cram_r_Copy(0,tmp2,vl)


*     **** read in vnl 3d blocks ****

ccccccccccccccccccc
      if (psp_type.eq.7) then
cccccccccccccccccccc
c Relativistic Non-Local PPot.
c spin-orbit included...
cccccccccccccccccccc
        vnl_tag  = cpsp_projector_alloc(brillioun_nbrillq(),
     >                                  nprj,2*npack1)
        vnl2_tag = cpsp_projector_alloc(brillioun_nbrillq(),
     >                                  nprj,2*npack1)
!*** XXXX double check this INPUT ****
        do nb=1,brillioun_nbrillioun() 
          call K1dB_ktoqp(nb,nbq,pk)
          do n=1,nprj
            call C3dB_c_Read(1,5,tmp2,tmp,-1,pk)
            if (pk.eq.taskid_k) then
               call Cram_c_pack(nbq,tmp2)
               call cpsp_projector_add(vnl_tag,nbq,n,tmp2)
            end if
          end do
          do n=1,nprj
            call C3dB_c_read(1,5,tmp2,tmp,-1,pk)
            if (pk.eq.taskid_k) then
               call Cram_c_pack(nbq,tmp2)
               call cpsp_projector_add(vnl2_tag,nbq,n,tmp2)
            end if
          end do
        end do

      else
        vnl_tag = cpsp_projector_alloc(brillioun_nbrillq(),
     >                                 nprj,npack1)
        do nb=1,brillioun_nbrillioun() 
          call K1dB_ktoqp(nb,nbq,pk)
          do n=1,nprj

            call C3dB_r_read(1,5,tmp2,tmp,-1,pk,.true.)

            if (pk.eq.taskid_k) then
               call Cram_r_pack(nbq,tmp2)
               call cpsp_projector_add(vnl_tag,nbq,n,tmp2)
            endif
          end do
        end do
*     **** read in v_spin_orbit 3d blocks ****
        if (psp_type.eq.1) then
          vnl2_tag = cpsp_projector_alloc(brillioun_nbrillq(),
     >                                    nprj,2*npack1)
          do nb=1,brillioun_nbrillioun()
            call K1dB_ktoqp(nb,nbq,pk)
            do n=1,nprj
              call C3dB_c_read(1,5,tmp2,tmp,-1,pk)
              if (pk.eq.taskid_k) then
                 call Cram_c_pack(nbq,tmp2)
                 call cpsp_projector_add(vnl2_tag,nbq,n,tmp2)
              endif
            end do
          end do
        end if
      end if
*     **** read in semicore density block ****
      if (semicore) then
         call C3dB_r_read(1,5,tmp2,tmp,-1,-1,.true.)
         call Cram_r_pack(0,tmp2)
         call Cram_r_Copy(0,tmp2,ncore(1))

         call C3dB_r_read(1,5,tmp2,tmp,-1,-1,.true.)
         call Cram_r_pack(0,tmp2)
         call Cram_r_Copy(0,tmp2,ncore(1+2*npack0))

         call C3dB_r_read(1,5,tmp2,tmp,-1,-1,.true.)
         call Cram_r_pack(0,tmp2)
         call Cram_r_Copy(0,tmp2,ncore(1+3*npack0))

         call C3dB_r_read(1,5,tmp2,tmp,-1,-1,.true.)
         call Cram_r_pack(0,tmp2)
         call Cram_r_Copy(0,tmp2,ncore(1+4*npack0))
      end if

*     *** close fname binary file ***
      if (taskid.eq.MASTER) then
         call closefile(5)
      end if

      ierr = 0
      return
      end


*     ***********************************
*     *					*
*     *	 	  cpsp_readall  	*
*     *					*
*     ***********************************

      subroutine cpsp_readall()
      implicit none
      
#include "bafdecls.fh"
#include "cpsp_common.fh"
#include "c_semicore_common.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "util.fh"

*     **** local variables ****
      integer ngp(3),version,nfft3d,npack0,npack1
      integer ia,l
      real*8 unita(3,3)
      character*12 boundry
      integer tmp(2),tmp2(2),ierr
      logical value,found,correct_box,mprint
      character*5  element
      character*50 fname

*     **** parallel i/o variable ****
      integer MASTER,taskid
      parameter(MASTER=0)

*     **** external functions ****
      logical      nwpw_filefind,control_spin_orbit,control_print
      integer      control_ngrid
      real*8       control_unita
      character*12 control_boundry
      character*4  ion_atom
      external     nwpw_filefind,control_spin_orbit,control_print
      external     control_ngrid
      external     control_unita
      external     control_boundry
      external     ion_atom
      

      call C3dB_nfft3d(1,nfft3d)
      call Cram_npack(0,npack0)
      call Cram_max_npack(npack1)
      call Parallel_taskid(taskid)
      !nbrill = brillioun_nbrillioun()

      mprint = (taskid.eq.MASTER).and.control_print(print_medium)

*     *** set semicore(0) *****
      log_mb(semicore(1)) = .false.

      value = BA_push_get(mt_dbl,(4*nfft3d),'tmp',tmp(2),tmp(1))
      if (.not. value) 
     > call errquit('cpsp_readall:out of stack memory',0,MA_ERR)

      value = BA_push_get(mt_dbl,(2*nfft3d),'tmp2',tmp2(2),tmp2(1))
      if (.not. value) 
     > call errquit('cpsp_readall:out of stack memory',0,MA_ERR)

      do_spin_orbit=.false.

*     **** read pseudopotentials ****
      do ia=1,npsp

*      **** define formatted psp name ****
       element = '     '
       element = ion_atom(ia)
       l = index(element,' ') - 1
       fname = element(1:l)//'.cpp'

       found = .false.
       do while (.not.found)

         if (nwpw_filefind(fname)) then
            call cpsp_read(fname,comment(ia),
     >                  int_mb(psp_type(1)+ia-1),
     >                  version,
     >                  ngp,unita,
     >                  atom(ia),
     >                  dbl_mb(amass(1)+ia-1),
     >                  dbl_mb(zv(1)+ia-1),
     >                  int_mb(lmmax(1)+ia-1),
     >                  int_mb(lmax(1)+ia-1),
     >                  int_mb(locp(1)+ia-1),
     >                  int_mb(nmax(1)+ia-1),
     >                  dbl_mb(rc(1) + (ia-1)*(lmax_max+1)),
     >                  int_mb(nprj(1)+ia-1),
     >                  int_mb(n_projector(1)
     >                     + (ia-1)*jmmax_max),
     >                  int_mb(l_projector(1)
     >                     + (ia-1)*jmmax_max),
     >                  int_mb(m_projector(1)
     >                     + (ia-1)*jmmax_max),
     >                  dbl_mb(Gijl(1)
     >                     + (ia-1)*gij_stride),
     >                  dbl_mb(Kijl(1)
     >                     + (ia-1)*kij_stride),
     >                  nfft3d,npack0,npack1,lmmax_max,nmax_max,
     >                  dbl_mb(vl(1) + (ia-1)*npack0),
     >                  int_mb(vnl(1)   + ia-1),
     >                  int_mb(vnlso(1) + ia-1),
     >                  log_mb(semicore(1)+ia),
     >                  dbl_mb(rcore(1)+ia-1),
     >                  dbl_mb(ncore(1)+ (ia-1)*npack0*5),
     >                  dbl_mb(tmp(1)),dbl_mb(tmp2(1)),
     >                  ierr)

           if ((int_mb(psp_type(1)+(ia-1)).eq.1).and.
     >           control_spin_orbit()) then
              do_spin_orbit=.true.
              if (mprint) then
                write(*,*)"HGH spin orbit pseudopotential activated"
                call flush(6)
              end if
           end if
           if ((int_mb(psp_type(1)+(ia-1)).eq.7)) then
             if (mprint) then
               write(*,*)"RTC spin orbit pseudopotential activated"
               call flush(6)
             end if
           end if

*          **** set semicore(0) ****
           if (log_mb(semicore(1)+ia)) log_mb(semicore(1)) = .true.
           if (ierr.gt.1) go to 9000
*          **************************************************************
*          ***** logic for finding out if psp is correctly formatted ****
*          **************************************************************
           correct_box = .true.
           boundry = control_boundry()
           l =index(boundry,' ') - 1
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
              if (mprint) then
              write(luout,*)
              write(luout,*) 
     >        "pseudopotential is not correctly formatted:",fname
              end if
              if ((ierr.eq.0).and.(int_mb(nprj(1)+ia-1).gt.0)) then
                 call cpsp_projector_dealloc(int_mb(vnl(1)+ia-1))
                 if ((int_mb(psp_type(1)+ia-1).eq.7)
     >               .or.(do_spin_orbit)) then
                    call cpsp_projector_dealloc(int_mb(vnlso(1)+ia-1))
                 end if
              end if
           end if
           if (correct_box) found = .true.
           if (ierr.eq.1)   then
              found = .false.
              if (mprint) then
              write(luout,*)
              write(luout,*) 
     >        "pseudopotential is not correctly formatted---",
     >         "bad brillioun zone:",fname
              end if
           end if

         end if

*        **** generate formatted pseudopotential atom.cpp *****
         if (.not.found) then
             call cpsp_formatter_auto(ion_atom(ia))
         end if

       end do !***do while ****


      end do

 9000 value =           BA_pop_stack(tmp2(2))
      value = value.and.BA_pop_stack(tmp(2))
      if (.not. value) 
     > call errquit('cpsp_readall:error popping stack',0,MA_ERR)

*     **** done reading set nprj_max and prjtmp for nonlocal psp operator ****
      call cpsp_proj_init()

      return
      end




