#define TCGMSG

      subroutine cpsp_init()      
      implicit none

#include "mafdecls.fh"
#include "cpsp_common.fh"
#include "errquit.fh"

      !integer nfft3d
      integer npack0,npack1,nbrillioun
      logical value

*     **** external functions *****
      integer  ion_nkatm, brillioun_nbrillioun
      external ion_nkatm, brillioun_nbrillioun
 


      npsp = ion_nkatm()
      nbrillioun = brillioun_nbrillioun()

      call Cram_npack(0,npack0)
      call Cram_max_npack(npack1)
      !call C3dB_nfft3d(1,nfft3d)

      value = MA_alloc_get(mt_dcpl,(nmax_max*lmmax_max*npack1),
     >                    'prjtmp',prjtmp(2),prjtmp(1))

      value = value.and.
     >        MA_alloc_get(mt_dbl,(npsp*npack0),
     >                    'vl',vl(2),vl(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,
     >                    (nmax_max*lmmax_max*npsp*npack1*nbrillioun),
     >                    'vnl',vnl(2),vnl(1))
      value = value.and.
     >        MA_alloc_get(mt_dcpl,
     >                    (nmax_max*lmmax_max*npsp*npack1*nbrillioun),
     >                    'vnlso',vnlso(2),vnlso(1))

      value = value.and.
     >      MA_alloc_get(mt_dbl,(nmax_max*nmax_max*(lmax_max+1)*npsp),
     >                    'Gijl',Gijl(2),Gijl(1))

      value = value.and.
     >      MA_alloc_get(mt_dbl,(nmax_max*nmax_max*(lmax_max+1)*npsp),
     >                    'Kijl',Kijl(2),Kijl(1))

      value = value.and.
     >        MA_alloc_get(mt_int,(npsp),'nprj',nprj(2),nprj(1))
      value = value.and.
     >        MA_alloc_get(mt_int,(npsp*nmax_max*lmmax_max),
     >                     'n_projector',n_projector(2),n_projector(1))
      value = value.and.
     >        MA_alloc_get(mt_int,(npsp*nmax_max*lmmax_max),
     >                     'l_projector',l_projector(2),l_projector(1))
      value = value.and.
     >        MA_alloc_get(mt_int,(npsp*nmax_max*lmmax_max),
     >                     'm_projector',m_projector(2),m_projector(1))


      value = value.and.
     >        MA_alloc_get(mt_dbl,(npsp),'zv',zv(2),zv(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,(npsp),'amass',amass(2),amass(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,(npsp*(lmax_max+1)),'rc',rc(2),rc(1))
      value = value.and.
     >        MA_alloc_get(mt_int,(npsp),'lmmax',lmmax(2),lmmax(1))
      value = value.and.
     >        MA_alloc_get(mt_int,(npsp),'lmax',lmax(2),lmax(1))
      value = value.and.
     >        MA_alloc_get(mt_int,(npsp),'locp',locp(2),locp(1))
      value = value.and.
     >        MA_alloc_get(mt_int,(npsp),'nmax',nmax(2),nmax(1))
      value = value.and.
     >        MA_alloc_get(mt_int,(npsp),
     >                     'psp_type',psp_type(2),psp_type(1))


      if (.not. value) 
     >   call errquit('cpsp_init:out of heap memory',0, MA_ERR)

      call dcopy(npsp*npack0,          0.0d0,0,dbl_mb(vl(1)), 1)
      call dcopy(lmmax_max*npsp*npack1*nbrillioun,
     >           0.0d0,0,dbl_mb(vnl(1)),1)
      call dcopy(npsp,                 0.0d0,0,dbl_mb(zv(1)),1)
      call dcopy(npsp,                 0.0d0,0,dbl_mb(amass(1)),1)
      call dcopy(npsp*(lmax_max+1),    0.0d0,0,dbl_mb(rc(1)),1)

*     **** allocate semicore data ****
      call c_semicore_init()

      return
      end

      subroutine cpsp_end()      
      implicit none

#include "mafdecls.fh"
#include "errquit.fh"
#include "cpsp_common.fh"


      logical value

*     **** external functions ****

*     **** deallocate semicore data ****
      call c_semicore_end()

      value = MA_free_heap(prjtmp(2))
      value = value.and.MA_free_heap(vl(2))
      value = value.and.MA_free_heap(vnlso(2))
      value = value.and.MA_free_heap(vnl(2))
      value = value.and.MA_free_heap(Kijl(2))
      value = value.and.MA_free_heap(Gijl(2))
      value = value.and.MA_free_heap(nprj(2))
      value = value.and.MA_free_heap(n_projector(2))
      value = value.and.MA_free_heap(l_projector(2))
      value = value.and.MA_free_heap(m_projector(2))
      value = value.and.MA_free_heap(zv(2))
      value = value.and.MA_free_heap(amass(2))
      value = value.and.MA_free_heap(rc(2))
      value = value.and.MA_free_heap(lmmax(2))
      value = value.and.MA_free_heap(lmax(2))
      value = value.and.MA_free_heap(locp(2))
      value = value.and.MA_free_heap(nmax(2))
      value = value.and.MA_free_heap(psp_type(2))
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

#include "mafdecls.fh"
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

#include "mafdecls.fh"
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

#include "mafdecls.fh"
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

#include "mafdecls.fh"
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

#include "mafdecls.fh"
#include "cpsp_common.fh"

      cpsp_nprj  = int_mb(nprj(1)+ia-1)
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

#include "mafdecls.fh"
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

#include "mafdecls.fh"
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

#include "mafdecls.fh"
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

#include "mafdecls.fh"
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

      value = MA_push_get(mt_dcpl,npack0,'exi', exi(2), exi(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,nfft3d,'vtmp',vtmp(2),vtmp(1))
      if (.not. value) call errquit('cpsp_v_local: pushing stack',0,
     &       MA_ERR)

      if (move) then
        value = MA_push_get(mt_dbl, nfft3d,'xtmp',xtmp(2),xtmp(1))
        value = value.and.
     >          MA_push_get(mt_dbl, nfft3d,'Gx',Gx(2),Gx(1))
        value = value.and.
     >          MA_push_get(mt_dbl, nfft3d,'Gy',Gy(2),Gy(1))
        value = value.and.
     >          MA_push_get(mt_dbl, nfft3d,'Gz',Gz(2),Gz(1))
        if (.not. value) call errquit('cpsp_v_local: pushing stack',1,
     &       MA_ERR)
        G(1)  = c_G_indx(1)
        G(2)  = c_G_indx(2)
        G(3)  = c_G_indx(3)

*       **** define Gx,Gy and Gz in packed space ****
        call C3dB_r_Copy(1,dbl_mb(G(1)),dbl_mb(Gx(1)))
        call C3dB_r_Copy(1,dbl_mb(G(2)),dbl_mb(Gy(1)))
        call C3dB_r_Copy(1,dbl_mb(G(3)),dbl_mb(Gz(1)))
        call Cram_r_pack(0,dbl_mb(Gx(1)))
        call Cram_r_pack(0,dbl_mb(Gy(1)))
        call Cram_r_pack(0,dbl_mb(Gz(1)))
      end if

      call dcopy((2*npack0),0.0d0,0,vl_out,1)
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
        call Cram_cc_Sum(0,vl_out,dcpl_mb(vtmp(1)),vl_out)


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
        value = MA_pop_stack(Gz(2))
        value = value.and.MA_pop_stack(Gy(2))
        value = value.and.MA_pop_stack(Gx(2))
        value = value.and.MA_pop_stack(xtmp(2))
        if (.not. value) call errquit('cpsp_v_local: popping stack',0,
     &       MA_ERR)
      end if
      value = MA_pop_stack(vtmp(2))
      value = value.and.MA_pop_stack(exi(2))
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

#include "mafdecls.fh"
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
      value = MA_push_get(mt_dcpl,npack0,'exi', exi(2), exi(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,nfft3d,'vtmp',vtmp(2),vtmp(1))
      value = value.and.
     >        MA_push_get(mt_dbl, nfft3d,'xtmp',xtmp(2),xtmp(1))
      value = value.and.
     >        MA_push_get(mt_dbl, nfft3d,'Gx',Gx(2),Gx(1))
      value = value.and.
     >        MA_push_get(mt_dbl, nfft3d,'Gy',Gy(2),Gy(1))
      value = value.and.
     >        MA_push_get(mt_dbl, nfft3d,'Gz',Gz(2),Gz(1))
      if (.not. value) call errquit('cpsp_f_vlocal:pushing stack',0,
     &       MA_ERR)
      G(1)  = c_G_indx(1)
      G(2)  = c_G_indx(2)
      G(3)  = c_G_indx(3)

*     **** define Gx,Gy and Gz in packed space ****
      call C3dB_r_Copy(1,dbl_mb(G(1)),dbl_mb(Gx(1)))
      call C3dB_r_Copy(1,dbl_mb(G(2)),dbl_mb(Gy(1)))
      call C3dB_r_Copy(1,dbl_mb(G(3)),dbl_mb(Gz(1)))
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
      value = MA_pop_stack(Gz(2))
      value = value.and.MA_pop_stack(Gy(2))
      value = value.and.MA_pop_stack(Gx(2))
      value = value.and.MA_pop_stack(xtmp(2))
      value = value.and.MA_pop_stack(vtmp(2))
      value = value.and.MA_pop_stack(exi(2))
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

      subroutine cpsp_v_nonlocal(ispin,ne,
     >                           psi1,psi2,move,fion)
      implicit none
      integer    ispin,ne(2)
      complex*16 psi1(*)
      complex*16 psi2(*)
      logical move
      real*8 fion(3,*)

#include "mafdecls.fh"
#include "cpsp_common.fh"
#include "errquit.fh"

*     *** local variables ***
      complex*16 one,mone
      parameter  (one=(1.0d0,0.0d0), mone=(-1.0d0,0.0d0))

      integer nfft3d,G(3),npack1,npack
      integer i,ii,ia,l,n,nn,nb,neall,nbrill
      integer k,shift,l_prj,nproj
      real*8  omega,weight,scal
      complex*16 ctmp,cxr
      integer exi(2),xtmp(2),zsw1(2),zsw2(2),sum(2)
      integer Gx(2),Gy(2),Gz(2)
      logical value,sd_function

*     **** external functions ****
      logical  is_sORd
      integer  ion_nion,ion_katm,c_G_indx,brillioun_nbrillioun
      real*8   lattice_omega,brillioun_weight
      external is_sORd
      external ion_nion,ion_katm,c_G_indx,brillioun_nbrillioun
      external lattice_omega,brillioun_weight

      call nwpw_timing_start(6)
      
*     **** allocate local memory ****
      nn = ne(1)+ne(2)
      neall = ne(1)+ne(2)
      nbrill = brillioun_nbrillioun()
      call C3dB_nfft3d(1,nfft3d)
      call Cram_max_npack(npack1)

      value = MA_push_get(mt_dcpl,npack1,'exi', exi(2), exi(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,nn*nmax_max*lmmax_max,
     >                    'zsw1',zsw1(2),zsw1(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,nn*nmax_max*lmmax_max,
     >                    'zsw2',zsw2(2),zsw2(1))
      if (.not.value) call errquit('cpsp_v_nonlocal:pushing stack',0,
     &       MA_ERR)

      if (move) then
        value = MA_push_get(mt_dbl, nfft3d,'xtmp',xtmp(2),xtmp(1))
        value = value.and.
     >          MA_push_get(mt_dbl, nfft3d,'Gx',Gx(2),Gx(1))
        value = value.and.
     >          MA_push_get(mt_dbl, nfft3d,'Gy',Gy(2),Gy(1))
        value = value.and.
     >          MA_push_get(mt_dbl, nfft3d,'Gz',Gz(2),Gz(1))
        value = value.and.
     >        MA_push_get(mt_dbl,3*nn,'sum',sum(2),sum(1))
        if (.not.value) call errquit('cpsp_v_nonlocal:pushing stack',1,
     &       MA_ERR)
        G(1)  = c_G_indx(1)
        G(2)  = c_G_indx(2)
        G(3)  = c_G_indx(3)

*       **** define Gx,Gy and Gz in packed space ****
        call C3dB_r_Copy(1,dbl_mb(G(1)),dbl_mb(Gx(1)))
        call C3dB_r_Copy(1,dbl_mb(G(2)),dbl_mb(Gy(1)))
        call C3dB_r_Copy(1,dbl_mb(G(3)),dbl_mb(Gz(1)))
      end if

      omega = lattice_omega()
      scal  = 1.0d0/omega

      do ii=1,ion_nion()
        ia=ion_katm(ii)

        nproj = int_mb(nprj(1)+ia-1)
        if (nproj.gt.0) then


        do nb=1,nbrill
           call Cram_npack(nb,npack)

*       **** structure factor pseudopotential ****
           call cstrfac_pack(nb,ii,dcpl_mb(exi(1)))
           call cstrfac_k(ii,nb,cxr)
           call Cram_c_ZMul(nb,cxr,dcpl_mb(exi(1)),dcpl_mb(exi(1)))


*          **** generate zsw1's and projectors ****
           do l=1,nproj

              shift = vnl(1)+(l-1) *npack1
     >                      +(nb-1)*npack1*nmax_max*lmmax_max
     >                      +(ia-1)*npack1*nmax_max*lmmax_max*nbrill
              l_prj = int_mb(l_projector(1)+(l-1)
     >                                     +(ia-1)*(nmax_max*lmmax_max))

              !sd_function = .not.and(l_prj,1)
#ifdef GCC4
              k = iand(l_prj,1)
#else
              k = and(l_prj,1)
#endif
              sd_function = (k.eq.0)

*             **** phase factor does not matter therefore ****
*             **** (-i)^l is the same as (i)^l in the     ****
*             **** Rayleigh scattering formula            ****

c             *** current function is s or d ****
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

*             **** compute nnXnproj matrix zsw1 = <psi1|prj> ****
              call Cram_cc_inzdot(nb,nn,
     >                      psi1(1+(nb-1)*nn*npack1),
     >                      dcpl_mb(prjtmp(1)+(l-1)*npack1),
     >                      dcpl_mb(zsw1(1)+(l-1)*nn))
           end do !**l**
           call C3dB_Vector_SumAll((2*nn*nproj),dcpl_mb(zsw1(1)))
         

*          **** zsw2 = Gijl*zsw1 ******
           call Multiply_Gijl_zsw1(nn,
     >                         nproj,
     >                         int_mb(nmax(1)+ia-1),
     >                         int_mb(lmax(1)+ia-1),
     >                         int_mb(n_projector(1)
     >                                + (ia-1)*(nmax_max*lmmax_max)),
     >                         int_mb(l_projector(1)
     >                                + (ia-1)*(nmax_max*lmmax_max)),
     >                         int_mb(m_projector(1)
     >                                + (ia-1)*(nmax_max*lmmax_max)),
     >                         dbl_mb(Gijl(1)
     >                         + (ia-1)*nmax_max*nmax_max*(lmax_max+1)),
     >                         dcpl_mb(zsw1(1)),
     >                         dcpl_mb(zsw2(1)))


*          **** do Kleinman-Bylander Multiplication ****
           call dscal(2*nn*nproj,scal,dcpl_mb(zsw2(1)),1)
           call ZGEMM('N','C',npack,nn,nproj,
     >                mone,
     >                dcpl_mb(prjtmp(1)), npack1,
     >                dcpl_mb(zsw2(1)),   nn,
     >                one,
     >                psi2(1+(nb-1)*nn*npack1), npack1)


           if (move) then
              weight = brillioun_weight(nb)
              if (ispin.eq.1) call dscal(2*nn*nproj,dcpl_mb(zsw2(1)),1)
          
              do l=1,nproj
          
                 do n=1,nn

                    do i=1,npack
                    ctmp = psi1(i+(n-1)*npack1
     >                           +(nb-1)*neall*npack1)
     >                   *dconjg(dcpl_mb(prjtmp(1)+(l-1)*npack1+i-1))
     >                   *dconjg(dcpl_mb(zsw2(1)+(l-1)*nn+n-1))
                    dbl_mb(xtmp(1)+i-1) = dimag(ctmp)
                    end do

*                   **** define Gx,Gy and Gz in packed space ****
                    call C3dB_r_Copy(1,dbl_mb(G(1)),dbl_mb(Gx(1)))
                    call C3dB_r_Copy(1,dbl_mb(G(2)),dbl_mb(Gy(1)))
                    call C3dB_r_Copy(1,dbl_mb(G(3)),dbl_mb(Gz(1)))
                    call Cram_r_pack(nb,dbl_mb(Gx(1)))
                    call Cram_r_pack(nb,dbl_mb(Gy(1)))
                    call Cram_r_pack(nb,dbl_mb(Gz(1)))
                    call Cram_rr_idot(nb,dbl_mb(Gx(1)),dbl_mb(xtmp(1)),
     >                                dbl_mb(sum(1)+3*(n-1)))
                    call Cram_rr_idot(nb,dbl_mb(Gy(1)),dbl_mb(xtmp(1)),
     >                                dbl_mb(sum(1)+1+3*(n-1)))
                    call Cram_rr_idot(nb,dbl_mb(Gz(1)),dbl_mb(xtmp(1)),
     >                                dbl_mb(sum(1)+2+3*(n-1)))
                 end do !**n**

                 call C3dB_Vector_Sumall(3*(nn),dbl_mb(sum(1)))
  
                 do n=1,nn
                    fion(1,ii) = fion(1,ii) 
     >                         + 2.0d0*weight
     >                                *dbl_mb(sum(1)+3*(n-1))
                    fion(2,ii) = fion(2,ii) 
     >                         + 2.0d0*weight
     >                                *dbl_mb(sum(1)+1+3*(n-1))
                    fion(3,ii) = fion(3,ii) 
     >                         + 2.0d0*weight
     >                                *dbl_mb(sum(1)+2+3*(n-1))
                 end do !** nn **

              end do !** l **
           end if !** move **
       
        end do !** nb **
        end if !** nproj>0**

      end do !**ii**

      if (move) then
        value = MA_pop_stack(sum(2))
        value = value.and.MA_pop_stack(Gz(2))
        value = value.and.MA_pop_stack(Gy(2))
        value = value.and.MA_pop_stack(Gx(2))
        value = value.and.MA_pop_stack(xtmp(2))
        if (.not.value) call errquit('cpsp_v_nonlocal:popping stack',2,
     &       MA_ERR)
      end if

      value =           MA_pop_stack(zsw2(2))
      value = value.and.MA_pop_stack(zsw1(2))
      value = value.and.MA_pop_stack(exi(2))
      if (.not.value) call errquit('cpsp_v_nonlocal:popping stack',3,
     &       MA_ERR)

      call nwpw_timing_end(6)

      return 
      end



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

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Create the G(i,j,L)*<psi|(LS)|prj(j,L)> array
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine Multiply_Kijl_SO(nn,nprj,nmax,lmax,
     >                             n_prj,l_prj,m_prj,
     >                             G,
     >                             zsw1,zsw2)
      implicit none
      integer nn
      integer nprj,nmax,lmax,nh
      integer n_prj(nprj)
      integer l_prj(nprj)
      integer m_prj(nprj)
      real*8  G(nmax,nmax,0:lmax),xx,dl2,dm2,dlm
      complex*16 zsw1(nn,nprj)
      complex*16 zsw2(nn,nprj)

      !**** local variables ****
      integer a,b,na,nb,la,lb,ma,mb

      nh=nn/2
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
              xx = 0.5d0*G(nb,na,la)*dble(ma)  
              call daxpy(nn,xx,zsw1(1,a),1,zsw2(1,b),1)
              call daxpy(nn,(-xx),zsw1(nh+1,a),1,zsw2(nh+1,b),1)
            end if

         end do
      end do
      do b=1,nprj
         lb = l_prj(b)
         mb = m_prj(b)
         do a=1,nprj
            la = l_prj(a)
            ma = m_prj(a)
            if ((la.eq.lb).and.(ma.eq.(mb+1)).and.(mb.ne.lb)) then
              dl2=dble(la*(la+1))
              dm2=dble(ma*mb)
              dlm=dsqrt(dl2-dm2)
              na = n_prj(a)
              nb = n_prj(b)
              xx = 0.5d0*G(nb,na,la)*dlm
              call daxpy(nn,xx,zsw1(1,a),1,zsw2(nh+1,b),1)
            end if
         end do
      end do
      do b=1,nprj
         lb = l_prj(b)
         mb = m_prj(b)
         do a=1,nprj
            la = l_prj(a)
            ma = m_prj(a)
            if ((la.eq.lb).and.(ma.eq.(mb-1)).and.(mb.ne.(-lb))) then
              dl2=dble(la*(la+1))
              dm2=dble(ma*mb)
              dlm=dsqrt(dl2-dm2)
              na = n_prj(a)
              nb = n_prj(b)
              xx = 0.5d0*G(nb,na,la)*dlm
              call daxpy(nn,xx,zsw1(nh+1,a),1,zsw2(1,b),1)
            end if
         end do
      end do
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Create the G(i,j,L)*<psi|(LS)|prj(j,L)> array
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine Multiply_Kijl_SO_x(nn,nprj,nmax,lmax,
     >                             n_prj,l_prj,m_prj,
     >                             G,
     >                             zsw1u,zsw1d,zsw2u,zsw2d)
      implicit none
      integer nn
      integer nprj,nmax,lmax,nh
      integer n_prj(nprj)
      integer l_prj(nprj)
      integer m_prj(nprj)
      real*8  G(nmax,nmax,0:lmax),xx,dl2,dm2,dlm
      complex*16 zsw1u(nprj),zsw1d(nprj)
      complex*16 zsw2u(nprj),zsw2d(nprj)

      !**** local variables ****
      integer a,b,na,nb,la,lb,ma,mb

      call dcopy(2*nn*nprj,0.0d0,0,zsw2u,1)
      call dcopy(2*nn*nprj,0.0d0,0,zsw2d,1)
      do b=1,nprj
         lb = l_prj(b)
         mb = m_prj(b)

         do a=1,nprj
            la = l_prj(a)
            ma = m_prj(a)

            if ((la.eq.lb).and.(ma.eq.mb)) then
              na = n_prj(a)
              nb = n_prj(b)
              xx = 0.5d0*G(nb,na,la)*dble(ma)  
              call daxpy(nn,xx,zsw1u(a),1,zsw2u(b),1)
              call daxpy(nn,(-xx),zsw1d(a),1,zsw2d(b),1)
            end if

         end do
      end do
      do b=1,nprj
         lb = l_prj(b)
         mb = m_prj(b)
         do a=1,nprj
            la = l_prj(a)
            ma = m_prj(a)
            if ((la.eq.lb).and.(ma.eq.(mb+1)).and.(mb.ne.lb)) then
              dl2=dble(la*(la+1))
              dm2=dble(ma*mb)
              dlm=dsqrt(dl2-dm2)
              na = n_prj(a)
              nb = n_prj(b)
              xx = 0.5d0*G(nb,na,la)*dlm
              call daxpy(nn,xx,zsw1u(a),1,zsw2d(b),1)
            end if
         end do
      end do
      do b=1,nprj
         lb = l_prj(b)
         mb = m_prj(b)
         do a=1,nprj
            la = l_prj(a)
            ma = m_prj(a)
            if ((la.eq.lb).and.(ma.eq.(mb-1)).and.(mb.ne.(-lb))) then
              dl2=dble(la*(la+1))
              dm2=dble(ma*mb)
              dlm=dsqrt(dl2-dm2)
              na = n_prj(a)
              nb = n_prj(b)
              xx = 0.5d0*G(nb,na,la)*dlm
              call daxpy(nn,xx,zsw1d(a),1,zsw2u(b),1)
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

      subroutine cpsp_v_nonlocal_orb(nb,psi1,psi2)
      implicit none
      integer    nb
      complex*16 psi1(*)
      complex*16 psi2(*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "cpsp_common.fh"


*     *** local variables ***
      complex*16 one,mone
      parameter  (one=(1.0d0,0.0d0), mone=(-1.0d0,0.0d0))

      integer nfft3d,npack1,npack,nbrill
      integer ii,ia,l
      integer k,shift,l_prj,nproj
      real*8  omega,scal
      complex*16 cxr
      integer exi(2),zsw1(2),zsw2(2)
      logical value,sd_function

*     **** external functions ****
      logical  is_sORd,cpsi_spin_orbit
      integer  ion_nion,ion_katm,brillioun_nbrillioun
      real*8   lattice_omega
      external is_sORd,cpsi_spin_orbit
      external ion_nion,ion_katm,brillioun_nbrillioun
      external lattice_omega

      if (cpsi_spin_orbit()) then
        call cpsp_v_nonlocal_orb_2com(nb,psi1,psi2)
        return
      end if

      call nwpw_timing_start(6)
      
*     **** allocate local memory ****
      call C3dB_nfft3d(1,nfft3d)
      call Cram_max_npack(npack1)
      nbrill = brillioun_nbrillioun()

      value = MA_push_get(mt_dcpl,npack1,'exi', exi(2), exi(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,nmax_max*lmmax_max,
     >                    'zsw1',zsw1(2),zsw1(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,nmax_max*lmmax_max,
     >                    'zsw2',zsw2(2),zsw2(1))
      if (.not.value) 
     > call errquit('cpsp_v_nonlocal_orb:pushing stack',0,MA_ERR)


      omega = lattice_omega()
      scal  = 1.0d0/omega

      do ii=1,ion_nion()
         ia=ion_katm(ii)
         nproj = int_mb(nprj(1)+ia-1)

         if (nproj.gt.0) then

*           **** structure factor ****
            call Cram_npack(nb,npack)
            call cstrfac_pack(nb,ii,dcpl_mb(exi(1)))
            call cstrfac_k(ii,nb,cxr)
            call Cram_c_ZMul(nb,cxr,dcpl_mb(exi(1)),dcpl_mb(exi(1)))

            do l=1,nproj

              shift = vnl(1)+(l-1) *npack1
     >                      +(nb-1)*npack1*nmax_max*lmmax_max
     >                      +(ia-1)*npack1*nmax_max*lmmax_max*nbrill
              l_prj = int_mb(l_projector(1)+(l-1)
     >                                     +(ia-1)*(nmax_max*lmmax_max))

              !sd_function = .not.and(l_prj,1)
#ifdef GCC4
              k = iand(l_prj,1)
#else
              k = and(l_prj,1)
#endif
              sd_function = (k.eq.0)

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
     >                      psi1,
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
     >                                + (ia-1)*(nmax_max*lmmax_max)),
     >                         int_mb(l_projector(1)
     >                                + (ia-1)*(nmax_max*lmmax_max)),
     >                         int_mb(m_projector(1)
     >                                + (ia-1)*(nmax_max*lmmax_max)),
     >                         dbl_mb(Gijl(1)
     >                         + (ia-1)*nmax_max*nmax_max*(lmax_max+1)),
     >                         dcpl_mb(zsw1(1)),
     >                         dcpl_mb(zsw2(1)))

*           **** do Kleinman-Bylander Multiplication ****
            call dscal(2*nproj,scal,dcpl_mb(zsw2(1)),1)
            call ZGEMM('N','C',npack1,1,nproj,
     >                 mone,
     >                 dcpl_mb(prjtmp(1)), npack1,
     >                 dcpl_mb(zsw2(1)),   1,
     >                 one,
     >                 psi2, npack1)

         end if !** nproj>0 **
      end do !** ii***

      value =           MA_pop_stack(zsw2(2))
      value = value.and.MA_pop_stack(zsw1(2))
      value = value.and.MA_pop_stack(exi(2))
      if (.not.value) 
     > call errquit('cpsp_v_nonlocal_orb:popping stack',3,MA_ERR)

      call nwpw_timing_end(6)

      return 
      end


*     *******************************************
*     *				  		*
*     *	 	 cpsp_v_nonlocal_orb_2com  	*
*     *						*
*     *******************************************

      subroutine cpsp_v_nonlocal_orb_2com(nb,psi1,psi2)
      implicit none
      integer    nb
      complex*16 psi1(*)
      complex*16 psi2(*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "cpsp_common.fh"


*     *** local variables ***
      complex*16 one,mone
      parameter  (one=(1.0d0,0.0d0), mone=(-1.0d0,0.0d0))

      integer nfft3d,npack1,npack,nbrill
      integer ii,ia,l
      integer k,shift,l_prj,nproj,shifts,ne1
      real*8  omega,scal
      complex*16 cxr
      integer exi(2),zsw1u(2),zsw2u(2),zsw1d(2),zsw2d(2)
      logical value,sd_function

*     **** external functions ****
      logical  is_sORd
      integer  ion_nion,ion_katm,brillioun_nbrillioun
      integer  cpsi_ne
      real*8   lattice_omega
      external is_sORd
      external ion_nion,ion_katm,brillioun_nbrillioun
      external lattice_omega,cpsi_ne

      call nwpw_timing_start(6)
      
*     **** allocate local memory ****
      call C3dB_nfft3d(1,nfft3d)
      call Cram_max_npack(npack1)
      nbrill = brillioun_nbrillioun()

      value = MA_push_get(mt_dcpl,npack1,'exi', exi(2), exi(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,nmax_max*lmmax_max,
     >                    'zsw1u',zsw1u(2),zsw1u(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,nmax_max*lmmax_max,
     >                    'zsw2u',zsw2u(2),zsw2u(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,nmax_max*lmmax_max,
     >                    'zsw1d',zsw1d(2),zsw1d(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,nmax_max*lmmax_max,
     >                    'zsw2d',zsw2d(2),zsw2d(1))
      if (.not.value) 
     > call errquit('cpsp_v_nonlocal_orb2com:pushing stack',0,MA_ERR)

      ne1=cpsi_ne(1)
      shifts = npack1*ne1
      omega = lattice_omega()
      scal  = 1.0d0/omega

      do ii=1,ion_nion()
         ia=ion_katm(ii)
         nproj = int_mb(nprj(1)+ia-1)

         if (nproj.gt.0) then

*           **** structure factor ****
            call Cram_npack(nb,npack)
            call cstrfac_pack(nb,ii,dcpl_mb(exi(1)))
            call cstrfac_k(ii,nb,cxr)
            call Cram_c_ZMul(nb,cxr,dcpl_mb(exi(1)),dcpl_mb(exi(1)))

            do l=1,nproj

              shift = vnl(1)+(l-1) *npack1
     >                      +(nb-1)*npack1*nmax_max*lmmax_max
     >                      +(ia-1)*npack1*nmax_max*lmmax_max*nbrill
              l_prj = int_mb(l_projector(1)+(l-1)
     >                                     +(ia-1)*(nmax_max*lmmax_max))

              !sd_function = .not.and(l_prj,1)
#ifdef GCC4
              k = iand(l_prj,1)
#else
              k = and(l_prj,1)
#endif
              sd_function = (k.eq.0)

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
     >                      psi1,
     >                      dcpl_mb(prjtmp(1)+(l-1)*npack1),
     >                      dcpl_mb(zsw1u(1)+(l-1)))
              call Cram_cc_izdot(nb,
     >                      psi1(shifts+1),
     >                      dcpl_mb(prjtmp(1)+(l-1)*npack1),
     >                      dcpl_mb(zsw1d(1)+(l-1)))
            end do !**l**
            call C3dB_Vector_SumAll((2*nproj),dcpl_mb(zsw1u(1)))
            call C3dB_Vector_SumAll((2*nproj),dcpl_mb(zsw1d(1)))

*           **** zsw2 = Gijl*zsw1 ******
            call Multiply_Gijl_zsw1(1,
     >                         nproj,
     >                         int_mb(nmax(1)+ia-1),
     >                         int_mb(lmax(1)+ia-1),
     >                         int_mb(n_projector(1)
     >                                + (ia-1)*(nmax_max*lmmax_max)),
     >                         int_mb(l_projector(1)
     >                                + (ia-1)*(nmax_max*lmmax_max)),
     >                         int_mb(m_projector(1)
     >                                + (ia-1)*(nmax_max*lmmax_max)),
     >                         dbl_mb(Gijl(1)
     >                         + (ia-1)*nmax_max*nmax_max*(lmax_max+1)),
     >                         dcpl_mb(zsw1u(1)),
     >                         dcpl_mb(zsw2u(1)))
            call Multiply_Gijl_zsw1(1,
     >                         nproj,
     >                         int_mb(nmax(1)+ia-1),
     >                         int_mb(lmax(1)+ia-1),
     >                         int_mb(n_projector(1)
     >                                + (ia-1)*(nmax_max*lmmax_max)),
     >                         int_mb(l_projector(1)
     >                                + (ia-1)*(nmax_max*lmmax_max)),
     >                         int_mb(m_projector(1)
     >                                + (ia-1)*(nmax_max*lmmax_max)),
     >                         dbl_mb(Gijl(1)
     >                         + (ia-1)*nmax_max*nmax_max*(lmax_max+1)),
     >                         dcpl_mb(zsw1d(1)),
     >                         dcpl_mb(zsw2d(1)))

*           **** do Kleinman-Bylander Multiplication ****
            call dscal(2*nproj,scal,dcpl_mb(zsw2u(1)),1)
            call dscal(2*nproj,scal,dcpl_mb(zsw2d(1)),1)
            call ZGEMM('N','C',npack1,1,nproj,
     >                 mone,
     >                 dcpl_mb(prjtmp(1)), npack1,
     >                 dcpl_mb(zsw2u(1)),   1,
     >                 one,
     >                 psi2, npack1)
            call ZGEMM('N','C',npack1,1,nproj,
     >                 mone,
     >                 dcpl_mb(prjtmp(1)), npack1,
     >                 dcpl_mb(zsw2d(1)),   1,
     >                 one,
     >                 psi2(1+shifts), npack1)

         end if !** nproj>0 **
      end do !** ii***

      value =           MA_pop_stack(zsw2d(2))
      value = value.and.MA_pop_stack(zsw1d(2))
      value = value.and.MA_pop_stack(zsw2u(2))
      value = value.and.MA_pop_stack(zsw1u(2))
      value = value.and.MA_pop_stack(exi(2))
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

      subroutine cpsp_f_vnonlocal(ispin,ne,psi1,fion)
      implicit none
      integer    ispin,ne(2)
      complex*16 psi1(*)
      real*8 fion(3,*)

#include "mafdecls.fh"
#include "cpsp_common.fh"
#include "errquit.fh"


*     *** local variables ***
      integer nfft3d,G(3),npack1,npack,nbrill,neall,shift
      integer i,ii,ia,k,l,n,nn,nb,l_prj,nproj
      real*8  omega,weight,scal
      complex*16 ctmp,cxr
      integer exi(2),xtmp(2),zsw1(2),zsw2(2),sum(2)
      integer Gx(2),Gy(2),Gz(2)
      logical value,sd_function

*     **** external functions ****
      logical  is_sORd
      integer  ion_nion,ion_katm,c_G_indx,brillioun_nbrillioun
      real*8   lattice_omega,brillioun_weight
      external is_sORd
      external ion_nion,ion_katm,c_G_indx,brillioun_nbrillioun
      external lattice_omega,brillioun_weight


      call nwpw_timing_start(6)
      
*     **** allocate local memory ****
      nn = ne(1)+ne(2)
      neall = ne(1)+ne(2)
      nbrill = brillioun_nbrillioun()
      call C3dB_nfft3d(1,nfft3d)
      call Cram_max_npack(npack1)

      value = MA_push_get(mt_dcpl,npack1,'exi', exi(2), exi(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,nn*nmax_max*lmmax_max,
     >                    'zsw1',zsw1(2),zsw1(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,nn*nmax_max*lmmax_max,
     >                    'zsw2',zsw2(2),zsw2(1))
      if (.not.value) call errquit('cpsp_v_nonlocal:pushing stack',0,
     &       MA_ERR)

      value = MA_push_get(mt_dbl, nfft3d,'xtmp',xtmp(2),xtmp(1))
      value = value.and.
     >        MA_push_get(mt_dbl, nfft3d,'Gx',Gx(2),Gx(1))
      value = value.and.
     >        MA_push_get(mt_dbl, nfft3d,'Gy',Gy(2),Gy(1))
      value = value.and.
     >        MA_push_get(mt_dbl, nfft3d,'Gz',Gz(2),Gz(1))

      value = value.and.
     >      MA_push_get(mt_dbl,3*nn,'sum',sum(2),sum(1))
      if (.not.value) call errquit('cpsp_v_nonlocal:pushing stack',1,
     &       MA_ERR)
      G(1)  = c_G_indx(1)
      G(2)  = c_G_indx(2)
      G(3)  = c_G_indx(3)

      omega = lattice_omega()
      scal  = 1.0d0/omega

      do nb = 1,nbrill
        call Cram_npack(nb,npack)
        weight = brillioun_weight(nb)

*       **** define Gx,Gy and Gz in packed space ****
        call C3dB_r_Copy(1,dbl_mb(G(1)),dbl_mb(Gx(1)))
        call C3dB_r_Copy(1,dbl_mb(G(2)),dbl_mb(Gy(1)))
        call C3dB_r_Copy(1,dbl_mb(G(3)),dbl_mb(Gz(1)))
        call Cram_r_pack(nb,dbl_mb(Gx(1)))
        call Cram_r_pack(nb,dbl_mb(Gy(1)))
        call Cram_r_pack(nb,dbl_mb(Gz(1)))

        do ii=1,ion_nion()
           ia=ion_katm(ii)
           nproj = int_mb(nprj(1)+ia-1)
           if (nproj.gt.0) then

*             **** structure factor ****
              call cstrfac_pack(nb,ii,dcpl_mb(exi(1)))
              call cstrfac_k(ii,nb,cxr)
              call Cram_c_ZMul(nb,cxr,dcpl_mb(exi(1)),dcpl_mb(exi(1)))

              do l=1,nproj
                 shift = vnl(1)
     >                      +(l-1) *npack1
     >                      +(nb-1)*npack1*nmax_max*lmmax_max
     >                      +(ia-1)*npack1*nmax_max*lmmax_max*nbrill
                 l_prj = int_mb(l_projector(1)+(l-1)
     >                          +(ia-1)*(nmax_max*lmmax_max))
              
     >       
                 !sd_function = .not.and(l_prj,1)
#ifdef GCC4
                 k = iand(l_prj,1)
#else
                 k = and(l_prj,1)
#endif
                 sd_function = (k.eq.0)

*                **** phase factor does not matter therefore ****
*                **** (-i)^l is the same as (i)^l in the     ****
*                **** Rayleigh scattering formula            ****
*                *** current function is s or d ****
                 if (sd_function) then
                    call Cram_rc_Mul(nb,dbl_mb(shift),
     >                                  dcpl_mb(exi(1)),
     >                                  dcpl_mb(prjtmp(1)+(l-1)*npack1))
*                *** current function is p or f ****
                 else
                    call Cram_irc_Mul(nb,dbl_mb(shift),
     >                                   dcpl_mb(exi(1)),
     >                                  dcpl_mb(prjtmp(1)+(l-1)*npack1))
                 end if

*                **** compute nnXnproj matrix zsw1 = <psi1|prj> ****
                 call Cram_cc_inzdot(nb,nn,
     >                         psi1(1+(nb-1)*neall*npack1),
     >                         dcpl_mb(prjtmp(1)+(l-1)*npack1),
     >                         dcpl_mb(zsw1(1)+(l-1)*nn))
              end do
              call C3dB_Vector_SumAll((2*nn*nproj),dcpl_mb(zsw1(1)))


*             **** do kleinman-bylander multiplication ****
*             **** sw2 = Gijl*sw1 ******
              call Multiply_Gijl_zsw1(nn,
     >                          nproj,
     >                          int_mb(nmax(1)+ia-1),
     >                          int_mb(lmax(1)+ia-1),
     >                          int_mb(n_projector(1)
     >                                 + (ia-1)*(nmax_max*lmmax_max)),
     >                          int_mb(l_projector(1)
     >                                 + (ia-1)*(nmax_max*lmmax_max)),
     >                          int_mb(m_projector(1)
     >                                 + (ia-1)*(nmax_max*lmmax_max)),
     >                          dbl_mb(Gijl(1)
     >                         + (ia-1)*nmax_max*nmax_max*(lmax_max+1)),
     >                          dcpl_mb(zsw1(1)),
     >                          dcpl_mb(zsw2(1)))

              call dscal(2*nn*nproj,scal,dcpl_mb(zsw2(1)),1)
              if (ispin.eq.1) call dscal(2*nn*nproj,dcpl_mb(zsw2(1)),1)

              do l=1,nproj

                 do n=1,nn
                    do i=1,npack
                       ctmp = psi1(i+(n-1)*npack1
     >                              +(nb-1)*neall*npack1)
     >                      *dconjg(dcpl_mb(prjtmp(1)+(l-1)*npack1+i-1))
     >                      *dconjg(dcpl_mb(zsw2(1)+(l-1)*nn+n-1))
                       dbl_mb(xtmp(1)+i-1) = dimag(ctmp)
                    end do
 
                    call Cram_rr_idot(nb,dbl_mb(Gx(1)),dbl_mb(xtmp(1)),
     >                                   dbl_mb(sum(1)+3*(n-1)))
                    call Cram_rr_idot(nb,dbl_mb(Gy(1)),dbl_mb(xtmp(1)),
     >                                dbl_mb(sum(1)+1+3*(n-1)))
                    call Cram_rr_idot(nb,dbl_mb(Gz(1)),dbl_mb(xtmp(1)),
     >                                dbl_mb(sum(1)+2+3*(n-1)))
                 end do
                 call C3dB_Vector_Sumall(3*(nn),dbl_mb(sum(1)))
  
                 do n=1,nn
                    fion(1,ii) = fion(1,ii) 
     >                         + 2.0d0*weight*dbl_mb(sum(1)+3*(n-1))
                    fion(2,ii) = fion(2,ii) 
     >                         + 2.0d0*weight*dbl_mb(sum(1)+1+3*(n-1))
                    fion(3,ii) = fion(3,ii) 
     >                         + 2.0d0*weight*dbl_mb(sum(1)+2+3*(n-1))
                 end do
        
              end do !** l **

           end if !** nproj>0) ***
        end do !** ii **

      end do !** nb ***

      value =           MA_pop_stack(sum(2))
      value = value.and.MA_pop_stack(Gz(2))
      value = value.and.MA_pop_stack(Gy(2))
      value = value.and.MA_pop_stack(Gx(2))
      value = value.and.MA_pop_stack(xtmp(2))
      if (.not.value) call errquit('cpsp_f_vnonlocal:popping stack',2,
     &       MA_ERR)

      value =           MA_pop_stack(zsw2(2))
      value = value.and.MA_pop_stack(zsw1(2))
      value = value.and.MA_pop_stack(exi(2))
      if (.not.value) call errquit('cpsp_f_vnonlocal:popping stack',3,
     &       MA_ERR)

      call nwpw_timing_end(6)
      

      return 
      end


*     ***********************************
*     *					*
*     *	 	 cpsp_v_spin_orbit  	*
*     *					*
*     ***********************************

      subroutine cpsp_v_spin_orbit(ispin,ne,
     >                           psi1,psi2,move,fion)
      implicit none
      integer    ispin,ne(2)
      complex*16 psi1(*)
      complex*16 psi2(*)
      logical move
      real*8 fion(3,*)

#include "mafdecls.fh"
#include "cpsp_common.fh"
#include "errquit.fh"

*     *** local variables ***
      complex*16 one,mone
      parameter  (one=(1.0d0,0.0d0), mone=(-1.0d0,0.0d0))

      integer nfft3d,G(3),npack1,npack
      integer i,ii,ia,l,n,nn,nb,neall,nbrill
      integer k,shift,l_prj,nproj
      real*8  omega,weight,scal
      complex*16 ctmp,cxr
      integer exi(2),xtmp(2),zsw1(2),zsw2(2),sum(2)
      integer Gx(2),Gy(2),Gz(2)
      logical value,sd_function

*     **** external functions ****
      logical  is_sORd
      integer  ion_nion,ion_katm,c_G_indx,brillioun_nbrillioun
      real*8   lattice_omega,brillioun_weight
      external is_sORd
      external ion_nion,ion_katm,c_G_indx,brillioun_nbrillioun
      external lattice_omega,brillioun_weight

      if (.not.do_spin_orbit) return
      call nwpw_timing_start(6)
      
*     **** allocate local memory ****
      nn = ne(1)+ne(2)
      neall = ne(1)+ne(2)
      nbrill = brillioun_nbrillioun()
      call C3dB_nfft3d(1,nfft3d)
      call Cram_max_npack(npack1)

      value = MA_push_get(mt_dcpl,npack1,'exi', exi(2), exi(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,nn*nmax_max*lmmax_max,
     >                    'zsw1',zsw1(2),zsw1(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,nn*nmax_max*lmmax_max,
     >                    'zsw2',zsw2(2),zsw2(1))
      if (.not.value) call errquit('cpsp_v_nonlocal:pushing stack',0,
     &       MA_ERR)

      if (move) then
        value = MA_push_get(mt_dbl, nfft3d,'xtmp',xtmp(2),xtmp(1))
        value = value.and.
     >          MA_push_get(mt_dbl, nfft3d,'Gx',Gx(2),Gx(1))
        value = value.and.
     >          MA_push_get(mt_dbl, nfft3d,'Gy',Gy(2),Gy(1))
        value = value.and.
     >          MA_push_get(mt_dbl, nfft3d,'Gz',Gz(2),Gz(1))
        value = value.and.
     >        MA_push_get(mt_dbl,3*nn,'sum',sum(2),sum(1))
        if (.not.value) call errquit('cpsp_v_nonlocal:pushing stack',1,
     &       MA_ERR)
        G(1)  = c_G_indx(1)
        G(2)  = c_G_indx(2)
        G(3)  = c_G_indx(3)

*       **** define Gx,Gy and Gz in packed space ****
        call C3dB_r_Copy(1,dbl_mb(G(1)),dbl_mb(Gx(1)))
        call C3dB_r_Copy(1,dbl_mb(G(2)),dbl_mb(Gy(1)))
        call C3dB_r_Copy(1,dbl_mb(G(3)),dbl_mb(Gz(1)))
      end if

      omega = lattice_omega()
      scal  = 1.0d0/omega

      do ii=1,ion_nion()
        ia=ion_katm(ii)

        nproj = int_mb(nprj(1)+ia-1)
        if (nproj.gt.0) then


        do nb=1,nbrill
           call Cram_npack(nb,npack)

*       **** structure factor pseudopotential ****
           call cstrfac_pack(nb,ii,dcpl_mb(exi(1)))
           call cstrfac_k(ii,nb,cxr)
           call Cram_c_ZMul(nb,cxr,dcpl_mb(exi(1)),dcpl_mb(exi(1)))


*          **** generate zsw1's and projectors ****
           do 105 l=1,nproj

              shift = vnlso(1)+(l-1) *npack1
     >                      +(nb-1)*npack1*nmax_max*lmmax_max
     >                      +(ia-1)*npack1*nmax_max*lmmax_max*nbrill
              l_prj = int_mb(l_projector(1)+(l-1)
     >                                     +(ia-1)*(nmax_max*lmmax_max))
              if (l_prj.eq.0) then
                call dcopy(npack1*2,0.0d0,0,
     >               dcpl_mb(prjtmp(1)+(l-1)*npack1),1)
                goto 105
              end if
                sd_function = ((l_prj.eq.0).or.(l_prj.eq.2))

*             **** phase factor does not matter therefore ****
*             **** (-i)^l is the same as (i)^l in the     ****
*             **** Rayleigh scattering formula            ****

c             *** current function is s or d ****
                if (sd_function) then
                  call Cram_cc_Mul(nb,dcpl_mb(shift),
     >                               dcpl_mb(exi(1)),
     >                               dcpl_mb(prjtmp(1)+(l-1)*npack1))
*             *** current function is p or f ****
                else
                  call Cram_icc_Mul(nb,dcpl_mb(shift),
     >                                dcpl_mb(exi(1)),
     >                                dcpl_mb(prjtmp(1)+(l-1)*npack1))
                end if

*             **** compute nnXnproj matrix zsw1 = <psi1|prj> ****
                call Cram_cc_inzdot(nb,nn,
     >                      psi1(1+(nb-1)*nn*npack1),
     >                      dcpl_mb(prjtmp(1)+(l-1)*npack1),
     >                      dcpl_mb(zsw1(1)+(l-1)*nn))
           
 105       continue 
           call C3dB_Vector_SumAll((2*nn*nproj),dcpl_mb(zsw1(1)))
         

*          **** zsw2 = Kijl*zsw1 ******
           call Multiply_Kijl_SO(nn,
     >                         nproj,
     >                         int_mb(nmax(1)+ia-1),
     >                         int_mb(lmax(1)+ia-1),
     >                         int_mb(n_projector(1)
     >                                + (ia-1)*(nmax_max*lmmax_max)),
     >                         int_mb(l_projector(1)
     >                                + (ia-1)*(nmax_max*lmmax_max)),
     >                         int_mb(m_projector(1)
     >                                + (ia-1)*(nmax_max*lmmax_max)),
     >                         dbl_mb(Kijl(1)
     >                         + (ia-1)*nmax_max*nmax_max*(lmax_max+1)),
     >                         dcpl_mb(zsw1(1)),
     >                         dcpl_mb(zsw2(1)))


*          **** do Kleinman-Bylander Multiplication ****
           call dscal(2*nn*nproj,scal,dcpl_mb(zsw2(1)),1)
           call ZGEMM('N','C',npack,nn,nproj,
     >                mone,
     >                dcpl_mb(prjtmp(1)), npack1,
     >                dcpl_mb(zsw2(1)),   nn,
     >                one,
     >                psi2(1+(nb-1)*nn*npack1), npack1)


           if (move) then
              weight = brillioun_weight(nb)
              if (ispin.eq.1) call dscal(2*nn*nproj,dcpl_mb(zsw2(1)),1)
          
              do l=1,nproj
          
                 do n=1,nn

                    do i=1,npack
                    ctmp = psi1(i+(n-1)*npack1
     >                           +(nb-1)*neall*npack1)
     >                   *dconjg(dcpl_mb(prjtmp(1)+(l-1)*npack1+i-1))
     >                   *dconjg(dcpl_mb(zsw2(1)+(l-1)*nn+n-1))
                    dbl_mb(xtmp(1)+i-1) = dimag(ctmp)
                    end do

*                   **** define Gx,Gy and Gz in packed space ****
                    call C3dB_r_Copy(1,dbl_mb(G(1)),dbl_mb(Gx(1)))
                    call C3dB_r_Copy(1,dbl_mb(G(2)),dbl_mb(Gy(1)))
                    call C3dB_r_Copy(1,dbl_mb(G(3)),dbl_mb(Gz(1)))
                    call Cram_r_pack(nb,dbl_mb(Gx(1)))
                    call Cram_r_pack(nb,dbl_mb(Gy(1)))
                    call Cram_r_pack(nb,dbl_mb(Gz(1)))
                    call Cram_rr_idot(nb,dbl_mb(Gx(1)),dbl_mb(xtmp(1)),
     >                                dbl_mb(sum(1)+3*(n-1)))
                    call Cram_rr_idot(nb,dbl_mb(Gy(1)),dbl_mb(xtmp(1)),
     >                                dbl_mb(sum(1)+1+3*(n-1)))
                    call Cram_rr_idot(nb,dbl_mb(Gz(1)),dbl_mb(xtmp(1)),
     >                                dbl_mb(sum(1)+2+3*(n-1)))
                 end do !**n**

                 call C3dB_Vector_Sumall(3*(nn),dbl_mb(sum(1)))
  
                 do n=1,nn
                    fion(1,ii) = fion(1,ii) 
     >                         + 2.0d0*weight
     >                                *dbl_mb(sum(1)+3*(n-1))
                    fion(2,ii) = fion(2,ii) 
     >                         + 2.0d0*weight
     >                                *dbl_mb(sum(1)+1+3*(n-1))
                    fion(3,ii) = fion(3,ii) 
     >                         + 2.0d0*weight
     >                                *dbl_mb(sum(1)+2+3*(n-1))
                 end do !** nn **

              end do !** l **
           end if !** move **
       
        end do !** nb **
        end if !** nproj>0**

      end do !**ii**

      if (move) then
        value = MA_pop_stack(sum(2))
        value = value.and.MA_pop_stack(Gz(2))
        value = value.and.MA_pop_stack(Gy(2))
        value = value.and.MA_pop_stack(Gx(2))
        value = value.and.MA_pop_stack(xtmp(2))
        if (.not.value) call errquit('cpsp_v_nonlocal:popping stack',2,
     &       MA_ERR)
      end if

      value =           MA_pop_stack(zsw2(2))
      value = value.and.MA_pop_stack(zsw1(2))
      value = value.and.MA_pop_stack(exi(2))
      if (.not.value) call errquit('cpsp_v_nonlocal:popping stack',3,
     &       MA_ERR)

      call nwpw_timing_end(6)

      return 
      end


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
     >                       vl,vnl,vnlsox,
     >                       semicore,rcore,ncore,
     >                       tmp,tmp2,ztmp,ztmp2,
     >                       ierr)
      implicit none 
      character*(*) comment
      character*20 fname
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
      real*8 vnl(npack1,lmmax_max*nmax_max,*)
      complex*16 vnlsox(npack1,lmmax_max*nmax_max,*)
      logical semicore
      real*8  rcore
      real*8  ncore(*)
      real*8     tmp(*)
      real*8     tmp2(*)
      complex*16  ztmp(*)
      complex*16  ztmp2(*)
      integer ierr

#ifdef TCGMSG
#include "tcgmsg.fh"
#include "msgtypesf.h"
#endif

*    *** local variables ***
      integer MASTER,taskid
      parameter(MASTER=0)
      integer n,l,nb
      integer msglen
      integer iatom(2)
      character*255 full_filename
      real*8 kv(3)
      integer nbrillioun

      integer  brillioun_nbrillioun
      real*8   brillioun_k
      external brillioun_nbrillioun
      external brillioun_k


      call Parallel_taskid(taskid)
  
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
         call dread(5,rc,lmax+1)

         call iread(5,nprj,1)
         if (nprj.gt.0) then
         call iread(5,n_projector,nprj)
         call iread(5,l_projector,nprj)
         call iread(5,m_projector,nprj)
         call dread(5,Gijl,nmax*nmax*(lmax+1)) !** number of matrix elements = nmax*nmax*(lmax+1) **
         if (psp_type.eq.1) then
         call dread(5,Kijl,nmax*nmax*(lmax+1)) !** number of matrix elements = nmax*nmax*(lmax+1) **
         end if
         end if

         call dread(5,rcore,1)
         call iread(5,nbrillioun,1)
         ierr = 0
         if (nbrillioun.eq.brillioun_nbrillioun()) then
            do nb=1,nbrillioun
               call dread(5,kv,3)
               if ((brillioun_k(1,nb).ne.kv(1)).or.
     >             (brillioun_k(2,nb).ne.kv(2)).or.
     >             (brillioun_k(3,nb).ne.kv(3))) ierr = 1
            end do
         else
            ierr = 1
         end if
      end if

      msglen = 1
      call BRDCST(9+MSGINT,ierr,mitob(msglen),MASTER)
      if (ierr.ne.0) then
         return
      end if

*     **** send header data to all processors ****
#ifdef TCGMSG
      msglen = 1
      call BRDCST(9+MSGINT,version,mitob(msglen),MASTER)
      call BRDCST(9+MSGINT,psp_type,mitob(msglen),MASTER)
      msglen = 3
      call BRDCST(9+MSGINT,nfft,mitob(msglen),MASTER)
      msglen = 9
      call BRDCST(9+MSGDBL,unita,mdtob(msglen),MASTER)

      iatom(1) = ichar(atom(1:1))
      iatom(2) = ichar(atom(2:2))
      msglen = 2
      call BRDCST(9+MSGCHR,iatom,mitob(msglen),MASTER)
      atom(1:1) = char(iatom(1))
      atom(2:2) = char(iatom(2))

      msglen = 1
      call BRDCST(9+MSGDBL,amass,mdtob(msglen),MASTER)
      msglen = 1
      call BRDCST(9+MSGDBL,zv,mdtob(msglen),MASTER)
      msglen = 1
      call BRDCST(9+MSGINT,lmax,mitob(msglen),MASTER)
      call BRDCST(9+MSGINT,locp,mitob(msglen),MASTER)
      call BRDCST(9+MSGINT,nmax,mitob(msglen),MASTER)
      call BRDCST(9+MSGINT,nprj,mitob(msglen),MASTER)
      lmmax=(lmax+1)**2 - (2*locp+1)

      msglen=lmax+1
      call BRDCST(9+MSGDBL,rc,mdtob(msglen),MASTER)


      msglen=nprj
      call BRDCST(9+MSGINT,n_projector,mitob(msglen),MASTER)
      call BRDCST(9+MSGINT,l_projector,mitob(msglen),MASTER)
      call BRDCST(9+MSGINT,m_projector,mitob(msglen),MASTER)

      msglen=nmax*nmax*(lmax+1)
      call BRDCST(9+MSGDBL,Gijl,mdtob(msglen),MASTER)
      call BRDCST(9+MSGDBL,Kijl,mdtob(msglen),MASTER)

      msglen=1
      call BRDCST(9+MSGDBL,rcore,mdtob(msglen),MASTER)

#endif

*     **** determine semicore value ****
      if (rcore.gt.0.0d0) then
         semicore = .true.
      else
         semicore = .false.
      end if


 
*     *** read in vl 3d block ***
      call C3dB_r_read(1,5,tmp2,tmp)
      call Cram_r_pack(0,tmp2)
      call Cram_r_Copy(0,tmp2,vl)

*     **** read in vnl 3d blocks ****
      do nb=1,brillioun_nbrillioun()
      do n=1,nprj
         call C3dB_r_read(1,5,tmp2,tmp)
         call Cram_r_pack(nb,tmp2)
         call Cram_r_Copy(nb,tmp2,vnl(1,n,nb))
      end do
      end do
*     **** read in v_spin_orbit 3d blocks ****
      if (psp_type.eq.1) then
      do nb=1,brillioun_nbrillioun()
      do n=1,nprj
         call C3dB_c_read(1,5,ztmp2,ztmp)
         call Cram_c_pack(nb,ztmp2)
         call Cram_c_Copy(nb,ztmp2,vnlsox(1,n,nb))
      end do
      end do
      end if
*     **** read in semicore density block ****
      if (semicore) then
         call C3dB_r_read(1,5,tmp2,tmp)
         call Cram_r_pack(0,tmp2)
         call Cram_r_Copy(0,tmp2,ncore(1))

         call C3dB_r_read(1,5,tmp2,tmp)
         call Cram_r_pack(0,tmp2)
         call Cram_r_Copy(0,tmp2,ncore(1+2*npack0))

         call C3dB_r_read(1,5,tmp2,tmp)
         call Cram_r_pack(0,tmp2)
         call Cram_r_Copy(0,tmp2,ncore(1+3*npack0))

         call C3dB_r_read(1,5,tmp2,tmp)
         call Cram_r_pack(0,tmp2)
         call Cram_r_Copy(0,tmp2,ncore(1+4*npack0))
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
*     *	 	  cpsp_readall  	*
*     *					*
*     ***********************************

      subroutine cpsp_readall()
      implicit none
      
#include "mafdecls.fh"
#include "cpsp_common.fh"
#include "c_semicore_common.fh"
#include "stdio.fh"
#include "errquit.fh"



*     **** local variables ****
      integer ngp(3),version,nfft3d,npack0,npack1,nbrill
      integer ia,l
      real*8 unita(3,3)
      character*12 boundry
      integer tmp(2),tmp2(2),ztmp(2),ztmp2(2),ierr
      logical value,found,correct_box
      character*4  element
      character*20 fname

*     **** parallel i/o variable ****
      integer MASTER,taskid
      parameter(MASTER=0)

*     **** external functions ****
      logical      nwpw_filefind,control_spin_orbit
      integer      control_ngrid,brillioun_nbrillioun
      real*8       control_unita
      character*12 control_boundry
      character*2  ion_atom
      external     nwpw_filefind,control_spin_orbit
      external     control_ngrid,brillioun_nbrillioun
      external     control_unita
      external     control_boundry
      external     ion_atom
      

      call C3dB_nfft3d(1,nfft3d)
      call Cram_npack(0,npack0)
      call Cram_max_npack(npack1)
      call Parallel_taskid(taskid)
      nbrill = brillioun_nbrillioun()

*     *** set semicore(0) *****
      log_mb(semicore(1)) = .false.

      value = MA_push_get(mt_dbl,(2*nfft3d),'tmp',tmp(2),tmp(1))
      if (.not. value) call errquit('out of stack memory',0,
     &       MA_ERR)

      value = MA_push_get(mt_dbl,(nfft3d),'tmp2',tmp2(2),tmp2(1))
      if (.not. value) call errquit('out of stack memory',0,
     &       MA_ERR)

      value = MA_push_get(mt_dcpl,(2*nfft3d),'ztmp',ztmp(2),ztmp(1))
      if (.not. value) call errquit('out of stack memory',0,
     &       MA_ERR)

      value = MA_push_get(mt_dcpl,(nfft3d),'ztmp2',ztmp2(2),ztmp2(1))
      if (.not. value) call errquit('out of stack memory',0,
     &       MA_ERR)

      do_spin_orbit=.false.
*     **** read pseudopotentials ****
      do ia=1,npsp

*      **** define formatted psp name ****
       element = '    '
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
     >                         + (ia-1)*(nmax_max*lmmax_max)),
     >                  int_mb(l_projector(1)
     >                         + (ia-1)*(nmax_max*lmmax_max)),
     >                  int_mb(m_projector(1)
     >                         + (ia-1)*(nmax_max*lmmax_max)),
     >                  dbl_mb(Gijl(1)
     >                         + (ia-1)*nmax_max*nmax_max*(lmax_max+1)),
     >                  dbl_mb(Kijl(1)
     >                         + (ia-1)*nmax_max*nmax_max*(lmax_max+1)),
     >                  nfft3d,npack0,npack1,lmmax_max,nmax_max,
     >                  dbl_mb(vl(1) + (ia-1)*npack0),
     >                  dbl_mb(vnl(1)
     >                        +(ia-1)*npack1*nmax_max*lmmax_max*nbrill),
     >                  dcpl_mb(vnlso(1)
     >                        +(ia-1)*npack1*nmax_max*lmmax_max*nbrill),
     >                  log_mb(semicore(1)+ia),
     >                  dbl_mb(rcore(1)+ia-1),
     >                  dbl_mb(ncore(1)+ (ia-1)*npack0*5),
     >                  dbl_mb(tmp(1)),dbl_mb(tmp2(1)),
     >                  dcpl_mb(ztmp(1)),dcpl_mb(ztmp2(1)),
     >                  ierr)

           if ((int_mb(psp_type(1)+(ia-1)).eq.1).and.
     >           control_spin_orbit()) then
              do_spin_orbit=.true.
              write(*,*)"SPIN ORBIT PPOT ENABLED!"
             call flush(6)
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
              if (taskid.eq.MASTER) then
              write(luout,*)
              write(luout,*) 
     >        "pseudopotential is not correctly formatted:",fname
              end if
           end if
           if (correct_box) found = .true.
           if (ierr.eq.1)   then
              found = .false.
              if (taskid.eq.MASTER) then
              write(luout,*)
              write(luout,*) 
     >        "pseudopotential is not correctly formatted-",
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
 9000 value = MA_pop_stack(ztmp2(2))
      value = MA_pop_stack(ztmp(2)) 
      value = MA_pop_stack(tmp2(2))
      value = MA_pop_stack(tmp(2))

      return
      end

*     *******************************************
*     *				  		*
*     *	 	 cpsp_v_spin_orbit_orb          *
*     *						*
*     *******************************************

      subroutine cpsp_v_spin_orbit_orb(nb,psi1,psi2)
      implicit none
      integer    nb
      complex*16 psi1(*)
      complex*16 psi2(*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "cpsp_common.fh"


*     *** local variables ***
      complex*16 one,mone
      parameter  (one=(1.0d0,0.0d0), mone=(-1.0d0,0.0d0))

      integer nfft3d,npack1,npack,nbrill
      integer ii,ia,l
      integer k,shift,l_prj,nproj,shifts,ne1
      real*8  omega,scal
      complex*16 cxr
      integer exi(2),zsw1u(2),zsw2u(2),zsw1d(2),zsw2d(2)
      logical value,sd_function

*     **** external functions ****
      logical  is_sORd
      integer  ion_nion,ion_katm,brillioun_nbrillioun
      integer  cpsi_ne
      real*8   lattice_omega
      external is_sORd
      external ion_nion,ion_katm,brillioun_nbrillioun
      external lattice_omega,cpsi_ne

      call nwpw_timing_start(6)
      
*     **** allocate local memory ****
      call C3dB_nfft3d(1,nfft3d)
      call Cram_max_npack(npack1)
      nbrill = brillioun_nbrillioun()

      value = MA_push_get(mt_dcpl,npack1,'exi', exi(2), exi(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,nmax_max*lmmax_max,
     >                    'zsw1u',zsw1u(2),zsw1u(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,nmax_max*lmmax_max,
     >                    'zsw2u',zsw2u(2),zsw2u(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,nmax_max*lmmax_max,
     >                    'zsw1d',zsw1d(2),zsw1d(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,nmax_max*lmmax_max,
     >                    'zsw2d',zsw2d(2),zsw2d(1))
      if (.not.value) 
     > call errquit('cpsp_v_spin_orbit_orb2com:pushing stack',0,MA_ERR)

      ne1=cpsi_ne(1)
      shifts = npack1*ne1
      omega = lattice_omega()
      scal  = 1.0d0/omega

      do ii=1,ion_nion()
         ia=ion_katm(ii)
         nproj = int_mb(nprj(1)+ia-1)

         if (nproj.gt.0) then

*           **** structure factor ****
            call Cram_npack(nb,npack)
            call cstrfac_pack(nb,ii,dcpl_mb(exi(1)))
            call cstrfac_k(ii,nb,cxr)
            call Cram_c_ZMul(nb,cxr,dcpl_mb(exi(1)),dcpl_mb(exi(1)))

            do l=1,nproj

              shift = vnlso(1)+(l-1) *npack1
     >                      +(nb-1)*npack1*nmax_max*lmmax_max
     >                      +(ia-1)*npack1*nmax_max*lmmax_max*nbrill
              l_prj = int_mb(l_projector(1)+(l-1)
     >                                     +(ia-1)*(nmax_max*lmmax_max))

              !sd_function = .not.and(l_prj,1)
#ifdef GCC4
              k = iand(l_prj,1)
#else
              k = and(l_prj,1)
#endif
              sd_function = (k.eq.0)

*             **** phase factor does not matter therefore ****
*             **** (-i)^l is the same as (i)^l in the     ****
*             **** Rayleigh scattering formula            ****
*             *** current function is s or d ****
              if (sd_function) then
                 call Cram_cc_Mul(nb,dcpl_mb(shift),
     >                               dcpl_mb(exi(1)),
     >                               dcpl_mb(prjtmp(1)+(l-1)*npack1))
*             *** current function is p or f ****
              else
                 call Cram_icc_Mul(nb,dcpl_mb(shift),
     >                                dcpl_mb(exi(1)),
     >                                dcpl_mb(prjtmp(1)+(l-1)*npack1))
              end if

*             **** compute 1Xnproj matrix zsw1 = <psi1|prj> ****
              call Cram_cc_izdot(nb,
     >                      psi1,
     >                      dcpl_mb(prjtmp(1)+(l-1)*npack1),
     >                      dcpl_mb(zsw1u(1)+(l-1)))
              call Cram_cc_izdot(nb,
     >                      psi1(shifts+1),
     >                      dcpl_mb(prjtmp(1)+(l-1)*npack1),
     >                      dcpl_mb(zsw1d(1)+(l-1)))
            end do !**l**
            call C3dB_Vector_SumAll((2*nproj),dcpl_mb(zsw1u(1)))
            call C3dB_Vector_SumAll((2*nproj),dcpl_mb(zsw1d(1)))

*           **** zsw2 = Gijl*zsw1 ******
            call Multiply_Kijl_SO_x(1,
     >                         nproj,
     >                         int_mb(nmax(1)+ia-1),
     >                         int_mb(lmax(1)+ia-1),
     >                         int_mb(n_projector(1)
     >                                + (ia-1)*(nmax_max*lmmax_max)),
     >                         int_mb(l_projector(1)
     >                                + (ia-1)*(nmax_max*lmmax_max)),
     >                         int_mb(m_projector(1)
     >                                + (ia-1)*(nmax_max*lmmax_max)),
     >                         dbl_mb(Kijl(1)
     >                         + (ia-1)*nmax_max*nmax_max*(lmax_max+1)),
     >                         dcpl_mb(zsw1u(1)),dcpl_mb(zsw1d(1)),
     >                         dcpl_mb(zsw2u(1)),dcpl_mb(zsw2d(1)))
*           **** do Kleinman-Bylander Multiplication ****
            call dscal(2*nproj,scal,dcpl_mb(zsw2u(1)),1)
            call dscal(2*nproj,scal,dcpl_mb(zsw2d(1)),1)
            call ZGEMM('N','C',npack1,1,nproj,
     >                 mone,
     >                 dcpl_mb(prjtmp(1)), npack1,
     >                 dcpl_mb(zsw2u(1)),   1,
     >                 one,
     >                 psi2, npack1)
            call ZGEMM('N','C',npack1,1,nproj,
     >                 mone,
     >                 dcpl_mb(prjtmp(1)), npack1,
     >                 dcpl_mb(zsw2d(1)),   1,
     >                 one,
     >                 psi2(1+shifts), npack1)

         end if !** nproj>0 **
      end do !** ii***

      value =           MA_pop_stack(zsw2d(2))
      value = value.and.MA_pop_stack(zsw1d(2))
      value = value.and.MA_pop_stack(zsw2u(2))
      value = value.and.MA_pop_stack(zsw1u(2))
      value = value.and.MA_pop_stack(exi(2))
      if (.not.value) 
     > call errquit('cpsp_v_spin_orbit_orb:popping stack',3,MA_ERR)

      call nwpw_timing_end(6)

      return 
      end



