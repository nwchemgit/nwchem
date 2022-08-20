*
* $Id$
*

*     ***********************************
*     *					*
*     *	 	  psp_U_init  		*
*     *					*
*     ***********************************

      subroutine psp_U_init()      
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"
#include "psp.fh"

      integer taskid,MASTER
      parameter (MASTER=0)

      integer nion,i,j,ii
      logical value

*     **** external functions *****
      logical  control_psputerm
      external control_psputerm
      integer  ion_nion,control_pspnuterms
      external ion_nion,control_pspnuterms
 

*     **** setup psputerm structure - used for adding Hubbard U term ****
      pspnuterms = control_pspnuterms()
      psputerm   = control_psputerm().and.(pspnuterms.gt.0)
      if (psputerm) then
         nion  = ion_nion()
         value = BA_alloc_get(mt_int,pspnuterms,'psputerm_l',
     >                psputerm_l(2),psputerm_l(1))

         value = value.and.
     >   BA_alloc_get(mt_dbl,pspnuterms,'psputerm_uscale',
     >                psputerm_uscale(2),psputerm_uscale(1))
         value = value.and.
     >   BA_alloc_get(mt_dbl,pspnuterms,'psputerm_jscale',
     >                psputerm_jscale(2),psputerm_jscale(1))
         value = value.and.
     >   BA_alloc_get(mt_log,nion*pspnuterms,'psputerm_ions',
     >                psputerm_ions(2),psputerm_ions(1))

         value = value.and.
     >   BA_alloc_get(mt_int,pspnuterms,'psputerm_vstart',
     >                psputerm_vstart(2),psputerm_vstart(1))

         if (.not. value)
     >   call errquit('psp_U_init:out of heap memory',0, MA_ERR)

         call control_set_psputerm(nion,pspnuterms,
     >                             int_mb(psputerm_l(1)),
     >                             dbl_mb(psputerm_uscale(1)),
     >                             dbl_mb(psputerm_jscale(1)),
     >                             log_mb(psputerm_ions(1)))

         psputerm_vcount = 0
         do i=1,pspnuterms
            int_mb(psputerm_vstart(1)+i-1) = psputerm_vcount
            psputerm_vcount = psputerm_vcount 
     >                      + (2*int_mb(psputerm_l(1)+i-1)+1)**4
         end do

         j = 0
         do ii=1,nion
         do i=1,pspnuterms
         if (log_mb(psputerm_ions(1)+ii-1+(i-1)*nion)) then
            j = j + 1
         end if
         end do
         end do
         value = BA_alloc_get(mt_int,j,'psputerm_ld_start',
     >             psputerm_ld_start(2),psputerm_ld_start(1))
         if (.not. value)
     >   call errquit('psp_U_init:out of heap memory',1, MA_ERR)


         j = 0
         psputerm_ld_count = 0
         do ii=1,nion
            do i=1,pspnuterms
               if (log_mb(psputerm_ions(1)+ii-1+(i-1)*nion)) then
                  int_mb(psputerm_ld_start(1)+j) = psputerm_ld_count
                  psputerm_ld_count = psputerm_ld_count 
     >                   + 2*(2*int_mb(psputerm_l(1)+i-1)+1)**2
                  j = j + 1
               end if
            end do
         end do

         value = 
     >   BA_alloc_get(mt_dbl,psputerm_vcount,'psputerm_vmmmm',
     >                psputerm_vmmmm(2),psputerm_vmmmm(1))
         value = value.and.
     >   BA_alloc_get(mt_dbl,psputerm_ld_count,'psputerm_ld',
     >                psputerm_ld(2),psputerm_ld(1))
         if (.not. value)
     >   call errquit('psp_U_init:out of heap memory',2, MA_ERR)

         do i=1,pspnuterms
             call psputerm_gen_vmmmm(int_mb(psputerm_l(1)+i-1),
     >                               2*int_mb(psputerm_l(1)+i-1)+1,
     >                               dbl_mb(psputerm_uscale(1)+i-1),
     >                               dbl_mb(psputerm_jscale(1)+i-1),
     >                               dbl_mb(psputerm_vmmmm(1)
     >                               +int_mb(psputerm_vstart(1)+i-1)))
         end do

      end if

      return
      end

*     ***********************************
*     *					*
*     *	 	  psp_U_end  		*
*     *					*
*     ***********************************

      subroutine psp_U_end()      
      implicit none

#include "errquit.fh"
#include "bafdecls.fh"
#include "psp.fh"

*     **** version4 common block ****
      logical value

      value = .true.
      if (psputerm) then
         value = value.and.BA_free_heap(psputerm_l(2))
         value = value.and.BA_free_heap(psputerm_uscale(2))
         value = value.and.BA_free_heap(psputerm_jscale(2))
         value = value.and.BA_free_heap(psputerm_ions(2))
         value = value.and.BA_free_heap(psputerm_vstart(2))
         value = value.and.BA_free_heap(psputerm_vmmmm(2))
         value = value.and.BA_free_heap(psputerm_ld_start(2))
         value = value.and.BA_free_heap(psputerm_ld(2))
      end if

      if (.not.value) 
     >  call errquit('psp_U_end:error freeing heap memory',0,MA_ERR)

      return
      end


*     ***********************************
*     *                                 *
*     *            psp_U_psputerm       *
*     *                                 *
*     ***********************************
      logical function psp_U_psputerm()
      implicit none

#include "psp.fh"

      psp_U_psputerm = psputerm
      return
      end



*     ***********************************
*     *					*
*     *	 	   psp_U_v_nonlocal	*
*     *					*
*     ***********************************

*    This routine computes the Kleinman-Bylander non-local 
* pseudopotential projection.
*
*  Note - This routine was restructured 5-13-2002 to improve
*         parallel efficiency.
*
      subroutine psp_U_v_nonlocal(ispin,ne,psi1,psi2,move,fion,
     >                            fractional,occ,ld_update)
      implicit none
      integer    ispin,ne(2)
      complex*16 psi1(*)
      complex*16 psi2(*)
      logical move
      real*8 fion(3,*)
      logical fractional
      real*8 occ(*)
      logical ld_update

#include "bafdecls.fh"
#include "psp.fh"
#include "errquit.fh"


*     *** local variables ***
      integer G(3),npack1,nion,nu
      integer i,j,ii,ia,l,n,nn
      integer k,shift,l_prj,nproj,Gijl_indx
      real*8  omega,scal,ff(3)
      complex*16 ctmp
      integer exi(2),xtmp(2),sw1(2),sw2(2),sw3(2),asum(2)
      logical value,sd_function
      real*8 vmm(150)
      integer ld_ptr

*     **** external functions ****
      logical  is_sORd
      integer  ion_nion,ion_katm,Pack_G_indx
      integer  psi_data_get_ptr,psi_data_get_chnk
      real*8   lattice_omega
      external is_sORd
      external ion_nion,ion_katm,Pack_G_indx
      external psi_data_get_ptr,psi_data_get_chnk
      external lattice_omega

      call nwpw_timing_start(6) 

*     **** allocate local memory ****
      nion = ion_nion()
      nn = ne(1)+ne(2)
      call Pack_npack(1,npack1)

      value = BA_push_get(mt_dcpl,npack1,'exi', exi(2), exi(1))
      value = value.and.
     >   BA_push_get(mt_dbl,nn*nprj_max*nprj_mult,'sw1',sw1(2),sw1(1))
      value = value.and.
     >   BA_push_get(mt_dbl,nn*nprj_max*nprj_mult,'sw2',sw2(2),sw2(1))
      value = value.and.
     >   BA_push_get(mt_dbl,nn*nprj_max*nprj_mult,'sw3',sw3(2),sw3(1))
      if (.not.value) 
     >  call errquit('v_nonlocal: out of stack',0, MA_ERR)

      if (move) then
       value = value.and.BA_push_get(mt_dbl,npack1,
     >                               'xtmp',xtmp(2),xtmp(1))
       value = value.and.BA_push_get(mt_dbl,3*nn,'asum',asum(2),asum(1))
       if (.not. value) 
     >  call errquit('v_nonlocal:out of stack memory',1,MA_ERR)

       G(1)  = Pack_G_indx(1,1)
       G(2)  = Pack_G_indx(1,2)
       G(3)  = Pack_G_indx(1,3)
      end if

      omega = lattice_omega()
      scal = 1.0d0/omega

      j = 0
      psputerm_edftu = 0.0d0
      psputerm_pdftu = 0.0d0
      do ii=1,nion
        value = .false.
        do nu=1,pspnuterms
           value=value.or.log_mb(psputerm_ions(1)+ii-1+(nu-1)*nion)
        end do
        if (value) then
           ia=ion_katm(ii)
           nproj = int_mb(nprj(1)+ia-1)

*          **** structure factor and local pseudopotential ****
           call strfac_pack(1,ii,dcpl_mb(exi(1)))

*          **** generate sw1's and projectors ****
           do l=1,nproj
              shift = psi_data_get_ptr(int_mb(vnl(1)+ia-1),l)
              l_prj = int_mb(l_projector(1)+(l-1) 
     >                                  + (ia-1)*(nmax_max*lmmax_max))

#ifdef GCC4
              k = iand(l_prj,1)
#else
              k = and(l_prj,1)
#endif
              sd_function = (k.eq.0)

*             *** current function is s or d ****
              if (sd_function) then
                 call Pack_tc_Mul(1,dbl_mb(shift),
     >                              dcpl_mb(exi(1)),
     >                              dcpl_mb(prjtmp(1)+(l-1)*npack1))

*             *** current function is p or f ****
              else
                 call Pack_tc_iMul(1,dbl_mb(shift),
     >                              dcpl_mb(exi(1)),
     >                              dcpl_mb(prjtmp(1)+(l-1)*npack1))

              end if
              call Pack_cc_indot(1,nn,
     >                      psi1,
     >                      dcpl_mb(prjtmp(1)+(l-1)*npack1),
     >                      dbl_mb(sw1(1)+(l-1)*nn))
           end do
           call D3dB_Vector_SumAll((nn*nproj*nprj_mult),dbl_mb(sw1(1)))

*          **** sw2 = Gijl*sw1 ******
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

           call dcopy(nn*nproj,0.0d0,0,dbl_mb(sw3(1)),1)
           do nu=1,pspnuterms
           if (log_mb(psputerm_ions(1)+ii-1+(nu-1)*nion)) then
              ld_ptr =psputerm_ld(1)+int_mb(psputerm_ld_start(1)+j)
              j = j + 1
              if (ld_update) then
                 call psputerm_gen_l_density(ispin,ne,nproj,
     >                        int_mb(l_projector(1)
     >                               +(ia-1)*(nmax_max*lmmax_max)),
     >                        int_mb(m_projector(1)
     >                              +(ia-1)*(nmax_max*lmmax_max)),
     >                        dbl_mb(sw2(1)),
     >                        int_mb(psputerm_l(1)+nu-1),
     >                        2*int_mb(psputerm_l(1)+nu-1)+1,
     >                        dbl_mb(ld_ptr))
              end if
              call psputerm_gen_vmm(ispin,
     >                     int_mb(psputerm_l(1)+nu-1),
     >                     2*int_mb(psputerm_l(1)+nu-1)+1,
     >                     dbl_mb(psputerm_uscale(1)+nu-1),
     >                     dbl_mb(psputerm_jscale(1)+nu-1),
     >                     dbl_mb(psputerm_vmmmm(1)
     >                           +int_mb(psputerm_vstart(1)+nu-1)),
     >                     dbl_mb(ld_ptr),vmm,
     >                     psputerm_edftu,psputerm_pdftu)
              call psputerm_add_upotential(ispin,ne,nproj,
     >                     int_mb(l_projector(1)
     >                            +(ia-1)*(nmax_max*lmmax_max)),
     >                     int_mb(m_projector(1)
     >                           +(ia-1)*(nmax_max*lmmax_max)),
     >                     int_mb(psputerm_l(1)+nu-1),
     >                     2*int_mb(psputerm_l(1)+nu-1)+1,
     >                     vmm,dbl_mb(sw2(1)),dbl_mb(sw3(1)))
           end if
           end do

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
     >                         dbl_mb(sw3(1)),
     >                         dbl_mb(sw2(1)))
        
*          **** do Kleinman-Bylander Multiplication ****
           call dscal(nn*int_mb(nprj(1)+ia-1)*nprj_mult,
     >                scal,dbl_mb(sw2(1)),1)

          call DGEMM('N','T',2*npack1,nn,int_mb(nprj(1)+ia-1)*nprj_mult,
     >               (1.0d0),
     >               dcpl_mb(prjtmp(1)), 2*npack1,
     >               dbl_mb(sw2(1)),     nn,
     >               (1.0d0),
     >               psi2,               2*npack1)


        if (move) then
        do l=1,nproj
             do n=1,nn
                if (ispin.eq.1) 
     >            dbl_mb(sw2(1)+n-1+(l-1)*nn)
     >            =2.0d0*dbl_mb(sw2(1)+n-1+(l-1)*nn) !// change

#ifndef CRAY
!DIR$ ivdep
#endif
                do i=1,npack1
                   ctmp = psi1(i+(n-1)*npack1)
     >                  *dconjg(dcpl_mb(prjtmp(1)+(l-1)*npack1 + i-1))
                   dbl_mb(xtmp(1)+i-1) = dimag(ctmp)
                end do
                call Pack_tt_idot(1,dbl_mb(G(1)),dbl_mb(xtmp(1)),
     >                            dbl_mb(asum(1)+3*(n-1)))
                call Pack_tt_idot(1,dbl_mb(G(2)),dbl_mb(xtmp(1)),
     >                            dbl_mb(asum(1)+1+3*(n-1)))
                call Pack_tt_idot(1,dbl_mb(G(3)),dbl_mb(xtmp(1)),
     >                            dbl_mb(asum(1)+2+3*(n-1)))

             end do

             call D3dB_Vector_SumAll(3*(nn),dbl_mb(asum(1)))

             !**** fractional weighting ****
             if (fractional) then
             do n=1,nn
              call Dneall_qton(n,i)
              dbl_mb(asum(1)  +3*(n-1))=dbl_mb(asum(1)  +3*(n-1))*occ(i)
              dbl_mb(asum(1)+1+3*(n-1))=dbl_mb(asum(1)+1+3*(n-1))*occ(i)
              dbl_mb(asum(1)+2+3*(n-1))=dbl_mb(asum(1)+2+3*(n-1))*occ(i)
             end do
             end if
  
             ff(1) = 0.0d0
             ff(2) = 0.0d0
             ff(3) = 0.0d0
             do n=1,nn
                ff(1) = ff(1) + 2.0d0*dbl_mb(sw2(1)+n-1+(l-1)*nn) !// change
     >                               *dbl_mb(asum(1)+  3*(n-1))
                ff(2) = ff(2) + 2.0d0*dbl_mb(sw2(1)+n-1+(l-1)*nn) !// change
     >                               *dbl_mb(asum(1)+1+3*(n-1))
                ff(3) = ff(3) + 2.0d0*dbl_mb(sw2(1)+n-1+(l-1)*nn) !// change
     >                               *dbl_mb(asum(1)+2+3*(n-1))
             end do
             call D1dB_Vector_SumAll(3,ff)
             fion(1,ii) = fion(1,ii)  - ff(1)
             fion(2,ii) = fion(2,ii)  - ff(2)
             fion(3,ii) = fion(3,ii)  - ff(3)
       
        end do !** l **
        end if !** move **


        end if !** if nproj>0 **
      end do !** ii **

      value = .true.
      if (move) then
      value = value.and.BA_pop_stack(asum(2))
      value = value.and.BA_pop_stack(xtmp(2))
      end if
      value = value.and.BA_pop_stack(sw3(2))
      value = value.and.BA_pop_stack(sw2(2))
      value = value.and.BA_pop_stack(sw1(2))
      value = value.and.BA_pop_stack(exi(2))
      if (.not.value) call errquit('v_nonlocal: popping stack',3,
     &       MA_ERR)
      call nwpw_timing_end(6)

      return 
      end


*     ***********************************
*     *					*
*     *	      f_psp_U_v_nonlocal	*
*     *					*
*     ***********************************

*    This routine computes the Kleinman-Bylander non-local 
* pseudopotential projection.
*
*  Note - This routine was restructured 5-13-2002 to improve
*         parallel efficiency.
*
      subroutine f_psp_U_v_nonlocal(ispin,ne,psi1,fion,
     >                            fractional,occ,ld_update)
      implicit none
      integer    ispin,ne(2)
      complex*16 psi1(*)
      real*8 fion(3,*)
      logical fractional
      real*8 occ(*)
      logical ld_update

#include "bafdecls.fh"
#include "psp.fh"
#include "errquit.fh"

*     *** local variables ***
      integer G(3),npack1,nion,nu
      integer i,j,ii,ia,l,n,nn
      integer k,shift,l_prj,nproj,Gijl_indx
      real*8  omega,scal,ff(3)
      complex*16 ctmp
      integer exi(2),xtmp(2),sw1(2),sw2(2),sw3(2),asum(2)
      logical value,sd_function
      real*8 vmm(150)
      integer ld_ptr

*     **** external functions ****
      logical  is_sORd
      integer  ion_nion,ion_katm,Pack_G_indx
      integer  psi_data_get_ptr,psi_data_get_chnk
      real*8   lattice_omega
      external is_sORd
      external ion_nion,ion_katm,Pack_G_indx
      external psi_data_get_ptr,psi_data_get_chnk
      external lattice_omega

      call nwpw_timing_start(6) 

*     **** allocate local memory ****
      nion = ion_nion()
      nn = ne(1)+ne(2)
      call Pack_npack(1,npack1)

      value = BA_push_get(mt_dcpl,npack1,'exi', exi(2), exi(1))
      value = value.and.
     >   BA_push_get(mt_dbl,nn*nprj_max*nprj_mult,'sw1',sw1(2),sw1(1))
      value = value.and.
     >   BA_push_get(mt_dbl,nn*nprj_max*nprj_mult,'sw2',sw2(2),sw2(1))
      value = value.and.
     >   BA_push_get(mt_dbl,nn*nprj_max*nprj_mult,'sw3',sw3(2),sw3(1))
      if (.not.value) 
     >  call errquit('v_nonlocal: out of stack',0, MA_ERR)

      value = value.and.BA_push_get(mt_dbl,npack1,
     >                               'xtmp',xtmp(2),xtmp(1))
      value = value.and.BA_push_get(mt_dbl,3*nn,'asum',asum(2),asum(1))
      if (.not. value) 
     > call errquit('v_nonlocal:out of stack memory',1,MA_ERR)

      G(1)  = Pack_G_indx(1,1)
      G(2)  = Pack_G_indx(1,2)
      G(3)  = Pack_G_indx(1,3)

      omega = lattice_omega()
      scal = 1.0d0/omega

      j = 0
      psputerm_edftu = 0.0d0
      psputerm_pdftu = 0.0d0
      do ii=1,nion
        value = .false.
        do nu=1,pspnuterms
           value=value.or.log_mb(psputerm_ions(1)+ii-1+(nu-1)*nion)
        end do
        if (value) then
           ia=ion_katm(ii)
           nproj = int_mb(nprj(1)+ia-1)

*          **** structure factor and local pseudopotential ****
           call strfac_pack(1,ii,dcpl_mb(exi(1)))

*          **** generate sw1's and projectors ****
           do l=1,nproj
              shift = psi_data_get_ptr(int_mb(vnl(1)+ia-1),l)
              l_prj = int_mb(l_projector(1)+(l-1) 
     >                                  + (ia-1)*(nmax_max*lmmax_max))

#ifdef GCC4
              k = iand(l_prj,1)
#else
              k = and(l_prj,1)
#endif
              sd_function = (k.eq.0)

*             *** current function is s or d ****
              if (sd_function) then
                 call Pack_tc_Mul(1,dbl_mb(shift),
     >                              dcpl_mb(exi(1)),
     >                              dcpl_mb(prjtmp(1)+(l-1)*npack1))

*             *** current function is p or f ****
              else
                 call Pack_tc_iMul(1,dbl_mb(shift),
     >                              dcpl_mb(exi(1)),
     >                              dcpl_mb(prjtmp(1)+(l-1)*npack1))

              end if
              call Pack_cc_indot(1,nn,
     >                      psi1,
     >                      dcpl_mb(prjtmp(1)+(l-1)*npack1),
     >                      dbl_mb(sw1(1)+(l-1)*nn))
           end do
           call D3dB_Vector_SumAll((nn*nproj*nprj_mult),dbl_mb(sw1(1)))

*          **** sw2 = Gijl*sw1 ******
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

           call dcopy(nn*nproj,0.0d0,0,dbl_mb(sw3(1)),1)
           do nu=1,pspnuterms
           if (log_mb(psputerm_ions(1)+ii-1+(nu-1)*nion)) then
              ld_ptr =psputerm_ld(1)+int_mb(psputerm_ld_start(1)+j)
              j = j + 1
              if (ld_update) then
                 call psputerm_gen_l_density(ispin,ne,nproj,
     >                        int_mb(l_projector(1)
     >                               +(ia-1)*(nmax_max*lmmax_max)),
     >                        int_mb(m_projector(1)
     >                              +(ia-1)*(nmax_max*lmmax_max)),
     >                        dbl_mb(sw2(1)),
     >                        int_mb(psputerm_l(1)+nu-1),
     >                        2*int_mb(psputerm_l(1)+nu-1)+1,
     >                        dbl_mb(ld_ptr))
              end if
              call psputerm_gen_vmm(ispin,
     >                     int_mb(psputerm_l(1)+nu-1),
     >                     2*int_mb(psputerm_l(1)+nu-1)+1,
     >                     dbl_mb(psputerm_uscale(1)+nu-1),
     >                     dbl_mb(psputerm_jscale(1)+nu-1),
     >                     dbl_mb(psputerm_vmmmm(1)
     >                           +int_mb(psputerm_vstart(1)+nu-1)),
     >                     dbl_mb(ld_ptr),vmm,
     >                     psputerm_edftu,psputerm_pdftu)
              call psputerm_add_upotential(ispin,ne,nproj,
     >                     int_mb(l_projector(1)
     >                            +(ia-1)*(nmax_max*lmmax_max)),
     >                     int_mb(m_projector(1)
     >                           +(ia-1)*(nmax_max*lmmax_max)),
     >                     int_mb(psputerm_l(1)+nu-1),
     >                     2*int_mb(psputerm_l(1)+nu-1)+1,
     >                     vmm,dbl_mb(sw2(1)),dbl_mb(sw3(1)))
           end if
           end do

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
     >                         dbl_mb(sw3(1)),
     >                         dbl_mb(sw2(1)))
        
*          **** do Kleinman-Bylander Multiplication ****
           call dscal(nn*int_mb(nprj(1)+ia-1)*nprj_mult,
     >                scal,dbl_mb(sw2(1)),1)


        do l=1,nproj
             do n=1,nn
                if (ispin.eq.1) 
     >            dbl_mb(sw2(1)+n-1+(l-1)*nn)
     >            =2.0d0*dbl_mb(sw2(1)+n-1+(l-1)*nn) !// change

#ifndef CRAY
!DIR$ ivdep
#endif
                do i=1,npack1
                   ctmp = psi1(i+(n-1)*npack1)
     >                  *dconjg(dcpl_mb(prjtmp(1)+(l-1)*npack1 + i-1))
                   dbl_mb(xtmp(1)+i-1) = dimag(ctmp)
                end do
                call Pack_tt_idot(1,dbl_mb(G(1)),dbl_mb(xtmp(1)),
     >                            dbl_mb(asum(1)+3*(n-1)))
                call Pack_tt_idot(1,dbl_mb(G(2)),dbl_mb(xtmp(1)),
     >                            dbl_mb(asum(1)+1+3*(n-1)))
                call Pack_tt_idot(1,dbl_mb(G(3)),dbl_mb(xtmp(1)),
     >                            dbl_mb(asum(1)+2+3*(n-1)))

             end do

             call D3dB_Vector_SumAll(3*(nn),dbl_mb(asum(1)))

             !**** fractional weighting ****
             if (fractional) then
             do n=1,nn
              call Dneall_qton(n,i)
              dbl_mb(asum(1)  +3*(n-1))=dbl_mb(asum(1)  +3*(n-1))*occ(i)
              dbl_mb(asum(1)+1+3*(n-1))=dbl_mb(asum(1)+1+3*(n-1))*occ(i)
              dbl_mb(asum(1)+2+3*(n-1))=dbl_mb(asum(1)+2+3*(n-1))*occ(i)
             end do
             end if
  
             ff(1) = 0.0d0
             ff(2) = 0.0d0
             ff(3) = 0.0d0
             do n=1,nn
                ff(1) = ff(1) + 2.0d0*dbl_mb(sw2(1)+n-1+(l-1)*nn) !// change
     >                               *dbl_mb(asum(1)+  3*(n-1))
                ff(2) = ff(2) + 2.0d0*dbl_mb(sw2(1)+n-1+(l-1)*nn) !// change
     >                               *dbl_mb(asum(1)+1+3*(n-1))
                ff(3) = ff(3) + 2.0d0*dbl_mb(sw2(1)+n-1+(l-1)*nn) !// change
     >                               *dbl_mb(asum(1)+2+3*(n-1))
             end do
             call D1dB_Vector_SumAll(3,ff)
             fion(1,ii) = fion(1,ii)  - ff(1)
             fion(2,ii) = fion(2,ii)  - ff(2)
             fion(3,ii) = fion(3,ii)  - ff(3)
       
        end do !** l **


        end if !** if nproj>0 **
      end do !** ii **

      value = .true.
      value = value.and.BA_pop_stack(asum(2))
      value = value.and.BA_pop_stack(xtmp(2))
      value = value.and.BA_pop_stack(sw3(2))
      value = value.and.BA_pop_stack(sw2(2))
      value = value.and.BA_pop_stack(sw1(2))
      value = value.and.BA_pop_stack(exi(2))
      if (.not.value) call errquit('v_nonlocal: popping stack',3,
     &       MA_ERR)
      call nwpw_timing_end(6)

      return 
      end





*     ****************************************************
*     *                                                  *
*     *                psputerm_gen_vmmmm                *
*     *                                                  *
*     ****************************************************
      subroutine psputerm_gen_vmmmm(l,lm,U,J,vmmmm)
      implicit none
      integer l,lm
      real*8 U,J
      real*8 vmmmm(lm,lm,lm,lm)

*     **** local variables ****
      integer m1,m2,m3,m4,k,q,kk
      real*8 F(4),fourpi,aterm

*     **** external functions ****
      real*8   nwpw_gaunt
      external nwpw_gaunt

      fourpi = 16.0d0*datan(1.0d0)

      if (l.eq.0) then
        F(1) = U
        F(2) = 0.0d0
        F(3) = 0.0d0
        F(4) = 0.0d0
      end if
      if (l.eq.1) then
        F(1) = U
        F(2) = 5*J
        F(3) = 0.0d0
        F(4) = 0.0d0
      end if
      if (l.eq.2) then
        F(1) = U
        F(2) = 112.0d0/13.0d0*J
        F(3) =   70.d0/13.0d0*J
        F(4) = 0.0d0
      end if
      if (l.eq.3) then
        F(1) = U
        F(2) = 4125.0d0/346.0d0*J
        F(3) = 5511.0d0/692.0d0*J
        F(4) = 8151.0d0/13840d0*J
      end if

      call nwpw_gaunt_init(.false.,3)
      do m1=-l,l
      do m3=-l,l
      do m2=-l,l
      do m4=-l,l
         vmmmm(m1+l+1,m2+l+1,m3+l+1,m4+l+1) = 0.0d0

         do kk=1,4
            k = 2*(kk-1)
            aterm = 0.0d0
            do q=-k,k
               aterm = aterm + nwpw_gaunt(.true.,k,q,l,m1,l,m2)
     >                        *nwpw_gaunt(.true.,k,q,l,m4,l,m3)
            end do
            aterm = aterm*fourpi/dble(2*k+1)

            vmmmm(m1+l+1,m2+l+1,m3+l+1,m4+l+1) 
     >         = vmmmm(m1+l+1,m2+l+1,m3+l+1,m4+l+1)
     >         +  F(kk)*aterm
         end do
      end do
      end do
      end do
      end do
      call nwpw_gaunt_end()

      return
      end 

*     ****************************************************
*     *                                                  *
*     *                psputerm_gen_vmm                  *
*     *                                                  *
*     ****************************************************
      subroutine psputerm_gen_vmm(ispin,l,lm,U,J,vmmmm,l_density,vmm,
     >                            edftu,pdftu)
      implicit none
      integer ispin,l,lm
      real*8 U,J
      real*8 vmmmm(lm,lm,lm,lm)
      real*8 l_density(lm,lm,ispin)
      real*8 vmm(lm,lm,ispin)
      real*8 edftu,pdftu

*     *** local variables ****
      integer ms,m1,m2,m3,m4,n1,n2,n3,n4
      real*8  density_total,density(2),edftu1,pdftu1

      do ms=1,ispin
         density(ms) = 0.0d0
         do n1=1,lm
            density(ms) = density(ms) + l_density(n1,n1,ms)
         end do
      end do
      density_total = density(1)+density(ispin)
  
      edftu1 = 0.0d0
      pdftu1 = 0.0d0

      do ms=1,ispin
         do m1=-l,l
         do m2=-l,l
            n1=m1+l+1
            n2=m2+l+1
            vmm(n1,n2,ms)=0.0d0
            do m3=-l,l
            do m4=-l,l
               n3 = m3+l+1
               n4 = m4+l+1
            
               vmm(n1,n2,ms) = vmm(n1,n2,ms)
     >                       + (vmmmm(n1,n2,n3,n4)-vmmmm(n1,n4,n3,n2))
     >                        *l_density(n3,n4,ms)
               vmm(n1,n2,ms) = vmm(n1,n2,ms)
     >                       + vmmmm(n1,n2,n3,n4)
     >                        *l_density(n3,n4,ispin+1-ms)
            end do
            end do
            edftu1 = edftu1 + 0.5d0*vmm(n1,n2,ms)*l_density(n1,n2,ms)
            if(m1==m2)  then
               vmm(n1,n2,ms)=vmm(n1,n2,ms)-U*(density_total-0.5d0)
               vmm(n1,n2,ms)=vmm(n1,n2,ms)+J*(density(ms)  -0.5d0)
            end if
            pdftu1 = pdftu1 + vmm(n1,n2,ms)*l_density(n1,n2,ms)
         end do
         end do
         edftu1 = edftu1 + 0.5d0*J*density(ms)*(density(ms)-1.0d0)
      end do
      edftu1 = edftu1 - 0.5d0*U*density_total*(density_total-1.0d0)

      if (ispin.eq.1) then
         edftu1 = edftu1 + edftu1
         pdftu1 = pdftu1 + pdftu1
      end if
      edftu = edftu + edftu1
      pdftu = pdftu + pdftu1


      return 
      end 
      

*     ****************************************************
*     *                                                  *
*     *                psputerm_gen_l_density            *
*     *                                                  *
*     ****************************************************
      subroutine psputerm_gen_l_density(ispin,ne,nprj,l_prj,m_prj,sw2,
     >                                  l,lm,l_density)
      implicit none
      integer ispin,ne(2),nprj
      integer l_prj(nprj),m_prj(nprj)
      real*8  sw2(ne(1)+ne(2),nprj)
      integer l,lm
      real*8  l_density(lm,lm,ispin)

*     **** local variables ****
      real*8  omega,trace
      integer a,b,c,s,k,ma,mb,ms,n,n1(2),n2(2)

*     **** external functions ****
      real*8   lattice_omega
      external lattice_omega

      omega = lattice_omega()
      n1(1) = 1
      n1(2) = ne(1)+1
      n2(1) = ne(1)
      n2(2) = ne(1)+ne(2)

      call dcopy(lm*lm*ispin,0.0d0,0,l_density,1)
      do a=1,nprj
      do b=1,nprj
      if ((l_prj(a).eq.l).and.(l_prj(b).eq.l)) then
         mb = m_prj(b)+l+1
         ma = m_prj(a)+l+1
         do ms=1,ispin
         do n=n1(ms),n2(ms)
            l_density(ma,mb,ms) = l_density(ma,mb,ms)
     >                          + sw2(n,a)*sw2(n,b)/omega
         end do
         end do
      end if
      end do
      end do
      call D1dB_Vector_SumAll(lm*lm*ispin,l_density)

      return
      end 

*     ****************************************************
*     *                                                  *
*     *                psputerm_add_upotential           *
*     *                                                  *
*     ****************************************************
      subroutine psputerm_add_upotential(ispin,ne,nprj,l_prj,m_prj,
     >                                   l,lm,vmm,sw2,sw3)
      implicit none
      integer ispin,ne(2),nprj
      integer l_prj(nprj),m_prj(nprj)
      integer l,lm
      real*8  vmm(lm,lm,ispin)
      real*8  sw2(ne(1)+ne(2),nprj)
      real*8  sw3(ne(1)+ne(2),nprj)

*     **** local variables ****
      integer a,b,ms,i,ma,mb,n1(2),n2(2)
      real*8 omega

*     **** external functions ****
      real*8   lattice_omega
      external lattice_omega
  
      omega = lattice_omega()
      n1(1) = 1
      n1(2) = ne(1)+1
      n2(1) = ne(1)
      n2(2) = ne(1)+ne(2)

      do ms =1,ispin
      do i=n1(ms),n2(ms)
         do a=1,nprj
         do b=1,nprj
            if ((l_prj(a).eq.l).and.(l_prj(b).eq.l)) then
               mb = m_prj(b)+l+1
               ma = m_prj(a)+l+1
               sw3(i,a) = sw3(i,a)
     >                  - sw2(i,b)*vmm(ma,mb,ms)
            end if
         end do
         end do
      end do
      end do

      return 
      end

*     ****************************************************
*     *                                                  *
*     *                psp_U_psputerm_energy             *
*     *                                                  *
*     ****************************************************
      subroutine psp_U_psputerm_energy(edftu,pdftu)
      implicit none
      real*8 edftu,pdftu

#include "psp.fh"
      
      edftu = psputerm_edftu
      pdftu = psputerm_pdftu
      return
      end 
