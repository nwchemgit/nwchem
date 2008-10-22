*
* $Id: paw.F,v 1.3 2008-10-22 23:56:45 bylaska Exp $
*

*     ***********************************
*     *					*
*     *	 	   psp_overlap_S	*
*     *					*
*     ***********************************

*    This routine computes the paw overlap S operator to psi1
*      psi2 = S*psi1
*
      subroutine psp_overlap_S(ispin,ne,psi1,psi2)
      implicit none
      integer    ispin,ne(2)
      complex*16 psi1(*)
      complex*16 psi2(*)

#include "mafdecls.fh"
#include "psp.fh"
#include "errquit.fh"

*     *** local variables ***
      integer npack1,nion
      integer ii,ia,l,nn
      integer k,shift,l_prj,nproj,Gijl_indx
      real*8  omega,scal,scalsqr
      integer exi(2),sw1(2),sw2(2)
      logical value,sd_function

*     **** external functions ****
      integer  ion_nion,ion_katm
      integer  psi_data_get_ptr
      real*8   lattice_omega
      external ion_nion,ion_katm
      external psi_data_get_ptr
      external lattice_omega

      call nwpw_timing_start(6) 

*     **** allocate local memory ****
      nion = ion_nion()
      nn = ne(1)+ne(2)
      call Pack_npack(1,npack1)

      value = MA_push_get(mt_dcpl,npack1,'exi', exi(2), exi(1))
      value = value.and.
     >        MA_push_get(mt_dbl,nn*nprj_max,'sw1',sw1(2),sw1(1))
      value = value.and.
     >        MA_push_get(mt_dbl,nn*nprj_max,'sw2',sw2(2),sw2(1))
      if (.not.value) 
     >  call errquit('psp_overlap_S: out of stack',0,MA_ERR)

      omega = lattice_omega()
      scal = 1.0d0/(omega)
      scalsqr = scal*scal

      call dcopy(2*npack1*nn,psi1,1,psi2,1)
      do ii=1,nion
        ia=ion_katm(ii)

        nproj   = int_mb(nprj(1)+ia-1)

        if ((int_mb(psp_type(1)+ia-1).eq.4).and.(nproj.gt.0)) then

*       **** structure factor and local pseudopotential ****
        call strfac_pack(1,ii,dcpl_mb(exi(1)))

*       **** generate sw1's and projectors ****
        do l=1,nproj

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


*          **** phase factor does not matter therefore ****
*          **** (-i)^l is the same as (i)^l in the     ****
*          **** Rayleigh scattering formula            ****

*          *** current function is s or d ****
           if (sd_function) then
              call Pack_tc_Mul(1,dbl_mb(shift),
     >                           dcpl_mb(exi(1)),
     >                           dcpl_mb(prjtmp(1)+(l-1)*npack1))
*          *** current function is p or f ****
           else
              call Pack_tc_iMul(1,dbl_mb(shift),
     >                           dcpl_mb(exi(1)),
     >                           dcpl_mb(prjtmp(1)+(l-1)*npack1))
           end if
           call Pack_cc_indot(1,nn,
     >                      psi1,
     >                      dcpl_mb(prjtmp(1)+(l-1)*npack1),
     >                      dbl_mb(sw1(1)+(l-1)*nn))

        end do
        call D3dB_Vector_SumAll((nn*nproj),dbl_mb(sw1(1)))

*       **** sw2 = Sijl*sw1 ******
        Gijl_indx = psi_data_get_ptr(int_mb(Gijl(1)+ia-1),2)
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
        
*       **** do Kleinman-Bylander Multiplication ****
        !scal = 1.0d0/(omega)
        call dscal(nn*int_mb(nprj(1)+ia-1),scal,dbl_mb(sw2(1)),1)
        call DGEMM('N','T',2*npack1,nn,int_mb(nprj(1)+ia-1),
     >             (1.0d0),
     >             dcpl_mb(prjtmp(1)), 2*npack1,
     >             dbl_mb(sw2(1)),     nn,
     >             (1.0d0),
     >             psi2,               2*npack1)


        end if !** nproj>0 **
      end do !** ii **

      value =           MA_pop_stack(sw2(2))
      value = value.and.MA_pop_stack(sw1(2))
      value = value.and.MA_pop_stack(exi(2))
      if (.not.value) call errquit('psp_overlap_S: popping stack',3,
     &       MA_ERR)
      call nwpw_timing_end(6)
      return 
      end

*     ***********************************
*     *					*
*     *	 	   psp_overlap		*
*     *					*
*     ***********************************

*    This routine computes the paw overlap S operator to psi1
*      psi2 = S*psi1
*
      subroutine psp_overlap(ispin,ne,psi1,S)
      implicit none
      integer    ispin,ne(2)
      complex*16 psi1(*)
      real*8     S(*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "psp.fh"

*     *** local variables ***
      integer npack1,nion
      integer ii,ia,l,nn,ms,shifts,shiftsw
      integer k,shift,l_prj,nproj,Gijl_indx
      real*8  omega,scal,scalsqr
      integer exi(2),sw1(2),sw2(2)
      logical value,sd_function

*     **** external functions ****
      integer  ion_nion,ion_katm
      integer  psi_data_get_ptr
      real*8   lattice_omega
      external ion_nion,ion_katm
      external psi_data_get_ptr
      external lattice_omega

      call nwpw_timing_start(6) 

*     **** S = transpose(psi)*psi ****
      call Dneall_ffm_sym_Multiply(0,psi1,psi1,npack1,S)

*     **** allocate local memory ****
      nion = ion_nion()
      nn = ne(1)+ne(2)
      call Pack_npack(1,npack1)

      value = MA_push_get(mt_dcpl,npack1,'exi', exi(2), exi(1))
      value = value.and.
     >        MA_push_get(mt_dbl,nn*nprj_max,'sw1',sw1(2),sw1(1))
      value = value.and.
     >        MA_push_get(mt_dbl,nn*nprj_max,'sw2',sw2(2),sw2(1))
      if (.not.value) 
     >  call errquit('psp_overlap: out of stack',0,MA_ERR)

      omega = lattice_omega()
      scal    = 1.0d0/(omega)
      scalsqr = scal*scal

      do ii=1,nion

        ia    = ion_katm(ii)
        nproj = int_mb(nprj(1)+ia-1)

        if ((int_mb(psp_type(1)+ia-1).eq.4).and.(nproj.gt.0)) then

*       **** structure factor and local pseudopotential ****
        call strfac_pack(1,ii,dcpl_mb(exi(1)))

*       **** generate sw1's and projectors ****
        do l=1,nproj

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


*          **** phase factor does not matter therefore ****
*          **** (-i)^l is the same as (i)^l in the     ****
*          **** Rayleigh scattering formula            ****

*          *** current function is s or d ****
           if (sd_function) then
              call Pack_tc_Mul(1,dbl_mb(shift),
     >                           dcpl_mb(exi(1)),
     >                           dcpl_mb(prjtmp(1)+(l-1)*npack1))
*          *** current function is p or f ****
           else
              call Pack_tc_iMul(1,dbl_mb(shift),
     >                           dcpl_mb(exi(1)),
     >                           dcpl_mb(prjtmp(1)+(l-1)*npack1))
           end if
           call Pack_cc_indot(1,nn,
     >                      psi1,
     >                      dcpl_mb(prjtmp(1)+(l-1)*npack1),
     >                      dbl_mb(sw1(1)+(l-1)*nn))
        end do
        call D3dB_Vector_SumAll((nn*nproj),dbl_mb(sw1(1)))

*       **** sw2 = Sijl*sw1 ******
        Gijl_indx = psi_data_get_ptr(int_mb(Gijl(1)+ia-1),2)
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

*       *** routine needs to be parallelized over orbitals ****
*       **** S = S + sw1*transpose(sw2) ****
        do ms=1,ispin
          shifts  = 1+(ms-1)*ne(1)*ne(1)
          shiftsw =   (ms-1)*ne(1)
         write(*,*) "into DGEMM ",ms
          call DGEMM('N','T',
     >              ne(ms),ne(ms),int_mb(nprj(1)+ia-1),
     >              (scal),
     >              dbl_mb(sw1(1)+shiftsw), nn,
     >              dbl_mb(sw2(1)+shiftsw), nn,
     >              (1.0d0),
     >              S(shifts), ne(ms))
        end do
        

        end if !** nproj>0 **
      end do !** ii **

      value =           MA_pop_stack(sw2(2))
      value = value.and.MA_pop_stack(sw1(2))
      value = value.and.MA_pop_stack(exi(2))
      if (.not.value) call errquit('psp_overlap: popping stack',3,
     &       MA_ERR)
      call nwpw_timing_end(6)
      return 
      end



*     ***********************************
*     *					*
*     *	 	   psp_overlap_orb	*
*     *					*
*     ***********************************

*    This routine computes the paw overlap S operator to psi1
*      psi2 = S*psi1
*
      subroutine psp_overlap_orb(n,psi1,S)
      implicit none
      integer    n
      complex*16 psi1(*)
      real*8     S(*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "psp.fh"

*     *** local variables ***
      integer npack1,nion
      integer ii,ia,l
      integer k,shift,l_prj,nproj,Gijl_indx
      real*8  omega,scal,scalsqr
      integer exi(2),sw1(2),sw2(2)
      logical value,sd_function

*     **** external functions ****
      integer  ion_nion,ion_katm
      integer  psi_data_get_ptr
      real*8   lattice_omega
      external ion_nion,ion_katm
      external psi_data_get_ptr
      external lattice_omega

      call nwpw_timing_start(6) 

*     **** S = transpose(psi)*psi ****
      call Pack_ccm_sym_dot(1,n,psi1,psi1,S)


*     **** allocate local memory ****
      nion = ion_nion()
      call Pack_npack(1,npack1)

      value = MA_push_get(mt_dcpl,npack1,'exi', exi(2), exi(1))
      value = value.and.
     >        MA_push_get(mt_dbl,n*nprj_max,'sw1',sw1(2),sw1(1))
      value = value.and.
     >        MA_push_get(mt_dbl,n*nprj_max,'sw2',sw2(2),sw2(1))
      if (.not.value) 
     >  call errquit('psp_overlap_orb: out of stack',0,MA_ERR)

      omega = lattice_omega()
      scal    = 1.0d0/(omega)
      scalsqr = scal*scal


      do ii=1,nion
        ia    = ion_katm(ii)
        nproj = int_mb(nprj(1)+ia-1)

        if ((int_mb(psp_type(1)+ia-1).eq.4).and.(nproj.gt.0)) then

*       **** structure factor and local pseudopotential ****
        call strfac_pack(1,ii,dcpl_mb(exi(1)))

*       **** generate sw1's and projectors ****
        do l=1,nproj

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


*          **** phase factor does not matter therefore ****
*          **** (-i)^l is the same as (i)^l in the     ****
*          **** Rayleigh scattering formula            ****

*          *** current function is s or d ****
           if (sd_function) then
              call Pack_tc_Mul(1,dbl_mb(shift),
     >                           dcpl_mb(exi(1)),
     >                           dcpl_mb(prjtmp(1)+(l-1)*npack1))
*          *** current function is p or f ****
           else
              call Pack_tc_iMul(1,dbl_mb(shift),
     >                           dcpl_mb(exi(1)),
     >                           dcpl_mb(prjtmp(1)+(l-1)*npack1))
           end if
           call Pack_cc_indot(1,n,
     >                      psi1,
     >                      dcpl_mb(prjtmp(1)+(l-1)*npack1),
     >                      dbl_mb(sw1(1)+(l-1)*n))

        end do
        call D3dB_Vector_SumAll((n*nproj),dbl_mb(sw1(1)))


*       **** sw2 = Sijl*sw1 ******
        Gijl_indx = psi_data_get_ptr(int_mb(Gijl(1)+ia-1),2)
        call Multiply_Gijl_sw1(n,
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

*       *** routine needs to be parallelized over orbitals ****
*       **** S = S + sw1*transpose(sw2) ****
        call DGEMM('N','T',n,n,int_mb(nprj(1)+ia-1),
     >              (scal),
     >              dbl_mb(sw1(1)), n,
     >              dbl_mb(sw2(1)), n,
     >              (1.0d0),
     >              S, n)

        end if !** nproj>0 **
      end do !** ii **

      value =           MA_pop_stack(sw2(2))
      value = value.and.MA_pop_stack(sw1(2))
      value = value.and.MA_pop_stack(exi(2))
      if (.not.value) call errquit('psp_overlap_orb: popping stack',3,
     &       MA_ERR)
      call nwpw_timing_end(6)
      return 
      end




*     ***********************************
*     *					*
*     *	        psp_kinetic_core   	*
*     *					*
*     ***********************************
*
*    This routine returns the paw kinetic energy for the core density
*
      real*8 function psp_kinetic_core()
      implicit none

#include "mafdecls.fh"
#include "psp.fh"

*     *** local variables ***
      integer ii,ia
      real*8  ecore

*     **** external functions ****
      integer  ion_nion,ion_katm
      external ion_nion,ion_katm

      ecore = 0.0d0
      do ii=1,ion_nion()
        ia    = ion_katm(ii)
        if (int_mb(psp_type(1)+ia-1).eq.4) then
           ecore = ecore + dbl_mb(core_kin(1)+ia-1)
        end if
      end do

      psp_kinetic_core = ecore
      return 
      end




*     ***********************************
*     *					*
*     *	        psp_kinetic_atom	*
*     *					*
*     ***********************************

*    This routine computes the paw atomic kinetic energy
*
      real*8 function psp_kinetic_atom(ispin,ne,psi1)
      implicit none
      integer    ispin,ne(2)
      complex*16 psi1(*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "psp.fh"

*     *** local variables ***
      integer npack1,nion
      integer ii,ia,l,nn
      integer k,shift,l_prj,nproj,Gijl_indx
      real*8  omega,scal,scalsqr
      integer exi(2),sw1(2),sw2(2)
      logical value,sd_function

      real*8  kinetic_atom

*     **** external functions ****
      integer  ion_nion,ion_katm
      integer  psi_data_get_ptr
      real*8   lattice_omega
      external ion_nion,ion_katm
      external psi_data_get_ptr
      external lattice_omega

      kinetic_atom = 0.0d0

      nn = ne(1)+ne(2)

*     **** allocate local memory ****
      nion = ion_nion()
      call Pack_npack(1,npack1)

      value = MA_push_get(mt_dcpl,npack1,'exi', exi(2), exi(1))
      value = value.and.
     >        MA_push_get(mt_dbl,nn*nprj_max,'sw1',sw1(2),sw1(1))
      value = value.and.
     >        MA_push_get(mt_dbl,nn*nprj_max,'sw2',sw2(2),sw2(1))
      if (.not.value) 
     >  call errquit('psp_overlap_orb: out of stack',0,MA_ERR)

      omega = lattice_omega()
      scal    = 1.0d0/(omega)
      scalsqr = scal*scal


      do ii=1,nion
        ia    = ion_katm(ii)
        nproj = int_mb(nprj(1)+ia-1)

        if ((int_mb(psp_type(1)+ia-1).eq.4).and.(nproj.gt.0)) then

*       **** structure factor and local pseudopotential ****
        call strfac_pack(1,ii,dcpl_mb(exi(1)))

*       **** generate sw1's and projectors ****
        do l=1,nproj

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


*          **** phase factor does not matter therefore ****
*          **** (-i)^l is the same as (i)^l in the     ****
*          **** Rayleigh scattering formula            ****

*          *** current function is s or d ****
           if (sd_function) then
              call Pack_tc_Mul(1,dbl_mb(shift),
     >                           dcpl_mb(exi(1)),
     >                           dcpl_mb(prjtmp(1)+(l-1)*npack1))
*          *** current function is p or f ****
           else
              call Pack_tc_iMul(1,dbl_mb(shift),
     >                           dcpl_mb(exi(1)),
     >                           dcpl_mb(prjtmp(1)+(l-1)*npack1))
           end if
           call Pack_cc_indot(1,nn,
     >                      psi1,
     >                      dcpl_mb(prjtmp(1)+(l-1)*npack1),
     >                      dbl_mb(sw1(1)+(l-1)*nn))

        end do
        call D3dB_Vector_SumAll((nn*nproj),dbl_mb(sw1(1)))


*       **** sw2 = Tijl*sw1 ******
        Gijl_indx = psi_data_get_ptr(int_mb(Gijl(1)+ia-1),3)
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

*       **** keatom = transpose(sw1)*sw2) ****
        do l=0,(nn*nproj-1)
         kinetic_atom = kinetic_atom+dbl_mb(sw1(1)+l)*dbl_mb(sw2(1)+l)
        end do

        end if !** nproj>0 **
      end do !** ii **

      value =           MA_pop_stack(sw2(2))
      value = value.and.MA_pop_stack(sw1(2))
      value = value.and.MA_pop_stack(exi(2))
      if (.not.value) call errquit('psp_overlap_orb: popping stack',3,
     &       MA_ERR)

      if (ispin.eq.1) kinetic_atom = kinetic_atom+kinetic_atom
      kinetic_atom = kinetic_atom*scal

      psp_kinetic_atom = kinetic_atom
      return 
      end



