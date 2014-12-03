*
* $Id$
*

*     ***********************************
*     *					*
*     *	 	   psp_ld_trace		*
*     *					*
*     ***********************************

*    This routine returns the angular momentum ll local density of states of atom ii. It
* uses the pseudopotential projectors to define the local density of states.
*
*
      real*8 function psp_ld_trace(ii,ll,ispin,ne,psi1)
      implicit none
      integer    ii,ll
      integer    ispin,ne(2)
      complex*16 psi1(*)

#include "bafdecls.fh"
#include "psp.fh"
#include "errquit.fh"

*     *** local variables ***
      integer npack1
      integer i,j,ia,l,n,nn
      integer k,shift,l_prj,nproj,Gijl_indx
      integer exi(2),sw1(2),sw2(2)
      logical value,sd_function
      real*8 ld

*     **** external functions ****
      logical  is_sORd
      integer  ion_katm
      integer  psi_data_get_ptr,psi_data_get_chnk
      real*8   psp_ld_gen_trace
      external is_sORd
      external ion_katm
      external psi_data_get_ptr,psi_data_get_chnk
      external psp_ld_gen_trace

      call nwpw_timing_start(6) 

*     **** allocate local memory ****
      nn = ne(1)+ne(2)
      call Pack_npack(1,npack1)

      value = BA_push_get(mt_dcpl,npack1,'exi', exi(2), exi(1))
      value = value.and.
     >   BA_push_get(mt_dbl,nn*nprj_max*nprj_mult,'sw1',sw1(2),sw1(1))
      value = value.and.
     >   BA_push_get(mt_dbl,nn*nprj_max*nprj_mult,'sw2',sw2(2),sw2(1))
      if (.not.value) 
     >  call errquit('psp_ld: out of stack',0, MA_ERR)

      call dcopy(nn*nprj_max*nprj_mult,0.0d0,0,dbl_mb(sw1(1)),1)

      ld = 0.0d0

      ia=ion_katm(ii)
      nproj = int_mb(nprj(1)+ia-1)

*     **** structure factor and local pseudopotential ****
      call strfac_pack(1,ii,dcpl_mb(exi(1)))

*     **** generate sw1's and projectors ****
      do l=1,nproj
         shift = psi_data_get_ptr(int_mb(vnl(1)+ia-1),l)
         l_prj = int_mb(l_projector(1)+(l-1) 
     >                             + (ia-1)*(nmax_max*lmmax_max))

#ifdef GCC4
         k = iand(l_prj,1)
#else
         k = and(l_prj,1)
#endif
         sd_function = (k.eq.0)

*             *** current function is s or d ****
         if (sd_function) then
            call Pack_tc_Mul(1,dbl_mb(shift),
     >                         dcpl_mb(exi(1)),
     >                         dcpl_mb(prjtmp(1)+(l-1)*npack1))

*        *** current function is p or f ****
         else
            call Pack_tc_iMul(1,dbl_mb(shift),
     >                         dcpl_mb(exi(1)),
     >                         dcpl_mb(prjtmp(1)+(l-1)*npack1))

         end if
         call Pack_cc_indot(1,nn,
     >                 psi1,
     >                 dcpl_mb(prjtmp(1)+(l-1)*npack1),
     >                 dbl_mb(sw1(1)+(l-1)*nn))
      end do
      call D3dB_Vector_SumAll((nn*nproj*nprj_mult),dbl_mb(sw1(1)))

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

      ld = psp_ld_gen_trace(ispin,ne,nproj,
     >                      int_mb(l_projector(1)
     >                               +(ia-1)*(nmax_max*lmmax_max)),
     >                        int_mb(m_projector(1)
     >                              +(ia-1)*(nmax_max*lmmax_max)),
     >                        dbl_mb(sw2(1)),
     >                        ll)

      call D1dB_SumAll(ld)

      value = .true.
      value = value.and.BA_pop_stack(sw2(2))
      value = value.and.BA_pop_stack(sw1(2))
      value = value.and.BA_pop_stack(exi(2))
      if (.not.value) call errquit('psp_ld_trace: popping stack',3,
     &       MA_ERR)
      call nwpw_timing_end(6)

      psp_ld_trace = ld
      return 
      end

*     ***********************************
*     *					*
*     *	    psp_ld_trace_gradient	*
*     *					*
*     ***********************************
      subroutine psp_ld_trace_gradient(ii,ll,ispin,ne,psi1,
     >                                 fractional,occ,
     >                                 dv,dpsi1,move,fion)
      implicit none
      integer    ii,ll
      integer    ispin,ne(2)
      complex*16 psi1(*)
      logical    fractional
      real*8     occ(*)
      real*8     dv
      complex*16 dpsi1(*)  !*** dld/dpsi1*  ***
      logical    move
      real*8     fion(3,*)   !*** dld/dRii    ***

#include "bafdecls.fh"
#include "psp.fh"
#include "errquit.fh"

*     *** local variables ***
      integer npack1,G(3)
      integer i,j,ia,l,n,nn
      integer k,shift,l_prj,nproj,Gijl_indx
      integer exi(2),xtmp(2),sw1(2),sw2(2),sw3(2),asum(2)
      logical value,sd_function
      complex*16 ctmp
      real*8 scal,omega,ff(3)

*     **** external functions ****
      logical  is_sORd
      integer  ion_katm,Pack_G_indx
      integer  psi_data_get_ptr,psi_data_get_chnk
      real*8   psp_ld_gen_trace,lattice_omega
      external is_sORd
      external ion_katm,Pack_G_indx
      external psi_data_get_ptr,psi_data_get_chnk
      external psp_ld_gen_trace,lattice_omega

      call nwpw_timing_start(6)

*     **** allocate local memory ****
      nn = ne(1)+ne(2)
      call Pack_npack(1,npack1)

      value = BA_push_get(mt_dcpl,npack1,'exi', exi(2), exi(1))
      value = value.and.
     >   BA_push_get(mt_dbl,nn*nprj_max*nprj_mult,'sw1',sw1(2),sw1(1))
      value = value.and.
     >   BA_push_get(mt_dbl,nn*nprj_max*nprj_mult,'sw2',sw2(2),sw2(1))
      value = value.and.
     >   BA_push_get(mt_dbl,nn*nprj_max*nprj_mult,'sw3',sw3(2),sw3(1))
      value = value.and.BA_push_get(mt_dbl,npack1,
     >                               'xtmp',xtmp(2),xtmp(1))
      value = value.and.BA_push_get(mt_dbl,3*nn,'asum',asum(2),asum(1))
      if (.not.value)
     >  call errquit('psp_ld: out of stack',0, MA_ERR)

      G(1)  = Pack_G_indx(1,1)
      G(2)  = Pack_G_indx(1,2)
      G(3)  = Pack_G_indx(1,3)


      call dcopy(nn*nprj_max*nprj_mult,0.0d0,0,dbl_mb(sw1(1)),1)

      omega = lattice_omega()
      scal = 1.0d0/omega

      ia=ion_katm(ii)
      nproj = int_mb(nprj(1)+ia-1)

*     **** structure factor and local pseudopotential ****
      call strfac_pack(1,ii,dcpl_mb(exi(1)))

*     **** generate sw1's and projectors ****
      do l=1,nproj
         shift = psi_data_get_ptr(int_mb(vnl(1)+ia-1),l)
         l_prj = int_mb(l_projector(1)+(l-1)
     >                             + (ia-1)*(nmax_max*lmmax_max))

#ifdef GCC4
         k = iand(l_prj,1)
#else
         k = and(l_prj,1)
#endif
         sd_function = (k.eq.0)

*        *** current function is s or d ****
         if (sd_function) then
            call Pack_tc_Mul(1,dbl_mb(shift),
     >                         dcpl_mb(exi(1)),
     >                         dcpl_mb(prjtmp(1)+(l-1)*npack1))

*        *** current function is p or f ****
         else
            call Pack_tc_iMul(1,dbl_mb(shift),
     >                         dcpl_mb(exi(1)),
     >                         dcpl_mb(prjtmp(1)+(l-1)*npack1))

         end if
         call Pack_cc_indot(1,nn,
     >                      psi1,
     >                      dcpl_mb(prjtmp(1)+(l-1)*npack1),
     >                      dbl_mb(sw1(1)+(l-1)*nn))
      end do
      call Parallel_Vector_SumAll((nn*nproj*nprj_mult),dbl_mb(sw1(1)))

*     **** sw2 = Gijl*sw1 ******
      Gijl_indx = psi_data_get_ptr(int_mb(Gijl(1)+ia-1),1)
      call Multiply_Gijl_sw1(nn,
     >                    nproj,
     >                    int_mb(nmax(1)+ia-1),
     >                    int_mb(lmax(1)+ia-1),
     >                    int_mb(n_projector(1)
     >                           + (ia-1)*(nmax_max*lmmax_max)),
     >                    int_mb(l_projector(1)
     >                           + (ia-1)*(nmax_max*lmmax_max)),
     >                    int_mb(m_projector(1)
     >                           + (ia-1)*(nmax_max*lmmax_max)),
     >                    dbl_mb(Gijl_indx),
     >                    dbl_mb(sw1(1)),
     >                    dbl_mb(sw2(1)))

      call psp_ld_gen_sw3(ispin,ne,nproj,
     >                    int_mb(l_projector(1)
     >                           +(ia-1)*(nmax_max*lmmax_max)),
     >                    int_mb(m_projector(1)
     >                          +(ia-1)*(nmax_max*lmmax_max)),
     >                    ll,dbl_mb(sw2(1)),dbl_mb(sw3(1)))

      call Multiply_Gijl_sw1(nn,
     >                    nproj,
     >                    int_mb(nmax(1)+ia-1),
     >                    int_mb(lmax(1)+ia-1),
     >                    int_mb(n_projector(1)
     >                           + (ia-1)*(nmax_max*lmmax_max)),
     >                    int_mb(l_projector(1)
     >                           + (ia-1)*(nmax_max*lmmax_max)),
     >                    int_mb(m_projector(1)
     >                           + (ia-1)*(nmax_max*lmmax_max)),
     >                    dbl_mb(Gijl_indx),
     >                    dbl_mb(sw3(1)),
     >                    dbl_mb(sw2(1)))

*     **** do Kleinman-Bylander Multiplication ****
      call dscal(nn*int_mb(nprj(1)+ia-1)*nprj_mult,
     >           dv*scal,dbl_mb(sw2(1)),1)

      call DGEMM('N','T',2*npack1,nn,int_mb(nprj(1)+ia-1)*nprj_mult,
     >           (1.0d0),
     >           dcpl_mb(prjtmp(1)), 2*npack1,
     >           dbl_mb(sw2(1)),     nn,
     >           (1.0d0),
     >           dpsi1,              2*npack1)


      if (move) then
      do l=1,nproj
         do n=1,nn
            if (ispin.eq.1) 
     >         dbl_mb(sw2(1)+n-1+(l-1)*nn)
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
     >                           *dbl_mb(asum(1)+  3*(n-1))
            ff(2) = ff(2) + 2.0d0*dbl_mb(sw2(1)+n-1+(l-1)*nn) !// change
     >                           *dbl_mb(asum(1)+1+3*(n-1))
            ff(3) = ff(3) + 2.0d0*dbl_mb(sw2(1)+n-1+(l-1)*nn) !// change
     >                           *dbl_mb(asum(1)+2+3*(n-1))
         end do
         call D1dB_Vector_SumAll(3,ff)
         fion(1,ii) = fion(1,ii)  - ff(1)
         fion(2,ii) = fion(2,ii)  - ff(2)
         fion(3,ii) = fion(3,ii)  - ff(3)

      end do !** l **
      end if

      value = .true.
      value = value.and.BA_pop_stack(asum(2))
      value = value.and.BA_pop_stack(xtmp(2))
      value = value.and.BA_pop_stack(sw3(2))
      value = value.and.BA_pop_stack(sw2(2))
      value = value.and.BA_pop_stack(sw1(2))
      value = value.and.BA_pop_stack(exi(2))
      if (.not.value) 
     >   call errquit('psp_ld_trace_gradient: popping stack',3,
     &                MA_ERR)
      call nwpw_timing_end(6)
      return 
      end


*     ****************************************************
*     *                                                  *
*     *                psp_ld_gen_trace                  *
*     *                                                  *
*     ****************************************************
      real*8 function psp_ld_gen_trace(ispin,ne,
     >                                     nprj,l_prj,m_prj,sw2,l)
      implicit none
      integer ispin,ne(2),nprj
      integer l_prj(nprj),m_prj(nprj)
      real*8  sw2(ne(1)+ne(2),nprj)
      integer l

*     **** local variables ****
      real*8  omega,ld
      integer a,ma,ms,n,n1(2),n2(2)

*     **** external functions ****
      real*8   lattice_omega
      external lattice_omega

      omega = lattice_omega()
      n1(1) = 1
      n1(2) = ne(1)+1
      n2(1) = ne(1)
      n2(2) = ne(1)+ne(2)

      ld = 0.0d0
      do a=1,nprj
      if ((l_prj(a).eq.l)) then
         ma = m_prj(a)+l+1
         do ms=1,ispin
         do n=n1(ms),n2(ms)
            ld = ld + sw2(n,a)*sw2(n,a)/omega
         end do
         end do
      end if
      end do
      if (ispin.eq.1) ld = ld + ld

      psp_ld_gen_trace = ld
      return
      end 




*     ***********************************
*     *					*
*     *	 	   psp_ld_matrix	*
*     *					*
*     ***********************************

*    This routine returns the ll local density of states of atom ii. It
* uses the pseudopotential projectors to define the local density of states.
*
*
      subroutine psp_ld_matrix(ii,ll,ispin,ne,psi1,ld_density)
      implicit none
      integer    ii,ll
      integer    ispin,ne(2)
      complex*16 psi1(*)
      real*8     ld_density(*)

#include "bafdecls.fh"
#include "psp.fh"
#include "errquit.fh"


*     *** local variables ***
      integer npack1
      integer i,j,ia,l,n,nn
      integer k,shift,l_prj,nproj,Gijl_indx
      integer exi(2),xtmp(2),sw1(2),sw2(2),asum(2)
      logical value,sd_function
      real*8 ld

*     **** external functions ****
      logical  is_sORd
      integer  ion_katm
      integer  psi_data_get_ptr,psi_data_get_chnk
      external is_sORd
      external ion_katm
      external psi_data_get_ptr,psi_data_get_chnk

      call nwpw_timing_start(6) 

*     **** allocate local memory ****
      nn = ne(1)+ne(2)
      call Pack_npack(1,npack1)

      value = BA_push_get(mt_dcpl,npack1,'exi', exi(2), exi(1))
      value = value.and.
     >   BA_push_get(mt_dbl,nn*nprj_max*nprj_mult,'sw1',sw1(2),sw1(1))
      value = value.and.
     >   BA_push_get(mt_dbl,nn*nprj_max*nprj_mult,'sw2',sw2(2),sw2(1))
      if (.not.value) 
     >  call errquit('psp_ld_matrix: out of stack',0, MA_ERR)

      call dcopy(nn*nprj_max*nprj_mult,0.0d0,0,dbl_mb(sw1(1)),1)

      ld = 0.0d0

      ia=ion_katm(ii)
      nproj = int_mb(nprj(1)+ia-1)

*     **** structure factor and local pseudopotential ****
      call strfac_pack(1,ii,dcpl_mb(exi(1)))

*     **** generate sw1's and projectors ****
      do l=1,nproj
         shift = psi_data_get_ptr(int_mb(vnl(1)+ia-1),l)
         l_prj = int_mb(l_projector(1)+(l-1) 
     >                             + (ia-1)*(nmax_max*lmmax_max))

#ifdef GCC4
         k = iand(l_prj,1)
#else
         k = and(l_prj,1)
#endif
         sd_function = (k.eq.0)

*             *** current function is s or d ****
         if (sd_function) then
            call Pack_tc_Mul(1,dbl_mb(shift),
     >                         dcpl_mb(exi(1)),
     >                         dcpl_mb(prjtmp(1)+(l-1)*npack1))

*        *** current function is p or f ****
         else
            call Pack_tc_iMul(1,dbl_mb(shift),
     >                         dcpl_mb(exi(1)),
     >                         dcpl_mb(prjtmp(1)+(l-1)*npack1))

         end if
         call Pack_cc_indot(1,nn,
     >                 psi1,
     >                 dcpl_mb(prjtmp(1)+(l-1)*npack1),
     >                 dbl_mb(sw1(1)+(l-1)*nn))
      end do
      call D3dB_Vector_SumAll((nn*nproj*nprj_mult),dbl_mb(sw1(1)))

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

      call psp_ld_gen_matrix(ispin,ne,nproj,
     >                      int_mb(l_projector(1)
     >                               +(ia-1)*(nmax_max*lmmax_max)),
     >                        int_mb(m_projector(1)
     >                              +(ia-1)*(nmax_max*lmmax_max)),
     >                        dbl_mb(sw2(1)),
     >                        ll,ld_density)

      call D1dB_Vector_SumAll(ispin*(2*ll+1)**2,ld_density)

      value = .true.
      value = value.and.BA_pop_stack(sw2(2))
      value = value.and.BA_pop_stack(sw1(2))
      value = value.and.BA_pop_stack(exi(2))
      if (.not.value) call errquit('psp_ld_matrix: popping stack',3,
     &       MA_ERR)
      call nwpw_timing_end(6)

      return 
      end





*     ****************************************************
*     *                                                  *
*     *                psp_ld_gen_matrix                 *
*     *                                                  *
*     ****************************************************
      subroutine psp_ld_gen_matrix(ispin,ne,nprj,l_prj,m_prj,sw2,
     >                             l,l_density)
      implicit none
      integer ispin,ne(2),nprj
      integer l_prj(nprj),m_prj(nprj)
      real*8  sw2(ne(1)+ne(2),nprj)
      integer l
      real*8  l_density(2*l+1,2*l+1,ispin)

*     **** local variables ****
      integer a,b,ma,mb,ms,n,n1(2),n2(2),lm
      real*8  omega

*     **** external functions ****
      real*8   lattice_omega
      external lattice_omega

      lm = 2*l+1
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
      if (ispin.eq.1) call dscal(ispin*lm*lm,2.0d0,l_density,1)

      return
      end




*     ****************************************************
*     *                                                  *
*     *                psp_ld_gen_sw3                    *
*     *                                                  *
*     ****************************************************
      subroutine psp_ld_gen_sw3(ispin,ne,nprj,l_prj,m_prj,
     >                          l,sw2,sw3)
      implicit none
      integer ispin,ne(2),nprj
      integer l_prj(nprj),m_prj(nprj)
      integer l
c      real*8  vmm(2*l+1,2*l+1,ispin)
      real*8  sw2(ne(1)+ne(2),nprj)
      real*8  sw3(ne(1)+ne(2),nprj)

*     **** local variables ****
      integer a,ms,i,lm
  
      call dcopy((ne(1)+ne(2))*nprj,0.0d0,0,sw3,1)
      do a=1,nprj
         if ((l_prj(a).eq.l)) then
            do i=1,ne(1)+ne(2)
               sw3(i,a) =  sw2(i,a)
            end do
         end if
      end do

      return
      end

