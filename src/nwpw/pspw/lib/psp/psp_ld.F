*
* $Id: psp_ld.F 22647 2012-07-28 00:03:24Z bylaska $
*


*     ***********************************
*     *					*
*     *	 	   psp_ld		*
*     *					*
*     ***********************************

*    This routine returns the ll local density of states of atom ii. It
* uses the pseudopotential projectors to define the local density of states.
*
*
      real*8 function psp_ld(ii,ll,ispin,ne,psi1)
      implicit none
      integer    ii,ll
      integer    ispin,ne(2)
      complex*16 psi1(*)

#include "mafdecls.fh"
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
      real*8   psp_ld_gen_l_density
      external is_sORd
      external ion_katm
      external psi_data_get_ptr,psi_data_get_chnk
      external psp_ld_gen_l_density

      call nwpw_timing_start(6) 

*     **** allocate local memory ****
      nn = ne(1)+ne(2)
      call Pack_npack(1,npack1)

      value = MA_push_get(mt_dcpl,npack1,'exi', exi(2), exi(1))
      value = value.and.
     >   MA_push_get(mt_dbl,nn*nprj_max*nprj_mult,'sw1',sw1(2),sw1(1))
      value = value.and.
     >   MA_push_get(mt_dbl,nn*nprj_max*nprj_mult,'sw2',sw2(2),sw2(1))
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

      ld = psp_ld_gen_l_density(ispin,ne,nproj,
     >                        int_mb(l_projector(1)
     >                               +(ia-1)*(nmax_max*lmmax_max)),
     >                        int_mb(m_projector(1)
     >                              +(ia-1)*(nmax_max*lmmax_max)),
     >                        dbl_mb(sw2(1)),
     >                        ll)

      call D1dB_SumAll(ld)

      value = .true.
      value = value.and.MA_pop_stack(sw2(2))
      value = value.and.MA_pop_stack(sw1(2))
      value = value.and.MA_pop_stack(exi(2))
      if (.not.value) call errquit('psp_ld: popping stack',3,
     &       MA_ERR)
      call nwpw_timing_end(6)

      psp_ld = ld
      return 
      end



*     ****************************************************
*     *                                                  *
*     *                psp_ld_gen_l_density              *
*     *                                                  *
*     ****************************************************
      real*8 function psp_ld_gen_l_density(ispin,ne,
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

      psp_ld_gen_l_density = ld
      return
      end 
