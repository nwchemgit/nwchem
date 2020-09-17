c
c $Id$
c


********************************************************************
*
*	pspw_APC module: used to generate derived atomic point charges
*         	         from a plane-wave density.
*
*   The algorithms used in this module are based on the work of
*   P.E. Blochl, J. Chem. Phys. vol. 103, page 7422 (1995).
*
*	Data structure
*
********************************************************************

*     ***********************************
*     *					*
*     *		pspw_init_APC		*
*     *					*
*     ***********************************
      subroutine pspw_init_APC(rtdb)
      implicit none
      integer rtdb

#include "bafdecls.fh"
#include "btdb.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "pspw_APC.fh"

*     **** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)

      logical value,doapc,oprint
      integer i,j,k,n,npack0,nfft3d
      integer nx,ny,nxh
      integer zero,pzero,qzero
      integer tmp1(2),G(3)
      real*8  coef,x,gg,fourpi

*     **** external functions ****
      logical  control_Mulliken,control_print
      logical  nwpw_born_on,nwpw_cosmo2_on,nwpw_cdft_on
      integer  ion_nion_qm,G_indx
      real*8   lattice_omega
      external control_Mulliken,control_print
      external nwpw_born_on,nwpw_cosmo2_on,nwpw_cdft_on
      external ion_nion_qm,G_indx
      external lattice_omega

      if (.not.btdb_get(rtdb,'nwpw:APC',mt_log,1,doapc))
     >  doapc = control_Mulliken().or.
     >          nwpw_born_on().or.
     >          nwpw_cosmo2_on().or.
     >          nwpw_cdft_on()

*     **** read in nga from rtdb ***
      value = btdb_get(rtdb,'nwpw_APC:nga',mt_int,1,nga)

*     **** return and set nga to zero if nothing in rtdb and not mulliken****
      if ((.not.value).or.(nga.le.0)) then
         if (doapc) then
            nga = 3
         else
            nga = 0
            return
         end if
      end if
      ngs = nga*ion_nion_qm()
      call Parallel_taskid(taskid)
      oprint = (taskid.eq.MASTER).and.control_print(print_debug)

      call Pack_npack(0,npack0)

*     **** Allocate memory needed for pspw_APC ****
      value = BA_alloc_get(mt_dbl,ngs*ngs,'AAPC',A(2),A(1))
      value = value.and.
     >          BA_alloc_get(mt_dbl,ngs*ngs,'AmAPC',Am(2),Am(1))
      value = value.and.
     >          BA_alloc_get(mt_dbl,ngs,'bAPC',b(2),b(1))
      value = value.and.
     >          BA_alloc_get(mt_dbl,ngs,'qAPC',q(2),q(1))
      value = value.and.
     >          BA_alloc_get(mt_dbl,nga,'gammaAPC',gamma(2),gamma(1))
      value = value.and.
     >          BA_alloc_get(mt_dbl,npack0,'wAPC',w(2),w(1))
      value = value.and.
     >          BA_alloc_get(mt_dbl,npack0*nga,'gausAPC',
     >                       gaus(2),gaus(1))

*     **** EXIT IF OUT OF MEMORY ****
      if (.not. value)
     >   call errquit('pspw_init_APC:out of heap memory',0, MA_ERR)

*     **** read in Gc and gamma from rtdb ***
      if (.not.btdb_get(rtdb,'nwpw_APC:Gc',mt_dbl,1,Gc)) Gc = 2.5d0
      if (.not.btdb_get(rtdb,'nwpw_APC:gamma',mt_dbl,nga,
     >   dbl_mb(gamma(1)))) then
         if (nga.eq.3) then
            dbl_mb(gamma(1))   = 0.6d0
            dbl_mb(gamma(1)+1) = 0.9d0
            dbl_mb(gamma(1)+2) = 1.35d0
         else
            call errquit('pspw_init_APC:error reading rtdb',0,RTDB_ERR)
         end if
      end if


*     ********************************
*     **** define weight function ****
*     ********************************

      fourpi = 4.0d0*(4.0d0*datan(1.0d0))
      call D3dB_nfft3d(1,nfft3d)
      G(1)= G_indx(1)
      G(2)= G_indx(2)
      G(3)= G_indx(3)

*     ***** find the G==0 point in the lattice *****
      i=0
      j=0
      k=0
      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      nxh = nx/2
c     call D3dB_ktoqp(1,k+1,qzero,pzero)
c     zero = (qzero-1)*(nxh+1)*ny
c    >     + j*(nxh+1)
c    >     + i+1
      call D3dB_ijktoindexp(1,i+1,j+1,k+1,zero,pzero)

      value = BA_push_get(mt_dbl,nfft3d,'tmp1',tmp1(2),tmp1(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

      do i = 1,nfft3d
         gg  = ( dbl_mb(G(1)+i-1)*dbl_mb(G(1)+i-1)
     >         + dbl_mb(G(2)+i-1)*dbl_mb(G(2)+i-1)
     >         + dbl_mb(G(3)+i-1)*dbl_mb(G(3)+i-1))

         dbl_mb(tmp1(1)+i-1) = 0.0d0
         if ((gg.gt.1.0d-6).and.(gg.lt.(Gc*Gc))) then
            dbl_mb(tmp1(1)+i-1) = fourpi*(gg-Gc*Gc)**2/(gg*Gc*Gc)
         end if
      end do
      call Pack_t_pack(0,dbl_mb(tmp1(1)))
      call Pack_t_Copy(0,dbl_mb(tmp1(1)),dbl_mb(w(1)))
      value = BA_pop_stack(tmp1(2))
      if (.not. value) call errquit('popping of stack memory',0, MA_ERR)


*     *************************************
*     **** define Gaussians in G-space ****
*     *************************************

      value = BA_push_get(mt_dbl,nfft3d,'tmp1',tmp1(2),tmp1(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)


      coef = 1.0d0/lattice_omega()
      do n=1,nga
         x = dbl_mb(gamma(1)+(n-1))
         x = x*x/4.0d0
         do i = 1,nfft3d
            gg  = ( dbl_mb(G(1)+i-1)*dbl_mb(G(1)+i-1)
     >            + dbl_mb(G(2)+i-1)*dbl_mb(G(2)+i-1)
     >            + dbl_mb(G(3)+i-1)*dbl_mb(G(3)+i-1))
            dbl_mb(tmp1(1)+i-1) = coef*exp(-x*gg)
         end do
         call Pack_t_pack(0,dbl_mb(tmp1(1)))
         call Pack_t_Copy(0,dbl_mb(tmp1(1)),
     >                      dbl_mb(gaus(1)+(n-1)*npack0))
      end do

      value = BA_pop_stack(tmp1(2))
      if (.not. value) call errquit('popping of stack memory',0, MA_ERR)

*     **** DEBUG - write out information *****
      if (oprint) then
         write(luout,*) 
         write(luout,*) 'initializing pspw_APC data structure'
         write(luout,*) '------------------------------------'
         write(luout,*) 'nga, ngs:',nga,ngs
         write(luout,*) 'Gc      :',Gc
         do i=1,nga
           write(luout,*) 'APC gamma:',i, dbl_mb(gamma(1)+i-1)
         end do
      end if

      return
      end

*     ***********************************
*     *					*
*     *		pspw_end_APC		*
*     *					*
*     ***********************************
      subroutine pspw_end_APC()
      implicit none
#include "errquit.fh"

#include "bafdecls.fh"
#include "pspw_APC.fh"

      logical value

      
      if (nga.gt.0) then
        value = BA_free_heap(A(2))
        value = value.and.BA_free_heap(Am(2))
        value = value.and.BA_free_heap(b(2))
        value = value.and.BA_free_heap(q(2))
        value = value.and.BA_free_heap(gamma(2))
        value = value.and.BA_free_heap(w(2))
        value = value.and.BA_free_heap(gaus(2))
         if (.not. value) call errquit('error freeing heap',0, MA_ERR)
        nga = 0
        ngs = 0
      end if

      return
      end

      subroutine pspw_gen_db_APC_d(ispin,ne,dng,nion,dbdR)
      implicit none
      integer ispin,ne(2)
      complex*16 dng(*)
      integer nion
      real*8 dbdR(nion,*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "pspw_APC.fh"

      return
      end


*     ***********************************
*     *					*
*     *		pspw_gen_APC		*
*     *					*
*     ***********************************
      subroutine pspw_gen_APC(ispin,ne,dng)
      implicit none
      integer ispin,ne(2)
      complex*16 dng(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "pspw_APC.fh"

*     ***** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)

      logical value
      integer i,j,ii,iii,jj,jjj,indx,indxt
      integer nfft3d,npack0
      real*8 omega,N,sum,sum1

c      integer G(3)
      integer exi(2),exj(2)
      integer gaus_i(2),gaus_j(2),ipiv(2)
      integer AAA(2),work(2),lwork,rank,info

*     ***** external functions ****
      integer  ion_nion_qm,Pack_G_indx
      external ion_nion_qm,Pack_G_indx
      real*8   lattice_omega
      external lattice_omega

      if (nga.gt.0) then
         omega = lattice_omega()
         call Parallel_taskid(taskid)
         call D3dB_nfft3d(1,nfft3d)
         call Pack_npack(0,npack0)

c         G(1)  = Pack_G_indx(0,1)
c         G(2)  = Pack_G_indx(0,2)
c         G(3)  = Pack_G_indx(0,3)

*         **** get memory from stack ****
         value = BA_push_get(mt_dcpl,nfft3d,'ttexi',exi(2),exi(1))
         value = value.and.
     >           BA_push_get(mt_dcpl,nfft3d,'ttexj',exj(2),exj(1))
         value = value.and.
     >           BA_push_get(mt_dcpl,npack0,'ghaus_i1',
     >                       gaus_i(2),gaus_i(1))
         value = value.and.
     >           BA_push_get(mt_dcpl,npack0,'ghaus_j2',
     >                                       gaus_j(2),gaus_j(1))
         if (.not. value) call errquit('out of stack memory',0, MA_ERR)


*        ***** calculate N = dng(G=0)*omega *****
         N = dble(ne(1)+ne(ispin))

*        ***** calculate c_i = omega*gaus_i(G=0) = 1.0d0 *****

*        ***** calculate b_i = omega*Sum(G) w(G)*Re(dcongj(dng(G))*gaus_i(G)) ****
         do ii=1,ion_nion_qm()

*           **** structure factor  ****
            call strfac(ii,dcpl_mb(exi(1)))
            call Pack_c_pack(0,dcpl_mb(exi(1)))

            do iii=1,nga
               i = iii + (ii-1)*nga

*              *** gaus_i(G))***
               call Pack_tc_Mul(0,dbl_mb(gaus(1)+npack0*(iii-1)),
     >                          dcpl_mb(exi(1)),
     >                          dcpl_mb(gaus_i(1)))

*              *** w(G)*gaus_i(G))***
c               call Pack_tc_Mul(0,dbl_mb(w(1)),
c     >                          dcpl_mb(gaus_i(1)),
c     >                          dcpl_mb(gaus_i(1)))
               call Pack_tc_Mul2(0,dbl_mb(w(1)),dcpl_mb(gaus_i(1)))

*              *** omega*Sum(G) w(G)*Re(dcongj(dng(G))*gaus_i(G))***
               call Pack_cc_dot(0,dng,
     >                            dcpl_mb(gaus_i(1)),
     >                            sum)
!$OMP MASTER
               dbl_mb(b(1)+i-1) = sum*omega
!$OMP END MASTER


c               call Pack_cct_iconjgMulb(0,
c     >                            dng,
c     >                            dcpl_mb(gaus_i(1)),
c     >                            dbl_mb(xtmp(1)))
c               call Pack_t_dscal1(0,omega,dbl_mb(xtmp(1)))
c
c               call Pack_tt_dot(0,dbl_mb(G(1)),
c     >                           dbl_mb(xtmp(1)),
c     >                           dbl_mb(dbdR(1)+3*(i-1)))
c               call Pack_tt_dot(0,dbl_mb(G(2)),
c     >                           dbl_mb(xtmp(1)),
c     >                           dbl_mb(dbdR(1)+3*(i-1)+1))
c               call Pack_tt_dot(0,
c     >                           dbl_mb(G(3)),
c     >                           dbl_mb(xtmp(1)),
c     >                           dbl_mb(dbdR(1)+3*(i-1)+2))
            end do
          end do



*        ***** calculate A_ij = omega*Sum(G) w(G)*dcongj(gaus_i(G))*gaus_j(G)) ****
         do ii=1,ion_nion_qm()
*           **** structure factor  ****
            call strfac(ii,dcpl_mb(exi(1)))
            call Pack_c_pack(0,dcpl_mb(exi(1)))

             do jj=ii,ion_nion_qm()
*              **** structure factor  ****
               call strfac(jj,dcpl_mb(exj(1)))
               call Pack_c_pack(0,dcpl_mb(exj(1)))


               do iii=1,nga
*                 *** gaus_i(G))***
                  call Pack_tc_Mul(0,dbl_mb(gaus(1)+npack0*(iii-1)),
     >                               dcpl_mb(exi(1)),
     >                               dcpl_mb(gaus_i(1)))

*                 *** w(G)*gaus_i(G))***
c                  call Pack_tc_Mul(0,dbl_mb(w(1)),
c     >                               dcpl_mb(gaus_i(1)),
c     >                               dcpl_mb(gaus_i(1)))
                  call Pack_tc_Mul2(0,dbl_mb(w(1)),dcpl_mb(gaus_i(1)))

                  do jjj=1,nga
*                   *** gaus_j(G))***
                     call Pack_tc_Mul(0,dbl_mb(gaus(1)+npack0*(jjj-1)),
     >                                  dcpl_mb(exj(1)),
     >                                  dcpl_mb(gaus_j(1)))

*                   *** omega*Sum(G) w(G)*gaus_i(G)*gaus_j(G))***
                    call Pack_cc_dot(0,dcpl_mb(gaus_i(1)),
     >                                 dcpl_mb(gaus_j(1)),
     >                                 sum)
                    sum = sum*omega
                    i = iii + (ii-1)*nga
                    j = jjj + (jj-1)*nga

                    indx  = (i-1) + (j-1)*ngs
                    indxt = (j-1) + (i-1)*ngs
!$OMP MASTER
                    dbl_mb(A(1)+indx)  = sum
                    dbl_mb(A(1)+indxt) = sum
!$OMP END MASTER
               end do
             end do
           end do
         end do
         value = BA_pop_stack(gaus_j(2))
         value = value.and.BA_pop_stack(gaus_i(2))
         value = value.and.BA_pop_stack(exj(2))
         value = value.and.BA_pop_stack(exi(2))
         if (.not. value) call errquit('popping of stack memory',0,
     &       MA_ERR)

*        **** perform matrix operations in serial ****
         !call dcopy(ngs*ngs,0.0d0,0,dbl_mb(Am(1)),1)
         !call dcopy(ngs,0.0d0,0,dbl_mb(q(1)),1)
         call Parallel_shared_vector_zero(.true.,ngs*ngs,dbl_mb(Am(1)))
         call Parallel_shared_vector_zero(.true.,ngs,dbl_mb(q(1)))

         if (taskid.eq.MASTER) then

           lwork = 5*ngs*ngs
           !value = BA_push_get(mt_int,ngs,'ipivAPC',ipiv(2),ipiv(1))
           value = BA_push_get(mt_dbl,lwork,'work',work(2),work(1))
           value = value.and.
     >             BA_push_get(mt_dbl,ngs*ngs,'AAAAPC',AAA(2),AAA(1))
           if (.not. value) call errquit('out of stack memory',0,
     &       MA_ERR)

*          **** calculate Am_ij ****
           !call dcopy(ngs*ngs,dbl_mb(A(1)),1,dbl_mb(AAA(1)),1)
           !call dcopy(ngs*ngs,0.0d0,0,dbl_mb(Am(1)),1)
           call Parallel_shared_vector_copy(.true.,ngs*ngs,
     >                                      dbl_mb(A(1)),dbl_mb(AAA(1)))
           call Parallel_shared_vector_zero(.true.,ngs*ngs,
     >                                      dbl_mb(Am(1)))
!$OMP MASTER
           do i=1,ngs
              indx = i + (i-1)*ngs
              dbl_mb(Am(1)+indx-1) = 1.0d0
              !int_mb(ipiv(1)+i-1) = 0
           end do
c           call DGESV(ngs,ngs,dbl_mb(AAA(1)), ngs,
c     >                        int_mb(ipiv(1)),
c     >                        dbl_mb(Am(1)),ngs,
c     >                        j)

           call DGELSS(ngs,ngs,ngs,dbl_mb(AAA(1)),ngs,
     >                        dbl_mb(Am(1)),ngs,
     >                        dbl_mb(q(1)),1.0d-9,
     >                        rank,
     >                        dbl_mb(work(1)),lwork,
     >                        info)
!$OMP END MASTER
            value = BA_pop_stack(AAA(2))
            value = value.and.BA_pop_stack(work(2))
            !value = value.and.BA_pop_stack(ipiv(2))
            if (.not. value) call errquit('popping stack memory',0,
     &       MA_ERR)

           !call dcopy(ngs,0.0d0,0,dbl_mb(q(1)),1)
           call Parallel_shared_vector_zero(.true.,ngs,dbl_mb(q(1)))
*          **** calculate q_i ****
           sum  = 0.0d0
           sum1 = 0.0d0
           do j=1,ngs
           do i=1,ngs
              indx = (i-1) + (j-1)*ngs
              sum  = sum  + dbl_mb(Am(1)+indx)*dbl_mb(b(1)+j-1)
              sum1 = sum1 + dbl_mb(Am(1)+indx)
           end do
           end do
           sum = (sum-N)/sum1
   
           do i=1,ngs
              sum1 = 0.0d0
              do j=1,ngs
                indx = (i-1) + (j-1)*ngs
                sum1 = sum1 
     >               + dbl_mb(Am(1)+indx)*(dbl_mb(b(1)+j-1)-sum)
              end do
!$OMP MASTER
              dbl_mb(q(1)+i-1) = sum1
!$OMP END MASTER
           end do
         end if
*        **** synchronization ****
         call D3dB_Vector_SumAll(ngs,dbl_mb(q(1)))
         call D3dB_Vector_SumAll(ngs*ngs,dbl_mb(Am(1)))

      end if

      return
      end

*     ***********************************
*     *					*
*     *		pspw_dngen_APC		*
*     *					*
*     ***********************************
      subroutine pspw_dngen_APC(ispin,ne,dn)
      implicit none
      integer ispin,ne(2)
      real*8 dn(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "pspw_APC.fh"

*     ***** local variables ****
      logical value
      integer nx,ny,nz,n2ft3d,nfft3d,npack0
      real*8 scal1

      integer dng(2),tmp(2)

      if (nga.gt.0) then

      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      scal1 = 1.0d0/dble(nx*ny*nz)

      call Pack_npack(0,npack0)
      call D3dB_nfft3d(1,nfft3d)
      n2ft3d = 2*nfft3d

*     **** get memory from stack ****
      value = BA_push_get(mt_dcpl,npack0,'dngakjs',dng(2),dng(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,nfft3d,'tmpxkjz',tmp(2),tmp(1))
      if (.not. value) call errquit('out of stack memory',0,
     &       MA_ERR)

*     **** generate dng ****
      call D3dB_rr_Sum(1,dn(1),dn(1+(ispin-1)*n2ft3d),dcpl_mb(tmp(1)))
c      call D3dB_r_SMul(1,scal1,dcpl_mb(tmp(1)),dcpl_mb(tmp(1)))
      call D3dB_r_SMul1(1,scal1,dcpl_mb(tmp(1)))
      call D3dB_rc_fft3f(1,dcpl_mb(tmp(1)))
      call Pack_c_pack(0,dcpl_mb(tmp(1)))
      call Pack_c_Copy(0,dcpl_mb(tmp(1)),dcpl_mb(dng(1)))
      value = BA_pop_stack(tmp(2))
      if (.not. value) call errquit('popping of stack memory',0,
     &       MA_ERR)

*     **** generate APC *****
      call pspw_gen_APC(ispin,ne,dcpl_mb(dng(1)))


      value = value.and.BA_pop_stack(dng(2))
      if (.not. value) call errquit('popping of stack memory',0,
     &       MA_ERR)

      end if

      return
      end



*     ***********************************
*     *                                 *
*     *         pspw_sumAm_APC          *
*     *                                 *
*     ***********************************
      real*8 function pspw_sumAm_APC(ngs,Am)
      implicit none
      integer ngs
      real*8 Am(ngs,ngs)

*     **** local variables ***
      integer i,j
      real*8  sum0

!$OMP MASTER
      sum0 = 0.0d0
      do i=1,ngs
         do j=1,ngs
            sum0 = sum0 + Am(i,j)
         end do
      end do
!$OMP END MASTER
!$OMP BARRIER

      pspw_sumAm_APC = sum0
      return
      end

*     ***********************************
*     *                                 *
*     *         pspw_Amtimesu_APC       *
*     *                                 *
*     ***********************************
      real*8 function pspw_Amtimesu_APC(ngs,Am,u)
      implicit none
      integer ngs
      real*8 Am(ngs,ngs)
      real*8 u(ngs)

*     **** local variables ***
      integer i,j
      real*8  sum0

!$OMP MASTER
      sum0 = 0.0d0
      do i=1,ngs
         do j=1,ngs
            sum0 = sum0 + Am(i,j)*u(j)
         end do
      end do
!$OMP END MASTER
!$OMP BARRIER

      pspw_Amtimesu_APC = sum0
      return
      end

*     ***********************************
*     *                                 *
*     *         pspw_Vfac_APC           *
*     *                                 *
*     ***********************************
      real*8 function pspw_Vfac_APC(ngs,Am,u,i)
      implicit none
      integer ngs
      real*8 Am(ngs,ngs)
      real*8 u(ngs)
      integer i

*     **** local variables ****
      integer j
      real*8 sum0

!$OMP MASTER
      sum0 = 0.0d0
      do j=1,ngs
         sum0 = sum0 + Am(i,j)*u(j)
      end do
!$OMP END MASTER
!$OMP BARRIER

      pspw_Vfac_APC = sum0
      return
      end


*     ***********************************
*     *					*
*     *		pspw_VQ_APC		*
*     *					*
*     ***********************************
*
*   This routine calculates dE/drho(G) where E is a function of model q. In order 
* to use this routine Am = inverse A must be calculated
* for the current geometry.
*
*   Entry - nion: number of qm atoms
*         - u: dE/dq(ii) - derivative of E wrt to model charges q(ii)
*   Exit - VQ(G) = dE/dq(ii)*dq(ii)/drho(G)
*
*   Note - pspw_gen_APC needs to be called to Am=inv(A) before 
*          this routine is called.
*
      subroutine pspw_VQ_APC(nion,u,VQ)
      implicit none
      integer    nion
      real*8     u(*)
      complex*16 VQ(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "pspw_APC.fh"

*     ***** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)

      logical value
      integer i,j,iii,ii,tid,nthreads
      integer npack0
      real*8  sumAm,afac,sumAmU,omega

      integer exi(2)
      integer gaus_i(2)

*     ***** external functions ****
      real*8   lattice_omega,pspw_Vfac_APC
      external lattice_omega,pspw_Vfac_APC
      real*8   pspw_sumAm_APC,pspw_Amtimesu_APC
      external pspw_sumAm_APC,pspw_Amtimesu_APC
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads

      if (nga.gt.0) then
         omega = lattice_omega()
         call Parallel_taskid(taskid)
         call Pack_npack(0,npack0)
         tid = Parallel_threadid()
         nthreads = Parallel_nthreads()

*        **** get memory from stack ****
         value = BA_push_get(mt_dcpl,npack0,'ttexi',exi(2),exi(1))
         value = value.and.
     >           BA_push_get(mt_dcpl,npack0,'ghaus_i1',
     >                       gaus_i(2),gaus_i(1))
         if (.not. value) 
     >      call errquit('pspw_VQ_APC:push stack',0, MA_ERR)

         sumAm  = pspw_sumAm_APC(ngs,dbl_mb(Am(1)))
         sumAmU = pspw_Amtimesu_APC(ngs,dbl_mb(Am(1)),u)

         do i=tid+1,ngs,nthreads
            u(i) = u(i) - (sumAmU/SumAm)
         end do


*        ***** calculate VQ(G) ****
         !***call Parallel_shared_vector_zero(.true.,2*npack0,VQ)
         do ii=1,nion

*           **** structure factor  ****
            call strfac_pack(0,ii,dcpl_mb(exi(1)))

            do iii=1,nga
               i = iii + (ii-1)*nga

               afac = omega*pspw_Vfac_APC(ngs,dbl_mb(Am(1)),u,i)
               !afac = dsqrt(omega)*pspw_Vfac_APC(ngs,dbl_mb(Am(1)),u,i)
               !afac = pspw_Vfac_APC(ngs,dbl_mb(Am(1)),u,i)

*              *** gaus_i(G))***
               call Pack_tc_Mul(0,dbl_mb(gaus(1)+npack0*(iii-1)),
     >                       dcpl_mb(exi(1)),
     >                       dcpl_mb(gaus_i(1)))

*              *** VQ(G) += afac*w(G)*gaus_i(G))***
               call Pack_tc_aMulAdd(0,afac,
     >                              dbl_mb(w(1)),
     >                              dcpl_mb(gaus_i(1)),
     >                              VQ)
            end do
         end do
         !call Pack_c_addzero(0,omega*sumAmU/sumAm,VQ)
         call Pack_c_addzero(0,sumAmU/sumAm,VQ)

         value =           BA_pop_stack(gaus_i(2))
         value = value.and.BA_pop_stack(exi(2))
         if (.not. value) 
     >      call errquit('pspw_VQ_APC:pop stack',0,MA_ERR)


      end if
      return
      end


*     ***********************************
*     *					*
*     *		pspw_cosmo_V0_APC	*
*     *					*
*     ***********************************
*
      subroutine pspw_cosmo_V0_APC(vcosmo)
      implicit none
      complex*16 vcosmo(*)

#include "bafdecls.fh"
#include "pspw_APC.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value
      integer taskid,np,tid,nthreads
      integer j,ii,k,nion,nion_q
      real*8 sb,qj,x,y,z,r
      integer u(2),qcoord_ptr

*     **** external functions ****
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads
      integer  ion_nion_qm,nwpw_cosmo_nq
      external ion_nion_qm,nwpw_cosmo_nq
      integer  nwpw_cosmo_qcoord_ptr
      external nwpw_cosmo_qcoord_ptr
      real*8   nwpw_cosmo_screen_qc0,ion_rion
      external nwpw_cosmo_screen_qc0,ion_rion

      call Parallel_taskid(taskid)
      call Parallel_np(np)
      tid = Parallel_threadid()
      nthreads = Parallel_nthreads()
      nion     = ion_nion_qm()
      nion_q   = nwpw_cosmo_nq()
      qcoord_ptr = nwpw_cosmo_qcoord_ptr()

*     **** get memory from stack ****
      value = BA_push_get(mt_dbl,ngs,'u',u(2),u(1))
      if (.not.value) 
     >   call errquit('pspw_cosmo_V0_APC:push stack',0, MA_ERR)

*     **** calculate u = dEQelcq/dq ****
      call Parallel_shared_vector_zero(.true.,ngs,dbl_mb(u(1)))
      do j=tid+1,nion_q,nthreads
         qj = nwpw_cosmo_screen_qc0(j)
         do ii=taskid+1,nion,np
            x = dbl_mb(qcoord_ptr+3*(j-1))  -ion_rion(1,ii)
            y = dbl_mb(qcoord_ptr+3*(j-1)+1)-ion_rion(2,ii)
            z = dbl_mb(qcoord_ptr+3*(j-1)+2)-ion_rion(3,ii)
            r = dsqrt(x*x + y*y + z*z)
            sb = -qj/r
            do k=1,nga
                 dbl_mb(u(1)+(ii-1)*nga+k-1)
     >         = dbl_mb(u(1)+(ii-1)*nga+k-1) + sb
            end do
         end do
      end do
      call Parallel_Vector_SumAll(ngs,dbl_mb(u(1)))

      call pspw_VQ_APC(nion,dbl_mb(u(1)),vcosmo)

      value = BA_pop_stack(u(2))
      if (.not.value)
     >   call errquit('pspw_cosmo_V0_APC:pop stack',0, MA_ERR)

      return
      end 


*     ***********************************
*     *					*
*     *		pspw_cosmo_V_APC	*
*     *					*
*     ***********************************
*
      subroutine pspw_cosmo_V_APC(ispin,ne,dng,vcosmo,ecosmo,pcosmo)
      implicit none
      integer ispin,ne(2)
      complex*16 dng(*), vcosmo(*)
      real*8    ecosmo,pcosmo

#include "bafdecls.fh"
#include "pspw_APC.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value
      integer nthreads,tid,np,taskid
      integer j,ii,ia,k
      integer nion_q,nion,npack0
      real*8 x,y,z,r,sb,qj
      real*8 suma,sumb

      integer BQ_ptr,u(2),qcoord_ptr,vtmp(2)
      
*     **** external functions ****
      integer  Parallel_nthreads,Parallel_threadid
      external Parallel_nthreads,Parallel_threadid
      integer  nwpw_cosmo_BQ_ptr,nwpw_cosmo_nq,ion_nion_qm,ion_katm
      external nwpw_cosmo_BQ_ptr,nwpw_cosmo_nq,ion_nion_qm,ion_katm
      real*8   pspw_getQtot_APC,ion_rion,psp_zv,pspw_getQ_APC
      external pspw_getQtot_APC,ion_rion,psp_zv,pspw_getQ_APC
      integer  nwpw_cosmo_qcoord_ptr
      external nwpw_cosmo_qcoord_ptr
      real*8   nwpw_cosmo_screen_qc0
      external nwpw_cosmo_screen_qc0
      real*8   nwpw_cosmo_EQelcq
      external nwpw_cosmo_EQelcq

      call Parallel_taskid(taskid)
      call Parallel_np(np)
      tid = Parallel_threadid()
      nthreads = Parallel_nthreads()
      nion     = ion_nion_qm()
      nion_q   = nwpw_cosmo_nq()
      call Pack_npack(0,npack0)

      value = BA_push_get(mt_dcpl,npack0,'vtmp',vtmp(2),vtmp(1))
      if (.not.value)
     >   call errquit('pspw_cosmo_V_APC:push stack',0,MA_ERR)

      qcoord_ptr = nwpw_cosmo_qcoord_ptr()
      BQ_ptr     = nwpw_cosmo_BQ_ptr()
      call Parallel_shared_vector_zero(.true.,nion_q,dbl_mb(BQ_ptr))
      call Parallel_shared_vector_zero(.true.,2*npack0,dcpl_mb(vtmp(1)))
      


*     **** generate APC charges ****
      call pspw_gen_APC(ispin,ne,dng)

      !**** q-Qion + q-Qelc BQ contributions ****
      do j=tid+1,nion_q,nthreads
         do ii=taskid+1,nion,np
            ia=ion_katm(ii)
            x = dbl_mb(qcoord_ptr+3*(j-1))  -ion_rion(1,ii)
            y = dbl_mb(qcoord_ptr+3*(j-1)+1)-ion_rion(2,ii)
            z = dbl_mb(qcoord_ptr+3*(j-1)+2)-ion_rion(3,ii)
            r = dsqrt(x*x + y*y +z*z)
            sb = (psp_zv(ia)+pspw_getQtot_APC(ii))/r
            dbl_mb(BQ_ptr+j-1)= dbl_mb(BQ_ptr+j-1) + sb
         end do
      end do
      call Parallel_Vector_SumAll(nion_q,dbl_mb(BQ_ptr))

*     **** solve for cosmo q charges ****
      call nwpw_cosmo_solve_q()


*     **** generate cosmo potential using current APC charges ****
*     **** generate APC potential ****
      call pspw_cosmo_V0_APC(dcpl_mb(vtmp(1)))
      call Pack_cc_daxpy(0,1.0d0,dcpl_mb(vtmp(1)),vcosmo)

      call Pack_cc_dot(0,dng,dcpl_mb(vtmp(1)),pcosmo)
      ecosmo = nwpw_cosmo_EQelcq()

      value = value.and.BA_pop_stack(vtmp(2))
      if (.not.value)
     >   call errquit('pspw_cosmo_V_APC:pop stack',0,MA_ERR)
      return
      end


*     ***********************************
*     *					*
*     *		pspw_born_V0_APC	*
*     *					*
*     ***********************************
*
      subroutine pspw_born_V0_APC(vborn)
      implicit none
      complex*16 vborn(*)

#include "bafdecls.fh"
#include "pspw_APC.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value
      integer taskid,np,tid,nthreads
      integer j,ii,k,nion,nion_q
      real*8 sb,qj,x,y,z,r
      integer u(2),uion(2),qion(2)

*     **** external functions ****
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads
      integer  ion_nion
      external ion_nion
      real*8   pspw_getQtot_APC
      external pspw_getQtot_APC

      call Parallel_taskid(taskid)
      call Parallel_np(np)
      tid = Parallel_threadid()
      nthreads = Parallel_nthreads()
      nion     = ion_nion()

*     **** get memory from stack ****
      value = BA_push_get(mt_dbl,nion,'qion',qion(2),qion(1))
      value = value.and.
     >        BA_push_get(mt_dbl,nion,'uion',uion(2),uion(1))
      value = value.and.
     >        BA_push_get(mt_dbl,ngs,'u',u(2),u(1))
      if (.not.value) 
     >   call errquit('pspw_born_V_APC:push stack',0, MA_ERR)

      do ii=tid+1,nion,nthreads
         dbl_mb(qion(1)+ii-1) = pspw_getQtot_APC(ii)
      end do

*     **** calculate u = dEQelcq/dq ****
      call nwpw_born_dVdq(nion,dbl_mb(qion(1)),
     >                         dbl_mb(uion(1)))
      call Parallel_shared_vector_zero(.true.,ngs,dbl_mb(u(1)))
      do ii=tid+1,nion,nthreads
         do k=1,nga
              dbl_mb(u(1)+(ii-1)*nga+k-1)
     >      = dbl_mb(u(1)+(ii-1)*nga+k-1) + dbl_mb(uion(1)+ii-1)
         end do
      end do

      call pspw_VQ_APC(nion,dbl_mb(u(1)),vborn)

      value = BA_pop_stack(u(2))
      value = value.and.BA_pop_stack(uion(2))
      value = value.and.BA_pop_stack(qion(2))
      if (.not.value)
     >   call errquit('pspw_born_V0_APC:pop stack',0, MA_ERR)

      return
      end 


*     ***********************************
*     *                                 *
*     *         pspw_born_V_APC         *
*     *                                 *
*     ***********************************
*
      subroutine pspw_born_V_APC(ispin,ne,dng,vborn)
      implicit none
      integer ispin,ne(2)
      complex*16 dng(*), vborn(*)

#include "bafdecls.fh"
#include "pspw_APC.fh"
#include "errquit.fh"

*     **** generate APC charges ****
      call pspw_gen_APC(ispin,ne,dng)

*     **** generate APC potential ****
      call pspw_born_V0_APC(vborn)
      return
      end 


*     ***********************************
*     *                                 *
*     *         pspw_cdft_V_APC         *
*     *                                 *
*     ***********************************
*
      subroutine pspw_cdft_V_APC(ispin,ne,dng,vcdft,ecdft,pcdft)
      implicit none
      integer ispin,ne(2)
      complex*16 dng(*),vcdft(*)
      real*8     ecdft,pcdft

#include "bafdecls.fh"
#include "pspw_APC.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value
      integer taskid,np,tid,nthreads
      integer nion,i,k,ii
      integer u(2),vtmp(2),npack0
      real*8 elocal,elocal1


*     **** external functions ****
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads
      integer  nwpw_cdft_u_ptr,nwpw_cdft_q_ptr
      external nwpw_cdft_u_ptr,nwpw_cdft_q_ptr
      integer  ion_nion,ion_katm_qm
      external ion_nion,ion_katm_qm
      real*8   pspw_getQtot_APC,psp_zv,nwpw_cdft_energy
      external pspw_getQtot_APC,psp_zv,nwpw_cdft_energy
      real*8   lattice_omega
      external lattice_omega

      call Parallel_taskid(taskid)
      call Parallel_np(np)
      tid = Parallel_threadid()
      nthreads = Parallel_nthreads()
      nion     = ion_nion()
      call Pack_npack(0,npack0)


*     **** get memory from stack ****
      value = BA_push_get(mt_dcpl,npack0,'vtmp',vtmp(2),vtmp(1))
      value = BA_push_get(mt_dbl,ngs,'u',u(2),u(1))
      if (.not.value)
     >   call errquit('pspw_cdft_V_APC:push stack',0, MA_ERR)

      call Parallel_shared_vector_zero(.true.,ngs,dbl_mb(u(1)))
      do ii=tid+1,nion,nthreads
         do k=1,nga
            dbl_mb(u(1)+(ii-1)*nga+k-1) = dbl_mb(u(1)+(ii-1)*nga+k-1)
     >                                  - dbl_mb(nwpw_cdft_u_ptr()+ii-1)
         end do
      end do

*     **** generate APC charges ****
      call pspw_gen_APC(ispin,ne,dng)
      do ii=1,nion
         dbl_mb(nwpw_cdft_q_ptr()+ii-1) = pspw_getQtot_APC(ii) 
     >                                  + psp_zv(ion_katm_qm(ii))
      end do


      call Pack_cc_dot(0,dng,vcdft,elocal)

*     **** generate APC potential ****
      call Parallel_shared_vector_zero(.true.,2*npack0,dcpl_mb(vtmp(1)))
      call pspw_VQ_APC(nion,dbl_mb(u(1)),dcpl_mb(vtmp(1)))
      call Pack_cc_daxpy(0,1.0d0,dcpl_mb(vtmp(1)),vcdft)

      call Pack_cc_dot(0,dng,dcpl_mb(vtmp(1)),pcdft)
      ecdft = nwpw_cdft_energy()


      !call Pack_cc_dot(0,dng,vcdft,elocal1)
c      write(*,*) "ELOCAL           =", elocal
c      write(*,*) "PLOCAL+Pcdft     =", elocal1
c      write(*,*) "dng*vcdft, diff  =", pcdft, (elocal1-elocal)
c      write(*,*) "ECDFT (point Q)  = ",ecdft
c      write(*,*) "      (kcal/mol) = ",ecdft*27.2116d0*23.06d0
c      do ii=1,nion
c         write(*,*) " - ii,Q=",ii, dbl_mb(nwpw_cdft_q_ptr()+ii-1) 
c      end do
c      write(*,*)


*     **** pop stack ****
      value = BA_pop_stack(u(2))
      value = value.and.BA_pop_stack(vtmp(2))
      if (.not.value)
     >   call errquit('pspw_cdft_V_APC:pop stack',0, MA_ERR)

      return
      end

      ! - calculate dbi/dR_I
      ! - calculate dAij/dR_I



*     ***********************************
*     *					*
*     *		pspw_V_APC_on	        *
*     *					*
*     ***********************************
      logical function pspw_V_APC_on()
      implicit none

*     **** external functions ****
      logical  nwpw_cosmo2_on,nwpw_born_on,nwpw_cdft_on
      external nwpw_cosmo2_on,nwpw_born_on,nwpw_cdft_on

      pspw_V_APC_on =  nwpw_cosmo2_on().or.
     >                 nwpw_born_on().or.
     >                 nwpw_cdft_on()
      return
      end


*     ***********************************
*     *					*
*     *		pspw_V_APC	        *
*     *					*
*     ***********************************
      subroutine pspw_V_APC(ispin,ne,dng,vapc,Eapc,Papc)
      implicit none
      integer ispin,ne(2)
      complex*16 dng(*), vapc(*)
      real*8 Eapc,Papc

      real*8 Eapc_APC,Papc_APC
      common /e_apc_common/ Eapc_APC,Papc_APC
      
*     **** external functions ****
      logical  nwpw_cosmo2_on,nwpw_born_on,nwpw_cdft_on
      external nwpw_cosmo2_on,nwpw_born_on,nwpw_cdft_on

      Eapc = 0.0d0
      Papc = 0.0d0
      if (nwpw_cosmo2_on()) then
         call pspw_cosmo_V_APC(ispin,ne,dng,vapc,Eapc,Papc)
      else if (nwpw_born_on()) then
         call pspw_born_V_APC(ispin,ne,dng,vapc)
      else if (nwpw_cdft_on()) then
         call pspw_cdft_V_APC(ispin,ne,dng,vapc,Eapc,Papc)
      end if
      Eapc_APC = Eapc
      Papc_APC = Papc

      return
      end


*     ***********************************
*     *					*
*     *		pspw_E_APC	        *
*     *					*
*     ***********************************
      subroutine pspw_E_APC(Eapc,Papc)
      implicit none
      real*8 Eapc,Papc

      real*8 Eapc_APC,Papc_APC
      common /e_apc_common/ Eapc_APC,Papc_APC

      Eapc = Eapc_APC
      Papc = Papc_APC
      return
      end







*     ***********************************
*     *					*
*     *		pspw_getQ_APC		*
*     *					*
*     ***********************************
      real*8 function pspw_getQ_APC(ii,n)
      implicit none
      integer ii,n

#include "bafdecls.fh"
#include "pspw_APC.fh"

*     **** local variables ****
      integer i
      real*8 qq

      qq = 0.0d0

      if (nga.gt.0) then
         i = n + (ii-1)*nga
         qq = dbl_mb(q(1)+i-1)
      end if
      pspw_getQ_APC = qq

      return
      end

*     ***********************************
*     *					*
*     *		pspw_getQtot_APC	*
*     *					*
*     ***********************************
      real*8 function pspw_getQtot_APC(ii)
      implicit none
      integer ii

#include "bafdecls.fh"
#include "pspw_APC.fh"

*     **** local variables ****
      integer i,n
      real*8 qq

      qq = 0.0d0

      if (nga.gt.0) then
         do n=1,nga
            i = n + (ii-1)*nga
            qq = qq - dbl_mb(q(1)+i-1)
         end do
      end if
      pspw_getQtot_APC = qq
      return
      end



*     ***********************************
*     *                                 *
*     *         pspw_shortprint_APC     *
*     *                                 *
*     ***********************************
      subroutine pspw_shortprint_APC(unit)
      implicit none
      integer unit

#include "bafdecls.fh"
#include "pspw_APC.fh"
#include "util.fh"

      integer  ion_nion,nwpw_cdft_u_ptr,nwpw_cdft_q_ptr
      external ion_nion,nwpw_cdft_u_ptr,nwpw_cdft_q_ptr

*     *** local variables ***
      integer ii

      write(unit,*)
      write(unit,*) "APC Potential:"
      write(unit,'(10F14.9)') (dbl_mb(nwpw_cdft_u_ptr()+ii-1),
     >                         ii=1,ion_nion())
      write(unit,*)
      write(unit,*) "APC Point Charges:"
      write(unit,'(10F14.9)') (dbl_mb(nwpw_cdft_q_ptr()+ii-1),
     >                         ii=1,ion_nion())

      return
      end



*     ***********************************
*     *					*
*     *		pspw_print_APC		*
*     *					*
*     ***********************************
      subroutine pspw_print_APC(unit)
      implicit none
      integer unit

#include "bafdecls.fh"
#include "pspw_APC.fh"
#include "util.fh"

*     **** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)
      integer i,j,indx
      real*8 sum,sumi,sume
      logical oprint

*     **** external functions ****
      character*4  ion_atom
      integer      ion_nion_qm,ion_katm_qm
      real*8       psp_zv
      external     ion_atom
      external     ion_nion_qm,ion_katm_qm
      external     psp_zv
      logical      control_print
      external     control_print
      
      if (nga.gt.0) then
      call Parallel_taskid(taskid)
      oprint = (taskid.eq.MASTER).and.control_print(print_medium)

      if (oprint) then

         WRITE(unit,*)
         WRITE(unit,*)
     >  '*************************************************************'
         WRITE(unit,*)
     >  '**                                                         **'
         WRITE(unit,*)
     >  '**          PSPW Atomic Point Charge (APC) Analysis        **'
         WRITE(unit,*)
     >  '**                                                         **'
         WRITE(unit,*)
     >  '**   Point charge analysis based on paper by P.E. Blochl   **'
         WRITE(unit,*)
     >  '**         (J. Chem. Phys. vol 103, page 7422, 1995)       **'
         WRITE(unit,*)
     >  '**                                                         **'
         WRITE(unit,*)
     >  '*************************************************************'


*        **** write out information *****
         write(unit,*) 
         write(unit,*) 'pspw_APC data structure'
         write(unit,*) '-----------------------'
         write(unit,*) 'nga, ngs:',nga,ngs
         write(unit,*) 'Gc      :',Gc
         do i=1,nga
           write(unit,*) 'APC gamma:',i, dbl_mb(gamma(1)+i-1)
         end do

         write(unit,*) 
         WRITE(unit,*) 'charge analysis on each atom'
         WRITE(unit,*) '----------------------------'
         sume=0.0d0
         sumi=0.0d0
         write(unit,100) 'no','atom','Qelc','Qion','Qtotal'
         write(unit,105) 
         do j=1,ion_nion_qm()
           sum = 0.0d0
           do i=1,nga
              indx = (i-1) + (j-1)*nga
              sum = sum + dbl_mb(q(1)+indx)
           end do
           sume = sume - sum
           sumi = sumi + psp_zv(ion_katm_qm(j))
           write(unit,110) j,ion_atom(ion_katm_qm(j)),
     >                  (-sum),
     >                  (psp_zv(ion_katm_qm(j))),
     >                  (psp_zv(ion_katm_qm(j))-sum)
         end do
         write(unit,120) sume,sumi,(sume+sumi)

         write(unit,*) 
         write(unit,*) 
         WRITE(unit,*) 'gaussian coefficients of model density'
         WRITE(unit,*) '--------------------------------------'
         write(unit,200) 'no','atom',
     >                   'g=',0.0d0,
     >                    ('g=',dbl_mb(gamma(1)+i-1),i=1,nga)
         write(unit,205) ('-------',i=0,nga)
         do j=1,ion_nion_qm()
           write(unit,210) j,ion_atom(ion_katm_qm(j)),
     >                     psp_zv(ion_katm_qm(j)),
     >                     (-dbl_mb(q(1)+i-1+(j-1)*nga),i=1,nga)
         end do


      end if

      call nwpw_born_print(nga,ion_nion_qm(),dbl_mb(q(1)))

      end if
      return
  100 format(/2x,A4,A6,3A12)
  105 format(4x,'--',2x,'----',
     >       5x,'-------',
     >       5x,'-------',
     >       5x,'-------')
  110 format(2x,I4,A6,3F12.3)
  120 format(2x,'   Total Q',3F12.3)
  200 format(/2x,A4,A6,20(5x,A2,F5.3))
  205 format(4x,'--',2x,'----',
     >       20(5x,A7))
  210 format(2x,I4,A6,20F12.3)
      end

