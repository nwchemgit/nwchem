c
c $Id$
c


********************************************************************
*
*	band_APC module: used to generate derived atomic point charges
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
*     *		band_init_APC		*
*     *					*
*     ***********************************
      subroutine band_init_APC(rtdb)
      implicit none
      integer rtdb

#include "bafdecls.fh"
#include "btdb.fh"
#include "errquit.fh"
#include "band_APC.fh"

*     **** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)

      logical value,doapc
      integer i,j,k,n,npack0,nfft3d
      !integer nx,ny,nxh
      integer zero,pzero,qzero
      integer tmp1(2),G(3)
      real*8  coef,x,gg,fourpi

*     **** external functions ****
      logical  control_Mulliken
      integer  ion_nion_qm, c_G_indx
      real*8   lattice_omega
      external control_Mulliken
      external ion_nion_qm, c_G_indx
      external lattice_omega

      if (.not.btdb_get(rtdb,'nwpw:APC',mt_log,1,doapc))
     >  doapc = control_Mulliken()

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
      call Cram_npack(0,npack0)


*     **** Allocate memory needed for nwpw_APC ****
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
     >   call errquit('nwpw_init_APC:out of heap memory',0, MA_ERR)

*     **** read in Gc and gamma from rtdb ***
      if (.not.btdb_get(rtdb,'nwpw_APC:Gc',mt_dbl,1,Gc)) Gc = 2.5d0
      if (.not.btdb_get(rtdb,'nwpw_APC:gamma',mt_dbl,nga,
     >   dbl_mb(gamma(1)))) then
         if (nga.eq.3) then
            dbl_mb(gamma(1))   = 0.6d0
            dbl_mb(gamma(1)+1) = 0.9d0
            dbl_mb(gamma(1)+2) = 1.35d0
         else
            call errquit('nwpw:init_APC:error reading rtdb',0,RTDB_ERR)
         end if
      end if


*     ********************************
*     **** define weight function ****
*     ********************************

      fourpi = 4.0d0*(4.0d0*datan(1.0d0))
      call C3dB_nfft3d(1,nfft3d)
      G(1)= c_G_indx(1)
      G(2)= c_G_indx(2)
      G(3)= c_G_indx(3)

*     ***** find the G==0 point in the lattice *****
      i=0
      j=0
      k=0
      !call C3dB_nx(1,nx)
      !call C3dB_ny(1,ny)
      !nxh = nx/2
c     call D3dB_ktoqp(1,k+1,qzero,pzero)
c     zero = (qzero-1)*(nxh+1)*ny
c    >     + j*(nxh+1)
c    >     + i+1
      call C3dB_ijktoindexp(1,i+1,j+1,k+1,zero,pzero)

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
      call Cram_r_pack(0,dbl_mb(tmp1(1)))
      call Cram_r_Copy(0,dbl_mb(tmp1(1)),dbl_mb(w(1)))
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
         call Cram_r_pack(0,dbl_mb(tmp1(1)))
         call Cram_r_Copy(0,dbl_mb(tmp1(1)),
     >                      dbl_mb(gaus(1)+(n-1)*npack0))
      end do

      value = BA_pop_stack(tmp1(2))
      if (.not. value) call errquit('popping of stack memory',0, MA_ERR)

*     **** DEBUG - write out information *****
      if (taskid.eq.MASTER) then
         write(*,*) 
         write(*,*) 'initializing band_APC data structure'
         write(*,*) '------------------------------------'
         write(*,*) 'nga, ngs:',nga,ngs
         write(*,*) 'Gc      :',Gc
         do i=1,nga
           write(*,*) 'APC gamma:',i, dbl_mb(gamma(1)+i-1)
         end do
      end if

      return
      end

*     ***********************************
*     *					*
*     *		band_end_APC		*
*     *					*
*     ***********************************
      subroutine band_end_APC()
      implicit none
#include "errquit.fh"

#include "bafdecls.fh"
#include "band_APC.fh"

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



*     ***********************************
*     *					*
*     *		band_gen_APC		*
*     *					*
*     ***********************************
      subroutine band_gen_APC(ispin,ne,dng)
      implicit none
      integer ispin,ne(2)
      complex*16 dng(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "band_APC.fh"

*     ***** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)

      logical value
      integer i,j,ii,iii,jj,jjj,indx,indxt
      integer nfft3d,npack0
      real*8 omega,N,sum,sum1

      integer exi(2),exj(2)
      integer gaus_i(2),gaus_j(2),ipiv(2)
      integer AAA(2),work(2),lwork,rank,info

*     ***** external functions ****
      integer  ion_nion_qm
      real*8   lattice_omega
      external ion_nion_qm
      external lattice_omega

      if (nga.gt.0) then
         omega = lattice_omega()
         call Parallel_taskid(taskid)
         call C3dB_nfft3d(1,nfft3d)
         call Cram_npack(0,npack0)

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
            call cstrfac(ii,dcpl_mb(exi(1)))
            call Cram_c_pack(0,dcpl_mb(exi(1)))

            do iii=1,nga
               i = iii + (ii-1)*nga

*              *** gaus_i(G))***
               call Cram_rc_Mul(0,dbl_mb(gaus(1)+npack0*(iii-1)),
     >                          dcpl_mb(exi(1)),
     >                          dcpl_mb(gaus_i(1)))

*              *** w(G)*gaus_i(G))***
c               call Pack_tc_Mul(0,dbl_mb(w(1)),
c     >                          dcpl_mb(gaus_i(1)),
c     >                          dcpl_mb(gaus_i(1)))
               call Cram_rc_Mul2(0,dbl_mb(w(1)),dcpl_mb(gaus_i(1)))

*              *** omega*Sum(G) w(G)*Re(dcongj(dng(G))*gaus_i(G))***
               call Cram_cc_dot(0,dng,
     >                            dcpl_mb(gaus_i(1)),
     >                            sum)
               dbl_mb(b(1)+i-1) = sum*omega
            end do
          end do



*        ***** calculate A_ij = omega*Sum(G) w(G)*dcongj(gaus_i(G))*gaus_j(G)) ****
         do ii=1,ion_nion_qm()
*           **** structure factor  ****
            call cstrfac(ii,dcpl_mb(exi(1)))
            call Cram_c_pack(0,dcpl_mb(exi(1)))

             do jj=ii,ion_nion_qm()
*              **** structure factor  ****
               call cstrfac(jj,dcpl_mb(exj(1)))
               call Cram_c_pack(0,dcpl_mb(exj(1)))


               do iii=1,nga
*                 *** gaus_i(G))***
                  call Cram_rc_Mul(0,dbl_mb(gaus(1)+npack0*(iii-1)),
     >                               dcpl_mb(exi(1)),
     >                               dcpl_mb(gaus_i(1)))

*                 *** w(G)*gaus_i(G))***
c                  call Pack_tc_Mul(0,dbl_mb(w(1)),
c     >                               dcpl_mb(gaus_i(1)),
c     >                               dcpl_mb(gaus_i(1)))
                  call Cram_rc_Mul2(0,dbl_mb(w(1)),dcpl_mb(gaus_i(1)))

                  do jjj=1,nga
*                   *** gaus_j(G))***
                     call Cram_rc_Mul(0,dbl_mb(gaus(1)+npack0*(jjj-1)),
     >                                  dcpl_mb(exj(1)),
     >                                  dcpl_mb(gaus_j(1)))

*                   *** omega*Sum(G) w(G)*gaus_i(G)*gaus_j(G))***
                    call Cram_cc_dot(0,dcpl_mb(gaus_i(1)),
     >                                 dcpl_mb(gaus_j(1)),
     >                                 sum)
                    sum = sum*omega
                    i = iii + (ii-1)*nga
                    j = jjj + (jj-1)*nga

                    indx  = (i-1) + (j-1)*ngs
                    indxt = (j-1) + (i-1)*ngs
                    dbl_mb(A(1)+indx)  = sum
                    dbl_mb(A(1)+indxt) = sum
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
         call dcopy(ngs*ngs,0.0d0,0,dbl_mb(Am(1)),1)
         call dcopy(ngs,0.0d0,0,dbl_mb(q(1)),1)
         if (taskid.eq.MASTER) then

           lwork = 5*ngs*ngs
           !value = BA_push_get(mt_int,ngs,'ipivAPC',ipiv(2),ipiv(1))
           value = BA_push_get(mt_dbl,lwork,'work',work(2),work(1))
           value = value.and.
     >             BA_push_get(mt_dbl,ngs*ngs,'AAAAPC',AAA(2),AAA(1))
           if (.not. value) call errquit('out of stack memory',0,
     &       MA_ERR)

*          **** calculate Am_ij ****
           call dcopy(ngs*ngs,dbl_mb(A(1)),1,dbl_mb(AAA(1)),1)
           call dcopy(ngs*ngs,0.0d0,0,dbl_mb(Am(1)),1)
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
            value = BA_pop_stack(AAA(2))
            value = value.and.BA_pop_stack(work(2))
            !value = value.and.BA_pop_stack(ipiv(2))
            if (.not. value) call errquit('popping stack memory',0,
     &       MA_ERR)

           call dcopy(ngs,0.0d0,0,dbl_mb(q(1)),1)
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
              dbl_mb(q(1)+i-1) = sum1
           end do
         end if
*        **** synchronization ****
         call C3dB_Vector_SumAll(ngs,dbl_mb(q(1)))
         call C3dB_Vector_SumAll(ngs*ngs,dbl_mb(Am(1)))

      end if

      return
      end

*     ***********************************
*     *					*
*     *		band_dngen_APC		*
*     *					*
*     ***********************************
      subroutine band_dngen_APC(ispin,ne,dn)
      implicit none
      integer ispin,ne(2)
      real*8 dn(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "band_APC.fh"

*     ***** local variables ****
      logical value
      integer nx,ny,nz,n2ft3d,nfft3d,npack0
      real*8 scal1

      integer dng(2),tmp(2)

      if (nga.gt.0) then

      call C3dB_nx(1,nx)
      call C3dB_ny(1,ny)
      call C3dB_nz(1,nz)
      scal1 = 1.0d0/dble(nx*ny*nz)

      call Cram_npack(0,npack0)
      call C3dB_nfft3d(1,nfft3d)
      n2ft3d = 2*nfft3d

*     **** get memory from stack ****
      value = BA_push_get(mt_dcpl,npack0,'dngakjs',dng(2),dng(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,nfft3d,'tmpxkjz',tmp(2),tmp(1))
      if (.not. value) call errquit('out of stack memory',0,
     &       MA_ERR)

*     **** generate dng ****
      call C3dB_rr_Sum(1,dn(1),dn(1+(ispin-1)*n2ft3d),dcpl_mb(tmp(1)))
c      call D3dB_r_SMul(1,scal1,dcpl_mb(tmp(1)),dcpl_mb(tmp(1)))
      call C3dB_r_SMul1(1,scal1,dcpl_mb(tmp(1)))
      call C3dB_rc_fft3f(1,dcpl_mb(tmp(1)))
      call Cram_c_pack(0,dcpl_mb(tmp(1)))
      call Cram_c_Copy(0,dcpl_mb(tmp(1)),dcpl_mb(dng(1)))
      value = BA_pop_stack(tmp(2))
      if (.not. value) call errquit('popping of stack memory',0,
     &       MA_ERR)

*     **** generate APC *****
      call band_gen_APC(ispin,ne,dcpl_mb(dng(1)))


      value = value.and.BA_pop_stack(dng(2))
      if (.not. value) call errquit('popping of stack memory',0,
     &       MA_ERR)

      end if

      return
      end

*     ***********************************
*     *					*
*     *		band_getQ_APC		*
*     *					*
*     ***********************************
      real*8 function band_getQ_APC(ii,n)
      implicit none
      integer ii,n

#include "bafdecls.fh"
#include "band_APC.fh"

*     **** local variables ****
      integer i
      real*8 qq

      qq = 0.0d0

      if (nga.gt.0) then
         i = n + (ii-1)*nga
         qq = dbl_mb(q(1)+i-1)
      end if
      band_getQ_APC = qq

      return
      end



*     ***********************************
*     *					*
*     *		band_print_APC		*
*     *					*
*     ***********************************
      subroutine band_print_APC(unit)
      implicit none
      integer unit

#include "bafdecls.fh"
#include "band_APC.fh"

*     **** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)
      integer i,j,indx
      real*8 sum,sumi,sume

*     **** external functions ****
      character*4  ion_atom
      integer      ion_nion_qm,ion_katm_qm
      real*8       cpsp_zv
      external     ion_atom
      external     ion_nion_qm,ion_katm_qm
      external     cpsp_zv
      
      if (nga.gt.0) then
      call Parallel_taskid(taskid)

      if (taskid.eq.MASTER) then

         WRITE(unit,*)
         WRITE(unit,*)
     >  '*************************************************************'
         WRITE(unit,*)
     >  '**                                                         **'
         WRITE(unit,*)
     >  '**          BAND Atomic Point Charge (APC) Analysis        **'
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
         write(unit,*) 'band_APC data structure'
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
           sumi = sumi + cpsp_zv(ion_katm_qm(j))
           write(unit,110) j,ion_atom(ion_katm_qm(j)),
     >                  (-sum),
     >                  (cpsp_zv(ion_katm_qm(j))),
     >                  (cpsp_zv(ion_katm_qm(j))-sum)
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
     >                     cpsp_zv(ion_katm_qm(j)),
     >                     (-dbl_mb(q(1)+i-1+(j-1)*nga),i=1,nga)
         end do


      end if

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

