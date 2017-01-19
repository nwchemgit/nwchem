
*     ************************************************
*     *                                              *
*     *        pspw_hfx_localize_init                *
*     *                                              *
*     ************************************************
      subroutine pspw_hfx_localize_init()
      implicit none

#include "bafdecls.fh"
#include "pspw_hfx.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value
      integer neq(2),ne(2)
      logical  Dneall_m_allocate
      external Dneall_m_allocate

      if (localize_on) then
         call nwpw_timing_start(52)
         call Dneall_neq(neq)
         ne(1) = norbs(1)
         ne(2) = norbs(2)
         findcolumns = .true.

         value = BA_alloc_get(mt_dbl,3*(ne(1)+ne(2)),
     >                        'rwav',rwav(2),rwav(1))
         value = value.and.
     >           BA_alloc_get(mt_int,2*(ne(1)+ne(2)),
     >                        'dmindx',dmindx(2),dmindx(1))
         value = value.and.
     >           BA_alloc_get(mt_dbl,n2ft3d*(neq(1)+neq(2)),
     >                        'psiloc_r',psiloc_r(2),psiloc_r(1))
         value = BA_alloc_get(mt_dbl,(ne(1)*ne(1)+ne(2)*ne(2)),
     >                'switchr',switchr(2),switchr(1))
         value = value.and.Dneall_m_allocate(0,rotmat21)
         value = value.and.Dneall_m_allocate(0,absmat)
         value = BA_alloc_get(mt_log,(ne(1)*ne(1)+ne(2)*ne(2)),
     >                'closeenough',closeenough(2),closeenough(1))
         value = BA_alloc_get(mt_int,3*(ne(1)*ne(1)+ne(2)*ne(2)),
     >                'centeroverlap',centeroverlap(2),centeroverlap(1))
         value = BA_alloc_get(mt_dbl,3*(ne(1)*ne(1)+ne(2)*ne(2)),
     >                'roverlap',roverlap(2),roverlap(1))
         if (.not.value)
     >      call errquit('pspw_hfx_localize_init: out of heap',3,MA_ERR)
         call nwpw_timing_end(52)

         if (localize2_on) call pspw_hfx_localize2_init()
      end if

      return
      end 


*     ************************************************
*     *                                              *
*     *        pspw_hfx_localize2_init               *
*     *                                              *
*     ************************************************
      subroutine pspw_hfx_localize2_init()
      implicit none

#include "bafdecls.fh"
#include "pspw_hfx.fh"
#include "errquit.fh"

*     **** local variables ****
      integer MASTER,taskid,tmp_len
      parameter (MASTER=0,tmp_len=140)

      logical oprint
      integer i,j,info
      real*8 wrk(6,6),bmat(3,3)
      real*8 b(3),ixmat(3,6)
      real*8 scal,xx,yy,zz,tmp1(tmp_len)


*     **** external functions ****
      real*8   lattice_omega,lattice_unitg
      external lattice_omega,lattice_unitg

      call Parallel_taskid(taskid)
      oprint = (taskid.eq.MASTER)
      oprint = .false.

      call Pack_npack(1,npack1)


c     *** Silvestrelli G1 ***
      ixmat(1,1)=1.0d0
      ixmat(2,1)=0.0d0
      ixmat(3,1)=0.0d0

c     *** Silvestrelli G4 ***
      ixmat(1,2)=1.0d0
      ixmat(2,2)=1.0d0
c      ixmat(3,2)=0.0d0

c     *** Silvestrelli G5 ***
      ixmat(1,3)=1.0d0
      ixmat(2,3)=0.0d0
      ixmat(3,3)=1.0d0

c     *** Silvestrelli G2 ***
      ixmat(1,4)=0.0d0
      ixmat(2,4)=1.0d0
      ixmat(3,4)=0.0d0

c     *** Silvestrelli G6 ***
      ixmat(1,5)=0.0d0
      ixmat(2,5)=1.0d0
      ixmat(3,5)=1.0d0

c     *** Silvestrelli G3 ***
      ixmat(1,6)=0.0d0
      ixmat(2,6)=0.0d0
      ixmat(3,6)=1.0d0

      do i=1,3
         bmat(i,1)=lattice_unitg(1,i)
         bmat(i,2)=lattice_unitg(2,i)
         bmat(i,3)=lattice_unitg(3,i)
      end do

      do i=1,6
         xx=0.0d0
         yy=0.0d0
         zz=0.0d0
         do j=1,3
           xx=xx+bmat(j,1)*ixmat(j,i)
           yy=yy+bmat(j,2)*ixmat(j,i)
           zz=zz+bmat(j,3)*ixmat(j,i)
         end do
         bv(1,i)=xx
         bv(2,i)=yy
         bv(3,i)=zz
      end do
      do i=1,6
         wrk(1,i)=bv(1,i)*bv(1,i)
         wrk(2,i)=bv(1,i)*bv(2,i)
         wrk(3,i)=bv(1,i)*bv(3,i)
         wrk(4,i)=bv(2,i)*bv(2,i)
         wrk(5,i)=bv(2,i)*bv(3,i)
         wrk(6,i)=bv(3,i)*bv(3,i)
         wts(i)=0.0d0
      end do

*     *** scal=(2*pi/L)**2 ***
      scal = lattice_omega()**(1.0d0/3.0d0)
      scal = 8.0*datan(1.0d0)/scal
      scal = scal*scal
      wts(1)=scal
      wts(4)=scal
      wts(6)=scal
      call dgels('N',6,6,1,wrk,6,wts,6,tmp1,tmp_len,info)
      if (info.ne.0) then
        write(*,*)"Illegal argument in call to dgels"
        call flush(6)
      end if
      rank=0
      do i=1,6
         if (dabs(wts(i)).gt.1.e-6) then
           rank=rank+1
           wrk(1,rank)=bv(1,i)
           wrk(2,rank)=bv(2,i)
           wrk(3,rank)=bv(3,i)
           wrk(4,rank)=wts(i)
         end if
      end do
      do i=1,rank
         bv(1,i)=wrk(1,i)
         bv(2,i)=wrk(2,i)
         bv(3,i)=wrk(3,i)
         wts(i)=wrk(4,i)
      end do

      if (oprint) then
         write(*,1260)
         write(*,1261) rank
         do i=1,rank
            write(*,1262) i,bv(1,i),bv(2,i),bv(3,i),wts(i)
         end do
      end if

      return
 1260 FORMAT(/5x,' Silvestrelli reciprocal lattice vectors:')
 1261 FORMAT(5x,' rank = ',I1)
 1262 FORMAT(5x,' bs',I1,' =<',3f8.3,' > weight=',f8.3)
      end


*     ************************************************
*     *                                              *
*     *          pspw_hfx_localize_end               *
*     *                                              *
*     ************************************************
      subroutine pspw_hfx_localize_end()
      implicit none

#include "bafdecls.fh"
#include "pspw_hfx.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value
      logical  Dneall_m_free
      external Dneall_m_free

      if (localize_on) then
         call nwpw_timing_start(52)
         value = BA_free_heap(rwav(2))
         value = value.and.BA_free_heap(dmindx(2))
         value = value.and.BA_free_heap(psiloc_r(2))
         value = value.and.BA_free_heap(switchr(2))
         value = value.and.Dneall_m_free(rotmat21)
         value = value.and.Dneall_m_free(absmat)
         value = value.and.BA_free_heap(closeenough(2))
         value = value.and.BA_free_heap(centeroverlap(2))
         value = value.and.BA_free_heap(roverlap(2))
         if (.not.value)
     >    call errquit('pspw_hfx_localize_end:freeing heap',3,MA_ERR)
         call nwpw_timing_end(52)
      end if
      return
      end




*     ************************************************
*     *                                              *
*     *            pspw_hfx_localize2_n              *
*     *                                              *
*     ************************************************
      subroutine pspw_hfx_localize2_n(npsi,psi_r,psi1,psi2,psi3)
      implicit none
      integer npsi
      real*8 psi_r(*)
      complex*16 psi1(*),psi2(*),psi3(*)

#include "bafdecls.fh"
#include "pspw_hfx.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value
      integer i,ms,mshift,msize,neq(2),ne(2)
      integer X(2,6),A(2)

*     **** external functions ****
      logical  Dneall_w_push_get,Dneall_w_pop_stack
      external Dneall_w_push_get,Dneall_w_pop_stack

      if (hfx_on) then
      if (localize2_on) then

*        **** push stack memory ****
         value = .true.
         do i=1,6
            value = value.and.Dneall_w_push_get(1,X(1,i))
         end do
         value = value.and.Dneall_w_push_get(0,A)
         if (.not. value)
     >      call errquit('pspw_hfx_localize_n:push stack',1,MA_ERR)

         call Dneall_neq(neq)
         ne(1) = norbs(1)
         ne(2) = norbs(2)
         call Pack_npack(1,npack1)
         call Dneall_m_size(1,msize)

         do ms=1,ispin
            mshift = (ms-1)*msize
            do i=1,rank
               call silvestrelli_overlap(
     >                     bv(1,i),ms,ne,neq,
     >                     psi_r,
     >                     dbl_mb(psiloc_r(1)),
     >                     dcpl_mb(X(1,i)))
            end do
            call silvestrelli_iterate(maxiter,ms,ne(ms),rank,wts,
     >                       dcpl_mb(X(1,1)),dcpl_mb(X(1,2)),
     >                       dcpl_mb(X(1,3)),dcpl_mb(X(1,4)),
     >                       dcpl_mb(X(1,5)),dcpl_mb(X(1,6)),
     >                       dcpl_mb(A(1)))

            call Dneall_wtom_Real(ms,dcpl_mb(A(1)),
     >                            dbl_mb(rotmat21(1)+mshift))
         end do

*        **** rotate psi to produce Wannier orbitals ****
         call Dneall_mm_transpose(0,dbl_mb(rotmat21(1)),dcpl_mb(A(1)))

         call dcopy(n2ft3d*(neq(1)+neq(2)),
     >              psi_r,1,
     >              dbl_mb(psiloc_r(1)),1)
         call Dneall_gmg_Multiply(0,
     >                  dbl_mb(psiloc_r(1)),n2ft3d,
     >                  dcpl_mb(A(1)),1.0d0,
     >                  psi_r,
     >                  0.0d0)

         if (npsi.gt.0) then
            call dcopy(2*npack1*(neq(1)+neq(2)),
     >                 psi1,1,
     >                 dbl_mb(psiloc_r(1)),1)
            call Dneall_fmf_Multiply(0,
     >                     dbl_mb(psiloc_r(1)),npack1,
     >                     dcpl_mb(A(1)),1.0d0,
     >                     psi1,
     >                     0.0d0)
         end if
         if (npsi.gt.1) then
            call dcopy(2*npack1*(neq(1)+neq(2)),
     >                 psi2,1,
     >                 dbl_mb(psiloc_r(1)),1)
            call Dneall_fmf_Multiply(0,
     >                     dbl_mb(psiloc_r(1)),npack1,
     >                     dcpl_mb(A(1)),1.0d0,
     >                     psi2,
     >                     0.0d0)
         end if
         if (npsi.gt.2) then
            call dcopy(2*npack1*(neq(1)+neq(2)),
     >                 psi3,1,
     >                 dbl_mb(psiloc_r(1)),1)
            call Dneall_fmf_Multiply(0,
     >                     dbl_mb(psiloc_r(1)),npack1,
     >                     dcpl_mb(A(1)),1.0d0,
     >                     psi3,
     >                     0.0d0)
         end if

*        **** pop stack memory ****
         value = Dneall_w_pop_stack(A)
         do i=6,1,-1
            value = value.and.Dneall_w_pop_stack(X(1,i))
         end do
         if (.not. value)
     >   call errquit('pspw_hfx_localize_n:pop stack',2,MA_ERR)


      end if
      end if

      return
      end

*     ************************************************
*     *                                              *
*     *          pspw_hfx_find_centers               *
*     *                                              *
*     ************************************************
      subroutine pspw_hfx_find_centers(psi_r,rcent)
      implicit none
      integer ispin,ne(2),neq(2)
      real*8  psi_r(*)
      real*8  rcent(3,*)

      return
      end


*     ************************************************
*     *                                              *
*     *        pspw_hfx_localize_start               *
*     *                                              *
*     ************************************************
      subroutine pspw_hfx_localize_start(psi_r)
      implicit none
      real*8  psi_r(*)

#include "bafdecls.fh"
#include "pspw_hfx.fh"
#include "errquit.fh"
#include "stdio.fh"

*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)

      logical value,oprint
      integer neq(2),ii,nx,ny,nz,q,ms,n,n1(2),n2(2),ne(2)
      integer i,j,k,center(3)
      real*8  scal1,a(3,3),rcenter(3)
      real*8  cdx1,cdx2,cdx3,cdy1,cdy2,cdy3,t1,t2,t3,x,y,z,twopi
      real*8  xx,yy,zz,dx,dy,dz,rr
      integer tmp1(2),xsgrid(2),xs(6)


*     **** external functions ****
      integer  ion_nearest_index
      external ion_nearest_index
      real*8   lattice_unita
      external lattice_unita
      character*2 ion_aname_nocarat
      external    ion_aname_nocarat


      if (localize_on) then
         call nwpw_timing_start(52)
         call Parallel_taskid(taskid)
         oprint = (taskid.eq.MASTER)
         oprint = .false.


         call Dneall_neq(neq)
         ne(1) = norbs(1)
         ne(2) = norbs(2)
         call D3dB_nx(1,nx)
         call D3dB_ny(1,ny)
         call D3dB_nz(1,nz)
         scal1 = 1.0d0/dble(nx*ny*nz)
         twopi = 8.0d0*datan(1.0d0)

*        ***** localize matrix ****
         if (localize2_on) then
*           **** use psi_r - assumed that it is localized ***
            call dcopy(n2ft3d*(neq(1)+neq(2)),
     >                 psi_r,1,
     >                 dbl_mb(psiloc_r(1)),1)
         else
               
*           **** generate localize psi_r ****
            do i=1,neq(1)+neq(2)
               call D3dB_r_Zero_Ends(1,psi_r(1+(i-1)*n2ft3d))
            end do
            call psi_dmatrix_maxallcolumns(ispin,ne,neq,n2ft3d,
     >                                     psi_r,dbl_mb(psiloc_r(1)),
     >                                     int_mb(dmindx(1)))
            call psi_dmatrix_rotate(ispin,ne,neq,n2ft3d,
     >                              int_mb(dmindx(1)),
     >                              psi_r,dbl_mb(psiloc_r(1)))

            do i=1,neq(1)+neq(2)
               call D3dB_r_Zero_Ends(1,dbl_mb(psiloc_r(1)+(i-1)*n2ft3d))
            end do

*           **** generate < psiloc_r|psi_r> matrix  ****
            call Dneall_ggm_Multiply(0,dbl_mb(psiloc_r(1)),psi_r,
     >                            scal1,
     >                            n2ft3d,dbl_mb(rotmat21(1)),0.0d0)

         end if



*        **************************************************
*        **** Find centers of localized orbitals, rwav ****
*        **************************************************
         value = BA_push_get(mt_dbl,n2ft3d,"tmp1",tmp1(2),tmp1(1))
         value = value.and.
     >           BA_push_get(mt_dbl,6*n2ft3d,"xsgrid",
     >                       xsgrid(2),xsgrid(1))
         if (.not.value)
     >   call errquit('pspw_hfx_localize_start:push stack',1,MA_ERR)
         a(1,1) = lattice_unita(1,1)
         a(2,1) = lattice_unita(2,1)
         a(3,1) = lattice_unita(3,1)
         a(1,2) = lattice_unita(1,2)
         a(2,2) = lattice_unita(2,2)
         a(3,2) = lattice_unita(3,2)
         a(1,3) = lattice_unita(1,3)
         a(2,3) = lattice_unita(2,3)
         a(3,3) = lattice_unita(3,3)

         call lattice_p_grid(dbl_mb(xsgrid(1)))

         xs(1) = xsgrid(1)
         xs(2) = xsgrid(1)+n2ft3d
         xs(3) = xsgrid(1)+2*n2ft3d
         xs(4) = xsgrid(1)+3*n2ft3d
         xs(5) = xsgrid(1)+4*n2ft3d
         xs(6) = xsgrid(1)+5*n2ft3d

         n1(1) = 1 
         n1(2) = neq(1)+1
         n2(1) = neq(1)
         n2(2) = neq(1)+neq(2)
         call dcopy(3*(ne(1)+ne(2)),0.0d0,0,dbl_mb(rwav(1)),1)
         do ms=1,ispin
            do q=n1(ms),n2(ms)
               call Dneall_qton(q,n)
               call D3dB_rr_Sqr(1,dbl_mb(psiloc_r(1)+(q-1)*n2ft3d),
     >                         dbl_mb(tmp1(1)))

               call D3dB_rr_dot(1,dbl_mb(xs(1)),dbl_mb(tmp1(1)),cdx1)
               call D3dB_rr_dot(1,dbl_mb(xs(2)),dbl_mb(tmp1(1)),cdy1)

               call D3dB_rr_dot(1,dbl_mb(xs(3)),dbl_mb(tmp1(1)),cdx2)
               call D3dB_rr_dot(1,dbl_mb(xs(4)),dbl_mb(tmp1(1)),cdy2)

               call D3dB_rr_dot(1,dbl_mb(xs(5)),dbl_mb(tmp1(1)),cdx3)
               call D3dB_rr_dot(1,dbl_mb(xs(6)),dbl_mb(tmp1(1)),cdy3)
               t1 = datan2(cdy1*scal1,cdx1*scal1)/twopi
               t2 = datan2(cdy2*scal1,cdx2*scal1)/twopi
               t3 = datan2(cdy3*scal1,cdx3*scal1)/twopi
               x  = a(1,1)*t1 + a(1,2)*t2 + a(1,3)*t3
               y  = a(2,1)*t1 + a(2,2)*t2 + a(2,3)*t3
               z  = a(3,1)*t1 + a(3,2)*t2 + a(3,3)*t3
               !dbl_mb(rwav(1)+3*((ms-1)*ne(1)+n-1))   = x
               !dbl_mb(rwav(1)+3*((ms-1)*ne(1)+n-1)+1) = y
               !dbl_mb(rwav(1)+3*((ms-1)*ne(1)+n-1)+2) = z
               dbl_mb(rwav(1)+3*(n-1))   = x
               dbl_mb(rwav(1)+3*(n-1)+1) = y
               dbl_mb(rwav(1)+3*(n-1)+2) = z

            end do
         end do
         call D1dB_Vector_SumAll(3*(ne(1)+ne(2)),dbl_mb(rwav(1)))

         value =           BA_pop_stack(xsgrid(2))
         value = value.and.BA_pop_stack(tmp1(2))
         if (.not.value)
     >    call errquit('pspw_hfx_localize_start:pop stack',2,MA_ERR)

*        ***********************************************************************
*        **** compute the overlap centers between localized orbital centers ****
*        ***********************************************************************
         do ms=1,ispin
            do n=1,ne(ms)
            do q=1,n
               call pspw_hfx_localize_calc_ovlpcntr(ms,n,q,a,
     >                  dbl_mb(rwav(1)+3*(ms-1)*ne(1)),rcenter)
               nx = 3*((q-1)+(n-1)*ne(ms) + (ms-1)*ne(1)*ne(1))
               ny = 3*((n-1)+(q-1)*ne(ms) + (ms-1)*ne(1)*ne(1))

               call lattice_center0_xyz_to_ijk(
     >                      rcenter(1),rcenter(2),rcenter(3),
     >                      center(1),center(2),center(3))

               dbl_mb(roverlap(1)+nx)   = rcenter(1)
               dbl_mb(roverlap(1)+nx+1) = rcenter(2)
               dbl_mb(roverlap(1)+nx+2) = rcenter(3)
               dbl_mb(roverlap(1)+ny)   = rcenter(1)
               dbl_mb(roverlap(1)+ny+1) = rcenter(2)
               dbl_mb(roverlap(1)+ny+2) = rcenter(3)

               int_mb(centeroverlap(1)+nx)   = center(1)
               int_mb(centeroverlap(1)+nx+1) = center(2)
               int_mb(centeroverlap(1)+nx+2) = center(3)
               int_mb(centeroverlap(1)+ny)   = center(1)
               int_mb(centeroverlap(1)+ny+1) = center(2)
               int_mb(centeroverlap(1)+ny+2) = center(3)
              
            end do
            end do
         end do


*        *************************************************
*        **** prune interactions using overlap matrix ****
*        *************************************************
         if (prune_algorithm.eq.0) then

*           **** compute overlap matrix of localized orbitals ****
            call Dneall_ggm_AMultiply(0,dbl_mb(psiloc_r(1)),
     >                                  dbl_mb(psiloc_r(1)),
     >                                  scal1,
     >                                  n2ft3d,dbl_mb(absmat(1)),0.0d0)
            do ms=1,ispin
               do n=1,ne(ms)
               do q=1,n
                  nx = (q-1)+(n-1)*ne(ms) + (ms-1)*ne(1)*ne(1)
                  ny = (n-1)+(q-1)*ne(ms) + (ms-1)*ne(1)*ne(1)
                  cdx1 = dbl_mb(absmat(1)+nx)

                  ! rin_thresh  = 5.0d-3
                  ! rout_thresh = 1.0d-2
                  if (cdx1.gt.rin_thresh) then
                     value = .true.
                     if (cdx1.lt.rout_thresh) then
                        cdy1  = (cdx1**2-(rin_thresh)**2)**2
     >                   *(3*(rout_thresh)**2-2*cdx1**2-(rin_thresh)**2)
     >                   /((rout_thresh)**2-(rin_thresh)**2)**3
                     else
                        cdy1  = 1.0d0
                     end if
                  else
                     value = .false.
                     cdy1  = 0.0d0
                  end if
                  log_mb(closeenough(1)+nx) = value
                  log_mb(closeenough(1)+ny) = value
                  dbl_mb(switchr(1)+nx) = cdy1
                  dbl_mb(switchr(1)+ny) = cdy1
               end do
               end do
            end do


*        *****************************************************************************
*        **** prune interactions using distance between localized orbital centers ****
*        *****************************************************************************
         else

            do ms=1,ispin
               do n=1,ne(ms)
               do q=1,n
                  x = dbl_mb(rwav(1)+3*((ms-1)*ne(1)+n-1))
                  y = dbl_mb(rwav(1)+3*((ms-1)*ne(1)+n-1)+1)
                  z = dbl_mb(rwav(1)+3*((ms-1)*ne(1)+n-1)+2)

                  t1= dbl_mb(rwav(1)+3*((ms-1)*ne(1)+q-1))
                  t2= dbl_mb(rwav(1)+3*((ms-1)*ne(1)+q-1)+1)
                  t3= dbl_mb(rwav(1)+3*((ms-1)*ne(1)+q-1)+2)
                  xx = x-t1
                  yy = y-t2
                  zz = z-t3
                  cdx1 = 9.99e9
                  do k = -1,1
                  do j = -1,1
                  do i = -1,1
                     dx = xx+a(1,1)*i+a(1,2)*j+a(1,3)*k
                     dy = yy+a(2,1)*i+a(2,2)*j+a(2,3)*k
                     dz = zz+a(3,1)*i+a(3,2)*j+a(3,3)*k
                     rr = dsqrt(dx*dx + dy*dy + dz*dz)
                     if (rr.lt.cdx1) cdx1 = rr
                  end do
                  end do
                  end do
                  if (cdx1.lt.rout_thresh) then
                     value = .true.
                     if (cdx1.le.rin_thresh) then
                        cdy1 = 1.0d0
                     else
                        cdy1 = (rout_thresh**2 - cdx1**2)**2
     >                       *(rout_thresh**2+2*cdx1**2-3*rin_thresh**2)
     >                       /(rout_thresh**2 - rin_thresh**2)**3
                     end if
                  else
                     value = .false.
                     cdy1  = 0.0d0
                  end if
                  nx = (q-1)+(n-1)*ne(ms) + (ms-1)*ne(1)*ne(1)
                  ny = (n-1)+(q-1)*ne(ms) + (ms-1)*ne(1)*ne(1)
                  log_mb(closeenough(1)+nx) = value
                  log_mb(closeenough(1)+ny) = value
                  dbl_mb(switchr(1)+nx) = cdy1
                  dbl_mb(switchr(1)+ny) = cdy1
               end do
               end do
            end do


*           **** print out centers ****
            if (oprint) then
               write(luout,1382)
               do ms=1,ispin
                  do n=1,ne(ms)
                     x = dbl_mb(rwav(1)+3*((ms-1)*ne(1)+n-1))
                     y = dbl_mb(rwav(1)+3*((ms-1)*ne(1)+n-1)+1)
                     z = dbl_mb(rwav(1)+3*((ms-1)*ne(1)+n-1)+2)
                     ii = ion_nearest_index(x,y,z)
                     if (ms.eq.1) then
                        write(luout,1390) n,x,y,z,
     >                                    ion_aname_nocarat(ii),ii
                     end if
                     if (ms.eq.2) then
                        write(luout,1391) n+ne(1),x,y,z,
     >                                    ion_aname_nocarat(ii),ii
                     end if
                  end do
               end do
            end if
         end if

      call nwpw_timing_end(52)
      end if
      
      return
 1382 FORMAT(/' position of localized orbital centers:')
 1390 FORMAT(5X,'psi',I5,' (',3F11.5,
     >       ' ) - spin up orbital, nearest ion=',A3,I5)
 1391 FORMAT(5X,'psi',I5,' (',3F11.5,
     >       ' ) - spin down orbital, nearest ion=',A3,I5)

      end


*     ************************************************
*     *                                              *
*     *        pspw_hfx_localize_stop                *
*     *                                              *
*     ************************************************
      subroutine pspw_hfx_localize_stop(Hpsi_r)
      implicit none
      real*8  Hpsi_r(*)

#include "bafdecls.fh"
#include "pspw_hfx.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value
      integer neq(2),ne(2)
      integer tmat1(2),tmat2(2),nx,ny,nz,ms,msize
      real*8  scal1

      logical  Dneall_m_push_get,Dneall_m_pop_stack
      external Dneall_m_push_get,Dneall_m_pop_stack


      if (localize_on.or.localize2_on) then
         call nwpw_timing_start(52)

*        **** allocate temp matrices from stack ***
         value =           Dneall_m_push_get(0,tmat2)
         value = value.and.Dneall_m_push_get(0,tmat1)
         if (.not.value)
     >     call errquit('pspw_hfx_localize_stop:push stack',1,MA_ERR)

         call Dneall_m_size(0,msize)
         call Dneall_neq(neq)
         call D3dB_nx(1,nx)
         call D3dB_ny(1,ny)
         call D3dB_nz(1,nz)
         scal1 = 1.0d0/dble(nx*ny*nz)
         ne(1) = norbs(1)
         ne(2) = norbs(2)

*        **** Make Hpsi Hermitian ****
         call Dneall_ggm_Multiply(0,dbl_mb(psiloc_r(1)),Hpsi_r,
     >                            scal1,
     >                            n2ft3d,dbl_mb(tmat1(1)),0.0d0)
         call Dneall_mm_transpose(0,dbl_mb(tmat1(1)),dbl_mb(tmat2(1)))
         call daxpy(msize,-1.0d0,dbl_mb(tmat1(1)),1,dbl_mb(tmat2(1)),1)
         call dscal(msize,0.5d0,dbl_mb(tmat2(1)),1)

         call Dneall_gmg_Multiply(0,dbl_mb(psiloc_r(1)),n2ft3d,
     >                              dbl_mb(tmat2(1)), 1.0d0,
     >                              Hpsi_r,           1.0d0)

*        **** rotate Hermitian Hpsi if not localize2 ****
         if (.not.localize2_on) then
            call Dneall_gmg_Multiply(0,Hpsi_r,n2ft3d,
     >                              dbl_mb(rotmat21(1)),1.0d0,
     >                              dbl_mb(psiloc_r(1)),0.0d0)
            call dcopy(n2ft3d*(neq(1)+neq(2)),dbl_mb(psiloc_r(1)),1,
     >                 Hpsi_r,1)
         end if


*        **** pop memory from stack ****
         value =           Dneall_m_pop_stack(tmat1)
         value = value.and.Dneall_m_pop_stack(tmat2)
         if (.not.value)
     >     call errquit('pspw_hfx_localize_stop:pop stack',2,MA_ERR)

         call nwpw_timing_end(52)
      end if

      return
      end

*     ************************************************
*     *                                              *
*     *        pspw_hfx_localize_set_findcolumns     *
*     *                                              *
*     ************************************************
      subroutine pspw_hfx_localize_set_findcolumns(fc)
      implicit none
      logical fc

#include "bafdecls.fh"
#include "pspw_hfx.fh"
#include "errquit.fh"

      findcolumns = fc
      return
      end
 
*     ************************************************
*     *                                              *
*     *       pspw_hfx_localize_calc_ovlpcntr        *
*     *                                              *
*     ************************************************
      subroutine pspw_hfx_localize_calc_ovlpcntr(ms,i,j,a,rwav,rcenter)
      implicit none
      integer ms
      integer i,j
      real*8  a(3,3),rwav(3,*)
      real*8  rcenter(3)

*     **** local variables ****
      integer i1,i2,i3,ii1,ii2,ii3
      real*8 x1,y1,z1,x2,y2,z2,dx,dy,dz,xc,yc,zc,r2min,r2

      rcenter(1) = 0.0d0
      rcenter(2) = 0.0d0
      rcenter(3) = 0.0d0

      x1 = rwav(1,i)
      y1 = rwav(2,i)
      z1 = rwav(3,i)

      x2= rwav(1,j)
      y2= rwav(2,j)
      z2= rwav(3,j)
         
      ii1 = 0
      ii2 = 0
      ii3 = 0
      r2min = 99.9d99
      do i3=-1,1 
      do i2=-1,1 
      do i1=-1,1 
         dx = (x1-x2+a(1,1)*i1+a(1,2)*i2+a(1,3)*i3)
         dy = (y1-y2+a(2,1)*i1+a(2,2)*i2+a(2,3)*i3)
         dz = (z1-z2+a(3,1)*i1+a(3,2)*i2+a(3,3)*i3)
         r2 = dx*dx + dy*dy + dz*dz
         if (r2.lt.r2min) then
            r2min = r2
            ii1 = i1
            ii2 = i2
            ii3 = i3
         end if
      end do
      end do
      end do
      xc = 0.5d0*(x1 + x2 + a(1,1)*ii1 + a(1,2)*ii2 + a(1,3)*ii3)
      yc = 0.5d0*(y1 + y2 + a(2,1)*ii1 + a(2,2)*ii2 + a(2,3)*ii3)
      zc = 0.5d0*(z1 + z2 + a(3,1)*ii1 + a(3,2)*ii2 + a(3,3)*ii3)
      call lattice_min_difference(xc,yc,zc)
      rcenter(1) = xc
      rcenter(2) = yc
      rcenter(3) = zc

      return
      end


*     ************************************************
*     *                                              *
*     *        pspw_hfx_localize_closeenough         *
*     *                                              *
*     ************************************************

      logical function pspw_hfx_localize_closeenough(i,j)
      implicit none
      integer i,j

#include "bafdecls.fh"
#include "pspw_hfx.fh"
#include "errquit.fh"
      
*     **** local variables ****
      integer indx,ne(2)

      if (localize_on) then
         ne(1) = norbs(1)
         ne(2) = norbs(2)
         indx = 0
         if ((i.gt.ne(1)).and.(j.gt.ne(1))) then
            indx = (i-ne(1)-1) + (j-ne(1)-1)*ne(2) + ne(1)*ne(1)
         else if ((i.le.ne(1)).and.(j.le.ne(1))) then
            indx = (i-1) + (j-1)*ne(1) 
         else
           call errquit('pspw_hfx_localize_closeenough: failed',0,UERR)
         end if
         pspw_hfx_localize_closeenough = log_mb(closeenough(1)+indx)
      else
         pspw_hfx_localize_closeenough = .true.
      end if
      return
      end

*     ************************************************
*     *                                              *
*     *        pspw_hfx_localize_switchr             *
*     *                                              *
*     ************************************************

      real*8 function pspw_hfx_localize_switchr(i,j)
      implicit none
      integer i,j

#include "bafdecls.fh"
#include "pspw_hfx.fh"
#include "errquit.fh"

*     **** local variables ****
      integer indx,ne(2)

      if (localize_on) then
         ne(1) = norbs(1)
         ne(2) = norbs(2)
         indx = 0
         if ((i.gt.ne(1)).and.(j.gt.ne(1))) then
            indx = (i-ne(1)-1) + (j-ne(1)-1)*ne(2) + ne(1)*ne(1)
         else if ((i.le.ne(1)).and.(j.le.ne(1))) then
            indx = (i-1) + (j-1)*ne(1)
         else
           call errquit('pspw_hfx_localize_switchr: failed',0,UERR)
         end if
         pspw_hfx_localize_switchr = dbl_mb(switchr(1)+indx)
      else
         pspw_hfx_localize_switchr = 1.0d0
      end if
      return
      end

*     ************************************************
*     *                                              *
*     *        pspw_hfx_localize_rcenter_ovlp        *
*     *                                              *
*     ************************************************

      subroutine pspw_hfx_localize_rcenter_ovlp(i,j,rcenter)
      implicit none
      integer i,j
      real*8 rcenter(3)

#include "bafdecls.fh"
#include "pspw_hfx.fh"
#include "errquit.fh"

*     **** local variables ****
      integer indx,ne(2)

      rcenter(1) = 0.0d0
      rcenter(2) = 0.0d0
      rcenter(3) = 0.0d0
      if (localize_on) then
         ne(1) = norbs(1)
         ne(2) = norbs(2)
         indx = 0
         if ((i.gt.ne(1)).and.(j.gt.ne(1))) then
            indx = (i-ne(1)-1) + (j-ne(1)-1)*ne(2) + ne(1)*ne(1)
         else if ((i.le.ne(1)).and.(j.le.ne(1))) then
            indx = (i-1) + (j-1)*ne(1)
         else
           call errquit('pspw_hfx_localize_rcenter_ovlp:failed',0,UERR)
         end if
         rcenter(1) = dbl_mb(roverlap(1)+3*indx)
         rcenter(2) = dbl_mb(roverlap(1)+3*indx+1)
         rcenter(3) = dbl_mb(roverlap(1)+3*indx+2)
      end if

      return
      end


*     ************************************************
*     *                                              *
*     *        pspw_hfx_localize_center_ovlp         *
*     *                                              *
*     ************************************************
      subroutine pspw_hfx_localize_center_ovlp(i,j,center)
      implicit none
      integer i,j
      integer center(3)

#include "bafdecls.fh"
#include "pspw_hfx.fh"
#include "errquit.fh"

*     **** local variables ****
      integer indx,ne(2)

      center(1) = 0
      center(2) = 0
      center(3) = 0
      if (localize_on) then
         ne(1) = norbs(1)
         ne(2) = norbs(2)
         indx = 0
         if ((i.gt.ne(1)).and.(j.gt.ne(1))) then
            indx = (i-ne(1)-1) + (j-ne(1)-1)*ne(2) + ne(1)*ne(1)
         else if ((i.le.ne(1)).and.(j.le.ne(1))) then
            indx = (i-1) + (j-1)*ne(1)
         else
           call errquit('pspw_hfx_localize_center_ovlp:failed',0,UERR)
         end if
         center(1) = int_mb(centeroverlap(1)+3*indx)
         center(2) = int_mb(centeroverlap(1)+3*indx+1)
         center(3) = int_mb(centeroverlap(1)+3*indx+2)
      end if

      return
      end




c*     ************************************************
c*     *                                              *
c*     *        pspw_potential_hfx_localize           *
c*     *                                              *
c*     ************************************************
c
c      subroutine pspw_potential_hfx_localize(psi_r,Hpsi_r)
c      implicit none
c      real*8 psi_r(*)
c      real*8 Hpsi_r(*)
c
c#include "bafdecls.fh"
c#include "pspw_hfx.fh"
c#include "errquit.fh"
c
c      !*** variables to use from pspw_hfx.fh ***
c      !* ispin - 1-restricted, 2-unrestricted
c      !* norbs(2) - number of up orbitals and down orbitals
c      !* n2ft3d   - number of real space grid points ~ (nfft1+2)*nfft2*nfft3/nprocs
c      !*******************************************
c
c      !*** the exchange correction has to be added (or appended) to Hpsi_r not just set
c
c
c      !**** need to set these variables ****
c      ehfx = 0.0
c      phfx = 0.0  !* phfx = 2*ehfx
c
c
c      write(*,*) "Hello from pspw_potential_localized_hfx",Hpsi_r(11)
c
c
c      return
c      end
