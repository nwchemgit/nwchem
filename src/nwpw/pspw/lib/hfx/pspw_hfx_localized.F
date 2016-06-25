
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
         call Dneall_neq(neq)
         ne(1) = norbs(1)
         ne(2) = norbs(2)

         value = BA_alloc_get(mt_dbl,3*(ne(1)+ne(2)),
     >                        'rwav',rwav(2),rwav(1))
         value = value.and.
     >           BA_alloc_get(mt_dbl,n2ft3d*(neq(1)+neq(2)),
     >                        'psiloc_r',psiloc_r(2),psiloc_r(1))
         value = BA_alloc_get(mt_log,(ne(1)*ne(1)+ne(2)*ne(2)),
     >                'closeenough',closeenough(2),closeenough(1))
         value = value.and.Dneall_m_allocate(0,rotmat)
         if (.not.value)
     >      call errquit('pspw_hfx_localize_init: out of heap',3,MA_ERR)
      end if

      return
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
         value = BA_free_heap(rwav(2))
         value = value.and.BA_free_heap(psiloc_r(2))
         value = value.and.BA_free_heap(closeenough(2))
         value = value.and.Dneall_m_free(rotmat)
         if (.not.value)
     >    call errquit('pspw_hfx_localize_end:freeing heap',3,MA_ERR)
      end if
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
      real*8  scal1,a(3,3)
      real*8  cdx1,cdx2,cdx3,cdy1,cdy2,cdy3,t1,t2,t3,x,y,z,twopi
      integer tmp1(2),xsgrid(2),xs(6)

*     **** external functions ****
      integer  ion_nearest_index
      external ion_nearest_index
      real*8   lattice_unita
      external lattice_unita
      character*2 ion_aname_nocarat
      external    ion_aname_nocarat

      if (localize_on) then
         call Parallel_taskid(taskid)
         oprint = (taskid.eq.MASTER)
         oprint = .false.

         value = BA_push_get(mt_dbl,n2ft3d,"tmp1",tmp1(2),tmp1(1))
         value = value.and.
     >        BA_push_get(mt_dbl,6*n2ft3d,"xsgrid",xsgrid(2),xsgrid(1))
         if (.not.value)
     >    call errquit('pspw_hfx_localize_start:push stack',1,MA_ERR)

         call Dneall_neq(neq)
         ne(1) = norbs(1)
         ne(2) = norbs(2)
         call D3dB_nx(1,nx)
         call D3dB_ny(1,ny)
         call D3dB_nz(1,nz)
         scal1 = 1.0d0/dble(nx*ny*nz)
         twopi = 8.0d0*datan(1.0d0)

*        **** generate localize psi_r ****
         call psi_dmatrix_localize(ispin,neq,n2ft3d,
     >                             psi_r,dbl_mb(psiloc_r(1)))

*        **** generate < psiloc_r|psi_r< matrix  ****
         call Dneall_ggm_Multiply(0,dbl_mb(psiloc_r(1)),psi_r,scal1,
     >                           n2ft3d,dbl_mb(rotmat(1)),0.0d0)

*        **** Find centers of localized orbitals, rwav ****
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
         call dcopy(3*(ne(1)+ne(2)),0.0d0,0,dbl_mb(Rwav(1)),1)
         do ms=1,ispin
            do q=n1(ms),n2(ms)
               call Dneall_qton(q,n)
               call D3dB_rr_Sqr(1,dbl_mb(psiloc_r(1)+(q-1)*n2ft3d),
     >                            dbl_mb(tmp1(1)))

               call D3dB_rr_dot(1,dbl_mb(xs(1)),dbl_mb(tmp1(1)),cdx1)
               call D3dB_rr_dot(1,dbl_mb(xs(2)),dbl_mb(tmp1(1)),cdy1)

               call D3dB_rr_dot(1,dbl_mb(xs(3)),dbl_mb(tmp1(1)),cdx2)
               call D3dB_rr_dot(1,dbl_mb(xs(4)),dbl_mb(tmp1(1)),cdy2)

               call D3dB_rr_dot(1,dbl_mb(xs(5)),dbl_mb(tmp1(1)),cdx3)
               call D3dB_rr_dot(1,dbl_mb(xs(6)),dbl_mb(tmp1(1)),cdy3)
               t1 = datan2(cdy1,cdx1)/twopi
               t2 = datan2(cdy2,cdx2)/twopi
               t3 = datan2(cdy3,cdx3)/twopi
               x  = a(1,1)*t1 + a(1,2)*t2 + a(1,3)*t3
               y  = a(2,1)*t1 + a(2,2)*t2 + a(2,3)*t3
               z  = a(3,1)*t1 + a(3,2)*t2 + a(3,3)*t3
               dbl_mb(rwav(1)+3*((ms-1)*ne(1)+n-1))   = x
               dbl_mb(rwav(1)+3*((ms-1)*ne(1)+n-1)+1) = y
               dbl_mb(rwav(1)+3*((ms-1)*ne(1)+n-1)+2) = z

            end do
         end do
         call D1dB_Vector_SumAll(3*(ne(1)+ne(2)),dbl_mb(Rwav(1)))

         value =           BA_pop_stack(xsgrid(2))
         value = value.and.BA_pop_stack(tmp1(2))
         if (.not.value)
     >    call errquit('pspw_hfx_localize_start:pop stack',2,MA_ERR)

         do ms=1,ispin
            do n=1,ne(ms)
            do q=1,n
               x = dbl_mb(rwav(1)+3*((ms-1)*ne(1)+n-1))
               y = dbl_mb(rwav(1)+3*((ms-1)*ne(1)+n-1)+1)
               z = dbl_mb(rwav(1)+3*((ms-1)*ne(1)+n-1)+2)

               t1= dbl_mb(rwav(1)+3*((ms-1)*ne(1)+q-1))
               t2= dbl_mb(rwav(1)+3*((ms-1)*ne(1)+q-1)+1)
               t3= dbl_mb(rwav(1)+3*((ms-1)*ne(1)+q-1)+2)
               cdx1 = dsqrt((x-t1)**2 + (y-t2)**2 + (z-t3)**2)
               if (cdx1.le.rwav_thresh) then
                  value = .true.
               else
                  value = .false.
               end if

               nx = (q-1)+(n-1)*ne(ms) + (ms-1)*ne(1)*ne(1)
               ny = (n-1)+(q-1)*ne(ms) + (ms-1)*ne(1)*ne(1)
               log_mb(closeenough(1)+nx) = value
               log_mb(closeenough(1)+ny) = value
            end do
            end do
         end do

*        **** print out centers ****
         if (oprint) then
            write(luout,1382)
            do ms=1,ispin
               do n=1,ne(ms)
                  x = dbl_mb(rwav(1)+3*((ms-1)*ne(1)+n-1))
                  y = dbl_mb(rwav(1)+3*((ms-1)*ne(1)+n-1)+1)
                  z = dbl_mb(rwav(1)+3*((ms-1)*ne(1)+n-1)+2)
                  ii = ion_nearest_index(x,y,z)
                  if (ms.eq.1) then
                     write(luout,1390) n,x,y,z,ion_aname_nocarat(ii),ii
                  end if
                  if (ms.eq.2) then
                     write(luout,1391) n+ne(1),x,y,z,
     >                              ion_aname_nocarat(ii),ii
                  end if
               end do
            end do
         end if

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
      integer neq(2)

      if (localize_on) then
         call Dneall_neq(neq)
         call Dneall_gmg_Multiply(0,Hpsi_r,n2ft3d,
     >                           dbl_mb(rotmat(1)),  1.0d0,
     >                           dbl_mb(psiloc_r(1)),0.0d0)
         call dcopy(n2ft3d*(neq(1)+neq(2)),dbl_mb(psiloc_r(1)),1,
     >              Hpsi_r,1)
      end if

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
