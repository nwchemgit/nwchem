
*     ************************************
*     *                    	         *
*     *       pspw_atom_analysis         *
*     *                                  *
*     ************************************
      subroutine pspw_atom_analysis(rtdb,ispin,n2ft3d,dn)
      implicit none
      integer rtdb
      integer ispin
      integer n2ft3d
      real*8 dn(n2ft3d,ispin)

#include "bafdecls.fh"
#include "btdb.fh"
#include "errquit.fh"

*     ***** local variables ****
      integer ma_type,nactive_atoms,h_actlist,l_actlist
      integer ii
      real*8 rcpp(2)


*     **** pspw:atom_analysis ****
      if (rtdb_ma_get(rtdb, 'pspw:atom_analysis', ma_type,
     >        nactive_atoms, h_actlist)) then


         if (.not.BA_get_index(h_actlist,l_actlist))
     >   call errquit(
     >     'pspw_Lubin_Water_analysis: ma_get_index failed',911, MA_ERR)

         if (.not.btdb_get(rtdb,'pspw:atom_analysis_rcut_pp', 
     >                    mt_dbl,2,rcpp)) then
           rcpp(1) = 8.0d0
           rcpp(2) = 8.0d0
         end if

         do ii=1,nactive_atoms
           call pspw_sub_atom_analysis(ispin,n2ft3d,dn,
     >                                 int_mb(l_actlist+ii-1),rcpp)
         end do

*       **** free heap ****
        if (.not. BA_free_heap(h_actlist))
     >   call errquit('h_actlist:error freeing heap memory',0, MA_ERR)

      end if

      return
      end

*     ************************************
*     *                    	         *
*     *    pspw_sub_atom_analysis        *
*     *                         	 *
*     ************************************
      subroutine pspw_sub_atom_analysis(ispin,n2ft3d,dn,indx,rcpp)
      implicit none
      integer ispin
      integer n2ft3d
      real*8 dn(n2ft3d,ispin)
      integer indx
      real*8 rcpp(2)

#include "bafdecls.fh"
#include "errquit.fh"


*     **** Parameters ****
      integer taskid,MASTER
      parameter (MASTER=0)

*     **** local variables ****
      logical value
      integer k,ii
      integer nx,ny,nz
      real*8 ncharge,x
      real*8 dv,rcut,pp,pp2

*     **** MA stack memory ****
      integer rgrid(2),rgx(2),fcut(2)

 
*     **** external functions ****
      real*8   lattice_omega,ion_rion
      external lattice_omega,ion_rion

      call Parallel_taskid(taskid)

*     **** allocate memory ****
      value = BA_push_get(mt_dbl,3*n2ft3d,'rgrid', rgrid(2), rgrid(1))
      value = value.and.
     >        BA_push_get(mt_dbl,n2ft3d,'rgx', rgx(2), rgx(1))
      value = value.and.
     >        BA_push_get(mt_dbl, n2ft3d,'fcut',fcut(2),fcut(1))     
      if (.not. value) 
     >   call errquit('atom_analysis: out of stack memory',0, MA_ERR)

      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      dv=lattice_omega()/dble(nx*ny*nz)
      call lattice_r_grid(dbl_mb(rgrid(1)))

*     **** define fcut =  *****
      rcut = rcpp(1)
      pp   = rcpp(2)
      pp2  = pp+2.0d0
      do k=1,n2ft3d
         x = (dbl_mb(rgrid(1)+3*(k-1))   - ion_rion(1,indx))**2
     >     + (dbl_mb(rgrid(1)+3*(k-1)+1) - ion_rion(2,indx))**2
     >     + (dbl_mb(rgrid(1)+3*(k-1)+2) - ion_rion(3,indx))**2
         x=dsqrt(x)
         dbl_mb(fcut(1)+k-1)=(1.0d0-(1.0d0-dexp(-(x/rcut)**pp2))**pp)
      end do


*     **** calculate ncharge ****
      call D3dB_r_Zero(1,dbl_mb(rgx(1)))
      call D3dB_rr_Sum(1,dn(1,1),dn(1,ispin),dbl_mb(rgx(1)))
      call D3dB_rr_dot(1,dbl_mb(rgx(1)),dbl_mb(fcut(1)),ncharge)
      ncharge = ncharge*dv

      if (taskid.eq.MASTER) then
         WRITE(*,100)
         WRITE(*,110)
         WRITE(*,111) rcut
         WRITE(*,112) pp
         WRITE(*,115)
         WRITE(*,120) ion_rion(1,indx),
     >                ion_rion(2,indx),
     >                ion_rion(3,indx)
         WRITE(*,131) ncharge
      end if

*     **** pop stack memory ****
      value = value.and.BA_pop_stack(fcut(2))
      value = value.and.BA_pop_stack(rgx(2))
      value = value.and.BA_pop_stack(rgrid(2))
      if (.not. value) 
     >   call errquit('atom analysis: cannot pop stack memory',0,
     >       MA_ERR)

      return

  100 format(//'Atom Density Analysis')
  110 format( '----------------------')
  111 format('Sphere radius     :',F10.4)
  112 format('Boundary sharpness:',F10.4)
  115 format(/'Atom position:') 
  120 format('RO           ','  (',F10.4,',',F10.4,',',F10.4,' )')
  131 format('Qe = ',F10.4,' au')  

      end



c $Id$
