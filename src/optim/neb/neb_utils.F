*
* $Id$
*

*     ***********************************************
*     *                                             *
*     *             neb_path_energy                 *
*     *                                             *
*     ***********************************************
      subroutine neb_path_energy(bead_list,dpath,epath)
      implicit none
      character*(*) bead_list
      real*8 dpath
      real*8 epath

#include "mafdecls.fh"

*     **** local variables ****
      logical value
      integer i,index,index_p,nbeads,nion,ng
      integer e(2),c(2),r(2)
      real*8  dist
      real*8  amat(3,3)
     
*     **** external functions ****
      integer  size_bead_list,nion_bead_list
      real*8   ddot
      external size_bead_list,nion_bead_list
      external ddot

      nbeads = size_bead_list(bead_list)
      nion   = nion_bead_list(bead_list,1)
      ng     = 3*nion*nbeads

      call amatrix_get_bead_list(bead_list,1,amat)

      value  = MA_alloc_get(mt_dbl,nbeads,'e',e(2),e(1))
      value  = value.and.
     >         MA_alloc_get(mt_dbl,3*nion,'r',r(2),r(1))
      value  = value.and.
     >         MA_alloc_get(mt_dbl,ng,'c',c(2),c(1))
      if (.not.value) 
     > call errquit('neb_path_energy failed - increase stack',0,0)
    
      call neb_energies_get(bead_list,dbl_mb(e(1)))
      call neb_coords_get(bead_list,dbl_mb(c(1)))
      
*     **** calculate the path length and path energy ****
      dpath = 0.0d0
      epath = 0.0d0
      do i=1,(nbeads-1)
         index   = (i-1)*3*nion + 1
         index_p = (i  )*3*nion + 1
         call neb_cdiff(nion,dbl_mb(c(1)+index_p-1), 
     >                       dbl_mb(c(1)+index-1),
     >                       amat,dbl_mb(r(1)))
c         call dcopy(3*nion,dbl_mb(c(1)+index_p-1),1,
c     >                     dbl_mb(r(1)),1)
c         call daxpy(3*nion,(-1.0d0),
c     >              dbl_mb(c(1)+index-1),1,
c     >              dbl_mb(r(1)),   1)
         dist = dsqrt(ddot(3*nion,dbl_mb(r(1)),1,
     >                            dbl_mb(r(1)),1))
         dpath = dpath + dist
         epath = epath + 0.5d0*dbl_mb(e(1)+i-1)*dist
         epath = epath + 0.5d0*dbl_mb(e(1)+i  )*dist
      end do
      value = value.and.MA_free_heap(c(2))
      value = value.and.MA_free_heap(r(2))
      value = value.and.MA_free_heap(e(2))      
      if (.not.value) call errquit('neb_path_energy failed',1,0)

      return
      end

*     ***********************************************
*     *                                             *
*     *             neb_gradient_get                *
*     *                                             *
*     ***********************************************
      subroutine neb_gradient_get(bead_list,kbeads,c,e,t,g,
     >                            freeze1,freezen)
      implicit none
      character*(*) bead_list
      real*8 kbeads
      real*8 c(*)
      real*8 e(*)
      real*8 t(*)
      real*8 g(*)
      logical freeze1,freezen

*     **** local variables ****
      integer i,index,nbeads,nion,istart,iend
      real*8  k,norm
      real*8  amat(3,3)
     
*     **** external functions ****
      integer  size_bead_list,nion_bead_list
      real*8   ddot
      external size_bead_list,nion_bead_list
      external ddot


      k = kbeads

      nbeads = size_bead_list(bead_list)
      nion   = nion_bead_list(bead_list,1)

      call amatrix_get_bead_list(bead_list,1,amat)

*     *** get the tangent ****
      call neb_tangent(nbeads,nion,c,e,t,amat)

*     **** normal gradient |g> = |g> - |t><t|g> ****
      call dcopy(nbeads*3*nion,0.0d0,0,g,1)
      !do i=2,(nbeads-1)
      istart = 2
      iend = nbeads-1
      if (.not.freeze1) istart = 1 
      if (.not.freezen) iend    = nbeads
      do i=istart,iend
         index = (i-1)*3*nion + 1
         call gradient_get_bead_list(bead_list,i,g(index))
         norm = ddot(3*nion,t(index),1,g(index),1)
         call daxpy(3*nion,(-norm),
     >              t(index),1,
     >              g(index),1)
      end do

*     *** add to normal gradient the tangent spring gradient *** 
      call neb_add_spring_gradient(nbeads,nion,c,t,k,g,amat)
      
      return
      end

      subroutine neb_gradient_get0(bead_list,g)
      implicit none
      character*(*) bead_list
      real*8 g(*)

*     **** local variables ****
      integer i,index,nbeads,nion

*     **** external functions ****
      integer  size_bead_list,nion_bead_list
      external size_bead_list,nion_bead_list


      nbeads = size_bead_list(bead_list)
      nion   = nion_bead_list(bead_list,1)

      call dcopy(nbeads*3*nion,0.0d0,0,g,1)
      do i=2,(nbeads-1)
         index = (i-1)*3*nion + 1
         call gradient_get_bead_list(bead_list,i,g(index))
      end do
      return
      end


      subroutine neb_project_gradient(nion,nbeads,t,g)
      implicit none
      integer nion,nbeads
      real*8 t(*),g(*)

*     **** local variables ****
      integer i,index
      real*8 norm

*     **** external functions ****
      real*8   ddot
      external ddot

      do i=2,(nbeads-1)
         index = (i-1)*3*nion+1
         norm = ddot(3*nion,t(index),1,g(index),1)
         call daxpy(3*nion,(-norm),
     >              t(index),1,
     >              g(index),1)
      end do

      return
      end
     

*     ***********************************************
*     *                                             *
*     *                  neb_cdiff                  *
*     *                                             *
*     ***********************************************
  
      subroutine neb_cdiff(nion,c2,c1,amat,c2mc1)
      implicit none
      integer nion
      real*8 c1(3,*),c2(3,*)
      real*8 amat(3,3)
      real*8 c2mc1(3,*)

*     **** local variables ****
      logical isperiodic
      integer ii,i1,i2,i3
      real*8 dx,dy,dz,dx0,dy0,dz0,dxf,dyf,dzf,rr,rrmin
      
      rr = 0.0d0
      do i2=1,3
      do i1=1,3
         rr = rr + dabs(amat(i1,i2))
      end do
      end do
      isperiodic = (rr.gt.3.001d0)

      do ii=1,nion
         dx0 = c2(1,ii) - c1(1,ii)
         dy0 = c2(2,ii) - c1(2,ii)
         dz0 = c2(3,ii) - c1(3,ii)
         dxf = dx0
         dyf = dy0
         dzf = dz0
         if (isperiodic) then
            rrmin = dx0*dx0 + dy0*dy0 + dz0*dz0
            do i3=-2,2
            do i2=-2,2
            do i1=-2,2
               dx = dx0 + amat(1,1)*i1 + amat(1,2)*i2 + amat(1,3)*i3
               dy = dy0 + amat(2,1)*i1 + amat(2,2)*i2 + amat(2,3)*i3
               dz = dz0 + amat(3,1)*i1 + amat(3,2)*i2 + amat(3,3)*i3
               rr = dx*dx + dy*dy + dz*dz
               if (rr.lt.rrmin) then
                  rrmin = rr
                  dxf = dx
                  dyf = dy
                  dzf = dz
               end if
            end do
            end do
            end do
         end if

         c2mc1(1,ii) = dxf
         c2mc1(2,ii) = dyf
         c2mc1(3,ii) = dzf
      end do

      return
      end

*     ***********************************************
*     *                                             *
*     *             neb_add_spring_gradient         *
*     *                                             *
*     ***********************************************
      subroutine neb_add_spring_gradient(nbeads,nion,c,t,k,gs,amat)
      implicit none
      integer nbeads,nion
      real*8 c(*)
      real*8 t(*)
      real*8 k
      real*8 gs(*)
      real*8 amat(3,3)

*     **** local variables ****
      integer i,index,index_m,index_p
      integer rp,rm
      real*8  normm,normp

*     **** external functions ***
      real*8   ddot
      external ddot

      rm = 1
      rp = (nbeads-1)*3*nion + 1
      
      do i=2,(nbeads-1)
         index   = (i-1)*3*nion + 1
         index_m = (i-2)*3*nion + 1
         index_p = (i  )*3*nion + 1
         call neb_cdiff(nion,c(index),c(index_m),amat,t(rm))
         call neb_cdiff(nion,c(index_p),c(index),amat,t(rp))
c         call dcopy(3*nion,c(index),1,t(rm),1)
c         call daxpy(3*nion,(-1.0d0),
c     >               c(index_m),1,
c     >               t(rm),     1)
c         call dcopy(3*nion,c(index_p),1,t(rp),1)
c         call daxpy(3*nion,(-1.0d0),
c     >               c(index),1,
c     >               t(rp),   1)
         normm = ddot(3*nion,t(rm),1,t(rm),1)
         normp = ddot(3*nion,t(rp),1,t(rp),1)
         normp = -k*(dsqrt(normp) - dsqrt(normm))
         call daxpy(3*nion,normp,
     >              t(index), 1,
     >              gs(index),1)
      end do
      call dcopy(3*nion,0.0d0,0,t(rm),1)
      call dcopy(3*nion,0.0d0,0,t(rp),1)

      return
      end

*     ***********************************************
*     *                                             *
*     *             neb_tangent                     *
*     *                                             *
*     ***********************************************
      subroutine neb_tangent(nbeads,nion,c,e,t,amat)
      implicit none
      integer nbeads,nion
      real*8 c(*)
      real*8 e(*)
      real*8 t(*)
      real*8 amat(3,3)

*     **** local variables ****
      integer i,index,index_m,index_p
      integer rp,rm
      real*8  norm,dVmax,dVmin

*     **** external functions ***
      real*8   ddot
      external ddot

      rm = 1
      rp = (nbeads-1)*3*nion + 1
      
      do i=2,(nbeads-1)
         index   = (i-1)*3*nion + 1
         index_m = (i-2)*3*nion + 1
         index_p = (i  )*3*nion + 1
         call neb_cdiff(nion,c(index),c(index_m),amat,t(rm))
         call neb_cdiff(nion,c(index_p),c(index),amat,t(rp))

c         call dcopy(3*nion,c(index),1,t(rm),1)
c         call daxpy(3*nion,(-1.0d0),
c     >               c(index_m),1,
c     >               t(rm),     1)
c         call dcopy(3*nion,c(index_p),1,t(rp),1)
c         call daxpy(3*nion,(-1.0d0),
c     >               c(index),1,
c     >               t(rp),   1)

         if   ( ((e(i+1)-e(i)).gt.1.0e-2).and.
     >          ((e(i)-e(i-1)).gt.1.0d-2) ) then
            call dcopy(3*nion,t(rp),1,t(index),1)
         else if ( ((e(i-1)-e(i)).gt.1.0d-2) .and.
     >             ((e(i)-e(i+1)).gt.1.0d-2) ) then
            call dcopy(3*nion,t(rm),1,t(index),1)
         else
            
             if ( (dabs(e(i+1)-e(i))-dabs(e(i-1)-e(i)))
     >           .gt.1.0d-3) then
               dVmax = dabs(e(i+1)-e(i))
               dVmin = dabs(e(i-1)-e(i))
             else
               dVmax = dabs(e(i-1)-e(i))
               dVmin = dabs(e(i+1)-e(i))
             end if

             if ((e(i+1)-e(i-1)).gt.1.0d-2) then
               call dscal(3*nion,dVmax,t(rp),1)
               call dscal(3*nion,dVmin,t(rm),1)
             else
               call dscal(3*nion,dVmin,t(rp),1)
               call dscal(3*nion,dVmax,t(rm),1)
             end if
             call dcopy(3*nion,t(rp),1,t(index),1)
             call daxpy(3*nion,(1.0d0),
     >                  t(rm),1,
     >                  t(index),1)
         end if

*        *** normalize tangent ***
         norm = ddot(3*nion,t(index),1,t(index),1)
         norm = 1.0d0/dsqrt(norm)
         call dscal(3*nion,norm,t(index),1)

      end do
      call dcopy(3*nion,0.0d0,0,t(rm),1)
      call dcopy(3*nion,0.0d0,0,t(rp),1)

      return
      end

*     ***********************************************
*     *                                             *
*     *             neb_gradient_get1               *
*     *                                             *
*     ***********************************************

*     *** RRR SUBRT TO GET AND ADJUST GRADIENT ***
*     *** RRR CHANGE TO INCLUDE CLIBING IMAGE ***
*     *** RRR CHANGE COMBINE PROJECTION and sforce
*     *** RRR addition into one subroutine

      subroutine neb_gradient_get1(bead_list,kbeads,c,e,t,g)


      implicit none
      character*(*) bead_list
      real*8 kbeads
      real*8 c(*)
      real*8 e(*)
      real*8 t(*)
      real*8 g(*)

#include "rtdb.fh"
#include "mafdecls.fh"
#include "geom.fh"
#include "global.fh"

*     **** local variables ****
      integer i,index,index_m, index_p
      integer nbeads, nion, rm, rp
      real*8  k,norm, normm, normp
      real*8  sprf, sprnorm
      real*8  amat(3,3)
     
*     **** external functions ****
      integer  size_bead_list,nion_bead_list
      real*8   ddot
      external size_bead_list,nion_bead_list
      external ddot

      k = kbeads
      nbeads = size_bead_list(bead_list)
      nion   = nion_bead_list(bead_list,1)
      call amatrix_get_bead_list(bead_list,1,amat)

      rm = 1
      rp = (nbeads-1)*3*nion + 1


*     *** get the tangent ****
      call neb_tangent(nbeads,nion,c,e,t,amat)


*    *** RRR project out tangent PES force

*     **** normal gradient |g> = |g> - |t><t|g> ****
      call dcopy(nbeads*3*nion,0.0d0,0,g,1)

*    *** RRR loop over ions to be moved
      do i=2,(nbeads-1)

*    *** RRR set index numbers
         index = (i-1)*3*nion + 1
         index_m = (i-2)*3*nion + 1
         index_p = (i  )*3*nion + 1

*    *** RRR get original gradients
         call gradient_get_bead_list(bead_list,i,g(index))

*    *** RRR dot product of 
*    *** RRR original force and tangent vector
*    *** RRR to get magnitude of gradient
*    *** RRR  parallel to tangent direction
         norm = ddot(3*nion,t(index),1,g(index),1)

*    *** RRR check relative energy of current bead
*    *** RRR if both nneighbors lower, then
*    *** RRR invert parallel force component
*    *** RRR to climb uphill along MEP
         if ((e(i).gt.e(i-1)).and.(e(i).gt.e(i+1))) then
                  call daxpy(3*nion,2*(-norm),
     >              t(index),1,
     >              g(index),1)
         else
*    *** MV for other images use regular projection ***
             call daxpy(3*nion,-norm,
     >       t(index),1,
     >       g(index),1)
          

*    *** RRR calculate spring forces
            call neb_cdiff(nion,c(index),c(index_m),amat,t(rm))
            call neb_cdiff(nion,c(index_p),c(index),amat,t(rp))
c                 call dcopy(3*nion,c(index),1,t(rm),1)
c                 call daxpy(3*nion,(-1.0d0),
c     >               c(index_m),1,
c     >               t(rm),     1)
c                 call dcopy(3*nion,c(index_p),1,t(rp),1)
c                 call daxpy(3*nion,(-1.0d0),
c     >               c(index),1,
c     >               t(rp),   1)

                 normm = ddot(3*nion,t(rm),1,t(rm),1)
                 normp = ddot(3*nion,t(rp),1,t(rp),1)
                 sprf = -k*(dsqrt(normp) - dsqrt(normm))
                 sprnorm = sprf - norm
                 call daxpy(3*nion,sprnorm,
     >              t(index), 1,
     >              g(index),1)
         endif
        
      end do

*   *** RRR END OF FORCE MODIFICATION
      call dcopy(3*nion,0.0d0,0,t(rm),1)
      call dcopy(3*nion,0.0d0,0,t(rp),1)

      return
      end

*     ***********************************************
*     *                                             *
*     *             neb_tangent1                    *
*     *                                             *
*     ***********************************************

*    *** RRR GET TANGENT TO BEAD PATH
*    *** RRR  DEFINE TANGENT AS THE VECTOR
*    *** RRR   TO THE NEIGHBOR BEAD THAT IS
*    *** RRR   HIGHER IN ENERGY
*    *** RRR  IE COMPARE e(i-1) to e(e+1)
*    *** RRR  if e(i-1) > e(e+1) then
*    *** RRR  tan goes from e(i+1) to e(e-1)
    

      subroutine neb_tangent1(nbeads,nion,c,e,t,amat)
      implicit none
      integer nbeads,nion
      real*8 c(*)
      real*8 e(*)
      real*8 t(*)
      real*8 amat(3,3)

*     **** local variables ****
      integer i,index,index_m,index_p
      integer rp,rm
      real*8  norm,dVmax,dVmin

*     **** external functions ***
      real*8   ddot
      external ddot

      rm = 1
      rp = (nbeads-1)*3*nion + 1
      
      do i=2,(nbeads-1)
         index   = (i-1)*3*nion + 1
         index_m = (i-2)*3*nion + 1
         index_p = (i  )*3*nion + 1
         call neb_cdiff(nion,c(index),c(index_m),amat,t(rm))
         call neb_cdiff(nion,c(index_p),c(index),amat,t(rp))

c         call dcopy(3*nion,c(index),1,t(rm),1)
c         call daxpy(3*nion,(-1.0d0),
c     >               c(index_m),1,
c     >               t(rm),     1)
c         call dcopy(3*nion,c(index_p),1,t(rp),1)
c         call daxpy(3*nion,(-1.0d0),
c     >               c(index),1,
c     >               t(rp),   1)

         if      ((e(i+1).gt.e(i)).and.(e(i).gt.e(i-1))) then
            call dcopy(3*nion,t(rp),1,t(index),1)
         else if ((e(i-1).gt.e(i)).and.(e(i).gt.e(i+1))) then
            call dcopy(3*nion,t(rm),1,t(index),1)
         else
            
             if (dabs(e(i+1)-e(i)).gt.dabs(e(i-1)-e(i))) then
               dVmax = dabs(e(i+1)-e(i))
               dVmin = dabs(e(i-1)-e(i))
             else
               dVmax = dabs(e(i-1)-e(i))
               dVmin = dabs(e(i+1)-e(i))
             end if

             if (e(i+1).gt.e(i-1)) then
               call dscal(3*nion,dVmax,t(rp),1)
               call dscal(3*nion,dVmin,t(rm),1)
             else
               call dscal(3*nion,dVmin,t(rp),1)
               call dscal(3*nion,dVmax,t(rm),1)
             end if
             call dcopy(3*nion,t(rp),1,t(index),1)
             call daxpy(3*nion,(1.0d0),
     >                  t(rm),1,
     >                  t(index),1)
         end if

*        *** normalize tangent ***
         norm = ddot(3*nion,t(index),1,t(index),1)
         norm = 1.0d0/dsqrt(norm)
         call dscal(3*nion,norm,t(index),1)

      end do
      call dcopy(3*nion,0.0d0,0,t(rm),1)
      call dcopy(3*nion,0.0d0,0,t(rp),1)

      return
      end

*     ***********************************************
*     *                                             *
*     *             neb_energies_get                *
*     *                                             *
*     ***********************************************

      subroutine neb_energies_get(bead_list,e)
      implicit none
      character*(*) bead_list
      real*8 e(*)

*     **** local variables ****
      integer i,nbeads

*     **** external functions ****
      integer  size_bead_list
      real*8   energy_bead_list
      external size_bead_list
      external energy_bead_list

      nbeads = size_bead_list(bead_list)
    
      do i=1,nbeads
         e(i) = energy_bead_list(bead_list,i)
      end do

      return
      end

*     ***********************************************
*     *                                             *
*     *             neb_coords_get                  *
*     *                                             *
*     ***********************************************
      subroutine neb_coords_get(bead_list,c)
      implicit none
      character*(*) bead_list
      real*8 c(*)

*     **** local variables ****
      integer i,index,nbeads,nion
      real*8 amat(3,3)

*     **** external functions ****
      integer  size_bead_list,nion_bead_list
      external size_bead_list,nion_bead_list

      nbeads = size_bead_list(bead_list)
      nion   = nion_bead_list(bead_list,1)

      do i=1,(nbeads)
        index = (i-1)*3*nion+1
        call coords_get_bead_list(bead_list,i,c(index))
      end do

      return
      end

*     ***********************************************
*     *                                             *
*     *             neb_coords_set                  *
*     *                                             *
*     ***********************************************
      subroutine neb_coords_set(bead_list,c)
      implicit none
      character*(*) bead_list
      real*8 c(*)

*     **** local variables ****
      integer i,index,nbeads,nion

*     **** external functions ****
      integer  size_bead_list,nion_bead_list
      external size_bead_list,nion_bead_list

      nbeads = size_bead_list(bead_list)
      nion   = nion_bead_list(bead_list,1)

      do i=1,nbeads
        index = (i-1)*3*nion+1
        call coords_set_bead_list(bead_list,i,c(index))
      end do
      return
      end

*     ***********************************************
*     *                                             *
*     *             neb_masses_get                  *
*     *                                             *
*     ***********************************************
      subroutine  neb_masses_get(rtdb,bead_list,m)
      implicit none
      integer rtdb
      character*(*) bead_list
      real*8 m(*)

#include "geom.fh"

*     **** local variables ****
      logical value
      integer geom,nion,bead_listlen
      character*255 geom_name
      integer geomlen
*     **** external functions ****
      integer     inp_strlen
      character*7 bead_index_name
      external    inp_strlen
      external    bead_index_name


      bead_listlen = inp_strlen(bead_list)
      !geom_name   = 'bead'//bead_index_name(1)//':geom'
      geom_name   = bead_list(1:bead_listlen)
     >             //bead_index_name(1)//':geom'
      geomlen     = inp_strlen(geom_name)

      value = geom_create(geom,'neb_tmp')
      value = value.and.geom_rtdb_load(rtdb,geom,
     >                                 geom_name(1:geomlen))
      value = value.and.geom_ncent(geom,nion)
      value = value.and.geom_masses_get(geom,nion,m)
      value = value.and.geom_destroy(geom)
      if (.not.value) call errquit('neb_masses failed',0,0)

      call dscal(nion,1822.89d0,m,1)
      return
      end

*     ***********************************************
*     *                                             *
*     *             neb_initialize                  *
*     *                                             *
*     ***********************************************

      subroutine neb_initialize(rtdb, bead_list)
      implicit none
      integer rtdb
      character*(*) bead_list 

#include "nwc_const.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "inp.fh"
c     
c     This routine initializes the common /coptopt/ and
c     also creates and returns the geometry handle
c     
      integer nbeads
      character*80 neb_movecs
      logical custom_path

      if (.not.rtdb_cget(rtdb,'neb:movecs',1,neb_movecs)) then
         call util_file_prefix('movecs',neb_movecs)
      end if

      if (.not.rtdb_get(rtdb,'neb:custom_path',mt_log,1,custom_path)) 
     >    custom_path = .false.
 
      if (.not.rtdb_get(rtdb,'neb:nbeads',mt_int,1,nbeads)) then
         nbeads = 5
      end if

*     **** create bead_list *** 
      call init_bead_list(rtdb,bead_list,neb_movecs)
      if(custom_path) then
        call neb_initial_path_custom(rtdb,bead_list,nbeads)
      else
        call neb_initial_path(rtdb,bead_list,nbeads)
      end if

      return
      end

*     ***********************************************
*     *                                             *
*     *             neb_initial_path                *
*     *                                             *
*     ***********************************************

      subroutine neb_initial_path(rtdb,bead_list,nbeads)
      implicit none
      integer rtdb
      character*(*) bead_list
      integer nbeads

#include "geom.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "stdio.fh"

*     **** local variables ****
      logical value,impose,oprint,hasmiddle
      integer i,geom,geomlen,movecslen,nion,nfit,i2,i3,pathguess
      integer bead_listlen
      integer r1(2),r2(2),r3(2),rmid(2),ifit(2),wfit(2),rcoords(2)
      real*8  t,rms1,rms2
      character*255 geom_name,movecs_name

*     **** external functions ****
      integer     inp_strlen
      character*7 bead_index_name
      external    inp_strlen
      external    bead_index_name

      oprint = (ga_nodeid().eq.0)
      bead_listlen = inp_strlen(bead_list)

      if (.not. rtdb_get(rtdb,'neb:pathguess',mt_int,1,pathguess))
     $      pathguess = 2

      value = geom_create(geom,'neb_tmp')

      !*** check for neb_start, otherwise just read geometry ****
      if (.not.geom_rtdb_load(rtdb,geom,'neb_start')) then
         value = value.and.geom_rtdb_load(rtdb,geom,'geometry')
      end if

      value = value.and.geom_ncent(geom,nion)
      value = value.and.geom_destroy(geom)
      if (.not.value) call errquit('neb_initial_path failed',0,0)

      value = value.and.MA_push_get(mt_dbl, (3*nion), 'r1',r1(2),r1(1))
      value = value.and.MA_push_get(mt_dbl, (3*nion), 'r2',r2(2),r2(1))
      value = value.and.MA_push_get(mt_dbl, (3*nion), 'r3',r3(2),r3(1))
      value = value.and.MA_push_get(mt_dbl, (3*nion),
     >                              'rmid',rmid(2),rmid(1))
      if (.not.value) call errquit('neb_initial_path failed',1,0)

      value = value.and.geom_create(geom,'neb_tmp')

      !*** try neb_end, then endgeom ****
      if (.not.geom_rtdb_load(rtdb,geom,'neb_end')) then
         value = value.and.geom_rtdb_load(rtdb,geom,'endgeom')
      end if
      value = value.and.geom_cart_coords_get(geom,dbl_mb(r2(1)))

      !*** try neb_start, then geometry****
      if (.not.geom_rtdb_load(rtdb,geom,'neb_start')) then
         value = value.and.geom_rtdb_load(rtdb,geom,'geometry')
      end if
      value = value.and.geom_cart_coords_get(geom,dbl_mb(r1(1)))
      if (.not.value) call errquit('neb_initial_path failed',2,0)

      if (.not.rtdb_get(rtdb,'neb:hasmiddle',mt_log,1,hasmiddle))
     >   hasmiddle = .false.

      if (hasmiddle) then
         if (.not.geom_rtdb_load(rtdb,geom,'neb_middle')) then
            if (.not.geom_rtdb_load(rtdb,geom,'midgeom')) then
               hasmiddle = .false.
            end if
         else
            if (.not.geom_cart_coords_get(geom,dbl_mb(rmid(1))))
     >      call errquit('neb_initial_path failed',2,0)
         end if
      end if

      if (oprint) then
         write(luout,*) 
     >   " - Generating initial path by linear interpolation"
         write(luout,'(A,I4)') "    + number images = ",nbeads
         if (hasmiddle) then
            write(luout,'(3A)')    
     >      "    + neb_start (geometry) geometry -->",
     >      " neb_middle (midgeom) geometry -->",
     >      " neb_end (endgeom) geometry"
         else
            write(luout,'(A)')    
     >      "    + neb_start (geometry) geometry -->",
     >      " neb_end (endgeom) geometry"
         end if
      end if



      if (.not.rtdb_get(rtdb,'neb:impose',mt_log,1,impose))
     >   impose = .false.

      if (impose) then
         value =  MA_push_get(mt_int,(2*nion),'ifit',ifit(2),ifit(1))
         value = value.and.
     >            MA_push_get(mt_dbl,(nion),'wfit',wfit(2),wfit(1))
         if(.not.value) call errquit('neb_initial_path failed',3,MA_ERR)

         if (hasmiddle) then
            call neb_impose(nion,dbl_mb(r1(1)),dbl_mb(rmid(1)),
     >                   nfit,int_mb(ifit(1)),dbl_mb(wfit(1)),rms1,rms2)
            if (oprint) then
            write(luout,*) 
     >      " - Imposing neb_mid (midgeom) geometry onto",
     >      " neb_start (geometry) geometry"
            write(luout,'(A,F10.6)') "    + initial rmsq = ",rms1
            write(luout,'(A,F10.6)') "    + imposed rmsq = ",rms2
            end if

            call neb_impose(nion,dbl_mb(rmid(1)),dbl_mb(r2(1)),
     >                   nfit,int_mb(ifit(1)),dbl_mb(wfit(1)),rms1,rms2)
            if (oprint) then
            write(luout,*) 
     >      " - Imposing neb_end geometry onto neb_mid geometry"
            write(luout,'(A,F10.6)') "    + initial rmsq = ",rms1
            write(luout,'(A,F10.6)') "    + imposed rmsq = ",rms2
            end if
         else
            call neb_impose(nion,dbl_mb(r1(1)),dbl_mb(r2(1)),
     >                   nfit,int_mb(ifit(1)),dbl_mb(wfit(1)),rms1,rms2)
            if (oprint) then
            write(luout,*) 
     >      " - Imposing neb_end (endgeom) geometry",
     >      " onto neb_start (geometry) geometry"
            write(luout,'(A,F10.6)') "    + initial rmsq = ",rms1
            write(luout,'(A,F10.6)') "    + imposed rmsq = ",rms2
            end if
         end if

         value =           MA_pop_stack(wfit(2))
         value = value.and.MA_pop_stack(ifit(2))
         if(.not.value) call errquit('neb_initial_path failed',4,MA_ERR)
      end if



      if(.not.ma_push_get(mt_dbl,3*nion*nbeads,'rcoords',
     >                    rcoords(2), rcoords(1)))
     >     call errquit('neb_initialize: memory',3*nion*nbeads,MA_ERR)



      if (hasmiddle) then
         i2 = nbeads/2+1

        call dcopy(3*nion,dbl_mb(r1(1)),1,dbl_mb(rcoords(1)),1) 
        call dcopy(3*nion,dbl_mb(rmid(1)),1,
     >            dbl_mb(rcoords(1)+3*nion*(i2-1)),1) 
        call dcopy(3*nion,dbl_mb(r2(1)),1,
     >            dbl_mb(rcoords(1)+3*nion*(nbeads-1)),1) 

        call zts_guess(nion,nbeads,i2,dbl_mb(rcoords(1)),pathguess)

        do i=1,nbeads
          call dcopy(3*nion,dbl_mb(rcoords(1)+3*nion*(i-1)),1,
     >                      dbl_mb(r3(1)),1)

          !geom_name   = 'bead'//bead_index_name(i)//':geom'
          !movecs_name = 'bead'//bead_index_name(i)//'.movecs'
          geom_name   = bead_list(1:bead_listlen)
     >                  //bead_index_name(i)//':geom'
          movecs_name = bead_list(1:bead_listlen)
     >                  //bead_index_name(i)//'.movecs'
          geomlen     = inp_strlen(geom_name)
          movecslen   = inp_strlen(movecs_name)
          value = value.and.geom_cart_coords_set(geom,dbl_mb(r3(1)))
          value = value.and.geom_rtdb_store(rtdb,geom,
     >                                    geom_name(1:geomlen))
          call add_bead_list(bead_list,
     >                     movecs_name(1:movecslen),
     >                     geom_name(1:geomlen))
        end do


c        do i=1,i2
c          t = (i-1)/dble(i2-1)
c          call dcopy(3*nion,dbl_mb(r1(1)),1,dbl_mb(r3(1)),1)
c          call dscal(3*nion,(1.0d0-t),dbl_mb(r3(1)),1)
c          call daxpy(3*nion,t,dbl_mb(rmid(1)),1,dbl_mb(r3(1)),1)
c
c          geom_name   = 'bead'//bead_index_name(i)//':geom'
c          movecs_name = 'bead'//bead_index_name(i)//'.movecs'
c          geomlen     = inp_strlen(geom_name)
c          movecslen   = inp_strlen(movecs_name)
c          value = value.and.geom_cart_coords_set(geom,dbl_mb(r3(1)))
c          value = value.and.geom_rtdb_store(rtdb,geom,
c     >                                    geom_name(1:geomlen))
c          call add_bead_list(bead_list,
c     >                     movecs_name(1:movecslen),
c     >                     geom_name(1:geomlen))
c        end do
c
c        i3 = nbeads-i2
c        do i=1,i3
c          t = i/dble(i3)
c          call dcopy(3*nion,dbl_mb(rmid(1)),1,dbl_mb(r3(1)),1)
c          call dscal(3*nion,(1.0d0-t),dbl_mb(r3(1)),1)
c          call daxpy(3*nion,t,dbl_mb(r2(1)),1,dbl_mb(r3(1)),1)
c
c          geom_name   = 'bead'//bead_index_name(i+i2)//':geom'
c          movecs_name = 'bead'//bead_index_name(i+i2)//'.movecs'
c          geomlen     = inp_strlen(geom_name)
c          movecslen   = inp_strlen(movecs_name)
c          value = value.and.geom_cart_coords_set(geom,dbl_mb(r3(1)))
c          value = value.and.geom_rtdb_store(rtdb,geom,
c     >                                    geom_name(1:geomlen))
c          call add_bead_list(bead_list,
c     >                     movecs_name(1:movecslen),
c     >                     geom_name(1:geomlen))
c        end do

*     **** linear interpolation with Robinson Checking ****
      else
      
        call dcopy(3*nion,dbl_mb(r1(1)),1,dbl_mb(rcoords(1)),1) 
        call dcopy(3*nion,dbl_mb(r2(1)),1,
     >            dbl_mb(rcoords(1)+3*nion*(nbeads-1)),1) 

        call zts_guessall(nion,nbeads,dbl_mb(rcoords(1)),geom)

        do i=1,nbeads
c          t = (i-1)/dble(nbeads-1)
c
c          call dcopy(3*nion,dbl_mb(r1(1)),1,dbl_mb(r3(1)),1)
c          call dscal(3*nion,(1.0d0-t),dbl_mb(r3(1)),1)
c          call daxpy(3*nion,t,dbl_mb(r2(1)),1,dbl_mb(r3(1)),1)

          call dcopy(3*nion,dbl_mb(rcoords(1)+3*nion*(i-1)),1,
     >                      dbl_mb(r3(1)),1)
 
          !geom_name   = 'bead'//bead_index_name(i)//':geom'
          !movecs_name = 'bead'//bead_index_name(i)//'.movecs'
          geom_name   = bead_list(1:bead_listlen)
     >                  //bead_index_name(i)//':geom'
          movecs_name = bead_list(1:bead_listlen)
     >                  //bead_index_name(i)//'.movecs'
          geomlen     = inp_strlen(geom_name)
          movecslen   = inp_strlen(movecs_name)
          value = value.and.geom_cart_coords_set(geom,dbl_mb(r3(1)))
          value = value.and.geom_rtdb_store(rtdb,geom,
     >                                    geom_name(1:geomlen))

          call add_bead_list(bead_list,
     >                     movecs_name(1:movecslen),
     >                     geom_name(1:geomlen))

        end do
      end if
      value = value.and.MA_pop_stack(rcoords(2))
      value = value.and.geom_destroy(geom)
      value = value.and.MA_pop_stack(rmid(2))
      value = value.and.MA_pop_stack(r3(2))
      value = value.and.MA_pop_stack(r2(2))
      value = value.and.MA_pop_stack(r1(2))
      if (.not.value) call errquit('neb_initial_path failed',3,0)


      return
      end

*     ***********************************************
*     *                                             *
*     *             neb_initial_path_custom         *
*     *                                             *
*     ***********************************************
      subroutine neb_initial_path_custom(rtdb,bead_list,nbeads)
      implicit none
      integer rtdb
      character*(*) bead_list
      integer nbeads

#include "mafdecls.fh"
#include "geom.fh"

*     **** local variables ****
      logical value
      integer i,geomlen,movecslen,bead_listlen
      character*255 geom_name,movecs_name

*     **** external functions ****
      integer     inp_strlen
      character*7 bead_index_name
      external    inp_strlen
      external    bead_index_name

      bead_listlen = inp_strlen(bead_list)
      do i=1,nbeads
        !geom_name   = 'bead'//bead_index_name(i)//':geom'
        !movecs_name = 'bead'//bead_index_name(i)//'.movecs'
        geom_name   = bead_list(1:bead_listlen)
     >                //bead_index_name(i)//':geom'
        movecs_name = bead_list(1:bead_listlen)
     >                //bead_index_name(i)//'.movecs'
        geomlen     = inp_strlen(geom_name)
        movecslen   = inp_strlen(movecs_name)

        call add_bead_list(bead_list,
     >                     movecs_name(1:movecslen),
     >                     geom_name(1:geomlen))

      end do

      return
      end

*     ***********************************************
*     *                                             *
*     *             neb_verlet_update               *
*     *                                             *
*     ***********************************************

      subroutine neb_verlet_update(ng,c0,c1,v1,dti,g1)
        integer ng
        double precision c0(*)
        double precision c1(*)
        double precision v1(*)
        double precision dti(*)
        double precision g1(*)

        integer i

*        *** c1 <- 2*c1 - c0 + dti*g          ***
*        ***    <- c1 + ((c1-c0)/t)*t + dti*g ***

*        *** RRR loop over ions *** 
*        *** NEED: CONSTRAINED DYNAMICS (FROZEN ATOMS)
         do i=1,ng
           v1(i) = c1(i)-c0(i)
           if (v1(i)*g1(i).gt.0.0d0) v1(i) = 0.0d0
         end do
         call dcopy(ng,c1(1),1,c0(1),1)
         do i=1,ng
            c1(i) = 
     >              c1(i)
     >            + v1(i)  
     >            + dti(i)*g1(i)
         end do

      end

      subroutine neb_cg_direction(ng,g0,g1,s)
      integer ng
      double precision g0(*)
      double precision g1(*)
      double precision s(*)

      integer i
      double precision gamma1
      double precision gamma2
      double precision sn

c     *** choosing Polac-Ribiere coeff ***        
      gamma1=0.0d0
      gamma2=0.0d0
      do i=1,ng
         gamma1 = gamma1 + (g1(i)-g0(i))*g1(i)
         gamma2 = gamma2 + g0(i)*g0(i)
      end do

      do i=1,ng
         s(i) = -g1(i) + s(i)*gamma1/gamma2
      end do
        
      return
      end

      subroutine neb_move(ng,dt,c0,c1,s)
      integer ng
      double precision dt
      double precision c0(*)
      double precision c1(*)
      double precision s(*)

      integer i

      do i=1,ng
          c1(i) = c0(i)+dt*s(i)
      end do

      return
      end

      subroutine neb_calc_convergence(ng,g1,c0,c1,Gmax,Grms,Xmax,Xrms)
      implicit none
      integer ng
      double precision g1(*)
      double precision c0(*)
      double precision c1(*)
      double precision Gmax,Grms,Xmax,Xrms

#include "global.fh"

      !***** local variables ****
      logical oprint
      integer ii,imax
      double precision dx


      oprint = ga_nodeid() .eq. 0

      Gmax = 0.0d0
      imax = 0
      do ii=1,ng
        if (dabs(g1(ii)).gt.Gmax) then
            Gmax = dabs(g1(ii))
            imax = ii
        end if
      end do
      if (oprint) write(*,*) "neb: imax,Gmax=",imax,Gmax

      Grms = 0.0d0
      do ii=1,ng
        Grms = Grms + g1(ii)*g1(ii)
      end do
      Grms = dsqrt(Grms/dble(ng))

      Xmax = 0.0d0
      do ii=1,ng
        dx = dabs(c1(ii)-c0(ii))
        if (dabs(g1(ii)).gt.Xmax) Xmax = dx
      end do

      Xrms = 0.0d0
      do ii=1,ng
        dx = dabs(c1(ii)-c0(ii))
        Xrms = Xrms + dx*dx
      end do
      Xrms = dsqrt(Xrms/dble(ng))


      return
      end


*     *****************************************
*     *                                       *
*     *             neb_lmbfgs                *
*     *                                       *
*     *****************************************
      subroutine neb_lmbfgs(n,m,x,g,hg)
      implicit none
      integer n,m
      real*8 x(n,m)
      real*8 g(n,m)
      real*8 hg(n)

*     **** local variables ****
      integer k
      real*8 rho(25),alpha(25),beta(25)
      real*8 tmp

*     **** external functions ****
      real*8   ddot
      external ddot

*     **** compute rho(k) = 1/y(:,k)' * s(:,k) ****
      do k=1,m-1
         tmp =       ddot(n,x(1,k+1),1,g(1,k+1),1)
         tmp = tmp - ddot(n,x(1,k+1),1,g(1,k),  1)
         tmp = tmp - ddot(n,x(1,k),  1,g(1,k+1),1)
         tmp = tmp + ddot(n,x(1,k),  1,g(1,k),  1)
         if (dabs(tmp).gt.1.0d-9) then
             rho(k) = 1.0d0/tmp
         else
             rho(k) = 0.0d0
         end if
      end do

      call dcopy(n,g(1,m),1,hg,1)

      do k = (m-1),1,-1
        alpha(k) = rho(k)
     >            *(ddot(n,x(1,k+1),1,hg,1) - ddot(n,x(1,k),1,hg,1))

        call daxpy(n,(-alpha(k)),g(1,k+1),1,hg,1)
        call daxpy(n,( alpha(k)),g(1,k),  1,hg,1)
      end do

      do k = 1,(m-1)
        beta(k) = rho(k)
     >           *(ddot(n,g(1,k+1),1,hg,1) - ddot(n,g(1,k),1,hg,1))

        call daxpy(n,(alpha(k)-beta(k)),x(1,k+1),1,hg,1)
        call daxpy(n,(beta(k)-alpha(k)),x(1,k),  1,hg,1)
      end do

      return
      end


*     ***********************************************
*     *                                             *
*     *             neb_resize_path                 *
*     *                                             *
*     ***********************************************
      subroutine neb_resize_path(rtdb,bead_list,nbeads1,nbeads2)
      implicit none
      integer rtdb
      character*(*) bead_list
      integer nbeads1
      integer nbeads2

#include "mafdecls.fh"
#include "geom.fh"

*     **** local variables ****
      logical value
      integer i,geom,geomlen,movecslen,nion,bead_listlen
      integer j1,j2,shift
      integer r1(2),r2(2),r3(2),c(2)
      real*8  t,t1,t2,t3
      character*255 geom_name,movecs_name

*     **** external functions ****
      integer     inp_strlen
      character*7 bead_index_name
      external    inp_strlen
      external    bead_index_name

      bead_listlen = inp_strlen(bead_list)

      value = geom_create(geom,'neb_tmp')
      if (.not.geom_rtdb_load(rtdb,geom,'neb_start')) then
          value = value.and.geom_rtdb_load(rtdb,geom,'geometry')
      end if
      value = value.and.geom_ncent(geom,nion)
      value = value.and.geom_destroy(geom)
      if (.not.value) call errquit('neb_resize_path failed',0,0)


      value = value.and.MA_push_get(mt_dbl,(3*nion*nbeads1),
     >                              'c',c(2),c(1))
      value = value.and.MA_push_get(mt_dbl,(3*nion),'r1',r1(2),r1(1))
      value = value.and.MA_push_get(mt_dbl,(3*nion),'r2',r2(2),r2(1))
      value = value.and.MA_push_get(mt_dbl,(3*nion),'r3',r3(2),r3(1))
      if (.not.value) call errquit('neb_resize_path failed',1,0)

      value = value.and.geom_create(geom,'neb_tmp')
      if (.not.geom_rtdb_load(rtdb,geom,'neb_start')) then
         value = value.and.geom_rtdb_load(rtdb,geom,'geometry')
      end if
      value = value.and.geom_cart_coords_get(geom,dbl_mb(r3(1)))
      if (.not.value) call errquit('neb_resize_path failed',2,0)

      do i=1,nbeads1
         shift = (i-1)*3*nion
         call coords_get_bead_list(bead_list,i,dbl_mb(c(1)+shift))
      end do
      call reset_bead_list(bead_list)

      do i=1,nbeads2
        t = (i-1)/dble(nbeads2-1)

        j1 = t*(nbeads1-1) + 1
        j2 = j1+1
        t1 = (j1-1)/dble(nbeads1-1)
        t2 = (j2-1)/dble(nbeads1-1)
        t3 = (t-t1)/(t2-t1)

        if (j2.gt.nbeads1) then
           t3 = 0.0d0
           j2=nbeads1
        end if

        shift = (j1-1)*3*nion
        call dcopy(3*nion,dbl_mb(c(1)+shift),1,dbl_mb(r1(1)),1)

        shift = (j2-1)*3*nion
        call dcopy(3*nion,dbl_mb(c(1)+shift),1,dbl_mb(r2(1)),1)

        call dcopy(3*nion,dbl_mb(r1(1)),1,dbl_mb(r3(1)),1)
        call dscal(3*nion,(1.0d0-t3),dbl_mb(r3(1)),1)
        call daxpy(3*nion,t3,dbl_mb(r2(1)),1,dbl_mb(r3(1)),1)

        !geom_name   = 'bead'//bead_index_name(i)//':geom'
        !movecs_name = 'bead'//bead_index_name(i)//'.movecs'
        geom_name   = bead_list(1:bead_listlen)
     >                //bead_index_name(i)//':geom'
        movecs_name = bead_list(1:bead_listlen)
     >                //bead_index_name(i)//'.movecs'
        geomlen     = inp_strlen(geom_name)
        movecslen   = inp_strlen(movecs_name)
        value = value.and.geom_cart_coords_set(geom,dbl_mb(r3(1)))
        value = value.and.geom_rtdb_store(rtdb,geom,
     >                                    geom_name(1:geomlen))

        call add_bead_list(bead_list,
     >                     movecs_name(1:movecslen),
     >                     geom_name(1:geomlen))

      end do
      value = value.and.geom_destroy(geom)
      value = value.and.MA_pop_stack(r3(2))
      value = value.and.MA_pop_stack(r2(2))
      value = value.and.MA_pop_stack(r1(2))
      value = value.and.MA_pop_stack(c(2))
      if (.not.value) call errquit('neb_new_path failed',3,0)
      return
      end


*     ***********************************************
*     *                                             *
*     *             neb_impose                      *
*     *                                             *
*     ***********************************************
c     subroutine neb_impose  --  superimpose two coordinate sets
c
c     This routine performs the least squares best superposition
c     of two atomic coordinate sets via a quaternion method;
c     upon return, the first coordinate set is unchanged while
c     the second set is translated and rotated to give best fit;
c     the final root mean square fit is returned in "rmsvalue"
c
      subroutine neb_impose(nion,rion1,rion2,nfit,ifit,wfit,rms1,rms2)
      implicit none
      integer nion
      real*8 rion1(3,*),rion2(3,*)
      integer nfit,ifit(2,*)
      real*8  wfit(*)
      real*8 rms1,rms2

*     **** local variables ****
      integer i
      real*8 xmid,ymid,zmid

*     **** external functions ****
      real*8   neb_rmsfit
      external neb_rmsfit

      nfit = nion
      do i=1,nfit
         ifit(1,i) = i
         ifit(2,i) = i
         wfit(i)   = 1.0d0
      end do
      rms1 = neb_rmsfit(nfit,ifit,wfit,rion1,rion2)
c
c     superimpose the centroids of active atom pairs
c
      call neb_center(nion,rion1,rion2,
     >                      nfit,ifit,wfit,
     >                      xmid,ymid,zmid)

c
c     use a quaternion method to achieve the superposition
c
      call neb_quatfit(nion,rion1,rion2,nfit,ifit,wfit)
c
c     translate both coordinate sets so as to return
c     the first set to its original position
c
      do i = 1, nion
         rion1(1,i) = rion1(1,i) + xmid
         rion1(2,i) = rion1(2,i) + ymid
         rion1(3,i) = rion1(3,i) + zmid
      end do
      do i = 1, nion
         rion2(1,i) = rion2(1,i) + xmid
         rion2(2,i) = rion2(2,i) + ymid
         rion2(3,i) = rion2(3,i) + zmid
      end do
      rms2 = neb_rmsfit(nfit,ifit,wfit,rion1,rion2)

      return
      end

*     ***********************************************
*     *                                             *
*     *             neb_rmsfit                      *
*     *                                             *
*     ***********************************************
c     function neb_rmsfit  --  rms deviation for paired atoms
c
c     This routine computes the rms fit of two coordinate sets
c
      real*8 function neb_rmsfit(nfit,ifit,wfit,rion1,rion2)
      implicit none
      integer nfit,ifit(2,*)
      real*8 wfit(*)
      real*8 rion1(3,*),rion2(3,*)

      integer i,i1,i2
      real*8 rmsterm,rmsfit
      real*8 xr,yr,zr,dist2
      real*8 weight,norm
c
c     compute the rms fit over superimposed atom pairs
c
      rmsfit = 0.0d0
      norm = 0.0d0
      do i = 1, nfit
         i1 = ifit(1,i)
         i2 = ifit(2,i)
         weight = wfit(i)
         xr = rion1(1,i1) - rion2(1,i2)
         yr = rion1(2,i1) - rion2(2,i2)
         zr = rion1(3,i1) - rion2(3,i2)
         dist2 = xr**2 + yr**2 + zr**2
         norm = norm + weight
         rmsterm = dist2 * weight
         rmsfit = rmsfit + rmsterm
      end do
      neb_rmsfit = sqrt(rmsfit/norm)
      return
      end


*     ***********************************************
*     *                                             *
*     *             neb_center                      *
*     *                                             *
*     ***********************************************
c     subroutine neb_center  --  superimpose structure centroids 
c
c     This routine moves the weighted centroid of each coordinate
c     set to the origin during least squares superposition

      subroutine neb_center(nion,rion1,rion2,
     >                      nfit,ifit,wfit,
     >                      xmid,ymid,zmid)
      implicit none
      integer nion
      real*8 rion1(3,*),rion2(3,*)
      integer nfit,ifit(2,*)
      real*8  wfit(*)
      real*8 xmid,ymid,zmid

*     **** local variables ****
      integer i,k
      real*8 weight,norm
c
c
c     find the weighted centroid of the second
c     structure and translate it to the origin
c
      xmid = 0.0d0
      ymid = 0.0d0
      zmid = 0.0d0
      norm = 0.0d0
      do i = 1, nfit
         k = ifit(2,i)
         weight = wfit(i)
         xmid = xmid + rion2(1,k)*weight
         ymid = ymid + rion2(2,k)*weight
         zmid = zmid + rion2(3,k)*weight
         norm = norm + weight
      end do

      xmid = xmid / norm
      ymid = ymid / norm
      zmid = zmid / norm
      do i = 1, nion
         rion2(1,i) = rion2(1,i) - xmid
         rion2(2,i) = rion2(2,i) - ymid
         rion2(3,i) = rion2(3,i) - zmid
      end do
c
c     now repeat for the first structure, note
c     that this centroid position gets returned
c
      xmid = 0.0d0
      ymid = 0.0d0
      zmid = 0.0d0
      norm = 0.0d0
      do i = 1, nfit
         k = ifit(1,i)
         weight = wfit(i)
         xmid = xmid + rion1(1,k)*weight
         ymid = ymid + rion1(2,k)*weight
         zmid = zmid + rion1(3,k)*weight
         norm = norm + weight
      end do

      xmid = xmid / norm
      ymid = ymid / norm
      zmid = zmid / norm
      do i = 1, nion
         rion1(1,i) = rion1(1,i) - xmid
         rion1(2,i) = rion1(2,i) - ymid
         rion1(3,i) = rion1(3,i) - zmid
      end do

      return
      end


*     ***********************************************
*     *                                             *
*     *             neb_quatfit                     *
*     *                                             *
*     ***********************************************
c     subroutine quatfit  --  quaternion superposition of coords
c
c     This routine uses a quaternion-based method to achieve the best
c     fit superposition of two sets of coordinates
c
c     literature reference:
c
c     S. J. Kearsley, "An Algorithm for the Simultaneous Superposition
c     of a Structural Series", Journal of Computational Chemistry,
c     11, 1187-1192 (1990)
c
c     adapted from an original program written by David J. Heisterberg,
c     Ohio Supercomputer Center, Columbus, OH
c
      subroutine neb_quatfit(nion,rion1,rion2,nfit,ifit,wfit)
      implicit none
      integer nion
      real*8 rion1(3,*),rion2(3,*)
      integer nfit
      integer ifit(2,*)
      real*8 wfit(*)


      integer i,i1,i2,n1,n2
      real*8  weight,xrot,drot,zrot
      real*8  xxyx,xxyy,xxyz,xyyx,xyyy
      real*8  xyyz,xzyx,xzyy,xzyz
      real*8  rot(3,3),temp1(4),temp2(4)
      real*8  q(4),d(4),c(4,4),v(4,4)
c
c     build the upper triangle of the quadratic form matrix
c
      xxyx = 0.0d0
      xxyy = 0.0d0
      xxyz = 0.0d0
      xyyx = 0.0d0
      xyyy = 0.0d0
      xyyz = 0.0d0
      xzyx = 0.0d0
      xzyy = 0.0d0
      xzyz = 0.0d0
      do i = 1, nfit
         i1 = ifit(1,i)
         i2 = ifit(2,i)
         weight = wfit(i)
         xxyx = xxyx + weight*rion1(1,i1)*rion2(1,i2)
         xxyy = xxyy + weight*rion1(2,i1)*rion2(1,i2)
         xxyz = xxyz + weight*rion1(3,i1)*rion2(1,i2)
         xyyx = xyyx + weight*rion1(1,i1)*rion2(2,i2)
         xyyy = xyyy + weight*rion1(2,i1)*rion2(2,i2)
         xyyz = xyyz + weight*rion1(3,i1)*rion2(2,i2)
         xzyx = xzyx + weight*rion1(1,i1)*rion2(3,i2)
         xzyy = xzyy + weight*rion1(2,i1)*rion2(3,i2)
         xzyz = xzyz + weight*rion1(3,i1)*rion2(3,i2)
      end do
      c(1,1) = xxyx + xyyy + xzyz
      c(1,2) = xzyy - xyyz
      c(2,2) = xxyx - xyyy - xzyz
      c(1,3) = xxyz - xzyx
      c(2,3) = xxyy + xyyx
      c(3,3) = xyyy - xzyz - xxyx
      c(1,4) = xyyx - xxyy
      c(2,4) = xzyx + xxyz
      c(3,4) = xyyz + xzyy
      c(4,4) = xzyz - xxyx - xyyy
c
c     diagonalize the quadratic form matrix
c
      call neb_jacobi4(4,4,c,d,v,temp1,temp2)
c
c     extract the desired quaternion
c
      q(1) = v(1,4)
      q(2) = v(2,4)
      q(3) = v(3,4)
      q(4) = v(4,4)
c
c     assemble rotation matrix that superimposes the molecules
c
      rot(1,1) = q(1)**2 + q(2)**2 - q(3)**2 - q(4)**2
      rot(2,1) = 2.0d0 * (q(2) * q(3) - q(1) * q(4))
      rot(3,1) = 2.0d0 * (q(2) * q(4) + q(1) * q(3))
      rot(1,2) = 2.0d0 * (q(3) * q(2) + q(1) * q(4))
      rot(2,2) = q(1)**2 - q(2)**2 + q(3)**2 - q(4)**2
      rot(3,2) = 2.0d0 * (q(3) * q(4) - q(1) * q(2))
      rot(1,3) = 2.0d0 * (q(4) * q(2) - q(1) * q(3))
      rot(2,3) = 2.0d0 * (q(4) * q(3) + q(1) * q(2))
      rot(3,3) = q(1)**2 - q(2)**2 - q(3)**2 + q(4)**2
c
c     rotate second molecule to best fit with first molecule
c
      do i=1,nion
        xrot=rion2(1,i)*rot(1,1)+rion2(2,i)*rot(1,2)+rion2(3,i)*rot(1,3)
        drot=rion2(1,i)*rot(2,1)+rion2(2,i)*rot(2,2)+rion2(3,i)*rot(2,3)
        zrot=rion2(1,i)*rot(3,1)+rion2(2,i)*rot(3,2)+rion2(3,i)*rot(3,3)
        rion2(1,i) = xrot
        rion2(2,i) = drot
        rion2(3,i) = zrot
      end do

      return
      end


*     ***********************************************
*     *                                             *
*     *             neb_jacobi4                     *
*     *                                             *
*     ***********************************************
c     subroutine neb_jacobi4  --  jacobi matrix diagonalization
c
c
c     This routine performs a matrix diagonalization of a real
c     symmetric matrix by the method of Jacobi rotations
c
c     n    logical dimension of the matrix to be diagonalized
c     np   physical dimension of the matrix storage area
c     a    input with the matrix to be diagonalized; only
c             the upper triangle and diagonal are required
c     d    returned with the eigenvalues in ascending order
c     v    returned with the eigenvectors of the matrix
c     b    temporary work vector
c     z    temporary work vector
c
      subroutine neb_jacobi4(n,np,a,d,v,b,z)
      implicit none
      integer i,j,k,ip,iq,n,np,nrot,maxrot
      real*8  sm,tresh,s,c,t,theta,tau,h,g,p
      real*8  a(np,np),d(np),v(np,np),b(np),z(np)
c
c
c     setup and initialization
c
      maxrot = 100
      nrot = 0
      do ip = 1, n
         do iq = 1, n
            v(ip,iq) = 0.0d0
         end do
         v(ip,ip) = 1.0d0
      end do
      do ip = 1, n
         b(ip) = a(ip,ip)
         d(ip) = b(ip)
         z(ip) = 0.0d0
      end do
c
c     perform the jacobi rotations
c
      do i = 1, maxrot
         sm = 0.0d0
         do ip = 1, n-1
            do iq = ip+1, n
               sm = sm + abs(a(ip,iq))
            end do
         end do
         if (sm .eq. 0.0d0)  goto 10
         if (i .lt. 4) then
            tresh = 0.2d0*sm / n**2
         else
            tresh = 0.0d0
         end if
         do ip = 1, n-1
            do iq = ip+1, n
               g = 100.0d0 * abs(a(ip,iq))
               if (i.gt.4 .and. abs(d(ip))+g.eq.abs(d(ip))
     &                    .and. abs(d(iq))+g.eq.abs(d(iq))) then
                  a(ip,iq) = 0.0d0
               else if (abs(a(ip,iq)) .gt. tresh) then
                  h = d(iq) - d(ip)
                  if (abs(h)+g .eq. abs(h)) then
                     t = a(ip,iq) / h
                  else
                     theta = 0.5d0*h / a(ip,iq)
                     t = 1.0d0 / (abs(theta)+sqrt(1.0d0+theta**2))
                     if (theta .lt. 0.0d0)  t = -t
                  end if
                  c = 1.0d0 / sqrt(1.0d0+t**2)
                  s = t * c
                  tau = s / (1.0d0+c)
                  h = t * a(ip,iq)
                  z(ip) = z(ip) - h
                  z(iq) = z(iq) + h
                  d(ip) = d(ip) - h
                  d(iq) = d(iq) + h
                  a(ip,iq) = 0.0d0
                  do j = 1, ip-1
                     g = a(j,ip)
                     h = a(j,iq)
                     a(j,ip) = g - s*(h+g*tau)
                     a(j,iq) = h + s*(g-h*tau)
                  end do
                  do j = ip+1, iq-1
                     g = a(ip,j)
                     h = a(j,iq)
                     a(ip,j) = g - s*(h+g*tau)
                     a(j,iq) = h + s*(g-h*tau)
                  end do
                  do j = iq+1, n
                     g = a(ip,j)
                     h = a(iq,j)
                     a(ip,j) = g - s*(h+g*tau)
                     a(iq,j) = h + s*(g-h*tau)
                  end do
                  do j = 1, n
                     g = v(j,ip)
                     h = v(j,iq)
                     v(j,ip) = g - s*(h+g*tau)
                     v(j,iq) = h + s*(g-h*tau)
                  end do
                  nrot = nrot + 1
               end if
            end do
         end do
         do ip = 1, n
            b(ip) = b(ip) + z(ip)
            d(ip) = b(ip)
            z(ip) = 0.0d0
         end do
      end do
c
c     print warning if not converged
c
   10 continue
      if (nrot.eq.maxrot) then
         write(*,20)
   20    format (/,' JACOBI4 -- Matrix Diagonalization not Converged')
      end if
c
c     sort the eigenvalues and vectors
c
      do i = 1, n-1
         k = i
         p = d(i)
         do j = i+1, n
            if (d(j) .lt. p) then
               k = j
               p = d(j)
            end if
         end do
         if (k .ne. i) then
            d(k) = d(i)
            d(i) = p
            do j = 1, n
               p = v(j,i)
               v(j,i) = v(j,k)
               v(j,k) = p
            end do
         end if
      end do
      return
      end




*     ********************************
*     *                              *
*     *      neb_linesearch_init     *
*     *                              *
*     ********************************
      subroutine neb_linesearch_init()
      implicit none

*     **** neblinesearch_counter common block ****
      integer counter
      common / neblinesearch_counter / counter

      counter = 0

      return
      end

*     ********************************
*     *                              *
*     *      neb_linesearch_count    *
*     *                              *
*     ********************************
      integer function neb_linesearch_count()
      implicit none

*     **** neblinesearch_counter common block ****
      integer counter
      common / neblinesearch_counter / counter
 
      neb_linesearch_count = counter
      return
      end

*     ********************************
*     *                              *
*     *         neb_linesearch	     *
*     *                              *
*     ********************************

      real*8 function neb_linesearch(bead_list,kbeads,t0,f0,deltat,
     >                           tolerance,tmin,deltaE,
     >                           stoptype)
      implicit none
      character*(*) bead_list
      real*8 kbeads
      real*8 t0,f0
      real*8 deltat
      real*8   tolerance
      real*8   tmin
      real*8   deltaE
      integer  stoptype


*     **** local variables ****
      integer maxiter,iteration
      parameter (maxiter=8)
      logical secant,notfinished
      integer indx(3)
      real*8  t(3)
      real*8  f(3)
      real*8  t_last, f_last
      real*8  t_first,f_first
      real*8 up,down,fmin,deltaf

*     **** neblinesearch_counter common block ****
      integer counter
      common / neblinesearch_counter / counter

*     **** external functions ****
      real*8   neb_line_energy
      external neb_line_energy
   
      counter = counter + 1

      secant = .true.

      t(1) = t0
      f(1)  = f0
      t_last = t(1)
      f_last = f(1)
      t_first = t(1)
      f_first = f(1)

      f(2)  =  neb_line_energy(bead_list,kbeads,t(1)+deltat,0)
      
      iteration = 1
*     **** make sure that f2 < f1 ******
      do while ((f(2).gt.f(1)).and.(iteration.le.maxiter))
        deltat = 0.5d0*deltat
        f(2)  =  neb_line_energy(bead_list,kbeads,t(1)+deltat,0)
        iteration = iteration + 1
      end do
      t(2) = t(1) + deltat
      t_last = t(2)
      f_last = f(2)


*     **** use secant method to generate f(3) *****
      deltat = -f(1)*(t(2)-t(1))/(f(2)-f(1))
      t(3)   = t(1) + deltat
      f(3)   =  neb_line_energy(bead_list,kbeads,t(3),0)
      iteration = iteration + 1
      t_last = t(3)
      f_last = f(3)

*     **** sort the function values ****
      call neb_Order_Values(3,f,indx)

      deltaf = f(indx(2)) - f(indx(1))


      if (stoptype.eq.1) then
        notfinished = (dabs(deltaf).gt.tolerance)
     >                .and.(iteration.le.maxiter)
      else
        notfinished = (dabs(f(indx(1))/f_first).gt.tolerance)
     >                .and.(iteration.le.maxiter)
      end if


      do while (notfinished) 
      

*       **** use secant interpolation to generate tmin ***
        if (secant) then

          call util_flush(6)
          deltat = -f(indx(1))
     >           *(t(indx(2))-t(indx(1)))
     >           /(f(indx(2))-f(indx(1)))
          tmin = t(indx(1)) + deltat
          call util_flush(6)
          fmin  =  neb_line_energy(bead_list,kbeads,tmin,0)
          iteration = iteration + 1
          t_last = tmin
          f_last = fmin

*         **** finish using secant method ****
          if (fmin.ge.f(indx(1))) then
            secant = .false.
            if (fmin.lt.f(indx(3))) then
              t(indx(3))  = tmin
              f(indx(3))  = fmin
              call neb_Order_Values(3,f,indx)
            end if
          end if


        end if 

*       **** use quadradic interpolation to generate tmin ***
        if (.not.secant) then
          up  = (t(2)*t(2) - t(3)*t(3))*f(1)
     >        + (t(3)*t(3) - t(1)*t(1))*f(2)
     >        + (t(1)*t(1) - t(2)*t(2))*f(3)
          down = (t(2) - t(3))*f(1)
     >         + (t(3) - t(1))*f(2)
     >         + (t(1) - t(2))*f(3)

*         **** check added by E.Apra ****
          if(abs(down).gt.tolerance**2) then
             tmin = 0.5d0*up/down
             fmin  =  neb_line_energy(bead_list,kbeads,tmin,0)
             iteration = iteration + 1
             t_last = tmin
             f_last = fmin

*         **** parabola fit failed - exit loop ****
          else
             tmin=t(indx(3))
             fmin=f(indx(3))+tolerance
             iteration = maxiter+1
          endif

        end if


*       **** tolerance check and replacement ****
        if (fmin.lt.f(indx(3))) then
           t(indx(3))  = tmin
           f(indx(3))  = fmin
           call neb_Order_Values(3,f,indx)
           deltaf = f(indx(2)) - f(indx(1))
        else
           deltaf=0.0d0
        end if

        if (stoptype.eq.1) then
          notfinished = (dabs(deltaf).gt.tolerance)
     >                .and.(iteration.le.maxiter)
        else
          notfinished = (dabs(f(indx(1))/f_first).gt.tolerance)
     >                .and.(iteration.le.maxiter)
        end if
           
      end do

*     **** make sure that tmin is last functions call *****
      tmin = t(indx(1))
      fmin = f(indx(1))
      if (tmin.ne.t_last) 
     >   f_last = neb_line_energy(bead_list,kbeads,tmin,0)

      deltaE = (fmin-f_first)
  
      neb_linesearch = fmin
      return
      end

*     *****************************
*     *                           *
*     *     neb_Order_Values      *
*     *                           *
*     *****************************
*
*   this subroutine makes f(indx(1)) < f(indx(2)) < f(indx(3)) < ....
*   via a bubble sort
*   
*   Entry - n,f
*   Exit - returns indx
*

      subroutine neb_Order_Values(n,f,indx)
      implicit none
      integer n
      real*8 f(*)
      integer indx(*)

*     ****** local variables *****
      integer i,j,idum

      do i=1,n
         indx(i) = i
      end do
      do i=1,n-1
        do j=i+1,n
           if (f(indx(j)).lt.f(indx(i))) then
              idum    = indx(i)
              indx(i) = indx(j)
              indx(j) = idum
           end if
         end do
       end do
            
      return 
      end


*     *****************************
*     *                           *
*     *     neb_linesearch_start  *
*     *                           *
*     *****************************
      subroutine neb_linesearch_start(ng_in,c0_in,s_in,
     >                                m_in,cs_in,gs_in,
     >                                nbeads)
      implicit none
      integer ng_in,m_in
      real*8 c0_in(*)
      real*8 s_in(*)
      real*8 cs_in(*)
      real*8 gs_in(*)
      integer nbeads
      
#include "mafdecls.fh"

*     *** neb_linesearch ***
      integer cs(2),gs(2),c1(2),c0(2),s(2),g1(2),t1(2),e1(2),m,ng
      common / neblinesearch / cs,gs,c1,c0,s,g1,t1,e1,m,ng

*     **** local variables ****
      logical value

      m = m_in + 1
      ng = ng_in

      value  = MA_push_get(mt_dbl,m*ng,'lcs',cs(2),cs(1))
      value  = value.and.MA_push_get(mt_dbl,m*ng,'lgs',gs(2),gs(1))
      value  = value.and.MA_push_get(mt_dbl,ng,'lc1',c1(2),c1(1))
      value  = value.and.MA_push_get(mt_dbl,ng,'lc0',c0(2),c0(1))
      value  = value.and.MA_push_get(mt_dbl,ng,'ls',s(2),s(1))
      value  = value.and.MA_push_get(mt_dbl,ng,'lg1',g1(2),g1(1))
      value  = value.and.MA_push_get(mt_dbl,ng,'lt1',t1(2),t1(1))
      value  = value.and.MA_push_get(mt_dbl,nbeads,'le1',e1(2),e1(1))
      if (.not.value)
     >   call errquit('neb_linesearch_start:increase stack memory',1,0)

      call dcopy(m_in*ng,cs_in,1,dbl_mb(cs(1)),1)
      call dcopy(m_in*ng,gs_in,1,dbl_mb(gs(1)),1)
      call dcopy(ng,c0_in,1,dbl_mb(c0(1)),1)
      call dcopy(ng,s_in,1,dbl_mb(s(1)),1)

      return
      end

*     *****************************
*     *                           *
*     *     neb_linesearch_end    *
*     *                           *
*     *****************************
      subroutine neb_linesearch_end()
      implicit none

#include "mafdecls.fh"

*     *** neb_linesearch ***
      integer cs(2),gs(2),c1(2),c0(2),s(2),g1(2),t1(2),e1(2),m,ng
      common / neblinesearch / cs,gs,c1,c0,s,g1,t1,e1,m,ng

*     **** local variables ****
      logical value

      value =           MA_pop_stack(e1(2))
      value = value.and.MA_pop_stack(t1(2))
      value = value.and.MA_pop_stack(g1(2))
      value = value.and.MA_pop_stack(s(2))
      value = value.and.MA_pop_stack(c0(2))
      value = value.and.MA_pop_stack(c1(2))
      value = value.and.MA_pop_stack(gs(2))
      value = value.and.MA_pop_stack(cs(2))
      if (.not.value)
     >   call errquit('neb_linesearch_end: popping stack',1,0)
  
      return
      end

*     *****************************
*     *                           *
*     *     neb_line_energy	  *
*     *                           *
*     *****************************
      real*8 function neb_line_energy(bead_list,kbeads,alpha,opt)
      implicit none
      character*(*) bead_list
      real*8 kbeads,alpha
      integer opt

#include "mafdecls.fh"

*     *** neb_linesearch ***
      integer cs(2),gs(2),c1(2),c0(2),s(2),g1(2),t1(2),e1(2),m,ng
      common / neblinesearch / cs,gs,c1,c0,s,g1,t1,e1,m,ng

*     **** local variables ****
      integer shift
      real*8  ee

*     **** external functions ****
      logical  task_gradient
      external task_gradient
      real*8   ddot
      external ddot

       ee = 0.0d0
       call util_flush(6)
       call neb_move(ng,
     >               (-alpha),
     >               dbl_mb(c0(1)),
     >               dbl_mb(c1(1)),
     >               dbl_mb(s(1)))

       if (opt.eq.0) then
           call neb_coords_set(bead_list,dbl_mb(c1(1)))
           call runmid_bead_list(bead_list,task_gradient,.true.,.true.)
           call neb_energies_get(bead_list,dbl_mb(e1(1)))
           call neb_gradient_get(bead_list,kbeads,
     >                          dbl_mb(c1(1)),
     >                          dbl_mb(e1(1)),
     >                          dbl_mb(t1(1)),
     >                          dbl_mb(g1(1)))
           shift = (m-1)*ng
           call dcopy(ng,dbl_mb(c1(1)),1,dbl_mb(cs(1)+shift),1)
           call dcopy(ng,dbl_mb(g1(1)),1,dbl_mb(gs(1)+shift),1)
           call neb_lmbfgs(ng,m,
     >                     dbl_mb(cs(1)),
     >                     dbl_mb(gs(1)),
     >                     dbl_mb(s(1)))
           ee = ddot(ng,dbl_mb(g1(1)),1,dbl_mb(s(1)),1)
           call util_flush(6)
      end if


      neb_line_energy = ee
      return
      end
