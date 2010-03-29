
c     **********************************************
c     *                                            *
c     *            meta_initialize                 *
c     *                                            *
c     **********************************************
      subroutine meta_initialize(rtdb)
      implicit none
      integer rtdb

#include "mafdecls.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "metadynamics.fh"

*     **** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)
      logical value
      integer i,j
      real*8 x,y,z,r,rmax,dr,pi,da
      character*80 rtdb_name

*     **** external functions ****
      real*8      lattice_unita
      character*7 c_index_name
      external    lattice_unita
      external    c_index_name

      call Parallel_taskid(taskid)
      metaprintcount = 0
      rtdb_name = 'metadynamics_print'
      if (.not.rtdb_get(rtdb,rtdb_name,mt_int,1,maxmetaprintcount)) 
     >   maxmetaprintcount=100

      rtdb_name = 'metadynamics_nbond'
      if (.not.rtdb_get(rtdb,rtdb_name,mt_int,1,bnmeta)) bnmeta = 0

      rtdb_name = 'metadynamics_nangle'
      if (.not.rtdb_get(rtdb,rtdb_name,mt_int,1,anmeta)) anmeta = 0


      metafound = (bnmeta.gt.0).or.(anmeta.gt.0)
      if (metafound) then
         rtdb_name = 'metadynamics_update'
         if (.not.rtdb_get(rtdb,rtdb_name,mt_int,1,maxmetacount)) 
     >      maxmetacount = 1
         metacount = 0

         value = .true.
         if (bnmeta.gt.0) then
            bnxmeta = 501
            value = value.and.
     >              MA_alloc_get(mt_dbl,bnxmeta,'bxmeta',
     >                        bxmeta(2),bxmeta(1))
            value = value.and.
     >              MA_alloc_get(mt_dbl,bnmeta*bnxmeta,'bymeta',
     >                           bymeta(2),bymeta(1))
            value = value.and.
     >              MA_alloc_get(mt_dbl,bnmeta*bnxmeta,'by2meta',
     >                           by2meta(2),by2meta(1))

*           **** load old spline fits ****
            call dcopy(bnmeta*bnxmeta,0.0d0,0,dbl_mb(bymeta(1)),1)
            call dcopy(bnmeta*bnxmeta,0.0d0,0,dbl_mb(by2meta(1)),1)
            rtdb_name = 'metadynamics_bymeta'
            if (.not.rtdb_get(rtdb,rtdb_name,mt_dbl,bnmeta*bnxmeta,
     >                       dbl_mb(bymeta(1)))) then
               call dcopy(bnmeta*bnxmeta,0.0d0,0,dbl_mb(bymeta(1)),1)
            else
               if (taskid.eq.MASTER) 
     >            write(*,*) "  ... reading bond metadynamics data"
            end if
            rtdb_name = 'metadynamics_by2meta'
            if (.not.rtdb_get(rtdb,rtdb_name,mt_dbl,bnmeta*bnxmeta,
     >                       dbl_mb(by2meta(1)))) then
               call dcopy(bnmeta*bnxmeta,0.0d0,0,dbl_mb(by2meta(1)),1)
            else
               if (taskid.eq.MASTER) 
     >            write(*,*) "  ... reading more bond metadynamics data"
            end if

            value = value.and.
     >              MA_alloc_get(mt_dbl,2*bnmeta,'bwsigmameta',
     >                           bwsigmameta(2),bwsigmameta(1))
            value = value.and.
     >              MA_alloc_get(mt_int,2*bnmeta,'bindxmeta',
     >                           bindxmeta(2),bindxmeta(1))
            do j=1,bnmeta
               rtdb_name = 
     >         'metadynamics_bond_index:'//c_index_name(j)
                value = value.and.
     >             rtdb_get(rtdb,rtdb_name,mt_int,2,
     >                      int_mb(bindxmeta(1)+(j-1)*2))
               rtdb_name = 
     >        'metadynamics_bond_ksigma:'//c_index_name(j)
               value = value.and.
     >             rtdb_get(rtdb,rtdb_name,mt_dbl,2,
     >                      dbl_mb(bwsigmameta(1)+(j-1)*2))
            end do

            x = lattice_unita(1,1)
            y = lattice_unita(2,1)
            z = lattice_unita(3,1)
            r = dsqrt(x*x+y*y+z*z)
            rmax = r
            x = lattice_unita(1,2)
            y = lattice_unita(2,2)
            z = lattice_unita(3,2)
            r = dsqrt(x*x+y*y+z*z)
            if (r.gt.rmax) rmax = r
            x = lattice_unita(1,3)
            y = lattice_unita(2,3)
            z = lattice_unita(3,3)
            r = dsqrt(x*x+y*y+z*z)
            if (r.gt.rmax) rmax = r
            dr = rmax/dble(bnxmeta-1)
            do i=1,bnxmeta
               r = (i-1)*dr
               dbl_mb(bxmeta(1)+i-1) = r
            end do
         end if

         if (anmeta.gt.0) then
            anxmeta = 501
            value = value.and.
     >              MA_alloc_get(mt_dbl,anxmeta,'axmeta',
     >                           axmeta(2),axmeta(1))
            value = value.and.
     >              MA_alloc_get(mt_dbl,anmeta*anxmeta,'aymeta',
     >                           aymeta(2),aymeta(1))
            value = value.and.
     >              MA_alloc_get(mt_dbl,anmeta*anxmeta,'ay2meta',
     >                           ay2meta(2),ay2meta(1))

*           **** load old spline fits ****
            call dcopy(anmeta*anxmeta,0.0d0,0,dbl_mb(aymeta(1)),1)
            call dcopy(anmeta*anxmeta,0.0d0,0,dbl_mb(ay2meta(1)),1)
            rtdb_name = 'metadynamics_aymeta'
            if (.not.rtdb_get(rtdb,rtdb_name,mt_dbl,anmeta*anxmeta,
     >                       dbl_mb(aymeta(1)))) then
               call dcopy(anmeta*anxmeta,0.0d0,0,dbl_mb(aymeta(1)),1)
            else
               if (taskid.eq.MASTER) 
     >           write(*,*) "  ... reading angle metadynamics data"
            end if
            rtdb_name = 'metadynamics_ay2meta'
            if (.not.rtdb_get(rtdb,rtdb_name,mt_dbl,anmeta*anxmeta,
     >                       dbl_mb(ay2meta(1)))) then
               call dcopy(anmeta*anxmeta,0.0d0,0,dbl_mb(ay2meta(1)),1)
            else
               if (taskid.eq.MASTER) 
     >           write(*,*) "  ... reading more angle metadynamics data"
            end if


            value = value.and.
     >              MA_alloc_get(mt_dbl,2*anmeta,'awsigmameta',
     >                           awsigmameta(2),awsigmameta(1))
            value = value.and.
     >              MA_alloc_get(mt_int,3*anmeta,'aindxmeta',
     >                           aindxmeta(2),aindxmeta(1))
            do j=1,anmeta
               rtdb_name =
     >         'metadynamics_angle_index:'//c_index_name(j)
                value = value.and.
     >             rtdb_get(rtdb,rtdb_name,mt_int,3,
     >                      int_mb(aindxmeta(1)+(j-1)*3))
               rtdb_name = 
     >        'metadynamics_angle_ksigma:'//c_index_name(j)
               value = value.and.
     >             rtdb_get(rtdb,rtdb_name,mt_dbl,2,
     >                      dbl_mb(awsigmameta(1)+(j-1)*2))
            end do



            pi = 4.0d0*datan(1.0d0)
            da = pi/dble(anxmeta-1)
            do i=1,anxmeta
               r = (i-1)*da
               dbl_mb(axmeta(1)+i-1) = r
            end do
         end if
         value = value.and.
     >           MA_alloc_get(mt_dbl,501,'tmeta',
     >                        tmeta(2),tmeta(1))
         if (.not.value)
     >    call errquit('cannot allocate heap memory for metadynamics',0,
     >       MA_ERR)

*     **** write out header info ****
      call Parallel_taskid(taskid)
      if (taskid.eq.MASTER) then
         write(*,*) 
         write(*,*)  "MetaDynamics parameters:"
         do i=1,bnmeta
             write(*,'(A,5x,2I4,4x,A,F11.6,4x,A,F11.6)')
     >       "   - Bond Parameters =  ",
     >                  (int_mb(bindxmeta(1)+2*(i-1)+j-1),j=1,2),
     >                  'w=',
     >                  dbl_mb(bwsigmameta(1)+(i-1)*2),
     >                  'sigma=',
     >                  dbl_mb(bwsigmameta(1)+(i-1)*2+1)

         end do
         do i=1,anmeta
             write(*,'(A,3I4,4x,A,F11.6,4x,A,F11.6)')
     >       "   - Angle Parameters =  ",
     >                  (int_mb(aindxmeta(1)+3*(i-1)+j-1),j=1,3),
     >                  'w=',
     >                  dbl_mb(awsigmameta(1)+2*(i-1)),
     >                  'sigma=',
     >                  dbl_mb(awsigmameta(1)+2*(i-1)+1)

         end do

      end if
      
      end if
      return
      end

c     **********************************************
c     *                                            *
c     *            meta_finalize                   *
c     *                                            *
c     **********************************************
      subroutine meta_finalize(rtdb)
      implicit none
      integer rtdb

#include "mafdecls.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "metadynamics.fh"

*     **** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)

      logical value
      integer i,j,yshift,k
      real*8 r
      character*80 filename,rtdb_name
      character*255 full_filename

*     **** external functions ****
      character*7 c_index_name
      external    c_index_name

      call Parallel_taskid(taskid)

      if (metafound) then

*        **** print out potentials ******
         call meta_print_potentials()

         value = .true.
         if (bnmeta.gt.0) then

*           **** write out current spline fits ****
            rtdb_name = 'metadynamics_bymeta'
            value = value.and.rtdb_put(rtdb,rtdb_name,mt_dbl,
     >                        bnmeta*bnxmeta,dbl_mb(bymeta(1)))
            rtdb_name = 'metadynamics_by2meta'
            value = value.and.rtdb_put(rtdb,rtdb_name,mt_dbl,
     >                        bnmeta*bnxmeta,dbl_mb(by2meta(1)))

            value = value.and.MA_free_heap(bxmeta(2))
            value = value.and.MA_free_heap(bymeta(2))
            value = value.and.MA_free_heap(by2meta(2))
            value = value.and.MA_free_heap(bwsigmameta(2))
            value = value.and.MA_free_heap(bindxmeta(2))
         end if
         if (anmeta.gt.0) then

*           **** write out current spline fits ****
            rtdb_name = 'metadynamics_aymeta'
            value = value.and.rtdb_put(rtdb,rtdb_name,mt_dbl,
     >                        anmeta*anxmeta,dbl_mb(aymeta(1)))
            rtdb_name = 'metadynamics_ay2meta'
            value = value.and.rtdb_put(rtdb,rtdb_name,mt_dbl,
     >                        anmeta*anxmeta,dbl_mb(ay2meta(1)))

            value = value.and.MA_free_heap(axmeta(2))
            value = value.and.MA_free_heap(aymeta(2))
            value = value.and.MA_free_heap(ay2meta(2))
            value = value.and.MA_free_heap(awsigmameta(2))
            value = value.and.MA_free_heap(aindxmeta(2))
         end if
         value = value.and.MA_free_heap(tmeta(2))
      end if

      return
      end

c     **********************************************
c     *                                            *
c     *            meta_collective                 *
c     *                                            *
c     **********************************************
      real*8 function meta_collective(ntype,indx)
      implicit none
      integer ntype,indx(*)

*     **** local variables ****
      real*8 x1,y1,z1,x2,y2,z2,x3,y3,z3,dx,dy,dz,r,r1,r3
      real*8 dx1,dy1,dz1,dx3,dy3,dz3,theta

*     **** external functions ****
      real*8   ion_rion
      external ion_rion

      r = 0.0d0
c     **** bond distance ***
      if (ntype.eq.1) then
         x1 = ion_rion(1,indx(1))
         y1 = ion_rion(2,indx(1))
         z1 = ion_rion(3,indx(1))
         x2 = ion_rion(1,indx(2))
         y2 = ion_rion(2,indx(2))
         z2 = ion_rion(3,indx(2))
         dx = x1-x2
         dy = y1-y2
         dz = z1-z2
         call lattice_min_difference(dx,dy,dz)
         r  = dsqrt(dx**2 + dy**2 + dz**2)

c     **** bond angle ***
      else if (ntype.eq.2) then
         x1 = ion_rion(1,indx(1))
         y1 = ion_rion(2,indx(1))
         z1 = ion_rion(3,indx(1))
         x2 = ion_rion(1,indx(2))
         y2 = ion_rion(2,indx(2))
         z2 = ion_rion(3,indx(2))
         x3 = ion_rion(1,indx(3))
         y3 = ion_rion(2,indx(3))
         z3 = ion_rion(3,indx(3))
         dx1 = x1-x2
         dy1 = y1-y2
         dz1 = z1-z2
         call lattice_min_difference(dx1,dy1,dz1)
         r1  = dsqrt(dx1**2 + dy1**2 + dz1**2)
         dx3 = x3-x2
         dy3 = y3-y2
         dz3 = z3-z2
         call lattice_min_difference(dx3,dy3,dz3)
         r3  = dsqrt(dx3**2 + dy3**2 + dz3**2)
         theta = (dx1*dx3 + dy1*dy3 + dz1*dz3)/(r1*r3)
         if (theta.gt.1.0d0)  theta = 1.0d0
         if (theta.lt.-1.0d0) theta = -1.0d0
         theta = dacos(theta)

c     **** bond dihedral ****

      end if

c     **** coordination number ****

      meta_collective = r
      return
      end

c     **********************************************
c     *                                            *
c     *            meta_collective_force           *
c     *                                            *
c     **********************************************

      subroutine meta_collective_force(ntype,indx,dv,fmeta)
      implicit none
      integer ntype,indx(*)
      real*8 dv,fmeta(3,*)

*     **** local variables ****
      real*8 x1,y1,z1,r1,vx1,vx2,vy1,vy2,vz1,vz2
      real*8 x2,y2,z2,aa,a11,a12,a22,r
      real*8 x3,y3,z3,r3,ctheta,stheta,denom
      real*8 dx,dy,dz,dx1,dy1,dz1,dx3,dy3,dz3

*     **** external functions ****
      real*8   ion_rion
      external ion_rion


c     **** bond distance ***
      if (ntype.eq.1) then
         x1 = ion_rion(1,indx(1))
         y1 = ion_rion(2,indx(1))
         z1 = ion_rion(3,indx(1))
         x2 = ion_rion(1,indx(2))
         y2 = ion_rion(2,indx(2))
         z2 = ion_rion(3,indx(2))
         dx = x1-x2
         dy = y1-y2
         dz = z1-z2
         call lattice_min_difference(dx,dy,dz)
         r  = dsqrt(dx**2 + dy**2 + dz**2)

         fmeta(1,indx(1)) = fmeta(1,indx(1)) - (dx/r)*dv
         fmeta(2,indx(1)) = fmeta(2,indx(1)) - (dy/r)*dv
         fmeta(3,indx(1)) = fmeta(3,indx(1)) - (dz/r)*dv
         fmeta(1,indx(2)) = fmeta(1,indx(2)) + (dx/r)*dv
         fmeta(2,indx(2)) = fmeta(2,indx(2)) + (dy/r)*dv
         fmeta(3,indx(2)) = fmeta(3,indx(2)) + (dz/r)*dv

c     **** bond angle ***
      else if (ntype.eq.2) then
         x1 = ion_rion(1,indx(1))
         y1 = ion_rion(2,indx(1))
         z1 = ion_rion(3,indx(1))
         x2 = ion_rion(1,indx(2))
         y2 = ion_rion(2,indx(2))
         z2 = ion_rion(3,indx(2))
         x3 = ion_rion(1,indx(3))
         y3 = ion_rion(2,indx(3))
         z3 = ion_rion(3,indx(3))
         dx1 = x1-x2
         dy1 = y1-y2
         dz1 = z1-z2
         call lattice_min_difference(dx1,dy1,dz1)
         r1  = dsqrt(dx1**2 + dy1**2 + dz1**2)
         dx3 = x3-x2
         dy3 = y3-y2
         dz3 = z3-z2
         call lattice_min_difference(dx3,dy3,dz3)
         r3  = dsqrt(dx3**2 + dy3**2 + dz3**2)
         denom = r1*r3
         if (denom.gt.1.0d-11) then
            ctheta = (dx1*dx3 + dy1*dy3 + dz1*dz3)/(denom)
            if (ctheta.gt.1.0d0)  ctheta =  1.0d0
            if (ctheta.lt.-1.0d0) ctheta = -1.0d0
            stheta = dsqrt(1.0d0-ctheta*ctheta)
            if (stheta.lt.0.001d0) stheta = 0.001d0
            stheta = 1.0d0/stheta

            aa  =  dv*stheta
            a11 =  aa*ctheta/r1
            a12 = -aa/(denom)
            a22 =  aa*ctheta/r3

            vx1 = a11*dx1 + a12*dx3
            vx2 = a22*dx3 + a12*dx1

            vy1 = a11*dy1 + a12*dy3
            vy2 = a22*dy3 + a12*dy1

            vz1 = a11*dz1 + a12*dz3
            vz2 = a22*dz3 + a12*dz1

            fmeta(1,indx(1)) = fmeta(1,indx(1)) - vx1
            fmeta(2,indx(1)) = fmeta(2,indx(1)) - vy1
            fmeta(3,indx(1)) = fmeta(3,indx(1)) - vz1

            fmeta(1,indx(2)) = fmeta(1,indx(2)) + vx1 + vx2
            fmeta(2,indx(2)) = fmeta(2,indx(2)) + vy1 + vy2
            fmeta(3,indx(2)) = fmeta(3,indx(2)) + vz1 + vz2

            fmeta(1,indx(3)) = fmeta(1,indx(3)) - vx2
            fmeta(2,indx(3)) = fmeta(2,indx(3)) - vy2
            fmeta(3,indx(3)) = fmeta(3,indx(3)) - vz2
         end if




c     **** bond dihedral ****

      end if



      return
      end

*     ***********************************************
*     *                                             *
*     *              meta_update                    *
*     *                                             *
*     ***********************************************
      subroutine meta_update()
      implicit none

#include "mafdecls.fh"
#include "errquit.fh"
#include "metadynamics.fh"

*     **** local variables ****
      integer i,j,yshift,y2shift
      real*8  r0,r,dr,w,sigma,yp1,ypn

*     **** external functions ****
      real*8   meta_collective
      external meta_collective

      if (metafound) then

         metacount = metacount + 1
         if (metacount.ge.maxmetacount) then
         metacount = 0

         if (bnmeta.gt.0) then
            dr = dbl_mb(bxmeta(1)+1) - dbl_mb(bxmeta(1))
            do j=1,bnmeta
               r0    = meta_collective(1,int_mb(bindxmeta(1)+(j-1)*2))
               w     = dbl_mb(bwsigmameta(1)+(j-1)*2)
               sigma = dbl_mb(bwsigmameta(1)+(j-1)*2+1)
               yshift  = bymeta(1) +(j-1)*bnxmeta
               y2shift = by2meta(1)+(j-1)*bnxmeta
               do i=1,bnxmeta
                  r = dbl_mb(bxmeta(1)+i-1)-r0
                  dbl_mb(yshift+i-1) = dbl_mb(yshift+i-1) 
     >                               + w*dexp(-(r/sigma)**2)
               end do

               !**** five point formula ***
               yp1 = ( -50.0d0*dbl_mb(yshift)
     >                + 96.0d0*dbl_mb(yshift+1)
     >                - 72.0d0*dbl_mb(yshift+2)
     >                + 32.0d0*dbl_mb(yshift+3)
     >                -  6.0d0*dbl_mb(yshift+4))/(24.0d0*dr)
               ypn = ( +50.0d0*dbl_mb(yshift+bnxmeta-5)
     >                - 96.0d0*dbl_mb(yshift+bnxmeta-4)
     >                + 72.0d0*dbl_mb(yshift+bnxmeta-3)
     >                - 32.0d0*dbl_mb(yshift+bnxmeta-2)
     >                +  6.0d0*dbl_mb(yshift+bnxmeta-1))/(24.0d0*dr)
               call nwpw_spline(dbl_mb(bxmeta(1)),
     >                 dbl_mb(yshift),bnxmeta,yp1,ypn,
     >                 dbl_mb(y2shift),dbl_mb(tmeta(1)))

            end do
         end if

         if (anmeta.gt.0) then
            dr = dbl_mb(axmeta(1)+1) - dbl_mb(axmeta(1))
            do j=1,anmeta
               r0    = meta_collective(2,int_mb(aindxmeta(1)+(j-1)*3))
               w     = dbl_mb(awsigmameta(1)+(j-1)*2)
               sigma = dbl_mb(awsigmameta(1)+(j-1)*2+1)
               yshift  = aymeta(1) +(j-1)*anxmeta
               y2shift = ay2meta(1)+(j-1)*anxmeta
               do i=1,anxmeta
                  r = dbl_mb(axmeta(1)+i-1)-r0
                  dbl_mb(yshift+i-1) = dbl_mb(yshift+i-1)
     >                               + w*dexp(-(r/sigma)**2)
               end do

               !**** five point formula ***
               yp1 = ( -50.0d0*dbl_mb(yshift)
     >                + 96.0d0*dbl_mb(yshift+1)
     >                - 72.0d0*dbl_mb(yshift+2)
     >                + 32.0d0*dbl_mb(yshift+3)
     >                -  6.0d0*dbl_mb(yshift+4))/(24.0d0*dr)
               ypn = ( +50.0d0*dbl_mb(yshift+anxmeta-5)
     >                - 96.0d0*dbl_mb(yshift+anxmeta-4)
     >                + 72.0d0*dbl_mb(yshift+anxmeta-3)
     >                - 32.0d0*dbl_mb(yshift+anxmeta-2)
     >                +  6.0d0*dbl_mb(yshift+anxmeta-1))/(24.0d0*dr)
               call nwpw_spline(dbl_mb(axmeta(1)),
     >                 dbl_mb(yshift),anxmeta,yp1,ypn,
     >                 dbl_mb(y2shift),dbl_mb(tmeta(1)))
            end do
         end if

         end if

         metaprintcount = metaprintcount + 1
         if (metaprintcount.ge.maxmetaprintcount) then
            metaprintcount = 0
            call meta_print_potentials()
         end if

      end if
 
      return
      end

*     ***********************************************
*     *                                             *
*     *              meta_force                     *
*     *                                             *
*     ***********************************************
      subroutine meta_force(fmeta)
      implicit none
      real*8 fmeta(3,*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "metadynamics.fh"

*     **** local variables ****
      integer j,nr,yshift,y2shift
      real*8 r,dr,sigma,w,dv

*     **** external functions ****
      real*8   nwpw_dsplint,meta_collective
      external nwpw_dsplint,meta_collective
      
      if (metafound) then

         if (bnmeta.gt.0) then
            dr = dbl_mb(bxmeta(1)+1) - dbl_mb(bxmeta(1))
            do j=1,bnmeta
               r     = meta_collective(1,int_mb(bindxmeta(1)+(j-1)*2))
               w     = dbl_mb(bwsigmameta(1)+(j-1)*2)
               sigma = dbl_mb(bwsigmameta(1)+(j-1)*2+1)
               yshift  = bymeta(1) +(j-1)*bnxmeta
               y2shift = by2meta(1)+(j-1)*bnxmeta
               nr = r/dr + 1
               dv = nwpw_dsplint(dbl_mb(bxmeta(1)),
     >                  dbl_mb(yshift),dbl_mb(y2shift),bnmeta,nr,r)
               call meta_collective_force(1,
     >                  int_mb(bindxmeta(1)+(j-1)*2),dv,fmeta)
            end do
         end if

         if (anmeta.gt.0) then
            dr = dbl_mb(axmeta(1)+1) - dbl_mb(axmeta(1))
            do j=1,anmeta
               r     = meta_collective(1,int_mb(aindxmeta(1)+(j-1)*3))
               w     = dbl_mb(awsigmameta(1)+(j-1)*2)
               sigma = dbl_mb(awsigmameta(1)+(j-1)*2+1)
               yshift  = aymeta(1) +(j-1)*anxmeta
               y2shift = ay2meta(1)+(j-1)*anxmeta
               nr = r/dr + 1
               dv = nwpw_dsplint(dbl_mb(axmeta(1)),
     >                  dbl_mb(yshift),dbl_mb(y2shift),anmeta,nr,r)
               call meta_collective_force(1,
     >                  int_mb(aindxmeta(1)+(j-1)*3),dv,fmeta)
            end do

         end if
      end if


      return
      end


*     ***********************************************
*     *                                             *
*     *              meta_energy                    *
*     *                                             *
*     ***********************************************
      real*8 function meta_energy()
      implicit none

#include "mafdecls.fh"
#include "errquit.fh"
#include "metadynamics.fh"

*     **** local variables ****
      integer j,nr,yshift,y2shift
      real*8 r,dr,sigma,w,v

*     **** external functions ****
      real*8   nwpw_splint,meta_collective
      external nwpw_splint,meta_collective

      v = 0.0d0
      if (metafound) then

         if (bnmeta.gt.0) then
            dr = dbl_mb(bxmeta(1)+1) - dbl_mb(bxmeta(1))
            do j=1,bnmeta
               r     = meta_collective(1,int_mb(bindxmeta(1)+(j-1)*2))
               w     = dbl_mb(bwsigmameta(1)+(j-1)*2)
               sigma = dbl_mb(bwsigmameta(1)+(j-1)*2+1)
               yshift  = bymeta(1) +(j-1)*bnxmeta
               y2shift = by2meta(1)+(j-1)*bnxmeta
               nr = r/dr + 1
               v = v + nwpw_splint(dbl_mb(bxmeta(1)),
     >                  dbl_mb(yshift),dbl_mb(y2shift),bnmeta,nr,r)
            end do
         end if
         if (anmeta.gt.0) then
            dr = dbl_mb(axmeta(1)+1) - dbl_mb(axmeta(1))
            do j=1,anmeta
               r     = meta_collective(1,int_mb(aindxmeta(1)+(j-1)*2))
               w     = dbl_mb(awsigmameta(1)+(j-1)*2)
               sigma = dbl_mb(awsigmameta(1)+(j-1)*2+1)
               yshift  = aymeta(1) +(j-1)*anxmeta
               y2shift = ay2meta(1)+(j-1)*anxmeta
               nr = r/dr + 1
               v = v + nwpw_splint(dbl_mb(axmeta(1)),
     >                  dbl_mb(yshift),dbl_mb(y2shift),anmeta,nr,r)
            end do
         end if

      end if

      meta_energy = v 
      return
      end

*     ***********************************************
*     *                                             *
*     *              meta_found                     *
*     *                                             *
*     ***********************************************
      logical function meta_found()
      implicit none

#include "metadynamics.fh"

      meta_found = metafound
      return
      end


c     **********************************************
c     *                                            *
c     *            meta_print_potentials           *
c     *                                            *
c     **********************************************
      subroutine meta_print_potentials()
      implicit none

#include "mafdecls.fh"
#include "errquit.fh"
#include "metadynamics.fh"

*     **** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)

      integer i,j,yshift,k
      real*8 r
      character*80 filename
      character*255 full_filename

*     **** external functions ****
      character*7 c_index_name
      external    c_index_name

      call Parallel_taskid(taskid)

      if (metafound) then

*        **** print out potentials ******
         if (bnmeta.gt.0) then
            if (taskid.eq.MASTER) then
            do j=1,bnmeta

               filename = "bmeta"//c_index_name(j)//".dat"
               call util_file_name_noprefix(filename,.false.,
     >                                      .false.,
     >                                      full_filename)
               open(unit=53,file=full_filename,form='formatted')
               yshift  = bymeta(1) +(j-1)*bnxmeta
               do i=1,bnxmeta
                   write(53,'(2F15.6)') 
     >                  dbl_mb(bxmeta(1)+i-1),dbl_mb(yshift+i-1)
               end do
               close(53)
            end do
            end if
         end if
         if (anmeta.gt.0) then
           if (taskid.eq.MASTER) then
            do j=1,anmeta

               filename = "ameta"//c_index_name(j)//".dat"
               call util_file_name_noprefix(filename,.false.,
     >                                      .false.,
     >                                      full_filename)
               open(unit=53,file=full_filename,form='formatted')
               yshift  = aymeta(1) +(j-1)*anxmeta
               do i=1,bnxmeta
                   write(53,'(2F15.6)')
     >                  dbl_mb(axmeta(1)+i-1),dbl_mb(yshift+i-1)
               end do
               close(53)
            end do
            end if
         end if
      end if

      return
      end
