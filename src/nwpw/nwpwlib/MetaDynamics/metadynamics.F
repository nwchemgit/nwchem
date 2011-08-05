
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
      integer i,j,n1,n2,is1,is2
      real*8 x,y,z,r,rmax,dr,pi,da,maxcoord
      character*80 rtdb_name
      character*4 celement1,celement2

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

      rtdb_name = 'metadynamics_ncoordnum'
      if (.not.rtdb_get(rtdb,rtdb_name,mt_int,1,cnmeta)) cnmeta = 0

      metafound = (bnmeta.gt.0).or.(anmeta.gt.0).or.(cnmeta.gt.0)
      if (metafound) then

         rtdb_name = 'metadynamics_update'
         if (.not.rtdb_get(rtdb,rtdb_name,mt_int,1,maxmetacount)) 
     >      maxmetacount = 1

         rtdb_name = 'metadynamics_metacount'
         if (.not.rtdb_get(rtdb,rtdb_name,mt_int,1,metacount)) 
     >      metacount = 0

         rtdb_name = 'metadynamics_print_shift'
         if (.not.rtdb_get(rtdb,rtdb_name,mt_int,1,metarayshift)) 
     >      metarayshift = 0

         rtdb_name = 'metadynamics_metaraycount'
         if (.not.rtdb_get(rtdb,rtdb_name,mt_int,1,metaraycount)) 
     >      metaraycount = 0

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




         if (cnmeta.gt.0) then
            cnxmeta = 501
            value = value.and.
     >              MA_alloc_get(mt_dbl,cnxmeta,'cxmeta',
     >                           cxmeta(2),cxmeta(1))
            value = value.and.
     >              MA_alloc_get(mt_dbl,cnmeta*cnxmeta,'cymeta',
     >                           cymeta(2),cymeta(1))
            value = value.and.
     >              MA_alloc_get(mt_dbl,cnmeta*cnxmeta,'cy2meta',
     >                           cy2meta(2),cy2meta(1))

*           **** load old spline fits ****
            call dcopy(cnmeta*cnxmeta,0.0d0,0,dbl_mb(cymeta(1)),1)
            call dcopy(cnmeta*cnxmeta,0.0d0,0,dbl_mb(cy2meta(1)),1)
            rtdb_name = 'metadynamics_cymeta'
            if (.not.rtdb_get(rtdb,rtdb_name,mt_dbl,cnmeta*cnxmeta,
     >                       dbl_mb(cymeta(1)))) then
               call dcopy(cnmeta*cnxmeta,0.0d0,0,dbl_mb(cymeta(1)),1)
            else
               if (taskid.eq.MASTER)
     >           write(*,*) 
     >           "  ... reading coordination number metadynamics data"
            end if
            rtdb_name = 'metadynamics_cy2meta'
            if (.not.rtdb_get(rtdb,rtdb_name,mt_dbl,cnmeta*cnxmeta,
     >                       dbl_mb(cy2meta(1)))) then
               call dcopy(cnmeta*cnxmeta,0.0d0,0,dbl_mb(cy2meta(1)),1)
            else
               if (taskid.eq.MASTER)
     >           write(*,*) 
     >        "  ... reading more coordination number metadynamics data"
            end if

            value = value.and.
     >              MA_alloc_get(mt_dbl,9*cnmeta,'cwsigmameta',
     >                           cwsigmameta(2),cwsigmameta(1))
            value = value.and.
     >              MA_alloc_get(mt_int,cnmeta,'cindxsize1',
     >                           cindxsize1(2),cindxsize1(1))
            value = value.and.
     >              MA_alloc_get(mt_int,cnmeta,'cindxsize2',
     >                           cindxsize2(2),cindxsize2(1))

            maxcoord = 0.0d0
            is1 = 0
            is2 = 0
            do j=1,cnmeta
              rtdb_name =
     >        'metadynamics_coordnum_indxsize1:'//c_index_name(j)
              value=value.and.rtdb_get(rtdb,rtdb_name,mt_int,1,
     >                                 int_mb(cindxsize1(1)+j-1))
              rtdb_name =
     >        'metadynamics_coordnum_indxsize2:'//c_index_name(j)
              value=value.and.rtdb_get(rtdb,rtdb_name,mt_int,1,
     >                                 int_mb(cindxsize2(1)+j-1))
              is1 = is1 + int_mb(cindxsize1(1)+j-1)
              is2 = is2 + int_mb(cindxsize2(1)+j-1)

              rtdb_name='metadynamics_coordnum_ksigma:'//c_index_name(j)
              if (.not.rtdb_get(rtdb,rtdb_name,mt_dbl,9,
     >                      dbl_mb(cwsigmameta(1)+(j-1)*9))) then
                 dbl_mb(cwsigmameta(1)+(j-1)*9)   = 6.0d0
                 dbl_mb(cwsigmameta(1)+(j-1)*9+1) = 18.0d0
                 dbl_mb(cwsigmameta(1)+(j-1)*9+2) = 3.0d0
                 dbl_mb(cwsigmameta(1)+(j-1)*9+3) = 0.00005d0
                 dbl_mb(cwsigmameta(1)+(j-1)*9+4) = 0.10d0
                 dbl_mb(cwsigmameta(1)+(j-1)*9+5) = 25.0d0
                 dbl_mb(cwsigmameta(1)+(j-1)*9+6) = -1.0d0
                 dbl_mb(cwsigmameta(1)+(j-1)*9+7) = -1.0d0
                 dbl_mb(cwsigmameta(1)+(j-1)*9+8) = -1.0d0
              end if
              if (dbl_mb(cwsigmameta(1)+(j-1)*9).le.0.0d0)
     >            dbl_mb(cwsigmameta(1)+(j-1)*9)   = 6.0d0
              if (dbl_mb(cwsigmameta(1)+(j-1)*9+1).le.0.0d0)
     >            dbl_mb(cwsigmameta(1)+(j-1)*9+1)   = 18.0d0
              if (dbl_mb(cwsigmameta(1)+(j-1)*9+2).le.0.0d0)
     >            dbl_mb(cwsigmameta(1)+(j-1)*9+2)   = 3.0d0 
              if (dbl_mb(cwsigmameta(1)+(j-1)*9+3).le.0.0d0)
     >            dbl_mb(cwsigmameta(1)+(j-1)*9+3)   = 0.00005d0
              if (dbl_mb(cwsigmameta(1)+(j-1)*9+4).le.0.0d0)
     >            dbl_mb(cwsigmameta(1)+(j-1)*9+4)   = 0.10d0
              if (dbl_mb(cwsigmameta(1)+(j-1)*9+5).le.0.0d0)
     >            dbl_mb(cwsigmameta(1)+(j-1)*9+5)   = 25.0d0

              if (dbl_mb(cwsigmameta(1)+(j-1)*9+5).gt.maxcoord) 
     >           maxcoord=dbl_mb(cwsigmameta(1)+(j-1)*9+5)
            end do

            value = value.and.
     >              MA_alloc_get(mt_int,is1*cnmeta,'cindxmeta1',
     >                           cindxmeta1(2),cindxmeta1(1))
            value = value.and.
     >              MA_alloc_get(mt_int,is2*cnmeta,'cindxmeta2',
     >                           cindxmeta2(2),cindxmeta2(1))
            is1 = 0
            is2 = 0
            do j=1,cnmeta
              n1 = int_mb(cindxsize1(1)+j-1)
              n2 = int_mb(cindxsize2(1)+j-1)
              rtdb_name ='metadynamics_coordnum_indx1:'//c_index_name(j)
              value=value.and.rtdb_get(rtdb,rtdb_name,mt_int,n1,
     >                                 int_mb(cindxmeta1(1)+is1))
              rtdb_name ='metadynamics_coordnum_indx2:'//c_index_name(j)
              value=value.and.rtdb_get(rtdb,rtdb_name,mt_int,n2,
     >                                 int_mb(cindxmeta2(1)+is2))
              is1 = is1 + n1
              is2 = is2 + n2

            end do

            da = maxcoord/dble(cnxmeta-1)
            do i=1,cnxmeta
               r = (i-1)*da
               dbl_mb(cxmeta(1)+i-1) = r
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
         is1=0
         is2=0
         do i=1,cnmeta
            n1 = int_mb(cindxsize1(1)+i-1)
            n2 = int_mb(cindxsize2(1)+i-1)
            write(*,1001) (int_mb(cindxmeta1(1)+is1+j-1),j=1,n1)
            write(*,1002) (int_mb(cindxmeta2(1)+is2+j-1),j=1,n2)

            write(*,'(A,/6x,A,F11.6,4x,A,F11.6,
     >                  /6x,A,F11.6,4x,A,F11.6,
     >                  /6x,A,F11.6,4x,A,F11.6
     >                  /6x,A,F11.6,F11.6)')
     >       "  - Coordination Number Parameters: ",
     >                  'n= ',
     >                  dbl_mb(cwsigmameta(1)+9*(i-1)+0),
     >                  'm=       ',
     >                  dbl_mb(cwsigmameta(1)+9*(i-1)+1),
     >                  'w= ',
     >                  dbl_mb(cwsigmameta(1)+9*(i-1)+3),
     >                  'sigma=   ',
     >                  dbl_mb(cwsigmameta(1)+9*(i-1)+4),
     >                  'r0=',
     >                  dbl_mb(cwsigmameta(1)+9*(i-1)+2),
     >                  'maxcoord=',
     >                  dbl_mb(cwsigmameta(1)+9*(i-1)+5),
     >                  'range=   ',
     >                  dbl_mb(cwsigmameta(1)+9*(i-1)+6),
     >                  dbl_mb(cwsigmameta(1)+9*(i-1)+7)
              if (dbl_mb(cwsigmameta(1)+9*(i-1)+8).lt.0) then
                 write(*,'(6x,A)') '- LJ function form'
              else
                 write(*,'(6x,A)') '- Sprik function form'
              end if
            is1 = is1 + n1
            is2 = is2 + n2
         end do
 1001 FORMAT(2x,"- Coorination Number (Index1) :",10I5)
 1002 FORMAT(2x,"- Coorination Number (Index2) :",10I5)

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
         metaraycount = metaraycount + metarayshift
         call meta_print_potentials(metaraycount)

         value=.true.
         rtdb_name='metadynamics_metacount'
         value= value.and.rtdb_put(rtdb,rtdb_name,mt_int,1,metacount)
         rtdb_name='metadynamics_metaraycount'
         value= value.and.rtdb_put(rtdb,rtdb_name,mt_int,1,metaraycount)

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
         if (cnmeta.gt.0) then

*           **** write out current spline fits ****
            rtdb_name = 'metadynamics_cymeta'
            value = value.and.rtdb_put(rtdb,rtdb_name,mt_dbl,
     >                        cnmeta*cnxmeta,dbl_mb(cymeta(1)))
            rtdb_name = 'metadynamics_cy2meta'
            value = value.and.rtdb_put(rtdb,rtdb_name,mt_dbl,
     >                        cnmeta*cnxmeta,dbl_mb(cy2meta(1)))

            value = value.and.MA_free_heap(cxmeta(2))
            value = value.and.MA_free_heap(cymeta(2))
            value = value.and.MA_free_heap(cy2meta(2))
            value = value.and.MA_free_heap(cwsigmameta(2))
            value = value.and.MA_free_heap(cindxmeta1(2))
            value = value.and.MA_free_heap(cindxmeta2(2))
            value = value.and.MA_free_heap(cindxsize1(2))
            value = value.and.MA_free_heap(cindxsize2(2))
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
      integer i,j
      real*8 x1,y1,z1,x2,y2,z2,x3,y3,z3,dx,dy,dz,r,r1,r3
      real*8 dx1,dy1,dz1,dx3,dy3,dz3,theta,n,m,r0,f

*     **** external functions ****
      real*8   ion_rion
      external ion_rion

      f = 0.0d0
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
         f  = r

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
         f = theta

c     **** bond dihedral ****

      end if

c     **** coordination number ****

      meta_collective = f
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
      integer i,j
      real*8 x1,y1,z1,r1,vx1,vx2,vy1,vy2,vz1,vz2
      real*8 x2,y2,z2,aa,a11,a12,a22,r,df,n,m,r0,rn,rm
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


c     **********************************************
c     *                                            *
c     *            meta_collective2                *
c     *                                            *
c     **********************************************
      real*8 function meta_collective2(ntype,n1,indx1,n2,indx2,param)
      implicit none
      integer ntype,n1,indx1(*),n2,indx2(*)
      real*8  param(*)

*     **** local variables ****
      integer i,j
      real*8 x1,y1,z1,x2,y2,z2,x3,y3,z3,dx,dy,dz,r,r1,r3
      real*8 dx1,dy1,dz1,dx3,dy3,dz3,theta,n,m,r0,f

*     **** external functions ****
      real*8   ion_rion
      external ion_rion

      f = 0.0d0
c     **** coordination number ****
      if (ntype.eq.4) then
         if (param(9).lt.0.0d0) then
            f = 0.0d0
            n  = param(1)
            m  = param(2)
            r0 = param(3)
            do i=1,n1
               do j=1,n2
                  x1 = ion_rion(1,indx1(i))
                  y1 = ion_rion(2,indx1(i))
                  z1 = ion_rion(3,indx1(i))
                  x2 = ion_rion(1,indx2(j))
                  y2 = ion_rion(2,indx2(j))
                  z2 = ion_rion(3,indx2(j))
                  dx1 = x1-x2
                  dy1 = y1-y2
                  dz1 = z1-z2
                  call lattice_min_difference(dx1,dy1,dz1)
                  r1  = dsqrt(dx1**2 + dy1**2 + dz1**2)
                  f = f + (1.0d0-(r1/r0)**n)/(1.0d0-(r1/r0)**m)
               end do
            end do
         else
            f = 0.0d0
            n  = param(1)
            m  = param(2)
            r0 = param(3)
            do i=1,n1
               do j=1,n2
                  x1 = ion_rion(1,indx1(i))
                  y1 = ion_rion(2,indx1(i))
                  z1 = ion_rion(3,indx1(i))
                  x2 = ion_rion(1,indx2(j))
                  y2 = ion_rion(2,indx2(j))
                  z2 = ion_rion(3,indx2(j))
                  dx1 = x1-x2
                  dy1 = y1-y2
                  dz1 = z1-z2
                  call lattice_min_difference(dx1,dy1,dz1)
                  r1  = dsqrt(dx1**2 + dy1**2 + dz1**2)
c                 f = f + (1.0d0-(r1/r0)**n)/(1.0d0-(r1/r0)**m)
                  f = f + 1.d0/(1.d0+dexp(n*(r1-r0)))
               end do
            end do
         end if

      end if

c     **** coordination number ****

      meta_collective2 = f
      return
      end

c     **********************************************
c     *                                            *
c     *           meta_collective2_force           *
c     *                                            *
c     **********************************************

      subroutine meta_collective2_force(ntype,
     >                                  n1,indx1,
     >                                  n2,indx2,param,
     >                                  dv,fmeta)
      implicit none
      integer ntype,n1,indx1(*),n2,indx2(*)
      real*8 param(*),dv,fmeta(3,*)

*     **** local variables ****
      integer i,j
      real*8 x1,y1,z1,r1,vx1,vx2,vy1,vy2,vz1,vz2
      real*8 x2,y2,z2,aa,a11,a12,a22,r,df,n,m,r0,rn,rm
      real*8 x3,y3,z3,r3,ctheta,stheta,denom
      real*8 dx,dy,dz,dx1,dy1,dz1,dx3,dy3,dz3

*     **** external functions ****
      real*8   ion_rion
      external ion_rion

c     **** coordination number ***
      if (ntype.eq.4) then
         if (param(9).lt.0.0d0) then
            n  = param(1)
            m  = param(2)
            r0 = param(3)
            do i=1,n1
               do j=1,n2
                  x1 = ion_rion(1,indx1(i))
                  y1 = ion_rion(2,indx1(i))
                  z1 = ion_rion(3,indx1(i))
                  x2 = ion_rion(1,indx2(j))
                  y2 = ion_rion(2,indx2(j))
                  z2 = ion_rion(3,indx2(j))
                  dx = x1-x2
                  dy = y1-y2
                  dz = z1-z2
                  call lattice_min_difference(dx,dy,dz)
                  r  = dsqrt(dx**2 + dy**2 + dz**2)
                  rn  = (1.0d0-(r/r0)**n)
                  rm  = (1.0d0-(r/r0)**m)
                  df = (-n*rm/r0*(r/r0)**(n-1) + m*rn/r0*(r/r0)**(m-1))
     >                 / (rm**2)
                  fmeta(1,indx1(i)) = fmeta(1,indx1(i)) - (dx/r)*df*dv
                  fmeta(2,indx1(i)) = fmeta(2,indx1(i)) - (dy/r)*df*dv
                  fmeta(3,indx1(i)) = fmeta(3,indx1(i)) - (dz/r)*df*dv
                  fmeta(1,indx2(j)) = fmeta(1,indx2(j)) + (dx/r)*df*dv
                  fmeta(2,indx2(j)) = fmeta(2,indx2(j)) + (dy/r)*df*dv
                  fmeta(3,indx2(j)) = fmeta(3,indx2(j)) + (dz/r)*df*dv
               end do
            end do
         else
            n  = param(1)
            m  = param(2)
            r0 = param(3)
            do i=1,n1
               do j=1,n2
                  x1 = ion_rion(1,indx1(i))
                  y1 = ion_rion(2,indx1(i))
                  z1 = ion_rion(3,indx1(i))
                  x2 = ion_rion(1,indx2(j))
                  y2 = ion_rion(2,indx2(j))
                  z2 = ion_rion(3,indx2(j))
                  dx = x1-x2
                  dy = y1-y2
                  dz = z1-z2
                  call lattice_min_difference(dx,dy,dz)
                  r  = dsqrt(dx**2 + dy**2 + dz**2)
c                 rn  = (1.0d0-(r/r0)**n)
c                 rm  = (1.0d0-(r/r0)**m)
c                 df = (-n*rm/r0*(r/r0)**(n-1) + m*rn/r0*(r/r0)**(m-1))
c    >                 / (rm**2)
                  rn  = 1.d0+dexp(n*(r-r0))
                  df = -n*(rn-1.d0)/(rn*rn)
                  fmeta(1,indx1(i)) = fmeta(1,indx1(i)) - (dx/r)*df*dv
                  fmeta(2,indx1(i)) = fmeta(2,indx1(i)) - (dy/r)*df*dv
                  fmeta(3,indx1(i)) = fmeta(3,indx1(i)) - (dz/r)*df*dv
                  fmeta(1,indx2(j)) = fmeta(1,indx2(j)) + (dx/r)*df*dv
                  fmeta(2,indx2(j)) = fmeta(2,indx2(j)) + (dy/r)*df*dv
                  fmeta(3,indx2(j)) = fmeta(3,indx2(j)) + (dz/r)*df*dv
               end do
            end do
         endif

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
      integer i,j,yshift,y2shift,is1,is2,n1,n2
      real*8  r0,r,dr,w,sigma,yp1,ypn,a,b

*     **** external functions ****
      real*8   meta_collective,meta_collective2
      external meta_collective,meta_collective2

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

         if (cnmeta.gt.0) then
            dr = dbl_mb(cxmeta(1)+1) - dbl_mb(cxmeta(1))
            is1 = 0
            is2 = 0
            do j=1,cnmeta
               n1 = int_mb(cindxsize1(1)+j-1)
               n2 = int_mb(cindxsize2(1)+j-1)
               r0    = meta_collective2(4,
     >                 n1,int_mb(cindxmeta1(1)+is1),
     >                 n2,int_mb(cindxmeta2(1)+is2),
     >                 dbl_mb(cwsigmameta(1)+(j-1)*9))
               w     = dbl_mb(cwsigmameta(1)+(j-1)*9+3)
               sigma = dbl_mb(cwsigmameta(1)+(j-1)*9+4)
               a     = dbl_mb(cwsigmameta(1)+(j-1)*9+6)
               b     = dbl_mb(cwsigmameta(1)+(j-1)*9+7)
               yshift  = cymeta(1) +(j-1)*cnxmeta
               y2shift = cy2meta(1)+(j-1)*cnxmeta
               do i=1,cnxmeta
                  r = dbl_mb(cxmeta(1)+i-1)-r0
                  dbl_mb(yshift+i-1) = dbl_mb(yshift+i-1)
     >                               + w*dexp(-(r/sigma)**2)

                  if ((r.lt.a).and.(a.gt.0)) then
                     dbl_mb(yshift+i-1) = dbl_mb(yshift+i-1) 
     >                                  + 2.0d0*w*dexp(-(r/sigma)**2)
                  end if
                  if ((r.gt.b).and.(b.gt.0)) then
                     dbl_mb(yshift+i-1) = dbl_mb(yshift+i-1) 
     >                                  + 2.0d0*w*dexp(-(r/sigma)**2)
                  end if
               end do

               !**** five point formula ***
               yp1 = ( -50.0d0*dbl_mb(yshift)
     >                + 96.0d0*dbl_mb(yshift+1)
     >                - 72.0d0*dbl_mb(yshift+2)
     >                + 32.0d0*dbl_mb(yshift+3)
     >                -  6.0d0*dbl_mb(yshift+4))/(24.0d0*dr)
               ypn = ( +50.0d0*dbl_mb(yshift+cnxmeta-5)
     >                - 96.0d0*dbl_mb(yshift+cnxmeta-4)
     >                + 72.0d0*dbl_mb(yshift+cnxmeta-3)
     >                - 32.0d0*dbl_mb(yshift+cnxmeta-2)
     >                +  6.0d0*dbl_mb(yshift+cnxmeta-1))/(24.0d0*dr)
               call nwpw_spline(dbl_mb(cxmeta(1)),
     >                 dbl_mb(yshift),cnxmeta,yp1,ypn,
     >                 dbl_mb(y2shift),dbl_mb(tmeta(1)))
               is1 = is1 + n1
               is2 = is2 + n2
            end do
         end if

         end if

         metaprintcount = metaprintcount + 1
         if (metaprintcount.ge.maxmetaprintcount) then
            metaprintcount = 0
            metaraycount   = metaraycount + metarayshift
            call meta_print_potentials(metaraycount)
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
      integer j,nr,yshift,y2shift,is1,is2,n1,n2
      real*8 r,dr,sigma,w,dv

*     **** external functions ****
      real*8   nwpw_dsplint,meta_collective,meta_collective2
      external nwpw_dsplint,meta_collective,meta_collective2
      
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
               r     = meta_collective(2,int_mb(aindxmeta(1)+(j-1)*3))
               w     = dbl_mb(awsigmameta(1)+(j-1)*2)
               sigma = dbl_mb(awsigmameta(1)+(j-1)*2+1)
               yshift  = aymeta(1) +(j-1)*anxmeta
               y2shift = ay2meta(1)+(j-1)*anxmeta
               nr = r/dr + 1
               dv = nwpw_dsplint(dbl_mb(axmeta(1)),
     >                  dbl_mb(yshift),dbl_mb(y2shift),anmeta,nr,r)
               call meta_collective_force(2,
     >                  int_mb(aindxmeta(1)+(j-1)*3),dv,fmeta)
            end do
         end if

         if (cnmeta.gt.0) then
            dr = dbl_mb(cxmeta(1)+1) - dbl_mb(cxmeta(1))
            is1 = 0
            is2 = 0
            do j=1,cnmeta
               n1 = int_mb(cindxsize1(1)+j-1)
               n2 = int_mb(cindxsize2(1)+j-1)
               r     = meta_collective2(4,
     >                 n1,int_mb(cindxmeta1(1)+is1),
     >                 n2,int_mb(cindxmeta2(1)+is2),
     >                 dbl_mb(cwsigmameta(1)+(j-1)*9))
               w     = dbl_mb(cwsigmameta(1)+(j-1)*9+3)
               sigma = dbl_mb(cwsigmameta(1)+(j-1)*9+4)
               yshift  = cymeta(1) +(j-1)*cnxmeta
               y2shift = cy2meta(1)+(j-1)*cnxmeta
               nr = r/dr + 1
               dv = nwpw_dsplint(dbl_mb(cxmeta(1)),
     >                  dbl_mb(yshift),dbl_mb(y2shift),cnmeta,nr,r)
               call meta_collective2_force(4,
     >                   n1,int_mb(cindxmeta1(1)+is1),
     >                   n2,int_mb(cindxmeta2(1)+is2),
     >                   dbl_mb(cwsigmameta(1)+(j-1)*9),
     >                   dv,fmeta)
               is1 = is1 + n1
               is2 = is2 + n2
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
      integer j,nr,yshift,y2shift,n1,n2,is1,is2
      real*8 r,dr,sigma,w,v

*     **** external functions ****
      real*8   nwpw_splint,meta_collective,meta_collective2
      external nwpw_splint,meta_collective,meta_collective2

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
               r     = meta_collective(2,int_mb(aindxmeta(1)+(j-1)*3))
               w     = dbl_mb(awsigmameta(1)+(j-1)*2)
               sigma = dbl_mb(awsigmameta(1)+(j-1)*2+1)
               yshift  = aymeta(1) +(j-1)*anxmeta
               y2shift = ay2meta(1)+(j-1)*anxmeta
               nr = r/dr + 1
               v = v + nwpw_splint(dbl_mb(axmeta(1)),
     >                  dbl_mb(yshift),dbl_mb(y2shift),anmeta,nr,r)
            end do
         end if
         if (cnmeta.gt.0) then
            dr = dbl_mb(cxmeta(1)+1) - dbl_mb(cxmeta(1))
            is1=0
            is2=0
            do j=1,cnmeta
               n1 = int_mb(cindxsize1(1)+j-1)
               n2 = int_mb(cindxsize2(1)+j-1)
               r     = meta_collective2(4,
     >                 n1,int_mb(cindxmeta1(1)+is1),
     >                 n2,int_mb(cindxmeta2(1)+is2),
     >                 dbl_mb(cwsigmameta(1)+(j-1)*9))
               w     = dbl_mb(cwsigmameta(1)+(j-1)*9+3)
               sigma = dbl_mb(cwsigmameta(1)+(j-1)*9+4)
               yshift  = cymeta(1) +(j-1)*cnxmeta
               y2shift = cy2meta(1)+(j-1)*cnxmeta
               nr = r/dr + 1
               v = v + nwpw_splint(dbl_mb(cxmeta(1)),
     >                  dbl_mb(yshift),dbl_mb(y2shift),cnmeta,nr,r)
               is1 = is1 + n1
               is2 = is2 + n2
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
      subroutine meta_print_potentials(icount)
      implicit none
      integer icount

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
               if (icount.gt.0) then
                  filename = 
     >            "bmeta"//c_index_name(j)//c_index_name(icount)//".dat"
               else
                  filename = "bmeta"//c_index_name(j)//".dat"
               end if

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
               if (icount.gt.0) then
                  filename = 
     >            "ameta"//c_index_name(j)//c_index_name(icount)//".dat"
               else
                  filename = "ameta"//c_index_name(j)//".dat"
               end if

               call util_file_name_noprefix(filename,.false.,
     >                                      .false.,
     >                                      full_filename)
               open(unit=53,file=full_filename,form='formatted')
               yshift  = aymeta(1) +(j-1)*anxmeta
               do i=1,anxmeta
                   write(53,'(2F15.6)')
     >                  dbl_mb(axmeta(1)+i-1),dbl_mb(yshift+i-1)
               end do
               close(53)
            end do
            end if
         end if

         if (cnmeta.gt.0) then
           if (taskid.eq.MASTER) then
            do j=1,cnmeta
               if (icount.gt.0) then
                  filename = 
     >            "cmeta"//c_index_name(j)//c_index_name(icount)//".dat"
               else
                  filename = "cmeta"//c_index_name(j)//".dat"
               end if
               call util_file_name_noprefix(filename,.false.,
     >                                      .false.,
     >                                      full_filename)
               open(unit=53,file=full_filename,form='formatted')
               yshift  = cymeta(1) +(j-1)*cnxmeta
               do i=1,cnxmeta
                   write(53,'(2F15.6)')
     >                  dbl_mb(cxmeta(1)+i-1),dbl_mb(yshift+i-1)
               end do
               close(53)
            end do
            end if
         end if


      end if

      return
      end
