      subroutine neb_init(rtdb)
      implicit none
      integer rtdb

#include "nwc_const.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "neb_data.fh"
#include "errquit.fh"
#include "bead_list.fh"
c     
      logical value
      logical oprint
      character*80 title,neb_movecs
      integer itest
      integer sl
      integer nion,ng
c     

      oprint = util_print('information', print_default)
     $     .and. (ga_nodeid() .eq. 0)
c
      if (oprint) then
         write(6,*)
         write(6,*)
         call util_print_centered(6,
     $        'NWChem Minimum Energy Pathway Program',
     $        40,.true.)
         write(6,*)
         write(6,*)
      endif
c
      if (rtdb_cget(rtdb,'title',1,title)) then
         if (oprint) then
            write(6,*)
            write(6,*)
            call util_print_centered(6, title, 40, .false.)
            write(6,*)
            write(6,*)
         endif
      endif

      if (.not.rtdb_get(rtdb,'neb:nbeads',mt_int,1,nbeads)) then
         nbeads = 5
         if (.not.rtdb_put(rtdb,'neb:nbeads',mt_int,1,nbeads))
     &      call errquit("neb_init: nbeads",911,INPUT_ERR)
      end if

      if (.not.rtdb_get(rtdb,'neb:timestep',mt_dbl,1,timestep)) then
         timestep = 10.0
         if (.not.rtdb_put(rtdb,'neb:timestep',mt_dbl,1,timestep))
     &      call errquit("neb_init: timestep",911,INPUT_ERR)
      end if


c
c     check if valid restart present
c
      if (.not.rtdb_cget(rtdb,'neb:bead_list',1,neb_list)) then
         neb_list="neb"
         sl=inp_strlen(neb_list)
         if (.not.rtdb_cput(rtdb,'neb:bead_list',1,neb_list(1:sl)))
     &      call errquit("neb_init: bead_list",911,INPUT_ERR)
            neb_movecs=neb_list(1:sl)//".tmpmovecs"
         call init_bead_list(rtdb,neb_list,neb_movecs)
         call neb_initial_path(rtdb,neb_list,nbeads)
      else
         itest=size_bead_list(neb_list)
         if(nbeads.ne.itest)
     &     call errquit("conflicting nbeads values",0,0)
      end if
 
      if (oprint) then
         write(6,1) nbeads,timestep
 1       format(
     $        ' number of images in path         (nbeads) = ', i4,/,
     $        ' time step                      (timestep) = ', f10.6)
      end if

c
c     get number of ions from bead #1
c
      nion   = nion_bead_list(neb_list,1)
      ng     = 3*nion*nbeads


*     **** allocate space for gradients and coordinates ****
      value  = ma_alloc_get(mt_dbl,nbeads,'e1',e1(2),e1(1))
      value  = value.and.
     >         ma_alloc_get(mt_dbl,ng,'g0',g0(2),g0(1))
      value  = value.and.
     >         ma_alloc_get(mt_dbl,ng,'g1',g1(2),g1(1))
      value  = value.and.
     >         ma_alloc_get(mt_dbl,ng,'s',s(2),s(1))
      value  = value.and.
     >         ma_alloc_get(mt_dbl,ng,'t1',t1(2),t1(1))
      value  = value.and.
     >         ma_alloc_get(mt_dbl,ng,'v1',v1(2),v1(1))
      value  = value.and.
     >         ma_alloc_get(mt_dbl,ng,'c0',c0(2),c0(1))
      value  = value.and.
     >         ma_alloc_get(mt_dbl,ng,'c1',c1(2),c1(1))
      value  = value.and.
     >         ma_alloc_get(mt_dbl,ng,'dti',dti(2),dti(1))
      if (.not.value) 
     >  call errquit('neb failed - increase stack memory',1,0)

      return
      end

      subroutine neb_initial_path(rtdb,bead_list,nbeads)
      implicit none
      integer rtdb
      character*(*) bead_list
      integer nbeads

#include "mafdecls.fh"
#include "geom.fh"

*     **** local variables ****
      logical value
      integer i,geom,geomlen,movecslen,nion
      integer r1(2),r2(2),r3(2)
      real*8  t
      character*255 geom_name,movecs_name

*     **** external functions ****
      integer     inp_strlen
      character*7 bead_index_name
      external    inp_strlen
      external    bead_index_name

      call util_print_centered(6, 
     & "Creating Initial NEB Path", 
     &  36,
     & .true.)

      value = geom_create(geom,'neb_tmp')
      value = value.and.geom_rtdb_load(rtdb,geom,'neb_start')
      value = value.and.geom_ncent(geom,nion)
      value = value.and.geom_destroy(geom)
      if (.not.value) call errquit('neb_initial_path failed',0,0)

      value = value.and.MA_push_get(mt_dbl, (3*nion), 'r1',r1(2),r1(1))
      value = value.and.MA_push_get(mt_dbl, (3*nion), 'r2',r2(2),r2(1))
      value = value.and.MA_push_get(mt_dbl, (3*nion), 'r3',r3(2),r3(1))
      if (.not.value) call errquit('neb_initial_path failed',1,0)

      value = value.and.geom_create(geom,'neb_tmp')
      value = value.and.geom_rtdb_load(rtdb,geom,'neb_end')
      value = value.and.geom_cart_coords_get(geom,dbl_mb(r2(1)))
      value = value.and.geom_rtdb_load(rtdb,geom,'neb_start')
      value = value.and.geom_cart_coords_get(geom,dbl_mb(r1(1)))
      if (.not.value) call errquit('neb_initial_path failed',2,0)

      do i=1,nbeads
        t = (i-1)/dble(nbeads-1)

        call dcopy(3*nion,dbl_mb(r1(1)),1,dbl_mb(r3(1)),1)
        call dscal(3*nion,(1.0d0-t),dbl_mb(r3(1)),1)
        call daxpy(3*nion,t,dbl_mb(r2(1)),1,dbl_mb(r3(1)),1)

        geom_name   = 'neb_bead'//bead_index_name(i)//':geom'
        movecs_name = 'neb_bead'//bead_index_name(i)//'.movecs'
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
      if (.not.value) call errquit('neb_initial_path failed',3,0)


      return
      end





