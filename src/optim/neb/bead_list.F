c
c     $Id$
c
c
c     This file contains a list of beads, where a bead is defined
c as a movecs_filename, geom_name, and a Energy.  This list is used
c in the NEB code.
c
c     Author           - Eric Bylaska 
c     Creation Date    - 11/11/2002
c
c
c     ***************************************************
c     *                                                 *
c     *                create_bead_list                 *
c     *                                                 *
c     ***************************************************

c  This routine initializes a bead_list labeled 'tag' on
c  the runtime database, rtdb.
c
c     Entry - rtdb: pointeer to runtime database
c             tag: name of created bead list
c             perm_movecs_name: name of movecs files used by call to task
c
c 
      subroutine init_bead_list(rtdb,tag,perm_movecs_name)
      implicit none
      integer       rtdb
      character*(*) tag,perm_movecs_name

#include "rtdb.fh"
#include "mafdecls.fh"

*     ***** local variables ****
      logical value
      integer size,taglen,movecslen
      character*255 rtdb_name

      integer bead_rtdb
      common /bead_list/ bead_rtdb

*     **** external functions ****
      integer  inp_strlen
      external inp_strlen

      bead_rtdb = rtdb

      taglen    = inp_strlen(tag)
      movecslen = inp_strlen(perm_movecs_name)
   
      size = 0
      rtdb_name = tag(1:taglen)//':size'
      value = rtdb_put(bead_rtdb,rtdb_name,mt_int,1,size) 

      rtdb_name = tag(1:taglen)//':perm_movecs'
      value = value.and.
     >        rtdb_cput(bead_rtdb,rtdb_name,1,
     >                  perm_movecs_name(1:movecslen))
     
      if (.not.value) call errquit('init_bead_list failed',0,0)

      return
      end

      subroutine set_rtdb_bead_list(rtdb)
      implicit none
      integer rtdb

      integer bead_rtdb
      common /bead_list/ bead_rtdb

      bead_rtdb = rtdb
      return
      end
      

      subroutine reset_bead_list(tag)
      implicit none
      integer       rtdb
      character*(*) tag

#include "rtdb.fh"
#include "mafdecls.fh"

*     ***** local variables ****
      logical value
      integer size,taglen
      character*255 rtdb_name

      integer bead_rtdb
      common /bead_list/ bead_rtdb

*     **** external functions ****
      integer  inp_strlen
      external inp_strlen

      taglen    = inp_strlen(tag)
      size = 0
      rtdb_name = tag(1:taglen)//':size'
      value = rtdb_put(bead_rtdb,rtdb_name,mt_int,1,size) 

      return
      end

c     ***************************************************
c     *                                                 *
c     *                add_bead_list                    *
c     *                                                 *
c     ***************************************************

c  This routine add a bead to the bead_list labeled 'tag'.
c
c     Entry - tag: name of created bead list
c             movecs_name: movecs filename of added bead
c             geom_name: geometry name of added bead
c
c 
      subroutine add_bead_list(tag,movecs_name,geom_name)
      implicit none
      character*(*) tag,movecs_name,geom_name

#include "rtdb.fh"
#include "mafdecls.fh"

*     ***** local variables ****
      logical value
      integer size,taglen,movecslen,geomlen
      character*255 rtdb_name
      real*8 stress(9)

      integer bead_rtdb
      common /bead_list/ bead_rtdb

*     **** external functions ****
      character*7 bead_index_name
      integer     inp_strlen
      external    bead_index_name
      external    inp_strlen


      taglen    = inp_strlen(tag)
      movecslen = inp_strlen(movecs_name)
      geomlen   = inp_strlen(geom_name)
   
      
      rtdb_name = tag(1:taglen)//':size'
      value = rtdb_get(bead_rtdb,rtdb_name,mt_int,1,size) 
      size = size+1
      value = value.and.
     >        rtdb_put(bead_rtdb,rtdb_name,mt_int,1,size)


      rtdb_name = tag(1:taglen)//bead_index_name(size)//':movecs_name'
      value = value.and.
     >        rtdb_cput(bead_rtdb,rtdb_name,1,movecs_name(1:movecslen))

      rtdb_name = tag(1:taglen)//bead_index_name(size)//':geom_name'
      value = value.and.
     >        rtdb_cput(bead_rtdb,rtdb_name,1,geom_name(1:geomlen))

      rtdb_name = tag(1:taglen)//bead_index_name(size)//':energy'
      value = value.and.
     >        rtdb_put(bead_rtdb,rtdb_name,mt_dbl,1,0.0d0)

      call dcopy(9,0.0d0,0,stress,1)
      rtdb_name = tag(1:taglen)//bead_index_name(size)//':stress'
      value = value.and.
     >        rtdb_put(bead_rtdb,rtdb_name,mt_dbl,9,stress)

      if (.not.value) call errquit('add_bead_list failed',0,0)

      return
      end



c     ***************************************************
c     *                                                 *
c     *                delete_bead_list                 *
c     *                                                 *
c     ***************************************************

c  This routine deletes bead,i, from the bead_list labeled 'tag'.
c
c     Entry - tag: name of created bead list
c             movecs_name: movecs filename of added bead
c             geom_name: geometry name of added bead
c
c 
      subroutine delete_bead_list(tag,i)
      implicit none
      character*(*) tag
      integer i

#include "rtdb.fh"
#include "mafdecls.fh"

*     ***** local variables ****
      logical value
      integer size,taglen,tmplen,j
      real*8  energy,stress(9)
      character*255 rtdb_name,tmp_name

      integer bead_rtdb
      common /bead_list/ bead_rtdb

*     **** external functions ****
      character*7 bead_index_name
      integer     inp_strlen
      external    bead_index_name
      external    inp_strlen


      taglen    = inp_strlen(tag)
      
      rtdb_name = tag(1:taglen)//':size'
      value = rtdb_get(bead_rtdb,rtdb_name,mt_int,1,size) 
      if (i.gt.size) return

      value = value.and.rtdb_get(bead_rtdb,rtdb_name,mt_int,1,size-1) 

      do j=i,size-1

        rtdb_name = tag(1:taglen)//bead_index_name(j+1)//':movecs_name'
c        if (.not.rtdb_cget(bead_rtdb,rtdb_name,1,tmp_name))
c     >     tmp_name = 'bead'//bead_index_name(j+1)//'.movecs'
        if (.not.rtdb_cget(bead_rtdb,rtdb_name,1,tmp_name))
     >     tmp_name = tag(1:taglen)//bead_index_name(j+1)//'.movecs'
        tmplen = inp_strlen(tmp_name)
        rtdb_name = tag(1:taglen)//bead_index_name(j)//':movecs_name'
        value = value.and.
     >          rtdb_cput(bead_rtdb,rtdb_name,1,tmp_name(1:tmplen))

        rtdb_name = tag(1:taglen)//bead_index_name(j+1)//':geom_name'
c        if (.not.rtdb_cget(bead_rtdb,rtdb_name,1,tmp_name))
c     >     tmp_name = 'bead'//bead_index_name(j+1)//':geom'
        if (.not.rtdb_cget(bead_rtdb,rtdb_name,1,tmp_name))
     >     tmp_name = tag(1:taglen)//bead_index_name(j+1)//':geom'
        tmplen = inp_strlen(tmp_name)
        rtdb_name = tag(1:taglen)//bead_index_name(j)//':geom_name'
        value = value.and.
     >          rtdb_cput(bead_rtdb,rtdb_name,1,tmp_name(1:tmplen))

        rtdb_name = tag(1:taglen)//bead_index_name(j+1)//':energy'
        value = value.and.
     >          rtdb_get(bead_rtdb,rtdb_name,mt_dbl,1,energy)
        rtdb_name = tag(1:taglen)//bead_index_name(j)//':energy'
        value = value.and.
     >          rtdb_put(bead_rtdb,rtdb_name,mt_dbl,1,energy)

        rtdb_name = tag(1:taglen)//bead_index_name(j+1)//':stress'
        value = value.and.
     >          rtdb_get(bead_rtdb,rtdb_name,mt_dbl,9,stress)
        rtdb_name = tag(1:taglen)//bead_index_name(j)//':stress'
        value = value.and.
     >          rtdb_put(bead_rtdb,rtdb_name,mt_dbl,9,stress)

      end do

      if (.not.value) call errquit('delete_bead_list failed',0,0)

      return
      end


c     ***************************************************
c     *                                                 *
c     *                  size_bead_list                 *
c     *                                                 *
c     ***************************************************

c  This routine returns the number of bead in the bead_list 
c labeled "tag".
c
c     Entry - tag: name of bead list
c     Exit  - returns number of beads
c 
      integer function size_bead_list(tag)
      implicit none
      character*(*) tag

#include "rtdb.fh"
#include "mafdecls.fh"

*     ***** local variables ****
      logical value
      integer size,taglen
      character*255 rtdb_name

      integer bead_rtdb
      common /bead_list/ bead_rtdb

*     **** external functions ****
      integer     inp_strlen
      external    inp_strlen

      taglen    = inp_strlen(tag)
      
      rtdb_name = tag(1:taglen)//':size'
      value = rtdb_get(bead_rtdb,rtdb_name,mt_int,1,size) 
      if (.not.value) call errquit('size_bead_list failed',0,0)

      size_bead_list = size
      return
      end


c     ***************************************************
c     *                                                 *
c     *                  nion_bead_list                 *
c     *                                                 *
c     ***************************************************

c  This routine returns the number of ions, nions, in bead,i, 
c from the bead_list labeled 'tag'.
c
c     Entry - tag: name of bead list
c             i:  bead number
c     Exit  - returns the number ions in bead, i
c 
      integer function nion_bead_list(tag,i)
      implicit none
      character*(*) tag
      integer       i

#include "rtdb.fh"
#include "mafdecls.fh"
#include "geom.fh"

*     ***** local variables ****
      logical value,ostress
      integer size,taglen,nion,geom,geomlen
      character*255 rtdb_name,geom_name

      integer bead_rtdb
      common /bead_list/ bead_rtdb

*     **** external functions ****
      logical     bead_includestress
      integer     inp_strlen,size_bead_list
      character*7 bead_index_name
      external    bead_includestress
      external    inp_strlen,size_bead_list
      external    bead_index_name

      taglen    = inp_strlen(tag)
      
      size = size_bead_list(tag)
      if (i.gt.size) then
        nion_bead_list = 0 
        return
      end if

      rtdb_name = tag(1:taglen)//bead_index_name(i)//':geom_name'
c      if (.not.rtdb_cget(bead_rtdb,rtdb_name,1,geom_name))
c     >   geom_name   = 'bead'//bead_index_name(i)//':geom'
      if (.not.rtdb_cget(bead_rtdb,rtdb_name,1,geom_name))
     >   geom_name   = tag(1:taglen)//bead_index_name(i)//':geom'
      geomlen = inp_strlen(geom_name)            

      value =           geom_create(geom,'bead_tmp') 
      value = value.and.geom_rtdb_load(bead_rtdb,geom,
     >                                 geom_name(1:geomlen))
      value = value.and.geom_ncent(geom,nion)
      value = value.and.geom_destroy(geom)
      if (.not.value) call errquit('nion_bead_list failed',0,0)

      if (bead_includestress()) nion = nion + 3
    
      nion_bead_list = nion
      return
      end


c     ***************************************************
c     *                                                 *
c     *                  run_bead_list                  *
c     *                                                 *
c     ***************************************************

c  This routine runs bead,i, from the bead_list labeled "tag".
c
c     Entry - tag: name of bead list
c             i: bead number
c
c 
      subroutine run_bead_list(tag,i,task_gradient_gen)
      implicit none
      character*(*) tag
      integer i

#include "rtdb.fh"
#include "mafdecls.fh"
#include "geom.fh"
#include "global.fh"
#include "errquit.fh"

*     ***** local variables ****
      logical value,status,ignore
      integer geom,gradient(2),nion,nelem
      integer size,taglen,permlen,movecslen,geomlen
      real*8  energy,stress(9)
      character*255 rtdb_name,perm_name,movecs_name,geom_name
      character*32 theory

      integer bead_rtdb
      common /bead_list/ bead_rtdb

*     **** external functions ****
      logical     task_gradient_gen,bead_includestress
      external    task_gradient_gen,bead_includestress
      character*7 bead_index_name
      integer     inp_strlen
      external    bead_index_name
      external    inp_strlen

      value     = .true.
      taglen    = inp_strlen(tag)

      rtdb_name = tag(1:taglen)//':perm_movecs'
      if (.not.rtdb_cget(bead_rtdb,rtdb_name,1,perm_name))
     >   call util_file_prefix('movecs',perm_name)
      call util_file_name_resolve(perm_name, .false.)
      permlen = inp_strlen(perm_name)

      rtdb_name = tag(1:taglen)//':size'
      if (.not.rtdb_get(bead_rtdb,rtdb_name,mt_int,1,size))
     >   size = 0
      if (i.gt.size) return

      rtdb_name = tag(1:taglen)//bead_index_name(i)//':movecs_name'
c      if (.not.rtdb_cget(bead_rtdb,rtdb_name,1,movecs_name))
c     >   movecs_name = 'bead'//bead_index_name(i)//'.movecs'
      if (.not.rtdb_cget(bead_rtdb,rtdb_name,1,movecs_name))
     >   movecs_name = tag(1:taglen)//bead_index_name(i)//'.movecs'

      call util_file_name_resolve(movecs_name, .false.)
      movecslen = inp_strlen(movecs_name)

      rtdb_name = tag(1:taglen)//bead_index_name(i)//':geom_name'
c      if (.not.rtdb_cget(bead_rtdb,rtdb_name,1,geom_name))
c     >   geom_name = 'bead'//bead_index_name(i)//':geom'
      if (.not.rtdb_cget(bead_rtdb,rtdb_name,1,geom_name))
     >   geom_name = tag(1:taglen)//bead_index_name(i)//':geom'
      geomlen = inp_strlen(geom_name)

      value = value.and.
     >        rtdb_cput(bead_rtdb,'geometry',1,geom_name(1:geomlen))
  
      if (.not.value) call errquit('run_bead_list failed',0,0)

      
*     *** copy bead movecs to taskmovecs ***
      if (ga_nodeid().eq.0) then
        inquire(file=movecs_name,exist=status)
        if (status) then
           call util_file_copy(movecs_name(1:movecslen),
     >                         perm_name(1:permlen))
        end if
      end if

*     *** run gradient task ***
      if (.not.rtdb_get(bead_rtdb,"task:ignore",mt_log,1,ignore))
     >  ignore = .false.

      if(.not.rtdb_put(bead_rtdb,"scf:converged",mt_log,1,.false.))
     >  call errquit('scf:converged put',0,0)
      if(.not.rtdb_put(bead_rtdb,"dft:converged",mt_log,1,.false.))
     >  call errquit('dft:converged put',0,0)

      if(.not.rtdb_put(bead_rtdb,"neb:ibead",mt_int,1,i))
     >  call errquit('neb:ibead put',0,0)
      if(ga_nodeid().eq.0) write(*,*) "neb: running bead ",i
      value=value.and.(task_gradient_gen(bead_rtdb).or.ignore)
      if ((.not.value))
     >   call errquit('run_bead_list failed',1,0)


*     *** copy taskmovecs to bead movecs ****
      if (ga_nodeid().eq.0) then
         inquire(file=perm_name,exist=status)
         if (status) then
            call util_file_copy(perm_name(1:permlen),
     >                          movecs_name(1:movecslen))
         end if
      end if

*     ***** set the energy *****
      value = value.and.
     >        rtdb_get(bead_rtdb,'task:energy',mt_dbl,1,energy)
      rtdb_name = tag(1:taglen)//bead_index_name(i)//':energy'
      value = value.and.
     >        rtdb_put(bead_rtdb,rtdb_name,mt_dbl,1,energy)
      if (.not.value) call errquit('run_bead_list failed',2,0)


*     ***** set the forces ******
      value = value.and.geom_create(geom,'geometry')      
      value = value.and.geom_rtdb_load(bead_rtdb,geom,'geometry')
      value = value.and.geom_ncent(geom,nion)
      if (.not.value) call errquit('run_bead_list failed',3,0)

      value = value.and.
     >        MA_push_get(mt_dbl,(3*nion),
     >                    'gradient',gradient(2),gradient(1))
      value = value.and.
     >         rtdb_get(bead_rtdb,'task:gradient',mt_dbl,(3*nion),
     >                  dbl_mb(gradient(1)))

      value = value.and.geom_vel_set(geom,dbl_mb(gradient(1)))
      value = value.and.MA_pop_stack(gradient(2))
      value = value.and.geom_rtdb_delete(bead_rtdb,'geometry')
      value = value.and.geom_rtdb_store(bead_rtdb,geom,'geometry')
      value = value.and.geom_destroy(geom)
      if (.not.value) call errquit('run_bead_list failed',4,0)


*     ***** set the stresses ******
      if (bead_includestress()) then

         if (.not.rtdb_cget(bead_rtdb, 'task:theory', 1, theory))
     >   call errquit('run_bead_list: stress theory not specified',
     >                0,RTDB_ERR)
         if (theory.eq.'pspw') then
          if (.not.rtdb_get(bead_rtdb, 'pspw:stress', mt_dbl, 9,stress))
     >      call errquit('run_bead_list: could not get stress',0,0)
         else if (theory.eq.'band') then
          if (.not.rtdb_get(bead_rtdb, 'band:stress', mt_dbl, 9,stress))
     >       call errquit('run_bead_list: could not get stress',0,0)
         else if (theory.eq.'paw') then
          if (.not.rtdb_get(bead_rtdb, 'paw:stress', mt_dbl, 9,stress))
     >       call errquit('run_bead_list: could not get stress',0,0)
         else
           call errquit('run_bead_list: no stress in theory',0,RTDB_ERR)
         end if

         rtdb_name = tag(1:taglen)//bead_index_name(i)//':stress'
         value = value.and.
     >        rtdb_put(bead_rtdb,rtdb_name,mt_dbl,9,stress)
         if (.not.value) call errquit('run_bead_list failed',2,0)
      end if

      if(ga_nodeid().eq.0) then 
       write(*,*) "neb: finished bead ",i
       write(*,*) "neb: final energy",energy
      end if

      return
      end


c     ***************************************************
c     *                                                 *
c     *                  runall_bead_list               *
c     *                                                 *
c     ***************************************************

c  This routine runs all beads from the bead_list labeled 'tag'.
c
c     Entry - tag: name of bead list
c
c 
      subroutine runall_bead_list(tag,task_gradient_gen)
      implicit none
      character*(*) tag
      logical     task_gradient_gen
      external    task_gradient_gen
     
*     *** local variables ***
      integer i,nbeads

*     *** external functions ***
      integer  size_bead_list
      external size_bead_list

      nbeads = size_bead_list(tag)
      do i=1,nbeads
        call run_bead_list(tag,i,task_gradient_gen)
      end do

      return
      end



c     ***************************************************
c     *                                                 *
c     *                  runmid_bead_list               *
c     *                                                 *
c     ***************************************************

c  This routine runs all beads between 2 and nbeads-1 
c from the bead_list labeled 'tag'.
c
c     Entry - tag: name of bead list
c
c 
      subroutine runmid_bead_list(tag,task_gradient_gen,
     >                            freeze1,freezen)
      implicit none
      character*(*) tag
      logical     task_gradient_gen
      external    task_gradient_gen
      logical freeze1,freezen
     
*     *** local variables ***
      integer i,nbeads,istart,iend

*     *** external functions ***
      integer  size_bead_list
      external size_bead_list

      nbeads = size_bead_list(tag)
      !do i=2,(nbeads-1)
      istart = 2
      iend   = nbeads-1
      if (.not.freeze1) istart = 1 
      if (.not.freezen) iend   = nbeads
      do i=istart,iend
        call run_bead_list(tag,i,task_gradient_gen)
      end do

      return
      end



c     ***************************************************
c     *                                                 *
c     *                  energy_bead_list               *
c     *                                                 *
c     ***************************************************

c  This routine returns the current energy of bead, i, 
c from the bead_list labeled 'tag'.
c
c     Entry - tag: name of bead list
c             i: bead number
c 
      real*8 function energy_bead_list(tag,i)
      implicit none
      character*(*) tag
      integer i
  
#include "rtdb.fh"
#include "mafdecls.fh"

*     ***** local variables ****
      logical value
      integer taglen
      real*8  energy
      character*255 rtdb_name

      integer bead_rtdb
      common /bead_list/ bead_rtdb

*     **** external functions ****
      integer     inp_strlen
      character*7 bead_index_name
      external    inp_strlen
      external    bead_index_name

      taglen    = inp_strlen(tag)
      
      rtdb_name = tag(1:taglen)//bead_index_name(i)//':energy'
      value = rtdb_get(bead_rtdb,rtdb_name,mt_dbl,1,energy) 
      if (.not.value) call errquit('energy_bead_list failed',0,0)

      energy_bead_list = energy
      return
      end



c     ***************************************************
c     *                                                 *
c     *                 coords_get_bead_list            *
c     *                                                 *
c     ***************************************************

c  This routine returns the ion coordinates, c, of the "i" bead
c from the "tag" bead list.
c
c     Entry - tag: name of bead list
c             i: bead number
c     Exit  - c: ion coordinates
c 
      subroutine coords_get_bead_list(tag,i,c)
      implicit none
      character*(*) tag
      integer i
      real*8 c(3,*)

#include "rtdb.fh"
#include "mafdecls.fh"
#include "geom.fh"

*     ***** local variables ****
      logical value
      integer geom,nion
      integer size,taglen,geomlen
      real*8  energy
      character*255 rtdb_name,geom_name

      integer bead_rtdb
      common /bead_list/ bead_rtdb

*     **** external functions ****
      logical     bead_includestress,geom_grad_cart_to_frac
      character*7 bead_index_name
      integer     inp_strlen,size_bead_list
      external    bead_includestress,geom_grad_cart_to_frac
      external    bead_index_name
      external    inp_strlen,size_bead_list

      size = size_bead_list(tag) 
      if (i.gt.size) return
    
      taglen    = inp_strlen(tag)

      rtdb_name = tag(1:taglen)//bead_index_name(i)//':geom_name'
c      if (.not.rtdb_cget(bead_rtdb,rtdb_name,1,geom_name))
c     >   geom_name = 'bead'//bead_index_name(i)//':geom'
      if (.not.rtdb_cget(bead_rtdb,rtdb_name,1,geom_name))
     >   geom_name = tag(1:taglen)//bead_index_name(i)//':geom'
      geomlen   = inp_strlen(geom_name)            

      value =           geom_create(geom,'bead_tmp') 
      value = value.and.geom_rtdb_load(bead_rtdb,geom,
     >                                 geom_name(1:geomlen))
      value = value.and.geom_ncent(geom,nion)
      value = value.and.geom_cart_coords_get(geom,c)

      if (bead_includestress()) then

*        **** put gradient into fractional ****
        if (.not. geom_grad_cart_to_frac(geom,c))
     $  call errquit('coords_get_bead_list: cart_to_frac?',0,0)

*        **** get stress part of gradient ***
         if (.not. geom_amatrix_get(geom,c(1,nion+1)))
     >   call errquit('coords_get_bead_list:failed to get amatrix',0,0)
      end if

      value = value.and.geom_destroy(geom)
      if (.not.value) call errquit('coords_get_bead_list failed',0,0)
    
      return
      end

c     ***************************************************
c     *                                                 *
c     *                 coords_set_bead_list            *
c     *                                                 *
c     ***************************************************

c  This routine sets the ion coordinates, c, of "i" bead
c from the "tag" bead list.
c
c     Entry - tag: name of bead list
c             i: bead number
c             c: ion coordinates
c
c 
      subroutine coords_set_bead_list(tag,i,c)
      implicit none
      character*(*) tag
      integer i
      real*8 c(3,*)

#include "rtdb.fh"
#include "mafdecls.fh"
#include "geom.fh"

*     ***** local variables ****
      logical value
      integer geom,nion
      integer size,taglen,geomlen
      real*8  energy
      character*255 rtdb_name,geom_name

      integer bead_rtdb
      common /bead_list/ bead_rtdb

*     **** external functions ****
      logical     bead_includestress,geom_amatrix_set
      character*7 bead_index_name
      integer     inp_strlen,size_bead_list
      external    bead_includestress,geom_amatrix_set
      external    bead_index_name
      external    inp_strlen,size_bead_list

      size = size_bead_list(tag) 
      if (i.gt.size) return
    
      taglen    = inp_strlen(tag)

      rtdb_name = tag(1:taglen)//bead_index_name(i)//':geom_name'
c      if (.not.rtdb_cget(bead_rtdb,rtdb_name,1,geom_name))
c     >   geom_name = 'bead'//bead_index_name(i)//':geom'
      if (.not.rtdb_cget(bead_rtdb,rtdb_name,1,geom_name))
     >   geom_name = tag(1:taglen)//bead_index_name(i)//':geom'
      geomlen   = inp_strlen(geom_name)            

      value =           geom_create(geom,'bead_tmp') 
      value = value.and.geom_rtdb_load(bead_rtdb,geom,
     >                                 geom_name(1:geomlen))
      value = value.and.geom_ncent(geom,nion)

      if (bead_includestress()) then

        if (.not. geom_frac_to_cart(geom, c))
     $  call errquit('coords_set_bead_list: frac_to_cart?',0,0)

        if (.not. geom_amatrix_set(geom,c(1,nion+1)))
     $  call errquit('coords_set_bead_list:failed to set amatrix',0,0)


      end if

      value = value.and.geom_cart_coords_set(geom,c)
      value = value.and.geom_rtdb_delete(bead_rtdb,geom_name(1:geomlen))
      value = value.and.geom_rtdb_store(bead_rtdb,geom,
     >                                  geom_name(1:geomlen))
      value = value.and.geom_destroy(geom)
      if (.not.value) call errquit('coords_set_bead_list failed',0,0)
    
      return
      end



c     **********************************************
c     *                                            *
c     *              systype_bead_list             *
c     *                                            *
c     **********************************************

c  This routine returns if the systype, i.e. the number of dimensions
c from the "tag" bead list. Uses bead index=1.
c
c     Entry - tag: name of bead list
c     Exit  -
c
      integer function systype_bead_list(tag)
      implicit none
      character*(*) tag

#include "rtdb.fh"
#include "mafdecls.fh"
#include "geom.fh"

*     ***** local variables ****
      integer i,isystype
      real*8 amat(3,3)
      logical value
      integer geom,nion
      integer size,taglen,geomlen
      real*8  energy
      character*255 rtdb_name,geom_name

      integer bead_rtdb
      common /bead_list/ bead_rtdb

*     **** external functions ****
      character*7 bead_index_name
      integer     inp_strlen,size_bead_list
      external    bead_index_name
      external    inp_strlen,size_bead_list

      isystype = 0
      i = 1
      size = size_bead_list(tag)
      if (i.gt.size) then
         systype_bead_list = isystype
         return
      end if

      taglen    = inp_strlen(tag)

      rtdb_name = tag(1:taglen)//bead_index_name(i)//':geom_name'
c      if (.not.rtdb_cget(bead_rtdb,rtdb_name,1,geom_name))
c     >   geom_name = 'bead'//bead_index_name(i)//':geom'
      if (.not.rtdb_cget(bead_rtdb,rtdb_name,1,geom_name))
     >   geom_name = tag(1:taglen)//bead_index_name(i)//':geom'
      geomlen   = inp_strlen(geom_name)

      value =           geom_create(geom,'bead_tmp')
      value = value.and.geom_rtdb_load(bead_rtdb,geom,
     >                                 geom_name(1:geomlen))

*     **** get systype ***
      if (.not.geom_systype_get(geom,isystype)) then
         isystype = 0
      end if

      value = value.and.geom_destroy(geom)
      if (.not.value) call errquit('systype_bead_list failed',0,0)

      systype_bead_list = isystype
      return
      end




c     ***************************************************
c     *                                                 *
c     *              amatrix_get_bead_list              *
c     *                                                 *
c     ***************************************************

c  This routine returns the amatrix, amat, of the "i" bead
c from the "tag" bead list.
c
c     Entry - tag: name of bead list
c             i: bead number
c     Exit  - amat: amatrix
c
      subroutine amatrix_get_bead_list(tag,i,amat)
      implicit none
      character*(*) tag
      integer i
      real*8 amat(3,3)

#include "rtdb.fh"
#include "mafdecls.fh"
#include "geom.fh"

*     ***** local variables ****
      logical value
      integer geom,nion
      integer size,taglen,geomlen
      real*8  energy
      character*255 rtdb_name,geom_name

      integer bead_rtdb
      common /bead_list/ bead_rtdb

*     **** external functions ****
      character*7 bead_index_name
      integer     inp_strlen,size_bead_list
      external    bead_index_name
      external    inp_strlen,size_bead_list

      size = size_bead_list(tag)
      if (i.gt.size) return

      taglen    = inp_strlen(tag)

      rtdb_name = tag(1:taglen)//bead_index_name(i)//':geom_name'
c      if (.not.rtdb_cget(bead_rtdb,rtdb_name,1,geom_name))
c     >   geom_name = 'bead'//bead_index_name(i)//':geom'
      if (.not.rtdb_cget(bead_rtdb,rtdb_name,1,geom_name))
     >   geom_name = tag(1:taglen)//bead_index_name(i)//':geom'
      geomlen   = inp_strlen(geom_name)

      value =           geom_create(geom,'bead_tmp')
      value = value.and.geom_rtdb_load(bead_rtdb,geom,
     >                                 geom_name(1:geomlen))

*     **** get amatrix ***
      if (.not.geom_amatrix_get(geom,amat)) then
         call dcopy(9,0.0d0,0,amat,1)
      end if

      value = value.and.geom_destroy(geom)
      if (.not.value) call errquit('amatrix_get_bead_list failed',0,0)

      return
      end




c     ***************************************************
c     *                                                 *
c     *               gradient_get_bead_list            *
c     *                                                 *
c     ***************************************************

c  This routine returns the ion gradients, c, of "i" bead from
c the "tag" bead list.
c
c     Entry - tag: name of bead list
c             i: bead number
c     Exit  - c: ion forces
c            
      subroutine gradient_get_bead_list(tag,i,c)
      implicit none
      character*(*) tag
      integer i
      real*8 c(3,*)

#include "rtdb.fh"
#include "mafdecls.fh"
#include "geom.fh"

*     ***** local variables ****
      logical value
      integer geom,nion
      integer size,taglen,geomlen
      real*8  energy
      character*255 rtdb_name,geom_name

      integer bead_rtdb
      common /bead_list/ bead_rtdb

*     **** external functions ****
      logical     bead_includestress
      character*7 bead_index_name
      integer     inp_strlen,size_bead_list
      external    bead_includestress
      external    bead_index_name
      external    inp_strlen,size_bead_list

      size = size_bead_list(tag) 
      if (i.gt.size) return
    
      taglen    = inp_strlen(tag)

      rtdb_name = tag(1:taglen)//bead_index_name(i)//':geom_name'
c      if (.not.rtdb_cget(bead_rtdb,rtdb_name,1,geom_name))
c     >   geom_name   = 'bead'//bead_index_name(i)//':geom'
      if (.not.rtdb_cget(bead_rtdb,rtdb_name,1,geom_name))
     >   geom_name   = tag(1:taglen)//bead_index_name(i)//':geom'
      geomlen   = inp_strlen(geom_name)            

      value =           geom_create(geom,'bead_tmp') 
      value = value.and.geom_rtdb_load(bead_rtdb,geom,
     >                                 geom_name(1:geomlen))
      value = value.and.geom_ncent(geom,nion)
      value = value.and.geom_vel_get(geom,c)
      value = value.and.geom_destroy(geom)
      if (.not.value) call errquit('gradient_get_bead_list failed',0,0)

      if (bead_includestress()) then
         rtdb_name = tag(1:taglen)//bead_index_name(i)//':stress'
         value = value.and.
     >           rtdb_get(bead_rtdb,rtdb_name,mt_dbl,9,c(1,nion+1))
      end if
    
      return
      end

c     ***************************************************
c     *                                                 *
c     *               gradient_set_bead_list            *
c     *                                                 *
c     ***************************************************

c  This routine sets the ion gradients, c, of "i" bead from the
c "tag" bead list'.
c
c     Entry - tag: name of bead list
c             i: bead number
c             c: ion forces
c
c 
      subroutine gradient_set_bead_list(tag,i,c)
      implicit none
      character*(*) tag
      integer i
      real*8 c(3,*)

#include "rtdb.fh"
#include "mafdecls.fh"
#include "geom.fh"

*     ***** local variables ****
      logical value
      integer geom,nion
      integer size,taglen,geomlen
      real*8  energy
      character*255 rtdb_name,geom_name

      integer bead_rtdb
      common /bead_list/ bead_rtdb

*     **** external functions ****
      logical     bead_includestress
      character*7 bead_index_name
      integer     inp_strlen,size_bead_list
      external    bead_includestress
      external    bead_index_name
      external    inp_strlen,size_bead_list

      size = size_bead_list(tag) 
      if (i.gt.size) return
    
      taglen    = inp_strlen(tag)

      rtdb_name = tag(1:taglen)//bead_index_name(i)//':geom_name'
c      if (.not.rtdb_cget(bead_rtdb,rtdb_name,1,geom_name))
c     >   geom_name = 'bead'//bead_index_name(i)//':geom'
      if (.not.rtdb_cget(bead_rtdb,rtdb_name,1,geom_name))
     >   geom_name = tag(1:taglen)//bead_index_name(i)//':geom'
      geomlen   = inp_strlen(geom_name)            

      value =           geom_create(geom,'bead_tmp') 
      value = value.and.geom_rtdb_load(bead_rtdb,
     >                                 geom,geom_name(1:geomlen))
      value = value.and.geom_ncent(geom,nion)
      value = value.and.geom_vel_set(geom,c)
      value = value.and.geom_rtdb_delete(bead_rtdb,
     >                                   geom_name(1:geomlen))
      value = value.and.geom_rtdb_store(bead_rtdb,
     >                                  geom,geom_name(1:geomlen))
      value = value.and.geom_destroy(geom)

      if (bead_includestress()) then
         rtdb_name = tag(1:taglen)//bead_index_name(i)//':stress'
         value = value.and.
     >           rtdb_put(bead_rtdb,rtdb_name,mt_dbl,9,c(1,nion+1))

      end if
      if (.not.value) call errquit('gradient_get_bead_list failed',0,0)
    
      return
      end





c     ***************************************************
c     *                                                 *
c     *            create_xyz_file_bead_list            *
c     *                                                 *
c     ***************************************************

c  This routine returns the current energy of bead, i, 
c from the bead_list labeled 'tag'.
c
c     Entry - tag: name of bead list
c
c 
      subroutine create_xyz_file_bead_list(luout,tag,header)
      implicit none
      integer luout
      character*(*) tag
      logical header
  
#include "rtdb.fh"
#include "mafdecls.fh"
#include "geom.fh"
#include "global.fh"

*     ***** local variables ****
      logical value
      integer taglen,i,ii,geom,geomlen,nbeads,nion
      real*8  energy,q,rxyz(3)
      character*255 rtdb_name,geom_name
      character*2   symbol
      character*16  t,name

      integer bead_rtdb
      common /bead_list/ bead_rtdb

*     **** external functions ****
      integer     inp_strlen,size_bead_list
      character*7 bead_index_name
      external    inp_strlen,size_bead_list
      external    bead_index_name

    
      if ((ga_nodeid().eq.0).and.(header)) then
         write(luout,*)
         write(luout,*) 'XYZ FILE for bead_list:',tag
         write(luout,*) '------------------------------------------'
      end if
      taglen    = inp_strlen(tag)
      nbeads    = size_bead_list(tag)
      value = geom_create(geom,'bead_tmp') 
      do i=1,nbeads  
         rtdb_name = tag(1:taglen)//bead_index_name(i)//':energy'
         value = value.and.rtdb_get(bead_rtdb,rtdb_name,
     >                              mt_dbl,1,energy)

         rtdb_name = tag(1:taglen)//bead_index_name(i)//':geom_name'
c         if (.not.rtdb_cget(bead_rtdb,rtdb_name,1,geom_name))
c     >      geom_name = 'bead'//bead_index_name(i)//':geom'
         if (.not.rtdb_cget(bead_rtdb,rtdb_name,1,geom_name))
     >      geom_name = tag(1:taglen)//bead_index_name(i)//':geom'
         geomlen = inp_strlen(geom_name)      
         value = value.and.geom_rtdb_load(bead_rtdb,
     >                             geom,geom_name(1:geomlen))
         value = value.and.geom_ncent(geom,nion)
 
         if (ga_nodeid().eq.0) then
         write(luout,*) nion
         write(luout,*) 'energy=',energy
         end if
         do ii=1,nion       
            value = value.and.geom_cent_get(geom,ii,t,rxyz,q)
            value = value.and. geom_tag_to_element(t,symbol,name,q)
            rxyz(1)= rxyz(1)*0.529177d0
            rxyz(2)= rxyz(2)*0.529177d0
            rxyz(3)= rxyz(3)*0.529177d0
            if (ga_nodeid().eq.0) then
            write(luout,'(A2,6x,3F12.6)') symbol,rxyz
            end if
         end do
      end do
      value = value.and.geom_destroy(geom)
 
      if ((ga_nodeid().eq.0).and.(header)) then
         write(luout,*) 
      end if

      if (.not.value) 
     >  call errquit('create_xyz_file_bead_list failed',0,0)

      return
      end


c     ***********************************************
c     *                                             *
c     *            bead_index_name                  *
c     *                                             *
c     ***********************************************
      character*7 function bead_index_name(i)
      integer i

      integer itmp,j0,j1,j2,j3,j4,j5
      character*7 name

      itmp = i

      j5 = itmp/100000
      itmp = itmp - j5*100000
      j4 = itmp/10000
      itmp = itmp - j4*10000
      j3 = itmp/1000
      itmp = itmp - j3*1000
      j2 = itmp/100
      itmp = itmp - j2*100
      j1 = itmp/10
      itmp = itmp - j1*10
      j0 = itmp/1

      name(1:1) = '_'
      name(2:2) = CHAR(j5+ICHAR('0'))
      name(3:3) = CHAR(j4+ICHAR('0'))
      name(4:4) = CHAR(j3+ICHAR('0'))
      name(5:5) = CHAR(j2+ICHAR('0'))
      name(6:6) = CHAR(j1+ICHAR('0'))
      name(7:7) = CHAR(j0+ICHAR('0'))
      bead_index_name = name
      return
      end

c     ***************************************************
c     *                                                 *
c     *            bead_includestress                   *
c     *                                                 *
c     ***************************************************
      logical function bead_includestress()
      implicit none

#include "rtdb.fh"
#include "mafdecls.fh"
#include "geom.fh"

*     ***** local variables ****
      logical ostress

      integer bead_rtdb
      common /bead_list/ bead_rtdb

      if (.not.rtdb_get(bead_rtdb,'includestress',mt_log,1,ostress))
     >   ostress = .false.

      bead_includestress = ostress
      return
      end 

