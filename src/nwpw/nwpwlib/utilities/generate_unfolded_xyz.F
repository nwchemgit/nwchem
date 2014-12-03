c
c     $Id$
c

*     *********************************
*     *                               *
*     *     generate_unfolded_xyz     *
*     *                               *
*     *********************************
      subroutine generate_unfolded_xyz(motionfilename,xyzfilename)
      implicit none
      character*(*) motionfilename
      character*(*) xyzfilename

#include "bafdecls.fh"
#include "stdio.fh"
#include "errquit.fh"

c     *** parameters ***
      real*8 DELR
      parameter (DELR=0.01d0)
 
      real*8    AACONV
      parameter (AACONV=0.529177d0)
      integer MASTER,taskid
      parameter (MASTER=0)

*     **** local variables ****
      logical value
      integer i,j,bin,natoms,natoms1,t,nframes,iii
      integer rion1(2),rion2(2),vion2(2),tags(2)
      real*8  time1,time2,volume,volume1
      real*8  x,y,z,r,gr,v1,v2,v3
      real*8 unita(3,3)
      character*255 full_filename
      character*2 tag2
      character*4 tag4

*     *** synchronization ***
      call ga_sync()

      call Parallel_taskid(taskid)

*     ***************************************
*     ****** PAIR DISTRIBUTION routine ******
*     ***************************************
      if (taskid.eq.MASTER) then

*     *** get the number of frames from MOTION ***
      call MOTION_nframes(motionfilename,nframes,natoms,volume,unita)

      write(luout,*)
      write(luout,*) 
     =  '===== Generating unfolded xyz trajectory from MOTION file ===='
      write(luout,*) 
      write(luout,*) "nframes = ",nframes
      write(luout,*) "natoms  = ",natoms
      write(luout,*) 
      write(luout,*) "MOTION filename      : ",motionfilename
      write(luout,*) "unfolded xyz filename: ",xyzfilename


      open(unit=18,file=xyzfilename,form='formatted')
      open(unit=19,file=motionfilename,form='formatted')

      value = BA_push_get(mt_dbl,3*natoms,'rion1',rion1(2),rion1(1))
      value = value.and.
     >        BA_push_get(mt_dbl,3*natoms,'rion2',rion2(2),rion2(1))
      value = value.and.
     >        BA_push_get(mt_dbl,3*natoms,'vion2',vion2(2),vion2(1))
      value = value.and.
     >        BA_push_get(mt_byte,2*natoms,'tags',tags(2),tags(1))
      if (.not. value) 
     >    call errquit(
     >         'generate_unfolded_xyz: out of stack memory',0, MA_ERR)
      
      call dcopy(3*natoms,0.0d0,0,dbl_mb(rion1(1)),1)
    
      time1 = 0.0
      do t=1,nframes
        
         time2    = time1
         read(19,*) time1,natoms1,volume1,unita
         do i=1,natoms1
            tag2(1:1) = ''
            tag2(2:2) = ''
            read(19,*) iii,tag2,tag4,x,y,z,v1,v2,v3
            byte_mb(tags(1)+2*(i-1))   = tag2(1:1)
            byte_mb(tags(1)+2*(i-1)+1) = tag2(2:2)
            dbl_mb(rion2(1)+3*(i-1))   = x
            dbl_mb(rion2(1)+3*(i-1)+1) = y
            dbl_mb(rion2(1)+3*(i-1)+2) = z          
            dbl_mb(vion2(1)+3*(i-1))   = v1
            dbl_mb(vion2(1)+3*(i-1)+1) = v2
            dbl_mb(vion2(1)+3*(i-1)+2) = v3          
         end do
         if (t.gt.1) 
     >      call nwpw_unfold(natoms1,unita,
     >                       dbl_mb(rion1(1)),dbl_mb(rion2(1)))

         write(18,110) natoms1
         do i=1,natoms1
            tag2(1:1) = byte_mb(tags(1)+2*(i-1))
            tag2(2:2) = byte_mb(tags(1)+2*(i-1)+1)
            write(18,111) tag2,
     >                    dbl_mb(rion2(1)+3*(i-1))  *AACONV,
     >                    dbl_mb(rion2(1)+3*(i-1)+1)*AACONV,
     >                    dbl_mb(rion2(1)+3*(i-1)+2)*AACONV,
     >                    dbl_mb(vion2(1)+3*(i-1))  *AACONV,
     >                    dbl_mb(vion2(1)+3*(i-1)+1)*AACONV,
     >                    dbl_mb(vion2(1)+3*(i-1)+2)*AACONV
         end do
         call util_flush(18)

         call dcopy(3*natoms1,dbl_mb(rion2(1)),1,dbl_mb(rion1(1)),1)
      end do
      close(unit=18)
      close(unit=19)

      value =           BA_pop_stack(tags(2))
      value = value.and.BA_pop_stack(vion2(2))
      value = value.and.BA_pop_stack(rion2(2))
      value = value.and.BA_pop_stack(rion1(2))
      if (.not. value) 
     >  call errquit(
     >       'generate_unfolded_xyz: cannot release stack memory',0,
     &       MA_ERR)
     
      end if

*     *** synchronization ***
      call ga_sync()

      return
  110 format(I5/)
  111 format(A2,5x, 6e12.4)
      end

