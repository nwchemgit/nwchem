c
c     $Id$
c
*     *********************************
*     *                               *
*     *          pair_dist            *
*     *                               *
*     *********************************
      subroutine pair_dist(filename)
      implicit none
#include "errquit.fh"
      character*(*) filename

#include "bafdecls.fh"
#include "stdio.fh"

c     *** parameters ***
      real*8 DELR
      parameter (DELR=0.01d0)
 
      integer MASTER,taskid
      parameter (MASTER=0)

*     **** local variables ****
      logical value
      integer MAXBIN
      integer i,j,bin,natoms,natoms1,t,nframes,iii
      integer rx(2),ry(2),rz(2),hst(2)
      real*8  time1,time2,volume,volume1
      real*8  x,y,z,r,gr,v1,v2,v3
      real*8 RLOWER,RUPPER,NIDEAL,CONST    
      real*8 grs,grm2,grm1,grp1,grp2
      character*255 full_filename
      character*2 tag2
      character*4 tag4

      real*8 unita(3,3)

*     *** synchronization ***
      call ga_sync()

      call util_file_name_noprefix('PAIR_DISTRIBUTION',.false.,
     >                             .false.,
     >                             full_filename)

      call Parallel_taskid(taskid)



*     ***************************************
*     ****** PAIR DISTRIBUTION routine ******
*     ***************************************
      if (taskid.eq.MASTER) then

      write(luout,*)
      write(luout,*) 
     =  '===== Calculating pair distribution function ====='
      write(luout,*)

*     *** get the number of frames from MOTION ***
      call MOTION_nframes(filename,nframes,natoms,volume,unita)

      open(unit=19,file=filename,form='formatted')


      MAXBIN = INT(0.5d0*(volume**(0.33333333d0))/DELR)
      CONST = (16.0d0*datan(1.0d0)/3.0d0)*(natoms/volume)

      value = BA_push_get(mt_int,MAXBIN,'hst89',hst(2),hst(1))
      value = value.and.
     >        BA_push_get(mt_dbl,natoms,'rx',rx(2),rx(1))
      value = value.and.
     >        BA_push_get(mt_dbl,natoms,'ry',ry(2),ry(1))
      value = value.and.
     >        BA_push_get(mt_dbl,natoms,'rz',rz(2),rz(1))
      if (.not. value) 
     >    call errquit('cpmd_properties: out of stack memory',0, MA_ERR)
      
     
      do bin=0,MAXBIN-1
         int_mb(hst(1)+bin) = 0
      end do
    
*     **** calculate the pair distribution function ****
      time1   = 0.0d0
      do t=1,nframes
        
         time2    = time1
         read(19,*) time1,natoms1,volume1,unita
         do i=1,natoms1
            read(19,*) iii,tag2,tag4,x,y,z,v1,v2,v3
            dbl_mb(rx(1)+(i-1)) = x
            dbl_mb(ry(1)+(i-1)) = y
            dbl_mb(rz(1)+(i-1)) = z          
         end do

         do i=1,natoms1-1
         do j=i+1,natoms1
           x = dbl_mb(rx(1)+(i-1)) - dbl_mb(rx(1)+(j-1))
           y = dbl_mb(ry(1)+(i-1)) - dbl_mb(ry(1)+(j-1))
           z = dbl_mb(rz(1)+(i-1)) - dbl_mb(rz(1)+(j-1))
           r = dsqrt(x*x+y*y+z*z)
           bin = INT(r/DELR) 
           if (bin.le.(MAXBIN-1)) 
     >        int_mb(hst(1)+bin) = int_mb(hst(1)+bin) + 2

         end do
         end do

      end do
      close(unit=19)

      write(luout,*) 
      write(luout,*) 'frames used:  ',nframes
      write(luout,*) 'natoms:       ',natoms
      write(luout,*) 'delr:         ',DELR
      write(luout,*) 'volume:       ',volume
      write(luout,*) 'maxbin:       ',MAXBIN
      write(luout,*) 'MOTION filename:     ',filename
      write(luout,*) 'g(r)   filename:     ',full_filename
      write(luout,*) 
      call util_flush(luout)

*     ***** write out g(r) in PAIR_DISTRIBUTION ****
c     open(unit=19,file='PAIR_DISTRIBUTION',form='formatted')
      open(unit=19,file=full_filename,form='formatted')

***     **** normalize hist(bin) to give g(r) ****
        do bin=2,MAXBIN-3
          RLOWER = dble(bin)*DELR
          RUPPER = RLOWER + DELR
          NIDEAL = CONST*(RUPPER**3 - RLOWER**3)
          r  = 0.5d0*(RLOWER+RUPPER)
          gr = dble(int_mb(hst(1)+bin))/dble(nframes)
     >                                  /dble(natoms)
     >                                  /NIDEAL
          RLOWER = dble(bin-2)*DELR
          RUPPER = RLOWER + DELR
          NIDEAL = CONST*(RUPPER**3 - RLOWER**3)
          grm2 = dble(int_mb(hst(1)+bin-2))/dble(nframes)
     >                                     /dble(natoms)
     >                                     /NIDEAL
          RLOWER = dble(bin-1)*DELR
          RUPPER = RLOWER + DELR
          NIDEAL = CONST*(RUPPER**3 - RLOWER**3)
          grm1 = dble(int_mb(hst(1)+bin-1))/dble(nframes)
     >                                     /dble(natoms)
     >                                     /NIDEAL
          RLOWER = dble(bin+1)*DELR
          RUPPER = RLOWER + DELR
          NIDEAL = CONST*(RUPPER**3 - RLOWER**3)
          grp1 = dble(int_mb(hst(1)+bin+1))/dble(nframes)
     >                                     /dble(natoms)
     >                                     /NIDEAL
          RLOWER = dble(bin+2)*DELR
          RUPPER = RLOWER + DELR
          NIDEAL = CONST*(RUPPER**3 - RLOWER**3)
          grp2 = dble(int_mb(hst(1)+bin+2))/dble(nframes)
     >                                     /dble(natoms)
     >                                     /NIDEAL

          grs = ( -3.0d0*grm2
     >          + 12.0d0*grm1
     >          + 17.0d0*gr
     >          + 12.0d0*grp1
     >          -  3.0d0*grp2)/35.0d0

          write(19,*) r,grs,gr    
        end do
      close(unit=19)

      value =           BA_pop_stack(rz(2))
      value = value.and.BA_pop_stack(ry(2))
      value = value.and.BA_pop_stack(rx(2))
      value = value.and.BA_pop_stack(hst(2))
      if (.not. value) 
     >  call errquit('cpmd_properties: cannot release stack memory',0,
     &       MA_ERR)
     
      end if

*     *** synchronization ***
      call ga_sync()

      return
      end

