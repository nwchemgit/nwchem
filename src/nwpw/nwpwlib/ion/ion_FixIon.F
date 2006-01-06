*     *************************
*     *                       *
*     *     ion_init_FixIon       *
*     *                       *
*     *************************
      subroutine ion_init_FixIon(rtdb)
      implicit none
#include "errquit.fh"
      integer rtdb

#include "mafdecls.fh"
#include "rtdb.fh"


*     **** FixIon common block ****
      integer count
      integer ion_number(2)
      common / FixIon_Block / ion_number,count


*     **** FixIon2 common block - inverse forces ****
      integer count2
      integer ion_number1(2),ion_number1r(2)
      common / FixIon2_Block / ion_number1,ion_number1r,count2


      integer taskid
      integer MASTER
      parameter(MASTER=0)

      integer i,j,nion,nconstraints
      integer ma_type,nactive_atoms,h_actlist,l_actlist
      logical found

*     **** external functions ****
      integer  ion_nionall
      external ion_nionall
      logical  control_zero_forces
      external control_zero_forces

*     *** initialize the  number of constraints ****
      count = 0

*     **** get the nions, including charges and qmmm ****
      nion = ion_nionall()

*     ************************************************
*     **** all ions constraints have been imposed ****
*     ************************************************
      if (control_zero_forces()) then
        nconstraints = nion

*       ***** allocate constraint list, ion_number, from heap ****
        if (.not. MA_alloc_get(mt_int,nconstraints,
     >     'ion_number',ion_number(2),ion_number(1)))
     >      call errquit('out of heap memory',0, MA_ERR)

*       ***** set the constraints ****
        do i=1,nion
           count = count+1
           int_mb(ion_number(1)+count-1) = i
        end do

*     ********************************************************
*     **** geometry:actlist constraints have been imposed ****
*     ********************************************************
      else if (rtdb_ma_get(rtdb, 'geometry:actlist', ma_type,
     >        nactive_atoms, h_actlist)) then

         if (.not.MA_get_index(h_actlist,l_actlist)) 
     >      call errquit(
     >       'ion_init_FixIon: ma_get_index failed for actlist',911,
     &       MA_ERR)


*       ***** count the number of constraints *****
        nconstraints = 0
        do i=1,nion
           found = .false.
           do j=1,nactive_atoms
              if (i .eq. int_mb(l_actlist+j-1)) found = .true.
           end do

*          **** constraint added ****
           if (.not. found) then
              nconstraints = nconstraints + 1
           end if
        end do

        if (nconstraints.gt.0) then

*          ***** allocate constraint list, ion_number, from heap ****
           if (.not. MA_alloc_get(mt_int,nconstraints,
     >       'ion_number',ion_number(2),ion_number(1)))
     >         call errquit('out of heap memory',0, MA_ERR)

        

*          ***** set the constraints ****
           do i=1,nion
              found = .false.
              do j=1,nactive_atoms
                 if (i .eq. int_mb(l_actlist+j-1)) found = .true.
              end do

*             **** constraint added ****
              if (.not. found) then
                 count = count+1
                 int_mb(ion_number(1)+count-1) = i
              end if
           end do
        end if

*       *** free heap ***
        if (.not. MA_free_heap(h_actlist)) 
     >   call errquit('h_actlist:error freeing heap memory',0, MA_ERR)

      end if


  102 continue
      if (count.ne.0) then
         call Parallel_taskid(taskid)
         if (taskid.eq.MASTER) then
          write(*,*)
          write(*,1003) (int_mb(ion_number(1)+i-1),i=1,count)
          write(*,*)
         end if
      end if




      if (rtdb_ma_get(rtdb, 'geometry:inversionlist', ma_type,
     >        nconstraints, h_actlist)) then

         if (.not.MA_get_index(h_actlist,l_actlist))
     >      call errquit(
     >       'ion_init_FixIon: ma_get_index failed for actlist',911,
     &       MA_ERR)



           count2 = nconstraints/2
*          ***** allocate constraint list, ion_number, from heap ****
           if (.not. MA_alloc_get(mt_int,count2,
     >       'ion_number1',ion_number1(2),ion_number1(1)))
     >         call errquit('out of heap memory',0, MA_ERR)
           if (.not. MA_alloc_get(mt_int,count2,
     >       'ion_number1r',ion_number1r(2),ion_number1r(1)))
     >         call errquit('out of heap memory',0, MA_ERR)

           do i=1,count2
              int_mb(ion_number1(1)+i-1)  = int_mb(l_actlist+2*i-1-1)
              int_mb(ion_number1r(1)+i-1) = int_mb(l_actlist+2*i-1)
           end do



*       *** free heap ***
        if (.not. MA_free_heap(h_actlist))
     >   call errquit('h_actlist:error freeing heap memory',0, MA_ERR)
      end if


      if (count2.ne.0) then
         call Parallel_taskid(taskid)
         if (taskid.eq.MASTER) then
          write(*,*) 
          write(*,1004) (int_mb(ion_number1(1)+i-1),
     >                   int_mb(ion_number1r(1)+i-1), i=1,count2)
          write(*,*)
         end if
      end if
   
*     *** initialize molecular constraints ****
      call molecular_init(rtdb)

      return
 1003 FORMAT(1x,"Fixed ion positions:",10I5)
 1004 FORMAT(1x,"Inversion constraints:",10I5)
      end


*     ***************************
*     *                         *
*     *         ion_FixIon      *
*     *                         *
*     ***************************

      subroutine ion_FixIon(force)
      implicit none
      real*8 force(3,*)

#include "mafdecls.fh"

*     **** FixIon common block ****
      integer count
      integer ion_number(2)
      common / FixIon_Block / ion_number,count

*     **** FixIon2 common block - inverse forces ****
      integer count2
      integer ion_number1(2),ion_number1r(2)
      common / FixIon2_Block / ion_number1,ion_number1r,count2



*     **** local variables ****
      integer i

      if (count.ne.0) then
         do i=1,count
            force(1,int_mb(ion_number(1)+i-1)) = 0.0d0
            force(2,int_mb(ion_number(1)+i-1)) = 0.0d0
            force(3,int_mb(ion_number(1)+i-1)) = 0.0d0
         end do
      end if


      if (count2.ne.0) then
         do i=1,count2
            force(1,int_mb(ion_number1r(1)+i-1)) =
     >     -force(1,int_mb(ion_number1(1)+i-1))
            force(2,int_mb(ion_number1r(1)+i-1)) =
     >     -force(2,int_mb(ion_number1(1)+i-1))
            force(3,int_mb(ion_number1r(1)+i-1)) = 
     >     -force(3,int_mb(ion_number1(1)+i-1))
         end do
      end if


      return
      end
*     ********************************
*     *                              *
*     *         ion_same_FixIon      *
*     *                              *
*     ********************************

*   make r2 = r1 for fixed ions
*
      subroutine ion_same_FixIon(r1,r2)
      implicit none
      real*8 r1(3,*)
      real*8 r2(3,*)

#include "mafdecls.fh"

*     **** FixIon common block ****
      integer count
      integer ion_number(2)
      common / FixIon_Block / ion_number,count

*     **** local variables ****
      integer i

      if (count.ne.0) then
         do i=1,count
            r2(1,int_mb(ion_number(1)+i-1)) 
     >      = r1(1,int_mb(ion_number(1)+i-1)) 
            r2(2,int_mb(ion_number(1)+i-1))
     >      = r1(2,int_mb(ion_number(1)+i-1)) 
            r2(3,int_mb(ion_number(1)+i-1))
     >      = r1(3,int_mb(ion_number(1)+i-1)) 
         end do
      end if

      return
      end



*     *************************
*     *                       *
*     *     ion_end_FixIon    *
*     *                       *
*     *************************
      subroutine ion_end_FixIon()
      implicit none
#include "errquit.fh"

#include "mafdecls.fh"

*     **** FixIon common block ****
      integer count
      integer ion_number(2)
      common / FixIon_Block / ion_number,count

*     **** FixIon2 common block - inverse forces ****
      integer count2
      integer ion_number1(2),ion_number1r(2)
      common / FixIon2_Block / ion_number1,ion_number1r,count2



      logical value

      if (count.gt.0) then
        value =  MA_free_heap(ion_number(2)) 
        if (.not. value)
     >  call errquit('ion_number(2):error freeing heap memory',0,
     &       MA_ERR)
      end if
      if (count2.gt.0) then
        value =            MA_free_heap(ion_number1(2)) 
        value =  value.and.MA_free_heap(ion_number1r(2)) 
        if (.not. value)
     >  call errquit('ion_number1(2):error freeing heap memory',0,
     &       MA_ERR)
      end if

*     *** remove molecular constraints ****
      call molecular_end()

      return
      end


*     ***************************
*     *                         *
*     *         ion_q_FixIon    *
*     *                         *
*     ***************************

      logical function ion_q_FixIon(i)
      implicit none
      integer i

#include "mafdecls.fh"

*     **** FixIon common block ****
      integer count
      integer ion_number(2)
      common / FixIon_Block / ion_number,count

*     **** local variables ***
      logical found
      integer j

      found = .false.
      do j=1,count
         if (i .eq. int_mb(ion_number(1)+j-1)) found = .true.
      end do

      ion_q_FixIon = found
      return
      end
