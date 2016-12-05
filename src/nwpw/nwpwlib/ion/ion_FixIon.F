*     *************************
*     *                       *
*     *     ion_init_FixIon   *
*     *                       *
*     *************************
      subroutine ion_init_FixIon(rtdb)
      implicit none
      integer rtdb

#include "bafdecls.fh"
#include "btdb.fh"
#include "errquit.fh"


*     **** FixIon common block ****
      integer count
      integer ion_number(2)
      common / FixIon_Block / ion_number,count

*     **** FixIonz common block ****
      integer countx,county,countz
      integer ion_numberx(2),ion_numbery(2),ion_numberz(2)
      common / FixIonz_Block / ion_numberx,ion_numbery,ion_numberz,
     >                         countx,county,countz


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
      integer  ion_nion
      external ion_nion
      logical  control_zero_forces
      external control_zero_forces

*     *** initialize the  number of constraints ****
      count  = 0
      count2 = 0
      countx = 0
      county = 0
      countz = 0

*     **** get the nions, including charges and qmmm ****
      nion = ion_nion()

*     ************************************************
*     **** all ions constraints have been imposed ****
*     ************************************************
      if (control_zero_forces()) then
        nconstraints = nion

*       ***** allocate constraint list, ion_number, from heap ****
        if (.not. BA_alloc_get(mt_int,nconstraints,
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

         if (.not.BA_get_index(h_actlist,l_actlist)) 
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
           if (.not. BA_alloc_get(mt_int,nconstraints,
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
        if (.not. BA_free_heap(h_actlist)) 
     >   call errquit('h_actlist:error freeing heap memory',0, MA_ERR)

*     ********************************************************
*     ****   nwpw:actlist constraints have been imposed   ****
*     ********************************************************
      else if (rtdb_ma_get(rtdb, 'nwpw:actlist', ma_type,
     >        nactive_atoms, h_actlist)) then

         if (.not.BA_get_index(h_actlist,l_actlist))
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
           if (.not. BA_alloc_get(mt_int,nconstraints,
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
        if (.not. BA_free_heap(h_actlist))
     >   call errquit('h_actlist:error freeing heap memory',0, MA_ERR)


      end if


  102 continue
      if (count.ne.0) then
         call ion_add_constraint(3*count)
         call Parallel_taskid(taskid)
         if (taskid.eq.MASTER) then
          write(*,*)
          write(*,1003) (int_mb(ion_number(1)+i-1),i=1,count)
          write(*,*)
         end if
      end if

*     ******************************************************************
*     **** nwpw:xfreezelist constraints have been imposed           ****
*     ******************************************************************
      if (rtdb_ma_get(rtdb, 'nwpw:xfreezelist', ma_type,
     >        nconstraints, h_actlist)) then

         if (.not.BA_get_index(h_actlist,l_actlist))
     >      call errquit(
     >       'ion_init_FixIon: ma_get_index failed for actlist',911,
     >       MA_ERR)

           countx = nconstraints
*          ***** allocate constraint list, ion_number, from heap ****
           if (.not. BA_alloc_get(mt_int,nconstraints,
     >       'ion_numberx',ion_numberx(2),ion_numberx(1)))
     >         call errquit('out of heap memory',0, MA_ERR)

           do i=1,countx
              int_mb(ion_numberx(1)+i-1)  = int_mb(l_actlist+i-1)
           end do

*        *** free heap ***
         if (.not. BA_free_heap(h_actlist))
     >   call errquit('h_actlist:error freeing heap memory',0, MA_ERR)
      end if

      if (countx.ne.0) then
         call ion_add_constraint(countx)
         call Parallel_taskid(taskid)
         if (taskid.eq.MASTER) then
          write(*,*)
          write(*,1005) (int_mb(ion_numberx(1)+i-1),i=1,countx)
          write(*,*)
         end if
      end if

*     ******************************************************************
*     **** nwpw:yfreezelist constraints have been imposed           ****
*     ******************************************************************
      if (rtdb_ma_get(rtdb, 'nwpw:yfreezelist', ma_type,
     >        nconstraints, h_actlist)) then

         if (.not.BA_get_index(h_actlist,l_actlist))
     >      call errquit(
     >       'ion_init_FixIon: ma_get_index failed for actlist',911,
     >       MA_ERR)

           county = nconstraints
*          ***** allocate constraint list, ion_number, from heap ****
           if (.not. BA_alloc_get(mt_int,nconstraints,
     >       'ion_numbery',ion_numbery(2),ion_numbery(1)))
     >         call errquit('out of heap memory',0, MA_ERR)

           do i=1,county
              int_mb(ion_numbery(1)+i-1)  = int_mb(l_actlist+i-1)
           end do

*        *** free heap ***
         if (.not. BA_free_heap(h_actlist))
     >   call errquit('h_actlist:error freeing heap memory',0, MA_ERR)
      end if

      if (county.ne.0) then
         call ion_add_constraint(county)
         call Parallel_taskid(taskid)
         if (taskid.eq.MASTER) then
          write(*,*)
          write(*,1006) (int_mb(ion_numbery(1)+i-1),i=1,county)
          write(*,*)
         end if
      end if



*     ******************************************************************
*     **** nwpw:zfreezelist constraints have been imposed           ****
*     ******************************************************************
      if (rtdb_ma_get(rtdb, 'nwpw:zfreezelist', ma_type,
     >        nconstraints, h_actlist)) then

         if (.not.BA_get_index(h_actlist,l_actlist))
     >      call errquit(
     >       'ion_init_FixIon: ma_get_index failed for actlist',911,
     >       MA_ERR)

           countz = nconstraints
*          ***** allocate constraint list, ion_number, from heap ****
           if (.not. BA_alloc_get(mt_int,nconstraints,
     >       'ion_numberz',ion_numberz(2),ion_numberz(1)))
     >         call errquit('out of heap memory',0, MA_ERR)

           do i=1,countz
              int_mb(ion_numberz(1)+i-1)  = int_mb(l_actlist+i-1)
           end do

*        *** free heap ***
         if (.not. BA_free_heap(h_actlist))
     >   call errquit('h_actlist:error freeing heap memory',0, MA_ERR)
      end if

      if (countz.ne.0) then
         call ion_add_constraint(countz)
         call Parallel_taskid(taskid)
         if (taskid.eq.MASTER) then
          write(*,*)
          write(*,1007) (int_mb(ion_numberz(1)+i-1),i=1,countz)
          write(*,*)
         end if
      end if




*     ******************************************************************
*     **** geometry:inversionlistlist constraints have been imposed ****
*     ******************************************************************
      if (rtdb_ma_get(rtdb, 'geometry:inversionlist', ma_type,
     >        nconstraints, h_actlist)) then

         if (.not.BA_get_index(h_actlist,l_actlist))
     >      call errquit(
     >       'ion_init_FixIon: ma_get_index failed for actlist',911,
     &       MA_ERR)



           count2 = nconstraints/2
*          ***** allocate constraint list, ion_number, from heap ****
           if (.not. BA_alloc_get(mt_int,count2,
     >       'ion_number1',ion_number1(2),ion_number1(1)))
     >         call errquit('out of heap memory',0, MA_ERR)
           if (.not. BA_alloc_get(mt_int,count2,
     >       'ion_number1r',ion_number1r(2),ion_number1r(1)))
     >         call errquit('out of heap memory',0, MA_ERR)

           do i=1,count2
              int_mb(ion_number1(1)+i-1)  = int_mb(l_actlist+2*i-1-1)
              int_mb(ion_number1r(1)+i-1) = int_mb(l_actlist+2*i-1)
           end do



*       *** free heap ***
        if (.not. BA_free_heap(h_actlist))
     >   call errquit('h_actlist:error freeing heap memory',0, MA_ERR)
      end if


      if (count2.ne.0) then
         call ion_add_constraint(3*count2)
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
 1005 FORMAT(1x,"Fixed ion x-positions:",10I5)
 1006 FORMAT(1x,"Fixed ion y-positions:",10I5)
 1007 FORMAT(1x,"Fixed ion z-positions:",10I5)
      end


*     ***************************
*     *                         *
*     *         ion_FixIon      *
*     *                         *
*     ***************************

      subroutine ion_FixIon(force)
      implicit none
      real*8 force(3,*)

#include "bafdecls.fh"

*     **** FixIon common block ****
      integer count
      integer ion_number(2)
      common / FixIon_Block / ion_number,count

*     **** FixIonz common block ****
      integer countx,county,countz
      integer ion_numberx(2),ion_numbery(2),ion_numberz(2)
      common / FixIonz_Block / ion_numberx,ion_numbery,ion_numberz,
     >                         countx,county,countz


*     **** FixIon2 common block - inverse forces ****
      integer count2
      integer ion_number1(2),ion_number1r(2)
      common / FixIon2_Block / ion_number1,ion_number1r,count2



*     **** local variables ****
      integer i

      if (count.ne.0) then
!$OMP DO
         do i=1,count
            force(1,int_mb(ion_number(1)+i-1)) = 0.0d0
            force(2,int_mb(ion_number(1)+i-1)) = 0.0d0
            force(3,int_mb(ion_number(1)+i-1)) = 0.0d0
         end do
!$OMP END DO
      end if

      if (countx.ne.0) then
!$OMP DO
         do i=1,countx
            force(1,int_mb(ion_numberx(1)+i-1)) = 0.0d0
         end do
!$OMP END DO
      end if

      if (county.ne.0) then
!$OMP DO
         do i=1,county
            force(2,int_mb(ion_numbery(1)+i-1)) = 0.0d0
         end do
!$OMP END DO
      end if

      if (countz.ne.0) then
!$OMP DO
         do i=1,countz
            force(3,int_mb(ion_numberz(1)+i-1)) = 0.0d0
         end do
!$OMP END DO
      end if

      if (count2.ne.0) then
!$OMP DO
         do i=1,count2
            force(1,int_mb(ion_number1r(1)+i-1)) =
     >     -force(1,int_mb(ion_number1(1)+i-1))
            force(2,int_mb(ion_number1r(1)+i-1)) =
     >     -force(2,int_mb(ion_number1(1)+i-1))
            force(3,int_mb(ion_number1r(1)+i-1)) = 
     >     -force(3,int_mb(ion_number1(1)+i-1))
         end do
!$OMP END DO
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

#include "bafdecls.fh"

*     **** FixIon common block ****
      integer count
      integer ion_number(2)
      common / FixIon_Block / ion_number,count

*     **** local variables ****
      integer i

      if (count.ne.0) then
!$OMP DO
         do i=1,count
            r2(1,int_mb(ion_number(1)+i-1)) 
     >      = r1(1,int_mb(ion_number(1)+i-1)) 
            r2(2,int_mb(ion_number(1)+i-1))
     >      = r1(2,int_mb(ion_number(1)+i-1)) 
            r2(3,int_mb(ion_number(1)+i-1))
     >      = r1(3,int_mb(ion_number(1)+i-1)) 
         end do
!$OMP END DO
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

#include "bafdecls.fh"

*     **** FixIon common block ****
      integer count
      integer ion_number(2)
      common / FixIon_Block / ion_number,count

*     **** FixIonz common block ****
      integer countx,county,countz
      integer ion_numberx(2),ion_numbery(2),ion_numberz(2)
      common / FixIonz_Block / ion_numberx,ion_numbery,ion_numberz,
     >                         countx,county,countz

*     **** FixIon2 common block - inverse forces ****
      integer count2
      integer ion_number1(2),ion_number1r(2)
      common / FixIon2_Block / ion_number1,ion_number1r,count2

      logical value

      if (count.gt.0) then
        value =  BA_free_heap(ion_number(2)) 
        if (.not. value)
     >  call errquit('ion_number(2):error freeing heap memory',0,
     &       MA_ERR)
      end if

      if (countx.gt.0) then
        value =  BA_free_heap(ion_numberx(2)) 
        if (.not. value)
     >  call errquit('ion_numberx(2):error freeing heap memory',0,
     >       MA_ERR)
      end if

      if (county.gt.0) then
        value =  BA_free_heap(ion_numbery(2)) 
        if (.not. value)
     >  call errquit('ion_numbery(2):error freeing heap memory',0,
     >       MA_ERR)
      end if

      if (countz.gt.0) then
        value =  BA_free_heap(ion_numberz(2)) 
        if (.not. value)
     >  call errquit('ion_numberz(2):error freeing heap memory',0,
     >       MA_ERR)
      end if


      if (count2.gt.0) then
        value =            BA_free_heap(ion_number1(2)) 
        value =  value.and.BA_free_heap(ion_number1r(2)) 
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

#include "bafdecls.fh"

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



*     ***************************
*     *                         *
*     *       ion_q_xFixIon     *
*     *                         *
*     ***************************
      logical function ion_q_xFixIon(i)
      implicit none
      integer i

#include "bafdecls.fh"

*     **** FixIonz common block ****
      integer countx,county,countz
      integer ion_numberx(2),ion_numbery(2),ion_numberz(2)
      common / FixIonz_Block / ion_numberx,ion_numbery,ion_numberz,
     >                         countx,county,countz

*     **** local variables ***
      logical found
      integer j

      found = .false.
      do j=1,countx
         if (i .eq. int_mb(ion_numberx(1)+j-1)) found = .true.
      end do

      ion_q_xFixIon = found
      return
      end

*     ***************************
*     *                         *
*     *       ion_q_yFixIon     *
*     *                         *
*     ***************************
      logical function ion_q_yFixIon(i)
      implicit none
      integer i

#include "bafdecls.fh"

*     **** FixIonz common block ****
      integer countx,county,countz
      integer ion_numberx(2),ion_numbery(2),ion_numberz(2)
      common / FixIonz_Block / ion_numberx,ion_numbery,ion_numberz,
     >                         countx,county,countz

*     **** local variables ***
      logical found
      integer j

      found = .false.
      do j=1,county
         if (i .eq. int_mb(ion_numbery(1)+j-1)) found = .true.
      end do

      ion_q_yFixIon = found
      return
      end


*     ***************************
*     *                         *
*     *       ion_q_zFixIon     *
*     *                         *
*     ***************************
      logical function ion_q_zFixIon(i)
      implicit none
      integer i

#include "bafdecls.fh"

*     **** FixIonz common block ****
      integer countx,county,countz
      integer ion_numberx(2),ion_numbery(2),ion_numberz(2)
      common / FixIonz_Block / ion_numberx,ion_numbery,ion_numberz,
     >                         countx,county,countz

*     **** local variables ***
      logical found
      integer j

      found = .false.
      do j=1,countz
         if (i .eq. int_mb(ion_numberz(1)+j-1)) found = .true.
      end do

      ion_q_zFixIon = found
      return
      end

*     ***************************
*     *                         *
*     *     ion_q_xyzFixIon     *
*     *                         *
*     ***************************
      logical function ion_q_xyzFixIon(i)
      implicit none
      integer i

#include "bafdecls.fh"

*     **** FixIonz common block ****
      integer countx,county,countz
      integer ion_numberx(2),ion_numbery(2),ion_numberz(2)
      common / FixIonz_Block / ion_numberx,ion_numbery,ion_numberz,
     >                         countx,county,countz

*     **** local variables ***
      logical found
      integer j

      found = .false.
      do j=1,countx
         if (i .eq. int_mb(ion_numberx(1)+j-1)) found = .true.
      end do
      do j=1,county
         if (i .eq. int_mb(ion_numbery(1)+j-1)) found = .true.
      end do
      do j=1,countz
         if (i .eq. int_mb(ion_numberz(1)+j-1)) found = .true.
      end do

      ion_q_xyzFixIon = found
      return
      end

*     *********************************
*     *                               *
*     *     ion_q_xyzFixIon_label     *
*     *                               *
*     *********************************
      character*14 function ion_q_xyzFixIon_label(i)
      implicit none
      integer i

#include "bafdecls.fh"

*     **** FixIonz common block ****
      integer countx,county,countz
      integer ion_numberx(2),ion_numbery(2),ion_numberz(2)
      common / FixIonz_Block / ion_numberx,ion_numbery,ion_numberz,
     >                         countx,county,countz

*     **** local variables ***
      logical foundx,foundy,foundz
      integer j
      character*14 label

      foundx = .false.
      foundy = .false.
      foundz = .false.
      do j=1,countx
         if (i .eq. int_mb(ion_numberx(1)+j-1)) foundx = .true.
      end do
      do j=1,county
         if (i .eq. int_mb(ion_numbery(1)+j-1)) foundy = .true.
      end do
      do j=1,countz
         if (i .eq. int_mb(ion_numberz(1)+j-1)) foundz = .true.
      end do
      
      label = "              "
      if      (foundx.and.foundy.and.foundz) then
         label = " - x,y,z fixed"
      else if (foundx.and.foundy) then
         label = " - x,y   fixed"
      else if (foundx.and.foundz) then
         label = " - x,z   fixed"
      else if (foundy.and.foundz) then
         label = " - y,z   fixed"
      else if (foundx) then
         label = " - x     fixed"
      else if (foundy) then
         label = " - y f   fixed"
      else if (foundz) then
         label = " - z     fixed"
      end if

      ion_q_xyzFixIon_label = label
      return
      end 




*     ***********************************
*     *				        *
*     *	    xFixIon_force_average	*
*     *	 		                *
*     ***********************************
      subroutine xFixIon_force_average(filename)
      implicit none
      character*(*) filename

#include "bafdecls.fh"
#include "errquit.fh"

*     **** FixIonz common block ****
      integer countx,county,countz
      integer ion_numberx(2),ion_numbery(2),ion_numberz(2)
      common / FixIonz_Block / ion_numberx,ion_numbery,ion_numberz,
     >                         countx,county,countz


*     **** local variables ****
      logical ok
      integer nion
      integer m,indx,indx2,indx3,ii,iii,nframe,l
      double precision E
      double precision a1x,a1y,a1z
      double precision a2x,a2y,a2z
      double precision a3x,a3y,a3z
      integer r1(2),f1(2),af1(2),af2(2),ac1(2)
      character*2 tag2
      character*4 tag4

      if (countx.gt.0) then
         l = index(filename,' ') - 1
         write(*,*)
         write(*,*) "Force Energy File:",filename(1:l)

         open(unit=35,file=filename,form='formatted',status='old')
         read(35,*) nion
         close(unit=35)

         ok =        BA_push_get(mt_dbl,3*nion,'r1',r1(2),r1(1))
         ok = ok.and.BA_push_get(mt_dbl,3*nion,'f1',f1(2),f1(1))
         ok = ok.and.BA_push_get(mt_dbl,countx,'af1',af1(2),af1(1))
         ok = ok.and.BA_push_get(mt_dbl,countx,'af2',af2(2),af2(1))
         ok = ok.and.BA_push_get(mt_dbl,countx,'ac1',ac1(2),ac1(1))
         if (.not.ok)
     >    call errquit('xFixIon_force_average: out of stack',1,MA_ERR)
         call dcopy(countx,0.0d0,0,dbl_mb(af1(1)),1)
         call dcopy(countx,0.0d0,0,dbl_mb(af2(1)),1)
         call dcopy(countx,0.0d0,0,dbl_mb(ac1(1)),1)

         open(unit=35,file=filename,form='formatted',status='old')
         nframe = 0
         do while (.true.)
            read(35,*,ERR=30,END=30) nion
            read(35,*,ERR=30,END=30) E
            read(35,*,ERR=30,END=30) a1x,a1y,a1z
            read(35,*,ERR=30,END=30) a2x,a2y,a2z
            read(35,*,ERR=30,END=30) a3x,a3y,a3z
            indx = 0
            do ii=1,nion
              read(35,*,ERR=30,END=30) iii,tag2,tag4,
     >                                 dbl_mb(r1(1)+indx),
     >                                 dbl_mb(r1(1)+indx+1),
     >                                 dbl_mb(r1(1)+indx+2),
     >                                 dbl_mb(f1(1)+indx),
     >                                 dbl_mb(f1(1)+indx+1),
     >                                 dbl_mb(f1(1)+indx+2)
              indx = indx + 3
            end do

            nframe = nframe + 1
            call xFixIon_sub_force_average(dbl_mb(r1(1)),
     >                                     dbl_mb(f1(1)),
     >                                     countx,
     >                                     int_mb(ion_numberx(1)), 
     >                                     dbl_mb(af1(1)),
     >                                     dbl_mb(af2(1)),
     >                                     dbl_mb(ac1(1)))
         end do

  30     continue
         close(unit=35)

         if (nframe.eq.0) then
            E = 0.0d0
         else
            E = 1.0d0/dble(nframe)
         end if
         call dscal(countx,E,dbl_mb(af1(1)),1)
         call dscal(countx,E,dbl_mb(af2(1)),1)
         call dscal(countx,E,dbl_mb(ac1(1)),1)

         write(*,110) nframe
         write(*,111) countx
         write(*,112) nion
  110    format(1x,'frames used           =',I8)
  111    format(1x,'number of constraints =',I8)
  112    format(1x,'number of atoms       =',I8)
         do m=1,countx
            write(*,120) int_mb(ion_numberx(1)+m-1),
     >                   dbl_mb(ac1(1)+m-1),
     >                   dbl_mb(af1(1)+m-1),
     >                   dbl_mb(af2(1)+m-1)-dbl_mb(af1(1)+m-1)**2
         end do
  120    format(1x,'atom=',I5,
     >          1x,'<x>=', E14.6,
     >          4x,'<Fx>=',E14.6,
     >          1x,'(<F**2>-<F>**2=',E10.3,')')

         ok =        BA_pop_stack(ac1(2))
         ok = ok.and.BA_pop_stack(af2(2))
         ok = ok.and.BA_pop_stack(af1(2))
         ok = ok.and.BA_pop_stack(f1(2))
         ok = ok.and.BA_pop_stack(r1(2))
         if (.not.ok)
     >    call errquit(
     >         'xFixIon_force_average: error popping stack',1,MA_ERR)

      end if

      return
      end

      subroutine xFixIon_sub_force_average(rion,fion,
     >                                     nindx,indx,
     >                                     af,af2,ac)
      implicit none
      real*8 rion(3,*)
      real*8 fion(3,*)
      integer nindx,indx(*)
      real*8 af(*)
      real*8 af2(*)
      real*8 ac(*)

*     **** local variables
      integer m
      real*8 ff

      do m=1,nindx
         ff = fion(1,indx(m)) 
         af(m)  = af(m)  + ff
         af2(m) = af2(m) + ff*ff
         ac(m)  = ac(m)  + rion(1,indx(m))
      end do

      return
      end


*     ***********************************
*     *				        *
*     *	    yFixIon_force_average	*
*     *	 		                *
*     ***********************************
      subroutine yFixIon_force_average(filename)
      implicit none
      character*(*) filename

#include "bafdecls.fh"
#include "errquit.fh"

*     **** FixIonz common block ****
      integer countx,county,countz
      integer ion_numberx(2),ion_numbery(2),ion_numberz(2)
      common / FixIonz_Block / ion_numberx,ion_numbery,ion_numberz,
     >                         countx,county,countz


*     **** local variables ****
      logical ok
      integer nion
      integer m,indx,indx2,indx3,ii,iii,nframe,l
      double precision E
      double precision a1x,a1y,a1z
      double precision a2x,a2y,a2z
      double precision a3x,a3y,a3z
      integer r1(2),f1(2),af1(2),af2(2),ac1(2)
      character*2 tag2
      character*4 tag4

      if (county.gt.0) then
         l = index(filename,' ') - 1
         write(*,*)
         write(*,*) "Force Energy File:",filename(1:l)

         open(unit=35,file=filename,form='formatted',status='old')
         read(35,*) nion
         close(unit=35)

         ok =        BA_push_get(mt_dbl,3*nion,'r1',r1(2),r1(1))
         ok = ok.and.BA_push_get(mt_dbl,3*nion,'f1',f1(2),f1(1))
         ok = ok.and.BA_push_get(mt_dbl,county,'af1',af1(2),af1(1))
         ok = ok.and.BA_push_get(mt_dbl,county,'af2',af2(2),af2(1))
         ok = ok.and.BA_push_get(mt_dbl,county,'ac1',ac1(2),ac1(1))
         if (.not.ok)
     >    call errquit('yFixIon_force_average: out of stack',1,MA_ERR)
         call dcopy(county,0.0d0,0,dbl_mb(af1(1)),1)
         call dcopy(county,0.0d0,0,dbl_mb(af2(1)),1)
         call dcopy(county,0.0d0,0,dbl_mb(ac1(1)),1)

         open(unit=35,file=filename,form='formatted',status='old')
         nframe = 0
         do while (.true.)
            read(35,*,ERR=30,END=30) nion
            read(35,*,ERR=30,END=30) E
            read(35,*,ERR=30,END=30) a1x,a1y,a1z
            read(35,*,ERR=30,END=30) a2x,a2y,a2z
            read(35,*,ERR=30,END=30) a3x,a3y,a3z
            indx = 0
            do ii=1,nion
              read(35,*,ERR=30,END=30) iii,tag2,tag4,
     >                                 dbl_mb(r1(1)+indx),
     >                                 dbl_mb(r1(1)+indx+1),
     >                                 dbl_mb(r1(1)+indx+2),
     >                                 dbl_mb(f1(1)+indx),
     >                                 dbl_mb(f1(1)+indx+1),
     >                                 dbl_mb(f1(1)+indx+2)
              indx = indx + 3
            end do

            nframe = nframe + 1
            call yFixIon_sub_force_average(dbl_mb(r1(1)),
     >                                     dbl_mb(f1(1)),
     >                                     county,
     >                                     int_mb(ion_numberx(1)), 
     >                                     dbl_mb(af1(1)),
     >                                     dbl_mb(af2(1)),
     >                                     dbl_mb(ac1(1)))
         end do

  30     continue
         close(unit=35)

         if (nframe.eq.0) then
            E = 0.0d0
         else
            E = 1.0d0/dble(nframe)
         end if
         call dscal(county,E,dbl_mb(af1(1)),1)
         call dscal(county,E,dbl_mb(af2(1)),1)
         call dscal(county,E,dbl_mb(ac1(1)),1)

         write(*,110) nframe
         write(*,111) county
         write(*,112) nion
  110    format(1x,'frames used           =',I8)
  111    format(1x,'number of constraints =',I8)
  112    format(1x,'number of atoms       =',I8)
         do m=1,county
            write(*,120) int_mb(ion_numbery(1)+m-1),
     >                   dbl_mb(ac1(1)+m-1),
     >                   dbl_mb(af1(1)+m-1),
     >                   dbl_mb(af2(1)+m-1)-dbl_mb(af1(1)+m-1)**2
         end do
  120    format(1x,'atom=',I5,
     >          1x,'<y>=', E14.6,
     >          4x,'<Fy>=',E14.6,
     >          1x,'(<F**2>-<F>**2=',E10.3,')')

         ok =        BA_pop_stack(ac1(2))
         ok = ok.and.BA_pop_stack(af2(2))
         ok = ok.and.BA_pop_stack(af1(2))
         ok = ok.and.BA_pop_stack(f1(2))
         ok = ok.and.BA_pop_stack(r1(2))
         if (.not.ok)
     >    call errquit(
     >         'xFixIon_force_average: error popping stack',1,MA_ERR)

      end if

      return
      end

      subroutine yFixIon_sub_force_average(rion,fion,
     >                                     nindx,indx,
     >                                     af,af2,ac)
      implicit none
      real*8 rion(3,*)
      real*8 fion(3,*)
      integer nindx,indx(*)
      real*8 af(*)
      real*8 af2(*)
      real*8 ac(*)

*     **** local variables
      integer m
      real*8 ff

      do m=1,nindx
         ff = fion(2,indx(m)) 
         af(m)  = af(m)  + ff
         af2(m) = af2(m) + ff*ff
         ac(m)  = ac(m)  + rion(2,indx(m))
      end do

      return
      end



*     ***********************************
*     *				        *
*     *	    zFixIon_force_average	*
*     *	 		                *
*     ***********************************
      subroutine ZFixIon_force_average(filename)
      implicit none
      character*(*) filename

#include "bafdecls.fh"
#include "errquit.fh"

*     **** FixIonz common block ****
      integer countx,county,countz
      integer ion_numberx(2),ion_numbery(2),ion_numberz(2)
      common / FixIonz_Block / ion_numberx,ion_numbery,ion_numberz,
     >                         countx,county,countz


*     **** local variables ****
      logical ok
      integer nion
      integer m,indx,indx2,indx3,ii,iii,nframe,l
      double precision E
      double precision a1x,a1y,a1z
      double precision a2x,a2y,a2z
      double precision a3x,a3y,a3z
      integer r1(2),f1(2),af1(2),af2(2),ac1(2)
      character*2 tag2
      character*4 tag4

      if (countz.gt.0) then
         l = index(filename,' ') - 1
         write(*,*)
         write(*,*) "Force Energy File:",filename(1:l)

         open(unit=35,file=filename,form='formatted',status='old')
         read(35,*) nion
         close(unit=35)

         ok =        BA_push_get(mt_dbl,3*nion,'r1',r1(2),r1(1))
         ok = ok.and.BA_push_get(mt_dbl,3*nion,'f1',f1(2),f1(1))
         ok = ok.and.BA_push_get(mt_dbl,countz,'af1',af1(2),af1(1))
         ok = ok.and.BA_push_get(mt_dbl,countz,'af2',af2(2),af2(1))
         ok = ok.and.BA_push_get(mt_dbl,countz,'ac1',ac1(2),ac1(1))
         if (.not.ok)
     >    call errquit('zFixIon_force_average: out of stack',1,MA_ERR)
         call dcopy(countz,0.0d0,0,dbl_mb(af1(1)),1)
         call dcopy(countz,0.0d0,0,dbl_mb(af2(1)),1)
         call dcopy(countz,0.0d0,0,dbl_mb(ac1(1)),1)

         open(unit=35,file=filename,form='formatted',status='old')
         nframe = 0
         do while (.true.)
            read(35,*,ERR=30,END=30) nion
            read(35,*,ERR=30,END=30) E
            read(35,*,ERR=30,END=30) a1x,a1y,a1z
            read(35,*,ERR=30,END=30) a2x,a2y,a2z
            read(35,*,ERR=30,END=30) a3x,a3y,a3z
            indx = 0
            do ii=1,nion
              read(35,*,ERR=30,END=30) iii,tag2,tag4,
     >                                 dbl_mb(r1(1)+indx),
     >                                 dbl_mb(r1(1)+indx+1),
     >                                 dbl_mb(r1(1)+indx+2),
     >                                 dbl_mb(f1(1)+indx),
     >                                 dbl_mb(f1(1)+indx+1),
     >                                 dbl_mb(f1(1)+indx+2)
              indx = indx + 3
            end do

            nframe = nframe + 1
            call zFixIon_sub_force_average(dbl_mb(r1(1)),
     >                                     dbl_mb(f1(1)),
     >                                     countz,
     >                                     int_mb(ion_numberz(1)), 
     >                                     dbl_mb(af1(1)),
     >                                     dbl_mb(af2(1)),
     >                                     dbl_mb(ac1(1)))
         end do

  30     continue
         close(unit=35)

         if (nframe.eq.0) then
            E = 0.0d0
         else
            E = 1.0d0/dble(nframe)
         end if
         call dscal(countz,E,dbl_mb(af1(1)),1)
         call dscal(countz,E,dbl_mb(af2(1)),1)
         call dscal(countz,E,dbl_mb(ac1(1)),1)

         write(*,110) nframe
         write(*,111) countz
         write(*,112) nion
  110    format(1x,'frames used           =',I8)
  111    format(1x,'number of constraints =',I8)
  112    format(1x,'number of atoms       =',I8)
         do m=1,countz
            write(*,120) int_mb(ion_numberz(1)+m-1),
     >                   dbl_mb(ac1(1)+m-1),
     >                   dbl_mb(af1(1)+m-1),
     >                   dbl_mb(af2(1)+m-1)-dbl_mb(af1(1)+m-1)**2
         end do
  120    format(1x,'atom=',I5,
     >          1x,'<z>=', E14.6,
     >          4x,'<Fz>=',E14.6,
     >          1x,'(<F**2>-<F>**2=',E10.3,')')

         ok =        BA_pop_stack(ac1(2))
         ok = ok.and.BA_pop_stack(af2(2))
         ok = ok.and.BA_pop_stack(af1(2))
         ok = ok.and.BA_pop_stack(f1(2))
         ok = ok.and.BA_pop_stack(r1(2))
         if (.not.ok)
     >    call errquit(
     >         'zFixIon_force_average: error popping stack',1,MA_ERR)

      end if

      return
      end

      subroutine zFixIon_sub_force_average(rion,fion,
     >                                     nindx,indx,
     >                                     af,af2,ac)
      implicit none
      real*8 rion(3,*)
      real*8 fion(3,*)
      integer nindx,indx(*)
      real*8 af(*)
      real*8 af2(*)
      real*8 ac(*)

*     **** local variables
      integer m
      real*8 ff

      do m=1,nindx
         ff = fion(3,indx(m)) 
         af(m)  = af(m)  + ff
         af2(m) = af2(m) + ff*ff
         ac(m)  = ac(m)  + rion(3,indx(m))
      end do

      return
      end



*     ***********************************
*     *				        *
*     *	    FixIon_force_average	*
*     *	 		                *
*     ***********************************
      subroutine FixIon_force_average(filename)
      implicit none
      character*(*) filename

#include "bafdecls.fh"
#include "errquit.fh"

*     **** FixIon common block ****
      integer count
      integer ion_number(2)
      common / FixIon_Block / ion_number,count


*     **** local variables ****
      logical ok
      integer nion
      integer m,indx,indx2,indx3,ii,iii,nframe,l
      double precision E
      double precision a1x,a1y,a1z
      double precision a2x,a2y,a2z
      double precision a3x,a3y,a3z
      integer r1(2),f1(2),af1(2),af2(2),ac1(2)
      character*2 tag2
      character*4 tag4

      if (count.gt.0) then
         l = index(filename,' ') - 1
         write(*,*)
         write(*,*) "Force Energy File:",filename(1:l)

         open(unit=35,file=filename,form='formatted',status='old')
         read(35,*) nion
         close(unit=35)

         ok =        BA_push_get(mt_dbl,3*nion,'r1',r1(2),r1(1))
         ok = ok.and.BA_push_get(mt_dbl,3*nion,'f1',f1(2),f1(1))
         ok = ok.and.BA_push_get(mt_dbl,3*count,'af1',af1(2),af1(1))
         ok = ok.and.BA_push_get(mt_dbl,3*count,'af2',af2(2),af2(1))
         ok = ok.and.BA_push_get(mt_dbl,3*count,'ac1',ac1(2),ac1(1))
         if (.not.ok)
     >    call errquit('zFixIon_force_average: out of stack',1,MA_ERR)
         call dcopy(3*count,0.0d0,0,dbl_mb(af1(1)),1)
         call dcopy(3*count,0.0d0,0,dbl_mb(af2(1)),1)
         call dcopy(3*count,0.0d0,0,dbl_mb(ac1(1)),1)

         open(unit=35,file=filename,form='formatted',status='old')
         nframe = 0
         do while (.true.)
            read(35,*,ERR=30,END=30) nion
            read(35,*,ERR=30,END=30) E
            read(35,*,ERR=30,END=30) a1x,a1y,a1z
            read(35,*,ERR=30,END=30) a2x,a2y,a2z
            read(35,*,ERR=30,END=30) a3x,a3y,a3z
            indx = 0
            do ii=1,nion
              read(35,*,ERR=30,END=30) iii,tag2,tag4,
     >                                 dbl_mb(r1(1)+indx),
     >                                 dbl_mb(r1(1)+indx+1),
     >                                 dbl_mb(r1(1)+indx+2),
     >                                 dbl_mb(f1(1)+indx),
     >                                 dbl_mb(f1(1)+indx+1),
     >                                 dbl_mb(f1(1)+indx+2)
              indx = indx + 3
            end do

            nframe = nframe + 1
            call FixIon_sub_force_average(dbl_mb(r1(1)),
     >                                     dbl_mb(f1(1)),
     >                                     count,
     >                                     int_mb(ion_number(1)), 
     >                                     dbl_mb(af1(1)),
     >                                     dbl_mb(af2(1)),
     >                                     dbl_mb(ac1(1)))
         end do

  30     continue
         close(unit=35)

         if (nframe.eq.0) then
            E = 0.0d0
         else
            E = 1.0d0/dble(nframe)
         end if
         call dscal(3*count,E,dbl_mb(af1(1)),1)
         call dscal(3*count,E,dbl_mb(af2(1)),1)
         call dscal(3*count,E,dbl_mb(ac1(1)),1)

         write(*,110) nframe
         write(*,111) count
         write(*,112) nion
  110    format(1x,'frames used           =',I8)
  111    format(1x,'number of constraints =',I8)
  112    format(1x,'number of atoms       =',I8)
         do m=0,count-1
            write(*,120) int_mb(ion_number(1)+m),
     >                   dbl_mb(ac1(1)+3*m),
     >                   dbl_mb(af1(1)+3*m),
     >                   dbl_mb(af2(1)+3*m)-dbl_mb(af1(1)+3*m)**2
            write(*,121) int_mb(ion_number(1)+m),
     >                   dbl_mb(ac1(1)+3*m+1),
     >                   dbl_mb(af1(1)+3*m+1),
     >                   dbl_mb(af2(1)+3*m+1)-dbl_mb(af1(1)+3*m+1)**2
            write(*,122) int_mb(ion_number(1)+m),
     >                   dbl_mb(ac1(1)+3*m+2),
     >                   dbl_mb(af1(1)+3*m+2),
     >                   dbl_mb(af2(1)+3*m+2)-dbl_mb(af1(1)+3*m+2)**2
         end do
  120    format(1x,'atom=',I5,
     >          1x,'<x>=', E14.6,
     >          4x,'<Fx>=',E14.6,
     >          1x,'(<Fx**2>-<Fx>**2=',E10.3,')')
  121    format(1x,'atom=',I5,
     >          1x,'<y>=', E14.6,
     >          4x,'<Fy>=',E14.6,
     >          1x,'(<Fy**2>-<Fy>**2=',E10.3,')')
  122    format(1x,'atom=',I5,
     >          1x,'<z>=', E14.6,
     >          4x,'<Fz>=',E14.6,
     >          1x,'(<Fz**2>-<Fz>**2=',E10.3,')')

         ok =        BA_pop_stack(ac1(2))
         ok = ok.and.BA_pop_stack(af2(2))
         ok = ok.and.BA_pop_stack(af1(2))
         ok = ok.and.BA_pop_stack(f1(2))
         ok = ok.and.BA_pop_stack(r1(2))
         if (.not.ok)
     >    call errquit(
     >         'FixIon_force_average: error popping stack',1,MA_ERR)

      end if

      return
      end

      subroutine FixIon_sub_force_average(rion,fion,
     >                                    nindx,indx,
     >                                    af,af2,ac)
      implicit none
      real*8 rion(3,*)
      real*8 fion(3,*)
      integer nindx,indx(*)
      real*8 af(*)
      real*8 af2(*)
      real*8 ac(*)

*     **** local variables
      integer m,mx,my,mz
      real*8 fx,fy,fz

      do m=1,nindx
         mx =3*(m-1)+1
         my =3*(m-1)+2
         mz =3*(m-1)+3
         fx = fion(1,indx(m)) 
         fy = fion(2,indx(m)) 
         fz = fion(3,indx(m)) 

         af(mx)  = af(mx)  + fx
         af(my)  = af(my)  + fy
         af(mz)  = af(mz)  + fz

         af2(mx) = af2(mx) + fx*fx
         af2(my) = af2(my) + fy*fy
         af2(mz) = af2(mz) + fz*fz

         ac(mx)  = ac(mx)  + rion(1,indx(m))
         ac(my)  = ac(my)  + rion(2,indx(m))
         ac(mz)  = ac(mz)  + rion(3,indx(m))
      end do

      return
      end





*     ***************************
*     *			        *
*     *	     FixIon_analysis    *
*     *				*
*     ***************************
      subroutine FixIon_analysis(rtdb)
      implicit none 
      integer   rtdb

#include "btdb.fh"
#include "errquit.fh"
#include "stdio.fh"

*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)

      logical value
      character*50  filename
      character*255 full_filename

*     **** external functions ****

      call ga_sync()

      call Parallel_taskid(taskid)


      call ion_init_FixIon(rtdb)
      value = btdb_parallel(.false.)
      if (taskid.eq.MASTER) then



         if (.not.btdb_cget(rtdb,'cpmd:fei_filename',1,filename))
     >     call util_file_prefix('fei',filename)

         call util_file_name_noprefix(filename,.false.,
     >                             .false.,
     >                             full_filename)

         call util_file_name_noprefix(filename,.false.,
     >                               .false.,
     >                               full_filename)

         write(luout,*)
         write(luout,*) '===== Performing FixIon analysis ===='
         write(luout,*)
         call FixIon_force_average(full_filename)

         write(luout,*)
         write(luout,*) '===== Performing z FixIon analysis ===='
         write(luout,*)
         call zFixIon_force_average(full_filename)
         


      end if
      call ga_sync()
      value = btdb_parallel(.true.)
      call ion_end_FixIon()
   

      return
      end



c $Id$
