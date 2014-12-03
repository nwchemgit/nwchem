c
c     $Id$
c
*     *********************************
*     *                               *
*     *          pair_distab          *
*     *                               *
*     *********************************
      subroutine pair_distab(mfilename)
      implicit none
      character*(*) mfilename

#include "errquit.fh"
#include "bafdecls.fh"
#include "stdio.fh"

c     *** parameters ***
      integer MASTER,taskid,np
      parameter (MASTER=0)

      real*8 DELR
      parameter (DELR=0.01d0)
 

*     **** local variables ****
      logical value,found
      integer MAXBIN
      integer i,j,bin,natoms,natoms0,t,nframes,iii,icount
      integer nindxA,nindxB,nkatm,iA,iB
      integer rx(2),ry(2),rz(2),hst(2),symb(2),indxA(2),indxB(2)
      real*8  time1,time2,volume,volume1
      real*8  x,y,z,r,gr,v1,v2,v3
      real*8 RLOWER,RUPPER,NIDEAL,CONST    
      real*8 grs,grm2,grm1,grp1,grp2,sum1,sum2
      integer k1,k2,n1,n2,n3
      real*8 unita(3,3)
      character*2 tag2
      character*2 tag4
      character*2 tagA,tagB
      character*255 full_filename,gfilename

*     ***************************************
*     ****** PAIR DISTRIBUTION routine ******
*     ***************************************

      call Parallel_np(np)
      call Parallel_taskid(taskid)

      if (taskid.eq.MASTER) then
         write(luout,*)
         write(luout,*) 
     =   '===== Calculating radial distribution functions ===='
         write(luout,*)

*        *** get the number of frames from MOTION ***
         call MOTION_nframes(mfilename,nframes,natoms,volume,unita)
      end if
      call Parallel_Brdcst_ivalues(MASTER,1,nframes)
      call Parallel_Brdcst_ivalues(MASTER,1,natoms)
      call Parallel_Brdcst_values(MASTER,1,volume)
      call Parallel_Brdcst_values(MASTER,9,unita)
      

      MAXBIN = INT(0.5d0*(volume**(0.33333333d0))/DELR)

      value = BA_push_get(mt_int,MAXBIN,'hst89',hst(2),hst(1))
      value = value.and.
     >        BA_push_get(mt_int,natoms,'indxA',indxA(2),indxA(1))
      value = value.and.
     >        BA_push_get(mt_int,natoms,'indxB',indxB(2),indxB(1))
      value = value.and.
     >        BA_push_get(mt_dbl,natoms,'rx',rx(2),rx(1))
      value = value.and.
     >        BA_push_get(mt_dbl,natoms,'ry',ry(2),ry(1))
      value = value.and.
     >        BA_push_get(mt_dbl,natoms,'rz',rz(2),rz(1))
      value = value.and.
     >        BA_push_get(mt_byte,2*natoms,'symb',symb(2),symb(1))
      if (.not. value) 
     >   call errquit('pair_distab: out of stack memory',0,MA_ERR)
      
      if (taskid.eq.MASTER) then
         iA = 1
         iB = 3
         gfilename
     >   ='                                                    '
         nkatm = 0
         open(unit=19,file=mfilename,form='formatted')
         read(19,*) time1,natoms0,volume1,unita
         do i=1,natoms0
            read(19,*) iii,tag2,tag4,x,y,z,v1,v2,v3

            found = .false.
            j = 0
            do while ((j.lt.nkatm).and.(.not.found))
               tagA(1:1) = byte_mb(symb(1)+2*j)
               tagA(2:2) = byte_mb(symb(1)+2*j+1)
               if (tag2.eq.tagA) found = .true.
               j = j + 1
            end do

            if (.not.found) then
               byte_mb(symb(1)+2*nkatm)   = tag2(1:1)
               byte_mb(symb(1)+2*nkatm+1) = tag2(2:2)
               nkatm = nkatm+1
               gfilename(iA:iB) = tag2//' '
               iA = iA + 3
               iB = iB + 3
            end if

         end do
         close(19)

         write(luout,*) 
         write(luout,*) 'frames used: ',nframes
         write(luout,*) 'natoms:      ',natoms
         write(luout,*) 'nkatoms:     ',nkatm
         write(luout,*) 'atoms:       ',gfilename
         write(luout,*) 'delr:        ',DELR
         write(luout,*) 'volume:      ',volume
         write(luout,1100) '    a1:      ',
     >                     unita(1,1),unita(2,1),unita(3,1)
         write(luout,1100) '    a2:      ',
     >                     unita(1,2),unita(2,2),unita(3,2)
         write(luout,1100) '    a3:      ',
     >                     unita(1,3),unita(2,3),unita(3,3)
         write(luout,*) 'maxbin:      ',MAXBIN
         iii = index(mfilename,' ')
         write(luout,*) 'MOTION filename:    ',mfilename(1:iii)
         write(luout,*) 
         call util_flush(luout)
      end if
 1100 FORMAT(A,3F8.3)

      call Parallel_Brdcst_ivalues(MASTER,1,nkatm)
      do iA = 1,nkatm
      do iB = 1,nkatm
         do bin=0,MAXBIN-1
            int_mb(hst(1)+bin) = 0
         end do

         if (taskid.eq.MASTER) then
            tagA(1:1) = byte_mb(symb(1)+2*iA-2)
            tagA(2:2) = byte_mb(symb(1)+2*iA-1)
            tagB(1:1) = byte_mb(symb(1)+2*iB-2)
            tagB(2:2) = byte_mb(symb(1)+2*iB-1)
            open(unit=19,file=mfilename,form='formatted')
         end if
     
*        **** calculate the pair distribution function ****
         time1 = 0.0d0
         do t=1,nframes
        
            if (taskid.eq.MASTER) then
               time2    = time1
               read(19,*) time1,natoms0,volume1,unita
               nindxA = 0
               nindxB = 0
               do i=1,natoms0
                  read(19,*) iii,tag2,tag4,x,y,z,v1,v2,v3
                  dbl_mb(rx(1)+(i-1)) = x
                  dbl_mb(ry(1)+(i-1)) = y
                  dbl_mb(rz(1)+(i-1)) = z          
                  if (tag2.eq.tagA) then
                     int_mb(indxA(1)+nindxA) = i
                     nindxA = nindxA + 1
                  end if
                  if (tag2.eq.tagB) then
                     int_mb(indxB(1)+nindxB) = i
                     nindxB = nindxB + 1
                  end if
               end do
            end if
            call Parallel_Brdcst_ivalues(MASTER,1,natoms0)
            call Parallel_Brdcst_ivalues(MASTER,1,nindxA)
            call Parallel_Brdcst_ivalues(MASTER,1,nindxB)
            call Parallel_Brdcst_ivalues(MASTER,nindxA,int_mb(indxA(1)))
            call Parallel_Brdcst_ivalues(MASTER,nindxB,int_mb(indxB(1)))
            call Parallel_Brdcst_values(MASTER,natoms0,dbl_mb(rx(1)))
            call Parallel_Brdcst_values(MASTER,natoms0,dbl_mb(ry(1)))
            call Parallel_Brdcst_values(MASTER,natoms0,dbl_mb(rz(1)))

            icount = 0
            do k1=1,nindxA
              i = int_mb(indxA(1)+k1-1)
              do k2=1,nindxB
                j = int_mb(indxB(1)+k2-1)
                if (icount.eq.taskid) then
                   do n1=-1,1
                   do n2=-1,1
                   do n3=-1,1
                     x = dbl_mb(rx(1)+(i-1)) - dbl_mb(rx(1)+(j-1))
     >                 + n1*unita(1,1) + n2*unita(1,2) + n3*unita(1,3)
                     y = dbl_mb(ry(1)+(i-1)) - dbl_mb(ry(1)+(j-1))
     >                 + n1*unita(2,1) + n2*unita(2,2) + n3*unita(2,3)
                     z = dbl_mb(rz(1)+(i-1)) - dbl_mb(rz(1)+(j-1))
     >                 + n1*unita(3,1) + n2*unita(3,2) + n3*unita(3,3)
                     r = dsqrt(x*x+y*y+z*z)
                     if (r.gt.DELR) then
                       bin = INT(r/DELR) 
                       if (bin.le.(MAXBIN-1)) 
     >                  int_mb(hst(1)+bin) = int_mb(hst(1)+bin) + 1
                     end if
                   end do
                   end do
                   end do
                end if
                icount = mod(icount+1,np)
              end do
            end do

         end do
         if (taskid.eq.MASTER) close(unit=19)

         call Parallel_Vector_ISumAll(MAXBIN,int_mb(hst(1)))

         if (taskid.eq.MASTER) then
*        ***** write out g(r) in PAIR_DISTRIBUTION ****
         i = 2
         j = 2
         if (tagA(2:2).eq.' ') i = 1
         if (tagB(2:2).eq.' ') j = 1
         gfilename = tagA(1:i)//tagB(1:j)//'.gr'
         call util_file_name_noprefix(gfilename,.false.,
     >                             .false.,
     >                             full_filename)
         open(unit=19,file=full_filename,form='formatted')
         iii = index(full_filename,' ')
         write(luout,*) 'creating gr filename:',full_filename(1:iii)

         gfilename = tagA(1:i)//tagB(1:j)//'.hist'
         call util_file_name_noprefix(gfilename,.false.,
     >                             .false.,
     >                             full_filename)
         open(unit=21,file=full_filename,form='formatted')
         iii = index(full_filename,' ')
         write(luout,*) 'creating hist filename:',full_filename(1:iii)
         write(luout,*) 
         call util_flush(luout)

***      **** normalize hist(bin) to give g(r) ****
         sum1=0.0d0
         sum2=0.0d0
         CONST =       (16.0d0*datan(1.0d0)/3.0d0)*(nindxB/volume)
         do bin=2,MAXBIN-3
           RLOWER = dble(bin)*DELR
           RUPPER = RLOWER + DELR
           NIDEAL = CONST*(RUPPER**3 - RLOWER**3)
           !NIDEAL = NIDEAL*nindxA*nindxB
           NIDEAL = NIDEAL*nindxA
           r  = 0.5d0*(RLOWER+RUPPER)
           gr = dble(int_mb(hst(1)+bin))/dble(nframes)

           RLOWER = dble(bin-2)*DELR
           RUPPER = RLOWER + DELR
           !NIDEAL = CONST*(RUPPER**3 - RLOWER**3)
           grm2 = dble(int_mb(hst(1)+bin-2))/dble(nframes)

           RLOWER = dble(bin-1)*DELR
           RUPPER = RLOWER + DELR
           !NIDEAL = CONST*(RUPPER**3 - RLOWER**3)
           grm1 = dble(int_mb(hst(1)+bin-1))/dble(nframes)

           RLOWER = dble(bin+1)*DELR
           RUPPER = RLOWER + DELR
           !NIDEAL = CONST*(RUPPER**3 - RLOWER**3)
           grp1 = dble(int_mb(hst(1)+bin+1))/dble(nframes)

           RLOWER = dble(bin+2)*DELR
           RUPPER = RLOWER + DELR
           !NIDEAL = CONST*(RUPPER**3 - RLOWER**3)
           grp2 = dble(int_mb(hst(1)+bin+2))/dble(nframes)

           grs = ( -3.0d0*grm2
     >           + 12.0d0*grm1
     >           + 17.0d0*gr
     >           + 12.0d0*grp1
     >           -  3.0d0*grp2)/35.0d0

           sum1 = sum1 + grs
           sum2 = sum2 + gr
           write(19,1200) 0.529177d0*r,grs/NIDEAL,sum1,gr/NIDEAL,sum2
           write(21,1201) 0.529177d0*r,grs,gr
        end do
        close(unit=19)
        close(unit=21)
      end if

      end do
      end do
 1200 FORMAT(5E14.5)
 1201 FORMAT(3E14.5)

      value =           BA_pop_stack(symb(2))
      value = value.and.BA_pop_stack(rz(2))
      value = value.and.BA_pop_stack(ry(2))
      value = value.and.BA_pop_stack(rx(2))
      value = value.and.BA_pop_stack(indxB(2))
      value = value.and.BA_pop_stack(indxA(2))
      value = value.and.BA_pop_stack(hst(2))
      if (.not. value) 
     >  call errquit('pair_distab:cannot release stack memory',0,
     &       MA_ERR)


*     *** synchronization ***
      call ga_sync()

      return
      end

