*     ***********************************
*     *                 				*
*     *           Balance_Init 		    *
*     *                 				*
*     ***********************************

      subroutine Balance_Init(nida,nidb,nidb_out)
      implicit none
      integer nida(0:1),nidb(0:1)
      integer nidb_out(0:1)

#include "mafdecls.fh"


*     **** Balance Common Block ****
      integer packet_size(2,0:1),indx_start(2,0:1)
      integer proc_to(2,0:1),proc_from(2,0:1)
      integer npacket(0:1)
      logical receiver(0:1),sender(0:1)
      common /Balance_Block/ packet_size,indx_start,
     >                       proc_to,proc_from,npacket,receiver,sender

*     **** local variables ****
      logical value
      integer nb,np,taskid
      integer nwave,nwave_out

      call Parallel_np(np)
      call Parallel_taskid(taskid)

*     **** allocate balance memory ****
      value =  MA_alloc_get(mt_int,np,
     >                   'psizea',packet_size(2,0),packet_size(1,0))
      value = value.and.
     >       MA_alloc_get(mt_int,np,
     >                   'psizeb',packet_size(2,1),packet_size(1,1))
      value = value.and.
     >       MA_alloc_get(mt_int,np,
     >                   'indxsa',indx_start(2,0),indx_start(1,0))
      value = value.and.
     >       MA_alloc_get(mt_int,np,
     >                   'indxsb',indx_start(2,1),indx_start(1,1))
      value = value.and.
     >       MA_alloc_get(mt_int,np,
     >                   'prctoa',proc_to(2,0),proc_to(1,0))
      value = value.and.
     >       MA_alloc_get(mt_int,np,
     >                   'prctob',proc_to(2,1),proc_to(1,1))
      value = value.and.
     >       MA_alloc_get(mt_int,np,
     >                   'prcfra',proc_from(2,0),proc_from(1,0))
      value = value.and.
     >       MA_alloc_get(mt_int,np,
     >                   'prcfrb',proc_from(2,1),proc_from(1,1))
      if (.not. value) 
     >   call errquit('Balance_init: out of heap memory',0)



      do nb=0,1
         nwave = nida(nb)+nidb(nb)
         call Balance_Init_a(nwave,np,taskid,nwave_out,
     >      npacket(nb),
     >      receiver(nb),sender(nb),
     >      int_mb(proc_to(1,nb)),     int_mb(proc_from(1,nb)),
     >      int_mb(packet_size(1,nb)), int_mb(indx_start(1,nb)))

         nidb_out(nb) = nidb(nb) + (nwave_out-nwave)

      end do

      return
      end



*     ***********************************
*     *                 				*
*     *           Balance_End 		    *
*     *                 				*
*     ***********************************

      subroutine Balance_End()
      implicit none

#include "mafdecls.fh"


*     **** Balance Common Block ****
      integer packet_size(2,0:1),indx_start(2,0:1)
      integer proc_to(2,0:1),proc_from(2,0:1)
      integer npacket(0:1)
      logical receiver(0:1),sender(0:1)
      common /Balance_Block/ packet_size,indx_start,
     >                       proc_to,proc_from,npacket,receiver,sender

      logical value

      value =           MA_free_heap(packet_size(2,0))
      value = value.and.MA_free_heap(packet_size(2,1))
      value = value.and.MA_free_heap(indx_start(2,0))
      value = value.and.MA_free_heap(indx_start(2,1))
      value = value.and.MA_free_heap(proc_to(2,0))
      value = value.and.MA_free_heap(proc_to(2,1))
      value = value.and.MA_free_heap(proc_from(2,0))
      value = value.and.MA_free_heap(proc_from(2,1))
      if (.not. value) call errquit('error freeing heap memory',0)

      return
      end

*     ***********************************
*     *                 				*
*     *           Balance_Init_a	    *
*     *                 				*
*     ***********************************
*    This routine defines the balance data structure

      subroutine Balance_Init_a(nwave,np,taskid,
     >                          nwave_out,
     >                          npacket,receiver,sender,
     >                          proc_to,proc_from,
     >                          packet_size,indx_start)
      implicit none
      integer nwave,np,taskid
      integer nwave_out

      integer npacket
      logical receiver,sender
      integer proc_to(*),proc_from(*)
      integer packet_size(*)
      integer indx_start(*)

#include "mafdecls.fh"

*     ***** local variables ****
      logical done,value
      integer i,j
      integer ave,short,long
      integer above,below

c      integer nwave2(0:(np-1))
c      integer indx(0:(np-1))
      integer nwave2(2),indx(2)     

*     **** allocate nwave2 and indx off the stack ****
      value = MA_push_get(mt_int,(np),
     >                     'nwave2',nwave2(2),nwave2(1))
      value = value.and.
     >        MA_push_get(mt_int,(np),
     >                     'indx',indx(2),indx(1))
      if (.not. value) call errquit('out of stack memory',0)
   
*     **** define nwave2 ****
      do i=0,np-1
c        nwave2(i) = 0
         int_mb(nwave2(1)+i) = 0
      end do
c     nwave2(taskid) = nwave
      int_mb(nwave2(1)+taskid) = nwave
c     call D3dB_Vector_ISumAll(np,nwave2)
      call D3dB_Vector_ISumAll(np,int_mb(nwave2(1)))

*     **** get the sorting index ****
c     call nwave2_sort(np,nwave2,indx)
      call nwave2_sort(np,int_mb(nwave2(1)),int_mb(indx(1)))

*     ***** get the average ****
      ave = 0
      do i=0,np-1
c        ave = ave + nwave2(i)
         ave = ave + int_mb(nwave2(1)+i)
      end do
      ave = ave/np

*     ***** get below ***
      below = -1
      do while (int_mb(nwave2(1) + int_mb(indx(1)+below+1)).lt.ave) 
        below = below + 1
      end do
         
*     ***** get above ***
      above = np
      do while (int_mb(nwave2(1) + int_mb(indx(1)+above-1)).gt.ave) 
        above = above - 1
      end do


      npacket  = 0
      receiver = .false.
      sender   = .false.

      if (np.gt.1) then
        i = 0 
        j = np-1
        done = .false.
        if (i .gt. below) done = .true.
        if (j .lt. above) done = .true.
        do while (.not. done)
           short = ave - int_mb(nwave2(1)+int_mb(indx(1)+i))
           long =  int_mb(nwave2(1)+int_mb(indx(1)+j)) - ave

           if (taskid.eq.int_mb(indx(1)+i)) then
              npacket = npacket + 1
              proc_from(npacket) = int_mb(indx(1)+j)
              receiver = .true.
           end if

           if (taskid.eq.int_mb(indx(1)+j)) then
              npacket = npacket + 1
              proc_to(npacket) = int_mb(indx(1)+i)
              sender   = .true.
           end if


           if (short.eq.long) then

             if (taskid.eq.int_mb(indx(1)+i)) then
                packet_size(npacket) = short
                indx_start(npacket)  = 
     >              int_mb(nwave2(1)+int_mb(indx(1)+i)) + 1
             end if

             if (taskid.eq.int_mb(indx(1)+j)) then
                packet_size(npacket) = long
                indx_start(npacket) = 
     >            int_mb(nwave2(1)+int_mb(indx(1)+j)) - long + 1
             end if

             int_mb(nwave2(1)+int_mb(indx(1)+i)) = 
     >         int_mb(nwave2(1)+int_mb(indx(1)+i)) + short
             int_mb(nwave2(1)+int_mb(indx(1)+j)) = 
     >         int_mb(nwave2(1)+int_mb(indx(1)+j)) - long
             i = i + 1
             j = j - 1


           else if (short.lt.long) then

             if (taskid.eq.int_mb(indx(1)+i)) then
               packet_size(npacket) = short
               indx_start(npacket) = 
     >            int_mb(nwave2(1)+int_mb(indx(1)+i)) + 1
             end if

             if (taskid.eq.int_mb(indx(1)+j)) then
               packet_size(npacket) = short
               indx_start(npacket) = 
     >            int_mb(nwave2(1)+int_mb(indx(1)+j)) - short + 1
             end if

             int_mb(nwave2(1)+int_mb(indx(1)+i)) = 
     >         int_mb(nwave2(1)+int_mb(indx(1)+i)) + short
             int_mb(nwave2(1)+int_mb(indx(1)+j)) = 
     >         int_mb(nwave2(1)+int_mb(indx(1)+j)) - short
             i = i + 1


           else if (short.gt.long) then
             if (taskid.eq.int_mb(indx(1)+i)) then
               packet_size(npacket) = long
               indx_start(npacket) = 
     >           int_mb(nwave2(1)+int_mb(indx(1)+i)) + 1
             end if

             if (taskid.eq.int_mb(indx(1)+j)) then
               packet_size(npacket) = long
               indx_start(npacket) = 
     >           int_mb(nwave2(1)+int_mb(indx(1)+j)) - long + 1
             end if

             int_mb(nwave2(1)+int_mb(indx(1)+i)) = 
     >          int_mb(nwave2(1)+int_mb(indx(1)+i)) + long
             int_mb(nwave2(1)+int_mb(indx(1)+j)) = 
     >          int_mb(nwave2(1)+int_mb(indx(1)+j)) - long
             j = j - 1

           end if

           if (i .gt. below) done = .true.
           if (j .lt. above) done = .true.

        end do

      end if

      nwave_out = int_mb(nwave2(1)+taskid)

      value =           MA_pop_stack(indx(2))
      value = value.and.MA_pop_stack(nwave2(2))
      if (.not. value) call errquit('error freeing stack memory',0)
  

      return
      end

      subroutine nwave2_sort(n,f,indx)
      integer n
      integer f(0:(n-1))
      integer indx(0:(n-1))

      integer i,j,idum
      do i=0,n-1
        indx(i) = i
      end do
      do i=0,(n-2)
      do j=i+1,(n-1)
        if (f(indx(j)).lt.f(indx(i))) then
              idum    = indx(i)
              indx(i) = indx(j)
              indx(j) = idum
           end if
      end do
      end do

      return
      end 


      subroutine Pack_c_balance(nb,A)
      implicit none
      integer nb
      complex*16 A(*)

#include "mafdecls.fh"


#include "tcgmsg.fh"
#include "msgtypesf.h"
      integer  rcv_len,rcv_proc

*     **** Balance Common Block ****
      integer packet_size(2,0:1),indx_start(2,0:1)
      integer proc_to(2,0:1),proc_from(2,0:1)
      integer npacket(0:1)
      logical receiver(0:1),sender(0:1)
      common /Balance_Block/ packet_size,indx_start,
     >                       proc_to,proc_from,npacket,receiver,sender

      integer j
      integer pto,pfrom,msglen,indx

      if (sender(nb)) then
         do j=1,npacket(nb)
            pto  = int_mb(proc_to(1,nb)+j-1)    
            msglen = int_mb(packet_size(1,nb)+j-1)
            indx = int_mb(indx_start(1,nb)+j-1)
c            send data....
            if (msglen.gt.0) then
               call SND(9+MSGDBL,
     >                  A(indx),
     >                  mdtob(2*msglen),pto,1)
            end if
  

         end do
      end if

      if (receiver(nb)) then
         do j=1,npacket(nb)
            pfrom = int_mb(proc_from(1,nb)+j-1)    
            msglen = int_mb(packet_size(1,nb)+j-1)
            indx = int_mb(indx_start(1,nb)+j-1)
c            recieve data....
            if (msglen.gt.0) then
               call RCV(9+MSGDBL,
     >                  A(indx),
     >                  mdtob(2*msglen),rcv_len,
     >                  pfrom,rcv_proc,1)
            end if

         end do
      end if

      return
      end




      subroutine Pack_t_balance(nb,A)
      implicit none
      integer nb
      real*8 A(*)

#include "mafdecls.fh"

#include "tcgmsg.fh"
#include "msgtypesf.h"
      integer  rcv_len,rcv_proc


*     **** Balance Common Block ****
      integer packet_size(2,0:1),indx_start(2,0:1)
      integer proc_to(2,0:1),proc_from(2,0:1)
      integer npacket(0:1)
      logical receiver(0:1),sender(0:1)
      common /Balance_Block/ packet_size,indx_start,
     >                       proc_to,proc_from,npacket,receiver,sender

      integer j
      integer pto,pfrom,msglen,indx


      if (sender(nb)) then
         do j=1,npacket(nb)
            pto  = int_mb(proc_to(1,nb)+j-1)    
            msglen = int_mb(packet_size(1,nb)+j-1)
            indx = int_mb(indx_start(1,nb)+j-1)
c            send data....
            if (msglen.gt.0) then
               call SND(9+MSGDBL,
     >                  A(indx),
     >                  mdtob(msglen),pto,1)
            end if
   

         end do
      end if

      if (receiver(nb)) then
         do j=1,npacket(nb)
            pfrom = int_mb(proc_from(1,nb)+j-1)    
            msglen = int_mb(packet_size(1,nb)+j-1)
            indx = int_mb(indx_start(1,nb)+j-1)
c            recieve data....
            if (msglen.gt.0) then
               call RCV(9+MSGDBL,
     >                  A(indx),
     >                  mdtob(msglen),rcv_len,
     >                  pfrom,rcv_proc,1)
            end if


         end do
      end if

      return
      end




      subroutine Pack_c_unbalance(nb,A)
      implicit none
      integer nb
      complex*16 A(*)

#include "mafdecls.fh"


#include "tcgmsg.fh"
#include "msgtypesf.h"
      integer  rcv_len,rcv_proc


*     **** Balance Common Block ****
      integer packet_size(2,0:1),indx_start(2,0:1)
      integer proc_to(2,0:1),proc_from(2,0:1)
      integer npacket(0:1)
      logical receiver(0:1),sender(0:1)
      common /Balance_Block/ packet_size,indx_start,
     >                       proc_to,proc_from,npacket,receiver,sender

      integer j
      integer pto,pfrom,msglen,indx


      if (sender(nb)) then
         do j=1,npacket(nb)
            pfrom  = int_mb(proc_to(1,nb)+j-1)    
            msglen = int_mb(packet_size(1,nb)+j-1)
            indx = int_mb(indx_start(1,nb)+j-1)
c            recieve data....
            if (msglen.gt.0) then
               call RCV(9+MSGDBL,
     >                  A(indx),
     >                  mdtob(2*msglen),rcv_len,
     >                  pfrom,rcv_proc,1)
            end if
  
         end do
      end if

      if (receiver(nb)) then
         do j=1,npacket(nb)
            pto    = int_mb(proc_from(1,nb)+j-1)    
            msglen = int_mb(packet_size(1,nb)+j-1)
            indx   = int_mb(indx_start(1,nb)+j-1)
c            send data....
            if (msglen.gt.0) then
               call SND(9+MSGDBL,
     >                  A(indx),
     >                  mdtob(2*msglen),pto,1)
            end if

         end do
      end if

      return
      end

