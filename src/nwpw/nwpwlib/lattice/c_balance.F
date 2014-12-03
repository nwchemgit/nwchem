*
* $Id$
*

*     ***********************************
*     *               			*
*     *         c_Balance_Init 	        *
*     *                 		*
*     ***********************************

      subroutine c_Balance_Init(maxsize0,nidb,nidb_out)
      implicit none
      integer maxsize0
      integer nidb(0:maxsize0-1)
      integer nidb_out(0:maxsize0-1)

#include "bafdecls.fh"
#include "errquit.fh"
#include "balance_common.fh"

*     **** local variables ****
      logical value
      integer nb,np,taskid
      integer nwave,nwave_out
      integer dum(2)


      maxsize = maxsize0
      call Parallel3d_np_i(np)
      call Parallel3d_taskid_i(taskid)

*     **** allocate balance memory ****
      value =  BA_alloc_get(mt_int,2*maxsize,
     >                      'psizea_list',
     >                      packet_size_list(2),
     >                      packet_size_list(1))
      value = value.and.
     >       BA_alloc_get(mt_int,2*maxsize,
     >                    'indxsa_list',
     >                    indx_start_list(2),
     >                    indx_start_list(1))
      value = value.and.
     >       BA_alloc_get(mt_int,2*maxsize,
     >                    'prctoa_list',
     >                    proc_to_list(2),
     >                    proc_to_list(1))
      value = value.and.
     >       BA_alloc_get(mt_int,2*maxsize,
     >                   'prcfra_list',
     >                    proc_from_list(2),
     >                    proc_from_list(1))

      value = value.and.
     >       BA_alloc_get(mt_int,maxsize,
     >                   'npacket_list',
     >                    npacket_list(2),
     >                    npacket_list(1))
      value = value.and.
     >       BA_alloc_get(mt_log,maxsize,
     >                   'receiver_list',
     >                    receiver_list(2),
     >                    receiver_list(1))
      value = value.and.
     >       BA_alloc_get(mt_log,maxsize,
     >                   'sender_list',
     >                    sender_list(2),
     >                    sender_list(1))

      do nb=0,maxsize-1

*        **** allocate balance memory ****
         value =  value.and.
     >          BA_alloc_get(mt_int,np,
     >                   'psizea',dum(2),dum(1))
         int_mb(packet_size_list(1)+2*nb  ) = dum(1)
         int_mb(packet_size_list(1)+2*nb+1) = dum(2)

         value = value.and.
     >       BA_alloc_get(mt_int,np,
     >                   'indxsa',dum(2),dum(1))
         int_mb(indx_start_list(1)+2*nb  ) = dum(1)
         int_mb(indx_start_list(1)+2*nb+1) = dum(2)

         value = value.and.
     >       BA_alloc_get(mt_int,np,
     >                   'prctoa',dum(2),dum(1))
         int_mb(proc_to_list(1)+2*nb  ) = dum(1)
         int_mb(proc_to_list(1)+2*nb+1) = dum(2)

         value = value.and.
     >       BA_alloc_get(mt_int,np,
     >                   'prcfra',dum(2),dum(1))
         int_mb(proc_from_list(1)+2*nb  ) = dum(1)
         int_mb(proc_from_list(1)+2*nb+1) = dum(2)

      end do

      if (.not. value) 
     >   call errquit('c_Balance_init: out of heap memory',0, MA_ERR)



      do nb=0,maxsize-1
         nwave = nidb(nb)
         call c_Balance_Init_a(nwave,np,taskid,nwave_out,
     >      int_mb(npacket_list(1) +nb),
     >      log_mb(receiver_list(1)+nb),
     >      log_mb(sender_list(1)  +nb),
     >      int_mb(int_mb(proc_to_list(1)    +2*nb)),     
     >      int_mb(int_mb(proc_from_list(1)  +2*nb)),
     >      int_mb(int_mb(packet_size_list(1)+2*nb)), 
     >      int_mb(int_mb(indx_start_list(1) +2*nb)))

         nidb_out(nb) = nidb(nb) + (nwave_out-nwave)
      end do

      return
      end 


*     ***********************************
*     *                 		*
*     *          c_Balance_End 		*
*     *                 		*
*     ***********************************

      subroutine c_Balance_End()
      implicit none

#include "errquit.fh"
#include "bafdecls.fh"
#include "balance_common.fh"


*     **** local variables ****
      logical value
      integer nb,dum2

      value = .true.
      do nb=0,maxsize-1

         dum2 = int_mb(packet_size_list(1)+2*nb+1)
         value = value.and.BA_free_heap(dum2)

         dum2 = int_mb(indx_start_list(1)+2*nb+1) 
         value = value.and.BA_free_heap(dum2)

         dum2 = int_mb(proc_to_list(1)+2*nb+1)
         value = value.and.BA_free_heap(dum2)

         dum2 = int_mb(proc_from_list(1)+2*nb+1)
         value = value.and.BA_free_heap(dum2)

      end do

      value = value.and.BA_free_heap(packet_size_list(2))
      value = value.and.BA_free_heap(indx_start_list(2))
      value = value.and.BA_free_heap(proc_to_list(2))
      value = value.and.BA_free_heap(proc_from_list(2))

      value = value.and.BA_free_heap(npacket_list(2))
      value = value.and.BA_free_heap(receiver_list(2))
      value = value.and.BA_free_heap(sender_list(2))
      if (.not. value) 
     > call errquit('c_Balance_end: error freeing heap memory',0,MA_ERR)

      return
      end

*     ***********************************
*     *                 		*
*     *           c_Balance_Init_a	*
*     *                 		*
*     ***********************************
*    This routine defines the balance data structure

      subroutine c_Balance_Init_a(nwave,np,taskid,
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

#include "bafdecls.fh"
#include "errquit.fh"

*     ***** local variables ****
      logical done,value
      integer i,j
      integer ave,short,long
      integer above,below

c      integer nwave2(0:(np-1))
c      integer indx(0:(np-1))
      integer nwave2(2),indx(2)     

*     **** allocate nwave2 and indx off the stack ****
      value = BA_push_get(mt_int,(np),
     >                     'nwave2',nwave2(2),nwave2(1))
      value = value.and.
     >        BA_push_get(mt_int,(np),
     >                     'indx',indx(2),indx(1))
      if (.not. value) 
     >  call errquit('c_Balance_init_a:out of stack memory',0,MA_ERR)
   
*     **** define nwave2 ****
      do i=0,np-1
         int_mb(nwave2(1)+i) = 0
      end do
      int_mb(nwave2(1)+taskid) = nwave
      call C3dB_Vector_ISumAll(np,int_mb(nwave2(1)))
     
*     **** get the sorting index ****
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

      value =           BA_pop_stack(indx(2))
      value = value.and.BA_pop_stack(nwave2(2))
      if (.not. value) 
     > call errquit('c_Balance_init_a:error freeing stack memory',0,
     >       MA_ERR)
  

      return
      end


*     ************************************
*     *                                  *
*     *         c_Balance_c_balance      *
*     *                                  *
*     ************************************

      subroutine c_Balance_c_balance(nb,A)
      implicit none
      integer nb
      complex*16 A(*)

#include "bafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "balance_common.fh"

*     **** local variables ****
      integer  rcv_len,rcv_proc
      integer j
      integer pto,pfrom,msglen,indx

*     **** external functions ****
      integer  Parallel3d_convert_taskid_i
      external Parallel3d_convert_taskid_i


      if (log_mb(sender_list(1)+nb)) then
         do j=1,int_mb(npacket_list(1)+nb)
            pto    = int_mb(int_mb(proc_to_list(1)    +2*nb)+j-1)    
            msglen = int_mb(int_mb(packet_size_list(1)+2*nb)+j-1)
            indx   = int_mb(int_mb(indx_start_list(1) +2*nb)+j-1)
c            send data....
            if (msglen.gt.0) then
               call SND(9+MSGDBL,
     >                  A(indx),
     >                  mdtob(2*msglen),
     >                  Parallel3d_convert_taskid_i(pto),
     >                  1)
            end if
  

         end do
      end if

      if (log_mb(receiver_list(1)+nb)) then
         do j=1,int_mb(npacket_list(1)+nb)
            pfrom  = int_mb(int_mb(proc_from_list(1)  +2*nb)+j-1)    
            msglen = int_mb(int_mb(packet_size_list(1)+2*nb)+j-1)
            indx   = int_mb(int_mb(indx_start_list(1) +2*nb)+j-1)
c            recieve data....
            if (msglen.gt.0) then
               call RCV(9+MSGDBL,
     >                  A(indx),
     >                  mdtob(2*msglen),rcv_len,
     >                  Parallel3d_convert_taskid_i(pfrom),
     >                  rcv_proc,1)
            end if

         end do
      end if

      return
      end


*     ************************************
*     *                                  *
*     *         c_Balance_t_balance      *
*     *                                  *
*     ************************************

      subroutine c_Balance_t_balance(nb,A)
      implicit none
      integer nb
      real*8 A(*)

#include "bafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "balance_common.fh"

*     **** local variables ****
      integer  rcv_len,rcv_proc
      integer j
      integer pto,pfrom,msglen,indx

*     **** external functions ****
      integer  Parallel3d_convert_taskid_i
      external Parallel3d_convert_taskid_i



      if (log_mb(sender_list(1)+nb)) then
         do j=1,int_mb(npacket_list(1)+nb)
            pto    = int_mb(int_mb(proc_to_list(1)    +2*nb)+j-1)    
            msglen = int_mb(int_mb(packet_size_list(1)+2*nb)+j-1)
            indx   = int_mb(int_mb(indx_start_list(1) +2*nb)+j-1)
c            send data....
            if (msglen.gt.0) then
               call SND(9+MSGDBL,
     >                  A(indx),
     >                  mdtob(msglen),
     >                  Parallel3d_convert_taskid_i(pto),
     >                  1)
            end if
   

         end do
      end if

      if (log_mb(receiver_list(1)+nb)) then
         do j=1,int_mb(npacket_list(1)+nb)
            pfrom  = int_mb(int_mb(proc_from_list(1)  +2*nb)+j-1)    
            msglen = int_mb(int_mb(packet_size_list(1)+2*nb)+j-1)
            indx   = int_mb(int_mb(indx_start_list(1) +2*nb)+j-1)
c            recieve data....
            if (msglen.gt.0) then
               call RCV(9+MSGDBL,
     >                  A(indx),
     >                  mdtob(msglen),rcv_len,
     >                  Parallel3d_convert_taskid_i(pfrom),
     >                  rcv_proc,1)
            end if


         end do
      end if

      return
      end


*     ************************************
*     *                                  *
*     *         c_Balance_c_unbalance    *
*     *                                  *
*     ************************************

      subroutine c_Balance_c_unbalance(nb,A)
      implicit none
      integer nb
      complex*16 A(*)

#include "bafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "balance_common.fh"

*     **** local variables ****
      integer  rcv_len,rcv_proc
      integer j
      integer pto,pfrom,msglen,indx

*     **** external functions ****
      integer  Parallel3d_convert_taskid_i
      external Parallel3d_convert_taskid_i



      if (log_mb(sender_list(1)+nb)) then
         do j=1,int_mb(npacket_list(1)+nb)
            pfrom  = int_mb(int_mb(proc_to_list(1)    +2*nb)+j-1)
            msglen = int_mb(int_mb(packet_size_list(1)+2*nb)+j-1)
            indx   = int_mb(int_mb(indx_start_list(1) +2*nb)+j-1)
c            recieve data....
            if (msglen.gt.0) then
               call RCV(9+MSGDBL,
     >                  A(indx),
     >                  mdtob(2*msglen),rcv_len,
     >                  Parallel3d_convert_taskid_i(pfrom),
     >                  rcv_proc,1)
            end if
  
         end do
      end if

      if (log_mb(receiver_list(1)+nb)) then
         do j=1,int_mb(npacket_list(1)+nb)
            pto    = int_mb(int_mb(proc_from_list(1)  +2*nb)+j-1)    
            msglen = int_mb(int_mb(packet_size_list(1)+2*nb)+j-1)
            indx   = int_mb(int_mb(indx_start_list(1) +2*nb)+j-1)
c            send data....
            if (msglen.gt.0) then
               call SND(9+MSGDBL,
     >                  A(indx),
     >                  mdtob(2*msglen),
     >                  Parallel3d_convert_taskid_i(pto),1)
            end if

         end do
      end if

      return
      end


*     ************************************
*     *                                  *
*     *         c_Balance_i_balance      *
*     *                                  *
*     ************************************

      subroutine c_Balance_i_balance(nb,A)
      implicit none
      integer nb
      integer A(*)


#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "bafdecls.fh"
#include "balance_common.fh"

*     **** local variables ****
      integer  rcv_len,rcv_proc
      integer j
      integer pto,pfrom,msglen,indx

*     **** external functions ****
      integer  Parallel3d_convert_taskid_i
      external Parallel3d_convert_taskid_i

      if (log_mb(sender_list(1)+nb)) then
         do j=1,int_mb(npacket_list(1)+nb)
            pto    = int_mb(int_mb(proc_to_list(1)    +2*nb)+j-1)    
            msglen = int_mb(int_mb(packet_size_list(1)+2*nb)+j-1)
            indx   = int_mb(int_mb(indx_start_list(1) +2*nb)+j-1)
c            send data....
            if (msglen.gt.0) then
               call SND(9+MSGINT,
     >                  A(indx),
     >                  mitob(msglen),
     >                  Parallel3d_convert_taskid_i(pto),1)
            end if
  

         end do
      end if

      if (log_mb(receiver_list(1)+nb)) then
         do j=1,int_mb(npacket_list(1)+nb)
            pfrom  = int_mb(int_mb(proc_from_list(1)  +2*nb)+j-1)    
            msglen = int_mb(int_mb(packet_size_list(1)+2*nb)+j-1)
            indx   = int_mb(int_mb(indx_start_list(1) +2*nb)+j-1)
c            recieve data....
            if (msglen.gt.0) then
               call RCV(9+MSGINT,
     >                  A(indx),
     >                  mitob(msglen),rcv_len,
     >                  Parallel3d_convert_taskid_i(pfrom),
     >                  rcv_proc,1)
            end if

         end do
      end if

      return
      end

