*
* $Id$
*

#define NBLOCKS 4


*     ***********************************
*     *					*
*     *	        Dne_init		*	
*     *					*
*     ***********************************

      subroutine Dne_init(ispin_in,ne_in,map_in)
      implicit none
      integer ispin_in
      integer ne_in(2)
      integer map_in


#include "bafdecls.fh"
#include "errquit.fh"
#include "Dne.fh"

*     **** local variables ****
      logical value
      integer ms,k,i,j,npack1,nida1,n2ft3d,nthr
      integer nework,ncqmax0
      integer ii,jj,icur,jcur,indx0,indx1,indx2,nn

*     **** external functions ****
      logical  control_mparallelized
      external control_mparallelized
      integer  Parallel2d_comm_i,Parallel2d_comm_j
      external Parallel2d_comm_i,Parallel2d_comm_j
      integer  control_mreplicate_size,Parallel_maxthreads
      external control_mreplicate_size,Parallel_maxthreads

      ispin   = ispin_in
      ne(1)   = ne_in(1)
      ne(2)   = ne_in(2)
      

      neq(1) = 0
      neq(2) = 0
      call D1dB_init(1,ne(1),map_in)
      call D1dB_nq(1,neq(1))
      if ((ispin.eq.2).and.(ne(2).gt.0)) then
        call D1dB_init(2,ne(2),map_in)
        call D1dB_nq(2,neq(2))
      end if

      nthr =  Parallel_maxthreads()
      call Parallel2d_np_i(np_i)
      call Parallel2d_np_j(np_j)
      call Parallel2d_taskid_i(taskid_i)
      call Parallel2d_taskid_j(taskid_j)
      comm_i = Parallel2d_comm_i()
      comm_j = Parallel2d_comm_j()
      parallelized = (np_j.gt.1)
      mparallelized = (parallelized.and.control_mparallelized())
      mreplicate = ((np_i*np_j).gt.control_mreplicate_size())

      if (parallelized) then
         value = .true.
         mcq(1) = 0
         mcq(2) = 0
         ncq(1) = 0
         ncq(2) = 0
         mcqmax(1) = 0
         mcqmax(2) = 0
         ncqmax(1) = 0
         ncqmax(2) = 0
         ncqmax0 = 0
         do ms=1,ispin
             value = value.and.
     >               BA_alloc_get(mt_int,np_i,'ma',ma(2,ms),ma(1,ms))
             value = value.and.
     >               BA_alloc_get(mt_int,np_i,'ma1',ma1(2,ms),ma1(1,ms))
             value = value.and.
     >               BA_alloc_get(mt_int,np_i,'ma2',ma2(2,ms),ma2(1,ms))
             value = value.and.
     >               BA_alloc_get(mt_int,np_i,'mc',mc(2,ms),mc(1,ms))
             value = value.and.
     >               BA_alloc_get(mt_int,np_j,'na',na(2,ms),na(1,ms))
             value = value.and.
     >               BA_alloc_get(mt_int,np_j,'nc',nc(2,ms),nc(1,ms))
             if (.not.value) then
               call errquit('Dne_init: out of heap memory',0,MA_ERR)
             end if

             call icopy(np_i,0,0,int_mb(ma(1,ms)),1)
             call icopy(np_i,0,0,int_mb(ma1(1,ms)),1)
             call icopy(np_i,0,0,int_mb(ma2(1,ms)),1)
             call icopy(np_i,0,0,int_mb(mc(1,ms)),1)
             call icopy(np_j,0,0,int_mb(na(1,ms)),1)
             call icopy(np_j,0,0,int_mb(nc(1,ms)),1)

             i = 0
             j = 0
             do k=1,ne(ms)
                int_mb(mc(1,ms)+i) = int_mb(mc(1,ms)+i) + 1

                int_mb(nc(1,ms)+j) = int_mb(nc(1,ms)+j) + 1
                int_mb(na(1,ms)+j) = int_mb(na(1,ms)+j) + 1
                i = mod(i+1,np_i)
                j = mod(j+1,np_j)
             end do

             call Pack_npack(1,npack1)
             call Pack_nida(1, nida1)
             call D3dB_n2ft3d(1,n2ft3d)
             int_mb(ma(1,ms) +taskid_i) = 2*npack1
             int_mb(ma1(1,ms)+taskid_i) = 2*nida1
             int_mb(ma2(1,ms)+taskid_i) = n2ft3d
             call D3dB_Vector_iSumAll(np_i,int_mb(ma(1,ms)))
             call D3dB_Vector_iSumAll(np_i,int_mb(ma1(1,ms)))
             call D3dB_Vector_iSumAll(np_i,int_mb(ma2(1,ms)))

             mcq(ms) = int_mb(mc(1,ms)+taskid_i)
             ncq(ms) = int_mb(nc(1,ms)+taskid_j)
             mcqmax(ms) = 0
             do i=0,np_i-1
                if (int_mb(mc(1,ms)+i).gt.mcqmax(ms)) 
     >             mcqmax(ms) = int_mb(mc(1,ms)+i) 
             end do
             ncqmax(ms) = 0
             do j=0,np_j-1
                if (int_mb(nc(1,ms)+j).gt.ncqmax(ms)) 
     >             ncqmax(ms) = int_mb(nc(1,ms)+j) 
             end do

             ncqmax0 = 0
             do j=0,np_j-1
                if (int_mb(nc(1,1)+j).gt.ncqmax0) 
     >             ncqmax0 = int_mb(nc(1,1)+j) 
             end do

             npack1_all = 0
             nida1_all  = 0
             n2ft3d_all = 0
             do i=0,np_i-1
                npack1_all = npack1_all + int_mb(ma(1,ms) +i)
                nida1_all  = nida1_all  + int_mb(ma1(1,ms)+i)
                n2ft3d_all = n2ft3d_all + int_mb(ma2(1,ms)+i)
             end do

         end do

         mdq = 0
         ndq = 0
         if ((ispin.eq.2).and.(ne(2).gt.0)) then
            value = BA_alloc_get(mt_int,np_i,'md',md(2),md(1))
            value = value.and.
     >              BA_alloc_get(mt_int,np_j,'nd',nd(2),nd(1))
            if (.not.value) then
               call errquit('Dne_init: out of heap memory',1,MA_ERR)
            end if
            call icopy(np_i,0,0,int_mb(md(1)),1)
            call icopy(np_j,0,0,int_mb(nd(1)),1)

            i = 0
            do k=1,ne(1)
               int_mb(md(1)+i) = int_mb(md(1)+i) + 1
               i = mod(i+1,np_i)
            end do

            j = 0
            do k=1,ne(2)
               int_mb(nd(1)+j) = int_mb(nd(1)+j) + 1
               j = mod(j+1,np_j)
            end do
            mdq = int_mb(md(1)+taskid_i)
            ndq = int_mb(nd(1)+taskid_j)
         end if


         value = value.and.BA_alloc_get(mt_dbl,
     >                                 2*2*64*int_mb(ma2(1,1)+taskid_i),
     >                                  'work1',work1(2),work1(1))

         nework = 2*64*int_mb(nc(1,1)+taskid_j)
         if (nework.lt.mcq(1)*ncq(1)) nework = mcq(1)*ncq(1)
         value = value.and.BA_alloc_get(mt_dbl,
     >                                  3*nework,
     >                                  'work2',work2(2),work2(1))

         value = value.and.BA_alloc_get(mt_dbl,
     >                              ne(1)*ncqmax0,
     >                              'bcolwork',bcolwork(2),bcolwork(1))
         value = value.and.BA_alloc_get(mt_dbl,
     >                              ne(1)*ncqmax0,
     >                              'bwork2',bwork2(2),bwork2(1))
         nework = int_mb(na(1,1))
         value = value.and.BA_alloc_get(mt_dbl,
     >                              2*nework*int_mb(ma2(1,1)+taskid_i),
     >                              'rwork1',rwork1(2),rwork1(1))
         value = value.and.BA_alloc_get(mt_dbl,
     >                              2*nework*int_mb(ma2(1,1)+taskid_i),
     >                              'rwork2',rwork2(2),rwork2(1))
         nework = mcq(1)*ncq(1) + mcq(2)*ncq(2)
         value = value.and.BA_alloc_get(mt_dbl,3*nework,
     >                              'mat_tmp',mat_tmp(2),mat_tmp(1))

         nework = mcq(1)*ncq(1) + mcq(2)*ncq(2)
         value = value.and.BA_alloc_get(mt_int,nework,
     >                                 'mindx0',mindx(2,0),mindx(1,0))
         nework = mcq(1)*ncq(1)
         value = value.and.BA_alloc_get(mt_int,nework,
     >                                 'mindx1',mindx(2,1),mindx(1,1))
         if (ispin.eq.2) then
         nework = mcq(2)*ncq(2)
         value = value.and.BA_alloc_get(mt_int,nework,
     >                                 'mindx2',mindx(2,2),mindx(1,2))
         endif
         if (.not.value) 
     >      call errquit('Dne_init: out of heap memory',1,MA_ERR)

         mall(0)  = ne(1)*ne(1) + ne(2)*ne(2)
         mall(1)  = ne(1)*ne(1)
         mall(2)  = ne(2)*ne(2)
         mpack(0) = mcq(1)*ncq(1) + mcq(2)*ncq(2)
         mpack(1) = mcq(1)*ncq(1)
         mpack(2) = mcq(2)*ncq(2)
         indx0=0
         indx1=0
         indx2=0
         jj   = 1
         jcur = 0
         do j=1,ne(1)
            ii   = 1
            icur = 0
            do i=1,ne(1)
               if ((icur.eq.taskid_i).and.(jcur.eq.taskid_j)) then
                  int_mb(mindx(1,0)+indx0) = i+(j-1)*ne(1)
                  int_mb(mindx(1,1)+indx1) = i+(j-1)*ne(1)
                  indx0 = indx0 + 1
                  indx1 = indx1 + 1
               end if
               ii = ii+1
               if (ii.gt.int_mb(mc(1,1)+icur)) then
                  icur = icur + 1
                  ii   = 1
               end if
            end do
            jj = jj+1
            if (jj.gt.int_mb(nc(1,1)+jcur)) then
               jcur = jcur + 1
               jj   = 1
            end if
         end do
         if (ispin.eq.2) then
         jj   = 1
         jcur = 0
         do j=1,ne(2)
            ii   = 1
            icur = 0
            do i=1,ne(2)
               if ((icur.eq.taskid_i).and.(jcur.eq.taskid_j)) then
                  int_mb(mindx(1,0)+indx0) = i+(j-1)*ne(2) + ne(1)*ne(1)
                  int_mb(mindx(1,2)+indx2) = i+(j-1)*ne(2)
                  indx0 = indx0 + 1
                  indx2 = indx2 + 1
               end if
               ii = ii+1
               if (ii.gt.int_mb(mc(1,2)+icur)) then
                  icur = icur + 1
                  ii   = 1
               end if
            end do
            jj = jj+1
            if (jj.gt.int_mb(nc(1,2)+jcur)) then
               jcur = jcur + 1
               jj   = 1
            end if
         end do
         end if

         if ((ispin.eq.2).and.(ne(2).gt.0)) then
            nework = mdq*ndq
            value = value.and.BA_alloc_get(mt_int,nework,
     >                                 'smindx',smindx(2),smindx(1))
            small  = ne(1)*ne(2)
            smpack = mdq*ndq
            indx2=0
            jj   = 1
            jcur = 0
            do j=1,ne(2)
               ii   = 1
               icur = 0
               do i=1,ne(1)
               if ((icur.eq.taskid_i).and.(jcur.eq.taskid_j)) then
                  int_mb(smindx(1)+indx2) = i+(j-1)*ne(1)
                  indx2 = indx2 + 1
               end if
               ii = ii+1
               if (ii.gt.int_mb(md(1)+icur)) then
                  icur = icur + 1
                  ii   = 1
               end if
               end do
               jj = jj+1
               if (jj.gt.int_mb(nd(1)+jcur)) then
                  jcur = jcur + 1
                  jj   = 1
               end if
            end do
         end if

      end if
      if (parallelized) then
         nn =  (mcqmax(1)*ncqmax(1)+mcqmax(2)*ncqmax(2))*nthr
      else
         nn = (ne(1)*ne(1)+ne(2)*ne(2))*nthr
      end if
      value = BA_alloc_get(mt_dbl,nn,'thrwork1',thrwork1(2),thrwork1(1))
      if (.not.value) call errquit('Dne_init:out of heap',9,MA_ERR)

      return
      end


*     ***********************************
*     *					*
*     *	          Dne_end   		*	
*     *					*
*     ***********************************

      subroutine Dne_end()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"
#include "Dne.fh"

*     ***** local variables ****
      logical value
      integer ms

      call D1dB_end(1)
      value = BA_free_heap(thrwork1(2))
      if (.not.value) 
     >   call errquit('Dne_end:error deallocating heap',0,MA_ERR)

      if ((ispin.eq.2).and.(ne(2).gt.0)) call D1dB_end(2)
      if (parallelized) then
         value = .true.
         do ms=1,ispin
            value = value.and.BA_free_heap(ma(2,ms))
            value = value.and.BA_free_heap(ma1(2,ms))
            value = value.and.BA_free_heap(ma2(2,ms))
            value = value.and.BA_free_heap(na(2,ms))
            value = value.and.BA_free_heap(mc(2,ms))
            value = value.and.BA_free_heap(nc(2,ms))
         end do
         if ((ispin.eq.2).and.(ne(2).gt.0)) then
            value = value.and.BA_free_heap(md(2))
            value = value.and.BA_free_heap(nd(2))
            value = value.and.BA_free_heap(smindx(2))
         end if
         value = value.and.BA_free_heap(work1(2))
         value = value.and.BA_free_heap(work2(2))
         value = value.and.BA_free_heap(bcolwork(2))
         value = value.and.BA_free_heap(bwork2(2))
         value = value.and.BA_free_heap(rwork1(2))
         value = value.and.BA_free_heap(rwork2(2))
         value = value.and.BA_free_heap(mat_tmp(2))
         value = value.and.BA_free_heap(mindx(2,0))
         value = value.and.BA_free_heap(mindx(2,1))
         if (ispin.eq.2) value = value.and.BA_free_heap(mindx(2,2))
         if (.not.value) then
           call errquit('Dne_end:error deallocating heap',1,MA_ERR)
         end if
      end if

      return
      end
         

*     ***********************************
*     *					*
*     *	          Dneall_ntoqp 		*	
*     *					*
*     ***********************************

      subroutine Dneall_ntoqp(n,q,p)      
      implicit none
      integer n,q,p

#include "Dne.fh"

      if (n.le.ne(1)) then
        call D1dB_ktoqp(1,n,q,p)
      else
        call D1dB_ktoqp(2,n-ne(1),q,p)
        q = q + neq(1)
      end if

      return
      end




*     ***********************************
*     *                                 *
*     *           Dneall_qton           *
*     *                                 *
*     ***********************************

      subroutine Dneall_qton(q,n)
      implicit none
      integer n,q

#include "Dne.fh"

      if (q.le.neq(1)) then
        call D1dB_qtok(1,q,n)
      else
        call D1dB_qtok(2,q-neq(1),n)
        n = n + ne(1)
      end if

      return
      end



*     ***********************************
*     *					*
*     *	        Dneall_neq		*	
*     *					*
*     ***********************************

      subroutine Dneall_neq(nqtmp)
      implicit none
      integer nqtmp(2)

#include "Dne.fh"

      nqtmp(1) = neq(1)
      nqtmp(2) = neq(2)
      return 
      end


*     ***********************************
*     *					*
*     *	        Dneall_ne 		*	
*     *					*
*     ***********************************

      subroutine Dneall_ne(netmp)
      implicit none
      integer netmp(2)

#include "Dne.fh"

      netmp(1) = ne(1)
      netmp(2) = ne(2)
      return 
      end

     
*     ***********************************
*     *                                 *
*     *         Dneall_ispin            *       
*     *                                 *
*     ***********************************
      
      subroutine Dneall_ispin(ispintmp)
      implicit none
      integer ispintmp
         
#include "Dne.fh"
           
      ispintmp = ispin
      return
      end




*     ***********************************
*     *                                 *
*     *         Dneall_na               *
*     *                                 *
*     ***********************************

      integer function Dneall_na(ms,jj)
      implicit none
      integer ms,jj

#include "bafdecls.fh"
#include "Dne.fh"

      Dneall_na = int_mb(na(1,ms)+jj)
      return
      end



*     ***********************************
*     *                                 *
*     *         Dneall_na_ptr           *
*     *                                 *
*     ***********************************

      integer function Dneall_na_ptr(ms)
      implicit none
      integer ms

#include "Dne.fh"

      Dneall_na_ptr = na(1,ms)
      return
      end


c
c
c
c     ****************************************
c     *                                      *
c     *        Dneall_m_size                 *
c     *                                      *
c     ****************************************
      subroutine Dneall_m_size(mb,size)
      implicit none
      integer mb
      integer size

#include "Dne.fh"
     
      if (mparallelized) then
         if (mb.eq.0) then
            size = mcq(1)*ncq(1) + mcq(2)*ncq(2)
         else
            size = mcq(mb)*ncq(mb)
         end if
      else
         if (mb.eq.0) then 
            size = ne(1)*ne(1) + ne(2)*ne(2)
         else
            size = ne(mb)*ne(mb)
         end if
      end if
     
      return
      end

c
c
c
c     ****************************************
c     *                                      *
c     *        Dneall_m_allocate_block       *
c     *                                      *
c     ****************************************
      logical function Dneall_m_allocate_block(mb,nb,hml)
      implicit none
      integer mb,nb
      integer hml(2)

#include "Dne.fh"
#include "bafdecls.fh"

      integer size
     
      if (mparallelized) then
         if (mb.eq.0) then
            size = mcq(1)*ncq(1) + mcq(2)*ncq(2)
         else
            size = mcq(mb)*ncq(mb)
         end if
      else
         if (mb.eq.0) then 
            size = ne(1)*ne(1) + ne(2)*ne(2)
         else
            size = ne(mb)*ne(mb)
         end if
      end if
     
      Dneall_m_allocate_block 
     > = BA_alloc_get(mt_dbl,nb*size,'hmlab',hml(2),hml(1))
      return
      end
c
c
c
c     ****************************************
c     *                                      *
c     *        Dneall_m_allocate             *
c     *                                      *
c     ****************************************
      logical function Dneall_m_allocate(mb,hml)
      implicit none
      integer mb
      integer hml(2)

#include "Dne.fh"
#include "bafdecls.fh"

      integer size
     
      if (mparallelized) then
         if (mb.eq.0) then
            size = mcq(1)*ncq(1) + mcq(2)*ncq(2)
         else
            size = mcq(mb)*ncq(mb)
         end if
      else
         if (mb.eq.0) then 
            size = ne(1)*ne(1) + ne(2)*ne(2)
         else
            size = ne(mb)*ne(mb)
         end if
      end if
     
      Dneall_m_allocate = BA_alloc_get(mt_dbl,size,'hmla',hml(2),hml(1))
      return
      end


c
c
c
c     ****************************************
c     *                                      *
c     *        Dneall_w_allocate             *
c     *                                      *
c     ****************************************
      logical function Dneall_w_allocate(mb,hml)
      implicit none
      integer mb
      integer hml(2)

#include "Dne.fh"
#include "bafdecls.fh"

      integer size

      if (mparallelized) then
         if (mb.eq.0) then
            size = mcq(1)*ncq(1) + mcq(2)*ncq(2)
         else
            size = mcq(mb)*ncq(mb)
         end if
      else
         if (mb.eq.0) then
            size = ne(1)*ne(1) + ne(2)*ne(2)
         else
            size = ne(mb)*ne(mb)
         end if
      end if

      Dneall_w_allocate=BA_alloc_get(mt_dcpl,size,'hmla',hml(2),hml(1))
      return
      end


c
c
c
c     ****************************************
c     *                                      *
c     *        Dneall_w_push_get             *
c     *                                      *
c     ****************************************
      logical function Dneall_w_push_get(mb,hml)
      implicit none
      integer mb
      integer hml(2)

#include "Dne.fh"
#include "bafdecls.fh"

      integer size

      if (mparallelized) then
         if (mb.eq.0) then
            size = mcq(1)*ncq(1) + mcq(2)*ncq(2)
         else 
            size = mcq(mb)*ncq(mb)
         end if
      else
         if (mb.eq.0) then
            size = ne(1)*ne(1) + ne(2)*ne(2)
         else
            size = ne(mb)*ne(mb)
         end if
      end if

      Dneall_w_push_get=BA_push_get(mt_dcpl,size,'hmla',hml(2),hml(1))
      return
      end


c
c     ****************************************
c     *                                      *
c     *        Dneall_m_free                 *
c     *                                      *
c     ****************************************
      logical function Dneall_m_free(hml)
      implicit none
      integer hml(2)

#include "bafdecls.fh"

      Dneall_m_free = BA_free_heap(hml(2))
      return
      end


c
c     ****************************************
c     *                                      *
c     *        Dneall_w_free                 *
c     *                                      *
c     ****************************************
      logical function Dneall_w_free(hml)
      implicit none
      integer hml(2)

#include "bafdecls.fh"

      Dneall_w_free = BA_free_heap(hml(2))
      return
      end

c
c     ****************************************
c     *                                      *
c     *        Dneall_m_push_get             *
c     *                                      *
c     ****************************************
      logical function Dneall_m_push_get(mb,hml)
      implicit none
      integer mb
      integer hml(2)

#include "Dne.fh"
#include "bafdecls.fh"

      integer size

      if (mparallelized) then
         if (mb.eq.0) then
            size = mcq(1)*ncq(1) + mcq(2)*ncq(2)
         else
            size = mcq(mb)*ncq(mb)
         end if
      else
         if (mb.eq.0) then
            size = ne(1)*ne(1) + ne(2)*ne(2)
         else
            size = ne(mb)*ne(mb)
         end if
      end if

      Dneall_m_push_get = BA_push_get(mt_dbl,size,'hmls',hml(2),hml(1))
      return
      end


c
c     ****************************************
c     *                                      *
c     *        Dneall_m_push_get_block       *
c     *                                      *
c     ****************************************
      logical function Dneall_m_push_get_block(mb,nb,hml)
      implicit none
      integer mb,nb
      integer hml(2)

#include "Dne.fh"
#include "bafdecls.fh"

      integer size

      if (mparallelized) then
         if (mb.eq.0) then
            size = mcq(1)*ncq(1) + mcq(2)*ncq(2)
         else
            size = mcq(mb)*ncq(mb)
         end if
      else
         if (mb.eq.0) then
            size = ne(1)*ne(1) + ne(2)*ne(2)
         else
            size = ne(mb)*ne(mb)
         end if
      end if

      Dneall_m_push_get_block 
     >  = BA_push_get(mt_dbl,nb*size,'hmlsb',hml(2),hml(1))
      return
      end


c     ****************************************
c     *                                      *
c     *        Dneall_m_pop_stack            *
c     *                                      *
c     ****************************************
      logical function Dneall_m_pop_stack(hml)
      implicit none
      integer hml(2)

#include "bafdecls.fh"

      Dneall_m_pop_stack = BA_pop_stack(hml(2))
      return
      end



c     ****************************************
c     *                                      *
c     *        Dneall_w_pop_stack            *
c     *                                      *
c     ****************************************
      logical function Dneall_w_pop_stack(hml)
      implicit none
      integer hml(2)

#include "bafdecls.fh"

      Dneall_w_pop_stack = BA_pop_stack(hml(2))
      return
      end


c
c
c
c     ****************************************
c     *                                      *
c     *        Dneall_sm_size                *
c     *                                      *
c     ****************************************
      subroutine Dneall_sm_size(size)
      implicit none
      integer size

#include "Dne.fh"

      if (mparallelized) then
         size = mdq*ndq
      else
         size = ne(1)*ne(2)
      end if
      return
      end

c     ****************************************
c     *                                      *
c     *        Dneall_sm_push_get            *
c     *                                      *
c     ****************************************
      logical function Dneall_sm_push_get(hml)
      implicit none
      integer hml(2)

#include "Dne.fh"
#include "bafdecls.fh"

      integer size

      if (mparallelized) then
         size = mdq*ndq
      else
         size = ne(1)*ne(2)
      end if

      Dneall_sm_push_get = BA_push_get(mt_dbl,size,'hmls',hml(2),hml(1))
      return
      end


c     ****************************************
c     *                                      *
c     *        Dneall_sm_pop_stack           *
c     *                                      *
c     ****************************************
      logical function Dneall_sm_pop_stack(hml)
      implicit none
      integer hml(2)

#include "bafdecls.fh"

      Dneall_sm_pop_stack = BA_pop_stack(hml(2))
      return
      end



c     ****************************************
c     *                                      *
c     *        Dneall_fmf_Multiply           *
c     *                                      *
c     ****************************************

*  uses rotation algorithm

      subroutine Dneall_fmf_Multiply(mb,Ain,npack1,
     >                                   hml,alpha,
     >                                   Aout,beta)
      implicit none
      integer    mb
      complex*16 Ain(*)
      integer    npack1
      real*8     hml(*)
      real*8     alpha
      complex*16 Aout(*)
      real*8     beta

#include "bafdecls.fh"
#include "Dne.fh"

*     **** local variables ****
      integer ms,ms1,ms2,n,shift,shift2,shift3,ishift2,ishift3

c      call Pack_npack(1,npack1)

      call nwpw_timing_start(16)
      if (parallelized) then
       if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = mcq(1)*ncq(1)
            ishift3 = ne(1)*ne(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
            ishift3 = 0
         end if

         if (mparallelized) then
         do ms=ms1,ms2
            shift  = 1 + (ms-1)*neq(1)*npack1
            shift2 = 1 + (ms-1)*ishift2
            call DMatrix_dgemm1_rot(npack1_all,ne(ms),ne(ms),
     >             alpha,
     >             Ain(shift),int_mb(ma(1,ms)+taskid_i),
     >                        int_mb(ma(1,ms)),
     >                        int_mb(na(1,ms)),
     >             hml(shift2),mcq(ms),
     >                         int_mb(mc(1,ms)),
     >                         int_mb(nc(1,ms)),
     >             beta,
     >             Aout(shift),int_mb(ma(1,ms)+taskid_i),
     >                         int_mb(ma(1,ms)),
     >                         int_mb(na(1,ms)),
     >             taskid_i,taskid_j,
     >             np_i,np_j,
     >             comm_i, comm_j,
     >             dbl_mb(bcolwork(1)),dbl_mb(bwork2(1)),
     >             dbl_mb(rwork1(1)),dbl_mb(rwork2(1)))
         end do
         else
c         call Dneall_m_scatter(mpack(mb),int_mb(mindx(1,mb)),
c     >                         hml,dbl_mb(mat_tmp(1)))
         do ms=ms1,ms2
            shift  = 1 + (ms-1)*neq(1)*npack1
c            shift2 =     (ms-1)*ishift2
c            call DMatrix_dgemm1_rot(npack1_all,ne(ms),ne(ms),
c     >             alpha,
c     >             Ain(shift),int_mb(ma(1,ms)+taskid_i),
c     >                        int_mb(ma(1,ms)),
c     >                        int_mb(na(1,ms)),
c     >             dbl_mb(mat_tmp(1)+shift2),mcq(ms),
c     >                         int_mb(mc(1,ms)),
c     >                         int_mb(nc(1,ms)),
c     >             beta,
c     >             Aout(shift),int_mb(ma(1,ms)+taskid_i),
c     >                         int_mb(ma(1,ms)),
c     >                         int_mb(na(1,ms)),
c     >             taskid_i,taskid_j,
c     >             np_i,np_j,
c     >             comm_i, comm_j,
c     >             dbl_mb(bcolwork(1)),dbl_mb(bwork2(1)),
c     >             dbl_mb(rwork1(1)),dbl_mb(rwork2(1)))
            shift3 = 1 + (ms-1)*ishift3
            call DMatrix_dgemm1_rot2(npack1_all,ne(ms),ne(ms),
     >             alpha,
     >             Ain(shift),int_mb(ma(1,ms)+taskid_i),
     >                        int_mb(ma(1,ms)),
     >                        int_mb(na(1,ms)),
     >             hml(shift3),mcq(ms),
     >                         int_mb(mc(1,ms)),
     >                         int_mb(nc(1,ms)),
     >             beta,
     >             Aout(shift),int_mb(ma(1,ms)+taskid_i),
     >                         int_mb(ma(1,ms)),
     >                         int_mb(na(1,ms)),
     >             taskid_i,taskid_j,
     >             np_i,np_j,
     >             comm_i, comm_j,
     >             dbl_mb(bcolwork(1)),dbl_mb(bwork2(1)),
     >             dbl_mb(rwork1(1)),dbl_mb(rwork2(1)))
         end do
         end if


      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = ne(1)*ne(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if

         do ms=ms1,ms2
            n     = ne(ms)
            if (n.le.0) go to 30
            shift  = 1 + (ms-1)*ne(1)*npack1
            shift2 = 1 + (ms-1)*ishift2
            call DGEMM_OMP('N','N',2*npack1,n,n,
     >                (alpha),
     >                Ain(shift),  2*npack1,
     >                hml(shift2),    n,
     >                (beta),
     >                Aout(shift),2*npack1)
   30       continue
         end do
      end if

      call nwpw_timing_end(16)
      return
      end


c     ****************************************
c     *                                      *
c     *        Dneall_fmf_Multiply_summa     *
c     *                                      *
c     ****************************************

*  uses summa algorithm

      subroutine Dneall_fmf_Multiply_summa(mb,Ain,npack1,
     >                                     hml,alpha,
     >                                     Aout,beta)
      implicit none
      integer    mb
      complex*16 Ain(*)
      integer    npack1
      real*8     hml(*)
      real*8     alpha
      complex*16 Aout(*)
      real*8     beta
           
#include "bafdecls.fh"
#include "Dne.fh"

*     **** local variables ****
      integer ms,ms1,ms2,n,shift,shift2,ishift2

c      call Pack_npack(1,npack1)

      call nwpw_timing_start(16)
      if (parallelized) then
       if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = mcq(1)*ncq(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if

         if (mparallelized) then
         do ms=ms1,ms2
            shift  = 1 + (ms-1)*neq(1)*npack1
            shift2 = 1 + (ms-1)*ishift2
            call DMatrix_dgemm1(npack1_all,ne(ms),ne(ms),64,
     >             alpha,
     >             Ain(shift),int_mb(ma(1,ms)+taskid_i),
     >                        int_mb(ma(1,ms)),
     >                        int_mb(na(1,ms)),
     >             hml(shift2),mcq(ms),
     >                         int_mb(mc(1,ms)),
     >                         int_mb(nc(1,ms)),
     >             beta,
     >             Aout(shift),int_mb(ma(1,ms)+taskid_i),
     >                         int_mb(ma(1,ms)),
     >                         int_mb(na(1,ms)),
     >             taskid_i,taskid_j,
     >             np_i,np_j,
     >             comm_i, comm_j,
     >             dbl_mb(work1(1)),dbl_mb(work2(1)))

         end do
         else
         call Dneall_m_scatter(mpack(mb),int_mb(mindx(1,mb)),
     >                         hml,dbl_mb(mat_tmp(1)))
         do ms=ms1,ms2
            shift  = 1 + (ms-1)*neq(1)*npack1
            shift2 =     (ms-1)*ishift2
            call DMatrix_dgemm1(npack1_all,ne(ms),ne(ms),64,
     >             alpha,
     >             Ain(shift),int_mb(ma(1,ms)+taskid_i),
     >                        int_mb(ma(1,ms)),
     >                        int_mb(na(1,ms)),
     >             dbl_mb(mat_tmp(1)+shift2),mcq(ms),
     >                         int_mb(mc(1,ms)),
     >                         int_mb(nc(1,ms)),
     >             beta,
     >             Aout(shift),int_mb(ma(1,ms)+taskid_i),
     >                         int_mb(ma(1,ms)),
     >                         int_mb(na(1,ms)),
     >             taskid_i,taskid_j,
     >             np_i,np_j,
     >             comm_i, comm_j,
     >             dbl_mb(work1(1)),dbl_mb(work2(1)))

         end do
         end if
         
      else
         if (mb.eq.0) then
            ms1 = 1    
            ms2 = ispin
            ishift2 = ne(1)*ne(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if

         do ms=ms1,ms2
            n     = ne(ms)
            if (n.le.0) go to 30
            shift  = 1 + (ms-1)*ne(1)*npack1
            shift2 = 1 + (ms-1)*ishift2
            call DGEMM_OMP('N','N',2*npack1,n,n,
     >                (alpha),
     >                Ain(shift),  2*npack1,
     >                hml(shift2),    n,
     >                (beta),
     >                Aout(shift),2*npack1)
   30       continue
         end do
      end if

      call nwpw_timing_end(16)
      return
      end


c     ****************************************
c     *                                      *
c     *        Dneall_ffm_sym_Multiply       *
c     *                                      *
c     ****************************************

      subroutine Dneall_ffm_sym_Multiply(mb,A1,A2,npack1,hml)
      implicit none
      integer    mb
      complex*16 A1(*),A2(*)
      integer    npack1
      real*8     hml(*)
           
#include "bafdecls.fh"
#include "Dne.fh"

*     **** local variables ****
      integer ms,ms1,ms2,n,shift,shift2,ishift2
      integer tid,nthr,nn,i
      !integer mm

*     **** external functions ****
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads

      tid  = Parallel_threadid()
      nthr = Parallel_nthreads()

      call nwpw_timing_start(15)
      if (parallelized) then
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = mcq(1)*ncq(1)
            !mm = mcq(1)*ncq(1) + mcq(2)*ncq(2)
            nn = mcqmax(1)*ncqmax(1)+mcqmax(2)*ncqmax(2)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
            !mm = mcq(mb)*ncq(mb) 
            nn = mcqmax(mb)*ncqmax(mb)
         end if

         if (mparallelized) then
         do ms=ms1,ms2
            if (ne(ms).le.0) go to 20
            shift  = 1 + (ms-1)*neq(1)*npack1
            shift2 = 1 + (ms-1)*ishift2
            call DMatrix_dgemm2(ne(ms),ne(ms),npack1_all,128,
     >             2.0d0,
     >             A1(shift),int_mb(ma(1,ms)+taskid_i), 
     >                       int_mb(ma(1,ms)),
     >                       int_mb(na(1,ms)),
     >             A2(shift),int_mb(ma(1,ms)+taskid_i), 
     >                       int_mb(ma(1,ms)),
     >                       int_mb(na(1,ms)),
     >             0.0d0,
     >             hml(shift2),int_mb(mc(1,ms)+taskid_i), 
     >                         int_mb(mc(1,ms)),
     >                         int_mb(nc(1,ms)),
     >             taskid_i,taskid_j,
     >             np_i,np_j,
     >             comm_i, comm_j,
     >             dbl_mb(work1(1)),dbl_mb(work2(1)))

            call DMatrix_dgemm2(ne(ms),ne(ms),nida1_all,128,
     >             -1.0d0,
     >             A1(shift),int_mb(ma(1,ms)+taskid_i), 
     >                       int_mb(ma1(1,ms)),
     >                       int_mb(na(1,ms)),
     >             A2(shift),int_mb(ma(1,ms)+taskid_i), 
     >                       int_mb(ma1(1,ms)),
     >                       int_mb(na(1,ms)),
     >             1.0d0,
     >             hml(shift2),int_mb(mc(1,ms)+taskid_i), 
     >                         int_mb(mc(1,ms)),
     >                         int_mb(nc(1,ms)),
     >             taskid_i,taskid_j,
     >             np_i,np_j,
     >             comm_i, comm_j,
     >             dbl_mb(work1(1)),dbl_mb(work2(1)))

  20        continue
         end do
         else
         do ms=ms1,ms2
            if (ne(ms).le.0) go to 21
            shift  = 1 + (ms-1)*neq(1)*npack1
            shift2 =     (ms-1)*ishift2
            call DMatrix_dgemm2c_omp(ne(ms),ne(ms),npack1_all,128,
     >              A1(shift),A2(shift),int_mb(ma(1,ms)+taskid_i),
     >                                  int_mb(ma(1,ms)),
     >                                  int_mb(ma1(1,ms)),
     >                                  int_mb(na(1,ms)),
     >              dbl_mb(mat_tmp(1)+shift2),
     >                                        int_mb(mc(1,ms)+taskid_i),
     >                                        int_mb(mc(1,ms)),
     >                                        int_mb(nc(1,ms)),
     >              taskid_i,taskid_j,
     >              np_i,np_j,
     >              comm_i, comm_j,
     >              tid,nthr,
     >              dbl_mb(work1(1)),dbl_mb(work2(1)),
     >              dbl_mb(thrwork1(1)+tid*nn))

c            call DMatrix_dgemm2(ne(ms),ne(ms),npack1_all,128,
c     >             2.0d0,
c     >             A1(shift),int_mb(ma(1,ms)+taskid_i),
c     >                       int_mb(ma(1,ms)),
c     >                       int_mb(na(1,ms)),
c     >             A2(shift),int_mb(ma(1,ms)+taskid_i),
c     >                       int_mb(ma(1,ms)),
c     >                       int_mb(na(1,ms)),
c     >             0.0d0,
c     >          dbl_mb(mat_tmp(1)+shift2+mm),int_mb(mc(1,ms)+taskid_i),
c     >                         int_mb(mc(1,ms)),
c     >                         int_mb(nc(1,ms)),
c     >             taskid_i,taskid_j,
c     >             np_i,np_j,
c     >             comm_i, comm_j,
c     >             dbl_mb(work1(1)),dbl_mb(work2(1)))
c
c            call DMatrix_dgemm2(ne(ms),ne(ms),nida1_all,128,
c     >             -1.0d0,
c     >             A1(shift),int_mb(ma(1,ms)+taskid_i),
c     >                       int_mb(ma1(1,ms)),
c     >                       int_mb(na(1,ms)),
c     >             A2(shift),int_mb(ma(1,ms)+taskid_i),
c     >                       int_mb(ma1(1,ms)),
c     >                       int_mb(na(1,ms)),
c     >             1.0d0,
c     >           dbl_mb(mat_tmp(1)+shift2+mm),int_mb(mc(1,ms)+taskid_i),
c     >                         int_mb(mc(1,ms)),
c     >                         int_mb(nc(1,ms)),
c     >             taskid_i,taskid_j,
c     >             np_i,np_j,
c     >             comm_i, comm_j,
c     >             dbl_mb(work1(1)),dbl_mb(work2(1)))
c            write(*,*) mcq(1),ncq(1),mm,taskid_i,taskid_j,tid,
c     >                "matmp=",( dbl_mb(mat_tmp(1)+shift2+mm+i-1) 
c     >                          - dbl_mb(mat_tmp(1)+shift2+i-1), 
c     >                            i=1,mm)

  21        continue
         end do
         call Dneall_m_gather(mall(mb),mpack(mb),int_mb(mindx(1,mb)),
     >                        dbl_mb(mat_tmp(1)),hml)

         end if
      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = ne(1)*ne(1)
            nn = ne(1)*ne(1) + ne(2)*ne(2)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
            nn = ne(mb)*ne(mb)
         end if

         do ms=ms1,ms2
            shift  = 1 + (ms-1)*ne(1)*npack1
            shift2 = 1 + (ms-1)*ishift2
            n     = ne(ms)
            if (n.le.0) go to 30

            call Pack_ccm_sym_dot_omp(1,n,
     >                        A1(shift),
     >                        A2(shift),
     >                        hml(shift2),dbl_mb(thrwork1(1)+tid*nn))
c            call Pack_ccm_sym_dot(1,n,
c     >                        A1(shift),
c     >                        A2(shift),
c     >                        hml(shift2))
  30        continue
         end do
      end if

      call nwpw_timing_end(15)
      return
      end


c     ****************************************
c     *                                      *
c     *        Dneall_ffm_Multiply          *
c     *                                      *
c     ****************************************

      subroutine Dneall_ffm_Multiply(mb,A1,A2,npack1,hml)
      implicit none
      integer    mb
      complex*16 A1(*),A2(*)
      integer    npack1
      real*8     hml(*)
           
#include "bafdecls.fh"
#include "Dne.fh"

*     **** local variables ****
      integer ms,ms1,ms2,n,shift2,ishift2,shift
      integer tid,nthr,nn

*     **** external functions ****
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads

      tid  = Parallel_threadid()
      nthr = Parallel_nthreads()

      call nwpw_timing_start(15)
      if (parallelized) then
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = mcq(1)*ncq(1)
            nn = mcqmax(1)*ncqmax(1)+mcqmax(2)*ncqmax(2)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
            nn = mcqmax(mb)*ncqmax(mb)
         end if

         if (mparallelized) then
         do ms=ms1,ms2
            if (ne(ms).le.0) go to 20
            shift  = 1 + (ms-1)*neq(1)*npack1
            shift2 = 1 + (ms-1)*ishift2
            call DMatrix_dgemm2(ne(ms),ne(ms),npack1_all,128,
     >             2.0d0,
     >             A1(shift),int_mb(ma(1,ms)+taskid_i),
     >                       int_mb(ma(1,ms)),
     >                       int_mb(na(1,ms)),
     >             A2(shift),int_mb(ma(1,ms)+taskid_i),
     >                       int_mb(ma(1,ms)),
     >                       int_mb(na(1,ms)),
     >             0.0d0,
     >             hml(shift2),int_mb(mc(1,ms)+taskid_i),
     >                         int_mb(mc(1,ms)),
     >                         int_mb(nc(1,ms)),
     >             taskid_i,taskid_j,
     >             np_i,np_j,
     >             comm_i, comm_j,
     >             dbl_mb(work1(1)),dbl_mb(work2(1)))

            call DMatrix_dgemm2(ne(ms),ne(ms),nida1_all,128,
     >             -1.0d0,
     >             A1(shift),int_mb(ma(1,ms)+taskid_i),
     >                       int_mb(ma1(1,ms)),
     >                       int_mb(na(1,ms)),
     >             A2(shift),int_mb(ma(1,ms)+taskid_i),
     >                       int_mb(ma1(1,ms)),
     >                       int_mb(na(1,ms)),
     >             1.0d0,
     >             hml(shift2),int_mb(mc(1,ms)+taskid_i),
     >                         int_mb(mc(1,ms)),
     >                         int_mb(nc(1,ms)),
     >             taskid_i,taskid_j,
     >             np_i,np_j,
     >             comm_i, comm_j,
     >             dbl_mb(work1(1)),dbl_mb(work2(1)))

  20        continue
         end do
         else
         do ms=ms1,ms2
            if (ne(ms).le.0) go to 21
            shift  = 1 + (ms-1)*neq(1)*npack1
            shift2 =     (ms-1)*ishift2
            call DMatrix_dgemm2c_omp(ne(ms),ne(ms),npack1_all,128,
     >              A1(shift),A2(shift),int_mb(ma(1,ms)+taskid_i),
     >                                  int_mb(ma(1,ms)),
     >                                  int_mb(ma1(1,ms)),
     >                                  int_mb(na(1,ms)),
     >              dbl_mb(mat_tmp(1)+shift2),
     >                                        int_mb(mc(1,ms)+taskid_i),
     >                                        int_mb(mc(1,ms)),
     >                                        int_mb(nc(1,ms)),
     >              taskid_i,taskid_j,
     >              np_i,np_j,
     >              comm_i, comm_j,
     >              tid,nthr,
     >              dbl_mb(work1(1)),dbl_mb(work2(1)),
     >              dbl_mb(thrwork1(1)+tid*nn))
c            call DMatrix_dgemm2(ne(ms),ne(ms),npack1_all,128,
c     >             2.0d0,
c     >             A1(shift),int_mb(ma(1,ms)+taskid_i),
c     >                       int_mb(ma(1,ms)),
c     >                       int_mb(na(1,ms)),
c     >             A2(shift),int_mb(ma(1,ms)+taskid_i),
c     >                       int_mb(ma(1,ms)),
c     >                       int_mb(na(1,ms)),
c     >             0.0d0,
c     >             dbl_mb(mat_tmp(1)+shift2),int_mb(mc(1,ms)+taskid_i),
c     >                         int_mb(mc(1,ms)),
c     >                         int_mb(nc(1,ms)),
c     >             taskid_i,taskid_j,
c     >             np_i,np_j,
c     >             comm_i, comm_j,
c     >             dbl_mb(work1(1)),dbl_mb(work2(1)))
c
c            call DMatrix_dgemm2(ne(ms),ne(ms),nida1_all,128,
c     >             -1.0d0,
c     >             A1(shift),int_mb(ma(1,ms)+taskid_i),
c     >                       int_mb(ma1(1,ms)),
c     >                       int_mb(na(1,ms)),
c     >             A2(shift),int_mb(ma(1,ms)+taskid_i),
c     >                       int_mb(ma1(1,ms)),
c     >                       int_mb(na(1,ms)),
c     >             1.0d0,
c     >             dbl_mb(mat_tmp(1)+shift2),int_mb(mc(1,ms)+taskid_i),
c     >                         int_mb(mc(1,ms)),
c     >                         int_mb(nc(1,ms)),
c     >             taskid_i,taskid_j,
c     >             np_i,np_j,
c     >             comm_i, comm_j,
c     >             dbl_mb(work1(1)),dbl_mb(work2(1)))

  21        continue
         end do
         call Dneall_m_gather(mall(mb),mpack(mb),int_mb(mindx(1,mb)),
     >                        dbl_mb(mat_tmp(1)),hml)
         end if
         
      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = ne(1)*ne(1)
            nn = ne(1)*ne(1) + ne(2)*ne(2)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
            nn = ne(mb)*ne(mb)
         end if

         do ms=ms1,ms2
            shift  = 1+(ms-1)*ne(1)*npack1
            shift2 = 1+(ms-1)*ishift2
            n     = ne(ms)
            if (n.le.0) go to 30

            call Pack_ccm_dot_omp(1,n,
     >                        A1(shift),
     >                        A2(shift),
     >                        hml(shift2),dbl_mb(thrwork1(1)+tid*nn))
  30        continue
         end do
      end if

      call nwpw_timing_end(15)
      return
      end




c     ****************************************
c     *                                      *
c     *        Dneall_ffsm_Multiply          *
c     *                                      *
c     ****************************************

      subroutine Dneall_ffsm_Multiply(A1,A2,npack1,hml)
      implicit none
      complex*16 A1(*),A2(*)
      integer    npack1
      real*8     hml(*)
           
#include "bafdecls.fh"
#include "Dne.fh"

*     **** local variables ****
      integer ms,ms1,ms2,n,shift2,ishift2,shift

      call nwpw_timing_start(15)
      if (parallelized) then

         if (mparallelized) then

            if (ne(2).le.0) go to 20
            shift  = 1 + neq(1)*npack1
            call DMatrix_dgemm2(ne(1),ne(2),npack1_all,128,
     >             2.0d0,
     >             A1,       int_mb(ma(1,1)+taskid_i),
     >                       int_mb(ma(1,1)),
     >                       int_mb(na(1,1)),
     >             A2(shift),int_mb(ma(1,2)+taskid_i),
     >                       int_mb(ma(1,2)),
     >                       int_mb(na(1,2)),
     >             0.0d0,
     >             hml,        int_mb(md(1)+taskid_i),
     >                         int_mb(md(1)),
     >                         int_mb(nd(1)),
     >             taskid_i,taskid_j,
     >             np_i,np_j,
     >             comm_i, comm_j,
     >             dbl_mb(work1(1)),dbl_mb(work2(1)))

            call DMatrix_dgemm2(ne(1),ne(2),nida1_all,128,
     >             -1.0d0,
     >             A1,       int_mb(ma(1,1)+taskid_i),
     >                       int_mb(ma1(1,1)),
     >                       int_mb(na(1,1)),
     >             A2(shift),int_mb(ma(1,2)+taskid_i),
     >                       int_mb(ma1(1,2)),
     >                       int_mb(na(1,2)),
     >             1.0d0,
     >             hml,        int_mb(md(1)+taskid_i),
     >                         int_mb(md(1)),
     >                         int_mb(nd(1)),
     >             taskid_i,taskid_j,
     >             np_i,np_j,
     >             comm_i, comm_j,
     >             dbl_mb(work1(1)),dbl_mb(work2(1)))

  20        continue

         else

            if (ne(2).le.0) go to 21
            shift  = 1 + neq(1)*npack1
            call DMatrix_dgemm2(ne(1),ne(2),npack1_all,128,
     >             2.0d0,
     >             A1,       int_mb(ma(1,1)+taskid_i),
     >                       int_mb(ma(1,1)),
     >                       int_mb(na(1,1)),
     >             A2(shift),int_mb(ma(1,2)+taskid_i),
     >                       int_mb(ma(1,2)),
     >                       int_mb(na(1,2)),
     >             0.0d0,
     >             dbl_mb(mat_tmp(1)),int_mb(md(1)+taskid_i),
     >                                int_mb(md(1)),
     >                                int_mb(nd(1)),
     >             taskid_i,taskid_j,
     >             np_i,np_j,
     >             comm_i, comm_j,
     >             dbl_mb(work1(1)),dbl_mb(work2(1)))

            call DMatrix_dgemm2(ne(1),ne(2),nida1_all,128,
     >             -1.0d0,
     >             A1,       int_mb(ma(1,1)+taskid_i),
     >                       int_mb(ma1(1,1)),
     >                       int_mb(na(1,1)),
     >             A2(shift),int_mb(ma(1,2)+taskid_i),
     >                       int_mb(ma1(1,2)),
     >                       int_mb(na(1,2)),
     >             1.0d0,
     >             dbl_mb(mat_tmp(1)),int_mb(md(1)+taskid_i),
     >                         int_mb(md(1)),
     >                         int_mb(nd(1)),
     >             taskid_i,taskid_j,
     >             np_i,np_j,
     >             comm_i, comm_j,
     >             dbl_mb(work1(1)),dbl_mb(work2(1)))

            call Dneall_m_gather(small,smpack,int_mb(smindx(1)),
     >                        dbl_mb(mat_tmp(1)),hml)

  21        continue
         end if
         
      else

         shift  = 1+ne(1)*npack1
         if (ne(2).le.0) go to 30

         call Pack_ccmn_dot(1,ne(1),ne(2),
     >                     A1,
     >                     A2(shift),
     >                     hml)
  30        continue
      end if

      call nwpw_timing_end(15)
      return
      end


c     ****************************************
c     *                                      *
c     *        Dneall_w_determinant          *
c     *                                      *
c     ****************************************

c     ****************************************
c     *                                      *
c     *     Dneall_w_eigenvaluesvectors      *
c     *                                      *
c     ****************************************

c  This routine compute the complex eigenvalues and left and right eigenvectors of a complex matrix.

      subroutine Dneall_w_eigenvaluesvectors(mb,hml,eig,U,V)
      implicit none
      integer    mb
      complex*16 hml(*),eig(*)
      complex*16 U(*),V(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "Dne.fh"

*     ***** local variables ****
      integer MASTER
      parameter (MASTER=0)
      logical value
      integer ms,ms1,ms2,shift1,shift2,ishift1,ishift2,ierr
      integer tmp1(2),tmp2(2),tmp3(2),tmp4(2),tmp5(2)

      call nwpw_timing_start(17)

      if (mparallelized) then
         write(*,*) " WARNING : NEEDS to be fixed for distributed m"
         call errquit(
     >        'WARNING Dneall_w_eigenvaluesvectors failed',0,MA_ERR)
      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift1 = ne(1)
            ishift2 = ne(1)*ne(1)
            call Parallel_shared_vector_zero(.true.,2*(ne(1)+ne(2)),eig)
         else
            ms1 = mb
            ms2 = mb
            ishift1 = 0
            ishift2 = 0
            call Parallel_shared_vector_zero(.true.,2*ne(mb),eig)
         end if

         value = BA_push_get(mt_dcpl,(ne(1)*ne(1)),
     >                       'tmp1',tmp1(2),tmp1(1))
c         value = value.and.
c     >           BA_push_get(mt_dcpl,(ne(1)*ne(1)),
c     >                       'tmp2',tmp2(2),tmp2(1))
c         value = value.and.
c     >           BA_push_get(mt_dcpl,(ne(1)*ne(1)),
c     >                       'tmp3',tmp3(2),tmp3(1))
         value = value.and.
     >           BA_push_get(mt_dcpl,(2*ne(1)*ne(1)),
     >                       'tmp4',tmp4(2),tmp4(1))
         value = value.and.
     >           BA_push_get(mt_dbl,(2*ne(1)),
     >                       'tmp5',tmp5(2),tmp5(1))
         if (.not. value)
     >      call errquit('Dneall_w_eigenvaluesvectors:out of stack',
     >                    1,MA_ERR)

*        ***** diagonalize the matrix *****
         do ms=ms1,ms2
            shift1 = 1+(ms-1)*ishift1
            shift2 = 1+(ms-1)*ishift2
            if (ne(ms).le.0) go to 30
            call Parallel_shared_vector_copy(.true.,2*ne(ms)*ne(ms),
     >                                 hml(shift1),dcpl_mb(tmp1(1)))
!$OMP MASTER
            call ZGEEV('V','V',ne(ms),
     >                 dcpl_mb(tmp1(1)),ne(ms),
     >                 eig(shift1),
     >                 U(shift2),ne(ms),
     >                 V(shift2),ne(ms),
     >                 dcpl_mb(tmp4(1)),2*ne(1)*ne(1),
     >                 dbl_mb(tmp5(1)),
     >                 ierr)
!$OMP END MASTER
!$OMP BARRIER
  30        continue
         end do

         value =           BA_pop_stack(tmp5(2))
         value = value.and.BA_pop_stack(tmp4(2))
         !value = value.and.BA_pop_stack(tmp3(2))
         !value = value.and.BA_pop_stack(tmp2(2))
         value = value.and.BA_pop_stack(tmp1(2))
         if (.not.value)
     >    call errquit(
     >         'error popping stack in Dneall_w_eigenvaluesvectors',
     >         0,MA_ERR)

      end if

      call nwpw_timing_end(17)
      return
      end



c     ****************************************
c     *                                      *
c     *        Dneall_w_eigenvalues          *
c     *                                      *
c     ****************************************

c  This routine compute the complex eigenvalues of a complex matrix.

      subroutine Dneall_w_eigenvalues(mb,hml,eig)
      implicit none
      integer    mb
      complex*16 hml(*),eig(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "Dne.fh"

*     ***** local variables ****
      integer MASTER
      parameter (MASTER=0)
      logical value
      integer ms,ms1,ms2,shift1,shift2,ishift1,ishift2,ierr
      integer tmp1(2),tmp2(2),tmp3(2),V(2),VV(2),Q(2),tu(2),ework(2)

      call nwpw_timing_start(17)

      if (mparallelized) then
         write(*,*) " WARNING : NEEDS to be fixed for distributed m"
         call errquit('WARNING Dneall_w_eigenvalues failed',0,MA_ERR)
      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift1 = ne(1)
            ishift2 = ne(1)*ne(1)
            call dcopy(2*(ne(1)+ne(2)),0.0d0,0,eig,1)
         else
            ms1 = mb
            ms2 = mb
            ishift1 = 0
            ishift2 = 0
            call dcopy(2*ne(mb),0.0d0,0,eig,1)
         end if

         value = BA_push_get(mt_dbl,(2*ne(1)*ne(1)),
     >                       'tmp1',tmp1(2),tmp1(1))
         value = value.and.
     >                 BA_push_get(mt_dbl,(4*ne(1)*ne(1)),
     >                       'tmp2',tmp2(2),tmp2(1))
         value = value.and.
     >           BA_push_get(mt_dbl,(2*ne(1)),
     >                       'tmp3',tmp3(2),tmp3(1))
         if (.not. value)
     >      call errquit('Dneall_w_eigenvalues:out of stack',
     >                    1,MA_ERR)

*        ***** diagonalize the matrix *****
         do ms=ms1,ms2
            shift1 = 1+(ms-1)*ishift1
            shift2 = 1+(ms-1)*ishift2
            if (ne(ms).le.0) go to 30
!$OMP MASTER
            call dcopy(2*ne(ms)*ne(ms),hml(shift1),1,dbl_mb(tmp1(1)),1)
            call ZGEEV('N','N',ne(ms),
     >                 dbl_mb(tmp1(1)),ne(ms),
     >                 eig(shift1),
     >                 dbl_mb(tmp2(1)),1,dbl_mb(tmp2(1)),1,
     >                 dbl_mb(tmp2(1)),2*ne(1)*ne(1),
     >                 dbl_mb(tmp3(1)),
     >                 ierr)
!$OMP END MASTER
!$OMP BARRIER
  30        continue
         end do

         value =           BA_pop_stack(tmp3(2))
         value = value.and.BA_pop_stack(tmp2(2))
         value = value.and.BA_pop_stack(tmp1(2))
         if (.not.value)
     >    call errquit('error popping stack in Dneall_w_eigenvalues',
     >                 0,MA_ERR)

      end if

      call nwpw_timing_end(17)
      return
      end 




c     ****************************************
c     *                                      *
c     *        Dneall_m_diagonalize          *
c     *                                      *
c     ****************************************

*   This routine diagonalizes hml = v*eig*v'.  On exit hml is replaced by v.

      subroutine Dneall_m_diagonalize(mb,hml,eig,assending)
      implicit none
      integer    mb
      real*8     hml(*),eig(*)
      logical    assending

#include "bafdecls.fh"           
#include "errquit.fh"
#include "Dne.fh"

*     ***** local variables ****
      integer MASTER
      parameter (MASTER=0)
      logical value
      integer ms,ms1,ms2,shift1,shift2,ishift1,ishift2,ierr
      integer tmp1(2),V(2),VV(2),Q(2),tu(2),ework(2)



      call nwpw_timing_start(17)
      if (mparallelized) then
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift1 = ne(1)
            ishift2 = mcq(1)*ncq(1)
            call Parallel_shared_vector_zero(.true.,(ne(1)+ne(2)),eig)
         else
            ms1 = mb
            ms2 = mb
            ishift1 = 0
            ishift2 = 0
            !call dcopy(ne(mb),0.0d0,0,eig,1)
            call Parallel_shared_vector_zero(.true.,ne(mb),eig)
         end if
         value = BA_push_get(mt_dbl,mcq(1)*ncq(1),'V',V(2),V(1))
         value = value.and.
     >           BA_push_get(mt_dbl,mcq(1)*ncq(1),'VV',VV(2),VV(1))
         value = value.and.
     >           BA_push_get(mt_dbl,mcq(1)*ncq(1),'Q',Q(2),Q(1))
         value = value.and.
     >           BA_push_get(mt_dbl,ne(1),'tu',tu(2),tu(1))
         value = value.and.
     >           BA_push_get(mt_dbl,ne(1),'ework',ework(2),ework(1))
         if (.not. value) 
     >      call errquit('Dneall_m_diagonalize:out of stack',
     >                    0,MA_ERR)
         do ms=ms1,ms2
            shift1 = 1+(ms-1)*ishift1
            shift2 = 1+(ms-1)*ishift2

      
            call nwpw_timing_start(22)
            !write(*,*) "into tredq"
            call DMatrix_tredq(ne(ms),
     >                hml(shift2),dbl_mb(Q(1)),
     >                mcq(ms),int_mb(mc(1,ms)),int_mb(nc(1,ms)),
     >                taskid_i,taskid_j,
     >                np_i,np_j,
     >                comm_i,comm_j,
     >                dbl_mb(work1(1)),dbl_mb(work2(1)),
     >                dbl_mb(V(1)),dbl_mb(VV(1)))
            call nwpw_timing_end(22)

            !write(*,*) "into get_diags"
            call nwpw_timing_start(23)
            call DMatrix_getdiags(ne(ms),
     >                eig(shift1),dbl_mb(tu(1)),
     >                hml(shift2),
     >                mcq(ms),int_mb(mc(1,ms)),int_mb(nc(1,ms)),
     >                taskid_i,taskid_j,
     >                np_i,np_j,
     >                comm_i,comm_j,
     >                dbl_mb(ework(1)))
            call dcopy(mcq(ms)*ncq(ms),dbl_mb(Q(1)),1,hml(shift2),1)
            call nwpw_timing_end(23)

            !write(*,*) "into tqliq"
            call nwpw_timing_start(24)
            call DMatrix_tqliq(ne(ms),
     >                eig(shift1),dbl_mb(tu(1)),
     >                hml(shift2),
     >                mcq(ms),int_mb(mc(1,ms)),int_mb(nc(1,ms)),
     >                taskid_i,taskid_j,
     >                np_i,np_j,
     >                comm_i,comm_j,
     >                dbl_mb(work1(1)),dbl_mb(work2(1)))
            call nwpw_timing_end(24)

            call nwpw_timing_start(25)
            if (.not.assending)
     >      call DMatrix_eigsrtq(ne(ms),
     >              eig(shift1),
     >              hml(shift2),
     >              mcq(ms),int_mb(mc(1,ms)),int_mb(nc(1,ms)),
     >              taskid_i,taskid_j,
     >              np_i,np_j,
     >              comm_i,comm_j,
     >              dbl_mb(work1(1)),dbl_mb(work2(1)))
            call nwpw_timing_end(25)

         end do

  
 

         value =           BA_pop_stack(ework(2))
         value = value.and.BA_pop_stack(tu(2))
         value = value.and.BA_pop_stack(Q(2))
         value = value.and.BA_pop_stack(VV(2))
         value = value.and.BA_pop_stack(V(2))
         if (.not. value) 
     >    call errquit('error popping stack in Dneall_m_diagonalize',
     >                 0,MA_ERR)


      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift1 = ne(1)
            ishift2 = ne(1)*ne(1)
            call Parallel_shared_vector_zero(.true.,(ne(1)+ne(2)),eig)
         else
            ms1 = mb
            ms2 = mb
            ishift1 = 0
            ishift2 = 0
            !call dcopy(ne(mb),0.0d0,0,eig,1)
            call Parallel_shared_vector_zero(.true.,ne(mb),eig)
         end if

         value = BA_push_get(mt_dbl,(2*ne(1)*ne(1)),
     >                       'tmp1',tmp1(2),tmp1(1))
         value = value.and.BA_push_get(mt_dbl,ne(1),
     >                                'ework',ework(2),ework(1))
         if (.not. value) 
     >      call errquit('Dneall_m_diagonalize:out of stack',
     >                    1,MA_ERR)

*        ***** diagonalize the matrix *****
         do ms=ms1,ms2
            shift1 = 1+(ms-1)*ishift1
            shift2 = 1+(ms-1)*ishift2
            if (ne(ms).le.0) go to 30

!$OMP MASTER
            call EIGEN_SET(ne(ms),hml(shift2),dbl_mb(ework(1)))
            call DSYEV('V','U',ne(ms),
     >                 hml(shift2),ne(ms), 
     >                 eig(shift1),
     >                 dbl_mb(tmp1(1)),2*ne(1)*ne(1),
     >                 ierr)
c            ierr = 0
c            call EIGEN_UNSORT(ne(ms),ne(ms),
c     >                        hml(shift2),eig(shift1),
c     >                        dbl_mb(tmp1(1)),ierr)
            if (ierr.ne.0) then
              write(*,*) "dnem_iag:dsyev err=",ierr,0
              call EIGEN_FAIL(ne(ms),
     >                        hml(shift2),
     >                        dbl_mb(ework(1)),
     >                        eig(shift1))
              !call errquit("dnem_iag:dsyev err=",ierr,0)
            end if
!$OMP END MASTER
!$OMP BARRIER

            if (.not.assending)
     >         call EIGSRT(eig(shift1),
     >                   hml(shift2),
     >                   ne(ms),ne(ms))

  30        continue
         end do
         if (mb.eq.0) then
            call Parallel_Brdcst_values(MASTER,ne(1)*ne(1)+ne(2)*ne(2),
     >                                  hml)
         else
            call Parallel_Brdcst_values(MASTER,ne(mb)*ne(mb),hml)
         end if

         value = BA_pop_stack(ework(2))
         value = value.and.BA_pop_stack(tmp1(2))
         if (.not. value) 
     >    call errquit('error popping stack in Dneall_m_diagonalize',
     >                 0,MA_ERR)

      end if

      call nwpw_timing_end(17)
      return
      end

      subroutine EIGEN_SET(N,HML,EW)
      implicit none
      integer N
      real*8 HML(N,N),EW(N)
      integer i
      do i=1,N
         EW(I) = HML(i,i)
      end do
      return
      end

      subroutine EIGEN_FAIL(N,HML,EW,EIG)
      implicit none
      integer N
      real*8 HML(N,N),EW(N),EIG(N)
      integer i
      do i=1,N
         EIG(i) = EW(i)
      end do
      call dcopy(N*N,0.0d0,0,HML,1)
      do i=1,N
          HML(i,i) = 1.0d0
      end do
      return
      end

c     ****************************************
c     *                                      *
c     *        Dneall_m_cholesky             *
c     *                                      *
c     ****************************************

*   This routine computes the cholesky decomposition of a postitive definate
*  matrix A = L*L'   On exit A is replaced by L.
      subroutine Dneall_m_cholesky(mb,A)
      implicit none
      integer mb
      real*8 A(*)

#include "bafdecls.fh"           
#include "errquit.fh"
#include "Dne.fh"

*     **** local variables ****
      integer i,j,ms,ms1,ms2,nn,ishift2,shift2,ierr

      if (mb.eq.0) then
         ms1 = 1
         ms2 = ispin
         ishift2 = ne(1)*ne(1)
      else
         ms1 = mb
         ms2 = mb
         ishift2 = 0
      end if

*     ***** WARNING : NEEDS to fixed for distributed m ****
      if (mparallelized) then
         write(*,*) " WARNING : NEEDS to be fixed for distributed m"
         call errquit('WARNING Dneall_m_cholesky failed',0,MA_ERR)
      else
         do ms=ms1,ms2
            shift2 = (ms-1)*ishift2+1
            call DPOTRF('L',ne(ms),A(shift2),ne(ms),ierr)
         end do
      end if

      return
      end




c     ****************************************
c     *                                      *
c     *        Dneall_mg_forwardsolve        *
c     *                                      *
c     ****************************************

*   This routine computes transpose(B) =  inv(L)*transpose(A)
*  On exit A is replaced by B.
      subroutine Dneall_mg_forwardsolve(mb,L,n2ft3d,A)
      implicit none
      integer mb,n2ft3d
      real*8 L(*),A(n2ft3d,*)

#include "bafdecls.fh"           
#include "errquit.fh"
#include "Dne.fh"

*     **** local variables ****
      !real*8 tmp
      integer i,j,k,ms,ms1,ms2,shift1,shift2,ishift1,ishift2,indx,indx0
      integer p,q,ii,jj
      integer tmp(2)

      if (mb.eq.0) then
         ms1 = 1
         ms2 = ispin
         ishift1 = ne(1)
         ishift2 = ne(1)*ne(1)
      else
         ms1 = mb
         ms2 = mb
         ishift1 = 0
         ishift2 = 0
      end if

      if (mparallelized) then
         write(*,*) " WARNING : NEEDS to fixed for distributed m"
         call errquit('WARNING Dneall_mg_forwardsolve failed',0,MA_ERR)
      else

         if (.not.BA_push_get(mt_dbl,n2ft3d,'tmp',tmp(2),tmp(1)))
     >      call errquit('Dneall_mg_forwardsolve:out of stack',0,MA_ERR)

         do ms=ms1,ms2
            shift1 = (ms-1)*ishift1
            shift2 = (ms-1)*ishift2
            do i=1,ne(ms)
               indx0 = i+(i-1)*ne(ms)+shift2
               !tmp = 0.0d0
               !call dcopy(n2ft3d,0.0d0,0,dbl_mb(tmp(1)),1)
               call Parallel_shared_vector_zero(.true.,n2ft3d,
     >                                          dbl_mb(tmp(1)))
               do j=1,i-1
                  jj = j + shift1
                  call Dneall_ntoqp(jj,q,p)
                  if (p.eq.taskid_j) then
                     indx = i+(j-1)*ne(ms)+shift2
                     !tmp = tmp + L(indx)*A(k,q+shift1)
                     !write(*,*) "i,j,L=",i,j,L(indx)
                     call daxpy_omp(n2ft3d,L(indx),A(1,q),1,
     >                              dbl_mb(tmp(1)),1)
                  end if
               end do !*j*
               !call D1dB_SumAll(tmp)
               !A(k,j) = (A(k,j)-tmp)/L(indx0)
               !write(*,*) "i,L(i,i)=",i,L(indx0)
               call D1dB_Vector_SumAll(n2ft3d,dbl_mb(tmp(1)))
               ii = i + shift1
               call Dneall_ntoqp(ii,q,p)
               if (p.eq.taskid_j) then
                 call daxpy_omp(n2ft3d,-1.0d0,dbl_mb(tmp(1)),1,A(1,q),1)
                 call dscal_omp(n2ft3d,(1.0d0/L(indx0)),A(1,q),1)
               end if
            end do !*i*
         end do !*ms*

         if (.not.BA_pop_stack(tmp(2)))
     >      call errquit('Dneall_mg_forwardsolve: popping stack',
     >                   0,MA_ERR)

      end if

      return
      end 


c     ****************************************
c     *                                      *
c     *        Dneall_m_SVD                 *
c     *                                      *
c     ****************************************

*
*    Computes the SVD decomposition of A = U*S*V', where A is an neall x neall matrix.

      subroutine Dneall_m_SVD(mb,A,U,S,V)
      implicit none
      integer    mb
      real*8 A(*),U(*)
      real*8 S(*),V(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "Dne.fh"

*     **** local variables ****
      integer ms,ms1,ms2,n,shift,tmp2(2)
      integer nj
      real*8 sum

*     ***** external functions ****
      real*8   ddot
      external ddot

      if (mb.eq.0) then
         ms1 = 1
         ms2 = ispin
         nj  = ne(1)+ne(2)
      else
         ms1 = mb
         ms2 = mb
         nj  = ne(mb)
      end if

*     **** generate V and Sigma^2 ****
      call Dneall_mmm_Multiply2(mb,A,A,V)
      call Dneall_m_diagonalize(mb,V,S,.false.)

*     **** generate U*Sigma ****
      call Dneall_mmm_Multiply(mb,A,V,1.0d0,U,0.0d0)

*     **** normalize U*sigma ****
      do ms=1,ispin
         shift = 1+(ms-1)*ne(1)*ne(1)
         do n=1,ne(ms)
            sum = ddot(ne(ms),U(shift),1,U(shift),1)
            sum = 1.0d0/dsqrt(sum)
            call dscal(ne(ms),sum,U(shift),1)
            shift = shift + ne(ms)
         end do
      end do

*     **** calculated sqrt(S^2) ****
      do n=1,nj
         if (S(n).lt.0.0d0) S(n) = dabs(S(n))
         S(n) = dsqrt(S(n))
      end do

      return
      end







c     ****************************************
c     *                                      *
c     *        Dneall_f_SVD                  *
c     *                                      *
c     ****************************************

*   This routine computes the SVD decomposition of A = U*S*V'
*  where A has a dimension of A(npack1,neall)

      subroutine Dneall_f_SVD(mb,A,U,npack1,S,V)
      implicit none
      integer    mb
      complex*16 A(*),U(*)
      integer    npack1      
      real*8     S(*),V(*)

#include "bafdecls.fh"           
#include "errquit.fh"
#include "Dne.fh"

*     **** local variables ****
      integer ms,ms1,ms2,n,shift,tmp2(2)
      integer nj   

      if (mb.eq.0) then
         ms1 = 1
         ms2 = ispin
         nj  = ne(1)+ne(2)
      else
         ms1 = mb
         ms2 = mb
         nj  = ne(mb)
      end if

*     **** generate V and Sigma^2 ****
      call Dneall_ffm_sym_Multiply(mb,A,A,npack1,V)
      call Dneall_m_diagonalize(mb,V,S,.false.)

*     **** generate U*Sigma ****
      call Dneall_fmf_Multiply(mb,A,npack1,V,1.0d0,U,0.0d0)

*     **** normalize U*sigma ****
      if (.not. BA_push_get(mt_dbl,neq(1),'tmp2',tmp2(2),tmp2(1)))
     >   call errquit('Dneall_f_SVD:out of stack memory',0,MA_ERR)

      do ms=ms1,ms2
         if (neq(ms).le.0) go to 30
         shift = 1+(ms-1)*neq(1)*npack1
         do n=1,neq(ms)
            call Pack_cc_idot(1,U(shift),U(shift),dbl_mb(tmp2(1)+n-1))
            shift = shift + npack1
         end do
         call D3dB_Vector_SumAll(neq(ms),dbl_mb(tmp2(1)))

!$OMP DO
         do n=1,neq(ms)
            dbl_mb(tmp2(1)+n-1) = 1.0d0/dsqrt(dbl_mb(tmp2(1)+n-1))
         end do
!$OMP END DO

         shift = 1+(ms-1)*neq(1)*npack1
         do n=1,neq(ms)
            call DSCAL_OMP(2*npack1,dbl_mb(tmp2(1)+n-1),U(shift),1)
            shift = shift + npack1
         end do

 30     continue
      end do

      if (.not.BA_pop_stack(tmp2(2))) 
     >  call errquit('Dneall_f_SVD: popping stack memory',0,MA_ERR)


*     **** calculated sqrt(S^2) ****
!$OMP DO
      do n=1,nj
         if (S(n).lt.0.0d0) S(n) = dabs(S(n))
         S(n) = dsqrt(S(n))
      end do
!$OMP END DO

      return
      end


c     ****************************************
c     *                                      *
c     *        Dneall_f_SVD_ASA              *
c     *                                      *
c     ****************************************

*   This routine computes the SVD decomposition of A = U*Sigma*V', such that U'SU = I
*  where A has a dimension of A(npack1,neall)
*
*

      subroutine Dneall_f_SVD_ASA1(mb,A,SA,USigma,npack1,Sigma,V)
      implicit none
      integer    mb
      complex*16 A(*),SA(*),USigma(*)
      integer    npack1
      real*8     Sigma(*),V(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "Dne.fh"

*     **** local variables ****
      integer ms,ms1,ms2,n,shift,tmp2(2)
      integer nj

      if (mb.eq.0) then
         ms1 = 1
         ms2 = ispin
         nj  = ne(1)+ne(2)
      else
         ms1 = mb
         ms2 = mb
         nj  = ne(mb)
      end if

*     **** generate V and Sigma^2 ****
      call Dneall_ffm_sym_Multiply(mb,A,SA,npack1,V)
      call Dneall_m_diagonalize(mb,V,Sigma,.false.)

*     **** generate U*Sigma ****
      call Dneall_fmf_Multiply(mb,A,npack1,V,1.0d0,USigma,0.0d0)


*     **** calculated sqrt(S^2) ****
      do n=1,nj
         if (Sigma(n).lt.0.0d0) Sigma(n) = dabs(Sigma(n))
         Sigma(n) = dsqrt(Sigma(n))
      end do

      return
      end

      subroutine Dneall_f_SVD_ASA2(mb,U,SU,npack1)
      implicit none
      integer    mb
      complex*16 U(*),SU(*)
      integer    npack1

#include "bafdecls.fh"
#include "errquit.fh"
#include "Dne.fh"

*     **** local variables ****
      integer ms,ms1,ms2,n,shift,tmp2(2)
      integer nj

      if (mb.eq.0) then
         ms1 = 1
         ms2 = ispin
         nj  = ne(1)+ne(2)
      else
         ms1 = mb
         ms2 = mb
         nj  = ne(mb)
      end if

*     **** normalize U*sigma ****
      if (.not. BA_push_get(mt_dbl,neq(1),'tmp2',tmp2(2),tmp2(1)))
     >   call errquit('Dneall_f_SVD_ASA2:out of stack memory',0,MA_ERR)
      call dcopy(neq(1),0.0d0,0,dbl_mb(tmp2(1)),1)

      do ms=ms1,ms2
         if (neq(ms).le.0) go to 30
         shift = 1+(ms-1)*neq(1)*npack1
         do n=1,neq(ms)
            call Pack_cc_idot(1,U(shift),SU(shift),dbl_mb(tmp2(1)+n-1))
            shift = shift + npack1
         end do
         call D3dB_Vector_SumAll(neq(ms),dbl_mb(tmp2(1)))


         do n=1,neq(ms)
            dbl_mb(tmp2(1)+n-1) = 1.0d0/dsqrt(dbl_mb(tmp2(1)+n-1))
         end do

         shift = 1+(ms-1)*neq(1)*npack1
         do n=1,neq(ms)
            call dscal(2*npack1,dbl_mb(tmp2(1)+n-1),U(shift),1)
            shift = shift + npack1
         end do

 30     continue
      end do

      if (.not.BA_pop_stack(tmp2(2)))
     >  call errquit('Dneall_f_SVD_ASA2: popping stack memory',0,MA_ERR)

      return
      end






c     ****************************************
c     *                                     *
c     *        Dneall_f_ortho               *
c     *                                      *
c     ****************************************

      subroutine Dneall_f_ortho(mb,U,npack1)
      integer mb
      complex*16 U(*)
      integer    npack1      

#include "bafdecls.fh"           
#include "errquit.fh"
#include "Dne.fh"

*     **** local variables ****
      integer taskid
      integer ms,ms1,ms2,n,shift,asize,j,k,indxk,indxj
      integer V(2),tmp2(2),S(2),A(2)
      real*8 sum1

      logical  Dneall_m_allocate,Dneall_m_free
      external Dneall_m_allocate,Dneall_m_free
 

      
      if (parallelized) then
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            asize = (neq(1)+neq(2))*npack1
         else
            ms1 = mb
            ms2 = mb
            asize = neq(mb)*npack1
         end if

      
         if (.not.BA_push_get(mt_dcpl,asize,'A',A(2),A(1)))
     >   call errquit('Dneall_f_ortho:out of stack memory',0,MA_ERR)
         if (.not.Dneall_m_allocate(mb,V))
     >   call errquit('Dneall_f_ortho:out of stack memory',1,MA_ERR)
         if (.not.BA_push_get(mt_dbl,(ne(1)+ne(2)),'S',S(2),S(1)))
     >   call errquit('Dneall_f_ortho:out of stack memory',2,MA_ERR)


*        **** generate V and Sigma^2 ****
         call Dneall_ffm_sym_Multiply(mb,U,U,npack1,dbl_mb(V(1)))
         call Dneall_m_diagonalize(mb,dbl_mb(V(1)),dbl_mb(S(1)),.false.)
         if (.not.BA_pop_stack(S(2)))
     >   call errquit('Dneall_f_ortho: popping stack memory',2,MA_ERR)

*        **** generate U*Sigma ****
         call dcopy(2*asize,U,1,dcpl_mb(A(1)),1)
         call Dneall_fmf_Multiply(mb,dcpl_mb(A(1)),npack1,
     >                         dbl_mb(V(1)),1.0d0,U,0.0d0)

*        **** deallocate tmp space ****
         if (.not.Dneall_m_free(V))
     >     call errquit('Dneall_f_ortho: popping stack memory',0,MA_ERR)
         if (.not.Ma_pop_stack(A(2)))
     >     call errquit('Dneall_f_ortho: popping stack memory',1,MA_ERR)


*        **** normalize U*sigma ****
         if (.not. BA_push_get(mt_dbl,neq(1),'tmp2',tmp2(2),tmp2(1)))
     >      call errquit('Dneall_f_ortho:out of stack memory',3,MA_ERR)

         do ms=ms1,ms2
            if (neq(ms).le.0) go to 30
            shift = 1+(ms-1)*neq(1)*npack1
            do n=1,neq(ms)
              call Pack_cc_idot(1,U(shift),U(shift),dbl_mb(tmp2(1)+n-1))
              shift = shift + npack1
            end do
            call D3dB_Vector_SumAll(neq(ms),dbl_mb(tmp2(1)))

            do n=1,neq(ms)
               dbl_mb(tmp2(1)+n-1) = 1.0d0/dsqrt(dbl_mb(tmp2(1)+n-1))
            end do

            shift = 1+(ms-1)*neq(1)*npack1
            do n=1,neq(ms)
               call dscal(2*npack1,dbl_mb(tmp2(1)+n-1),U(shift),1)
               shift = shift + npack1
            end do

 30        continue
         end do
         if (.not.BA_pop_stack(tmp2(2)))
     >     call errquit('Dneall_f_ortho: popping stack memory',3,MA_ERR)



c     **** not parallized ****
      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
         else
            ms1 = mb
            ms2 = mb
         end if
         shift = ne(1)*npack1

         !**** orthogonalize from the bottom -> up ****
         do ms=ms1,ms2
         do k=ne(ms),1,-1
            indxk = 1+(k-1)*npack1 + (ms-1)*shift
            call Pack_cc_dot(1,U(indxk),U(indxk),sum1)
            sum1 = 1.0d0/dsqrt(sum1)
c            call Pack_c_SMul(1,sum1,U(indxk),U(indxk))
            call Pack_c_SMul1(1,sum1,U(indxk))

            do j=k-1,1,-1
               indxj = 1+(j-1)*npack1 + (ms-1)*shift
               call Pack_cc_dot(1,U(indxk),U(indxj),sum1)
               sum1 = -sum1
               call Pack_cc_daxpy(1,sum1,U(indxk),U(indxj))
            end do
         end do
         end do


      end if

      return
      end 


c     ****************************************
c     *                                      *
c     *        Dneall_mm_transpose           *
c     *                                      *
c     ****************************************

      subroutine Dneall_mm_transpose(mb,M0,M1)
      implicit none
      integer mb
      real*8  M0(*),M1(*)


#include "bafdecls.fh"
#include "Dne.fh"

*     **** local variables ****
      integer ms,ms1,ms2,shift2,ishift2

      if (mparallelized) then
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = mcq(1)*ncq(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if

         do ms=ms1,ms2
           if (ne(ms).le.0) goto 20
           shift2 = 1 + (ms-1)*ishift2
           call DMatrix_mm_transpose(ne(ms),
     >                 M0(shift2),M1(shift2),
     >                 mcq(ms),int_mb(mc(1,ms)),int_mb(nc(1,ms)))
 20        continue
         end do

      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = ne(1)*ne(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if

         do ms=ms1,ms2
           if (ne(ms).le.0) goto 30
           shift2 = 1 + (ms-1)*ishift2
           call Dneall_mm_transpose_sub(ne(ms),
     >                 M0(shift2),M1(shift2))
 30        continue
         end do
      end if

      return
      end

      subroutine Dneall_mm_transpose_sub(n,A,B)
      implicit none
      integer n
      real*8 A(n,n)
      real*8 B(n,n)
      integer i,j
!$OMP DO private(i,j)
      do j=1,n
         do i=1,n
            B(i,j) = A(j,i)
         end do
      end do
!$OMP END DO
      return
      end


c     ****************************************
c     *                                      *
c     *        Dneall_mm_SCtimesVtrans       *
c     *                                      *
c     ****************************************

      subroutine Dneall_mm_SCtimesVtrans(mb,t,S,Vt,A,B,SA,SB)
      implicit none
      integer mb
      real*8 t
      real*8 S(*)
      real*8 Vt(*)
      real*8 A(*),B(*)
      real*8 SA(*),SB(*)

#include "bafdecls.fh"
#include "Dne.fh"

*     **** local variables ****
      integer ms,ms1,ms2,shift1,shift2,ishift1,ishift2,k,j,indx1,indx2
      integer nj

      if (mparallelized) then
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift1 = ne(1)
            ishift2 = mcq(1)*ncq(1)
            nj = ne(1)+ne(2)
         else
            ms1 = mb
            ms2 = mb
            ishift1 = 0
            ishift2 = 0
            nj = ne(mb)
         end if

!$OMP DO
         do j=1,nj
            SA(j) = dcos(S(j)*t)
            SB(j) = dsin(S(j)*t)
         end do
!$OMP END DO

         do ms=ms1,ms2
             shift1 = 1 + (ms-1)*ishift1
             shift2 = 1 + (ms-1)*ishift2
             call DMatrix_SASBtimesVtrans(ne(ms),
     >                SA(shift1),SB(shift1),
     >                Vt(shift2),A(shift2),B(shift2),
     >                mcq(ms),int_mb(mc(1,ms)),int_mb(nc(1,ms)),
     >                taskid_i,taskid_j)
         end do
      
      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift1 = ne(1)
            ishift2 = ne(1)*ne(1)
            nj = ne(1)+ne(2)
         else
            ms1 = mb
            ms2 = mb
            ishift1 = 0
            ishift2 = 0
            nj = ne(mb)
         end if

!$OMP DO private(j)
         do j=1,nj
            SA(j) = dcos(S(j)*t)
            SB(j) = dsin(S(j)*t)
         end do
!$OMP END DO

         do ms=ms1,ms2
           if (ne(ms).le.0) go to 30
           shift1 = 1 + (ms-1)*ishift1
           shift2 = 1 + (ms-1)*ishift2

           !indx2 = shift2
!$OMP DO private(j,k,indx1,indx2)
           do k=1,ne(ms)
              indx1 = shift1
              indx2 = shift2 + (k-1)*ne(ms)
              do j=1,ne(ms)
                 A(indx2) = SA(indx1)*Vt(indx2)
                 B(indx2) = SB(indx1)*Vt(indx2)
                 indx2 = indx2 + 1
                 indx1 = indx1 + 1
              end do
           end do
!$OMP END DO
   
 30        continue
         end do
      end if

      return
      end



c     ****************************************
c     *                                      *
c     *        Dneall_mm_SCtimesVtrans2      *
c     *                                      *
c     ****************************************

      subroutine Dneall_mm_SCtimesVtrans2(mb,t,S,Vt,A,B,SA,SB)
      implicit none
      integer mb
      real*8 t
      real*8 S(*)
      real*8 Vt(*)
      real*8 A(*),B(*)
      real*8 SA(*),SB(*)

#include "bafdecls.fh"
#include "Dne.fh"

*     **** local variables ****
      integer ms,ms1,ms2,shift1,shift2,ishift1,ishift2,k,j,indx1,indx2
      integer nj
 
      if (mparallelized) then
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift1 = ne(1)
            ishift2 = mcq(1)*ncq(1)
            nj = ne(1)+ne(2)
         else
            ms1 = mb
            ms2 = mb
            ishift1 = 0
            ishift2 = 0
            nj = ne(mb)
         end if

!$OMP DO
         do j=1,nj
            SA(j) = S(j)*dsin(S(j)*t)
            SB(j) = S(j)*dcos(S(j)*t)
         end do
!$OMP END DO

         do ms=ms1,ms2
             shift1 = 1 + (ms-1)*ishift1
             shift2 = 1 + (ms-1)*ishift2
             call DMatrix_SASBtimesVtrans(ne(ms),
     >                SA(shift1),SB(shift1),
     >                Vt(shift2),A(shift2),B(shift2),
     >                mcq(ms),int_mb(mc(1,ms)),int_mb(nc(1,ms)),
     >                taskid_i,taskid_j)
         end do


      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift1 = ne(1)
            ishift2 = ne(1)*ne(1)
            nj = ne(1)+ne(2)
         else
            ms1 = mb
            ms2 = mb
            ishift1 = 0
            ishift2 = 0
            nj = ne(mb)
         end if

!$OMP DO private(j)
         do j=1,nj
            SA(j) = S(j)*dsin(S(j)*t)
            SB(j) = S(j)*dcos(S(j)*t)
         end do
!$OMP END DO

         do ms=ms1,ms2
           if (ne(ms).le.0) go to 30

           shift1 = 1 + (ms-1)*ishift1
           shift2 = 1 + (ms-1)*ishift2

           indx2 = shift2
!$OMP DO private(j,k,indx1,indx2)
           do k=1,ne(ms)
              indx1 = shift1
              indx2 = shift2 + (k-1)*ne(ms)
              do j=1,ne(ms)
                 A(indx2) = SA(indx1)*Vt(indx2)
                 B(indx2) = SB(indx1)*Vt(indx2)
                 indx1 = indx1 + 1
                 indx2 = indx2 + 1
              end do
           end do
!$OMP END DO

 30        continue
         end do
      end if

      return
      end



c     ****************************************
c     *                                      *
c     *        Dneall_mm_SCtimesVtrans3      *
c     *                                      *
c     ****************************************

      subroutine Dneall_mm_SCtimesVtrans3(mb,t,S,Vt,A,B,SA,SB)
      implicit none
      integer mb
      real*8 t
      real*8 S(*)
      real*8 Vt(*)
      real*8 A(*),B(*)
      real*8 SA(*),SB(*)

#include "bafdecls.fh"
#include "Dne.fh"

*     **** local variables ****
      integer ms,ms1,ms2,shift1,shift2,ishift1,ishift2,k,j,indx1,indx2
      integer nj

      if (mparallelized) then
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift1 = ne(1)
            ishift2 = mcq(1)*ncq(1)
            nj = ne(1)+ne(2)
         else
            ms1 = mb
            ms2 = mb
            ishift1 = 0
            ishift2 = 0
            nj = ne(mb)
         end if

!$OMP DO private(j)
         do j=1,nj
            SA(j) = dsin(S(j)*t)
            SB(j) = 1.0d0-dcos(S(j)*t)
         end do
!$OMP END DO

         do ms=ms1,ms2
             shift1 = 1 + (ms-1)*ishift1
             shift2 = 1 + (ms-1)*ishift2
             call DMatrix_SASBtimesVtrans(ne(ms),
     >                SA(shift1),SB(shift1),
     >                Vt(shift2),A(shift2),B(shift2),
     >                mcq(ms),int_mb(mc(1,ms)),int_mb(nc(1,ms)),
     >                taskid_i,taskid_j)
         end do
         
      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift1 = ne(1)
            ishift2 = ne(1)*ne(1)
            nj = ne(1)+ne(2)
         else
            ms1 = mb
            ms2 = mb
            ishift1 = 0
            ishift2 = 0
            nj = ne(mb)
         end if

!$OMP DO private(j)
         do j=1,nj
            SA(j) = dsin(S(j)*t)
            SB(j) = 1.0d0-dcos(S(j)*t)
         end do
!$OMP END DO 

         do ms=ms1,ms2
           if (ne(ms).le.0) go to 30
           shift1 = 1 + (ms-1)*ishift1
           shift2 = 1 + (ms-1)*ishift2

           indx2 = shift2
!$OMP DO private(j,k,indx1,indx2)
           do k=1,ne(ms)
              indx1 = shift1
              indx2 = shift2 + (k-1)*ne(ms)
              do j=1,ne(ms)
                 A(indx2) = SA(indx1)*Vt(indx2)
                 B(indx2) = SB(indx1)*Vt(indx2)
                 indx1 = indx1 + 1
                 indx2 = indx2 + 1
              end do
           end do
!$OMP END DO

 30        continue
         end do
      end if

      return
      end


c     ****************************************
c     *                                      *
c     *        Dneall_mmm_Multiply2         *
c     *                                      *
c     ****************************************

      subroutine Dneall_mmm_Multiply2(mb,A,B,C)
      implicit none
      integer mb
      real*8 A(*),B(*),C(*)
           
#include "bafdecls.fh"
#include "Dne.fh"

*     **** local variables ****
      integer MASTER
      parameter (MASTER=0)
      integer ms,ms1,ms2,n,shift2,ishift2
  
      call nwpw_timing_start(18)
      if (mparallelized) then
        if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = mcq(1)*ncq(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if

         do ms=ms1,ms2
            if (ne(ms).le.0) go to 20
            shift2 = 1 + (ms-1)*ishift2
            call DMatrix_dgemm2(ne(ms),ne(ms),ne(ms),64,
     >             1.0d0,
     >             A(shift2),mcq(ms),
     >                       int_mb(mc(1,ms)),
     >                       int_mb(nc(1,ms)),
     >             B(shift2),mcq(ms),
     >                       int_mb(mc(1,ms)),
     >                       int_mb(nc(1,ms)),
     >             0.0d0,
     >             C(shift2),mcq(ms),
     >                       int_mb(mc(1,ms)),
     >                       int_mb(nc(1,ms)),
     >             taskid_i,taskid_j,
     >             np_i,np_j,
     >             comm_i, comm_j,
     >             dbl_mb(work1(1)),dbl_mb(work2(1)))

   20       continue
         end do


      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = ne(1)*ne(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if

         do ms=ms1,ms2
            n     = ne(ms)
            if (n.le.0) go to 30
            shift2 = 1 + (ms-1)*ishift2
            call DGEMM_OMP('T','N',n,n,n,1.0d0,
     >                A(shift2), n,
     >                B(shift2), n,
     >                0.0d0,
     >                C(shift2), n)
   30       continue
         end do
         if (mb.eq.0) then
            call Parallel_Brdcst_values(MASTER,
     >                   ne(1)*ne(1)+ne(2)*ne(2),C)
         else
            call Parallel_Brdcst_values(MASTER,ne(mb)*ne(mb),C)
         end if

      end if

      call nwpw_timing_end(18)
      return
      end


c     ****************************************
c     *                                      *
c     *        Dneall_mmm_Multiply3          *
c     *                                      *
c     ****************************************

      subroutine Dneall_mmm_Multiply3(mb,A,B,C)
      implicit none
      integer mb
      real*8 A(*),B(*),C(*)

#include "bafdecls.fh"
#include "Dne.fh"

*     **** local variables ****
      integer MASTER
      parameter (MASTER=0)
      integer ms,ms1,ms2,n,shift2,ishift2

      call nwpw_timing_start(18)
      if (mparallelized) then
        if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = mcq(1)*ncq(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if

         do ms=ms1,ms2
            if (ne(ms).le.0) go to 20
            shift2 = 1 + (ms-1)*ishift2
            call DMatrix_dgemm3(ne(ms),ne(ms),ne(ms),64,
     >             1.0d0,
     >             A(shift2),mcq(ms),
     >                       int_mb(mc(1,ms)),
     >                       int_mb(nc(1,ms)),
     >             B(shift2),mcq(ms),
     >                       int_mb(mc(1,ms)),
     >                       int_mb(nc(1,ms)),
     >             0.0d0,
     >             C(shift2),mcq(ms),
     >                       int_mb(mc(1,ms)),
     >                       int_mb(nc(1,ms)),
     >             taskid_i,taskid_j,
     >             np_i,np_j,
     >             comm_i, comm_j,
     >             dbl_mb(work1(1)),dbl_mb(work2(1)))

   20       continue
         end do

      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = ne(1)*ne(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if

         do ms=ms1,ms2
            n     = ne(ms)
            if (n.le.0) go to 30
            shift2 = 1 + (ms-1)*ishift2
            call DGEMM_OMP('N','T',n,n,n,1.0d0,
     >                A(shift2), n,
     >                B(shift2), n,
     >                0.0d0,
     >                C(shift2), n)
   30       continue
         end do
         if (mb.eq.0) then
            call Parallel_Brdcst_values(MASTER,
     >                   ne(1)*ne(1)+ne(2)*ne(2),C)
         else
            call Parallel_Brdcst_values(MASTER,ne(mb)*ne(mb),C)
         end if

      end if 
      call nwpw_timing_end(18)

      return
      end





c     ****************************************
c     *                                      *
c     *        Dneall_mmm_Multiply           *
c     *                                      *
c     ****************************************

      subroutine Dneall_mmm_Multiply(mb,A,B,alpha,C,beta)
      implicit none
      integer mb
      real*8 A(*),B(*),C(*)
      real*8 alpha,beta
           
#include "bafdecls.fh"
#include "Dne.fh"

*     **** local variables ****
      integer MASTER
      parameter (MASTER=0)
      integer ms,ms1,ms2,n,shift2,ishift2
  
      call nwpw_timing_start(18)
      if (mparallelized) then
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = mcq(1)*ncq(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if

         do ms=ms1,ms2
            if (ne(ms).le.0) go to 20
            shift2 = 1 + (ms-1)*ishift2
            call DMatrix_dgemm1(ne(ms),ne(ms),ne(ms),64,
     >             alpha,
     >             A(shift2),mcq(ms),
     >                       int_mb(mc(1,ms)),
     >                       int_mb(nc(1,ms)),
     >             B(shift2),mcq(ms),
     >                       int_mb(mc(1,ms)),
     >                       int_mb(nc(1,ms)),
     >             beta,
     >             C(shift2),mcq(ms),
     >                       int_mb(mc(1,ms)),
     >                       int_mb(nc(1,ms)),
     >             taskid_i,taskid_j,
     >             np_i,np_j,
     >             comm_i, comm_j,
     >             dbl_mb(work1(1)),dbl_mb(work2(1)))

   20       continue
         end do

      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = ne(1)*ne(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if

         do ms=ms1,ms2
            n     = ne(ms)
            if (n.le.0) go to 30
            shift2 = 1 + (ms-1)*ishift2
            call DGEMM_OMP('N','N',n,n,n,
     >                alpha,
     >                A(shift2), n,
     >                B(shift2), n,
     >                beta,
     >                C(shift2), n)
   30       continue
         end do
         if (mb.eq.0) then
            call Parallel_Brdcst_values(MASTER,
     >                   ne(1)*ne(1)+ne(2)*ne(2),C)
         else
            call Parallel_Brdcst_values(MASTER,ne(mb)*ne(mb),C)
         end if

      end if
      call nwpw_timing_end(18)

      return
      end

c     ****************************************
c     *                                      *
c     *        Dneall_mmm_Multiply0          *
c     *                                      *
c     ****************************************

      subroutine Dneall_mmm_Multiply0(mb,A,B,alpha,C,beta)
      implicit none
      integer mb
      real*8 A(*),B(*),C(*)
      real*8 alpha,beta

#include "bafdecls.fh"
#include "Dne.fh"

*     **** local variables ****
      integer MASTER
      parameter (MASTER=0)
      integer ms,ms1,ms2,n,shift2,ishift2

      call nwpw_timing_start(18)
      if (mparallelized) then
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = mcq(1)*ncq(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if

         do ms=ms1,ms2
            if (ne(ms).le.0) go to 20
            shift2 = 1 + (ms-1)*ishift2
            call DMatrix_dgemm1(ne(ms),ne(ms),ne(ms),64,
     >             alpha,
     >             A(shift2),mcq(ms),
     >                       int_mb(mc(1,ms)),
     >                       int_mb(nc(1,ms)),
     >             B(shift2),mcq(ms),
     >                       int_mb(mc(1,ms)),
     >                       int_mb(nc(1,ms)),
     >             beta,
     >             C(shift2),mcq(ms),
     >                       int_mb(mc(1,ms)),
     >                       int_mb(nc(1,ms)),
     >             taskid_i,taskid_j,
     >             np_i,np_j,
     >             comm_i, comm_j,
     >             dbl_mb(work1(1)),dbl_mb(work2(1)))

   20       continue
         end do

      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = ne(1)*ne(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if

         do ms=ms1,ms2
            n     = ne(ms)
            if (n.le.0) go to 30
            shift2 = 1 + (ms-1)*ishift2
            call DGEMM_OMP('N','N',n,n,n,
     >                alpha,
     >                A(shift2), n,
     >                B(shift2), n,
     >                beta,
     >                C(shift2), n)
   30       continue
         end do
      end if
      call nwpw_timing_end(18)

      return
      end

c      subroutine Dneall_mmm_Multiply_sub0(N,Nb,A,B,alpha,C,beta,Ctmp)
c      implicit none
c      integer N,Nb
c      real*8 A(N,N),B(N,N),alpha,
c      real*8 C(N,N),beta
c      real*8 Ctmp(N,N)
c
c*     **** local variables ****
c      integer  s,r,i,j,k,shifti,shiftj,shiftk
c      integer  taskid,np,ic,taskcount
c
c*     **** external functions ****
c      integer  Parallel2d_taskid_i,Parallel2d_taskid_j
c      external Parallel2d_taskid_i,Parallel2d_taskid_j
c      integer  Parallel2d_np_i,Parallel2d_np_j
c      external Parallel2d_np_i,Parallel2d_np_j
c
c      np_i = Parallel2d_np_i()
c      np_j = Parallel2d_np_j()
c      Nb = np_i*np_j
c      if (np_i*np_j).gt.
c
c      taskid_i = Parallel_taskid_i()
c      taskid_j = Parallel_taskid_i()
c      taskcount = 0
c      call Parallel_shared_vector_zero(.true.,N*N,Ctmp)
c      s = N/Nb
c      r = mod(N,Nb)
c      do k=1,Nb
c         shiftk = 1 + (k-1)*s + min(k,r)
c         nk = s
c         if (k.lt.r) nk = nk + 1
c
c         do j=1,Nb
c            shiftj = 1 + (j-1)*s + min(j,r)
c            nj = s
c            if (j.lt.r) nj = nj + 1
c
c            do i=1,Nb
c               shifti = 1 + (i-1)*s + min(i,r)
c               ni = s
c               if (i.lt.r) ni = ni + 1
c         
c               if (taskid.eq.taskcount) 
c     >             call DGEMM_OMP0('N','N',ni,nj,nk,
c     >                alpha,
c     >                A(shifti,shiftk), N,
c     >                B(shiftk,shiftj), N,
c     >                0.0d0,
c     >                Ctmp(shifti,shiftj), N)
c               taskcount = taskcount + 1
c            end do
c         end do
c!$OMP BARRIER
c      end do
c
c      call Parallela_Vector_SumAll(ic,N*N,Ctmp)
c
c      call dscal_omp(N*N,beta,C)
c      call daxpy_omp(N*N,1.0d0,Ctmp,1,C)
c
c
c      return
c      end 



c     ****************************************
c     *                                      *
c     *        Dneall_mmm_replicated         *
c     *                                      *
c     ****************************************
      subroutine  Dneall_mmm_replicated1(n,alpha,a,b,beta,c,ctmp)
      implicit none
      integer n
      real*8 alpha,beta
      real*8 a(n,n),b(n,n),c(n,n),ctmp(n,n)

      integer taskid,np
      integer i,j,k,icount
      call Parallel_np(np)
      call Parallel_taskid(taskid)
      call dcopy(n*n,0.0d0,0,ctmp,1)
      icount = 0
      do j=1,n
         do i=1,n
            if (mod(icount,np).eq.taskid) then
               do k=1,n
                  c(i,j) = c(i,j) + a(i,k)*b(k,j)
               end do
            end if
            icount = icount + 1
         end do
      end do
      call Parallel_Vector_SumAll(n*n,ctmp)
      do j=1,n
         do i=1,n
            c(i,j) = beta*c(i,j) + alpha*ctmp(i,j)
         end do
      end do
      return 
      end

      subroutine  Dneall_mmm_replicated2(n,alpha,a,b,beta,c,ctmp)
      implicit none
      integer n
      real*8 alpha,beta
      real*8 a(n,n),b(n,n),c(n,n),ctmp(n,n)

      integer taskid,np
      integer i,j,k,icount
      call Parallel_np(np)
      call Parallel_taskid(taskid)
      call dcopy(n*n,0.0d0,0,ctmp,1)
      icount = 0
      do j=1,n
         do i=1,n
            if (mod(icount,np).eq.taskid) then
               do k=1,n
                  c(i,j) = c(i,j) + a(k,i)*b(k,j)
               end do
            end if
            icount = icount + 1
         end do
      end do
      call Parallel_Vector_SumAll(n*n,ctmp)
      do j=1,n
         do i=1,n
            c(i,j) = beta*c(i,j) + alpha*ctmp(i,j)
         end do
      end do
      return
      end

      subroutine  Dneall_mmm_replicated3(n,alpha,a,b,beta,c,ctmp)
      implicit none
      integer n
      real*8 alpha,beta
      real*8 a(n,n),b(n,n),c(n,n),ctmp(n,n)

      integer taskid,np
      integer i,j,k,icount
      call Parallel_np(np)
      call Parallel_taskid(taskid)
      call dcopy(n*n,0.0d0,0,ctmp,1)
      icount = 0
      do j=1,n
         do i=1,n
            if (mod(icount,np).eq.taskid) then
               do k=1,n
                  c(i,j) = c(i,j) + a(i,k)*b(j,k)
               end do
            end if
            icount = icount + 1
         end do
      end do
      call Parallel_Vector_SumAll(n*n,ctmp)
      do j=1,n
         do i=1,n
            c(i,j) = beta*c(i,j) + alpha*ctmp(i,j)
         end do
      end do
      return
      end







c     ****************************************
c     *                                      *
c     *           Dneall_m_trace             *
c     *                                      *
c     ****************************************

      double precision function Dneall_m_trace(mb,M)
      implicit none
      integer mb
      real*8  M(*)

#include "bafdecls.fh"
#include "Dne.fh"

*     **** local variables ****
      integer ms,ms1,ms2,shift2,ishift2,indx,i
      real*8  sum

*     ***** external functions ****
      real*8   DMatrix_trace
      external DMatrix_trace

      if (mparallelized) then
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = mcq(1)*ncq(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if
         sum = 0.0d0
         do ms=ms1,ms2
            shift2 = 1 + (ms-1)*ishift2
            sum =  sum + DMatrix_trace(ne(ms),
     >                                 M(shift2),mcq(ms),
     >                                           int_mb(mc(1,ms)),
     >                                           int_mb(nc(1,ms)),
     >                                 taskid_i,taskid_j)
         end do

      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = ne(1)*ne(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if

         sum = 0.0d0
         do ms=ms1,ms2

           if (ne(ms).le.0) goto 30
           shift2 = (ms-1)*ishift2

           do i=1,ne(ms)
              indx = i + (i-1)*ne(ms) + shift2
              sum = sum + M(indx)
           end do

 30        continue
         end do
      end if

      Dneall_m_trace = sum
      return
      end



c     ****************************************
c     *                                      *
c     *           Dneall_m_sqr_trace         *
c     *                                      *
c     ****************************************
      double precision function Dneall_m_sqr_trace(mb,M)
      implicit none
      integer mb
      real*8  M(*)

#include "bafdecls.fh"
#include "Dne.fh"

*     **** local variables ****
      integer ms,ms1,ms2,shift2,ishift2,indx,i,j
      real*8  sum

*     ***** external functions ****
      real*8   DMatrix_sqr_trace
      external DMatrix_sqr_trace
      if (mparallelized) then
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = mcq(1)*ncq(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if
         sum = 0.0d0
         do ms=ms1,ms2
            shift2 = 1 + (ms-1)*ishift2
            sum =  sum + DMatrix_sqr_trace(ne(ms),
     >                                     M(shift2),mcq(ms),
     >                                     int_mb(mc(1,ms)),
     >                                     int_mb(nc(1,ms)),
     >                                     taskid_i,taskid_j)
         end do

      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = ne(1)*ne(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if

         sum = 0.0d0
         do ms=ms1,ms2

           if (ne(ms).le.0) goto 30
           shift2 = (ms-1)*ishift2

           do j=1,ne(ms)
           do i=1,ne(ms)
              indx = i + (j-1)*ne(ms) + shift2
              sum = sum + M(indx)**2
           end do
           end do

 30        continue
         end do
      end if

      Dneall_m_sqr_trace = sum
      return
      end






c     ****************************************
c     *                                      *
c     *           Dneall_m_scal             *
c     *                                      *
c     ****************************************

      subroutine Dneall_m_scal(mb,alpha,M)
      implicit none
      integer mb
      real*8  alpha
      real*8  M(*)

#include "Dne.fh"

*     **** local variables ****
      integer nn

      if (mparallelized) then
         if (mb.eq.0) then
            nn = mcq(1)*ncq(1)+mcq(2)*ncq(2)
         else
            nn = mcq(mb)*ncq(mb)
         end if
      else
         if (mb.eq.0) then
            nn = ne(1)*ne(1) + ne(2)*ne(2)
         else
            nn = ne(mb)*ne(mb)
         end if
      end if

      call DSCAL_OMP(nn,alpha,M,1)
      return
      end

c     ****************************************
c     *                                      *
c     *           Dneall_mm_sum              *
c     *                                      *
c     ****************************************

      subroutine Dneall_mm_sum(mb,M1,M2,ss)
      implicit none
      integer mb
      real*8  M1(*)
      real*8  M2(*)
      real*8  ss

#include "Dne.fh"

*     **** local variables ****
      integer nn

*     **** external functions ****
      real*8   ddot
      external ddot

      if (mparallelized) then
         if (mb.eq.0) then
            nn = mcq(1)*ncq(1)+mcq(2)*ncq(2)
         else
            nn = mcq(mb)*ncq(mb)
         end if
         ss = ddot(nn,M1,1,M2,1)
         call Parallel_SumAll(ss)
      else
         if (mb.eq.0) then
            nn = ne(1)*ne(1) + ne(2)*ne(2)
         else
            nn = ne(mb)*ne(mb)
         end if
         ss = ddot(nn,M1,1,M2,1)
      end if
      return
      end




c     ****************************************
c     *                                      *
c     *           Dneall_m_diag_scal        *
c     *                                      *
c     ****************************************

      subroutine Dneall_m_diag_scal(mb,diag,M)
      implicit none
      integer mb
      real*8  diag(*)
      real*8  M(*)

#include "Dne.fh"

*     **** local variables ****
      integer ms,ms1,ms2,shift1,shift2,ishift1,ishift2,indx,i

      if (mparallelized) then
         call errquit('Dneall_m_diag_scal not finished',0,0)
      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift1 = ne(1)
            ishift2 = ne(1)*ne(1)
         else
            ms1 = mb
            ms2 = mb
            ishift1 = 0
            ishift2 = 0
         end if

         do ms=ms1,ms2

            if (ne(ms).le.0) goto 30
            shift1 = (ms-1)*ishift1
            shift2 = (ms-1)*ishift2

            do i=1,ne(ms)
               indx = i + (i-1)*ne(ms) + shift2
               M(indx) = M(indx)*diag(i+shift1)
            end do

 30        continue
         end do
      end if

      return
      end

c     ****************************************
c     *                                      *
c     *        Dneall_m_diag_scal_inv       *
c     *                                      *
c     ****************************************

      subroutine Dneall_m_diag_scal_inv(mb,diag,M)
      implicit none
      integer mb
      real*8  diag(*)
      real*8  M(*)

#include "Dne.fh"

*     **** local variables ****
      integer ms,ms1,ms2,shift1,shift2,ishift1,ishift2,indx,i

      if (mparallelized) then
         call errquit('Dneall_m_diag_scal_inv not finished',0,0)
      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift1 = ne(1)
            ishift2 = ne(1)*ne(1)
         else
            ms1 = mb
            ms2 = mb
            ishift1 = 0
            ishift2 = 0
         end if

         do ms=ms1,ms2

            if (ne(ms).le.0) goto 30
            shift1 = (ms-1)*ishift1
            shift2 = (ms-1)*ishift2

            do i=1,ne(ms)
               indx = i + (i-1)*ne(ms) + shift2
               M(indx) = M(indx)/diag(i+shift1)
            end do

 30        continue
         end do
      end if
      return
      end





c     ****************************************
c     *                                      *
c     *        Dneall_m_scale_s22           *
c     *                                      *
c     ****************************************

      subroutine Dneall_m_scale_s22(mb,dte,s22)
      implicit none
      integer mb
      real*8 dte
      real*8 s22(*)


#include "bafdecls.fh"
#include "Dne.fh"

*     **** local variables ****
      integer ms,ms1,ms2,shift2,ishift2,k,j,indx,indxt

      if (mparallelized) then
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = mcq(1)*ncq(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if
         do ms=ms1,ms2
            shift2 = 1 + (ms-1)*ishift2
         
            call DMatrix_s22(ne(ms),dte,
     >                       s22(shift2),mcq(ms),
     >                                   int_mb(mc(1,ms)),
     >                                   int_mb(nc(1,ms)),
     >                       taskid_i,taskid_j)
         end do
      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = ne(1)*ne(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if

         do ms=ms1,ms2
           if (ne(ms).le.0) go to 30
           shift2 = (ms-1)*ishift2

!$OMP DO
           do k=1,ne(ms)
              indx = k + (k-1)*ne(ms) + shift2
              s22(indx) = (1.0d0 - s22(indx))*0.5d0/dte

              do j=k+1,ne(ms)
                 indx  = j + (k-1)*ne(ms) + shift2
                 indxt = k + (j-1)*ne(ms) + shift2

                 s22(indx)  = -s22(indx)*0.5d0/dte
                 s22(indxt) = s22(indx)
              end do
           end do
!$OMP END DO

 30        continue
         end do
      end if

      return
      end


c     ****************************************
c     *                                      *
c     *        Dneall_m_scale_s21           *
c     *                                      *
c     ****************************************

      subroutine Dneall_m_scale_s21(mb,dte,s21)
      implicit none
      integer mb
      real*8 dte
      real*8 s21(*)

#include "bafdecls.fh"
#include "Dne.fh"

*     **** local variables ****
      integer ms,ms1,ms2,shift2,ishift2,k,j,indx,indxt

      if (mparallelized) then
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = mcq(1)*ncq(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if
         do ms=ms1,ms2
            shift2 = 1 + (ms-1)*ishift2

            call DMatrix_s21(ne(ms),dte,
     >                       s21(shift2),mcq(ms),
     >                                   int_mb(mc(1,ms)),
     >                                   int_mb(nc(1,ms)),
     >                       taskid_i,taskid_j)
         end do
      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = ne(1)*ne(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if

         do ms=ms1,ms2
           if (ne(ms).le.0) go to 30
           shift2 = (ms-1)*ishift2

!$OMP DO
           do k=1,ne(ms)
              indx = k + (k-1)*ne(ms) + shift2
              s21(indx) = (1.0d0 - s21(indx))*0.5d0

              do j=k+1,ne(ms)
                 indx  = j + (k-1)*ne(ms) + shift2
                 indxt = k + (j-1)*ne(ms) + shift2

                 s21(indx)  = -s21(indx)*0.5d0
                 s21(indxt) = s21(indx)
              end do
           end do
!$OMP END DO

 30        continue
         end do
      end if

      return
      end


c     ****************************************
c     *                                      *
c     *        Dneall_m_scale_s11           *
c     *                                      *
c     ****************************************

      subroutine Dneall_m_scale_s11(mb,dte,s11)
      implicit none
      integer mb
      real*8 dte
      real*8 s11(*)

#include "bafdecls.fh"
#include "Dne.fh"

*     **** local variables ****
      integer ms,ms1,ms2,shift2,ishift2,k,j,indx,indxt

      if (mparallelized) then
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = mcq(1)*ncq(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if
         do ms=ms1,ms2
            shift2 = 1 + (ms-1)*ishift2

            call DMatrix_s11(ne(ms),dte,
     >                       s11(shift2),mcq(ms),
     >                                   int_mb(mc(1,ms)),
     >                                   int_mb(nc(1,ms)),
     >                       taskid_i,taskid_j)
         end do
      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = ne(1)*ne(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if

         do ms=ms1,ms2
           if (ne(ms).le.0) go to 30
           shift2 = (ms-1)*ishift2

!$OMP DO
           do k=1,ne(ms)
              indx = k + (k-1)*ne(ms) + shift2
              s11(indx) = -s11(indx)*0.5d0*dte

              do j=k+1,ne(ms)
                 indx  = j + (k-1)*ne(ms) + shift2
                 indxt = k + (j-1)*ne(ms) + shift2

                 s11(indx)  = -s11(indx)*0.5d0*dte
                 s11(indxt) = s11(indx)
              end do
           end do
!$OMP END DO

 30        continue
         end do
      end if
      return
      end

c     ****************************************
c     *                                      *
c     *        Dneall_m_Kiril_BTransform     *
c     *                                      *
c     ****************************************

      subroutine Dneall_m_Kiril_BTransform(mb,s12,s21)
      implicit none
      integer mb
      real*8 s12(*),s21(*)

#include "bafdecls.fh"
#include "Dne.fh"

*     **** local variables ****
      integer ms,ms1,ms2,shift2,ishift2,i,j,ii
      real*8 tmp1

      if (mparallelized) then
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = mcq(1)*ncq(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if
         do ms=ms1,ms2
            if (ne(ms).le.0) go to 29
            shift2 = 1 + (ms-1)*ishift2
            call DMatrix_Kiril_BTransform(ne(ms),
     >                       s12(shift2),s21(shift2),
     >                       mcq(ms),int_mb(mc(1,ms)),int_mb(nc(1,ms)),
     >                       taskid_i,taskid_j)
 29         continue
         end do

      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = ne(1)*ne(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if
         do ms=ms1,ms2
            if (ne(ms).le.0) go to 30
            shift2 = 1 + (ms-1)*ishift2
            do i=1,ne(ms)
               do j=1,(i-1)
                  ii = i + j*(i-1) + shift2
                  tmp1 = 0.5d0*(s12(ii)+s21(ii))
                  s12(ii) = tmp1
                  s21(ii) = tmp1
               end do
            end do
 30        continue
         end do

      end if

      return
      end


c     ****************************************
c     *                                      *
c     *        Dneall_m_dmax                 *
c     *                                      *
c     ****************************************

      double precision function Dneall_m_dmax(mb,A)
      implicit none
      integer mb
      real*8 A(*)

#include "bafdecls.fh"
#include "Dne.fh"

*     **** local variables ****
      integer ms,ms1,ms2,shift2,ishift2
      double precision adiff1, adiff2

      integer  idamax
      external idamax

      if (mparallelized) then
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = mcq(1)*ncq(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if
         adiff1 = 0.0d0
         adiff2 = 0.0d0
         do ms=ms1,ms2
           if (ne(ms).le.0) go to 20
            shift2 = 1 + (ms-1)*ishift2

            adiff1 = adiff2
            adiff2 = A(shift2-1+idamax(mcq(ms)*ncq(ms),A(shift2),1))
            adiff2 = dabs(adiff2)
            call DMatrix_MaxAll(adiff2)
            if (adiff2.gt.adiff1) adiff1 = adiff2
 20        continue
         end do

      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = ne(1)*ne(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if

         adiff1 = 0.0d0
         adiff2 = 0.0d0
         do ms=ms1,ms2
           if (ne(ms).le.0) go to 30
           shift2 = 1 + (ms-1)*ishift2

           adiff1 = adiff2      
           adiff2 = A(shift2-1+idamax(ne(ms)*ne(ms),A(shift2),1))
           adiff2 = dabs(adiff2)
           if (adiff2.gt.adiff1) adiff1 = adiff2
 30        continue
         end do
      end if

      Dneall_m_dmax = adiff1
      return
      end




c     ****************************************
c     *                                      *
c     *        Dneall_mm_Expand             *
c     *                                      *
c     ****************************************

      subroutine Dneall_mm_Expand(mb,A,A0)
      implicit none
      integer mb
      real*8 A(*),A0(*)
           
#include "Dne.fh"

*     **** local variables ****
      integer shift2,nn
      
      shift2 = 1
      if (mparallelized) then

         if (mb.eq.0) then
            nn     = mcq(1)*ncq(1) + mcq(2)*ncq(2)
            shift2 = 1
         else if (mb.eq.1) then
            nn     = mcq(1)*ncq(1)
            shift2 = 1 
         else if (mb.eq.2) then
            nn     = mcq(2)*ncq(2)
            shift2 = 1+mcq(1)*ncq(1)
         end if

      else
         if (mb.eq.0) then
            nn     = ne(1)*ne(1) + ne(2)*ne(2)
            shift2 = 1
         else if (mb.eq.1) then
            nn     = ne(1)*ne(1)
            shift2 = 1 
         else if (mb.eq.2) then
            nn     = ne(2)*ne(2)
            shift2 = 1+ne(1)*ne(1)
         end if

      end if

      !call dcopy(nn,A,1,A0(shift2),1)
      call Parallel_shared_vector_copy(.true.,nn,A,A0(shift2))
      return
      end


      subroutine Dneall_m_print(mb,A)
      implicit none
      integer mb
      real*8  A(*)


#include "Dne.fh"

*     **** local variables ****
      integer ms,ms1,ms2,shift2,ishift2,i,j
      integer taskid

      call Parallel_taskid(taskid)
      if (mparallelized) then
         if (mb.eq.0) then
            ms1=1 
            ms2=ispin
            ishift2 = mcq(1)*ncq(1)
         else
            ms1=mb
            ms2=mb
            ishift2 = 0
         end if
         do ms=ms1,ms2
            shift2 = (ms-1)*ishift2
            write(*,*)
            write(*,*) taskid,taskid_i,taskid_j,
     >                 "  Dneall Matrix print, spin=",ms
            do i=1,mcq(ms)
               write(*,'(A,2I2,A,10F10.6)') 
     >                    "taskids=",taskid_i,taskid_j," : ",
     >                    (A(i+(j-1)*mcq(ms)+shift2), j=1,ncq(ms))
            end do
         end do

      else
         if (mb.eq.0) then
            ms1=1 
            ms2=ispin
            ishift2 = ne(1)*ne(1)
         else
            ms1=mb
            ms2=mb
            ishift2 = 0
         end if

         do ms=ms1,ms2
            shift2 = (ms-1)*ishift2
            write(*,*)
            write(*,*) taskid,"  Dneall Matrix print, spin=",ms
            do i=1,ne(ms)
               write(*,'(A,I2,A,10F10.6)') "taskid=",taskid," : ",
     >                    (A(i+(j-1)*ne(ms)+shift2), j=1,ne(ms))
            end do
         end do
      end if

      return
      end



      real*8 function  Dneall_m_value(mb,ms,i,j,A)
      implicit none
      integer mb,ms,i,j
      real*8  A(*)

#include "bafdecls.fh"
#include "Dne.fh"

*     **** local variables ****
      integer ishift2,shift2
      real*8 w

*     **** external functions ****
      real*8   DMatrix_m_get_value
      external DMatrix_m_get_value
      
      w = 0.0d0
      if (ne(ms).gt.0)  then
         if (mparallelized) then

            if (mb.eq.0) then
               ishift2 = mcq(1)*ncq(1)
            else
               ishift2 = 0
            end if
            shift2 = 1 + (ms-1)*ishift2
            w = DMatrix_m_get_value(i,j,A(shift2),
     >                 mcq(ms),int_mb(mc(1,ms)),int_mb(nc(1,ms)))
         else
            shift2 = (ms-1)*ne(1)*ne(1)
            w = A(i+(j-1)*ne(ms)+shift2)
         end if
      end if

      Dneall_m_value = w
      return
      end



      subroutine Dneall_m_set_value(w,mb,ms,i,j,A)
      implicit none
      real*8 w
      integer mb,ms,i,j
      real*8  A(*)

#include "bafdecls.fh"
#include "Dne.fh"

*     **** local variables ****
      integer ishift2,shift2

      if (ne(ms).gt.0)  then
         if (mparallelized) then
            if (mb.eq.0) then
               ishift2 = mcq(1)*ncq(1)
            else
               ishift2 = 0
            end if
            shift2 = 1 + (ms-1)*ishift2
            call DMatrix_m_set_value(w,i,j,A(shift2),
     >                 mcq(ms),int_mb(mc(1,ms)),int_mb(nc(1,ms)))
         else
            shift2 = (ms-1)*ne(1)*ne(1)
            A(i+(j-1)*ne(ms)+shift2) = w
         end if
      end if

      return
      end




      subroutine Dneall_m_add_value(w,mb,ms,i,j,A)
      implicit none
      real*8 w
      integer mb,ms,i,j
      real*8  A(*)

#include "bafdecls.fh"
#include "Dne.fh"

*     **** local variables ****
      integer ishift2,shift2

      if (ne(ms).gt.0)  then
         if (mparallelized) then
            if (mb.eq.0) then
               ishift2 = mcq(1)*ncq(1)
            else
               ishift2 = 0
            end if
            shift2 = 1 + (ms-1)*ishift2
            call DMatrix_m_add_value(w,i,j,A(shift2),
     >                 mcq(ms),int_mb(mc(1,ms)),int_mb(nc(1,ms)))
         else
            if (mb.eq.0) then
               shift2 = (ms-1)*ne(1)*ne(1)
            else
               shift2 = 0
            end if
            A(i+(j-1)*ne(ms)+shift2) = A(i+(j-1)*ne(ms)+shift2) + w
         end if
      end if

      return
      end



      real*8 function  Dneall_sm_value(i,j,A)
      implicit none
      integer i,j
      real*8  A(*)

#include "bafdecls.fh"
#include "Dne.fh"

*     **** local variables ****
      real*8 w

*     **** external functions ****
      real*8   DMatrix_m_get_value
      external DMatrix_m_get_value

      w = 0.0d0
      if (ne(2).gt.0)  then
         if (mparallelized) then
            w = DMatrix_m_get_value(i,j,A,
     >                 mdq,int_mb(md(1)),int_mb(nd(1)))
         else
            w = A(i+(j-1)*ne(1))
         end if
      end if

      Dneall_sm_value = w
      return
      end


*     ***********************************
*     *                                 *
*     *       Dneall_m_add_sw1sw2       *
*     *                                 *
*     ***********************************
      subroutine Dneall_m_add_sw1sw2(mb,nprj,scal,sw1,sw2,S)
      implicit none
      integer mb,nprj
      real*8 scal
      real*8 sw1(*)
      real*8 sw2(*)
      real*8 S(*)

#include "bafdecls.fh"
#include "Dne.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value
      integer ms,q,n,nn,nnq,shift1,shift2,ishift1,ishift2
      integer i,j,ii,jj,pto,taskid,ms1,ms2,mshiftp,mshiftm
      integer sw1b(2),sw2b(2)
      real*8  w

*     **** external functions ****
      real*8   ddot
      external ddot


      if (parallelized) then
         if (mb.eq.0) then
            nn  = ne(1)  + ne(2)
            nnq = neq(1) + neq(2)
            mshiftp = 0
            mshiftm = 0
         else
            nn  = ne(mb)
            nnq = neq(mb)
            mshiftp = (mb-1)*neq(1)
            mshiftm = (mb-1)*ne(1)
         end if
         value = BA_push_get(mt_dbl,nn*nprj,'sw1b',sw1b(2),sw1b(1))
         value = value.and.
     >           BA_push_get(mt_dbl,nn*nprj,'sw2b',sw2b(2),sw2b(1))
         if (.not.value)
     >   call errquit('Dneall_m_add_sw1sw2:push stack',0,MA_ERR)

         !call dcopy(nn*nprj,0.0d0,0,dbl_mb(sw1b(1)),1)
         !call dcopy(nn*nprj,0.0d0,0,dbl_mb(sw2b(1)),1)
         call Parallel_shared_vector_zero(.false.,nn*nprj,
     >                                    dbl_mb(sw1b(1)))
         call Parallel_shared_vector_zero(.true.,nn*nprj,
     >                                    dbl_mb(sw2b(1)))
!$OMP MASTER
         do q=1,nnq
            call Dneall_qton(q+mshiftp,n)
            n = n-mshiftm
            !do prj = 1,nprj
            !   sw1b(n,prj) = sw1(q,prj)
            !   sw2b(n,prj) = sw2(q,prj)
            !end do
            call dcopy(nprj,sw1(q),nnq,dbl_mb(sw1b(1)+n-1),nn)
            call dcopy(nprj,sw2(q),nnq,dbl_mb(sw2b(1)+n-1),nn)
         end do
!$OMP END MASTER
!$OMP BARRIER
         call D1dB_SumAll(nn*nprj,dbl_mb(sw1b(1)))
         call D1dB_SumAll(nn*nprj,dbl_mb(sw2b(1)))

         if (mparallelized) then
            call Parallel_taskid(taskid)
            if (mb.eq.0) then
               ms1 = 1
               ms2 = ispin
               ishift1 = ne(1)
            else
               ms1 = mb
               ms2 = mb
               ishift1 = 0
            end if
            do ms=ms1,ms2
               shift1 =   (ms-1)*ishift1
               do j=1,ne(ms)
               do i=1,ne(ms)
                  call DMatrix_m_ijtoiijjp(i,j,ii,jj,pto,
     >                                     mcq(ms),
     >                                     int_mb(mc(1,ms)),
     >                                     int_mb(nc(1,ms)))
                  if (pto.eq.taskid) then
!$OMP MASTER
                     w = ddot(nprj,dbl_mb(sw1b(1)+i-1+shift1),nn,
     >                             dbl_mb(sw2b(1)+j-1+shift1),nn)
                     call Dneall_m_add_value(scal*w,mb,ms,i,j,S)
!$OMP END MASTER
!$OMP BARRIER
                  end if
               end do
               end do
            end do

         else
            if (mb.eq.0) then
               ms1 = 1
               ms2 = ispin
               ishift1 = ne(1)
               ishift2 = ne(1)*ne(1)
            else
               ms1 = mb
               ms2 = mb
               ishift1 = 0
               ishift2 = 0
            end if
            do ms=ms1,ms2
               if (ne(ms).gt.0) then
                  shift1 =   (ms-1)*ishift1
                  shift2 = 1+(ms-1)*ishift2
                  call DGEMM_OMP('N','T',ne(ms),ne(ms),nprj,
     >                       scal,
     >                       dbl_mb(sw1b(1)+shift1), nn,
     >                       dbl_mb(sw2b(1)+shift1), nn,
     >                       1.0d0,
     >                       S(shift2),ne(ms))
               end if
            end do
         end if
         value =           BA_pop_stack(sw2b(2))
         value = value.and.BA_pop_stack(sw1b(2))
         if (.not.value)
     >   call errquit('Dneall_m_add_sw1sw2:pop stack',1,MA_ERR)
      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift1 = ne(1)
            ishift2 = ne(1)*ne(1)
            nn = ne(1) + ne(2)
         else
            ms1 = mb
            ms2 = mb
            ishift1 = 0
            ishift2 = 0
            nn = ne(mb)
         end if
         do ms=ms1,ms2
            if (ne(ms).gt.0) then
               shift1 = 1+(ms-1)*ishift1
               shift2 = 1+(ms-1)*ishift2
               call DGEMM_OMP('N','T',
     >                    ne(ms),ne(ms),nprj,
     >                    scal,
     >                    sw1(shift1), nn,
     >                    sw2(shift1), nn,
     >                    1.0d0,
     >                    S(shift2), ne(ms))
            end if
         end do
      end if

      return
      end





c     ****************************************
c     *                                      *
c     *        Dneall_m_HmltimesSA           *
c     *                                      *
c     ****************************************

      subroutine Dneall_m_HmltimesSA(mb,A,SA)
      implicit none
      integer mb
      real*8 A(*)
      real*8 SA(*)

#include "bafdecls.fh"
#include "Dne.fh"

*     **** local variables ****
      integer ms,ms1,ms2,shift1,shift2,ishift1,ishift2,k,j,indx1,indx2

      if (mparallelized) then
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift1 = ne(1)
            ishift2 = mcq(1)*ncq(1)
         else
            ms1 = mb
            ms2 = mb
            ishift1 = 0
            ishift2 = 0
         end if

         do ms=ms1,ms2
             shift1 = 1 + (ms-1)*ishift1
             shift2 = 1 + (ms-1)*ishift2
             call DMatrix_HmltimesSA(ne(ms),
     >                SA(shift1),A(shift2),
     >                mcq(ms),int_mb(mc(1,ms)),int_mb(nc(1,ms)),
     >                taskid_i,taskid_j)
         end do
      
      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift1 = ne(1)
            ishift2 = ne(1)*ne(1)
         else
            ms1 = mb
            ms2 = mb
            ishift1 = 0
            ishift2 = 0
         end if

         do ms=ms1,ms2
           if (ne(ms).le.0) go to 30
           shift1 = 1 + (ms-1)*ishift1
           shift2 = 1 + (ms-1)*ishift2

           indx1 = shift1
           indx2 = shift2
           do k=1,ne(ms)
              do j=1,ne(ms)
                A(indx2) = A(indx2)*SA(indx1)
                indx2 = indx2 + 1
              end do
              indx1 = indx1 + 1
           end do
   
 30        continue
         end do
      end if

      return
      end







c     ****************************************
c     *                                      *
c     *        Dneall_m_HmldivideSA          *
c     *                                      *
c     ****************************************

      subroutine Dneall_m_HmldivideSA(mb,A,SA)
      implicit none
      integer mb
      real*8 A(*)
      real*8 SA(*)

#include "bafdecls.fh"
#include "Dne.fh"

*     **** local variables ****
      integer ms,ms1,ms2,shift1,shift2,ishift1,ishift2,k,j,indx1,indx2

      if (mparallelized) then
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift1 = ne(1)
            ishift2 = mcq(1)*ncq(1)
         else
            ms1 = mb
            ms2 = mb
            ishift1 = 0
            ishift2 = 0
         end if

         do ms=ms1,ms2
             shift1 = 1 + (ms-1)*ishift1
             shift2 = 1 + (ms-1)*ishift2
             call DMatrix_HmldivideSA(ne(ms),
     >                SA(shift1),A(shift2),
     >                mcq(ms),int_mb(mc(1,ms)),int_mb(nc(1,ms)),
     >                taskid_i,taskid_j)
         end do
      
      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift1 = ne(1)
            ishift2 = ne(1)*ne(1)
         else
            ms1 = mb
            ms2 = mb
            ishift1 = 0
            ishift2 = 0
         end if

         do ms=ms1,ms2
           if (ne(ms).le.0) go to 30
           shift1 = 1 + (ms-1)*ishift1
           shift2 = 1 + (ms-1)*ishift2

           indx1 = shift1
           indx2 = shift2
           do k=1,ne(ms)
              do j=1,ne(ms)
                A(indx2) = A(indx2)/SA(indx1)
                indx2 = indx2 + 1
              end do
              indx1 = indx1 + 1
           end do
   
 30        continue
         end do
      end if

      return
      end






c     ****************************************
c     *                                      *
c     *        Dneall_m_HmldivideDplusD      *
c     *                                      *
c     ****************************************

      subroutine Dneall_m_HmldivideDplusD(mb,A,D)
      implicit none
      integer mb
      real*8 A(*)
      real*8 D(*)

#include "bafdecls.fh"
#include "Dne.fh"

*     **** local variables ****
      integer ms,ms1,ms2,shift1,shift2,ishift1,ishift2,k,j
      integer indx1,indx2,indx3

      if (mparallelized) then
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift1 = ne(1)
            ishift2 = mcq(1)*ncq(1)
         else
            ms1 = mb
            ms2 = mb
            ishift1 = 0
            ishift2 = 0
         end if

         do ms=ms1,ms2
             shift1 = 1 + (ms-1)*ishift1
             shift2 = 1 + (ms-1)*ishift2
             call DMatrix_HmldivideDplusD(ne(ms),
     >                D(shift1),A(shift2),
     >                mcq(ms),int_mb(mc(1,ms)),int_mb(nc(1,ms)),
     >                taskid_i,taskid_j)
         end do


      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift1 = ne(1)
            ishift2 = ne(1)*ne(1)
         else
            ms1 = mb
            ms2 = mb
            ishift1 = 0
            ishift2 = 0
         end if

         do ms=ms1,ms2
           if (ne(ms).le.0) go to 30
           shift1 = 1 + (ms-1)*ishift1
           shift2 = 1 + (ms-1)*ishift2

           !indx1 = shift1
           !indx2 = shift2
!$OMP DO
           do k=1,ne(ms)
              !indx3 = shift1
              indx1 = shift1+(k-1)
              do j=1,ne(ms)
                indx2 = shift2+(k-1)*ne(ms)+(j-1)
                indx3 = shift1+(j-1)
                A(indx2) = A(indx2)/(D(indx1)+D(indx3))
                !indx2 = indx2 + 1
                !indx3 = indx3 + 1
              end do
              !indx1 = indx1 + 1
           end do
!$OMP END DO

 30        continue
         end do
      end if

      return
      end




c     ****************************************
c     *                                      *
c     *        Dneall_m_Hmlfweightscale      *
c     *                                      *
c     ****************************************

      subroutine Dneall_m_Hmlfweightscale(mb,A,fw)
      implicit none
      integer mb
      real*8 A(*)
      real*8 fw(*)

#include "bafdecls.fh"
#include "Dne.fh"

*     **** local variables ****
      integer ms,ms1,ms2,shift1,shift2,ishift1,ishift2,k,j
      integer indx1,indx2,indx3

      if (mparallelized) then
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift1 = ne(1)
            ishift2 = mcq(1)*ncq(1)
         else
            ms1 = mb
            ms2 = mb
            ishift1 = 0
            ishift2 = 0
         end if

         do ms=ms1,ms2
             shift1 = 1 + (ms-1)*ishift1
             shift2 = 1 + (ms-1)*ishift2
             call DMatrix_Hmlfweightscale(ne(ms),
     >                fw(shift1),A(shift2),
     >                mcq(ms),int_mb(mc(1,ms)),int_mb(nc(1,ms)),
     >                taskid_i,taskid_j)
         end do


      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift1 = ne(1)
            ishift2 = ne(1)*ne(1)
         else
            ms1 = mb
            ms2 = mb
            ishift1 = 0
            ishift2 = 0
         end if

         do ms=ms1,ms2
           if (ne(ms).le.0) go to 30
           shift1 = 1 + (ms-1)*ishift1
           shift2 = 1 + (ms-1)*ishift2

           indx1 = shift1
           indx2 = shift2
           do k=1,ne(ms)
              indx3 = shift1
              do j=1,ne(ms)
                if ((fw(indx1)+fw(indx3)).gt.1.0d-9)
     >             A(indx2) = A(indx2)
     >                      *(2.0d0*fw(indx1)/(fw(indx1)+fw(indx3)))
                indx2 = indx2 + 1
                indx3 = indx3 + 1
              end do
              indx1 = indx1 + 1
           end do

 30        continue
         end do
      end if

      return
      end



c     ****************************************
c     *                                      *
c     *        Dneall_m_eye                  *
c     *                                      *
c     ****************************************

      subroutine Dneall_m_eye(mb,A,alpha)
      implicit none
      integer mb
      real*8 A(*)
      real*8 alpha

#include "bafdecls.fh"
#include "Dne.fh"


*     **** local variables ****
      integer ms,ms1,ms2,shift2,ishift2,k
      integer indx2

      if (mparallelized) then
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = mcq(1)*ncq(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if

         do ms=ms1,ms2
             shift2 = 1 + (ms-1)*ishift2
             call DMatrix_eye(ne(ms),ne(ms),alpha,
     >                A(shift2),
     >                mcq(ms),int_mb(mc(1,ms)),int_mb(nc(1,ms)),
     >                taskid_i,taskid_j)
         end do

      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = ne(1)*ne(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if

         do ms=ms1,ms2
           if (ne(ms).le.0) go to 30
           shift2 = 1 + (ms-1)*ishift2
           !call dcopy(ne(ms)*ne(ms),0.0d0,0,A(shift2),1)
           call Parallel_shared_vector_zero(.true.,ne(ms)*ne(ms),
     >                                      A(shift2))
           indx2 = shift2
!$OMP DO
           do k=1,ne(ms)
              indx2 = shift2 + (k-1)*(1+ne(ms))
              A(indx2) = alpha
              !indx2    = indx2 + (1+ne(ms))
           end do
!$OMP END DO

 30        continue
         end do
      end if

      return
      end

c     ****************************************
c     *                                      *
c     *        Dneall_m_setdiag              *
c     *                                      *
c     ****************************************
      subroutine Dneall_m_setdiag(mb,A,alpha)
      implicit none
      integer mb
      real*8 A(*)
      real*8 alpha(*)

#include "bafdecls.fh"
#include "Dne.fh"


*     **** local variables ****
      integer ms,ms1,ms2,shift2,ishift2,k
      integer indx2,shift1

      if (mparallelized) then
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = mcq(1)*ncq(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if

         do ms=ms1,ms2
             shift1 = 1 + (ms-1)*ne(1)
             shift2 = 1 + (ms-1)*ishift2
             call DMatrix_setdiag(ne(ms),ne(ms),alpha(shift1),
     >                A(shift2),
     >                mcq(ms),int_mb(mc(1,ms)),int_mb(nc(1,ms)),
     >                taskid_i,taskid_j)
         end do

      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = ne(1)*ne(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if

         do ms=ms1,ms2
           if (ne(ms).le.0) go to 30
           shift1 = (ms-1)*ne(1)
           shift2 = 1 + (ms-1)*ishift2
           call dcopy(ne(ms)*ne(ms),0.0d0,0,A(shift2),1)
           indx2 = shift2
           do k=1,ne(ms)
              A(indx2) = alpha(k+shift1)
              indx2    = indx2 + (1+ne(ms))
           end do

 30        continue
         end do
      end if

      return
      end






c     ****************************************
c     *                                      *
c     *        Dneall_w_eye                  *
c     *                                      *
c     ****************************************

      subroutine Dneall_w_eye(mb,A,alpha)
      implicit none
      integer mb
      complex*16 A(*)
      complex*16 alpha

#include "bafdecls.fh"
#include "Dne.fh"


*     **** local variables ****
      integer ms,ms1,ms2,shift2,ishift2,k
      integer indx2

      if (mparallelized) then
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = mcq(1)*ncq(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if

         do ms=ms1,ms2
             shift2 = 1 + (ms-1)*ishift2
             call CMatrix_eye(ne(ms),ne(ms),alpha,
     >                A(shift2),
     >                mcq(ms),int_mb(mc(1,ms)),int_mb(nc(1,ms)),
     >                taskid_i,taskid_j)
         end do

      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = ne(1)*ne(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if

         do ms=ms1,ms2
           if (ne(ms).le.0) go to 30
           shift2 = 1 + (ms-1)*ishift2
           call dcopy(2*ne(ms)*ne(ms),0.0d0,0,A(shift2),1)
           indx2 = shift2
           do k=1,ne(ms)
              A(indx2) = alpha
              indx2    = indx2 + (1+ne(ms))
           end do

 30        continue
         end do
      end if

      return
      end






c     ****************************************
c     *                                      *
c     *        Dneall_f_Sortho               *
c     *                                      *
c     ****************************************

      subroutine Dneall_f_Sortho(mb,U,SU,npack1)
      integer mb
      complex*16 U(*),SU(*)
      integer    npack1      

#include "bafdecls.fh"           
#include "errquit.fh"
#include "Dne.fh"

*     **** local variables ****
      integer taskid
      integer ms,ms1,ms2,n,shift,j,k,indxk,indxj
      integer V(2),tmp2(2),S(2),A(2)
      integer jj,kk,jcur,kcur

      real*8 w
      common /Dneall_f_Sortho_tmp/ w

 
      if (parallelized) then
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
         else
            ms1 = mb
            ms2 = mb
         end if
         shift = neq(1)*npack1

         if (.not.BA_push_get(mt_dcpl,npack1,'A',A(2),A(1)))
     >   call errquit('Dneall_f_Sortho:out of stack memory',0,MA_ERR)

         do ms=ms1,ms2

            kcur = np_j-1
            kk   = int_mb(na(1,ms)+kcur) - 1

            do k =ne(ms),1,-1
               if (kcur.eq.taskid_j) then
                 indxk = 1 + kk*npack1    + (ms-1)*shift
                  call Pack_cc_dot(1,U(indxk),SU(indxk),w)
!$OMP MASTER
                  w        = 1.0d0/dsqrt(w)
!$OMP END MASTER
!$OMP BARRIER
                  call Pack_c_SMul1(1,w,U(indxk))
                  !call dcopy(2*npack1,U(indxk),1,dcpl_mb(A(1)),1)
                  call Parallel_shared_vector_copy(.true.,2*npack1,
     >                                        U(indxk),
     >                                        dcpl_mb(A(1)))
               end if
               if (kcur.gt.0) 
     >            call D1dB_Brdcst_values(kcur,2*npack1,dcpl_mb(A(1)))

c              *** set j = k+1 ***
               jj   = kk
               jcur = kcur

               jj = jj - 1
               if (jj.lt.0) then
                  jcur = jcur - 1
                  jj = int_mb(na(1,ms) + jcur) - 1
               end if

               do j=k-1,1,-1
                  if (jcur.eq.taskid_j) then
                     indxj = 1 + jj*npack1    + (ms-1)*shift
                     call Pack_cc_dot(1,dcpl_mb(A(1)),SU(indxj),w)
!$OMP MASTER
                     w = -w
!$OMP END MASTER
!$OMP BARRIER
                     call Pack_cc_daxpy(1,w,dcpl_mb(A(1)),U(indxj))
                  end if

                  jj = jj - 1
                  if (jj.lt.0) then
                     jcur = jcur - 1
                     jj = int_mb(na(1,ms) + jcur) - 1
                  end if

               end do

               kk = kk - 1
               if (kk.lt.0) then
                  kcur = kcur - 1
                  kk = int_mb(na(1,ms) + kcur) - 1
               end if

            end do
         end do

c        **** deallocate local memory ****
         if (.not.Ma_pop_stack(A(2)))
     >    call errquit('Dneall_f_Sortho: popping stack memory',1,MA_ERR)
     
          
c     **** not parallized ****
      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
         else
            ms1 = mb
            ms2 = mb
         end if
         shift = ne(1)*npack1

         !**** orthogonalize from the bottom -> up ****
         do ms=ms1,ms2
         do k=ne(ms),1,-1
            indxk = 1+(k-1)*npack1 + (ms-1)*shift
            call Pack_cc_dot(1,U(indxk),SU(indxk),w)
!$OMP MASTER
            w = 1.0d0/dsqrt(w)
!$OMP END MASTER
!$OMP BARRIER
c            call Pack_c_SMul(1,w,U(indxk),U(indxk))
            call Pack_c_SMul1(1,w,U(indxk))

            do j=k-1,1,-1
               indxj = 1+(j-1)*npack1 + (ms-1)*shift
               call Pack_cc_dot(1,U(indxk),SU(indxj),w)
!$OMP MASTER
               w = -w
!$OMP END MASTER
!$OMP BARRIER
               call Pack_cc_daxpy(1,w,U(indxk),U(indxj))
            end do
         end do
         end do

      end if

      return
      end 

      subroutine Dneall_m_gather(nall,npack,indx,A,B)
      implicit none
      integer nall,npack,indx(*)
      real*8 A(*),B(*)
      integer i
      !call dcopy(nall,0.0d0,0,B,1)
      call Parallel_shared_vector_zero(.true.,nall,B)
!$OMP DO
      do i=1,npack
         B(indx(i)) = A(i)
      end do  
!$OMP END DO
      call Parallel_Vector_SumAll(nall,B)
      return
      end

      subroutine DneAll_m_scatter(npack,indx,A,B)
      implicit none
      integer npack,indx(*)
      real*8 A(*),B(*)
      integer i
!$OMP DO
      do i=1,npack
         B(i) = A(indx(i))
      end do  
!$OMP END DO
      return
      end


c     ****************************************
c     *                                      *
c     *        Dneall_fm_QR                  *
c     *                                      *
c     ****************************************

      subroutine Dneall_fm_QR(mb,Q,npack1,R)
      implicit none
      integer mb
      complex*16 Q(*)
      integer    npack1      
      real*8     R(*)

#include "bafdecls.fh"           
#include "errquit.fh"
#include "Dne.fh"

c     **** local variables ****
      integer ms,ms1,ms2,shift,shift2,j,k,indxk,indxj,indxm
      integer jj,kk,rr,ss,jcur,kcur,rcur,scur
      integer tmp(2)
      real*8 w,RRM

      if (parallelized) then

c        **** allocate temporary memory ****
         if (.not.BA_push_get(mt_dcpl,npack1,'tmp',tmp(2),tmp(1)))
     >   call errquit('Dneall_fm_QR:out of stack memory',0,MA_ERR)

c        **** mparallized ****
         if (mparallelized) then

         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            shift2 = mcq(1)*ncq(1)
         else
            ms1 = mb
            ms2 = mb
            shift2 = 0
         end if
         shift = neq(1)*npack1

         call dcopy((mcq(1)*ncq(1)+mcq(2)*ncq(2)),0.0d0,0,R,1)
         do ms=ms1,ms2
            rcur = 0
            scur = 0
            rr   = 0
            ss   = 0
            kcur = 0
            kk   = 0
            do k=1,ne(ms)

               if (kcur.eq.taskid_j) then
                  indxk = 1 + kk*npack1    + (ms-1)*shift
                  call Pack_cc_dot(1,Q(indxk),Q(indxk),w)
                  RRM = dsqrt(w)
                  w        = 1.0d0/RRM
c                  call Pack_c_SMul(1,w,Q(indxk),Q(indxk))
                  call Pack_c_SMul1(1,w,Q(indxk))
                  call dcopy(2*npack1,Q(indxk),1,dcpl_mb(tmp(1)),1)
               end if

*              *** Brdcst Q(indxk) ***
               call D1dB_Brdcst_values(kcur,2*npack1,dcpl_mb(tmp(1)))
               call D1dB_Brdcst_values(kcur,1,RRM)

c              *** set j = k+1 ***
               jj   = kk
               jcur = kcur
               jj = jj + 1
               if (jj.ge.int_mb(na(1,ms) + jcur)) then
                  jcur = jcur + 1
                  jj = 0
               end if

               rr   = 0
               rcur = 0
               do j=1,k
                  rr = rr + 1
                  if (rr.ge.int_mb(mc(1,ms) + rcur)) then
                     rcur = rcur + 1
                     rr = 0
                  end if
               end do
               if ((rcur.eq.taskid_i).and.(scur.eq.taskid_j)) then
                  indxm = (rr+1) + ss*mcq(ms) + (ms-1)*shift2
                  R(indxm) = RRM
               end if

               do j=k+1,ne(ms)
                  if (jcur.eq.taskid_j) then
                     indxj = 1 + jj*npack1    + (ms-1)*shift
                     call Pack_cc_dot(1,dcpl_mb(tmp(1)),Q(indxj),w)
                     RRM = w
                     w = -w
                     call Pack_cc_daxpy(1,w,dcpl_mb(tmp(1)),Q(indxj))
                  end if
                  call D1dB_Brdcst_values(jcur,1,RRM)

                  if ((rcur.eq.taskid_i).and.(scur.eq.taskid_j)) then
                     indxm = (rr+1) + ss*mcq(ms) + (ms-1)*shift2
                     R(indxm) = RRM
                  end if

                  jj = jj + 1
                  if (jj.ge.int_mb(na(1,ms) + jcur)) then
                     jcur = jcur + 1
                     jj = 0
                  end if

                  rr = rr + 1
                  if (rr.ge.int_mb(mc(1,ms) + rcur)) then
                     rcur = rcur + 1
                     rr = 0
                  end if
               end do

               kk = kk + 1
               if (kk.ge.int_mb(na(1,ms) + kcur)) then
                  kcur = kcur + 1
                  kk = 0
               end if

               ss = ss + 1
               if (ss.ge.int_mb(nc(1,ms) + scur)) then
                  scur = scur + 1
                  ss = 0
               end if

            end do
         end do


c        **** not mparallized ****
         else

         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            shift2 = ne(1)*ne(1)
         else
            ms1 = mb
            ms2 = mb
            shift2 = 0
         end if
         shift = neq(1)*npack1

         call dcopy((ne(1)*ne(1)+ne(2)*ne(2)),0.0d0,0,R,1)
         do ms=ms1,ms2
            kcur = 0
            kk   = 0
            do k=1,ne(ms)

               if (kcur.eq.taskid_j) then
                  indxk = 1 + kk*npack1    + (ms-1)*shift
                  indxm = k + (k-1)*ne(ms) + (ms-1)*shift2
                  call Pack_cc_dot(1,Q(indxk),Q(indxk),w)
                  R(indxm) = dsqrt(w)
                  w        = 1.0d0/R(indxm)
c                  call Pack_c_SMul(1,w,Q(indxk),Q(indxk))
                  call Pack_c_SMul1(1,w,Q(indxk))
                  call dcopy(2*npack1,Q(indxk),1,dcpl_mb(tmp(1)),1)
               end if

c              *** Brdcst Q(indxk) ***
               call D1dB_Brdcst_values(kcur,2*npack1,dcpl_mb(tmp(1)))
         
c              *** set j = k+1 ***
               jj   = kk
               jcur = kcur
               jj = jj + 1
               if (jj.ge.int_mb(na(1,ms) + jcur)) then
                  jcur = jcur + 1
                  jj = 0
               end if
            
               do j=k+1,ne(ms)
                  if (jcur.eq.taskid_j) then
                     indxj = 1 + jj*npack1    + (ms-1)*shift
                     indxm = k + (j-1)*ne(ms) + (ms-1)*shift2
                     call Pack_cc_dot(1,dcpl_mb(tmp(1)),Q(indxj),w)
                     R(indxm) = w
                     w = -w
                     call Pack_cc_daxpy(1,w,dcpl_mb(tmp(1)),Q(indxj))
                  end if

                  jj = jj + 1
                  if (jj.ge.int_mb(na(1,ms) + jcur)) then
                     jcur = jcur + 1
                     jj = 0
                  end if
               end do

               kk = kk + 1
               if (kk.ge.int_mb(na(1,ms) + kcur)) then
                  kcur = kcur + 1
                  kk = 0
               end if

            end do
         end do
         call D1dB_Vector_SumAll((ne(1)*ne(1)+ne(2)*ne(2)),R)
         end if

c        **** deallocate local memory ****
         if (.not.BA_pop_stack(tmp(2)))
     >   call errquit('Dneall_fm_QR:popping stack memory',0,MA_ERR)

c     **** not parallized ****
      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            shift2 = ne(1)*ne(1)
         else
            ms1 = mb
            ms2 = mb
            shift2 = 0
         end if
         shift = ne(1)*npack1

c         **** modified Gram-Schmidt ****
         call dcopy((ne(1)*ne(1)+ne(2)*ne(2)),0.0d0,0,R,1)
         do ms=ms1,ms2
         do k=1,ne(ms)
            indxk = 1 + (k-1)*npack1 + (ms-1)*shift
            indxm = k + (k-1)*ne(ms) + (ms-1)*shift2
            call Pack_cc_dot(1,Q(indxk),Q(indxk),w)
            R(indxm) = dsqrt(w)
            w        = 1.0d0/R(indxm)
c            call Pack_c_SMul(1,w,Q(indxk),Q(indxk))
            call Pack_c_SMul1(1,w,Q(indxk))

            do j=k+1,ne(ms)
               indxj = 1 + (j-1)*npack1 + (ms-1)*shift
               indxm = k + (j-1)*ne(ms) + (ms-1)*shift2
               call Pack_cc_dot(1,Q(indxk),Q(indxj),w)
               R(indxm) = w
               w = -w
               call Pack_cc_daxpy(1,w,Q(indxk),Q(indxj))
            end do
         end do
         end do
      end if

      return
      end 

c
c     ****************************************
c     *                                      *
c     *        Dneall_4m_size                *
c     *                                      *
c     ****************************************
      subroutine Dneall_4m_size(mb,size)
      implicit none
      integer mb
      integer size

#include "Dne.fh"
     
      if (mparallelized) then
         if (mb.eq.0) then
            size = m2cq(1)*n2cq(1) + m2cq(2)*n2cq(2)
         else
            size = m2cq(mb)*n2cq(mb)
         end if
      else
         if (mb.eq.0) then 
            size = 4*(ne(1)*ne(1) + ne(2)*ne(2))
         else
            size = 4*ne(mb)*ne(mb)
         end if
      end if
     
      return
      end

c     ****************************************
c     *                                      *
c     *        Dneall_4m_allocate_block      *
c     *                                      *
c     ****************************************
      logical function Dneall_4m_allocate_block(mb,nb,hml)
      implicit none
      integer mb,nb
      integer hml(2)

#include "Dne.fh"
#include "bafdecls.fh"

      integer size
     
      if (mparallelized) then
         if (mb.eq.0) then
            size = m2cq(1)*n2cq(1) + m2cq(2)*n2cq(2)
         else
            size = m2cq(mb)*n2cq(mb)
         end if
      else
         if (mb.eq.0) then 
            size = 4*(ne(1)*ne(1) + ne(2)*ne(2))
         else
            size = 4*ne(mb)*ne(mb)
         end if
      end if
     
      Dneall_4m_allocate_block 
     > = BA_alloc_get(mt_dbl,nb*size,'hmlab',hml(2),hml(1))
      return
      end

c     ****************************************
c     *                                      *
c     *        Dneall_4m_allocate             *
c     *                                      *
c     ****************************************
      logical function Dneall_4m_allocate(mb,hml)
      implicit none
      integer mb
      integer hml(2)

#include "Dne.fh"
#include "bafdecls.fh"

      integer size
     
      if (mparallelized) then
         if (mb.eq.0) then
            size = m2cq(1)*n2cq(1) + m2cq(2)*n2cq(2)
         else
            size = m2cq(mb)*n2cq(mb)
         end if
      else
         if (mb.eq.0) then 
            size = 4*(ne(1)*ne(1) + ne(2)*ne(2))
         else
            size = 4*ne(mb)*ne(mb)
         end if
      end if
     
      Dneall_4m_allocate = BA_alloc_get(mt_dbl,size,'hma',hml(2),hml(1))
      return
      end

c     ****************************************
c     *                                      *
c     *        Dneall_4m_push_get            *
c     *                                      *
c     ****************************************
      logical function Dneall_4m_push_get(mb,hml)
      implicit none
      integer mb
      integer hml(2)

#include "Dne.fh"
#include "bafdecls.fh"

      integer size

      if (mparallelized) then
         if (mb.eq.0) then
            size = m2cq(1)*n2cq(1) + m2cq(2)*n2cq(2)
         else
            size = m2cq(mb)*n2cq(mb)
         end if
      else
         if (mb.eq.0) then
            size = 4*(ne(1)*ne(1) + ne(2)*ne(2))
         else
            size = 4*ne(mb)*ne(mb)
         end if
      end if

      Dneall_4m_push_get = BA_push_get(mt_dbl,size,'hmls',hml(2),hml(1))
      return
      end

c     ****************************************
c     *                                      *
c     *        Dneall_4m_push_get_block      *
c     *                                      *
c     ****************************************
      logical function Dneall_4m_push_get_block(mb,nb,hml)
      implicit none
      integer mb,nb
      integer hml(2)

#include "Dne.fh"
#include "bafdecls.fh"

      integer size

      if (mparallelized) then
         if (mb.eq.0) then
            size = m2cq(1)*n2cq(1) + m2cq(2)*n2cq(2)
         else
            size = m2cq(mb)*n2cq(mb)
         end if
      else
         if (mb.eq.0) then
            size = 4*(ne(1)*ne(1) + ne(2)*ne(2))
         else
            size = 4*ne(mb)*ne(mb)
         end if
      end if

      Dneall_4m_push_get_block 
     >  = BA_push_get(mt_dbl,nb*size,'hmlsb',hml(2),hml(1))
      return
      end

c     ****************************************
c     *                                      *
c     *        Dneall_AR_to_4m               *
c     *                                      *
c     ****************************************
      subroutine Dneall_AR_to_4m(mb,A,R,T)
      implicit none
      integer mb
      real*8 A(*),R(*),T(*)

#include "bafdecls.fh"           
#include "errquit.fh"
#include "Dne.fh"

*     **** local variables ****
      integer ms,ms1,ms2,ishift1,ishift2,shift1,shift2

      if (mparallelized) then
         call errquit('Dneall_AR_to_4m: mparallel not finished',0,0)
      else
        if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift1 = ne(1)*ne(1)
            ishift2 = 4*ne(1)*ne(1)
         else
            ms1 = mb
            ms2 = mb
            ishift1 = 0
            ishift2 = 0
         end if
         call dcopy(4*(ne(1)*ne(1)+ne(2)*ne(2)),0.0d0,0,T,1)
         do ms=ms1,ms2
             shift1 = 1+(ms-1)*ishift1
             shift2 = 1+(ms-1)*ishift2
             call Dneall_AR_to_4msub(ne(ms),
     >                               A(shift1),
     >                               R(shift1),
     >                               T(shift2))
         end do
      end if
      return
      end

      subroutine Dneall_AR_to_4msub(n,A,R,T)
      implicit none
      integer n
      real*8 A(n,n)
      real*8 R(n,n)
      real*8 T(2*n,2*n)

*     **** local variables ****
      integer i,j

*     **** copy A to upper-left of T ****
      do j=1,n
      do i=1,n
         T(i,j) = A(i,j)
      end do
      end do

*     **** copy R to lower-left of T ****
      do j=1,n
      do i=1,n
         T(i+n,j) = R(i,j)
      end do
      end do

*     **** copy -R^t to upper-right of T ****
      do j=1,n
      do i=1,n
         T(i,j+n) = -R(j,i)
      end do
      end do

      return
      end

c     ****************************************
c     *                                      *
c     *        Dneall_4m_to_MN               *
c     *                                      *
c     ****************************************
      subroutine Dneall_4m_to_MN(mb,R,M,N)
      implicit none
      integer mb
      real*8 R(*),M(*),N(*)

#include "bafdecls.fh"           
#include "errquit.fh"
#include "Dne.fh"

*     **** local variables ****
      integer ms,ms1,ms2,ishift1,ishift2,shift1,shift2

      if (mparallelized) then
         call errquit('Dneall_AR_to_4m: mparallel not finished',0,0)
      else
        if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift1 = ne(1)*ne(1)
            ishift2 = 4*ne(1)*ne(1)
         else
            ms1 = mb
            ms2 = mb
            ishift1 = 0
            ishift2 = 0
         end if
         do ms=ms1,ms2
            if (ne(ms).le.0) go to 30
            shift1 = 1 + (ms-1)*ishift1
            shift2 = 1 + (ms-1)*ishift2
            call Dneall_4m_MN_sub(ne(ms),R(shift2),M(shift1),N(shift1))
  30        continue
         end do
      end if

      return
      end
      subroutine Dneall_4m_MN_sub(n,R,X,Y)
      implicit none
      integer n
      real*8 R(2*n,2*n)
      real*8 X(n,n)
      real*8 Y(n,n)
      integer i,j
      do j=1,n
      do i=1,n
         X(i,j) = R(i,j)
         Y(i,j) = R(i+n,j)
      end do
      end do

      return
      end


c     ****************************************
c     *                                      *
c     *        Dneall_4m_FactorSkew          *
c     *                                      *
c     ****************************************

*  This routine factors a 4m skew matrix K such that
*  
*  K = U*Sigma*U^H, where U = (V+i*W)
*  
*     Entry - mb:  
*             K: 4mskew matrix - destroyed on output
*     Exit - V: real factor matrix 
*            W: imaginary factor matrix
*            Sigma: eigenvalues
*
      subroutine Dneall_4m_FactorSkew(mb,K,V,W,Sigma)
      implicit none
      integer    mb
      real*8     K(*),V(*),W(*),Sigma(*)

#include "bafdecls.fh"           
#include "errquit.fh"
#include "Dne.fh"

*     **** local variables ****
      integer ms,ms1,ms2,shift1,shift2,ishift1,ishift2

      call nwpw_timing_start(17)

      if (mparallelized) then
        call errquit('Dneall_4m_FactorSkew: mparallel not finished',0,0)
      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift1 = 2*ne(1)
            ishift2 = 4*ne(1)*ne(1)
         else
            ms1 = mb
            ms2 = mb
            ishift1 = 0
            ishift2 = 0
         end if
         do ms=ms1,ms2
            shift1 = 1 + (ms-1)*ishift1
            shift2 = 1 + (ms-1)*ishift2

            call Factor_Skew(2*ne(ms),
     >                       K(shift2),
     >                       V(shift2),
     >                       W(shift2),
     >                       Sigma(shift1))
         end do
      end if

      call nwpw_timing_end(17)
      return
      end


c     ****************************************
c     *                                      *
c     *        Dneall_4mm_transpose          *
c     *                                      *
c     ****************************************

      subroutine Dneall_4mm_transpose(mb,Min,Mout)
      implicit none
      integer mb
      real*8  Min(*),Mout(*)


#include "bafdecls.fh"
#include "Dne.fh"

*     **** local variables ****
      integer ms,ms1,ms2,shift2,ishift2

      if (mparallelized) then
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = m2cq(1)*n2cq(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if

         do ms=ms1,ms2
           if (ne(ms).le.0) goto 20
           shift2 = 1 + (ms-1)*ishift2
           call DMatrix_mm_transpose(2*ne(ms),
     >                 Min(shift2),Mout(shift2),
     >                 m2cq(ms),int_mb(m2c(1,ms)),int_mb(n2c(1,ms)))
 20        continue
         end do

      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = 4*ne(1)*ne(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if

         do ms=ms1,ms2
           if (ne(ms).le.0) goto 30
           shift2 = 1 + (ms-1)*ishift2
           call Dneall_mm_transpose_sub(2*ne(ms),
     >                 Min(shift2),Mout(shift2))
 30        continue
         end do
      end if

      return
      end


c     ****************************************
c     *                                      *
c     *        Dneall_4m_RotateSkew          *
c     *                                      *
c     ****************************************

      subroutine Dneall_4m_RotateSkew(mb,t,V,W,Sigma,R)
      implicit none
      integer mb
      real*8 t
      real*8 V(*),W(*),Sigma(*),R(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "Dne.fh"

*     **** local variables ****
      logical value
      integer ms,ms1,ms2,shift1,shift2,ishift1,ishift2,k,j,indx1,indx2
      integer nj
      integer A(2),B(2),SA(2),SB(2)

*     **** external functions ****
      logical  Dneall_4m_push_get,Dneall_m_pop_stack
      external Dneall_4m_push_get,Dneall_m_pop_stack

      if (mparallelized) then
        call errquit('Dneall_4mm_RotateSkew:mparallel not finished',0,0)
      
      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift1 = 2*ne(1)
            ishift2 = 4*ne(1)*ne(1)
            nj = 2*(ne(1)+ne(2))
         else
            ms1 = mb
            ms2 = mb
            ishift1 = 0
            ishift2 = 0
            nj = 2*ne(mb)
         end if
         value =           Dneall_4m_push_get(mb,A)
         value = value.and.Dneall_4m_push_get(mb,B)
         value = value.and.BA_push_get(mt_dbl,nj,'SA',SA(2),SA(1))
         value = value.and.BA_push_get(mt_dbl,nj,'SB',SB(2),SB(1))
         if (.not. value) 
     >   call errquit('Dneall_4mm_RotateSkew: pushing stack',0,MA_ERR)


         call Dneall_4mm_RotateSkew_sub1(nj,t,
     >                                   Sigma,
     >                                   dbl_mb(SA(1)),
     >                                   dbl_mb(SB(1)))
         do ms=ms1,ms2
           if (ne(ms).le.0) go to 30
           shift1 = (ms-1)*ishift1
           shift2 = (ms-1)*ishift2
           call Dneall_4mm_RotateSkew_sub2(2*ne(ms),
     >                                     dbl_mb(SA(1)+shift1),
     >                                     dbl_mb(SB(1)+shift1),
     >                                     V(shift2+1),
     >                                     W(shift2+1),
     >                                     dbl_mb(A(1)+shift2),
     >                                     dbl_mb(B(1)+shift2))
 30        continue
         end do
         call Dneall_4mmm_Multiply3(mb,V,dbl_mb(A(1)),1.0d0,R,0.0d0)
         call Dneall_4mmm_Multiply3(mb,W,dbl_mb(B(1)),1.0d0,R,1.0d0)

         value =           BA_pop_stack(SB(2))
         value = value.and.BA_pop_stack(SA(2))
         value = value.and.Dneall_m_pop_stack(B)
         value = value.and.Dneall_m_pop_stack(A)
         if (.not. value) 
     >   call errquit('Dneall_4mm_FactorSkew:popping stack',0,MA_ERR)

      end if

      return
      end

      subroutine Dneall_4mm_RotateSkew_sub1(N,t,Sigma,SA,SB)
      implicit none
      integer N
      real*8 t
      real*8 Sigma(N),SA(N),SB(N)

      integer i
      do i=1,N
         SA(i) = dcos(Sigma(i)*t)
         SB(i) = dsin(Sigma(i)*t)
      end do
      return
      end

      subroutine Dneall_4mm_RotateSkew_sub2(N,SA,SB,V,W,A,B)
      implicit none
      integer N
      real*8 SA(N),SB(N)
      real*8 V(N,N),W(N,N)
      real*8 A(N,N),B(N,N)

      integer i,j

      do j=1,N
         do i=1,N
            A(i,j) = V(i,j)*SA(j) + W(i,j)*SB(j)
            B(i,j) = W(i,j)*SA(j) - V(i,j)*SB(j)
         end do
      end do
      return
      end



c     ****************************************
c     *                                      *
c     *        Dneall_4mmm_Multiply2         *
c     *                                      *
c     ****************************************

      subroutine Dneall_4mmm_Multiply2(mb,A,B,C)
      implicit none
      integer mb
      real*8 A(*),B(*),C(*)
           
#include "bafdecls.fh"
#include "Dne.fh"

*     **** local variables ****
      integer MASTER
      parameter (MASTER=0)
      integer ms,ms1,ms2,n,shift2,ishift2
  
      if (mparallelized) then
        if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = m2cq(1)*n2cq(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if

         do ms=ms1,ms2
            n = 2*ne(ms)
            if (n.le.0) go to 20
            shift2 = 1 + (ms-1)*ishift2
            call DMatrix_dgemm2(n,n,n,64,
     >             1.0d0,
     >             A(shift2),m2cq(ms),
     >                       int_mb(m2c(1,ms)),
     >                       int_mb(n2c(1,ms)),
     >             B(shift2),m2cq(ms),
     >                       int_mb(m2c(1,ms)),
     >                       int_mb(n2c(1,ms)),
     >             0.0d0,
     >             C(shift2),m2cq(ms),
     >                       int_mb(m2c(1,ms)),
     >                       int_mb(n2c(1,ms)),
     >             taskid_i,taskid_j,
     >             np_i,np_j,
     >             comm_i, comm_j,
     >             dbl_mb(work1(1)),dbl_mb(work2(1)))

   20       continue
         end do


      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = 4*ne(1)*ne(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if

         do ms=ms1,ms2
            n     = 2*ne(ms)
            if (n.le.0) go to 30
            shift2 = 1 + (ms-1)*ishift2
            call DGEMM('T','N',n,n,n,1.0d0,
     >                A(shift2), n,
     >                B(shift2), n,
     >                0.0d0,
     >                C(shift2), n)
   30       continue
         end do
         if (mb.eq.0) then
            call Parallel_Brdcst_values(MASTER,
     >                   4*(ne(1)*ne(1)+ne(2)*ne(2)),C)
         else
            call Parallel_Brdcst_values(MASTER,4*ne(mb)*ne(mb),C)
         end if

      end if

      return
      end


c     ****************************************
c     *                                      *
c     *        Dneall_4mmm_Multiply3         *
c     *                                      *
c     ****************************************

      subroutine Dneall_4mmm_Multiply3(mb,A,B,alpha,C,beta)
      implicit none
      integer mb
      real*8 A(*),B(*),C(*)
      real*8 alpha,beta

#include "bafdecls.fh"
#include "Dne.fh"

*     **** local variables ****
      integer MASTER
      parameter (MASTER=0)
      integer ms,ms1,ms2,n,shift2,ishift2

      if (mparallelized) then
        if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = m2cq(1)*n2cq(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if

         do ms=ms1,ms2
            n = 2*ne(ms)
            if (n.le.0) go to 20
            shift2 = 1 + (ms-1)*ishift2
            call DMatrix_dgemm3(n,n,n,64,
     >             alpha,
     >             A(shift2),m2cq(ms),
     >                       int_mb(m2c(1,ms)),
     >                       int_mb(n2c(1,ms)),
     >             B(shift2),m2cq(ms),
     >                       int_mb(m2c(1,ms)),
     >                       int_mb(n2c(1,ms)),
     >             beta,
     >             C(shift2),m2cq(ms),
     >                       int_mb(m2c(1,ms)),
     >                       int_mb(n2c(1,ms)),
     >             taskid_i,taskid_j,
     >             np_i,np_j,
     >             comm_i, comm_j,
     >             dbl_mb(work1(1)),dbl_mb(work2(1)))

   20       continue
         end do

      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = 4*ne(1)*ne(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if

         do ms=ms1,ms2
            n     = 2*ne(ms)
            if (n.le.0) go to 30
            shift2 = 1 + (ms-1)*ishift2
            call DGEMM('N','T',n,n,n,alpha,
     >                A(shift2), n,
     >                B(shift2), n,
     >                beta,
     >                C(shift2), n)
   30       continue
         end do
         if (mb.eq.0) then
            call Parallel_Brdcst_values(MASTER,
     >                   4*(ne(1)*ne(1)+ne(2)*ne(2)),C)
         else
            call Parallel_Brdcst_values(MASTER,4*ne(mb)*ne(mb),C)
         end if

      end if 

      return
      end





c     ****************************************
c     *                                      *
c     *        Dneall_4mmm_Multiply          *
c     *                                      *
c     ****************************************

      subroutine Dneall_4mmm_Multiply(mb,A,B,alpha,C,beta)
      implicit none
      integer mb
      real*8 A(*),B(*),C(*)
      real*8 alpha,beta
           
#include "bafdecls.fh"
#include "Dne.fh"

*     **** local variables ****
      integer MASTER
      parameter (MASTER=0)
      integer ms,ms1,ms2,n,shift2,ishift2
  
      if (mparallelized) then
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = m2cq(1)*n2cq(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if

         do ms=ms1,ms2
            n = 2*ne(ms)
            if (n.le.0) go to 20
            shift2 = 1 + (ms-1)*ishift2
            call DMatrix_dgemm1(n,n,n,64,
     >             alpha,
     >             A(shift2),m2cq(ms),
     >                       int_mb(m2c(1,ms)),
     >                       int_mb(n2c(1,ms)),
     >             B(shift2),m2cq(ms),
     >                       int_mb(m2c(1,ms)),
     >                       int_mb(n2c(1,ms)),
     >             beta,
     >             C(shift2),m2cq(ms),
     >                       int_mb(m2c(1,ms)),
     >                       int_mb(n2c(1,ms)),
     >             taskid_i,taskid_j,
     >             np_i,np_j,
     >             comm_i, comm_j,
     >             dbl_mb(work1(1)),dbl_mb(work2(1)))

   20       continue
         end do

      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = 4*ne(1)*ne(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if

         do ms=ms1,ms2
            n     = 2*ne(ms)
            if (n.le.0) go to 30
            shift2 = 1 + (ms-1)*ishift2
            call DGEMM('N','N',n,n,n,
     >                alpha,
     >                A(shift2), n,
     >                B(shift2), n,
     >                beta,
     >                C(shift2), n)
   30       continue
         end do
         if (mb.eq.0) then
            call Parallel_Brdcst_values(MASTER,
     >                   4*(ne(1)*ne(1)+ne(2)*ne(2)),C)
         else
            call Parallel_Brdcst_values(MASTER,4*ne(mb)*ne(mb),C)
         end if

      end if

      return
      end






c     ****************************************
c     *                                      *
c     *           Dneall_4m_trace            *
c     *                                      *
c     ****************************************

      double precision function Dneall_4m_trace(mb,M)
      implicit none
      integer mb
      real*8  M(*)

#include "bafdecls.fh"
#include "Dne.fh"

*     **** local variables ****
      integer ms,ms1,ms2,shift2,ishift2,indx,i
      real*8  sum

*     ***** external functions ****
      real*8   DMatrix_trace
      external DMatrix_trace

      if (mparallelized) then
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = m2cq(1)*n2cq(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if
         sum = 0.0d0
         do ms=ms1,ms2
            shift2 = 1 + (ms-1)*ishift2
            sum =  sum + DMatrix_trace(2*ne(ms),
     >                                 M(shift2),m2cq(ms),
     >                                           int_mb(m2c(1,ms)),
     >                                           int_mb(n2c(1,ms)),
     >                                 taskid_i,taskid_j)
         end do

      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = 4*ne(1)*ne(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if

         sum = 0.0d0
         do ms=ms1,ms2

           if (ne(ms).le.0) goto 30
           shift2 = (ms-1)*ishift2

           do i=1,ne(ms)
              indx = i + (i-1)*2*ne(ms) + shift2
              sum = sum + M(indx)
           end do

 30        continue
         end do
      end if

      Dneall_4m_trace = sum
      return
      end


c     ****************************************
c     *                                      *
c     *           Dneall_4m_scal             *
c     *                                      *
c     ****************************************

      subroutine Dneall_4m_scal(mb,alpha,M)
      implicit none
      integer mb
      real*8  alpha
      real*8  M(*)

#include "Dne.fh"

*     **** local variables ****
      integer nn

      if (mparallelized) then
         if (mb.eq.0) then
            nn = m2cq(1)*n2cq(1)+m2cq(2)*n2cq(2)
         else
            nn = m2cq(mb)*n2cq(mb)
         end if
      else
         if (mb.eq.0) then
            nn = 4*(ne(1)*ne(1) + ne(2)*ne(2))
         else
            nn = 4*ne(mb)*ne(mb)
         end if
      end if

      call dscal(nn,alpha,M,1)
      return
      end



c     ****************************************
c     *                                      *
c     *           Dneall_4m_diag_scal        *
c     *                                      *
c     ****************************************

      subroutine Dneall_4m_diag_scal(mb,diag,M)
      implicit none
      integer mb
      real*8  diag(*)
      real*8  M(*)

#include "Dne.fh"

*     **** local variables ****
      integer ms,ms1,ms2,shift1,shift2,ishift1,ishift2,indx,i

      if (mparallelized) then
         call errquit('Dneall_4m_diag_scal not finished',0,0)
      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift1 = 2*ne(1)
            ishift2 = 4*ne(1)*ne(1)
         else
            ms1 = mb
            ms2 = mb
            ishift1 = 0
            ishift2 = 0
         end if

         do ms=ms1,ms2

            if (ne(ms).le.0) goto 30
            shift1 = (ms-1)*ishift1
            shift2 = (ms-1)*ishift2

            do i=1,2*ne(ms)
               indx = i + (i-1)*2*ne(ms) + shift2
               M(indx) = M(indx)*diag(i+shift1)
            end do

 30        continue
         end do
      end if

      return
      end

c     ****************************************
c     *                                      *
c     *        Dneall_4m_diag_scal_inv       *
c     *                                      *
c     ****************************************

      subroutine Dneall_4m_diag_scal_inv(mb,diag,M)
      implicit none
      integer mb
      real*8  diag(*)
      real*8  M(*)

#include "Dne.fh"

*     **** local variables ****
      integer ms,ms1,ms2,shift1,shift2,ishift1,ishift2,indx,i

      if (mparallelized) then
         call errquit('Dneall_4m_diag_scal_inv not finished',0,0)
      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift1 = 2*ne(1)
            ishift2 = 4*ne(1)*ne(1)
         else
            ms1 = mb
            ms2 = mb
            ishift1 = 0
            ishift2 = 0
         end if

         do ms=ms1,ms2

            if (ne(ms).le.0) goto 30
            shift1 = (ms-1)*ishift1
            shift2 = (ms-1)*ishift2

            do i=1,2*ne(ms)
               indx = i + (i-1)*2*ne(ms) + shift2
               M(indx) = M(indx)/diag(i+shift1)
            end do

 30        continue
         end do
      end if
      return
      end


      real*8 function  Dneall_4m_value(mb,ms,i,j,A)
      implicit none
      integer mb,ms,i,j
      real*8  A(*)

#include "bafdecls.fh"
#include "Dne.fh"

*     **** local variables ****
      integer ishift2,shift2
      real*8 w

*     **** external functions ****
      real*8   DMatrix_m_get_value
      external DMatrix_m_get_value
      
      w = 0.0d0
      if (ne(ms).gt.0)  then
         if (mparallelized) then

            if (mb.eq.0) then
               ishift2 = m2cq(1)*n2cq(1)
            else
               ishift2 = 0
            end if
            shift2 = 1 + (ms-1)*ishift2
            w = DMatrix_m_get_value(i,j,A(shift2),
     >                 m2cq(ms),int_mb(m2c(1,ms)),int_mb(n2c(1,ms)))
         else
            shift2 = (ms-1)*ne(1)*2
            w = A(i+(j-1)*2*ne(ms)+shift2)
         end if
      end if

      Dneall_4m_value = w
      return
      end


c     ****************************************
c     *                                      *
c     *        Dneall_4m_eye                  *
c     *                                      *
c     ****************************************

      subroutine Dneall_4m_eye(mb,A,alpha)
      implicit none
      integer mb
      real*8 A(*)
      real*8 alpha

#include "bafdecls.fh"
#include "Dne.fh"


*     **** local variables ****
      integer ms,ms1,ms2,shift2,ishift2,k
      integer indx2

      if (mparallelized) then
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = m2cq(1)*n2cq(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if

         do ms=ms1,ms2
             shift2 = 1 + (ms-1)*ishift2
             call DMatrix_eye(2*ne(ms),2*ne(ms),
     >                A(shift2),
     >                m2cq(ms),int_mb(m2c(1,ms)),int_mb(n2c(1,ms)),
     >                taskid_i,taskid_j)
         end do

      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = 4*ne(1)*ne(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if

         do ms=ms1,ms2
           if (ne(ms).le.0) go to 30
           shift2 = 1 + (ms-1)*ishift2
           call dcopy(4*ne(ms)*ne(ms),0.0d0,0,A(shift2),1)
           indx2 = shift2
           do k=1,2*ne(ms)
              A(indx2) = alpha
              indx2    = indx2 + (1+2*ne(ms))
           end do

 30        continue
         end do
      end if

      return
      end


c     ****************************************
c     *                                      *
c     *        Dneall_mmm_Multiply2ab        *
c     *                                      *
c     ****************************************

      subroutine Dneall_mmm_Multiply2ab(mb,A,B,alpha,C,beta)
      implicit none
      integer mb
      real*8 A(*),B(*),C(*)
      real*8 alpha,beta
           
#include "bafdecls.fh"
#include "Dne.fh"

*     **** local variables ****
      integer MASTER
      parameter (MASTER=0)
      integer ms,ms1,ms2,n,shift2,ishift2
  
      if (mparallelized) then
        if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = mcq(1)*ncq(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if

         do ms=ms1,ms2
            if (ne(ms).le.0) go to 20
            shift2 = 1 + (ms-1)*ishift2
            call DMatrix_dgemm2(ne(ms),ne(ms),ne(ms),64,
     >             alpha,
     >             A(shift2),mcq(ms),
     >                       int_mb(mc(1,ms)),
     >                       int_mb(nc(1,ms)),
     >             B(shift2),mcq(ms),
     >                       int_mb(mc(1,ms)),
     >                       int_mb(nc(1,ms)),
     >             beta,
     >             C(shift2),mcq(ms),
     >                       int_mb(mc(1,ms)),
     >                       int_mb(nc(1,ms)),
     >             taskid_i,taskid_j,
     >             np_i,np_j,
     >             comm_i, comm_j,
     >             dbl_mb(work1(1)),dbl_mb(work2(1)))

   20       continue
         end do


      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = ne(1)*ne(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if

         do ms=ms1,ms2
            n     = ne(ms)
            if (n.le.0) go to 30
            shift2 = 1 + (ms-1)*ishift2
            call DGEMM('T','N',n,n,n,alpha,
     >                A(shift2), n,
     >                B(shift2), n,
     >                beta,
     >                C(shift2), n)
   30       continue
         end do
         if (mb.eq.0) then
            call Parallel_Brdcst_values(MASTER,
     >                   ne(1)*ne(1)+ne(2)*ne(2),C)
         else
            call Parallel_Brdcst_values(MASTER,ne(mb)*ne(mb),C)
         end if

      end if

      return
      end

c     ****************************************
c     *                                      *
c     *        Dneall_wtom_Real              *
c     *                                      *
c     ****************************************
      subroutine Dneall_wtom_Real(mb,A,B)
      implicit none
      integer mb
      complex*16 A(*)
      real*8     B(*)
      integer i,nsize
      call Dneall_m_size(mb,nsize)
!$OMP DO
      do i=1,nsize
         B(i) = dble(A(i))
      end do
!$OMP END DO
      return
      end


c     ****************************************
c     *                                      *
c     *        Dneall_mmtow_Cmplx            *
c     *                                      *
c     ****************************************
      subroutine Dneall_mmtow_Cmplx(mb,Ar,Br,C)
      implicit none
      integer mb
      real*8 Ar(*),Br(*)
      complex*16 C(*)
      integer i,nsize
      call Dneall_m_size(mb,nsize)
!$OMP DO
      do i=1,nsize
         C(i) = dcmplx(Ar(i),Br(i))
      end do
!$OMP END DO
      return
      end


        
c     ****************************************
c     *                                      *
c     *        Dne_ffm_combo_sym_Multiply    *
c     *                                      *
c     ****************************************
c
c   warning - ms==0 not allowed

      subroutine Dne_ffm_combo_sym_Multiply(ms,A1,A2,npack1,hml3)
      implicit none
      complex*16 A1(*),A2(*)
      integer    ms,npack1
      real*8     hml3(*)
           
#include "bafdecls.fh"
#include "Dne.fh"

*     **** local variables ****
      integer n,shift
      integer tid,nthr,nn

*     **** external functions ****
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads

      call nwpw_timing_start(15)
      tid  = Parallel_threadid()
      nthr = Parallel_nthreads()

      if (parallelized) then
         nn = mcqmax(ms)*ncqmax(ms)

         if (mparallelized) then
            if (ne(ms).le.0) go to 20
            shift  = 1 + (ms-1)*neq(1)*npack1
            
            call DMatrix_combo_dgemm2(ne(ms),ne(ms),npack1_all,128,
     >             2.0d0,
     >             A1(shift),A2(shift),
     >                       int_mb(ma(1,ms)+taskid_i), 
     >                       int_mb(ma(1,ms)),
     >                       int_mb(na(1,ms)),
     >             0.0d0,
     >             hml3,        int_mb(mc(1,ms)+taskid_i), 
     >                         int_mb(nc(1,ms)+taskid_j),
     >                         int_mb(mc(1,ms)),
     >                         int_mb(nc(1,ms)),
     >             taskid_i,taskid_j,
     >             np_i,np_j,
     >             comm_i, comm_j,
     >             dbl_mb(work1(1)),dbl_mb(work2(1)))

            call DMatrix_combo_dgemm2(ne(ms),ne(ms),nida1_all,128,
     >             -1.0d0,
     >             A1(shift),A2(shift),
     >                       int_mb(ma(1,ms)+taskid_i), 
     >                       int_mb(ma1(1,ms)),
     >                       int_mb(na(1,ms)),
     >             1.0d0,
     >             hml3,        int_mb(mc(1,ms)+taskid_i), 
     >                         int_mb(nc(1,ms)+taskid_j),
     >                         int_mb(mc(1,ms)),
     >                         int_mb(nc(1,ms)),
     >             taskid_i,taskid_j,
     >             np_i,np_j,
     >             comm_i, comm_j,
     >             dbl_mb(work1(1)),dbl_mb(work2(1)))

  20        continue

         else

            if (ne(ms).le.0) go to 21
            shift  = 1 + (ms-1)*neq(1)*npack1
            
            call Parallel_shared_vector_zero(.false.,3*mall(ms),hml3)
c            call DMatrix_combo_dgemm2(ne(ms),ne(ms),npack1_all,128,
c     >             2.0d0,
c     >             A1(shift),A2(shift),
c     >                       int_mb(ma(1,ms)+taskid_i),
c     >                       int_mb(ma(1,ms)),
c     >                       int_mb(na(1,ms)),
c     >             0.0d0,
c     >             dbl_mb(mat_tmp(1)),
c     >                         int_mb(mc(1,ms)+taskid_i),
c     >                         int_mb(nc(1,ms)+taskid_j),
c     >                         int_mb(mc(1,ms)),
c     >                         int_mb(nc(1,ms)),
c     >             taskid_i,taskid_j,
c     >             np_i,np_j,
c     >             comm_i, comm_j,
c     >             dbl_mb(work1(1)),dbl_mb(work2(1)))
c
c            call DMatrix_combo_dgemm2(ne(ms),ne(ms),nida1_all,128,
c     >             -1.0d0,
c     >             A1(shift),A2(shift),
c     >                       int_mb(ma(1,ms)+taskid_i),
c     >                       int_mb(ma1(1,ms)),
c     >                       int_mb(na(1,ms)),
c     >             1.0d0,
c     >             dbl_mb(mat_tmp(1)),
c     >                         int_mb(mc(1,ms)+taskid_i),
c     >                         int_mb(nc(1,ms)+taskid_j),
c     >                         int_mb(mc(1,ms)),
c     >                         int_mb(nc(1,ms)),
c     >             taskid_i,taskid_j,
c     >             np_i,np_j,
c     >             comm_i, comm_j,
c     >             dbl_mb(work1(1)),dbl_mb(work2(1)))
            call DMatrix_combo_dgemm2c_omp(ne(ms),ne(ms),npack1_all,128,
     >             A1(shift),A2(shift),int_mb(ma(1,ms)+taskid_i),
     >                                 int_mb(ma(1,ms)),
     >                                 int_mb(ma1(1,ms)),
     >                                 int_mb(na(1,ms)),
     >             dbl_mb(mat_tmp(1)),int_mb(mc(1,ms)+taskid_i),
     >                                int_mb(nc(1,ms)+taskid_j),
     >                                int_mb(mc(1,ms)),
     >                                int_mb(nc(1,ms)),
     >             taskid_i,taskid_j,
     >             np_i,np_j,
     >             comm_i, comm_j,
     >             tid,nthr,
     >             dbl_mb(work1(1)),dbl_mb(work2(1)),
     >             dbl_mb(thrwork1(1)+tid*nn))

  21        continue
            call Dne_m_combo_gather(mall(ms),
     >                        mpack(ms),int_mb(mindx(1,ms)),
     >                        dbl_mb(mat_tmp(1)),hml3)

         end if
      else

         nn = ne(ms)*ne(ms)
         if (ne(ms).le.0) go to 30
         shift  = 1 + (ms-1)*ne(1)*npack1
c         call Pack_ccm_combo_sym_dot(1,ne(ms),
c     >                     A1(shift),
c     >                     A2(shift),
c     >                     hml3)
         call Pack_ccm_combo_sym_dot_omp(1,ne(ms),
     >                     A1(shift),
     >                     A2(shift),
     >                     hml3,dbl_mb(thrwork1(1)+tid*nn))

  30     continue
      end if
      call nwpw_timing_end(15)
      return
      end

      subroutine Dne_m_combo_gather(nall,npack,indx,A,B)
      implicit none
      integer nall,npack,indx(*)
      real*8 A(npack,*),B(nall,*)
      integer i
      !call dcopy(3*nall,0.0d0,0,B,1)
!$OMP DO
      do i=1,npack
         B(indx(i),1) = A(i,1)
         B(indx(i),2) = A(i,2)
         B(indx(i),3) = A(i,3)
      end do  
!$OMP END DO
      call Parallel_Vector_SumAll(3*nall,B)
      return
      end








c     ****************************************
c     *                                      *
c     *        Dneall_f_GramSchmidt          *
c     *                                      *
c     ****************************************

c   performs a modified gramSchmidt on Q
c
      subroutine Dneall_f_GramSchmidt(mb,Q,npack1)
      implicit none
      integer mb
      complex*16 Q(*)
      integer    npack1      

#include "bafdecls.fh"           
#include "errquit.fh"
#include "Dne.fh"

c     **** local variables ****
      integer ms,ms1,ms2,shift,j,k,indxk,indxj
      integer jj,kk,jcur,kcur
      integer tmp(2)

      real*8 w
      common /Dneall_f_Sortho_tmp/ w

      if (parallelized) then

c        **** allocate temporary memory ****
         if (.not.BA_push_get(mt_dcpl,npack1,'tmp',tmp(2),tmp(1)))
     >   call errquit('Dneall_f_Gramschmidt:out of stack',0,MA_ERR)

         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
         else
            ms1 = mb
            ms2 = mb
         end if
         shift = neq(1)*npack1

         do ms=ms1,ms2

            kcur = np_j-1
            kk   = int_mb(na(1,ms)+kcur) - 1 

            do k=ne(ms),1,-1

               if (kcur.eq.taskid_j) then
                  indxk = 1 + kk*npack1    + (ms-1)*shift
                  call Pack_cc_dot(1,Q(indxk),Q(indxk),w)
!$OMP SINGLE
                  w        = 1.0d0/dsqrt(w)
!$OMP END SINGLE
                  call Pack_c_SMul1(1,w,Q(indxk))
                  !call dcopy(2*npack1,Q(indxk),1,dcpl_mb(tmp(1)),1)
                  call Parallel_shared_vector_copy(.true.,
     >                                             2*npack1,Q(indxk),
     >                                             dcpl_mb(tmp(1)))
               end if

*              *** Brdcst Q(indxk) ***
               if (kcur.gt.0) 
     >            call D1dB_Brdcst_values(kcur,2*npack1,dcpl_mb(tmp(1)))

c              *** set j = k+1 ***
               jj   = kk
               jcur = kcur

               jj = jj - 1
               if (jj.lt.0) then
                  jcur = jcur - 1
                  jj = int_mb(na(1,ms) + jcur) - 1
               end if

               do j=k-1,1,-1
                  if (jcur.eq.taskid_j) then
                     indxj = 1 + jj*npack1    + (ms-1)*shift
                     call Pack_cc_dot(1,dcpl_mb(tmp(1)),Q(indxj),w)
!$OMP SINGLE
                     w = -w
!$OMP END SINGLE
                     call Pack_cc_daxpy(1,w,dcpl_mb(tmp(1)),Q(indxj))
                  end if

                  jj = jj - 1
                  if (jj.lt.0) then
                     jcur = jcur - 1
                     jj = int_mb(na(1,ms) + jcur) - 1
                  end if

               end do

               kk = kk - 1
               if (kk.lt.0) then
                  kcur = kcur - 1
                  kk = int_mb(na(1,ms) + kcur) - 1
               end if

            end do
         end do

c        **** deallocate local memory ****
         if (.not.BA_pop_stack(tmp(2)))
     >   call errquit('Dneall_f_GramSchmidt:popping stack',0,MA_ERR)



c     **** not parallized ****
      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
         else
            ms1 = mb
            ms2 = mb
         end if
         shift = ne(1)*npack1

c         **** modified Gram-Schmidt ****
         do ms=ms1,ms2
         do k=ne(ms),1,-1
            indxk = 1 + (k-1)*npack1 + (ms-1)*shift
            call Pack_cc_dot(1,Q(indxk),Q(indxk),w)
!$OMP MASTER
            w        = 1.0d0/dsqrt(w)
!$OMP END MASTER
!$OMP BARRIER
            call Pack_c_SMul1(1,w,Q(indxk))

            do j=k-1,1,-1
               indxj = 1 + (j-1)*npack1 + (ms-1)*shift
               call Pack_cc_dot(1,Q(indxk),Q(indxj),w)
!$OMP MASTER
               w = -w
!$OMP END MASTER
!$OMP BARRIER
               call Pack_cc_daxpy(1,w,Q(indxk),Q(indxj))
            end do
         end do
         end do
      end if
!$OMP BARRIER

      return
      end 




c     ****************************************
c     *                                      *
c     *        Dneall_mne_diagonalize        *
c     *                                      *
c     ****************************************

      subroutine Dneall_mne_diagonalize(mb,nein,hml,eig,assending)
      implicit none
      integer    mb,nein(2)
      real*8     hml(*),eig(*)
      logical    assending

#include "bafdecls.fh"           
#include "errquit.fh"
#include "Dne.fh"

*     ***** local variables ****
      integer MASTER
      parameter (MASTER=0)
      logical value
      integer ms,ms1,ms2,shift1,shift2,ishift1,ishift2,ierr
      integer tmp1(2),V(2),VV(2),Q(2),tu(2),ework(2)

      integer mcqin(NBLOCKS),ncqin(NBLOCKS)
      integer mcin(2,NBLOCKS),ncin(2,NBLOCKS)
      integer work1in(2),work2in(2)
      integer i,j,k,nework


      call nwpw_timing_start(17)
      if (mparallelized) then

*        ****************************************************
*        **** define mcin,ncin,mcqin,ncqin based on nein ****
*        ****************************************************
         value = .true.
         mcqin(1) = 0
         mcqin(2) = 0
         ncqin(1) = 0
         ncqin(2) = 0
         do ms=1,ispin
             value = value.and.
     >           BA_push_get(mt_int,np_i,'mcin',mcin(2,ms),mcin(1,ms))
             value = value.and.
     >           BA_push_get(mt_int,np_j,'ncin',ncin(2,ms),ncin(1,ms))
             if (.not.value) then
               call errquit(': out of heap memory',0,MA_ERR)
             end if
             call icopy(np_i,0,0,int_mb(mcin(1,ms)),1)
             call icopy(np_j,0,0,int_mb(ncin(1,ms)),1)
             i = 0
             j = 0
             do k=1,nein(ms)
                int_mb(mcin(1,ms)+i) = int_mb(mcin(1,ms)+i) + 1
                int_mb(ncin(1,ms)+j) = int_mb(ncin(1,ms)+j) + 1
                i = mod(i+1,np_i)
                j = mod(j+1,np_j)
             end do
             mcqin(ms) = int_mb(mcin(1,ms)+taskid_i)
             ncqin(ms) = int_mb(ncin(1,ms)+taskid_j)
         end do
*        ****************************************************
*        **** define mcin,ncin,mcqin,ncqin based on nein ****
*        ****************************************************
*        ***************************************************************
*        ****** allocate work1in and work2in ***************************
*        ***************************************************************
         value = value.and.BA_push_get(mt_dbl,
     >                                2*2*64*int_mb(mcin(1,1)+taskid_i),
     >                                  'work1in',work1in(2),work1in(1))
         nework = 2*64*int_mb(ncin(1,1)+taskid_j)
         if (nework.lt.mcqin(1)*ncqin(1)) nework = mcqin(1)*ncqin(1)
         value = value.and.BA_push_get(mt_dbl,
     >                                  3*nework,
     >                                  'work2in',work2in(2),work2in(1))
         if (.not.value) 
     >      call errquit(': out of heap memory',1,MA_ERR)
*        ***************************************************************
*        ****** allocate work1in and work2in ***************************
*        ***************************************************************


         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift1 = nein(1)
            ishift2 = mcqin(1)*ncqin(1)
            call dcopy((nein(1)+nein(2)),0.0d0,0,eig,1)
         else
            ms1 = mb
            ms2 = mb
            ishift1 = 0
            ishift2 = 0
            call dcopy(nein(mb),0.0d0,0,eig,1)
         end if
         value = BA_push_get(mt_dbl,mcqin(1)*ncqin(1),'V',V(2),V(1))
         value = value.and.
     >           BA_push_get(mt_dbl,mcqin(1)*ncqin(1),'VV',VV(2),VV(1))
         value = value.and.
     >           BA_push_get(mt_dbl,mcqin(1)*ncqin(1),'Q',Q(2),Q(1))
         value = value.and.
     >           BA_push_get(mt_dbl,nein(1),'tu',tu(2),tu(1))
         value = value.and.
     >           BA_push_get(mt_dbl,nein(1),'ework',ework(2),ework(1))
         if (.not. value) 
     >      call errquit('Dneall_m_diagonalize:out of stack',
     >                    0,MA_ERR)
         do ms=ms1,ms2
            shift1 = 1+(ms-1)*ishift1
            shift2 = 1+(ms-1)*ishift2

      
            call nwpw_timing_start(22)
            call DMatrix_tredq(nein(ms),
     >                hml(shift2),dbl_mb(Q(1)),
     >                mcqin(ms),int_mb(mcin(1,ms)),int_mb(ncin(1,ms)),
     >                taskid_i,taskid_j,
     >                np_i,np_j,
     >                comm_i,comm_j,
     >                dbl_mb(work1in(1)),dbl_mb(work2in(1)),
     >                dbl_mb(V(1)),dbl_mb(VV(1)))
            call nwpw_timing_end(22)

            call nwpw_timing_start(23)
            call DMatrix_getdiags(nein(ms),
     >                eig(shift1),dbl_mb(tu(1)),
     >                hml(shift2),
     >                mcqin(ms),int_mb(mcin(1,ms)),int_mb(ncin(1,ms)),
     >                taskid_i,taskid_j,
     >                np_i,np_j,
     >                comm_i,comm_j,
     >                dbl_mb(ework(1)))
            call dcopy(mcqin(ms)*ncqin(ms),dbl_mb(Q(1)),1,hml(shift2),1)
            call nwpw_timing_end(23)

            call nwpw_timing_start(24)
            call DMatrix_tqliq(nein(ms),
     >                eig(shift1),dbl_mb(tu(1)),
     >                hml(shift2),
     >                mcqin(ms),int_mb(mcin(1,ms)),int_mb(ncin(1,ms)),
     >                taskid_i,taskid_j,
     >                np_i,np_j,
     >                comm_i,comm_j,
     >                dbl_mb(work1in(1)),dbl_mb(work2in(1)))
            call nwpw_timing_end(24)

            call nwpw_timing_start(25)
            if (.not.assending)
     >      call DMatrix_eigsrtq(nein(ms),
     >              eig(shift1),
     >              hml(shift2),
     >              mcqin(ms),int_mb(mcin(1,ms)),int_mb(ncin(1,ms)),
     >              taskid_i,taskid_j,
     >              np_i,np_j,
     >              comm_i,comm_j,
     >              dbl_mb(work1in(1)),dbl_mb(work2in(1)))
            call nwpw_timing_end(25)

         end do

         value =           BA_pop_stack(ework(2))
         value = value.and.BA_pop_stack(tu(2))
         value = value.and.BA_pop_stack(Q(2))
         value = value.and.BA_pop_stack(VV(2))
         value = value.and.BA_pop_stack(V(2))
         if (.not. value) 
     >    call errquit('error popping stack in Dneall_mne_diagonalize',
     >                 0,MA_ERR)

*        ****************************************************
*        **** deallocate mcin,ncin,work1in,work2in       ****
*        ****************************************************
         value = .true.
         value = value.and.BA_pop_stack(work2in(2))
         value = value.and.BA_pop_stack(work1in(2))
         do ms=ispin,1,-1
            value = value.and.BA_pop_stack(ncin(2,ms))
            value = value.and.BA_pop_stack(mcin(2,ms))
         end do
        if (.not.value) 
     >     call errquit('Dne_end: error deallocating stack',
     >                  0,MA_ERR)
*        ****************************************************
*        **** deallocate mcin,ncin,work1in,work2in       ****
*        ****************************************************

      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift1 = nein(1)
            ishift2 = nein(1)*ne(1)
            call Parallel_shared_vector_zero(.true.,
     >                            (nein(1)+nein(2)),eig)
         else
            ms1 = mb
            ms2 = mb
            ishift1 = 0
            ishift2 = 0
            call Parallel_shared_vector_zero(.true.,nein(mb),eig)
         end if

         value = BA_push_get(mt_dbl,(2*nein(1)*nein(1)),
     >                       'tmp1',tmp1(2),tmp1(1))
         if (.not. value) 
     >      call errquit('Dneall_mne_diagonalize:out of stack',
     >                    1,MA_ERR)

*        ***** diagonalize the matrix *****
         do ms=ms1,ms2
            shift1 = 1+(ms-1)*ishift1
            shift2 = 1+(ms-1)*ishift2
            if (nein(ms).le.0) go to 30

            ierr = 0
            call DSYEV('V','U',nein(ms),
     >                 hml(shift2),nein(ms), 
     >                 eig(shift1),
     >                 dbl_mb(tmp1(1)),2*nein(1)*nein(1),
     >                 ierr)
            if (ierr.ne.0) then
               ierr = 0
               call EIGEN_UNSORT(nein(ms),nein(ms),
     >                           hml(shift2),eig(shift1),
     >                           dbl_mb(tmp1(1)),ierr)
c              call errquit('dneamdg: dsyev err.ne.0 ',ierr,0)
            end if

            if (.not.assending)
     >       call eigsrt(eig(shift1),
     >                   hml(shift2),
     >                   nein(ms),nein(ms))

  30       continue
         end do
         if (mb.eq.0) then
            call Parallel_Brdcst_values(MASTER,
     >                   nein(1)*nein(1)+nein(2)*nein(2),hml)
         else
            call Parallel_Brdcst_values(MASTER,nein(mb)*nein(mb),hml)
         end if

         value = BA_pop_stack(tmp1(2))
         if (.not. value) 
     >    call errquit('error popping stack in Dneall_mne_diagonalize',
     >                 0,MA_ERR)

      end if

      call nwpw_timing_end(17)
      return
      end


      real*8 function  Dneall_mne_value(mb,nein,ms,i,j,A)
      implicit none
      integer mb,ms,i,j,nein(2)
      real*8  A(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "Dne.fh"

*     **** local variables ****
      integer ishift2,shift2
      real*8 w

      logical value
      integer mcqin(2),ncqin(2)
      integer mcin(2,NBLOCKS),ncin(2,NBLOCKS)
      integer ii,jj,k

*     **** external functions ****
      real*8   DMatrix_m_get_value
      external DMatrix_m_get_value
      
      w = 0.0d0
      if (nein(ms).gt.0)  then
         if (mparallelized) then


         value = .true.
         mcqin(1) = 0
         mcqin(2) = 0
         ncqin(1) = 0
         ncqin(2) = 0
         do ms=1,ispin
             value = value.and.
     >           BA_push_get(mt_int,np_i,'mcin',mcin(2,ms),mcin(1,ms))
             value = value.and.
     >           BA_push_get(mt_int,np_j,'ncin',ncin(2,ms),ncin(1,ms))
             if (.not.value) then
               call errquit(': out of heap memory',0,MA_ERR)
             end if
             call icopy(np_i,0,0,int_mb(mcin(1,ms)),1)
             call icopy(np_j,0,0,int_mb(ncin(1,ms)),1)
             ii = 0
             jj = 0
             do k=1,nein(ms)
                int_mb(mcin(1,ms)+ii) = int_mb(mcin(1,ms)+ii) + 1
                int_mb(ncin(1,ms)+jj) = int_mb(ncin(1,ms)+jj) + 1
                ii = mod(jj+1,np_i)
                jj = mod(jj+1,np_j)
             end do
             mcqin(ms) = int_mb(mcin(1,ms)+taskid_i)
             ncqin(ms) = int_mb(ncin(1,ms)+taskid_j)
         end do


            if (mb.eq.0) then
               ishift2 = mcqin(1)*ncqin(1)
            else
               ishift2 = 0
            end if
            shift2 = 1 + (ms-1)*ishift2
            w = DMatrix_m_get_value(i,j,A(shift2),
     >                 mcqin(ms),int_mb(mcin(1,ms)),int_mb(ncin(1,ms)))


         do ms=ispin,1,-1
            value = value.and.BA_pop_stack(ncin(2,ms))
            value = value.and.BA_pop_stack(mcin(2,ms))
         end do
        if (.not.value) 
     >     call errquit('Dne_end: error deallocating stack',
     >                  0,MA_ERR)



         else
            if (mb.eq.0) then
               ishift2 = nein(1)*nein(1)
            else
               ishift2 = 0
            end if
            shift2 = (ms-1)*ishift2
            w = A(i+(j-1)*nein(ms)+shift2)
         end if
      end if

      Dneall_mne_value = w
      return
      end



      subroutine Dneall_mne_set_value(w,mb,nein,ms,i,j,A)
      implicit none
      real*8 w
      integer mb,ms,i,j,nein(2)
      real*8  A(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "Dne.fh"

*     **** local variables ****
      integer ishift2,shift2

      logical value
      integer mcqin(2),ncqin(2)
      integer mcin(2,NBLOCKS),ncin(2,NBLOCKS)
      integer ii,jj,k

      if (nein(ms).gt.0)  then
         if (mparallelized) then



         value = .true.
         mcqin(1) = 0
         mcqin(2) = 0
         ncqin(1) = 0
         ncqin(2) = 0
         do ms=1,ispin
             value = value.and.
     >           BA_push_get(mt_int,np_i,'mcin',mcin(2,ms),mcin(1,ms))
             value = value.and.
     >           BA_push_get(mt_int,np_j,'ncin',ncin(2,ms),ncin(1,ms))
             if (.not.value) then
               call errquit(': out of heap memory',0,MA_ERR)
             end if
             call icopy(np_i,0,0,int_mb(mcin(1,ms)),1)
             call icopy(np_j,0,0,int_mb(ncin(1,ms)),1)
             ii = 0
             jj = 0
             do k=1,nein(ms)
                int_mb(mcin(1,ms)+ii) = int_mb(mcin(1,ms)+ii) + 1
                int_mb(ncin(1,ms)+jj) = int_mb(ncin(1,ms)+jj) + 1
                ii = mod(ii+1,np_i)
                jj = mod(jj+1,np_j)
             end do
             mcqin(ms) = int_mb(mcin(1,ms)+taskid_i)
             ncqin(ms) = int_mb(ncin(1,ms)+taskid_j)
         end do




            if (mb.eq.0) then
               ishift2 = mcqin(1)*ncqin(1)
            else
               ishift2 = 0
            end if
            shift2 = 1 + (ms-1)*ishift2
            call DMatrix_m_set_value(w,i,j,A(shift2),
     >                 mcqin(ms),int_mb(mcin(1,ms)),int_mb(ncin(1,ms)))



         do ms=ispin,1,-1
            value = value.and.BA_pop_stack(ncin(2,ms))
            value = value.and.BA_pop_stack(mcin(2,ms))
         end do
        if (.not.value)
     >     call errquit('Dne_end: error deallocating stack',
     >                  0,MA_ERR)



         else
            if (mb.eq.0) then
               ishift2 = nein(1)*nein(1)
            else
               ishift2 = 0
            end if
            shift2 = (ms-1)*ishift2
             A(i+(j-1)*nein(ms)+shift2) = w
         end if
      end if

      return
      end



c
c
c     ****************************************
c     *                                      *
c     *        Dneall_mne_allocate_block     *
c     *                                      *
c     ****************************************
      logical function Dneall_mne_allocate_block(mb,nb,nein,hml)
      implicit none
      integer mb,nb,nein(2)
      integer hml(2)

#include "Dne.fh"
#include "errquit.fh"
#include "bafdecls.fh"

      integer size

      logical value
      integer mcqin(2),ncqin(2)
      integer mcin(2,NBLOCKS),ncin(2,NBLOCKS)
      integer ii,jj,k,ms

     
      if (mparallelized) then

         value = .true.
         mcqin(1) = 0
         mcqin(2) = 0
         ncqin(1) = 0
         ncqin(2) = 0
         do ms=1,ispin
             value = value.and.
     >           BA_push_get(mt_int,np_i,'mcin',mcin(2,ms),mcin(1,ms))
             value = value.and.
     >           BA_push_get(mt_int,np_j,'ncin',ncin(2,ms),ncin(1,ms))
             if (.not.value) 
     >         call errquit(': out of stackmemory',0,MA_ERR)
     
             call icopy(np_i,0,0,int_mb(mcin(1,ms)),1)
             call icopy(np_j,0,0,int_mb(ncin(1,ms)),1)
             ii = 0
             jj = 0
             do k=1,nein(ms)
                int_mb(mcin(1,ms)+ii) = int_mb(mcin(1,ms)+ii) + 1
                int_mb(ncin(1,ms)+jj) = int_mb(ncin(1,ms)+jj) + 1
                ii = mod(ii+1,np_i)
                jj = mod(jj+1,np_j)
             end do
             mcqin(ms) = int_mb(mcin(1,ms)+taskid_i)
             ncqin(ms) = int_mb(ncin(1,ms)+taskid_j)
         end do



         if (mb.eq.0) then
            size = mcqin(1)*ncqin(1) + mcqin(2)*ncqin(2)
         else
            size = mcqin(mb)*ncqin(mb)
         end if


         do ms=ispin,1,-1
            value = value.and.BA_pop_stack(ncin(2,ms))
            value = value.and.BA_pop_stack(mcin(2,ms))
         end do
        if (.not.value)
     >     call errquit('Dne_end: error deallocating stack',
     >                  0,MA_ERR)


      else
         if (mb.eq.0) then 
            size = nein(1)*nein(1) + nein(2)*nein(2)
         else
            size = nein(mb)*nein(mb)
         end if
      end if
     
      Dneall_mne_allocate_block 
     > = BA_alloc_get(mt_dbl,nb*size,'hmlab',hml(2),hml(1))
      return
      end




c
c
c     ****************************************
c     *                                      *
c     *        Dneall_mne_size                *
c     *                                      *
c     ****************************************
      integer function Dneall_mne_size(mb,nein)
      implicit none
      integer mb,nein(2)

#include "Dne.fh"
#include "errquit.fh"
#include "bafdecls.fh"

      integer size

      logical value
      integer mcqin(2),ncqin(2)
      integer mcin(2,NBLOCKS),ncin(2,NBLOCKS)
      integer ii,jj,k,ms

     
      if (mparallelized) then

         value = .true.
         mcqin(1) = 0
         mcqin(2) = 0
         ncqin(1) = 0
         ncqin(2) = 0
         do ms=1,ispin
             value = value.and.
     >           BA_push_get(mt_int,np_i,'mcin',mcin(2,ms),mcin(1,ms))
             value = value.and.
     >           BA_push_get(mt_int,np_j,'ncin',ncin(2,ms),ncin(1,ms))
             if (.not.value) 
     >         call errquit(': out of stackmemory',0,MA_ERR)
     
             call icopy(np_i,0,0,int_mb(mcin(1,ms)),1)
             call icopy(np_j,0,0,int_mb(ncin(1,ms)),1)
             ii = 0
             jj = 0
             do k=1,nein(ms)
                int_mb(mcin(1,ms)+ii) = int_mb(mcin(1,ms)+ii) + 1
                int_mb(ncin(1,ms)+jj) = int_mb(ncin(1,ms)+jj) + 1
                ii = mod(ii+1,np_i)
                jj = mod(jj+1,np_j)
             end do
             mcqin(ms) = int_mb(mcin(1,ms)+taskid_i)
             ncqin(ms) = int_mb(ncin(1,ms)+taskid_j)
         end do



         if (mb.eq.0) then
            size = mcqin(1)*ncqin(1) + mcqin(2)*ncqin(2)
         else
            size = mcqin(mb)*ncqin(mb)
         end if


         do ms=ispin,1,-1
            value = value.and.BA_pop_stack(ncin(2,ms))
            value = value.and.BA_pop_stack(mcin(2,ms))
         end do
        if (.not.value)
     >     call errquit('Dne_end: error deallocating stack',
     >                  0,MA_ERR)


      else
         if (mb.eq.0) then 
            size = nein(1)*nein(1) + nein(2)*nein(2)
         else
            size = nein(mb)*nein(mb)
         end if
      end if
     
      Dneall_mne_size = size
      return
      end



c     ****************************************
c     *                                      *
c     *        Dneall_gmg_Multiply           *
c     *                                      *
c     ****************************************

*  uses rotation algorithm

      subroutine Dneall_gmg_Multiply(mb,Ain,n2ft3d,
     >                                   hml,alpha,
     >                                   Aout,beta)
      implicit none
      integer    mb
      real*8     Ain(*)
      integer    n2ft3d
      real*8     hml(*)
      real*8     alpha
      real*8     Aout(*)
      real*8     beta

#include "bafdecls.fh"
#include "Dne.fh"

*     **** local variables ****
      integer ms,ms1,ms2,n,shift,shift2,shift3,ishift2,ishift3


      call nwpw_timing_start(16)
      if (parallelized) then
       if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = mcq(1)*ncq(1)
            ishift3 = ne(1)*ne(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
            ishift3 = 0
         end if

         if (mparallelized) then
            do ms=ms1,ms2
               shift  = 1 + (ms-1)*neq(1)*n2ft3d
               shift2 = 1 + (ms-1)*ishift2
               call DMatrix_dgemm1_rot(n2ft3d_all,ne(ms),ne(ms),
     >                alpha,
     >                Ain(shift),int_mb(ma2(1,ms)+taskid_i),
     >                           int_mb(ma2(1,ms)),
     >                           int_mb(na(1,ms)),
     >                hml(shift2),mcq(ms),
     >                            int_mb(mc(1,ms)),
     >                            int_mb(nc(1,ms)),
     >                beta,
     >                Aout(shift),int_mb(ma2(1,ms)+taskid_i),
     >                            int_mb(ma2(1,ms)),
     >                            int_mb(na(1,ms)),
     >                taskid_i,taskid_j,
     >                np_i,np_j,
     >                comm_i, comm_j,
     >                dbl_mb(bcolwork(1)),dbl_mb(bwork2(1)),
     >                dbl_mb(rwork1(1)),dbl_mb(rwork2(1)))
            end do
         else
            do ms=ms1,ms2
               shift  = 1 + (ms-1)*neq(1)*n2ft3d
               shift3 = 1 + (ms-1)*ishift3
               call DMatrix_dgemm1_rot2(n2ft3d_all,ne(ms),ne(ms),
     >                alpha,
     >                Ain(shift),int_mb(ma2(1,ms)+taskid_i),
     >                           int_mb(ma2(1,ms)),
     >                           int_mb(na(1,ms)),
     >                hml(shift3),mcq(ms),
     >                            int_mb(mc(1,ms)),
     >                            int_mb(nc(1,ms)),
     >                beta,
     >                Aout(shift),int_mb(ma2(1,ms)+taskid_i),
     >                            int_mb(ma2(1,ms)),
     >                            int_mb(na(1,ms)),
     >                taskid_i,taskid_j,
     >                np_i,np_j,
     >                comm_i, comm_j,
     >                dbl_mb(bcolwork(1)),dbl_mb(bwork2(1)),
     >                dbl_mb(rwork1(1)),dbl_mb(rwork2(1)))
            end do
         end if


      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = ne(1)*ne(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if

         do ms=ms1,ms2
            n     = ne(ms)
            if (n.le.0) go to 30
            shift  = 1 + (ms-1)*ne(1)*n2ft3d
            shift2 = 1 + (ms-1)*ishift2
            call DGEMM_OMP('N','N',n2ft3d,n,n,
     >                (alpha),
     >                Ain(shift),  n2ft3d,
     >                hml(shift2),    n,
     >                (beta),
     >                Aout(shift),n2ft3d)
   30       continue
         end do
      end if

      call nwpw_timing_end(16)
      return
      end



c     ****************************************
c     *                                      *
c     *        Dneall_ggm_sym_Multiply       *
c     *                                      *
c     ****************************************

      subroutine Dneall_ggm_sym_Multiply(mb,A1,A2,n2ft3d,hml)
      implicit none
      integer    mb
      real*8     A1(*),A2(*)
      integer    n2ft3d
      real*8     hml(*)
           
#include "bafdecls.fh"
#include "Dne.fh"

*     **** local variables ****
      integer ms,ms1,ms2,n,shift,shift2,ishift2

      call nwpw_timing_start(15)
      if (parallelized) then
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = mcq(1)*ncq(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if

         if (mparallelized) then
         do ms=ms1,ms2
            if (ne(ms).le.0) go to 20
            shift  = 1 + (ms-1)*neq(1)*n2ft3d
            shift2 = 1 + (ms-1)*ishift2
            call DMatrix_dgemm2(ne(ms),ne(ms),n2ft3d_all,128,
     >             1.0d0,
     >             A1(shift),int_mb(ma2(1,ms)+taskid_i), 
     >                       int_mb(ma2(1,ms)),
     >                       int_mb(na(1,ms)),
     >             A2(shift),int_mb(ma2(1,ms)+taskid_i), 
     >                       int_mb(ma2(1,ms)),
     >                       int_mb(na(1,ms)),
     >             0.0d0,
     >             hml(shift2),int_mb(mc(1,ms)+taskid_i), 
     >                         int_mb(mc(1,ms)),
     >                         int_mb(nc(1,ms)),
     >             taskid_i,taskid_j,
     >             np_i,np_j,
     >             comm_i, comm_j,
     >             dbl_mb(work1(1)),dbl_mb(work2(1)))

  20        continue
         end do
         else
         do ms=ms1,ms2
            if (ne(ms).le.0) go to 21
            shift  = 1 + (ms-1)*neq(1)*n2ft3d
            shift2 =     (ms-1)*ishift2
            call DMatrix_dgemm2(ne(ms),ne(ms),n2ft3d_all,128,
     >             1.0d0,
     >             A1(shift),int_mb(ma2(1,ms)+taskid_i),
     >                       int_mb(ma2(1,ms)),
     >                       int_mb(na(1,ms)),
     >             A2(shift),int_mb(ma2(1,ms)+taskid_i),
     >                       int_mb(ma2(1,ms)),
     >                       int_mb(na(1,ms)),
     >             0.0d0,
     >             dbl_mb(mat_tmp(1)+shift2),int_mb(mc(1,ms)+taskid_i),
     >                         int_mb(mc(1,ms)),
     >                         int_mb(nc(1,ms)),
     >             taskid_i,taskid_j,
     >             np_i,np_j,
     >             comm_i, comm_j,
     >             dbl_mb(work1(1)),dbl_mb(work2(1)))

  21        continue
         end do
         call Dneall_m_gather(mall(mb),mpack(mb),int_mb(mindx(1,mb)),
     >                        dbl_mb(mat_tmp(1)),hml)

         end if
      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = ne(1)*ne(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if

         do ms=ms1,ms2
            shift  = 1 + (ms-1)*ne(1)*n2ft3d
            shift2 = 1 + (ms-1)*ishift2
            n     = ne(ms)
            if (n.le.0) go to 30

            call D3dB_rrm_sym_dot(1,n,
     >                        A1(shift),
     >                        A2(shift),
     >                        hml(shift2))
  30        continue
         end do
      end if

      call nwpw_timing_end(15)
      return
      end




c     ****************************************
c     *                                      *
c     *        Dneall_ggm_Multiply           *
c     *                                      *
c     ****************************************

      subroutine Dneall_ggm_Multiply(mb,A1,A2,alpha,n2ft3d,hml,beta)
      implicit none
      integer    mb
      real*8     A1(*),A2(*),alpha
      integer    n2ft3d
      real*8     hml(*),beta
           
#include "bafdecls.fh"
#include "Dne.fh"

*     **** local variables ****
      integer ms,ms1,ms2,n,shift,shift2,ishift2,nn

      call nwpw_timing_start(15)
      if (parallelized) then
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = mcq(1)*ncq(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if

         if (mparallelized) then
         do ms=ms1,ms2
            if (ne(ms).le.0) go to 20
            shift  = 1 + (ms-1)*neq(1)*n2ft3d
            shift2 = 1 + (ms-1)*ishift2
            call DMatrix_dgemm2(ne(ms),ne(ms),n2ft3d_all,128,
     >             alpha,
     >             A1(shift),int_mb(ma2(1,ms)+taskid_i), 
     >                       int_mb(ma2(1,ms)),
     >                       int_mb(na(1,ms)),
     >             A2(shift),int_mb(ma2(1,ms)+taskid_i), 
     >                       int_mb(ma2(1,ms)),
     >                       int_mb(na(1,ms)),
     >             beta,
     >             hml(shift2),int_mb(mc(1,ms)+taskid_i), 
     >                         int_mb(mc(1,ms)),
     >                         int_mb(nc(1,ms)),
     >             taskid_i,taskid_j,
     >             np_i,np_j,
     >             comm_i, comm_j,
     >             dbl_mb(work1(1)),dbl_mb(work2(1)))

  20        continue
         end do
         else
         do ms=ms1,ms2
            if (ne(ms).le.0) go to 21
            shift  = 1 + (ms-1)*neq(1)*n2ft3d
            shift2 =     (ms-1)*ishift2
            call DMatrix_dgemm2(ne(ms),ne(ms),n2ft3d_all,128,
     >             alpha,
     >             A1(shift),int_mb(ma2(1,ms)+taskid_i),
     >                       int_mb(ma2(1,ms)),
     >                       int_mb(na(1,ms)),
     >             A2(shift),int_mb(ma2(1,ms)+taskid_i),
     >                       int_mb(ma2(1,ms)),
     >                       int_mb(na(1,ms)),
     >             beta,
     >             dbl_mb(mat_tmp(1)+shift2),int_mb(mc(1,ms)+taskid_i),
     >                         int_mb(mc(1,ms)),
     >                         int_mb(nc(1,ms)),
     >             taskid_i,taskid_j,
     >             np_i,np_j,
     >             comm_i, comm_j,
     >             dbl_mb(work1(1)),dbl_mb(work2(1)))

  21        continue
         end do
         call Dneall_m_gather(mall(mb),mpack(mb),int_mb(mindx(1,mb)),
     >                        dbl_mb(mat_tmp(1)),hml)

         end if
      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = ne(1)*ne(1)
            nn      = ne(1)*ne(1) + ne(2)*ne(2)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
            nn      = ne(mb)*ne(mb)
         end if

         do ms=ms1,ms2
            shift  = 1 + (ms-1)*ne(1)*n2ft3d
            shift2 = 1 + (ms-1)*ishift2
            n     = ne(ms)
            if (n.le.0) go to 30

            call DGEMM_OMP('T','N',n,n,n2ft3d,
     >                alpha,
     >                A1(shift), n2ft3d,
     >                A2(shift), n2ft3d,
     >                beta,
     >                hml(shift2),n)

  30        continue
         end do
         call D3dB_Vector_SumAll(nn,hml)

      end if

      call nwpw_timing_end(15)
      return
      end




c     ****************************************
c     *                                      *
c     *        Dneall_ggm_AMultiply          *
c     *                                      *
c     ****************************************

      subroutine Dneall_ggm_AMultiply(mb,A1,A2,alpha,n2ft3d,hml,beta)
      implicit none
      integer    mb
      real*8     A1(*),A2(*),alpha
      integer    n2ft3d
      real*8     hml(*),beta
           
#include "bafdecls.fh"
#include "Dne.fh"

*     **** local variables ****
      integer ms,ms1,ms2,n,shift,shift2,ishift2,nn

      call nwpw_timing_start(15)
      if (parallelized) then
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = mcq(1)*ncq(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if

         if (mparallelized) then
         do ms=ms1,ms2
            if (ne(ms).le.0) go to 20
            shift  = 1 + (ms-1)*neq(1)*n2ft3d
            shift2 = 1 + (ms-1)*ishift2
            call DMatrix_dgamm2(ne(ms),ne(ms),n2ft3d_all,128,
     >             alpha,
     >             A1(shift),int_mb(ma2(1,ms)+taskid_i), 
     >                       int_mb(ma2(1,ms)),
     >                       int_mb(na(1,ms)),
     >             A2(shift),int_mb(ma2(1,ms)+taskid_i), 
     >                       int_mb(ma2(1,ms)),
     >                       int_mb(na(1,ms)),
     >             beta,
     >             hml(shift2),int_mb(mc(1,ms)+taskid_i), 
     >                         int_mb(mc(1,ms)),
     >                         int_mb(nc(1,ms)),
     >             taskid_i,taskid_j,
     >             np_i,np_j,
     >             comm_i, comm_j,
     >             dbl_mb(work1(1)),dbl_mb(work2(1)))

  20        continue
         end do
         else
         do ms=ms1,ms2
            if (ne(ms).le.0) go to 21
            shift  = 1 + (ms-1)*neq(1)*n2ft3d
            shift2 =     (ms-1)*ishift2
            call DMatrix_dgamm2(ne(ms),ne(ms),n2ft3d_all,128,
     >             alpha,
     >             A1(shift),int_mb(ma2(1,ms)+taskid_i),
     >                       int_mb(ma2(1,ms)),
     >                       int_mb(na(1,ms)),
     >             A2(shift),int_mb(ma2(1,ms)+taskid_i),
     >                       int_mb(ma2(1,ms)),
     >                       int_mb(na(1,ms)),
     >             beta,
     >             dbl_mb(mat_tmp(1)+shift2),int_mb(mc(1,ms)+taskid_i),
     >                         int_mb(mc(1,ms)),
     >                         int_mb(nc(1,ms)),
     >             taskid_i,taskid_j,
     >             np_i,np_j,
     >             comm_i, comm_j,
     >             dbl_mb(work1(1)),dbl_mb(work2(1)))

  21        continue
         end do
         call Dneall_m_gather(mall(mb),mpack(mb),int_mb(mindx(1,mb)),
     >                        dbl_mb(mat_tmp(1)),hml)

         end if
      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = ne(1)*ne(1)
            nn      = ne(1)*ne(1) + ne(2)*ne(2)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
            nn      = ne(mb)*ne(mb)
         end if

         do ms=ms1,ms2
            shift  = 1 + (ms-1)*ne(1)*n2ft3d
            shift2 = 1 + (ms-1)*ishift2
            n     = ne(ms)
            if (n.le.0) go to 30

            call DGAMM_OMP('T','N',n,n,n2ft3d,
     >                alpha,
     >                A1(shift), n2ft3d,
     >                A2(shift), n2ft3d,
     >                beta,
     >                hml(shift2),n)

  30        continue
         end do
         call D3dB_Vector_SumAll(nn,hml)

      end if

      call nwpw_timing_end(15)
      return
      end





c     ****************************************
c     *                                      *
c     *        Dneall_www_Multiply          *
c     *                                      *
c     ****************************************

      subroutine Dneall_www_Multiply(mb,A,B,alpha,C,beta)
      implicit none
      integer mb
      complex*16 A(*),B(*),C(*)
      complex*16 alpha,beta
           
#include "bafdecls.fh"
#include "Dne.fh"

*     **** local variables ****
      integer ms,ms1,ms2,n,shift2,ishift2
  
      if (mparallelized) then
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = mcq(1)*ncq(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if

         do ms=ms1,ms2
            if (ne(ms).le.0) go to 20
            shift2 = 1 + (ms-1)*ishift2
            call CMatrix_zgemm1(ne(ms),ne(ms),ne(ms),64,
     >             alpha,
     >             A(shift2),mcq(ms),
     >                       int_mb(mc(1,ms)),
     >                       int_mb(nc(1,ms)),
     >             B(shift2),mcq(ms),
     >                       int_mb(mc(1,ms)),
     >                       int_mb(nc(1,ms)),
     >             beta,
     >             C(shift2),mcq(ms),
     >                       int_mb(mc(1,ms)),
     >                       int_mb(nc(1,ms)),
     >             taskid_i,taskid_j,
     >             np_i,np_j,
     >             comm_i, comm_j,
     >             dbl_mb(work1(1)),dbl_mb(work2(1)))

   20       continue
         end do

      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = ne(1)*ne(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if

         do ms=ms1,ms2
            n     = ne(ms)
            if (n.le.0) go to 30
            shift2 = 1 + (ms-1)*ishift2
            call ZGEMM('N','N',n,n,n,
     >                alpha,
     >                A(shift2), n,
     >                B(shift2), n,
     >                beta,
     >                C(shift2), n)
   30       continue
         end do
      end if

      return
      end

c     ****************************************
c     *                                      *
c     *           Dneall_ww_copy             *
c     *                                      *
c     ****************************************

      subroutine Dneall_ww_copy(mb,M1,M2)
      implicit none
      integer mb
      complex*16 M1(*),M2(*)

#include "Dne.fh" 

*     **** local variables ****
      integer nn

      if (mparallelized) then
         if (mb.eq.0) then
            nn = mcq(1)*ncq(1)+mcq(2)*ncq(2)
         else
            nn = mcq(mb)*ncq(mb)
         end if
      else
         if (mb.eq.0) then
            nn = ne(1)*ne(1) + ne(2)*ne(2)
         else
            nn = ne(mb)*ne(mb)
         end if 
      end if

      call dcopy(2*nn,M1,1,M2,1)
      return
      end


c     ****************************************
c     *                                      *
c     *           Dneall_mm_copy             *
c     *                                      *
c     ****************************************

      subroutine Dneall_mm_copy(mb,M1,M2)
      implicit none
      integer mb
      real*8 M1(*),M2(*)

#include "Dne.fh"

*     **** local variables ****
      integer nn

      if (mparallelized) then
         if (mb.eq.0) then
            nn = mcq(1)*ncq(1)+mcq(2)*ncq(2)
         else
            nn = mcq(mb)*ncq(mb)
         end if
      else
         if (mb.eq.0) then
            nn = ne(1)*ne(1) + ne(2)*ne(2)
         else
            nn = ne(mb)*ne(mb)
         end if
      end if

      !call dcopy(nn,M1,1,M2,1)
      call Parallel_shared_vector_copy(.true.,nn,M1,M2)
      return
      end

c     ****************************************
c     *                                      *
c     *           Dneall_mm_daxpy            *
c     *                                      *
c     ****************************************

      subroutine Dneall_mm_daxpy(mb,alpha,M1,M2)
      implicit none
      integer mb
      real*8 alpha
      real*8 M1(*),M2(*)

#include "Dne.fh"

*     **** local variables ****
      integer nn

      if (mparallelized) then
         if (mb.eq.0) then
            nn = mcq(1)*ncq(1)+mcq(2)*ncq(2)
         else
            nn = mcq(mb)*ncq(mb)
         end if
      else
         if (mb.eq.0) then
            nn = ne(1)*ne(1) + ne(2)*ne(2)
         else
            nn = ne(mb)*ne(mb)
         end if
      end if

      call daxpy_omp(nn,alpha,M1,1,M2,1)
      return
      end



c     ****************************************
c     *                                      *
c     *           Dneall_m_zero              *
c     *                                      *
c     ****************************************

      subroutine Dneall_m_zero(mb,M1)
      implicit none
      integer mb
      real*8 M1(*)

#include "Dne.fh"

*     **** local variables ****
      integer nn

      if (mparallelized) then
         if (mb.eq.0) then
            nn = mcq(1)*ncq(1)+mcq(2)*ncq(2)
         else
            nn = mcq(mb)*ncq(mb)
         end if
      else
         if (mb.eq.0) then
            nn = ne(1)*ne(1) + ne(2)*ne(2)
         else
            nn = ne(mb)*ne(mb)
         end if
      end if

      !call dcopy(nn,0.0d0,0,M1,1)
      call Parallel_shared_vector_zero(.true.,nn,M1)
      return
      end





c     ****************************************
c     *                                      *
c     *            Dneall_w_tracesqr         *
c     *                                      *
c     ****************************************

      double precision function Dneall_w_tracesqr(mb,M)
      implicit none
      integer mb
      complex*16  M(*)

#include "bafdecls.fh"
#include "Dne.fh"

*     **** local variables ****
      integer ms,ms1,ms2,shift2,ishift2,indx,i
      real*8  sum

*     ***** external functions ****
      real*8   CMatrix_tracesqr
      external CMatrix_tracesqr

      if (mparallelized) then
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = mcq(1)*ncq(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if
         sum = 0.0d0
         do ms=ms1,ms2
            shift2 = 1 + (ms-1)*ishift2
            sum =  sum + CMatrix_tracesqr(ne(ms),
     >                                 M(shift2),mcq(ms),
     >                                           int_mb(mc(1,ms)),
     >                                           int_mb(nc(1,ms)),
     >                                 taskid_i,taskid_j)
         end do

      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = ne(1)*ne(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if

         sum = 0.0d0
         do ms=ms1,ms2

           if (ne(ms).le.0) goto 30
           shift2 = (ms-1)*ishift2

           do i=1,ne(ms)
              indx = i + (i-1)*ne(ms) + shift2
              sum = sum + dble(M(indx))**2 + dimag(M(indx))**2
           end do

 30        continue
         end do
      end if

      Dneall_w_tracesqr = sum
      return
      end



c     ****************************************
c     *                                      *
c     *            Dneall_w_Max              *
c     *                                      *
c     ****************************************

      double precision function Dneall_w_max(mb,M)
      implicit none
      integer mb
      complex*16  M(*)

#include "bafdecls.fh"
#include "Dne.fh"

*     **** local variables ****
      integer ms,ms1,ms2,shift2,ishift2,indx,i,j
      real*8  sum,maxgrad

*     ***** external functions ****
      real*8   CMatrix_max
      external CMatrix_max

      maxgrad = 0.0d0
      if (mparallelized) then
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = mcq(1)*ncq(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if
         sum = 0.0d0
         do ms=ms1,ms2
            shift2 = 1 + (ms-1)*ishift2
            sum =  CMatrix_max(ne(ms),
     >                         M(shift2),mcq(ms),
     >                         int_mb(mc(1,ms)),
     >                         int_mb(nc(1,ms)),
     >                         taskid_i,taskid_j)
            if (sum.gt.maxgrad) maxgrad = sum 
         end do

      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = ne(1)*ne(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if

         do ms=ms1,ms2

           if (ne(ms).le.0) goto 30
           shift2 = (ms-1)*ishift2

           do j=1,ne(ms)
           do i=1,ne(ms)
              indx = i + (j-1)*ne(ms) + shift2
              sum =  dsqrt(dble(M(indx))**2 + dimag(M(indx))**2)
              if (sum.gt.maxgrad) maxgrad = sum 
           end do
           end do

 30        continue
         end do
      end if

      Dneall_w_max = maxgrad
      return
      end

