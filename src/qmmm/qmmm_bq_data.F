      subroutine qmmm_bq_data_init(irtdb)
      implicit none
c
#include "mafdecls.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "qmmm_bq_data.fh"
      integer irtdb
c
      nbqs=0
      nbqw=0
      nbq=nbqs+nbqw

       if (.not.rtdb_get(irtdb,'qmmm:bq_exclude',mt_int,1,bq_exclude))
     + call errquit('qmmm_data_init: qmmm:bq_exclude',bq_exclude,
     &       RTDB_ERR)


      end

      subroutine qmmm_bq_data_load()
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm_bq_data.fh"
#include "qmmm.fh"
#include "qmmm_utils.fh"
#include "qmmm_params.fh"
#include "global.fh"
      logical link
      integer i,k
      integer i_lb
      integer i_mbq,h_mbq
      integer nlink
      integer il
      character*32 pname

      pname = "qmmm_bq_data_load"

      if(qmmm_print_debug()) 
     >       write(*,*) "in",pname
      link = .true.

      call qmmm_bq_data_dealloc()
c
c      if(bq_exclude.eq.qmmm_all_atoms) return

      call mm_prune_bqzone()

      call mm_get_tot_nbqs(nbqs,link)
      call mm_get_tot_nbqw(nbqw)
      nbq=nbqs+nbqw
      
      if(qmmm_print_debug()) then
         write(*,*) "Total number of Bq charges",nbq
         write(*,*) "Number of solute Bq charges",nbqs
         write(*,*) "Number of solvent Bq charges",nbqw
      end if

      if(nbq.eq.0) return

      call qmmm_bq_data_alloc()

      if(nbqs.ne.0) then
        call mm_get_solute_bq_ind(nbqs,link,int_mb(i_ibq))
        call mm_get_solute_geom_bq(nbqs,link,
     >                       int_mb(i_ibq),
     >                       dbl_mb(i_cbq),
     >                       dbl_mb(i_qbq))

        if(bq_exclude.eq.qmmm_all_atoms) then
          do i=1,nbqs
             dbl_mb(i_qbq+i-1) = 0.0d0
          end do
        end if
      end if
c
c     number of link atoms
c     -------------------
      nlink = qmmm_get_nlink()

      if(nlink.gt.0 .and.
     >   bq_exclude.eq.qmmm_linkbond) then
c
c       mapping from global to local index for link atoms
c       -------------------------------------------------
        if(.not.ma_push_get(mt_int,nlink,'mbq',h_mbq,i_mbq))
     +   call errquit( pname//'Failed to allocate memory for mbq',
     +   nlink, MA_ERR)
cc
c        get pointer to index of link pairs
c        ----------------------------------
         i_lb = qmmm_get_i_lb()
c
c        map to a local index for bq atoms
c        ---------------------------------
         if(.not.qmmm_map(nbqs,int_mb(i_ibq),
     >                    nlink,int_mb(i_lb+nlink),
     >                    int_mb(i_mbq)))  
     >        call errquit( pname//'Failed to map bq atoms',0,0)

c
         do i=1,nlink
           il=int_mb(i_mbq+i-1)
           dbl_mb(i_qbq+il-1) = 0.0d0     
         end do
        if(.not.ma_pop_stack(h_mbq))
     +   call errquit( pname//'Failed to deallocate memory for mbq',
     +   nlink, MA_ERR)
      end if

      if(qmmm_print_debug()) 
     >       write(*,*) "finished setting up solute bq charges"
c
c     append solvent bq coord and charges from MM
c     -------------------------------------------    
      if(nbqw.ne.0) then
        call mm_get_solvent_ind_bq(nbqw,int_mb(i_ibq+nbqs))
        call mm_get_solvent_geom(nbqw,int_mb(i_ibq+nbqs),
     >                       dbl_mb(i_cbq+3*nbqs),
     >                       dbl_mb(i_qbq+nbqs))
      end if

      if(qmmm_print_debug()) 
     >       write(*,*) "finished setting up solvent bq charges"

      if(qmmm_print_debug()) then
        call util_print_centered(6,"Bq information",
     >              36, .true.)
        write(*,'(/,A,T20,A,/)') "global index","coordinates(au)"
        do i=1,nbqs+nbqw
         write(*,'(I5,T20,3F12.6)') int_mb(i_ibq+i-1),
     >              (dbl_mb(i_cbq+3*(i-1)+k-1),k=1,3)
        end do
      end if

      if(qmmm_print_debug()) 
     >       write(*,*) "out",pname
    
      end

      subroutine qmmm_bq_update()
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm_bq_data.fh"
#include "qmmm.fh"
#include "global.fh"
      logical link
      integer i,k
c
      link = .true.
      if(qmmm_print_debug()) then
         write(*,*) "Total number of Bq charges",nbq
         write(*,*) "Number of solute Bq charges",nbqs
         write(*,*) "Number of solvent Bq charges",nbqw
      end if

      if(nbq.eq.0) return

      if(nbqs.ne.0) then
        call mm_get_solute_geom_bq(nbqs,link,
     >                       int_mb(i_ibq),
     >                       dbl_mb(i_cbq),
     >                       dbl_mb(i_qbq))

      end if

c
c     append solvent bq coord and charges from MM
c     -------------------------------------------    
      if(nbqw.ne.0) then
        call mm_get_solvent_ind_bq(nbqw,int_mb(i_ibq+nbqs))
        call mm_get_solvent_geom(nbqw,int_mb(i_ibq+nbqs),
     >                       dbl_mb(i_cbq+3*nbqs),
     >                       dbl_mb(i_qbq+nbqs))
      end if

      if(qmmm_print_debug()) then
        call util_print_centered(6,"Bq information",
     >              36, .true.)
        write(*,'(/,A,T20,A,/)') "global index","coordinates(au)"
        do i=1,nbqs+nbqw
         write(*,'(I5,T20,3F12.6)') int_mb(i_ibq+i-1),
     >              (dbl_mb(i_cbq+3*(i-1)+k-1),k=1,3)
        end do
      end if

      if(qmmm_print_debug()) 
     >       write(*,*) "finshed updating up bq charges"
    
      end

      subroutine qmmm_bq_data_alloc()
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm_bq_data.fh"

      if(nbq.eq.0) return

      if(.not.ma_alloc_get(mt_dbl,3*nbq,'cbq',h_cbq,i_cbq))
     + call errquit('qmmm: Failed to allocate memory for cbq',3*nbq,
     +       MA_ERR)
      call dfill(3*nbq,0.d0,dbl_mb(i_cbq),1)
      if(.not.ma_alloc_get(mt_dbl,nbq,'qbq',h_qbq,i_qbq))
     + call errquit('qmmm: Failed to allocate memory for qbq',nbq,
     &       MA_ERR)
      call dfill(nbq,0.d0,dbl_mb(i_qbq),1)
      if(.not.ma_alloc_get(mt_int,nbq,'ibq',h_ibq,i_ibq))
     + call errquit('qmmm: Failed to allocate memory for ibq',nbq,
     &       MA_ERR)
      call ifill(nbq,0,int_mb(i_ibq),1)
      
      end

      subroutine qmmm_bq_data_dealloc()
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm_bq_data.fh"

      if(nbq.eq.0) return

      if(.not.ma_free_heap(h_ibq))
     + call errquit('qmmm:memory deallocation for ibq',nbq,
     &       MA_ERR)
      
      if(.not.ma_free_heap(h_qbq))
     + call errquit('qmmm:memory deallocation for qbq',nbq,
     &       MA_ERR)

      if(.not.ma_free_heap(h_cbq))
     + call errquit('qmmm: memory deallocation for cbq',3*nbq,
     +       MA_ERR)

       nbqs=0
       nbqw=0
       nbq=0

      end

      function qmmm_get_nbq()
      implicit none
#include "qmmm_bq_data.fh"

      integer qmmm_get_nbq

      qmmm_get_nbq = nbqs+nbqw

      end 

      function qmmm_get_nbqs()
      implicit none
#include "qmmm_bq_data.fh"

      integer qmmm_get_nbqs

      qmmm_get_nbqs = nbqs

      end 

      function qmmm_get_nbqw()
      implicit none
#include "qmmm_bq_data.fh"

      integer qmmm_get_nbqw

      qmmm_get_nbqw = nbqw

      end 

      function qmmm_get_h_qbq()
      implicit none
#include "qmmm_bq_data.fh"
      integer qmmm_get_h_qbq

      qmmm_get_h_qbq = h_qbq

      end 

      function qmmm_get_i_qbq()
      implicit none
#include "qmmm_bq_data.fh"
      integer qmmm_get_i_qbq

      qmmm_get_i_qbq = i_qbq

      end 

      function qmmm_get_h_cbq()
      implicit none
#include "qmmm_bq_data.fh"
      integer qmmm_get_h_cbq

      qmmm_get_h_cbq = h_cbq

      end 

      function qmmm_get_i_cbq()
      implicit none
#include "qmmm_bq_data.fh"
      integer qmmm_get_i_cbq

      qmmm_get_i_cbq = i_cbq

      end 

      function qmmm_get_i_ibq()
      implicit none
#include "qmmm_bq_data.fh"
      integer qmmm_get_i_ibq

      qmmm_get_i_ibq = i_ibq

      end 

      function qmmm_get_i_tbq()
      implicit none
#include "qmmm_bq_data.fh"
      integer qmmm_get_i_tbq

      qmmm_get_i_tbq = i_tbq

      end 

      function qmmm_get_h_cbqs()
      implicit none
#include "qmmm_bq_data.fh"
      integer qmmm_get_h_cbqs

      qmmm_get_h_cbqs = h_cbqs

      end 

      function qmmm_get_i_cbqs()
      implicit none
#include "qmmm_bq_data.fh"
      integer qmmm_get_i_cbqs

      qmmm_get_i_cbqs = i_cbqs

      end 

      function qmmm_get_h_cbqw()
      implicit none
#include "qmmm_bq_data.fh"
      integer qmmm_get_h_cbqw

      qmmm_get_h_cbqw = h_cbqw

      end 

      function qmmm_get_i_cbqw()
      implicit none
#include "qmmm_bq_data.fh"
      integer qmmm_get_i_cbqw

      qmmm_get_i_cbqw = i_cbqw

      end 

      function qmmm_get_i_qbqw()
      implicit none
#include "qmmm_bq_data.fh"
      integer qmmm_get_i_qbqw

      qmmm_get_i_qbqw = i_qbqw

      end 

      function qmmm_get_bq_exclude()
      implicit none
#include "qmmm_bq_data.fh"

      integer qmmm_get_bq_exclude

      qmmm_get_bq_exclude = bq_exclude

      end 


