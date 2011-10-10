      subroutine qmmm_bq_data_init(irtdb)
      implicit none
c
#include "mafdecls.fh"
#include "rtdb.fh"
#include "bq.fh"
#include "errquit.fh"
#include "qmmm_bq_data.fh"
      integer irtdb
c
      character*30 bq_expand
      character*30 pname

      pname = "qmmm_bq_data_init"
      nbqs=0
      nbqw=0
      nbq=nbqs+nbqw

       if (.not.rtdb_get(irtdb,'qmmm:bq_exclude',mt_int,1,bq_exclude))
     + call errquit(pname,bq_exclude,
     &       RTDB_ERR)

      if (.not. rtdb_get(irtdb,'qmmm:bq_update',mt_int,1,bq_update)) 
     +    bq_update = 0
      
      readbq = .false.
      if (rtdb_cget(irtdb,"qmmm:bq:load:file",1,bqfile_in)) then
        if(bqfile_in.ne."none") then
         readbq = .true.
         if(bqfile_in.eq."default") then
           call util_file_name("bq.ind",.false.,.false.,bqfile_in)
         else
           call util_file_name_resolve(bqfile_in,.false.)
         end if
        end if
      end if


      bqsave = .false.
      if (rtdb_cget(irtdb,'qmmm:bq:save:file',1,bqfile_out)) then
        if(bqfile_out.ne."none") then
         bqsave = .true.
         if(bqfile_out.eq."default") then
           call util_file_name("bq.ind",.false.,.false.,bqfile_out)
         else
           call util_file_name_resolve(bqfile_out,.false.)
         end if
        end if
      end if
c
      if (.not.rtdb_cget(irtdb,"qmmm:bq_expand",1,bq_expand)) 
     >      bq_expand = "none"
      if(bq_expand.eq."none") then
        allbqs = .false.
        allbqw = .false.
      else if(bq_expand.eq."all") then
        allbqs = .true.
        allbqw = .true.
      else if(bq_expand.eq."solute") then
        allbqs = .true.
        allbqw = .false.
      else if(bq_expand.eq."solvent") then
        allbqs = .false.
        allbqw = .true.
      else
       call errquit(pname,"unknown option for charge expand",0,0)
      end if
c
c      if (.not.rtdb_get(irtdb,'qmmm:bq-solute',mt_log,1,allbqs))
c     >     allbqs = .false.
cc
c      if (.not.rtdb_get(irtdb,'qmmm:bq-solvent',mt_log,1,allbqw))
c     >     allbqw = .false.

      if(.not.bq_create("qmmm charges",bq_handle))
     + call errquit(pname//'Failed bq_create',0,CALC_ERR) 

      call qmmm_bq_data_load()
c
c    force separate calculation of bq energy
c    ---------------------------------------
      if (.not. rtdb_put(irtdb,'dft:bq_energy',mt_dbl,1,0.0d0))
     + call errquit(pname//'setting dft:bq_energy',0,CALC_ERR)
c

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
#include "bq.fh"
#include "util.fh"
#include "rtdb.fh"
      integer i,k
      integer i_lb
      integer i_mbq,h_mbq
      integer nlink
      integer il,j
      double precision charge
      integer mwa
      character*32 pname

      pname = "qmmm_bq_data_load"

      if(qmmm_print_debug()) 
     >       write(*,*) "in",pname

      call qmmm_bq_data_dealloc()
c
c     reset bq flags for solute and solvent
c     -------------------------------------
      if(allbqs) call mm_activate_bqszone()
      if(allbqw) call mm_activate_bqwzone()
c
c     exclude bq's as specified by qmmm input
c     ---------------------------------------
      call mm_prune_bqzone()
c
c     get numbers of solute and solvent bq's
c     --------------------------------------
      if(readbq) then
         call qmmm_read_nbq()
      else
         call mm_get_tot_nbqs(nbqs)
         call mm_get_tot_nbqw(nbqw)
      end if
      call mm_get_mwa(mwa)
      nbq=nbqs+nbqw
      mbqw = nbqw/mwa

      if(nbq.eq.0) return

c     allocate bq structure
c     ---------------------
      call qmmm_bq_data_alloc()

c     get global index for bq charges
c     -------------------------------
      if(readbq) then
        call qmmm_read_ibq()
        call mm_reset_solute_bqzone(nbqs,int_mb(i_ibq))
        call mm_reset_solvent_bqzone(nbqw,int_mb(i_ibq+nbqs))
      else
        if(nbqs.ne.0) then
          call mm_get_solute_bq_ind(nbqs,int_mb(i_ibq))
        end if
        if(nbqw.ne.0) then
          call mm_get_solvent_ind_bq(nbqw,int_mb(i_ibq+nbqs))
        end if
      end if
c
c     get solute bq coord and charges from MM
c     ---------------------------------------
      if(nbqs.ne.0) then
        call mm_get_solute_geom_bq(nbqs,
     >                       int_mb(i_ibq),
     >                       dbl_mb(i_cbq),
     >                       dbl_mb(i_qbq))

      end if

      if(qmmm_print_debug()) 
     >       write(*,*) "finished setting up solute bq charges"
c
c     append solvent bq coord and charges from MM
c     -------------------------------------------    
      if(nbqw.ne.0) then
        call mm_get_solvent_geom_bq(nbqw,int_mb(i_ibq+nbqs),
     >                       dbl_mb(i_cbq+3*nbqs),
     >                       dbl_mb(i_qbq+nbqs))
      end if
  
      if(qmmm_print_debug()) 
     >       write(*,*) "finished setting up solvent bq charges"

c
      if(qmmm_master()) then
         write(*,*) "Total number of Bq charges",nbq
         write(*,*) "Number of solute Bq charges",nbqs
         write(*,*) "Number of solvent Bq charges",nbqw
      end if
c     
c     save bq index into the file
c     ---------------------------
      if(bqsave)  
     +  call qmmm_save_bq_index()
c
c     create and activate bq structure in bq module
c     ---------------------------------------------
      if(.not.bq_pset(bq_handle,nbq,h_qbq,h_cbq))
     + call errquit(pname//'Failed bq_pset',0,CALC_ERR)
      if(.not.bq_activate(bq_handle)) 
     + call errquit(pname//'Failed bq_activate',0,CALC_ERR)
      charge = 0.0d0
      do i=1,nbq
       charge = charge + dbl_mb(i_qbq+i-1)
      end do
      if(ga_nodeid().eq.0) then
       write(*,*) "Total Bq charge: ", charge
      end if

      if(util_print('bq_geom', print_high)) 
     + call bq_print_info(bq_handle)

       if(ga_nodeid().eq.0) 
     > call qmmm_print_pdb_bq(nbq,"-bq.pdb",
     >                    dbl_mb(i_cbq),
     >                    dbl_mb(i_qbq))

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

      subroutine qmmm_bq_data_reload()
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm_bq_data.fh"
#include "qmmm.fh"
      integer counter
      character*32 pname
      save    counter
      DATA counter /0/

      pname = "qmmm_bq_data_reload"

      if(qmmm_print_debug()) 
     >       write(*,*) "in",pname

      counter = counter + 1
      if(counter.eq.bq_update) then
        call  qmmm_bq_data_load()
        counter = 0
      end if

      if(qmmm_print_debug()) 
     >       write(*,*) "out",pname

      return
      end

      subroutine qmmm_bq_data_update_active()
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm_bq_data.fh"
#include "qmmm.fh"
#include "qmmm_utils.fh"
#include "qmmm_params.fh"
#include "global.fh"
#include "bq.fh"
#include "util.fh"
      integer i,k
      character*32 pname

      pname = "qmmm_bq_data_create_active"

      if(qmmm_print_debug()) 
     >       write(*,*) "in",pname
c
      do i=1,nbq
         log_mb(i_abq+i-1)= .false.
      end do
c
c     get active list of bqs
c     ------------------------
      if(nbqs.ne.0) 
     >   call qmmm_cons_get_actmaps(nbqs,int_mb(i_ibq),
     >                       log_mb(i_abq))

      if(nbqw.ne.0) 
     >  call  qmmm_cons_get_actmapw(nbqw,int_mb(i_ibq+nbqs),
     >                       log_mb(i_abq+nbqs))

      nbqa = 0 
      do i=1,nbq
         if(log_mb(i_abq+i-1)) then
           nbqa = nbqa + 1
         end if
      end do

      if(ga_nodeid().eq.0) then
       write(*,*) "Total number of active Bq charges ",nbqa
      end if

      if(qmmm_print_debug()) 
     >       write(*,*) "out",pname
    
      end

      subroutine qmmm_bq_data_reset_active()
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm_bq_data.fh"
#include "qmmm.fh"
#include "qmmm_utils.fh"
#include "qmmm_params.fh"
#include "global.fh"
#include "bq.fh"
#include "util.fh"
      integer i,k
      character*32 pname

      pname = "qmmm_bq_data_create_active"

      if(qmmm_print_debug()) 
     >       write(*,*) "in",pname
c
      do i=1,nbq
         log_mb(i_abq+i-1)= .true.
      end do
      nbqa = nbq
      do i=1,nbq
         if(log_mb(i_abq+i-1)) then
           nbqa = nbqa + 1
         end if
      end do

      if(ga_nodeid().eq.0) then
       write(*,*) "number of active bqs ",nbqa
      end if

      if(qmmm_print_debug()) 
     >       write(*,*) "out",pname
    
      end

      subroutine qmmm_bq_print_info()
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm_bq_data.fh"
#include "qmmm.fh"
#include "qmmm_utils.fh"
#include "qmmm_params.fh"
#include "global.fh"
#include "bq.fh"
#include "util.fh"
      integer i,k
      character*32 pname

      pname = "qmmm_bq_print_info"

      if(qmmm_print_debug()) 
     >       write(*,*) "in",pname

        call util_print_centered(6,"Bq information",
     >              36, .true.)
        write(*,'(/,A,T20,A,/)') "global index","coordinates(au)"
        do i=1,nbqs+nbqw
         write(*,'(I5,T20,3F12.6)') int_mb(i_ibq+i-1),
     >              (dbl_mb(i_cbq+3*(i-1)+k-1),k=1,3)
        end do

      if(qmmm_print_debug()) 
     >       write(*,*) "out",pname
    
      end

      subroutine qmmm_bq_coord_update()
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm_bq_data.fh"
#include "qmmm.fh"
#include "global.fh"
c
      integer i
      if(nbqs.ne.0) then
        call mm_get_solute_coord(nbqs,
     >                       int_mb(i_ibq),
     >                       dbl_mb(i_cbq))

      end if

c
c     append solvent bq coord and charges from MM
c     -------------------------------------------    
      if(nbqw.ne.0) then
        call mm_get_solvent_coord(nbqw,int_mb(i_ibq+nbqs),
     >                       dbl_mb(i_cbq+3*nbqs))
      end if

      if(qmmm_print_debug()) 
     >       write(*,*) "finshed updating up bq coords"
    
      end

      subroutine qmmm_bq_data_alloc()
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm_bq_data.fh"
      integer i
      if(nbq.eq.0) return

      if(.not.ma_alloc_get(mt_log,nbq,'abq',h_abq,i_abq))
     + call errquit('qmmm: Failed to allocate memory for abq',nbq,
     +       MA_ERR)
      do i=1,nbq
        log_mb(i_abq+i-1)=.true.
      end do

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

      if(.not.ma_alloc_get(mt_int,nbq,'i_ibqa',h_ibqa,i_ibqa))
     + call errquit('qmmm: Failed to allocate memory for i_ibqa',nbq,
     &       MA_ERR)
      call ifill(nbq,0,int_mb(i_ibqa),1)
      
      end

      subroutine qmmm_bq_data_dealloc()
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm_bq_data.fh"
#include "bq.fh"
      logical ignore

      if(nbq.eq.0) return

      if(.not.ma_free_heap(h_ibqa))
     + call errquit('qmmm:memory deallocation for ibqa',nbqa,
     &       MA_ERR)
  
      if(.not.ma_free_heap(h_ibq))
     + call errquit('qmmm:memory deallocation for ibq',nbq,
     &       MA_ERR)
      
      if(.not.ma_free_heap(h_qbq))
     + call errquit('qmmm:memory deallocation for qbq',nbq,
     &       MA_ERR)

      if(.not.ma_free_heap(h_cbq))
     + call errquit('qmmm: memory deallocation for cbq',3*nbq,
     +       MA_ERR)

      if(.not.ma_free_heap(h_abq))
     + call errquit('qmmm: memory deallocation for cbq',3*nbq,
     +       MA_ERR)

       nbqs=0
       nbqw=0
       nbq=0

       ignore = bq_destroy(bq_handle)


      end

      function qmmm_get_nbq()
      implicit none
#include "qmmm_bq_data.fh"

      integer qmmm_get_nbq

      qmmm_get_nbq = nbqs+nbqw

      end 

      function qmmm_get_nbqa()
      implicit none
#include "qmmm_bq_data.fh"

      integer qmmm_get_nbqa

      qmmm_get_nbqa = nbqa

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

      function qmmm_get_i_ibqa()
      implicit none
#include "qmmm_bq_data.fh"
      integer qmmm_get_i_ibqa

      qmmm_get_i_ibqa = i_ibqa

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

      subroutine qmmm_read_ibq()
      implicit none
#include "stdio.fh"
#include "errquit.fh"
#include "qmmm.fh"
#include "util.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
#include "qmmm_bq_data.fh"

c     local variables
      integer nf
      character*30 pname

      integer i
      integer fn
      integer anbq,anbqs,anbqw

      pname = "qmmm_read_ibq"

      if(qmmm_master()) then
        if(.not.qmmm_get_io_unit(fn)) 
     >       call errquit("cannot get file number",0,0)
        open(unit=fn,status="old",form="formatted",file=bqfile_in)
        read(fn,*) anbqs,anbqw,anbq
        do i=1,anbq
          read(fn,*) int_mb(i_ibq+i-1)
        end do
        close(fn)
        call util_print_centered(luout,
     >              "Bq index has been read from external file "//
     >               bqfile_in,
     >              36, .true.)
      end if
      call ga_brdcst(msg_qmmm_ind,int_mb(i_ibq),
     >      nbq*ma_sizeof(mt_int,1,mt_byte),0)

      call ga_sync()

      end

      subroutine qmmm_read_nbq()
      implicit none
#include "errquit.fh"
#include "qmmm.fh"
#include "util.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
#include "qmmm_bq_data.fh"

c     local variables
      integer nf

      integer i
      integer fn
 
      if(qmmm_master()) then
        if(.not.qmmm_get_io_unit(fn)) 
     >       call errquit("cannot get file number",0,0)
        open(unit=fn,status="old",form="formatted",file=bqfile_in)
        read(fn,*) nbqs,nbqw
        close(fn)
      end if

      call ga_brdcst(msg_qmmm_nbqs,nbqs,ma_sizeof(mt_int,1,mt_byte),0)
      call ga_brdcst(msg_qmmm_nbqw,nbqw,ma_sizeof(mt_int,1,mt_byte),0)


      end

      subroutine qmmm_save_bq_index()
      implicit none
#include "stdio.fh"
#include "errquit.fh"
#include "qmmm.fh"
#include "util.fh"
#include "global.fh"
#include "inp.fh"
#include "qmmm_bq_data.fh"
#include "mafdecls.fh"

c     local variables
      integer nf
      character*(nw_max_path_len) filename

      integer i
      integer fn
 
      if(qmmm_master()) then

      if(.not.qmmm_get_io_unit(fn))
     >       call errquit("cannot get file number",0,0)

      open(unit=fn,status="unknown",form="formatted",file=bqfile_out)

      write(fn,*) nbqs,nbqw,nbq
      do i=1,nbq
       write(fn,*) int_mb(i_ibq+i-1)
      end do

      call util_flush(fn)
      close(fn)

      call util_print_centered(luout,
     >              "Bq index has been saved to external file "//
     >               bqfile_out,
     >              36, .true.)

      end if

      call ga_sync()

      end

      subroutine qmmm_bq_push_grad(nt,
     >                       g)
      implicit none

#include "qmmm_bq_data.fh"
#include "qmmm_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"

      integer nt
      double precision g(3*nt)
c
      integer i_ibqg,h_ibqg
      integer i,j,k,ia
      integer ii,jj
      integer nbqsa,nbqwa,mwa
      character*32 pname

      pname = "qmmm_bq_push_grad"
      if(.not.ma_push_get(mt_int,nbqa,'ibqg',h_ibqg,i_ibqg))
     +   call errquit( pname//'Failed to allocate memory for ibqs',
     +   nbqa, MA_ERR)


      nbqsa = 0
      do i=1,nbqs
         if(log_mb(i_abq+i-1)) then
           nbqsa = nbqsa + 1
           int_mb(i_ibqg+nbqsa-1) = int_mb(i_ibq+i-1)
         end if
      end do

      call mm_get_mwa(mwa)
      nbqwa = 0
      do i=1,nbq-nbqs,mwa
         if(log_mb(i_abq+nbqs+i-1)) then
           nbqwa = nbqwa + 1
           int_mb(i_ibqg+nbqsa+nbqwa-1) = int_mb(i_ibq+nbqs+(i-1)/mwa)
         end if
      end do


      if(nbqsa.ne.0)   
     >  call mm_add_solute_force(nbqsa,
     >                            int_mb(i_ibqg),
     >                            g)


      if(nbqwa.ne.0)   
     > call mm_add_solvent_force(mwa*nbqwa,
     >                            int_mb(i_ibqg+nbqsa),
     >                            g(3*nbqsa+1))

     
      if(.not.ma_pop_stack(h_ibqg))
     +   call errquit( pname//'Failed to deallocate memory for ibqg',
     +   nbqa, MA_ERR)


      end

      subroutine qmmm_get_cqbqa(nt,
     >                       c,
     >                       q)
      implicit none

#include "qmmm_bq_data.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm_params.fh"

      integer nt
      double precision c(3,nt)
      double precision q(nt)
c
      integer i,j,k,ia
      integer ii,jj


      ia = 0
      do i=1,nbqs 
        if(log_mb(i_abq+i-1)) then
           ia = ia + 1
           q(ia) = dbl_mb(i_qbq+i-1)
           do k=1,3
             c(k,ia) = dbl_mb(i_cbq+(i-1)*3+k-1)
           end do
        end if
      end do

      do i=nbqs+1,nbq
        if(log_mb(i_abq+i-1)) then
           ia = ia + 1
           q(ia) = dbl_mb(i_qbq+i-1)
           do k=1,3
             c(k,ia) = dbl_mb(i_cbq+(i-1)*3+k-1)
           end do
        end if
      end do


      end

c $Id$
