c $I#d: cons.F,v 1.1 2004/01/28 01:30:59 marat Exp $
      function cons_active(rtdb)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "geom.fh"
#include "cons_params.fh"
#include "cons.fh"
#include "global.fh"
      logical cons_active
      integer rtdb
c
      character*255 namespace
c     set default namespace if none specified
c     ---------------------------------------
      if(.not. rtdb_cget(rtdb,"constraints" , 1, namespace)) then
         namespace=cons_default
      end if
c     check if the namespace existst at all
      cons_active = cons_inquire(namespace,rtdb)
c     check if constraints are enabled
      cons_active = cons_active .and. cons_status(namespace,rtdb)
      return
      end
c
      subroutine cons_inquire_fixed_atoms(status,rtdb)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
      logical status
      integer rtdb
c
      character*26 ma_date
      integer ma_type
      integer nact
c
c     check for existence of active list
c     ----------------------------------
      status=rtdb_get_info(rtdb, 'geometry:actlist', ma_type,
     $        nact, ma_date)
c
c     return false if list of zero length
c     -----------------------------------
      if(status) then
        if(nact.eq.0) status=.false.
      end if
      return
      end

      subroutine cons_save_fixed_atoms(namespace,rtdb)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "geom.fh"
#include "cons_params.fh"
#include "cons.fh"
#include "global.fh"
      character*(*) namespace
      integer rtdb
c
      logical ignore
      logical status
      integer ma_type
      integer i
      integer nat
      integer nact
      integer h_act,i_act
      integer h_tmp_id,i_tmp_id 
      logical geom_rtdb_ncent
      external geom_rtdb_ncent

c
c     check if fixed atoms exist
c     --------------------------
      if (rtdb_ma_get(rtdb, 'geometry:actlist', ma_type,
     $        nact, h_act)) then
            if (.not. ma_get_index(h_act, i_act))
     $           call errquit('cons_init: ma_get_index failed',h_act,
     &       MA_ERR)
      else
c       nothing to process
        return
      end if
c
c     get number of atoms
c     --------------------
      status=geom_rtdb_ncent(rtdb,"geometry",nat)
      if(.not.status) 
     & call errquit('cons_init: geom_create?',70, GEOM_ERR)
c
      if ( .not. ma_push_get( MT_INT, nat, 'h_tmp_id',
     &      h_tmp_id, i_tmp_id) ) call errquit(
     &      'cons_init: unable to allocate scratch space',
     &      nat, MA_ERR)

c
c     delete prior fixed atoms if any
c     -------------------------------
      call cons_delete_fixed_atoms(namespace,rtdb)

      do i=1,nat
         int_mb(i_tmp_id+i-1)=i
      end do
      call cons_add_fixed_atoms(namespace,nat,int_mb(i_tmp_id),rtdb)
      call cons_unfix_atoms(namespace,nact,int_mb(i_act),rtdb)

c     release memory
        if (.not. ma_pop_stack(h_tmp_id) ) call errquit(
     &    'cons_init: unable to free h_tmp_id',
     &    0, MA_ERR)
        if (.not.ma_free_heap(h_act)) call errquit(
     &    'cons_init: unable to free actlist',
     &    0, MA_ERR)

      return
      end

      subroutine cons_process_old_style(namespace,rtdb)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "geom.fh"
#include "cons_params.fh"
#include "cons.fh"
#include "global.fh"
      character*(*) namespace
      integer rtdb
c
      logical ignore
      logical new_style
      logical status
      integer ma_type
      integer i
      integer nat
      integer nact
      integer h_act,i_act
      integer h_tmp_id,i_tmp_id 
      logical geom_rtdb_ncent
      external geom_rtdb_ncent

c
c     check if fixed atoms were specified using old style
c     ---------------------------------------------------
      if (rtdb_ma_get(rtdb, 'geometry:actlist', ma_type,
     $        nact, h_act)) then
            if (.not. ma_get_index(h_act, i_act))
     $           call errquit('cons_init: ma_get_index failed',h_act,
     &       MA_ERR)
      else
c       nothing to process
        return
      end if
c
c     get number of atoms
c     --------------------
      status=geom_rtdb_ncent(rtdb,"geometry",nat)
      if(.not.status) 
     & call errquit('cons_init: geom_create?',70, GEOM_ERR)
c
      if ( .not. ma_push_get( MT_INT, nat, 'h_tmp_id',
     &      h_tmp_id, i_tmp_id) ) call errquit(
     &      'cons_init: unable to allocate scratch space',
     &      nat, MA_ERR)

c
c     check of fixed atoms were specified using new style
c     ---------------------------------------------------
      if(cons_fixed_atoms(namespace,rtdb)) then
        new_style = .true.
      else
        new_style = .false. 
      end if

      if(new_style) then
        call cons_unfix_atoms(namespace,nact,int_mb(i_act),rtdb)
      else
        do i=1,nat
          int_mb(i_tmp_id+i-1)=i
        end do
        call cons_add_fixed_atoms(namespace,nat,int_mb(i_tmp_id),rtdb)
        call cons_unfix_atoms(namespace,nact,int_mb(i_act),rtdb)
      end if

c     release memory
        if (.not. ma_pop_stack(h_tmp_id) ) call errquit(
     &    'cons_init: unable to free h_tmp_id',
     &    0, MA_ERR)
        if (.not.ma_free_heap(h_act)) call errquit(
     &    'cons_init: unable to free actlist',
     &    0, MA_ERR)

      return
      end

      subroutine cons_proces_old_style(rtdb)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "geom.fh"
#include "cons_params.fh"
#include "cons.fh"
#include "global.fh"
      integer rtdb
c
      logical old_style
      logical new_style
      logical status
      integer geom
      integer ma_type
      integer i
      integer nat
      integer nact
      integer h_act,i_act
      integer h_tmp_id,i_tmp_id 
      character*255 namespace
      logical oprint
      logical ocons
      write(*,*) "in cons init"
c
c     set default namespace if none specified
c     ---------------------------------------
      if(.not. rtdb_cget(rtdb,"constraints" , 1, namespace)) then
         namespace=cons_default
      end if
c
      call cons_create(namespace,rtdb)
c
      call util_print_push()
      call util_print_rtdb_load(rtdb,'cons')
      oprint = util_print('information', print_debug)
      oprint = oprint .and. cons_print()

      if(oprint) 
     >    write(*,*) "entered cons_init"
c
      ocons = .false.
c
c     initialize harmonic bobds
      if(cons_hbonds(namespace,rtdb)) then
      ocons = .true.
      call cons_set_hbond_nhb(0)
      call cons_load_hbonds(namespace,rtdb)
      end if
c
c     load geometry 
c     -------------
      if (.not. geom_create(geom, 'geometry'))
     &     call errquit('cons_init: geom_create?',70, GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     &     call errquit('cons_init: no geometry ',72, RTDB_ERR)
c
c     get number of atoms
c     --------------------
      status=geom_ncent(geom,nat)
      if(.not.status) 
     & call errquit('cons_init: geom_create?',70, GEOM_ERR)

c
c     check if fixed atoms were specified using old style
c     ---------------------------------------------------
      if (rtdb_ma_get(rtdb, 'geometry:actlist', ma_type,
     $        nact, h_act)) then
            if (.not. ma_get_index(h_act, i_act))
     $           call errquit('cons_init: ma_get_index failed',h_act,
     &       MA_ERR)
       if(cons_print()) then
         call banner(6,"Warning!!!                               ","*" ,
     >             .true.,.false.,.true. )
         call banner(6,"geometry:actlist is an obsolete option   ","*" ,
     >             .false.,.false.,.true. )
         call banner(6,"consider using fixed atoms constraints   ","*" ,
     >             .false.,.true.,.true. )
       end if
       old_style=.true.
      else
       old_style=.false.
      end if

c
c     check of fixed atoms were specified using new style
c     ---------------------------------------------------
      if(cons_fixed_atoms(namespace,rtdb)) then
      new_style = .true.
      else
      new_style = .false. 
      end if

      if(old_style) then
       if ( .not. ma_push_get( MT_INT, nat, 'h_tmp_id',
     &      h_tmp_id, i_tmp_id) ) call errquit(
     &      'cons_init: unable to allocate scratch space',
     &      nat, MA_ERR)
        if(new_style) then
         if(cons_print()) then
         call banner(6,"Warning!!!                            ","*" ,
     >           .true.,.false.,.true. )
         call banner(6,"mixing both geometry:actlist  and     ","*" ,
     >             .false.,.false.,.true. )
         call banner(6,"fixed atom constraints is not         ","*" ,
     >             .false.,.false.,.true. )
         call banner(6,"recommended. Check output below for   ","*" ,
     >             .false.,.false.,.true. )
         call banner(6,"for potential problems.               ","*" ,
     >             .false.,.true.,.true. )
          end if
          call cons_unfix_atoms(namespace,nact,int_mb(i_act),rtdb)
        else
          do i=1,nat
            int_mb(i_tmp_id+i-1)=i
          end do
          call cons_add_fixed_atoms(namespace,nat,int_mb(i_tmp_id),rtdb)
          call cons_unfix_atoms(namespace,nact,int_mb(i_act),rtdb)
        end if
        if (.not. ma_pop_stack(h_tmp_id) ) call errquit(
     &    'cons_init: unable to free h_tmp_id',
     &    0, MA_ERR)
        if (.not.ma_free_heap(h_act)) call errquit(
     &    'cons_init: unable to free actlist',
     &    0, MA_ERR)

      end if

      status=geom_destroy(geom)
      if(.not.status) 
     & call errquit('cons_init: geom_destroy?',70, GEOM_ERR)

      call cons_load_fixed_atoms(namespace,rtdb)

      call cons_print_gen_info(namespace,rtdb)
c
      if(oprint) 
     >    write(*,*) "finished cons_init"
c
      call util_print_pop()
      write(*,*) "out cons init"

      return
      end
c
      subroutine cons_init(rtdb)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "geom.fh"
#include "cons_params.fh"
#include "cons.fh"
#include "global.fh"
      integer rtdb
c
      character*255 namespace
      logical oprint
      logical ocons
c
      call cons_get_namespace(namespace,rtdb)
      if(.not.cons_inquire(namespace,rtdb)) then
          call errquit('cons_init',0, RTDB_ERR)
      end if
c
      call util_print_push()
      call util_print_rtdb_load(rtdb,'cons')
      oprint = util_print('information', print_debug)
      oprint = oprint .and. cons_print()

      if(oprint) 
     >    write(*,*) "entered cons_init"
c
      ocons = .false.
c
c     initialize harmonic bobds
      if(cons_hbonds(namespace,rtdb)) then
      ocons = .true.
      call cons_set_hbond_nhb(0)
      call cons_load_hbonds(namespace,rtdb)
      end if
c
c
c      call cons_process_old_style(namespace,rtdb)
      call cons_load_fixed_atoms(namespace,rtdb)

      call cons_print_gen_info(namespace,rtdb)
c
      if(oprint) 
     >    write(*,*) "finished cons_init"
c
      call util_print_pop()

      return
      end
c
      subroutine cons_end(rtdb)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "geom.fh"
#include "cons_params.fh"
#include "cons.fh"
#include "global.fh"
      integer rtdb
c
      call cons_print_summary(rtdb)
      call cons_unload_hbonds()
      call cons_unload_fixed_atoms(rtdb)
      return
      end
c
      subroutine cons_load_geom(rtdb)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "geom.fh"
#include "cons.fh"
      integer rtdb
c
      integer i_c,h_c
      integer geom
      integer nat
      logical status
c
c     load geometry 
c     -------------
      if (.not. geom_create(geom, 'geometry'))
     &     call errquit('cons_load_geom',0, GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     &     call errquit('cons_load_geom',0, RTDB_ERR)
c
c     get cart coordinates
c     --------------------
      status=geom_ncent(geom,nat)
      if(.not.status) 
     & call errquit('cons_init: geom_create?',70, GEOM_ERR)
      if ( .not. ma_alloc_get( MT_DBL, nat*3, 'cons coord',
     &      h_c, i_c) ) call errquit(
     &      'cons_load_geom: unable to allocate cons coord',
     &      nat*3, MA_ERR)
      status=geom_cart_coords_get(geom, dbl_mb(i_c))
      if(.not.status) 
     & call errquit('cons_init: geom_create?',70, GEOM_ERR)
 
      status=geom_destroy(geom)
      if(.not.status) 
     & call errquit('cons_init: geom_destroy?',70, GEOM_ERR)
      call cons_set_h_c(h_c)
c      call cons_set_i_c(i_c)

      return
      end
c
      subroutine cons_unload_geom()
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "cons.fh"

      integer h_c
      h_c = cons_get_h_c()
      if ( .not. ma_free_heap(h_c)) call errquit(
     &      'input_set: unable to deallocate cons coord',
     &      0, MA_ERR)

      return
      end

      subroutine cons_unload_fixed_atoms(rtdb)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "cons.fh"
#include "geom.fh"
c     
      integer rtdb
c
      logical ignore

      ignore= rtdb_delete(rtdb,"geometry:actlist")

      end
c
      subroutine cons_load_fixed_atoms(namespace,rtdb)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "cons.fh"
#include "geom.fh"
c     
      character*(*) namespace
      integer rtdb
c
      integer nfa
      integer i,j,id
      logical status
      character*255 tag_id
      character*255 tag_n
      integer h_tmp_id,i_tmp_id 
      integer h_tmp_id1,i_tmp_id1 
      integer h_active,i_active
      integer nat
      integer geom
c
      call cons_fxd_atoms_id_tag(namespace,tag_id)
      call cons_fxd_atoms_n_tag(namespace,tag_n)
c
c     load fixed atoms
c
      status = rtdb_get(rtdb,tag_n,
     >     mt_int,1,nfa) 
      if(.not.status) then
          nfa=0
      end if
      if(nfa.eq.0) then
        return
      end if
c
c     load geometry 
c     -------------
      if (.not. geom_create(geom, 'geometry'))
     &     call errquit('cons_init: geom_create?',70, GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     &     call errquit('cons_init: no geometry ',72, RTDB_ERR)
c
c     get number of atoms
c     --------------------
      status=geom_ncent(geom,nat)
      if(.not.status) 
     & call errquit('cons_init: geom_create?',70, GEOM_ERR)

      status=geom_destroy(geom)
      if(.not.status) 
     & call errquit('cons_init: geom_destroy?',70, GEOM_ERR)

      if ( .not. ma_push_get( MT_INT, nfa, 'h_tmp_id',
     &      h_tmp_id, i_tmp_id) ) call errquit(
     &      'cons_add_fixed_atoms: unable to allocate scratch space',
     &      nfa, MA_ERR)

      if ( .not. ma_push_get( MT_INT, nat-nfa, 'h_tmp_id1',
     &      h_tmp_id1, i_tmp_id1) ) call errquit(
     &      'cons_add_fixed_atoms: unable to allocate scratch space',
     &      nfa, MA_ERR)

      if ( .not. ma_push_get( MT_LOG, nat, 'h_active',
     &      h_active, i_active) ) call errquit(
     &      'cons_add_fixed_atoms: unable to allocate scratch space',
     &      nfa, MA_ERR)

        status = rtdb_get(rtdb,tag_id,
     >     mt_int,nfa,int_mb(i_tmp_id)) 
        if(.not.status) call errquit(
     >     'cons_load_fixed_atoms: unable to get prior fixed atoms',
     >      0, MA_ERR)
        
        do i=1,nat
          log_mb(i_active+i-1)=.true.
        end do

        do i=1,nfa
          id = int_mb(i_tmp_id+i-1)
          log_mb(i_active+id-1)=.false.
        end do


        j=0
        do i=1,nat
          if(log_mb(i_active+i-1)) then
            j=j+1
            int_mb(i_tmp_id1+j-1)=i  
          end if
        end do

      status = rtdb_put(rtdb,'geometry:actlist',
     >     mt_int,nat-nfa,int_mb(i_tmp_id1))
        if(.not.status) call errquit(
     >     'cons_load_fixed_atoms: unable to store fixed atoms',
     >      0, MA_ERR)

c      status = rtdb_put(rtdb,'cons:inactlist',
c     >     mt_int,nfa,int_mb(i_tmp_id))
c        if(.not.status) call errquit(
c     >     'cons_load_fixed_atoms: unable to store fixed atoms',
c     >      0, MA_ERR)
c

      if (.not. ma_pop_stack(h_active) ) call errquit(
     &    'cons_load_fixed_atoms:ma_pop_stack h_active',
     &    0, MA_ERR)
      if (.not. ma_pop_stack(h_tmp_id1) ) call errquit(
     &    'cons_load_fixed_atoms:ma_pop_stack h_tmp_id1',
     &    0, MA_ERR)

      if (.not. ma_pop_stack(h_tmp_id) ) call errquit(
     &    'cons_load_fixed_atoms:ma_pop_stack h_tmp_id',
     &    0, MA_ERR)



      end
c
      subroutine cons_add_egrad(rtdb,energy,gx)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "cons.fh"
c     
      integer rtdb
      double precision energy
      double precision gx(*)
      character*255 namespace
c
      call cons_get_namespace(namespace,rtdb)
      if(.not.cons_status(namespace,rtdb)) return
      if(.not.cons_hbonds(namespace,rtdb)) return
c
      call cons_load_geom(rtdb)
c
      call cons_add_hbond_egrad(rtdb,energy,gx)
c
      call cons_unload_geom()
c
      return
      end
c
      subroutine cons_add_energy(rtdb,energy)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "cons.fh"
c     
      integer rtdb
      double precision energy
      character*255 namespace
c
      call cons_get_namespace(namespace,rtdb)
      if(.not.cons_status(namespace,rtdb)) return
      if(.not.cons_hbonds(namespace,rtdb)) return
c
      call cons_load_geom(rtdb)
c
      call cons_add_hbond_energy(rtdb,energy)
c
      call cons_unload_geom()
      return
      end
c
      subroutine cons_print_summary(rtdb)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "util.fh"
#include "cons.fh"
c     
      integer rtdb
      double precision energy
c     local variables
c     ---------------
      logical status
      logical master
      integer i
      integer inb,nhb
      integer iat,jat
      double precision r
      double precision r0,k
      integer i_c,i_rhb,i_khb,i_ijhb
      double precision etot,e,f(3)
      character*255 namespace
c
      if(.not. rtdb_get(rtdb, 'cons:active', mt_log, 1, status)) then
        return
      else if (.not.status) then
        return
      end if

      master= ga_nodeid().eq.0
      if(master) then
        call cons_print_header()
        call util_print_centered(6,
     >     "Summary",
     >     1,.true.)

      end if

c
      call cons_load_geom(rtdb)
      call cons_get_namespace(namespace,rtdb)
      if(cons_fixed_atoms(namespace,rtdb)) then
        call cons_info_fixed_atoms(namespace,rtdb) 
      end if
c
      call cons_get_hbond_nhb(nhb)
      
      if(nhb.ne.0) then
        i_c    = cons_get_i_c()
        i_ijhb = cons_get_i_hbond_id()
        i_khb  = cons_get_i_hbond_k()
        i_rhb  = cons_get_i_hbond_r()
c
        etot=0.0d0
        do inb=1,nhb
         iat=int_mb(i_ijhb+2*(inb-1))
         jat=int_mb(i_ijhb+2*(inb-1)+1)
         r0 =dbl_mb(i_rhb+inb-1)
         k  =dbl_mb(i_khb+inb-1)
c
         call cons_spring_force(k,r0,
     >             dbl_mb(i_c+(iat-1)*3),
     >             dbl_mb(i_c+(jat-1)*3),
     >             r,e,f)

         if(master) then
          write(6,'(A30,T32,I4)')"spring #                    :",inb
          write(6,'(A30,T32,2I4,2(F6.3,2X))')
     >      "spring parameters (i,j,k,r0):",iat,jat,k,r0
          write(6,'(A30,T32,F6.3)')"spring length               :",r
          write(6,'(A30,T32,F6.3)')"spring energy               :",e
          write(6,'(A30,T32,3(F6.3,2X))')
     >          "spring forces               :",
     >          (f(i),i=1,3)
          write(6,*)"  "
         end if
         etot=etot+e
        end do
      else
       etot = 0
      end if

      if (.not. rtdb_get(rtdb,'task:energy', mt_dbl, 1, energy))
     $        call errquit('driver: could not get energy',0, RTDB_ERR)

      if(master) then
      write(6,'(A39,3X,F12.6)')
     >       " Total energy with constraints        :",
     >         energy
      write(6,'(A39,3X,F12.6)')
     >       " Energy contribution from constraints :",
     >                       etot
      write(6,'(A39,3X,F12.6)')
     >       " Total energy without constraints     :",
     >                       energy-etot
      write(6,*) "--------------------------------------"
      end if
      call cons_unload_geom()

      return
      end  
c
