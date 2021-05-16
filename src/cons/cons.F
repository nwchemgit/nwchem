c $I#d: cons.F,v 1.1 2004/01/28 01:30:59 marat Exp $
      function cons_active(rtdb)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "geom.fh"
#include "cons_data.fh"
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
#include "cons_data.fh"
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
#include "cons_data.fh"
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
#include "stdio.fh"
#include "cons_data.fh"
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
     >    write(luout,*) "entered cons_init"
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
         call banner(luout,"Warning!!!                               ",
     >              "*" ,
     >             .true.,.false.,.true. )
         call banner(luout,"geometry:actlist is an obsolete option   ",
     >              "*" ,
     >             .false.,.false.,.true. )
         call banner(luout,"consider using fixed atoms constraints   ",
     >              "*" ,
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
         call banner(luout,"Warning!!!                            ","*",
     >           .true.,.false.,.true. )
         call banner(luout,"mixing both geometry:actlist  and     ","*",
     >             .false.,.false.,.true. )
         call banner(luout,"fixed atom constraints is not         ","*",
     >             .false.,.false.,.true. )
         call banner(luout,"recommended. Check output below for   ","*",
     >             .false.,.false.,.true. )
         call banner(luout,"for potential problems.               ","*",
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
      status=geom_amatrix_get(geom,unita)
      status=status.and.geom_systype_get(geom,isystype)

      status=status.and.geom_destroy(geom)
      if(.not.status) 
     & call errquit('cons_init: geom_destroy?',70, GEOM_ERR)

      call cons_load_fixed_atoms(namespace,rtdb)

      call cons_print_gen_info(namespace,rtdb)
c
      if(oprint) 
     >    write(luout,*) "finished cons_init"
c
      call util_print_pop()

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
#include "stdio.fh"
#include "cons_data.fh"
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

c      if(oprint) 
c     >    write(*,*) "entered cons_init"
      if(oprint) 
     >    write(luout,*) "entered cons_init"
c
      ocons = .false.
c
c     initialize harmonic bonds
      if(cons_hbonds(namespace,rtdb)) then
      ocons = .true.
      call cons_set_hbond_nhb(0)
      call cons_load_hbonds(namespace,rtdb)
      end if
c
c
c     initialize dihedrals
      if(cons_hdihed(namespace,rtdb)) then
      ocons = .true.
      call cons_set_hdihed_n(0)
      call cons_load_hdihed(namespace,rtdb)
      end if
c
c
c     initialize bondings
      if(cons_hbondings(namespace,rtdb)) then
      ocons = .true.
      call cons_set_hbondings_nhc(0)
      call cons_load_hbondings(namespace,rtdb)
      end if
c
c
c     initialize pbondings
      if(cons_pbondings(namespace,rtdb)) then
      ocons = .true.
      call cons_set_pbondings_nhp(0)
      call cons_load_pbondings(namespace,rtdb)
      end if


c

c      call cons_process_old_style(namespace,rtdb)
      call cons_load_fixed_atoms(namespace,rtdb)

      call cons_print_gen_info(namespace,rtdb)
c
      if(oprint) 
     >    write(luout,*) "finished cons_init"
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
#include "cons_data.fh"
#include "cons.fh"
#include "global.fh"
      integer rtdb
c
      call cons_print_summary(rtdb)
      call cons_unload_hbonds()
      call cons_unload_hbondings()
      call cons_unload_pbondings()
      call cons_unload_hdihed()
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
      real*8 unita(3,3)
      integer isystype
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
      status=status.and.geom_amatrix_get(geom,unita)
      status=status.and.geom_systype_get(geom,isystype)
      if(.not.status) 
     & call errquit('cons_init: geom_create?',70, GEOM_ERR)
 
      status=geom_destroy(geom)
      if(.not.status) 
     & call errquit('cons_init: geom_destroy?',70, GEOM_ERR)
      call cons_set_h_c(h_c)
      call cons_set_unita(unita)
      call cons_set_isystype(isystype)
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
          if(id.gt.nat) call errquit(
     $         ' cons: fixed atom not existing ',id,INPUT_ERR)
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
      subroutine cons_add_egrad(rtdb,energy,n,gx)
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
      integer n
      double precision gx(n)
      character*255 namespace
      logical obonds,odihedrals,ocons,obondings,opbondings
      logical osim
c
      call cons_get_namespace(namespace,rtdb)
      if(.not.cons_status(namespace,rtdb)) return
      obonds = cons_hbonds(namespace,rtdb)
      obondings = cons_hbondings(namespace,rtdb)
      opbondings = cons_pbondings(namespace,rtdb)
      odihedrals = cons_hdihed(namespace,rtdb)
      ocons = obonds.or.odihedrals.or.obondings
      if(.not.ocons) return
c
      call cons_load_geom(rtdb)
c      
      if(.not.rtdb_get(rtdb, 'cons:simulate', mt_log, 1, osim)) then
        osim = .false.
      end if
c
      if(osim) then
        call dfill(0.0d0,gx,1,n)
        energy = 0.0d0
      end if
c
      if(obonds)  
     >    call cons_add_hbond_egrad(rtdb,energy,gx)
      if(obondings)  
     >    call cons_add_hbondings_egrad(rtdb,energy,gx)
      if(opbondings)  
     >    call cons_add_pbondings_egrad(rtdb,energy,gx)
      if(odihedrals)     
     >   call cons_add_hdihed_egrad(rtdb,energy,gx)
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
      logical obonds,odihedrals,obondings,opbondings,ocons
      logical osim
c
      call cons_get_namespace(namespace,rtdb)
      if(.not.cons_status(namespace,rtdb)) return
      obonds = cons_hbonds(namespace,rtdb)
      obondings = cons_hbondings(namespace,rtdb)
      opbondings = cons_pbondings(namespace,rtdb)
      odihedrals = cons_hdihed(namespace,rtdb)
      ocons = obonds.or.odihedrals.or.obondings
      if(.not.ocons) return
c
      call cons_load_geom(rtdb)
      if(.not.rtdb_get(rtdb, 'cons:simulate', mt_log, 1, osim)) then
        osim = .false.
      end if
      if(osim) energy = 0.0d0
c
      if(obonds)  
     >   call cons_add_hbond_energy(rtdb,energy)
      if(obondings)  
     >   call cons_add_hbondings_energy(rtdb,energy)
      if(opbondings)  
     >   call cons_add_pbondings_energy(rtdb,energy)
      if(odihedrals)     
     >   call cons_add_hdihed_energy(rtdb,energy)
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
#include "stdio.fh"
#include "cons.fh"
c     
      integer rtdb
      double precision energy
c     local variables
c     ---------------
      logical status
      logical master
      integer i,ii,n,j
      integer inb,nhb,nhc,nhp
      integer iat,jat
      double precision r
      double precision r0,k,K0,gcut0,gamma0,c0,gamm
      integer i_c,i_rhb,i_khb,i_ijhb
      integer i_n,i_indx,i_coef,i_k0,i_gcut0,i_gamma0,i_c0
      double precision etot,e,f(3)
      character*255 namespace
      logical  cons_active
      external cons_active
c
      if(.not.cons_active(rtdb)) then
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
      
      etot = 0.0d0
      if(nhb.ne.0) then
        i_c    = cons_get_i_c()
        i_ijhb = cons_get_i_hbond_id()
        i_khb  = cons_get_i_hbond_k()
        i_rhb  = cons_get_i_hbond_r()
c
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
          write(luout,'(A30,T32,I4)')"spring #                    :",inb
          write(luout,'(A30,T32,2I4,2(F6.3,2X))')
     >      "spring parameters (i,j,k,r0):",iat,jat,k,r0
          write(luout,'(A30,T32,F6.3)')"spring length               :",r
          write(luout,'(A30,T32,F6.3)')"spring energy               :",e
          write(luout,'(A30,T32,3(F6.3,2X))')
     >          "spring forces               :",
     >          (f(i),i=1,3)
          write(luout,*)"  "
         end if
         etot=etot+e
        end do
      end if

      call cons_get_hbondings_nhc(nhc)

      if(nhc.ne.0) then
         i_c      = cons_get_i_c()
         i_n      = cons_get_i_hbondings_n0()
         i_indx   = cons_get_i_hbondings_indx()
         i_coef   = cons_get_i_hbondings_coef()
         i_k0     = cons_get_i_hbondings_k0()
         i_gamma0 = cons_get_i_hbondings_gamma0()

         ii = 0
         do i=1,nhc
            n      = int_mb(i_n+i-1)
            K0     = dbl_mb(i_k0+i-1)
            gamma0 = dbl_mb(i_gamma0+i-1)

            call cons_bondings_energy(n,int_mb(i_indx+2*ii),
     >                               dbl_mb(i_coef+ii),
     >                               K0,gamma0,dbl_mb(i_c),gamm,e)
            if(master) then
               write(luout,'(A,I5)')" bondings spring #            :",i
               write(luout,'(A,2(F12.6,2X))')
     >        " spring parameters (K0,gamma0):",K0,gamma0
               write(luout,'(A,F12.6)')
     >         " gamma                        :",gamm
               write(luout,'(A,F12.6)')
     >         " spring energy                :",e
               write(luout,'(A)')" coefficient index1 index2    :"
               do j=1,n
                  write(luout,'(F12.6,2I7)') dbl_mb(i_coef+ii+j-1),
     >                                   int_mb(i_indx+2*ii+2*(j-1)),
     >                                   int_mb(i_indx+2*ii+2*(j-1)+1)
               end do
               write(luout,*)"  "
            end if

            ii = ii + n
            etot=etot+e
         end do
      end if


      call cons_get_pbondings_nhp(nhp)

      if(nhp.ne.0) then
         i_c      = cons_get_i_c()
         i_n      = cons_get_i_pbondings_n0()
         i_indx   = cons_get_i_pbondings_indx()
         i_coef   = cons_get_i_pbondings_coef()
         i_k0     = cons_get_i_pbondings_k0()
         i_gcut0  = cons_get_i_pbondings_gcut0()
         i_gamma0 = cons_get_i_pbondings_gamma0()
         i_c0     = cons_get_i_pbondings_c0()

         ii = 0
         do i=1,nhp
            n      = int_mb(i_n+i-1)
            K0     = dbl_mb(i_k0+i-1)
            gcut0  = dbl_mb(i_gcut0+i-1)
            gamma0 = dbl_mb(i_gamma0+i-1)
            c0     = dbl_mb(i_c0+i-1)

            call cons_pbondings_energy(n,int_mb(i_indx+2*ii),
     >                               dbl_mb(i_coef+ii),
     >                               K0,gcut0,gamma0,c0,
     >                               dbl_mb(i_c),gamm,e)
            if(master) then
               write(luout,'(A,I5)')" bondings penalty#            :",i
               write(luout,'(A,4(F12.6,2X))')
     >     " spring parameters (K0,gcut0,gamma0,c0):",K0,gcut0,gamma0,c0
               write(luout,'(A,F12.6)')
     >         " gamma                        :",gamm
               write(luout,'(A,F12.6)')
     >         " penalty energy               :",e
               write(luout,'(A)')" coefficient index1 index2    :"
               do j=1,n
                  write(luout,'(F12.6,2I7)') dbl_mb(i_coef+ii+j-1),
     >                                   int_mb(i_indx+2*ii+2*(j-1)),
     >                                   int_mb(i_indx+2*ii+2*(j-1)+1)
               end do
               write(luout,*)"  "
            end if

            ii = ii + n
            etot=etot+e
         end do
      end if




      if (.not. rtdb_get(rtdb,'task:energy', mt_dbl, 1, energy))
     $   energy = 0.0d0

      if ((master).and.(dabs(energy).gt.1.0d-9)) then
      write(luout,'(A39,3X,F12.6)')
     >       " Total energy with constraints        :",
     >         energy
      write(luout,'(A39,3X,F12.6)')
     >       " Energy contribution from constraints :",
     >                       etot
      write(luout,'(A39,3X,F12.6)')
     >       " Total energy without constraints     :",
     >                       energy-etot
      write(luout,*) "--------------------------------------"
      end if
      call cons_unload_geom()

      return
      end  
c
c $Id$
