c $I#d: cons.F,v 1.1 2004/01/28 01:30:59 marat Exp $
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
      logical old_style
      logical new_style
      logical status
      integer geom
      integer ma_type
      integer i
      integer nat
      integer nact
      integer n
      integer h_act,i_act
      integer h_tmp_id,i_tmp_id 
      integer h_tmp_act,i_tmp_act
      character(255) tag_id
      character(255) tag_n
      character(255) namespace
      logical oprint
c
      call util_print_push()
      call util_print_rtdb_load(rtdb,'cons')
      oprint = util_print('information', print_debug)
      oprint = oprint .and. cons_print()

      if(oprint) 
     >    write(*,*) "entered cons_init"
      if(.not. rtdb_cget(rtdb,"constraints" , 1, namespace)) then
         namespace=cons_default
      end if

      if(.not.cons_activate(namespace,rtdb)) then
        status = rtdb_put(rtdb, 'cons:active', mt_log, 1, .false.)
        if(.not.status)
     &  call errquit('cons_init',0,RTDB_ERR) 
        write(*,*) "cons_init:no constraints found "//namespace
        write(*,*) "cons_init: constraints will be deactivated"
        call errquit("cons_init:no constraints found "//namespace,
     >               0,RTDB_ERR) 
        return
      end if
c
c     initialize harmonic bobds
      call cons_set_hbond_nhb(0)
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

      if(.not.rtdb_get(rtdb,tag_n,mt_int,1,n)) then
         n = 0
      end if

      if (n.eq.0) then
        new_style=.false.
      else
        new_style=.true.
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
         call banner(6,"reccomended. Check output below for   ","*" ,
     >             .false.,.false.,.true. )
         call banner(6,"for potential problems.               ","*" ,
     >             .false.,.true.,.true. )
          end if
          call cons_unfix_atoms(nact,int_mb(i_act),rtdb)
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

      call cons_print_gen_info(rtdb)
c
      if(oprint) 
     >    write(*,*) "finished cons_init"
c
      call util_print_pop()

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
      integer i
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

      subroutine cons_unload_hbonds()
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "cons.fh"
      integer nhb,h_rhb,h_khb,h_ijhb

      h_ijhb = cons_get_h_hbond_id()
      h_khb  = cons_get_h_hbond_k()
      h_rhb  = cons_get_h_hbond_r()
c
c     unload harmonic constraints if any
c
      call cons_get_hbond_nhb(nhb)
      if(nhb.gt.0) then
      if (.not. ma_free_heap(h_rhb) ) call errquit(
     &    'cons_bond_input: unable to free h_rhb',
     &    0, MA_ERR)
      if (.not. ma_free_heap(h_khb) ) call errquit(
     &    'cons_bond_input: unable to free h_khb',
     &    0, MA_ERR)
      if (.not. ma_free_heap(h_ijhb) ) call errquit(
     &    'cons_bond_input: unable to free h_ijhb',
     &    0, MA_ERR)
      call cons_set_hbond_nhb(0)
      end if
      end
c
      subroutine cons_load_hbonds(namespace,rtdb)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "geom.fh"
c     
      character*(*) namespace
      integer rtdb
c
      logical status
      integer nhb
      integer i_rhb,i_khb,i_ijhb
      integer h_rhb,h_khb,h_ijhb
      character(255) tag_id
      character(255) tag_n
      character(255) tag_r
      character(255) tag_k
c
      call cons_unload_hbonds()
c
      call cons_hbond_id_tag(namespace,tag_id)
      call cons_hbond_n_tag(namespace,tag_n)
      call cons_hbond_k_tag(namespace,tag_k)
      call cons_hbond_r_tag(namespace,tag_r)
c
c     load harmonic constraints
c
      status=rtdb_get(rtdb,tag_n,mt_int,1,nhb)
c     
c     exit if no harm bonds
      if(.not.status .or. nhb.eq.0) then
c            call errquit(
c     >     'cons_load_hbonds: unable to get number of harm bonds:',
c     >      nhb, MA_ERR)
      return
      end if
 
      if ( .not. ma_alloc_get( MT_INT, 2*nhb, tag_id,
     &      h_ijhb, i_ijhb) ) call errquit(
     &      'cons_load_hbonds: unable to allocate cons scratch space',
     &      2*nhb, MA_ERR)
      if ( .not. ma_alloc_get( MT_DBL, nhb, tag_k,
     &      h_khb, i_khb) ) call errquit(
     &      'cons_load_hbonds: unable to allocate cons scratch space',
     &      nhb, MA_ERR)
      if ( .not. ma_alloc_get( MT_DBL, nhb, tag_r,
     &      h_rhb, i_rhb) ) call errquit(
     &      'cons_load_hbonds: unable to allocate cons scratch space',
     &      nhb, MA_ERR)

        status = rtdb_get(rtdb,tag_id,
     >     mt_int,2*nhb,int_mb(i_ijhb)) 
        status = status . and. rtdb_get(rtdb,tag_k,
     >     mt_dbl,nhb,dbl_mb(i_khb)) 
        status = status . and. rtdb_get(rtdb,tag_r,
     >     mt_dbl,nhb,dbl_mb(i_rhb)) 
        if(.not.status) call errquit(
     >     'cons_load_hbonds: unable to get ma index for harm bonds',
     >      nhb, MA_ERR)
       
       call cons_set_hbond_nhb(nhb)
       call cons_set_h_hbond_id(h_ijhb)
       call cons_set_h_hbond_k(h_khb)
       call cons_set_h_hbond_r(h_rhb)

      end
c
      subroutine cons_load_constraints(namespace,rtdb)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "geom.fh"
      integer rtdb
      character*(*) namespace
      
      call cons_set_namespace(namespace)

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
      character(255) tag_id
      character(255) tag_n
      integer h_tmp_id,i_tmp_id 
      integer h_tmp_id1,i_tmp_id1 
      integer h_active,i_active
      integer nat
      integer geom
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
        write(*,*) "no fixed atoms found"
        return
      end if

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
            int_mb(i_tmp_id1-j+1)=i  
          end if
        end do
        
      status = rtdb_put(rtdb,'geometry:actlist',
     >     mt_int,nat-nfa,int_mb(i_tmp_id1))
        if(.not.status) call errquit(
     >     'cons_load_fixed_atoms: unable to store fixed atoms',
     >      0, MA_ERR)

      if (.not. ma_pop_stack(h_active) ) call errquit(
     &    'cons_load_fixed_atoms:ma_pop_stack',
     &    0, MA_ERR)
      if (.not. ma_pop_stack(h_tmp_id1) ) call errquit(
     &    'cons_load_fixed_atoms:ma_pop_stack',
     &    0, MA_ERR)

      if (.not. ma_pop_stack(h_tmp_id) ) call errquit(
     &    'cons_load_fixed_atoms:ma_pop_stack',
     &    0, MA_ERR)



      end
c
      subroutine cons_delete_hbonds(namespace,rtdb)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
c     
      character*(*) namespace
      integer rtdb
c
      logical status
      character(255) tag_id
      character(255) tag_n
      character(255) tag_r
      character(255) tag_k

      call cons_hbond_id_tag(namespace,tag_id)
      call cons_hbond_n_tag(namespace,tag_n)
      call cons_hbond_k_tag(namespace,tag_k)
      call cons_hbond_r_tag(namespace,tag_r)

      status =              rtdb_delete(rtdb,tag_n)
      status = status .and. rtdb_delete(rtdb,tag_id)
      status = status .and. rtdb_delete(rtdb,tag_k)
      status = status .and. rtdb_delete(rtdb,tag_r)
      if(.not.status) call errquit(
     >     'input_set: unable to delete prior harm bonds',
     >      0, MA_ERR)
      end
c
      subroutine cons_delete_fixed_atoms(rtdb)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
c     
      integer rtdb
c
      logical status
      status = rtdb_delete(rtdb,"cons:fpos_i")
      if(.not.status) call errquit(
     >     'input_set: unable to delete prior fixed atoms',
     >      0, MA_ERR)
      end
c
      subroutine cons_add_egrad(rtdb,energy,gx)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      integer rtdb
      double precision energy
      double precision gx(*)
      logical status
      integer nhb
      character(255) tag_n
      character(255) namespace
c
      if(.not. rtdb_get(rtdb, 'cons:active', mt_log, 1, status)) then
        return
      else if (.not.status) then
        return
      end if
c
      call cons_load_geom(rtdb)
      call cons_get_namespace(namespace)
      call cons_load_hbonds(namespace,rtdb)
      call cons_get_hbond_nhb(nhb)
      if(nhb .gt. 0 ) then
         call cons_add_hbond_egrad(rtdb,energy,gx)
      end if
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
c     
      integer rtdb
      double precision energy
      logical status
      integer nhb
      character(255) tag_n
      character(255) namespace
c
      if(.not. rtdb_get(rtdb, 'cons:active', mt_log, 1, status)) then
        return
      else if (.not.status) then
        return
      end if

      call cons_load_geom(rtdb)
      call cons_get_namespace(namespace)
      call cons_load_hbonds(namespace,rtdb)
      call cons_get_hbond_nhb(nhb)
      if(nhb .gt. 0 ) then
           call cons_add_hbond_energy(rtdb,energy)
      end if
c
      call cons_unload_geom()
      return
      end
c
      subroutine cons_add_hbond_energy(rtdb,energy)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "cons.fh"
c     
      integer rtdb
      double precision energy
      logical status
      logical oprint
      integer i
      integer inb
      integer iat,jat
      double precision r
      double precision ciat(3),cjat(3)
      double precision r0,k
      double precision e,enb,enb1
      integer i_c,i_rhb,i_khb,i_ijhb
      integer nhb
c
      call util_print_push()
      call util_print_rtdb_load(rtdb,'cons')
      oprint = util_print('energy', print_default)
      oprint =oprint .and. (ga_nodeid().eq.0)
c
      call cons_get_hbond_nhb(nhb)
      i_c    = cons_get_i_c()
      i_ijhb = cons_get_i_hbond_id()
      i_khb  = cons_get_i_hbond_k()
      i_rhb  = cons_get_i_hbond_r()

      e=0.0d0
      do inb=1,nhb
       iat=int_mb(i_ijhb+2*(inb-1))
       jat=int_mb(i_ijhb+2*(inb-1)+1)
       r0 =dbl_mb(i_rhb+inb-1)
       k  =dbl_mb(i_khb+inb-1)

       call cons_spring_energy(k,r0,
     >         dbl_mb(i_c+(iat-1)*3),
     >         dbl_mb(i_c+(jat-1)*3),
     >         r,
     >         enb)

       if(oprint) then
        write(6,*)"adding spring # ",inb
        write(6,*)"  spring parameters (i,j,k,r0):",iat,jat,k,r0
        write(6,*)"  spring length and energy    :",r,enb
       end if
       e=e+enb
      end do
      if(rtdb_get(rtdb, 'cons:simulate', mt_log, 1, status)) then
        write(*,*) "cons energy simulation"
        energy=e
      else 
        energy = energy + e
      end if
c
      call util_print_pop()
      return
      end
c
      subroutine cons_add_hbond_egrad(rtdb,energy,gx)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "util.fh"
#include "cons.fh"
c     
      integer rtdb
      double precision gx(*)
      double precision energy
c     local variables
c     ---------------
      logical status
      logical oprint
      logical pdefault
      integer i
      integer inb,nhb
      integer iat,jat
      double precision r
      double precision ciat(3),cjat(3)
      double precision r0,k
      double precision e,f(3)
      integer i_c,i_rhb,i_khb,i_ijhb
c
      call util_print_push()
      call util_print_rtdb_load(rtdb,'cons')
      oprint = util_print('energy', print_default)
      oprint =oprint .and. (ga_nodeid().eq.0)
c
      call cons_get_hbond_nhb(nhb)
      i_c    = cons_get_i_c()
      i_ijhb = cons_get_i_hbond_id()
      i_khb  = cons_get_i_hbond_k()
      i_rhb  = cons_get_i_hbond_r()

c
      e=0.0d0
      do inb=1,nhb
       iat=int_mb(i_ijhb+2*(inb-1))
       jat=int_mb(i_ijhb+2*(inb-1)+1)
       r0 =dbl_mb(i_rhb+inb-1)
       k  =dbl_mb(i_khb+inb-1)
c
       call cons_spring_force(k,r0,
     >           dbl_mb(i_c+(iat-1)*3),
     >           dbl_mb(i_c+(jat-1)*3),
     >           r,e,f)

       if(oprint) then
        write(6,*)"adding spring # ",inb
        write(6,*)"  spring parameters (i,j,k,r0):",iat,jat,k,r0
        write(6,*)"  spring length and energy    :",r,e
        write(6,*)"  spring forces               :",(f(i),i=1,3)
       end if

      if(rtdb_get(rtdb, 'cons:simulate', mt_log, 1, status)) then
       do i=1,3
         gx((iat-1)*3+i)=f(i)
         gx((jat-1)*3+i)=-f(i)
       end do
       energy = e
      else 
       do i=1,3
         gx((iat-1)*3+i)=gx((iat-1)*3+i)+f(i)
         gx((jat-1)*3+i)=gx((jat-1)*3+i)-f(i)
       end do
       energy = energy+e
      end if

      end do
      call util_print_pop()

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
      double precision ciat(3),cjat(3)
      double precision r0,k
      integer i_c,i_rhb,i_khb,i_ijhb
      double precision etot,e,f(3)
      character(255) namespace
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
      call cons_get_namespace(namespace)
      call cons_load_hbonds(namespace,rtdb)
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
          write(6,*)"  spring #                    : ",inb
          write(6,*)"  spring parameters (i,j,k,r0):",iat,jat,k,r0
          write(6,*)"  spring length               :",r
          write(6,*)"  spring energy               :",e
          write(6,*)"  spring forces               :",(f(i),i=1,3)
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
      write(6,*)" Total energy with constraints        :",energy
      write(6,*)" Energy contribution from constraints :",etot
      write(6,*)" Total energy without constraints  :",energy-etot
      write(6,*)" -----------------------------------"
      end if
      call cons_unload_hbonds(rtdb)
      call cons_unload_geom()

      return
      end  
c
      subroutine cons_spring_force(k,r0,r1,r2,r,energy,f)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "cons_data.fh"
c     
      double precision k,r0
      double precision r1(3)
      double precision r2(3)
      double precision r
      double precision energy
      double precision f(3)
c
      integer i
c
       r=(r1(1)-r2(1))**2+
     >   (r1(2)-r2(2))**2+
     >   (r1(3)-r2(3))**2
       r=sqrt(r)
       energy=k*(r-r0)**2

       do i=1,3
         f(i)=2*k*(r-r0)*
     &     (r1(i)-r2(i))/r
       end do

      return
      end

      subroutine cons_spring_energy(k,r0,r1,r2,r,energy)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "cons_data.fh"
c     
      double precision k,r0
      double precision r1(3)
      double precision r2(3)
      double precision r
      double precision energy
c
      integer i
c
       r=(r1(1)-r2(1))**2+
     >   (r1(2)-r2(2))**2+
     >   (r1(3)-r2(3))**2
       r=sqrt(r)
       energy=k*(r-r0)**2
      return
      end

