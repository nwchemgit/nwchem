c
c $Id$
c
c
C> \ingroup task
C> @{
      function task_qmmm_gradient(rtdb)
      implicit none
c
#include "rtdb.fh"
#include "util.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm_params.fh"
#include "qmmm.fh"
c
      integer rtdb
      logical  task_qmmm_gradient
c
      logical  qmmm_energy_gradient
      external qmmm_energy_gradient
      logical  qmmm_neb_gradient
      external qmmm_neb_gradient
c     local variables
      character*32 interface
      logical save_rst
      logical do_print
      character*30 operation
      double precision energy
c
      if (.not. rtdb_cget(rtdb, 'task:operation', 1, operation))
     $     operation = 'energy'
c
      save_rst = .false.
      if (operation.eq.'optimize' .or.
     >    operation.eq.'saddle') save_rst = .true.
c
      do_print = .false.
      if (operation.eq.'gradient') do_print = .true.
c
      interface = qmmm_get_interface()
c
      task_qmmm_gradient = .true.
      if(operation.eq."neb") then
        task_qmmm_gradient = 
     >     qmmm_neb_gradient(rtdb)
        goto 20
      end if
      if (interface.eq.'mm') then
         call qmmm_geom_push_active(rtdb)
         call md_sp()
         if (.not. rtdb_get(rtdb,'md:energy',mt_dbl,1,energy))
     $     call errquit('qmmm: failed getting  md energy', 0, RTDB_ERR)
c
         if (.not. rtdb_put(rtdb,'qmmm:energy',mt_dbl,1,energy/cau2kj))
     $     call errquit('qmmm: failed put energy', 0, RTDB_ERR)
c
         call qmmm_geom_get_grad(rtdb)
         call qmmm_print_energy2(rtdb)

      else if (interface.eq.'qm') then
         call qmmm_geom_push_active(rtdb)
         call md_sp_qmmm()
         task_qmmm_gradient = qmmm_energy_gradient(rtdb,.true.)
         call qmmm_geom_get_grad(rtdb)
         call qmmm_energy_rtdb_push(rtdb)
         call qmmm_print_energy(rtdb)
         if(save_rst)
     >    call mm_write_restart()
         if(do_print)
     >     call mm_print_info()
      endif
c
 20   continue
c
      end
C> @}

      function qmmm_energy_gradient(rtdb,do_grad)
      implicit none
c
#include "rtdb.fh"
#include "util.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "nwc_const.fh"
#include "geom.fh"
#include "qmmm_params.fh"
#include "qmmm.fh"
#include "bq.fh"
c
      integer rtdb
      logical do_grad
      logical qmmm_energy_gradient
      logical status
      logical task_energy_doit
      external task_energy_doit
      logical task_gradient_doit
      external task_gradient_doit
      character*255 theory
      double precision qm_energy
      integer i_g
      integer nfg
      integer i_ifg
      integer nact
      logical ignore
      logical frozen_density,esp
      logical numerical
      double precision ebq
      character*30 pname

      pname = "qmmm_energy_gradient"

      if(qmmm_print_debug()) then
         if(do_grad) then
           write(*,*) pname,"doing gradient"
         else
           write(*,*) pname,"doing energy" 
         end if
      end if
      status = .true.

      if (.not. rtdb_get(rtdb, 'task:numerical', mt_log, 1,
     $     numerical)) numerical = .false.

      call qmmm_bq_data_update_active()

      nfg = qmmm_geom_get_nfg()
      i_ifg = qmmm_geom_get_i_ifg()
      i_g = qmmm_geom_get_i_g()

      if (.not. rtdb_get(rtdb,'qmmm:frozen_density',
     $                   mt_log,1,frozen_density))
     $     frozen_density=.false.
c
      if (.not. rtdb_get(rtdb,'qmmm:esp_density',
     $                   mt_log,1,esp))
     $     esp=.false.
      
      if(esp) frozen_density=.true.
c
c     disable frozen density/esp if any of qm atoms are active
c     ----------------------------------------------------
      if(qmmm_master()) 
     >   write(*,*) "frozen_density,esp",frozen_density,esp
      if(frozen_density.or.esp) then
        call qmmm_cons_get_nacts_qm(nact)
        if(nact.gt.0) then
          if(qmmm_master()) 
     >       write(*,*) "disabling frozen_density,esp "//
     >                  "due to active qm atoms"
           frozen_density = .false.
           esp = .false.
        end if
      end if
      
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit('task:energy: theory not specified',0, RTDB_ERR)

c
c     create full geometry
c     -------------------
      call qmmm_geom_create_full(rtdb)
c
c     update coordinates of classical charges
c     ---------------------------------------
      call qmmm_bq_coord_update()
c
c     create active atom list for qm atoms
c     ------------------------------------
      if(do_grad) then
        call qmmm_map_fixed_atoms(rtdb,nact)
      else
        nact = 0
      end if
c
      if(qmmm_print_debug()) then
        write(*,*) "number of active qm atoms",nact
      end if
c
      ignore = rtdb_delete(rtdb,'qmmm:bq_nuc_energy')
c
c
c     qm atoms gradients if any
c     -------------------------
      qm_energy=0.0d0
      call dfill(3*nfg,0,dbl_mb(i_g),1)
      if(nact.ne.0) then
        if(numerical)
     >  call qmmm_geom_init_aux(rtdb,"geometry")
        status = task_gradient_doit(rtdb,theory,qm_energy,dbl_mb(i_g))
        if(numerical)
     >  call qmmm_geom_aux_end(rtdb)
        if(.not.bq_nuc_energy(rtdb,ebq))
     >    call errquit('failed bq energy',0, RTDB_ERR)
        if (.not. rtdb_put(rtdb,'qmmm:bq_nuc_energy',mt_dbl,1,ebq))
     $     call errquit('qmmm: failed put bq energy', 0, RTDB_ERR)
        call qmmm_links_adjust_forces(nfg,
     &                       int_mb(i_ifg),
     &                       dbl_mb(i_g))

        call qmmm_cons_reaction(rtdb,.true.,qm_energy,dbl_mb(i_g))

        call mm_add_solute_force(nfg,
     &                     int_mb(i_ifg),
     &                     dbl_mb(i_g))

      else
        if(nfg.ne.0) then  
         if(frozen_density) then
           if(qmmm_print_debug()) then
              write(*,*) "frozen qm density"
           end if
         call qmmm_bq_energy(qm_energy,rtdb,esp)
         else
         status = task_energy_doit(rtdb,theory,qm_energy)
         call qmmm_cons_reaction(rtdb,.false.,qm_energy,dbl_mb(i_g))
         if(.not.bq_nuc_energy(rtdb,ebq))
     >     call errquit('failed bq energy',0, RTDB_ERR)
         if (.not. rtdb_put(rtdb,'qmmm:bq_nuc_energy',mt_dbl,1,ebq))
     $     call errquit('qmmm: failed put bq energy', 0, RTDB_ERR)
         end if
         end if
      end if
c
c     gradients on classical atoms if any
c     -----------------------------------
      if(do_grad.and.nfg.ne.0) 
     >  call qmmm_bq_gradient(rtdb,esp)
c
      if (.not. rtdb_put(rtdb,'qmmm:qm_energy',mt_dbl,1,qm_energy))
     $     call errquit('qmmm: failed put energy', 0, RTDB_ERR)
c
      call qmmm_geom_restore(rtdb)
c
      if(do_grad)
     > call qmmm_restore_fixed_atoms(rtdb)
c
      qmmm_energy_gradient=status
c
      call mm_print_system()
c
      end

      function qmmm_neb_gradient(irtdb)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "global.fh"
#include "stdio.fh"
#include "util.fh"
#include "qmmm.fh"
#include "inp.fh"
      integer irtdb
      logical qmmm_neb_gradient
c     
      logical  status
      logical  loadrst
c     
      character*32 pname
      character*255 prefix
      character*255 rstfile
      character*255 rstfile0
      character*255 neb_region
      integer ibead
c     
      logical   qmmm_mixed_optimization
      external  qmmm_mixed_optimization
      logical   qmmm_energy_gradient
      external  qmmm_energy_gradient

      pname = "qmmm_neb_gradient"
      if(qmmm_print_debug()) 
     > write(luout,*) "in ", pname
c 
      qmmm_neb_gradient = .true.
c
c     which bead we are doing here
c     ---------------------------
      if(.not.rtdb_get(irtdb,"neb:ibead",mt_int,1,ibead))
     >  call errquit('neb:ibead put',0,0)
c
c     construct restart file name  
c     ---------------------------
      call mm_system_prefix(prefix)
      write(rstfile,13) prefix(1:inp_strlen(prefix)),ibead
c      call util_file_name_resolve(rstfile,.false.)
c
c     get definition of neb region
c     ----------------------------
      if (.not. rtdb_cget(irtdb,"qmmm:neb_region",
     +                    1,neb_region)) then
       neb_region="qmlink"
       if (.not. rtdb_cput(irtdb,"qmmm:neb_region",
     +                     1,neb_region)) 
     +  call errquit(pname//'Failed neb region',
     + 1, RTDB_ERR)
      end if
c
c     load restart file from perm dir if any
c     --------------------------------------
      if(qmmm_file_exist(rstfile,.true.)) then
        if(ga_nodeid().eq.0) then
          write(*,*) "reloading restart file",rstfile
        end if
        call mm_reload_rst_default(irtdb,rstfile)
      end if
      call ga_sync()
c
c     propagate coord of neb region
c     -----------------------------
      call qmmm_geom_push_active(irtdb)
c     
c     do inital energy in case esp charges are needed
c     -----------------------------------------------
      call md_sp()
      status = qmmm_energy_gradient(irtdb,.false.)
c
c     now relax other regions as requested
c     ------------------------------------
      status = qmmm_mixed_optimization(irtdb)
c     
c     do final energy and gradient
c     ----------------------------
      call qmmm_cons_free_exclusive(neb_region)
      call md_sp()
      status = qmmm_energy_gradient(irtdb,.true.)
      call qmmm_geom_get_grad(irtdb)
      call qmmm_energy_rtdb_push(irtdb)
c
c     finally save restart file to perm
c     ---------------------------------
      call util_file_name_resolve(rstfile,.false.)
      if(ga_nodeid().eq.0) then
          write(*,*) "saving restart file",rstfile
      end if
      call mm_write_restart_named(rstfile,.false.)

      qmmm_neb_gradient = status

 13   format(a,'_neb',i3.3,'.rst')
      return
      end

      subroutine qmmm_bq_gradient(rtdb,esp)
      implicit none
c
#include "rtdb.fh"
#include "util.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "nwc_const.fh"
#include "geom.fh"
#include "qmmm_params.fh"
#include "qmmm.fh"
c
      integer rtdb
      logical esp
c 
      integer i_c,h_c
      integer i_q,h_q
      integer i_g,h_g
  
      integer na
      character*30 pname
      character*32 theory
      logical pspw_bq_gradient
      external pspw_bq_gradient


      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit(pname//'theory not specified',555,
     &       INPUT_ERR)

      pname = "qmmm_bq_gradient"

      na  = qmmm_get_nbqa() 
     
      if(na.eq.0) return

      if(.not.ma_push_get(MT_DBL, na, 'active bq charges',
     &      h_q, i_q) ) call errquit(
     &      pname//'unable to allocate heap space',
     &      na, MA_ERR)
      call dfill(na,0.0d0,dbl_mb(i_q),1)

      if(.not.ma_push_get(MT_DBL,3*na, 'active bq coord',
     &      h_c, i_c) ) call errquit(
     &      pname//'unable to allocate heap space',
     &      3*na, MA_ERR)
      call dfill(3*na,0.0d0,dbl_mb(i_c),1)

      if(.not.ma_push_get(MT_DBL,3*na, 'active bq forces',
     &      h_g, i_g) ) call errquit(
     &      pname//'unable to allocate heap space',
     &      3*na, MA_ERR)
      call dfill(3*na,0.0d0,dbl_mb(i_g),1)


      call qmmm_get_cqbqa(na,
     >                      dbl_mb(i_c),
     >                      dbl_mb(i_q))


      if(esp) then
      
        call qmmm_bq_esp_grad(rtdb,
     >                      na,
     >                      dbl_mb(i_q),
     >                      dbl_mb(i_c),
     >                      dbl_mb(i_g))

      else if (theory.eq."pspw") then

         if(.not.pspw_bq_gradient(rtdb,
     >                      na,
     >                      dbl_mb(i_c),
     >                      dbl_mb(i_q),
     >                      dbl_mb(i_g)))
     >   call errquit(
     &      pname//'pspw_bq_gradient',0,0)

         call dscal(3*na,-1.0d0,dbl_mb(i_g),1)
 
      else 

        call qmmm_bq_forces(rtdb,
     >                      na,
     >                      dbl_mb(i_c),
     >                      dbl_mb(i_q),
     >                      dbl_mb(i_g))


      end if


      call qmmm_bq_push_grad(na,dbl_mb(i_g))

      if(.not.ma_pop_stack(h_g))
     &      call errquit(
     &      pname//'unable to pop stack',0,MA_ERR)

      if(.not.ma_pop_stack(h_c))
     &      call errquit(
     &      pname//'unable to pop stack',0,MA_ERR)

      if(.not.ma_pop_stack(h_q))
     &      call errquit(
     &      pname//'unable to pop stack',0,MA_ERR)

      
      end

      subroutine qmmm_bq_energy(e,rtdb,esp)
      implicit none
c
#include "rtdb.fh"
#include "util.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "nwc_const.fh"
#include "geom.fh"
#include "qmmm_params.fh"
#include "qmmm.fh"
#include "bq.fh"
c
      integer rtdb
      logical esp
      double precision e
c 
      integer i_c,h_c
      integer i_q,h_q
      integer na 
      integer bq_handle
      logical ext_bq_found
      character*30 pname


      pname = "qmmm_bq_energy"

      if(qmmm_print_debug())
     >   write(*,*) "in"//pname

      na  = qmmm_get_nbqa() 
     
      if(na.eq.0) return

      if(.not.ma_push_get(MT_DBL, na, 'active bq charges',
     &      h_q, i_q) ) call errquit(
     &      pname//'unable to allocate heap space',
     &      na, MA_ERR)
      call dfill(na,0.0d0,dbl_mb(i_q),1)

      if(.not.ma_push_get(MT_DBL,3*na, 'active bq coord',
     &      h_c, i_c) ) call errquit(
     &      pname//'unable to allocate heap space',
     &      3*na, MA_ERR)
      call dfill(3*na,0.0d0,dbl_mb(i_c),1)

      call qmmm_get_cqbqa(na,
     >                      dbl_mb(i_c),
     >                      dbl_mb(i_q))



      if(esp) then 

        call qmmm_bq_esp_energy(rtdb,
     >                      na,
     >                      dbl_mb(i_q),
     >                      dbl_mb(i_c),
     >                      e)

      else
        if(bq_get_active(bq_handle)) then
          ext_bq_found = .true.
          if(.not.bq_deactivate(bq_handle)) call errquit(
     &        pname//'unable to deactivate bqs',0,0)

        else
         ext_bq_found = .false.
        end if

        call qmmm_bq_energy_gen(rtdb,
     >                      na,
     >                      dbl_mb(i_c),
     >                      dbl_mb(i_q),
     >                      e)

        if(ext_bq_found) then
          if(.not.bq_activate(bq_handle)) call errquit(
     &        pname//'unable to activate bqs',0,0)
        end if

      end if

      if(.not.ma_pop_stack(h_c))
     &      call errquit(
     &      pname//'unable to pop stack',0,MA_ERR)

      if(.not.ma_pop_stack(h_q))
     &      call errquit(
     &      pname//'unable to pop stack',0,MA_ERR)

      if(qmmm_print_debug())
     >   write(*,*) "out"//pname

      
      end

      subroutine qmmm_map_fixed_atoms(rtdb,nact1)
      implicit none
c
c this will map active atom list from auxiliary geometry in qm modules
c to a  geometry used in task gradient. 
c
#include "rtdb.fh"
#include "util.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "nwc_const.fh"
#include "geom.fh"
#include "qmmm_params.fh"
#include "qmmm.fh"
c
      integer rtdb
      integer nact1
c
      character*32 pname
      integer nact
      integer i,j
      integer i_cons,h_cons
      integer i_act,h_act
      integer ma_type
      integer ng,nfg
      integer i_ig,i_ifg
      integer jj0,jj,ii,ii0
      integer nlink
      integer i_lb
      integer nal,h_al,i_al
      integer h_am,i_am
      logical aflag
     
      pname = "qmmm_map_fixed_atoms"

      if(qmmm_print_debug())
     >   write(*,*) "in"//pname
c
c     pointers to global index for auxiliary (i_ig)
c     and gradient geometry (i_ifg)
c     --------------------------------------------
      ng   = qmmm_geom_get_ngs()
      i_ig = qmmm_geom_get_i_ig()

      nfg   = qmmm_geom_get_nfg()
      i_ifg = qmmm_geom_get_i_ifg()

      nlink  = qmmm_get_nlink()
      i_lb   = qmmm_get_i_lb()
     
c
c     get active list
c     --------------
      if (rtdb_ma_get(rtdb, 'geometry:actlist', ma_type,
     $        nact, h_cons)) then

        if (.not. ma_get_index(h_cons, i_cons))
     $           call errquit(pname,h_cons,
     &       MA_ERR)

        if (.not.rtdb_delete(rtdb, 'geometry:actlist'))
     $       call errquit(pname,0,
     &       RTDB_ERR)

        if (.not.rtdb_put(rtdb, 'qmmm:actlist',
     >       mt_int,nact,int_mb(i_cons)))
     $       call errquit(pname,0,
     &       RTDB_ERR)

      else
        nact = ng
        if(.not.ma_alloc_get(mt_int, nact, 'qmmm actlist',
     &      h_cons, i_cons) ) call errquit(
     &      'qmmm_data_alloc: unable to allocate heap space',
     &      nact, MA_ERR)
        do i=1,nact
           int_mb(i_cons+i-1) = i
        end do
      end if

      if(qmmm_print_debug())
     >   write(*,*) " 1 "//pname

      if(qmmm_print_debug())
     >   write(*,*) "ng,nfg",ng,nfg

      if(qmmm_print_debug())
     >   write(*,*) "i_ig",(int_mb(i_ig+i-1),i=1,ng)

      if(qmmm_print_debug())
     >   write(*,*) "i_ifg",(int_mb(i_ifg+i-1),i=1,nfg)

      if(qmmm_print_debug())
     >   write(*,*) "i_cons",(int_mb(i_cons+i-1),i=1,nact)

c
c     create new active list
c     ----------------------
      if(.not.ma_push_get(MT_INT, nfg, 'qmmm fixed atom list',
     &      h_act, i_act) ) call errquit(
     &      pname//' unable to allocate stack',
     &      nact, MA_ERR)

      if(.not.ma_push_get(MT_LOG, ng, 'tmp qmmm act atom list',
     &      h_am, i_am) ) call errquit(
     &      pname//' unable to allocate stack',
     &      ng, MA_ERR)

c
c     get fixed atoms as specified by mm
c     if using aux_geom we assume that all the atoms there
c     are active by construction
c     ----------------------------------------------------
c      if(qmmm_aux_geom()) then
c        do i=1,ng
c         log_mb(i_am+i-1)=.true.
c        end do
c      else
        call qmmm_cons_get_map(ng,int_mb(i_ig),
     >                       log_mb(i_am))
c      end if


c     --------------------------------------------
c     find total number (nal) of link atoms
c     that are bonded to active qm atoms
c     and store their global index in i_al
c     Note that link atoms are always indexed
c     by the global index of the classical atom 
c     --------------------------------------------
      if(nlink.ne.0) then

        if(.not.ma_push_get(MT_INT, nlink, 'qmmm tmp link',
     &        h_al, i_al) ) call errquit(
     &        pname//' unable to allocate stack',
     &        nact, MA_ERR)

        nal=0
        do i=1,nact
c        get global index next active atom in the auxilary geometry
         ii0 = int_mb(i_cons+i-1)
         ii  = int_mb(i_ig+ii0-1)
         aflag = log_mb(i_am+ii0-1)
         if(aflag) then
           do j=1,nlink
c            get index of qm atom bonded to a link atom
             jj  = int_mb(i_lb+j-1)
c            if qm atom bonded to a link atom is active
c            store corresponding index of the link atom
             if(aflag.and.(jj.eq.ii)) then
               nal=nal+1
               int_mb(i_al+nal-1)=int_mb(i_lb+nlink+j-1)
             end if
           end do 
         end if
        end do
       else
        nal = 0
       end if

      if(qmmm_print_debug())
     >   write(*,*) "i_al",nal,(int_mb(i_al+i-1),i=1,nal)
c     ------------------------------------------------

c     ------------------------------------------------
c     construct active atom index
c     note that if classical boundary atom is
c     active the corresponding link atom is automatically
c     flagged as active in the first do loop because it 
c     carries the same global index. The second do loop takes
c     care of the case when qm boundary atom is active using
c     link index constructed above
c     ------------------------------------------------
      nact1=0
      do i=1,nfg
        ii = int_mb(i_ifg+i-1)
        do j=1,nact
         jj0 = int_mb(i_cons+j-1)
         jj  = int_mb(i_ig+jj0-1)
         aflag = log_mb(i_am+jj0-1)
         if(aflag.and.(ii.eq.jj)) then
            nact1 = nact1 + 1
            int_mb(i_act+nact1-1) = i
            goto 1
          end if
        end do
        do j=1,nal
         jj  = int_mb(i_al+j-1)
         if(ii.eq.jj) then
            nact1 = nact1 + 1
            int_mb(i_act+nact1-1) = i
            goto 1
          end if
        end do
1       continue
      end do

      if(qmmm_print_debug())
     >   write(*,*) " 2 "//pname

      if(nact1.ne.0) then
        if (.not.rtdb_put(rtdb, 'geometry:actlist',
     >          mt_int,nact1,int_mb(i_act)))
     >         call errquit(pname,0,
     >         RTDB_ERR)
      end if

      if (.not.ma_free_heap(h_cons))
     $   call errquit(pname,h_cons,
     &       MA_ERR)


      if (.not. ma_chop_stack(h_act) ) call errquit(
     &    pname//' ma_pop_stack ',
     &    0, MA_ERR)

      if(qmmm_print_debug())
     >   write(*,*) "out"//pname

      end

      subroutine qmmm_restore_fixed_atoms(rtdb)
      implicit none
c
#include "rtdb.fh"
#include "util.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "nwc_const.fh"
#include "geom.fh"
#include "qmmm_params.fh"
#include "qmmm.fh"
c
      integer rtdb
c
      character*32 pname
      integer nact
      integer i_cons,h_cons
      integer ma_type
      logical ignore
     
      pname = "qmmm_restore_fixed_atoms"
     
      ignore = rtdb_delete(rtdb, 'geometry:actlist')
     
      if (rtdb_ma_get(rtdb, 'qmmm:actlist', ma_type,
     $        nact, h_cons)) then
            if (.not. ma_get_index(h_cons, i_cons))
     $           call errquit(pname,h_cons,
     &       MA_ERR)
      else
        return
      end if

      if (.not.rtdb_delete(rtdb, 'qmmm:actlist'))
     $       call errquit(pname,0,
     &       RTDB_ERR)

      if (.not.rtdb_put(rtdb, 'geometry:actlist',
     >       mt_int,nact,int_mb(i_cons)))
     $       call errquit(pname,0,
     &       RTDB_ERR)

      if (.not.ma_free_heap(h_cons))
     $   call errquit(pname,h_cons,
     &       MA_ERR)

      end
 
