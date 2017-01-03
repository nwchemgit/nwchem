c $Id: bq_data.F 23019 2012-10-30 00:59:12Z d3y133 $

      function  mmi_init(rtdb)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "mm_data.fh"
      logical mmi_init
      integer rtdb
c
      logical ignore
      logical scnb_adj
      character*30 operation
      character*255 crdparmfile
      character*255 field
      character*84 tag
      integer fn
      character*30 pname

      pname = "mm_init"
c     write(*,*) pname

      ignore = MA_set_hard_fail(.true.)
      ignore=MA_set_auto_verify(.true.)
 
      ignore = rtdb_cget(rtdb,'task:operation',1,operation)

      if(operation.eq."neb") then
        call mm_neb_init(rtdb)
      else
        call mm_coords_init(rtdb)
        call mm_bndparm_init(rtdb)
        call mm_bond_coords_init()  
        call mm_vdw_init(rtdb)
        call mm_vdw_coords_init()
c       call mm_bond_call_test()  
        call mm_links_init(rtdb) 

        scnb_adj = .true.
        if(.not.rtdb_get(rtdb,'mm:scnb_adjust',mt_log,1,scnb_adj))
     >    scnb_adj = .false.
c       adjust scnb if true
        if(scnb_adj) call mm_vdw_adj_scnb()

        call mm_geom_init(rtdb,"geometry")      
        call mm_bq_init(rtdb)
      end if

      mmi_init = .true.
      return
911   mmi_init = .false.
      return
      end

      function  mmi_end(rtdb)
      implicit none
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "mm_link_data.fh"
#include "mm_geom_data.fh"
      integer rtdb

      logical mmi_end

      logical task_energy
      external task_energy
c
      
      character*30 pname
      character*30 operation
      logical ignore
      logical push_geom
      logical oprint
      logical rtdb_mode
      integer i, k
      integer master
      character*16 tag
      character*50 filename
      character*255 full_filename

      pname = "mm_end"
c     write(*,*) pname

      master=0

c     TP: if operation = "optimize" run final single point energy
c     and  aux_geom = .false.
      ignore = rtdb_cget(rtdb,'task:operation',1,operation)
      push_geom = .false.
      if(operation.eq."optimize") push_geom = .true.
      if(operation.eq."saddle") push_geom = .true.
      
      if(push_geom) then
        aux_geom = .false.
        call mm_geom_push_active(rtdb)
        call mm_geom_create_full(rtdb)
        call mm_link_update_bq_coords(rtdb)
      end if

      rtdb_mode = rtdb_parallel(.false.)

c     Let master node handle for writing necessary outputs
      if(ga_nodeid().eq.master) then
        oprint = .true.
        if(.not.rtdb_get(rtdb,"qmmm:active_region_xyz",mt_log,1,oprint))
     >   oprint = .false.

     
        if(oprint) then
          call util_file_prefix('act.xyz',filename)
          call util_file_name_noprefix(filename,.false.,
     >                                 .false.,
     >                                  full_filename)
          call mm_write_active_region_xyz(24, full_filename)
        end if

c       Write out final crdparms file. 
c       If opeation=neb crdparms file for each bead is written 
c       in mm_add_egrad.
        if(operation.eq."optimize".or.operation.eq."saddle") then
          call mm_create_crdparm(rtdb) 
        end if
      end if !! if(ga_nodeid().eq.master)

      call ga_sync()

      rtdb_mode = rtdb_parallel(rtdb_mode)

      ignore = rtdb_delete(rtdb,'bq') 

c     Deallocation!!
      call mm_coords_deallocate()
      call mm_vdw_deallocate()
      call mm_vdw_coords_deallocate()
      call mm_bonded_deallocate()  
      call mm_bond_coords_end()
      call mm_links_end()
      call mm_geom_end()
      call mm_bq_end()      

      mmi_end = .true.
      return
911   mmi_end = .false.
      return
      end

      subroutine mm_test(n,t,c)
      implicit none

      integer n
      character*(16) t(n)
      double precision c(3,n)
      integer i

      do i=1,n
         write(6,*) t(i),c(1,i),c(2,i),c(3,i)
      end do
      
      end

      subroutine mm_open_file(filename,fn)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
      character*(*) filename
      integer       fn
c
      character*180 buffer
      character*180 message
      character*30 pname,atag
c
      logical util_io_unit
      external util_io_unit
c
      pname = "mm_open_file"
c
      if(.not.util_io_unit(80,90,fn))
     +  call errquit(pname//"cannot get io unit",0,0)
c     first try to open file in the run directory
      buffer = filename
      message = "opening file "//buffer
      open(unit=fn,file=buffer,status='old',form="formatted",ERR=10)
      goto 800
10    continue
c     now try perm directory
      call util_file_name_resolve(buffer, .false.)
      message = "opening file "//buffer
      open(unit=fn,file=buffer,status='old',form="formatted",ERR=911)
800   continue
c     write(luout,*) "Successfully "//trim(message)
c     write(luout,*)
      return
911   call errquit("error "//trim(message),0,
     >        -1)
      end

      subroutine mm_add_energy(rtdb,e)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "global.fh"
#include "rtdb.fh"
#include "mm_bq_data.fh"
#include "mm_geom_data.fh"
#include "bq.fh"

      integer rtdb
      double precision e

      double precision bq_el_energy
      double precision bq_nuc_energy0
      double precision eqm0
      double precision eqmtot
      double precision emm
      double precision ebq0
c     dummy variables to pass as an argument in mm_cons_reaction
      double precision g(3,nact)  
c
      logical rtdb_mode
      integer master
      character*30 pname
c
      pname = "mm_add_energy"
c     write(*,*) pname
      master=0

      if(.not.bq_nuc_energy(rtdb,ebq0))
     >  call errquit('failed bq energy',0,RTDB_ERR)

      call mm_geom_restore(rtdb)
      call mm_cons_reaction(rtdb,.false.,e,g)

      eqmtot = e
      bq_nuc_energy0 = ebq0

c     bq elec interaction energy
      if (.not. rtdb_get(rtdb,'dft:bq_energy',mt_dbl,1,bq_el_energy))
     $     bq_el_energy = 0.0d0

      eqm0 = eqmtot - bq_el_energy - bq_nuc_energy0

      rtdb_mode = rtdb_parallel(.false.)

      if(ga_nodeid().eq.master) then
        call mm_links_bq_add_energy(rtdb,e)
c       call mm_links_bq_scaled_add_energy(rtdb,e)
        call mm_bonded_add_energy(rtdb,e)
        call mm_vdw_add_energy(rtdb,e)
      end if

      if(ga_nodeid().eq.master) then

        emm = e - eqmtot

        if(.not. rtdb_put(rtdb,'qmmm:qm_en',mt_dbl,1,eqmtot))
     >  call errquit('qmmm: failed put qm energy',0,RTDB_ERR)

        if(.not. rtdb_put(rtdb,'qmmm:qm_int_en',mt_dbl,1,eqm0))
     >  call errquit('qmmm: failed put qm internal energy',0,RTDB_ERR)

        if(.not. rtdb_put(rtdb,'qmmm:bq_nuc',mt_dbl,1,bq_nuc_energy0))
     >  call errquit('qmmm: failed put bq_nuc energy',0,RTDB_ERR)

        if(.not. rtdb_put(rtdb,'qmmm:bq_el',mt_dbl,1,bq_el_energy))
     >  call errquit('qmmm: failed put bq_el energy',0,RTDB_ERR)

        if(.not. rtdb_put(rtdb,'qmmm:mm_en',mt_dbl,1,emm))
     >  call errquit('qmmm: failed put emm',0,RTDB_ERR)

        if(.not. rtdb_put(rtdb,'qmmm:tot_en',mt_dbl,1,e))
     >  call errquit('qmmm: failed put total energy',0,RTDB_ERR)

        call mm_print_energy(rtdb)
      end if

      call ga_sync()
      call ga_brdcst(msg_smd,e,
     >     ma_sizeof(mt_dbl,1,mt_byte),master)
      call ga_sync()
      rtdb_mode = rtdb_parallel(rtdb_mode)

      end

      subroutine mm_add_egrad(rtdb,e,n,g)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "global.fh"
#include "rtdb.fh"
#include "mm_coords_data.fh"
#include "mm_link_data.fh"
#include "mm_bq_data.fh"
#include "bq.fh"

      integer rtdb
      double precision e

      double precision bq_el_energy
      double precision bq_nuc_energy0
      double precision eqm0
      double precision eqmtot
      double precision emm
      double precision ebq0

      integer n
      double precision g(3,n)
c
      integer i,j
      character*30 pname
      integer master
      logical rtdb_mode
      integer l_act, k_act
      integer nactive
      character*30 operation
      logical ignore 

      pname = "mm_add_egrad"
c     write(*,*) pname

      ignore = rtdb_cget(rtdb,'task:operation',1,operation)

      master=0

      if(ga_nodeid().eq.master) then
        call mm_links_adjust_forces(n,int_mb(i_iqml),g)
      end if
      call ga_sync()

      call mm_link_ebq_add_grad(rtdb,n,g)

      if(.not.bq_nuc_energy(rtdb,ebq0))
     >  call errquit('failed bq energy',0,RTDB_ERR)

      bq_nuc_energy0 = ebq0

      call mm_geom_restore(rtdb)
      call mm_cons_reaction(rtdb,.true.,e,g)

      eqmtot = e

c     bq elec interaction energy
      if (.not. rtdb_get(rtdb,'dft:bq_energy',mt_dbl,1,bq_el_energy))
     $     bq_el_energy = 0.0d0

      eqm0 = eqmtot - bq_el_energy - bq_nuc_energy0

      rtdb_mode = rtdb_parallel(.false.)

      if(ga_nodeid().eq.master) then
C     lookup table and list of active atoms
      if (.not. ma_push_get(MT_LOG,n,'active atoms',l_act,k_act))
     $     call errquit('grad: could not allocate l_act',n, MA_ERR)
        call mm_links_bq_add_egrad(rtdb,e,n,g)
c       call mm_links_bq_scaled_add_egrad(rtdb,e,n,g)
        call mm_bonded_add_grad(rtdb,e,n,g)
        call mm_vdw_add_egrad(rtdb,e,n,g)

        call grad_active_atoms(rtdb, n, log_mb(k_act), nactive)

c       TP: callin zero_forces() for fixed atom 
c       defined in constraint module
        call zero_forces(g,log_mb(k_act),n)

        call mm_print_grad_tot(rtdb,n,g)

      if (.not.ma_pop_stack(l_act))
     $     call errquit('grad:ma free act',1, MA_ERR)
      end if


      if(ga_nodeid().eq.master) then

        emm = e - eqmtot

        if(.not. rtdb_put(rtdb,'qmmm:qm_en',mt_dbl,1,eqmtot))
     >  call errquit('qmmm: failed put qm energy',0,RTDB_ERR)

        if(.not. rtdb_put(rtdb,'qmmm:qm_int_en',mt_dbl,1,eqm0))
     >  call errquit('qmmm: failed put qm internal energy',0,RTDB_ERR)

        if(.not. rtdb_put(rtdb,'qmmm:bq_nuc',mt_dbl,1,bq_nuc_energy0))
     >  call errquit('qmmm: failed put bq_nuc energy',0,RTDB_ERR)

        if(.not. rtdb_put(rtdb,'qmmm:bq_el',mt_dbl,1,bq_el_energy))
     >  call errquit('qmmm: failed put bq_el energy',0,RTDB_ERR)

        if(.not. rtdb_put(rtdb,'qmmm:mm_en',mt_dbl,1,emm))
     >  call errquit('qmmm: failed put emm',0,RTDB_ERR)

        if(.not. rtdb_put(rtdb,'qmmm:tot_en',mt_dbl,1,e))
     >  call errquit('qmmm: failed put total energy',0,RTDB_ERR)

        call mm_print_energy(rtdb)
        call mm_create_xyz_file(rtdb)
C       TP: create intermediate crdparms file for each bead for next iter
        if(operation.eq."neb") call mm_create_crdparm(rtdb)
      end if
           
      call ga_sync()
      call ga_brdcst(msg_smd,e,
     >     ma_sizeof(mt_dbl,1,mt_byte),master)
      call ga_brdcst(msg_smd,g,
     >     3*n*ma_sizeof(mt_dbl,1,mt_byte),master)
      call ga_sync()

      rtdb_mode = rtdb_parallel(rtdb_mode)

c     call mm_create_xyz_file(rtdb)

      return
911   call errquit("error "//trim(pname),0,-1)
      end

      subroutine mm_bonded_add_energy(rtdb,e)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "global.fh"
#include "rtdb.fh"
      integer rtdb
      double precision e
      character*30 pname

      pname = "mm_bonded_add_energy"
c     write(*,*) pname
  
      call mm_bond_add_energy(rtdb,e)
      call mm_angl_add_energy(rtdb,e)
      call mm_dihe_add_energy(rtdb,e)

      end


      subroutine mm_bonded_add_grad(rtdb,e,n,g)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "global.fh"
#include "rtdb.fh"
      integer rtdb
      double precision e
      integer n
      double precision g(n)
      character*30 pname

      pname = "mm_bonded_add_grad"
c     write(*,*) pname
  
      call mm_bond_add_egrad(rtdb,e,n,g)
      call mm_angl_add_egrad(rtdb,e,n,g)
      call mm_dihe_add_egrad(rtdb,e,n,g)

      end

      subroutine mm_print_grad_tot(rtdb,n,g)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "mm_geom_data.fh"
#include "mm_link_data.fh"
#include "mm_coords_data.fh"
      integer rtdb
      integer n
      double precision g(3,n)

      logical status
      logical ignore
      integer geom
      integer i_c,h_c
      integer i_tag,h_tag
      integer i, j, ncent
      character*30 pname
      character*30 message
      character*16 tag
      double precision scale
      logical geom_cart_get1
      external geom_cart_get1

      pname="mm_print_grad_tot"
c     write(*,*) pname

      if(.not.ma_push_get(mt_dbl,3*nact,'c',h_c,i_c))
     & call errquit('mm: Failed to allocate memory for c',
     & 3*nact, MA_ERR)

      call dcopy(nact*3,dbl_mb(i_rqml),1,dbl_mb(i_c),1)
      call util_convert_units("angstrom","au",scale)
      call dscal(3*nact, scale,dbl_mb(i_c),1)

      write(6,1000) "QM + MM",
     $     'x','y','z','x','y','z'
      do i = 1, nact
         call mm_coords_tag_set(byte_mb(i_tqml+16*(i-1)),tag)
         write(6,2000) i, tag,
     $        (dbl_mb(i_c+3*(i-1)+j),j=0,2),
     $         g(1,i),g(2,i),g(3,i)
      enddo
      write(6,*)

 1000 format(/,/,25X,A,' ENERGY GRADIENTS',/,/,4X,'atom',15X,
     $     'coordinates',
     $     24X,'gradient',/,6X,2(1X,(3(10X,A1))))
 2000 format(1X,I3,1X,A4,2(1X,3(1X,F10.6)))
      write(6,*)
      call util_flush(6)

      message = "memory deallocation"
      if(.not.ma_pop_stack(h_c)) goto 911

      return

911   call errquit("error "//trim(message),0,-1)

      end 

      subroutine mm_print_energy(rtdb)
      implicit none
#include "util.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "rtdb.fh"
#include "msgids.fh"
#include "global.fh"
#include "stdio.fh"

      integer rtdb
      
      character*32 pname

      double precision bq_el_energy
      double precision bq_nuc_energy
      double precision bq_energy
      double precision eqm0
      double precision eqmtot
      double precision emm
      double precision etot

      double precision scale_energy

      pname = "mm_print_energy"    
c     write(*,*) pname
 
      if(ga_nodeid().ne.0) return

      if(.not. rtdb_get(rtdb,'qmmm:qm_en',mt_dbl,1,eqmtot))
     >  call errquit('qmmm: failed get qm energy',0,RTDB_ERR)

      if(.not. rtdb_get(rtdb,'qmmm:qm_int_en',mt_dbl,1,eqm0))
     >  call errquit('qmmm: failed get qm internal energy',0,RTDB_ERR)

      if(.not. rtdb_get(rtdb,'qmmm:bq_nuc',mt_dbl,1,bq_nuc_energy))
     >  call errquit('qmmm: failed get bq_nuc energy',0,RTDB_ERR)

      if(.not. rtdb_get(rtdb,'qmmm:bq_el',mt_dbl,1,bq_el_energy))
     >  call errquit('qmmm: failed get bq_el energy',0,RTDB_ERR)

      if(.not. rtdb_get(rtdb,'qmmm:mm_en',mt_dbl,1,emm))
     >  call errquit('qmmm: failed get emm',0,RTDB_ERR)

      if(.not. rtdb_get(rtdb,'qmmm:tot_en',mt_dbl,1,etot))
     >  call errquit('qmmm: failed get total energy',0,RTDB_ERR)

      bq_energy = bq_el_energy + bq_nuc_energy

      call util_convert_units("au","kcal",scale_energy)

      write(*,19) 
      write(*,21) "QM/MM Energy"
      write(*,19) 
      write(*,21) "quantum energy", eqmtot, eqmtot*scale_energy
c     write(*,21) "quantum energy adjusted", qm_energy-eatoms,
c    >            (qm_energy-eatoms)*cau2kj
      if(bq_energy.ne.0.0d0) then
        write(*,21) "quantum energy internal", eqm0,
     >               eqm0*scale_energy
        write(*,21) "Bq-nuclear energy", bq_nuc_energy, 
     >               bq_nuc_energy*scale_energy
        write(*,21) "Bq-electron energy", bq_el_energy, 
     >               bq_el_energy*scale_energy
      end if
      write(*,21) "classical energy", emm,
     >            emm*scale_energy
      write(*,21) "total qmmm energy", etot, etot*scale_energy
      write(*,19) 

      write(*,*) 
19    FORMAT("@",72("-"))
20    FORMAT("@",1X,A,T34,F18.9)
21    FORMAT("@",1X,A,T34,F18.9," (",E12.6,2X, "kcal/mol)")

      end

      subroutine mm_create_xyz_file(rtdb)
      implicit none
#include "util.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "rtdb.fh"
#include "msgids.fh"
#include "global.fh"
#include "stdio.fh"
      integer rtdb     
 
      character*50 filename
      character*255 full_filename
      character*30 operation
      logical ignore
      character*32 pname
      character*20 tag
      integer ibead
      logical wrt_neb_xyz
      logical wrt_opt_xyz

      pname = "mm_create_xyz_file"
c     write(*,*) pname

      ignore = rtdb_cget(rtdb,'task:operation',1,operation)

      wrt_neb_xyz = .false.
      wrt_opt_xyz = .false.

      if(operation.eq."neb") wrt_neb_xyz = .true.
      if(operation.eq."optimize") wrt_opt_xyz = .true.
      if(operation.eq."saddle") wrt_opt_xyz = .true.

      if(wrt_neb_xyz) then
        if(.not.rtdb_get(rtdb,'neb:ibead',mt_int,1,ibead))
     >    call errquit('neb:ibead get',0,-1)
        tag = ' '
        write(tag,10) ibead
        call util_file_prefix(tag,filename)
        call util_file_name_noprefix(filename,.false.,
     >                               .false.,
     >                               full_filename)

        open(unit=19,file=full_filename,form='formatted')
        call mm_neb_create_xyz_file(19,ibead)
        close(19)

      else if(wrt_opt_xyz) then
        call mm_geom_create_xyz_file(rtdb)
      else 
        return
      end if

10    format('neb_',i3.3,'.xyz')   

      end 

      subroutine mm_create_crdparm(rtdb) 
      implicit none
#include "util.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "rtdb.fh"
#include "msgids.fh"
#include "global.fh"
#include "stdio.fh"
      integer rtdb     
 
      character*50 filename
      character*255 full_filename
      character*30 operation
      logical ignore
      character*32 pname
      character*20 tag
      integer ibead
 
      pname = "mm_create_crdparm"
c     write(*,*) pname

      ignore = rtdb_cget(rtdb,'task:operation',1,operation)

      if(operation.eq."neb") then
        if(.not.rtdb_get(rtdb,'neb:ibead',mt_int,1,ibead))
     >    call errquit('neb:ibead get',0,-1)
        filename = ' '
        write(filename,10) ibead
        

      else if(operation.eq."optimize".or.operation.eq."saddle") then
        if(.not.rtdb_cget(rtdb,"mm:crdparms:load",1,
     &                   filename))
     &  call errquit("error "//"mm:crdparmfile",0,-1)

      else
        return        
      end if

      open(unit=19,file=filename,form='formatted')
      call mm_write_qmcoords(19)
      call mm_write_mmcoords(19)
      call mm_write_bond(19)
      call mm_write_angl(19)
      call mm_write_dihe(19)
      call mm_write_vdw(19)
      call mm_write_vdw14(19)
      close(19)

10    format('nwchem_',i3.3,'.crdparms')

      end

      subroutine mm_task_gradient(rtdb,theory,e,g,status)
#include "util.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "rtdb.fh"
#include "msgids.fh"
#include "global.fh"
#include "stdio.fh"

c     input variables
      integer rtdb
      character*32 theory
      double precision e
      double precision g(3,*)
      logical status

c     local variables
      logical ignore
      character*30 operation
   
      logical task_gradient_doit
      external task_gradient_doit

      ignore = rtdb_cget(rtdb,'task:operation',1,operation)
      if(operation.eq."neb")
     >   call mm_neb_update_ibead(rtdb)
      call mm_geom_push_active(rtdb)
      call mm_geom_create_full(rtdb)
      call mm_map_fixed_atoms(rtdb)
      call mm_link_update_bq_coords(rtdb)
      status = task_gradient_doit(rtdb,theory,e,g)
      call mm_restore_fixed_atoms(rtdb)

      end

      subroutine mm_task_energy(rtdb,theory,e,status)
#include "util.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "rtdb.fh"
#include "msgids.fh"
#include "global.fh"
#include "stdio.fh"

c     input variables
      integer rtdb
      character*32 theory
      double precision e
      logical status

c     local variables
      logical ignore
      character*30 operation

      logical task_energy_doit
      external task_energy_doit
  
      ignore = rtdb_cget(rtdb,'task:operation',1,operation)
      if(operation.eq."neb")
     >   call mm_neb_update_ibead(rtdb)
      call mm_geom_push_active(rtdb)
      call mm_geom_create_full(rtdb)
      call mm_link_update_bq_coords(rtdb)
      status = task_energy_doit(rtdb,theory,e)

      end
