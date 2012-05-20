      subroutine smd_init(operation,rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*(*) operation
      integer rtdb
c
      character*80 pname
      character*80 tag

      pname = "smd_init"
      tag = operation
      if(tag.eq."energy") then
        call smd_init_sp(rtdb)
      else if(tag.eq."dynamics") then
        call smd_init_dynamics(rtdb)
      else if(tag.eq."fragment") then
        call smd_init_fragment(rtdb)
      else
        call errquit(
     >       pname//'unknown operation'//tag,0, 0)
      end if


      return
      end

      subroutine smd_init_sp(rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      integer rtdb
c
      logical ostatus
      character*72 tag
      character*32 pname
      character*255 filename
      character*72 sp_system
      character*72 sp_atom
      character*72 sp_param
      character*72 sp_bond
      character*72 sp_shakelist
      character*72 sp_bondlist
      character*72 sp_lattice
      character*72 sp_type
      character*72 sp_pot
      character*72 sp_coord
      character*72 sp_charge
      character*72 sp_force
      character*72 sp_energy
      character*72 sp_exlist
      character*72 sp_vlist
      character*72 sp_mass
      character*72 sp_cutoff
      character*72 sp_temperature
      character*72 sp_rtdb
      character*72 sp_ewald
c
      logical olattice
      logical obond
      integer master
      logical odistrib
c
      odistrib = .false.
      master = 0
      pname = "smd_init_dynamics"
c
      sp_atom      = "atom"
      sp_param     = "parameters"
      sp_lattice   = "lattice"
      sp_type      = "type"
      sp_pot       = "potential"
      sp_coord     = "coordinates"
      sp_charge    = "charge"
      sp_force     = "force"
      sp_energy    = "energy"
      sp_exlist    = "excluded_list"
      sp_vlist     = "verlet_list"
      sp_shakelist = "shakelist"
      sp_bondlist  = "bondlist"
      sp_bond      = "bond"
      sp_mass  = "mass"
      sp_cutoff  = "cutoff"
      sp_temperature = "temperature"
      sp_rtdb = "rtdb"
      sp_ewald = "ewald"
c
c
c     CREATE SYSTEM
c     ------------
      tag = "atom"
      call smd_system_set_component(tag,sp_atom)

      tag = "parameters"
      call smd_system_set_component(tag,sp_param)

      tag = "type"
      call smd_system_set_component(tag,sp_type)

      tag = "lattice"
      call smd_system_set_component(tag,sp_lattice)

      tag = "potential"
      call smd_system_set_component(tag,sp_pot)

      tag = "coordinates"
      call smd_system_set_component(tag,sp_coord)

      tag = "force"
      call smd_system_set_component(tag,sp_force)

      tag = "charge"
      call smd_system_set_component(tag,sp_charge)

      tag = "excl_list"
      call smd_system_set_component(tag,sp_exlist)

      tag = "verlet_list"
      call smd_system_set_component(tag,sp_vlist)

      tag = "energy"
      call smd_system_set_component(tag,sp_energy)

      tag = "bond"
      call smd_system_set_component(tag,sp_bond)

      tag = "shakelist"
      call smd_system_set_component(tag,sp_shakelist)

      tag = "bondlist"
      call smd_system_set_component(tag,sp_bondlist)

      tag = "mass"
      call smd_system_set_component(tag,sp_mass)

      tag = "cutoff"
      call smd_system_set_component(tag,sp_cutoff)

      tag = "temperature"
      call smd_system_set_component(tag,sp_temperature)

      tag = "rtdb"
      call smd_system_set_component(tag,sp_rtdb)

      tag = "ewald"
      call smd_system_set_component(tag,sp_ewald)
c
c      write(*,*) "in "//pname
c     DEFINE ALL THE COMPONENTS
c     --------------------------------
c
      if(ga_nodeid().eq.master) then
      call smd_rtdb_init(.false.,rtdb)
      call smd_cutoff_init()
      call smd_temper_init()
      call smd_atom_init_system()
      call smd_param_init_system()
      call smd_lat_init_system()
      call smd_pot_init_system()
      call smd_type_init_system()
      call smd_energy_init_system()
      call smd_force_init_system()
      call smd_charge_init_system()
      call smd_coords_init_system()
      call smd_exlist_init_system()
      call smd_vlist_init_system()
      call smd_bond_init_system()
      call smd_shakelist_init_system()
      call smd_bondlist_init_system()
      call smd_mass_init_system()
c
c     finish system setup
c     -------------------
      call smd_coords_rebox()
      call smd_ewald_init()
      call smd_system_print(6)

      end if
      call ga_sync()
      call smd_data_brdcst_create(master)
      return
      end

      subroutine smd_init_fragment(rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      integer rtdb
c
      logical ostatus
      character*72 tag
      character*32 pname
      character*255 filename
      character*72 sp_system
      character*72 sp_atom
      character*72 sp_param
      character*72 sp_bond
      character*72 sp_shakelist
      character*72 sp_bondlist
      character*72 sp_lattice
      character*72 sp_type
      character*72 sp_pot
      character*72 sp_coord
      character*72 sp_charge
      character*72 sp_force
      character*72 sp_energy
      character*72 sp_exlist
      character*72 sp_vlist
      character*72 sp_mass
      character*72 sp_cutoff
      character*72 sp_temperature
      character*72 sp_rtdb
      character*72 sp_ewald
      character*72 sp_fragment
c
      logical olattice
      logical obond
      integer master
      logical odistrib
c
      odistrib = .true.
c
      master = 0
      pname = "smd_init_dynamics"
c
      sp_atom      = "atom"
      sp_param     = "parameters"
      sp_lattice   = "lattice"
      sp_type      = "type"
      sp_pot       = "potential"
      sp_coord     = "coordinates"
      sp_charge    = "charge"
      sp_fragment    = "fragment"
      sp_force     = "force"
      sp_energy    = "energy"
      sp_exlist    = "excluded_list"
      sp_vlist     = "verlet_list"
      sp_shakelist = "shakelist"
      sp_bondlist  = "bondlist"
      sp_bond      = "bond"
      sp_mass  = "mass"
      sp_cutoff  = "cutoff"
      sp_temperature = "temperature"
      sp_rtdb = "rtdb"
      sp_ewald = "ewald"
c
c
c     CREATE SYSTEM
c     ------------
      tag = "atom"
      call smd_system_set_component(tag,sp_atom)

      tag = "parameters"
      call smd_system_set_component(tag,sp_param)

      tag = "type"
      call smd_system_set_component(tag,sp_type)

      tag = "lattice"
      call smd_system_set_component(tag,sp_lattice)

      tag = "potential"
      call smd_system_set_component(tag,sp_pot)

      tag = "coordinates"
      call smd_system_set_component(tag,sp_coord)

      tag = "force"
      call smd_system_set_component(tag,sp_force)

      tag = "charge"
      call smd_system_set_component(tag,sp_charge)

      tag = "fragment"
      call smd_system_set_component(tag,sp_fragment)

      tag = "excl_list"
      call smd_system_set_component(tag,sp_exlist)

      tag = "verlet_list"
      call smd_system_set_component(tag,sp_vlist)

      tag = "energy"
      call smd_system_set_component(tag,sp_energy)

      tag = "bond"
      call smd_system_set_component(tag,sp_bond)

      tag = "shakelist"
      call smd_system_set_component(tag,sp_shakelist)

      tag = "bondlist"
      call smd_system_set_component(tag,sp_bondlist)

      tag = "mass"
      call smd_system_set_component(tag,sp_mass)

      tag = "cutoff"
      call smd_system_set_component(tag,sp_cutoff)

      tag = "temperature"
      call smd_system_set_component(tag,sp_temperature)

      tag = "rtdb"
      call smd_system_set_component(tag,sp_rtdb)

      tag = "ewald"
      call smd_system_set_component(tag,sp_ewald)
c
c      write(*,*) "in "//pname
c     DEFINE ALL THE COMPONENTS
c     --------------------------------
c
      if(ga_nodeid().eq.master) then
      call smd_rtdb_init(.false.,rtdb)
      call smd_cutoff_init()
      call smd_temper_init()
      call smd_atom_init_system()
      call smd_param_init_system()
      call smd_lat_init_system()
      call smd_pot_init_system()
      call smd_type_init_system()
      call smd_energy_init_system()
      call smd_force_init_system()
      call smd_charge_init_system()
      call smd_coords_init_system()
      call smd_exlist_init_system()
      call smd_vlist_init_system()
      call smd_bond_init_system()
      call smd_shakelist_init_system()
      call smd_bondlist_init_system()
      call smd_mass_init_system()
c
c     finish system setup
c     -------------------
      call smd_coords_rebox()
      call smd_ewald_init()
      call smd_system_print(6)

      end if
      call ga_sync()
      call smd_data_brdcst_create(master)
      call smd_fragment_init_system()
      return
      end

      subroutine smd_init_dynamics(rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      integer rtdb
c
      logical ostatus
      character*72 tag
      character*32 pname
      character*255 filename
      character*72 sp_system
      character*72 sp_atom
      character*72 sp_param
      character*72 sp_bond
      character*72 sp_shakelist
      character*72 sp_bondlist
      character*72 sp_lattice
      character*72 sp_type
      character*72 sp_pot
      character*72 sp_coord
      character*72 sp_charge
      character*72 sp_force
      character*72 sp_energy
      character*72 sp_exlist
      character*72 sp_vlist
      character*72 sp_velocity
      character*72 sp_mass
      character*72 sp_cutoff
      character*72 sp_temperature
      character*72 sp_rtdb
      character*72 sp_ewald
c
      logical olattice
      logical obond
      integer master
c
      logical odistrib
c
      odistrib = .false.
      master = 0
      pname = "smd_init_dynamics"
c
      sp_atom      = "atom"
      sp_param     = "parameters"
      sp_lattice   = "lattice"
      sp_type      = "type"
      sp_pot       = "potential"
      sp_coord     = "coordinates"
      sp_charge    = "charge"
      sp_force     = "force"
      sp_energy    = "energy"
      sp_exlist    = "excluded_list"
      sp_vlist     = "verlet_list"
      sp_shakelist = "shakelist"
      sp_bondlist  = "bondlist"
      sp_bond      = "bond"
      sp_velocity  = "velocity"
      sp_mass  = "mass"
      sp_cutoff  = "cutoff"
      sp_temperature = "temperature"
      sp_rtdb = "rtdb"
      sp_ewald = "ewald"
c
c
c     CREATE SYSTEM
c     ------------
      tag = "atom"
      call smd_system_set_component(tag,sp_atom)

      tag = "parameters"
      call smd_system_set_component(tag,sp_param)

      tag = "type"
      call smd_system_set_component(tag,sp_type)

      tag = "lattice"
      call smd_system_set_component(tag,sp_lattice)

      tag = "potential"
      call smd_system_set_component(tag,sp_pot)

      tag = "coordinates"
      call smd_system_set_component(tag,sp_coord)

      tag = "force"
      call smd_system_set_component(tag,sp_force)

      tag = "charge"
      call smd_system_set_component(tag,sp_charge)

      tag = "excl_list"
      call smd_system_set_component(tag,sp_exlist)

      tag = "verlet_list"
      call smd_system_set_component(tag,sp_vlist)

      tag = "energy"
      call smd_system_set_component(tag,sp_energy)

      tag = "bond"
      call smd_system_set_component(tag,sp_bond)

      tag = "shakelist"
      call smd_system_set_component(tag,sp_shakelist)

      tag = "bondlist"
      call smd_system_set_component(tag,sp_bondlist)

      tag = "velocity"
      call smd_system_set_component(tag,sp_velocity)

      tag = "mass"
      call smd_system_set_component(tag,sp_mass)

      tag = "cutoff"
      call smd_system_set_component(tag,sp_cutoff)

      tag = "temperature"
      call smd_system_set_component(tag,sp_temperature)

      tag = "rtdb"
      call smd_system_set_component(tag,sp_rtdb)

      tag = "ewald"
      call smd_system_set_component(tag,sp_ewald)
c
c      write(*,*) "in "//pname
c     DEFINE ALL THE COMPONENTS
c     --------------------------------
c
      if(ga_nodeid().eq.master) then
      call smd_rtdb_init(.false.,rtdb)
      call smd_cutoff_init()
      call smd_temper_init()
      call smd_atom_init_system()
      call smd_param_init_system()
      call smd_lat_init_system()
      call smd_pot_init_system()
      call smd_type_init_system()
      call smd_energy_init_system()
      call smd_force_init_system()
      call smd_coords_init_system()
      call smd_charge_init_system()
      call smd_exlist_init_system()
      call smd_vlist_init_system()
      call smd_bond_init_system()
      call smd_shakelist_init_system()
      call smd_bondlist_init_system()
      call smd_vel_init_system()
      call smd_mass_init_system()
c
c     finish system setup
c     -------------------
      call smd_coords_rebox()
      call smd_ewald_init()
      call smd_system_print(6)

      end if
      call ga_sync()
      call smd_data_brdcst_create(master)
      return
      end

c $Id$
