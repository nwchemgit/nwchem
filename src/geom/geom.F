      block data geom_data
      implicit none
#include "geomP.fh"
c     
      data ngeom_rtdb /0/
      data active /max_geom*.false./
      data symbols/
     $     'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne',
     $     'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', 'K ', 'Ca',
     $     'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
     $     'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y ', 'Zr',
     $     'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
     $     'Sb', 'Te', 'I ', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
     $     'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
     $     'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
     $     'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
     $     'Pa', 'U ', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',
     $     'Md', 'No', 'Lr'/
      data elements/
     $     'Hydrogen', 'Helium', 'Lithium', 'Beryllium', 'Boron',
     $     'Carbon', 'Nitrogen', 'Oxygen', 'Fluorine', 'Neon', 'Sodium',
     $     'Magnesium', 'Aluminium', 'Silicon', 'Phosphorous',
     $     'Sulphur', 'Chlorine', 'Argon', 'Potassium', 'Calcium',
     $     'Scandium', 'Titanium', 'Vanadium', 'Chromium', 'Manganese',
     $     'Iron', 'Cobalt', 'Nickel', 'Copper', 'Zinc', 'Gallium',
     $     'Germanium', 'Arsenic', 'Selenium', 'Bromine', 'Krypton',
     $     'Rubidium', 'Strontium', 'Yttrium', 'Zirconium', 'Niobium',
     $     'Molybdenum', 'Technetium', 'Ruthenium', 'Rhodium',
     $     'Palladium', 'Silver', 'Cadmium', 'Indium', 'Tin',
     $     'Antinomy', 'Tellurium', 'Iodine', 'Xenon', 'Caesium',
     $     'Barium', 'Lanthanum', 'Cerium', 'Praseodymium', 'Neodymium',
     $     'Promethium', 'Samarium', 'Europium', 'Gadolinium',
     $     'Terbium', 'Dysprosium', 'Holmium', 'Erbium', 'Thulium',
     $     'Ytterbium', 'Lutetium', 'Hafnium', 'Tantalum', 'Tungsten',
     $     'Rhenium', 'Osmium', 'Iridium', 'Platinum', 'Gold',
     $     'Mercury', 'Thallium', 'Lead', 'Bismuth', 'Polonium',
     $     'Astatine', 'Radon', 'Francium', 'Radium', 'Actinium',
     $     'Thorium', 'Protoactinium', 'Uranium', 'Neptunium',
     $     'Plutonium', 'Americium', 'Curium', 'Berkelium',
     $     'Californium', 'Einsteinium', 'Fermium', 'Mendelevium',
     $     'Nobelium', 'Lawrencium'/
c     
      end
      logical function geom_check_handle(geom, msg)
      implicit none
#include "geomP.fh"
c     
      integer geom              ! [input]
      character*(*) msg         ! [input]
c     
      geom_check_handle = geom.gt.0 .and. geom.le.max_geom
      if (geom_check_handle) geom_check_handle = geom_check_handle
     $     .and. active(geom)
c     
      if (.not. geom_check_handle) then
         write(6,*) msg,': geometry handle invalid ', geom
         call geom_err_info(msg)
      endif
c     
      end
      logical function geom_check_cent(geom, msg, icent)
      implicit none
#include "geomP.fh"
c     
      integer geom              ! [input]
      character*(*) msg         ! [input]
      integer icent             ! [input]
c     
      geom_check_cent = icent.gt.0 .and. icent.le.ncenter(geom)
      if (.not. geom_check_cent) then
         write(6,*) msg,': icent invalid ', icent,
     $        names(geom)(1:lenn(geom))
         call geom_err_info(msg)
         call geom_print(geom)
      endif
c     
      end
      logical function geom_rtdb_in(rtdb)
      implicit none
#include "geomP.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
c     
      integer rtdb              ! [input]
      integer geom
c     
c     load in info about known geometries ... this is more
c     for diagnostic and debugging purposes
c     
      geom_rtdb_in = .false.
      ngeom_rtdb = 0
      if (rtdb_get(rtdb, 'geometry:ngeom', mt_int, 1, ngeom_rtdb))
     $     then
         if (.not. rtdb_cget(rtdb, 'geometry:names', max_geom,
     $        names_rtdb)) then
            write(6,*) 'geom_rtdb_in: rtdb corrupt'
         else
            do geom = 1, ngeom_rtdb
               lenr(geom) = inp_strlen(names_rtdb(geom))
            enddo
            geom_rtdb_in = .true.
         endif
      endif
c     
      end
      logical function geom_rtdb_out(rtdb)
      implicit none
#include "geomP.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
c     
      integer rtdb              ! [input]
c     
c     output to rtdb info about known geometries
c     
      geom_rtdb_out  =
     $     rtdb_put(rtdb, 'geometry:ngeom', mt_int, 1, ngeom_rtdb)
     $     .and.
     $     rtdb_cput(rtdb, 'geometry:names', max_geom, names_rtdb)
      if (.not. geom_rtdb_out) 
     $     write(6,*) ' geom_rtdb_out: rtdb is corrupt '
c     
      end
      logical function geom_rtdb_add(rtdb, name)
      implicit none
#include "geomP.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
c     
      integer rtdb              ! [input]
      character*(*) name        ! [input]
      integer geom
      logical status
      integer ln
      logical geom_rtdb_in, geom_rtdb_out
      external geom_rtdb_in, geom_rtdb_out
c     
c     See if name is on the rtdb already
c     
      ln = inp_strlen(name)
      status = geom_rtdb_in(rtdb)
      geom_rtdb_add = .true.
      do geom = 1, ngeom_rtdb
         if (name(1:ln) .eq. names_rtdb(geom)(1:lenr(geom))) return
      enddo
c     
c     Name is not present ... add and rewrite info
c     
      if (ngeom_rtdb .eq. max_geom_rtdb) then
         write(6,*) ' geom_rtdb_add: too many geometries on rtdb ', name
         geom_rtdb_add = .false.
         return
      endif
      ngeom_rtdb = ngeom_rtdb + 1
      names_rtdb(ngeom_rtdb) = name
      lenr(ngeom_rtdb) = ln
c     
      if (.not. geom_rtdb_out(rtdb)) then
         write(6,*) ' geom_rtdb_add: rtdb error adding ', name(1:ln)
         geom_rtdb_add = .false.
         return
      endif
c     
      geom_rtdb_add = .true.
c     
      end
      subroutine geom_err_info(info)
      implicit none
#include "geomP.fh"
c     
      character*(*) info        ! [input]
      integer geom
      integer ngeom
c     
c     For internal use of the geom routines only: print out
c     info of known geometries to aid in diagnosing a problem
c     
      ngeom = 0
      do geom = 1, max_geom
         if (active(geom)) ngeom = ngeom + 1
      enddo
      write(6,1) info, ngeom
 1    format(' ',a,': open geometies: ',i2)
      ngeom = 0
      do geom = 1, max_geom
         if (active(geom)) then
            write(6,2) ngeom, info, names(geom)(1:lenn(geom)),
     $           trans(geom)(1:lent(geom))
 2          format(' ',i2,' ',a,': "',a, '" -> "', a,'"')
         endif
      enddo
      if (ngeom_rtdb .gt. 0) then
         write(6,3) info, ngeom_rtdb
 3       format(' ',a,': geometries in last accessed data base: ', i2)
         do geom = 1, ngeom_rtdb
            write(6,4) names_rtdb(geom)(1:lenr(geom))
 4          format(' ',a)
         enddo
      endif
c     
      end
      logical function geom_rtdb_load(rtdb, geom, name)
      implicit none
#include "rtdb.fh"
#include "geomP.fh"
#include "context.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"

c     
      integer rtdb              ! [input]
      integer geom              ! [input]
      character*(*) name        ! [input]
c     
      character*256 tmp
      integer k, nelem
      logical status
      logical geom_check_handle, geom_rtdb_in
      external geom_check_handle, geom_rtdb_in
c     
      geom_rtdb_load = geom_check_handle(geom, 'geom_rtdb_load')
      if (.not. geom_rtdb_load) return
      status = geom_rtdb_in(rtdb)
c     
c     translate the provided name
c     
      names(geom) = name
      lenn(geom) = inp_strlen(name)
      trans(geom) = 'junk'
      if (.not.context_rtdb_match(rtdb, name, trans(geom)))
     $     trans(geom) = name
      lent(geom) = inp_strlen(trans(geom))
c     
c     now get the info from the data base
c     
      tmp = 'geometry:'//trans(geom)(1:lent(geom))
      k = inp_strlen(tmp)+1
      status = .true.
c
      tmp(k:) = ':ncenter'
      status = status .and. 
     $     rtdb_get(rtdb, tmp, mt_int, 1, ncenter(geom))
      tmp(k:) = ':coords'
      status = status .and. 
     $     rtdb_get(rtdb, tmp, mt_dbl, max_cent*3, coords(1,1,geom))
      tmp(k:) = ':charges'
      status = status .and. 
     $     rtdb_get(rtdb, tmp, mt_dbl, max_cent, charge(1,geom))
      tmp(k:) = ':efield'
      status = status .and. 
     $     rtdb_get(rtdb, tmp, mt_dbl, 3, efield(1,geom))
      tmp(k:) = ':lattice vectors'
      status = status .and. 
     $     rtdb_get(rtdb, tmp, mt_dbl, 3, lattice_vectors(1,geom))
      tmp(k:) = ':lattice angles'
      status = status .and. 
     $     rtdb_get(rtdb, tmp, mt_dbl, 3, lattice_angles(1,geom))
      tmp(k:) = ':system type'
      status = status .and. 
     $     rtdb_get(rtdb, tmp, mt_int, 1, isystype(geom))
      tmp(k:) = ':no. unique centers'
      status = status .and. 
     $     rtdb_get(rtdb, tmp, mt_int, 1, ncenter_unique(geom))
      tmp(k:) = ':group number'
      status = status .and. 
     $     rtdb_get(rtdb, tmp, mt_int, 1, group_number(geom))
      tmp(k:) = ':setting number'
      status = status .and. 
     $     rtdb_get(rtdb, tmp, mt_int, 1, setting_number(geom))
      tmp(k:) = ':recip vectors'
      status = status .and. 
     $     rtdb_get(rtdb, tmp, mt_dbl, 3, recip_lat_vectors(1,geom))
      tmp(k:) = ':recip angles'
      status = status .and. 
     $     rtdb_get(rtdb, tmp, mt_dbl, 3, recip_lat_angles(1,geom))
      tmp(k:) = ':direct volume'
      status = status .and. 
     $     rtdb_get(rtdb, tmp, mt_dbl, 1, volume_direct(geom))
      tmp(k:) = ':unique centers'
      status = status .and. 
     $     rtdb_get(rtdb, tmp, mt_int, ncenter_unique(geom), 
     $     unique_cent(1,geom))
      tmp(k:) = ':tags'
      status = status .and. 
     $     rtdb_cget(rtdb, tmp, max_cent, tags(1,geom))
c
c--> get symmetry operators, number of operators and operator/atom
c    map from rtdb
c
      tmp(k:) = ':num_operators'
      status = status .and. 
     $     rtdb_get(rtdb, tmp, mt_int, 1, sym_num_ops(geom))
      tmp(k:) = ':operators'
      status = status .and. 
     $     rtdb_get(rtdb, tmp, mt_dbl, max_sym_ops*3*4,
     $     sym_ops(1,1,geom))
      tmp(k:) = ':map_atoms'
      status = status .and. 
     $     rtdb_ma_get(rtdb, tmp, mt_int, nelem,
     $     sym_center_map_handle(geom))
      if (nelem .ne. sym_num_ops(geom)*ncenter(geom)) call errquit
     $     ('geom_rtdb_load: invalid no. of element in sym map', nelem)
      if (.not. ma_get_index(sym_center_map_handle(geom),
     $     sym_center_map_index(geom)))call errquit
     $     ('geom_rtdb_load: bad ma handle', 0)
c     
      if (.not. status) then
         write(6,*) ' geom_rtdb_load: not found or rtdb corrupt: ',
     $        names(geom)(1:lenn(geom)), ' -> ',
     $        trans(geom)(1:lent(geom))
         call geom_err_info('geom_rtdb_load')
         geom_rtdb_load = .false.
         return
      endif
c     
c     Determine if external fields are applied
c     
      oefield(geom) = 
     $     ddot(3, efield(1,geom), 1, efield(1,geom), 1) .gt. 0.0d0
c
c     compute effective nuclear repulsion energy, dipole and
c     interaction with external fields
c
      call geom_compute_values(geom)
c
      active(geom) = .true.
      geom_rtdb_load = .true.
c     
      end
      subroutine geom_compute_values(geom)
      implicit none
#include "geomP.fh"
#include "util.fh"
      integer geom              ! [input]
c
c     compute effective nuclear repulsion energy, dipole and
c     interaction with external fields
c
      double precision e, r
      integer i, j
c
      e = 0.0d0
      ndipole(1,geom) = 0.0d0
      ndipole(2,geom) = 0.0d0
      ndipole(3,geom) = 0.0d0
c
c     compute nuclear dipole moment and usual nuclear repulsion energy
c
      do i = 1,ncenter(geom)
         do j = 1, 3
            ndipole(j,geom) = ndipole(j,geom) +
     $           charge(i,geom)*coords(j,i,geom)
         enddo
         do j = i+1, ncenter(geom)
            r = dsqrt(
     $           (coords(1,i,geom)-coords(1,j,geom))**2 + 
     $           (coords(2,i,geom)-coords(2,j,geom))**2 + 
     $           (coords(3,i,geom)-coords(3,j,geom))**2)
            e = e + charge(i,geom)*charge(j,geom)/r
         enddo
      enddo
c
c     add in interaction of nuclear dipole with external field
c
      e = e + ddot(3, ndipole(1,geom), 1, efield(1,geom), 1)
c
      erep(geom) = e
c
      end
      logical function geom_rtdb_store(rtdb, geom, name)
      implicit none
#include "geomP.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "context.fh"
#include "tcgmsg.fh"
#include "inp.fh"
c     
      integer rtdb              ! [input]
      character*(*) name        ! [input] ... if blank use current name
      integer geom              ! [input]
      logical geom_check_handle, geom_rtdb_add
      external geom_check_handle, geom_rtdb_add
      logical status
      character*256 tmp
      integer k
c     
      geom_rtdb_store =  geom_check_handle(geom, 'geom_rtdb_store')
      if (.not. geom_rtdb_store) return
      if (name .ne. ' ') then
         names(geom) = name
         lenn(geom) = inp_strlen(name)
      endif
c
c     If not process 0 return ... this is so that input routines
c     can be completely single threaded
c
      if (nodeid() .ne. 0) return
c     
c     try to translate the name
c     
      trans(geom) = 'junk'
      if (.not. context_rtdb_match(rtdb, name, trans(geom)))
     $     trans(geom) = name
      lent(geom) = inp_strlen(trans(geom))
c     
c     now put the info into the data base
c     
      tmp = 'geometry:'//trans(geom)(1:lent(geom))
      k = inp_strlen(tmp)+1
      status = .true.
c
      tmp(k:) = ':ncenter'
      status = status .and. 
     $     rtdb_put(rtdb, tmp, mt_int, 1, ncenter(geom))
      tmp(k:) = ':coords'
      status = status .and. 
     $     rtdb_put(rtdb, tmp, mt_dbl, ncenter(geom)*3,
     $     coords(1,1,geom))
      tmp(k:) = ':charges'
      status = status .and. 
     $     rtdb_put(rtdb, tmp, mt_dbl, ncenter(geom),
     $     charge(1,geom))
      tmp(k:) = ':efield'
      status = status .and. 
     $     rtdb_put(rtdb, tmp, mt_dbl, 3, efield(1,geom))
      tmp(k:) = ':lattice vectors'
      status = status .and. 
     $     rtdb_put(rtdb, tmp, mt_dbl, 3, lattice_vectors(1,geom))
      tmp(k:) = ':lattice angles'
      status = status .and. 
     $     rtdb_put(rtdb, tmp, mt_dbl, 3, lattice_angles(1,geom))
      tmp(k:) = ':system type'
      status = status .and. 
     $     rtdb_put(rtdb, tmp, mt_int, 1, isystype(geom))
      tmp(k:) = ':no. unique centers'
      status = status .and. 
     $     rtdb_put(rtdb, tmp, mt_int, 1, ncenter_unique(geom))
      tmp(k:) = ':group number'
      status = status .and. 
     $     rtdb_put(rtdb, tmp, mt_int, 1, group_number(geom))
      tmp(k:) = ':setting number'
      status = status .and. 
     $     rtdb_put(rtdb, tmp, mt_int, 1, setting_number(geom))
      tmp(k:) = ':recip vectors'
      status = status .and. 
     $     rtdb_put(rtdb, tmp, mt_dbl, 3, recip_lat_vectors(1,geom))
      tmp(k:) = ':recip angles'
      status = status .and. 
     $     rtdb_put(rtdb, tmp, mt_dbl, 3, recip_lat_angles(1,geom))
      tmp(k:) = ':direct volume'
      status = status .and. 
     $     rtdb_put(rtdb, tmp, mt_dbl, 1, volume_direct(geom))
      tmp(k:) = ':unique centers'
      status = status .and. 
     $     rtdb_put(rtdb, tmp, mt_int, ncenter_unique(geom), 
     $     unique_cent(1,geom))
      tmp(k:) = ':tags'
      status = status .and. 
     $     rtdb_cput(rtdb, tmp, ncenter(geom), tags(1,geom))
c
c--> put symmetry operators, number of operators and operator/atom
c    map out to rtdb
c
      tmp(k:) = ':num_operators'
      status = status .and. 
     $     rtdb_put(rtdb, tmp, mt_int, 1, sym_num_ops(geom))
      tmp(k:) = ':operators'
      status = status .and. 
     $     rtdb_put(rtdb, tmp, mt_dbl, max_sym_ops*3*4,
     $     sym_ops(1,1,geom))
      tmp(k:) = ':map_atoms'
      status = status .and. 
     $     rtdb_put(rtdb, tmp, mt_int, ncenter(geom)*sym_num_ops(geom),
     $     int_mb(sym_center_map_index(geom)))
c
c     insert translated name into list of known geometries
c
      status = status .and. geom_rtdb_add(rtdb, name)
c
c     check that all rtdb operations were successful
c
      if (.not. status) then
         write(6,*) ' geom_rtdb_store: write to rtdb failed',
     $        names(geom)(1:lenn(geom)), ' -> ',
     $        trans(geom)(1:lent(geom))
         call geom_err_info('geom_rtdb_store')
         geom_rtdb_store = .false.
         return
      endif
      geom_rtdb_store = .true.
c
      end
      logical function geom_rtdb_delete(rtdb, name)
      implicit none
#include "geomP.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
c
      integer rtdb              ! [input]
      character*(*) name        ! [input]
      character*256 translation, tmp
      integer lt, geom, geom2, k
      logical status
      logical geom_rtdb_in, geom_rtdb_out
      external geom_rtdb_in, geom_rtdb_out
c
c     try to translate the provided name
c
      if (.not. rtdb_cget(rtdb, name, 1, translation))
     $     translation = name
      lt = inp_strlen(translation)
c
c     locate name in list and remove
c     
      status = geom_rtdb_in(rtdb)
      do geom = 1, ngeom_rtdb
         if (names_rtdb(geom)(1:lenr(geom)) .eq. translation(1:lt))
     $        goto 10
      enddo
      write(6,*) ' geom_rtdb_delete: no such geometry ',
     $     name(1:inp_strlen(name)),' -> ', translation(1:lt)
 10   do geom2 = geom+1, ngeom_rtdb
         names_rtdb(geom2-1) = names_rtdb(geom2)
      enddo
      ngeom_rtdb = ngeom_rtdb - 1
      status = geom_rtdb_out(rtdb)
c
c     delete each entry assoicated with a geometry in the database
c
      tmp = 'geometry:'//translation(1:lt)
      k = inp_strlen(tmp)+1
c
      tmp(k:) = ':ncenter'
      status = status .and. rtdb_delete(rtdb, tmp)
      tmp(k:) = ':coords'
      status = status .and. rtdb_delete(rtdb, tmp)
      tmp(k:) = ':charges'
      status = status .and. rtdb_delete(rtdb, tmp)
      tmp(k:) = ':efield'
      status = status .and. rtdb_delete(rtdb, tmp)
      tmp(k:) = ':latvec'
      status = status .and. rtdb_delete(rtdb, tmp)
      tmp(k:) = ':tags'
      status = status .and. rtdb_delete(rtdb, tmp)
c
c     check status of all rtdb stores
c
      if (.not. status) then
         write(6,*) ' geom_rtdb_delete: rtdb corrupt ',translation(1:lt)
         call geom_err_info('geom_rtdb_delete')
         geom_rtdb_delete = .false.
         return
      endif
c
      geom_rtdb_delete = .true.
c
      end
      logical function geom_create(geom, name)
      implicit none
#include "geomP.fh"
#include "inp.fh"
c     
      integer geom              ! [output]
      character*(*) name        ! [input]
c
      integer i,j
c
c     Assign the next free slot for a geometry
c
      do geom = 1, max_geom
         if (.not. active(geom)) goto 10
      enddo
      write(6,1) name
 1    format(' geom_create: too many geoms trying to create ', a)
      call geom_err_info('geom_create')
      geom_create = .false.
      return
 10   continue
c
c     store info about the geometry
c
      names(geom) = name
      trans(geom) = ' '
      lenn(geom) = inp_strlen(name)
      ncenter(geom) = 0
      active(geom) = .true.
      geom_create = .true.
      oefield(geom) = .false.
      operiodic(geom) = .false.
      ncenter_unique(geom) = 0
      isystype(geom) = 0
      group_number(geom) = 1
      setting_number(geom) = 0
      sym_center_map_handle(geom) = -1
c
      do i = 1, 3
        lattice_vectors(i,geom) = 0
        lattice_angles(i,geom) = 0
        do j = 1, 3
           amatrix(j,i,geom) = 0.0d0
           amatrix_inv(j,i,geom) = 0.0d0
        enddo
        amatrix(i,i,geom) = 1.0d0
        amatrix_inv(i,i,geom) = 1.0d0
      enddo
c
      end
      logical function geom_destroy(geom)
      implicit none
#include "geomP.fh"
#include "mafdecls.fh"
c
      integer geom              ! [input]
      logical geom_check_handle, geom_check_cent
      external geom_check_handle, geom_check_cent
c     
      geom_destroy = geom_check_handle(geom, 'geom_destroy')
      if (.not. geom_destroy) return
c
      active(geom) = .false.
      geom_destroy = .true.
      if (sym_center_map_handle(geom) .ne. -1) then
         if (.not. ma_free_heap(sym_center_map_handle(geom)))
     $        call errquit('geom_destroy: free of atom map', 0)
      endif
c
      end
      logical function geom_cart_set(geom, ncent, t, c, q)
      implicit none
#include "geomP.fh"
c     
      integer geom              ! [input]
      integer ncent             ! [input]
      character*(*) t(ncent)    ! [input]
      double precision c(3, ncent) ! [input]
      double precision q(ncent) ! [input]
      logical geom_check_handle, geom_check_cent
      external geom_check_handle, geom_check_cent
      integer i
c
      geom_cart_set = geom_check_handle(geom, 'geom_cart_set')
      if (.not. geom_cart_set) return
c
      if (ncent.le.0 .or. ncent.gt.max_cent) then
         write(6,*) ' geom_cart_set: too many centers ',ncent,
     $        names(geom)(1:lenn(geom))
         return
      endif
c
      ncenter(geom) = ncent
      ncenter_unique(geom) = ncent
      do i = 1, ncent
         tags(i,geom) = t(i)
         charge(i,geom) = q(i)
         coords(1,i,geom) = c(1,i)
         coords(2,i,geom) = c(2,i)
         coords(3,i,geom) = c(3,i)
         unique_cent(i,geom) = i
      enddo
c
c
c     compute effective nuclear repulsion energy, dipole and
c     interaction with external fields
c
      call geom_compute_values(geom)
c
      end
      logical function geom_cart_get(geom, ncent, t, c, q)
      implicit none
#include "geomP.fh"
c     
      integer geom              ! [input]
      integer ncent             ! [output]
      character*(*) t(ncent)    ! [output]
      double precision c(3, ncent) ! [output]
      double precision q(ncent) ! [output]
      logical geom_check_handle, geom_check_cent
      external geom_check_handle, geom_check_cent
      integer i
c
      geom_cart_get = geom_check_handle(geom, 'geom_cart_get')
      if (.not. geom_cart_get) return
c
      ncent = ncenter(geom)
      do i = 1, ncent
         t(i) = tags(i,geom)
         q(i) = charge(i,geom)
         c(1,i) = coords(1,i,geom)
         c(2,i) = coords(2,i,geom)
         c(3,i) = coords(3,i,geom)
      enddo
c
      end
      logical function geom_cent_get(geom, icent, t, c, q)
      implicit none
#include "geomP.fh"
c     
      integer geom              ! [input]
      integer icent             ! [input]
      character*(*) t           ! [output]
      double precision c(3)     ! [output]
      double precision q        ! [output]
      logical geom_check_handle, geom_check_cent
      external geom_check_handle, geom_check_cent
c
      geom_cent_get = geom_check_handle(geom, 'geom_cent_get')
      if (.not. geom_cent_get) return
      geom_cent_get = geom_check_cent(geom, 'geom_cent_get', icent)
      if (.not. geom_cent_get) return

c
      t = tags(icent,geom)
      c(1) = coords(1,icent,geom)
      c(2) = coords(2,icent,geom)
      c(3) = coords(3,icent,geom)
      q = charge(icent,geom)
      geom_cent_get = .true.
c
      end
      logical function geom_cent_set(geom, icent, t, c, q)
      implicit none
#include "geomP.fh"
c     
      integer geom              ! [input]
      integer icent             ! [input]
      character*(*) t           ! [input]
      double precision c(3)     ! [input]
      double precision q        ! [input]
      logical geom_check_handle, geom_check_cent
      external geom_check_handle, geom_check_cent
c
      geom_cent_set = geom_check_handle(geom, 'geom_cent_set')
      if (.not. geom_cent_set) return
      geom_cent_set = geom_check_cent(geom, 'geom_cent_set', icent)
      if (.not. geom_cent_set) return
c
      tags(icent,geom) = t
      coords(1,icent,geom) = c(1)
      coords(2,icent,geom) = c(2)
      coords(3,icent,geom) = c(3)
      charge(icent,geom) = q
c
c     compute effective nuclear repulsion energy, dipole and
c     interaction with external fields
c
      call geom_compute_values(geom)
c
      end
      logical function geom_ncent(geom, ncent)
      implicit none
#include "geomP.fh"
c     
      integer geom              ! [input]
      integer ncent             ! [output]
      logical geom_check_handle, geom_check_cent
      external geom_check_handle, geom_check_cent
c
      geom_ncent = geom_check_handle(geom, 'geom_ncent')
      if (.not. geom_ncent) return
      ncent = ncenter(geom)
c      
      end
      logical function geom_cent_tag(geom, icent, tag)
      implicit none
#include "geomP.fh"
c     
      integer geom              ! [input]
      integer icent             ! [input]
      character*(*) tag         ! [output]
      logical geom_check_handle, geom_check_cent
      external geom_check_handle, geom_check_cent
c
      geom_cent_tag = geom_check_handle(geom, 'geom_cent_tag')
      if (.not. geom_cent_tag) return
      geom_cent_tag = geom_check_cent(geom, 'geom_cent_tag', icent)
      if (.not. geom_cent_tag) return
c
      tag = tags(icent,geom)
      geom_cent_tag = .true.
c
      end
      logical function geom_latvec_set(geom, vectors)
      implicit none
      integer geom              ! [input]
      double precision vectors(3,3) ! [input]
c
#include "geomP.fh"
c     
      call errquit('geom_latvec_set: not yet!', 0)
      geom_latvec_set = .false.
      end
      logical function geom_efield_set(geom, ef)
      implicit none
#include "geomP.fh"
c     
      integer geom              ! [input]
      double precision ef       ! [input]
c
      call errquit('geom_efield_set: not yet!', 0)
c     call geom_set_values(geom)
      geom_efield_set = .false.
      end
      logical function geom_efield_get(geom, ef)
      implicit none
#include "geomP.fh"
c     
      integer geom              ! [input]
      double precision ef       ! [output]
c
      call errquit('geom_get_efield: not yet!', 0)
      geom_efield_get = .false.
      end
      logical function geom_print(geom)
      implicit none
#include "geomP.fh"
#include "util.fh"
c
c     Basic printing of cartesian geometry ... needs support for
c     user defined units, internal coords, different formats, ...
c
      integer geom              ! [input]
      integer icent
      logical geom_check_handle, geom_check_cent
      external geom_check_handle, geom_check_cent
      integer i
c
      if (.not. geom_check_handle(geom, 'geom_print')) then
         geom_print = .false.
         return
      endif
c
      write(6,1) names(geom)(1:lenn(geom)), trans(geom)(1:lent(geom))
 1    format(//'                      Geometry (au) "',a,'" -> "',a,'"')
      write(6,2)
 2    format(  '                      -------------'/)
      write(6,3) 
 3    format('  No.       Tag          Charge          X',
     $     '              Y              Z'/
     $       ' ---- ---------------- ---------- --------------',
     $     ' -------------- --------------')
      do icent = 1, ncenter(geom)
         write(6,4) icent, tags(icent,geom), charge(icent,geom),
     $        coords(1,icent,geom),coords(2,icent,geom),
     $        coords(3,icent,geom)
 4       format(' ',i4,' ',a16,' ',f10.6,3f15.8)
      enddo
c
      write(6,41) erep(geom)
 41   format(/' Effective nuclear repulsion energy (au) ', f18.10/)
c
c     Nuclear dipole moment is usually zero ... only print it if it
c     is not.
c
      if (ddot(3,ndipole(1,geom),1,ndipole(1,geom),1) .gt. 1.0d-10) then
         write(6,91)
 91      format('             Nuclear Dipole moment (au) ')
         write(6,101)
 101     format('             --------------------------')
         write(6,7)
 7       format('        X                 Y               Z'/
     $        ' ---------------- ---------------- ----------------')
         write(6,8) (ndipole(i,geom), i=1,3)
 8       format(3(1x,f16.10))
         write(6,*)
      endif
c
      if (oefield(geom)) then
         write(6,9)
 9       format('         Electric Field (au) ')
         write(6,10)
 10      format('         -------------------')
         write(6,7)
         write(6,8) (efield(i,geom), i=1,3)
         write(6,*)
      endif
c
***      call sym_print_all(geom)
c
      geom_print = .true.
c
      end
      logical function geom_tag_to_element(tag, symbol, element, atn)
      implicit none
#include "inp.fh"
#include "geomP.fh"
      character*(*) tag         ! [input]
      character*(*) symbol      ! [output]
      character*(*) element     ! [output]
      integer atn               ! [output]
c
c     
c     attempt to figure out which element a tag refers to
c     and return the symbol, name and atomic no.
c
      integer lbuf, ind
      character*16 buf
      character*1 sym1(14)      ! 1 character atomic symbols+atomic no.s
      integer atn1(14)
      data sym1/'H','B','C','N','O','F','P','S','K','V','Y','I','W','U'/
      data atn1/ 1 , 5 , 6 , 7 , 8 , 9 , 15, 16, 19, 23, 39, 53, 74, 92/

      geom_tag_to_element = .false.
c
c     eliminate conventions that refer to centers used for
c     computation purposes .. just bq for now
c
      buf = tag
      lbuf = inp_strlen(buf)
      if (lbuf .eq. 0) return
c
      call inp_lcase(buf)
      if (buf(1:2) .eq. 'bq') then 
         element = 'point charge' ! Note that false is returned
         symbol  = 'bq'
         atn     = 0
         return
      endif
c
c     Attempt to match the first 4 characters of the
c     full names of the elements
c
      atn = 0
      if (lbuf .ge. 4) then
         if (inp_match(nelements,.false.,buf(1:4),elements,ind)) then
            symbol  = symbols(ind)
            element = elements(ind)
            atn     = ind
            geom_tag_to_element = .true.
            return
         endif
      endif
c
c     Failed ... attempt to match the first two characters
c     against two character element names
c
      if (buf(2:2) .ne. ' ') then
         if (inp_match(nelements,.false.,buf(1:2),symbols,ind)) then
            symbol  = symbols(ind)
            element = elements(ind)
            atn     = ind
            geom_tag_to_element = .true.
            return
         endif
      endif
c
c     Last ditch attempt ... match against 1 character symbols
c
      if (inp_match(14, .false., buf(1:1), sym1, ind)) then
         ind = atn1(ind)
         symbol  = symbols(ind)
         element = elements(ind)
         atn     = ind
         geom_tag_to_element = .true.
         return
      endif
c
c     Nothing matched
c
      symbol = ' '
      element = ' '
      atn = 0
      return
c
      end
      logical function geom_charge_center(geom)
      implicit none
#include "geomP.fh"
      integer geom              ! [input]
      integer i, k
      double precision qsum, shift(3)
      logical geom_check_handle
      external geom_check_handle
c
c     Adjust the cartesian coordinates so that the nuclear
c     dipole moment is zero ... i.e., the origin of the 
c     coordinate system is at the center of charge
c
      geom_charge_center = geom_check_handle(geom,'geom_charge_center')
      if (.not. geom_charge_center) return
      call geom_compute_values(geom)
      qsum = 0.0d0
      do i = 1, ncenter(geom)
         qsum = qsum + charge(i,geom)
      enddo
c
      if (qsum .eq. 0.0d0) return ! System is charge neutral
c
      do k = 1, 3
         shift(k) = ndipole(k,geom)/qsum
      enddo
      do i = 1, ncenter(geom)
         do k = 1, 3
            coords(k,i,geom) = coords(k,i,geom) - shift(k)
         enddo
      enddo
c     
      call geom_compute_values(geom)
c
      end
      logical function geom_nuc_rep_energy(geom, energy)
#include "geomP.fh"
      integer geom              ! [input]
      double precision energy   ! [output]
      logical geom_check_handle
      external geom_check_handle
c
c     return the effective nuclear repulsion energy etc.
c
      geom_nuc_rep_energy = geom_check_handle(geom, 'geom_nuc_rep_e')
      if (.not. geom_nuc_rep_energy) return
      energy = erep(geom)
c
      end
      logical function geom_nuc_charge(geom, total_charge)
#include "geomP.fh"
      integer geom              ! [input]
      double precision total_charge ! [output]
      logical geom_check_handle
      external geom_check_handle
      integer i
c
c     return the sum of the nuclear charges
c
      geom_nuc_charge = geom_check_handle(geom, 'geom_nuc_charge')
      if (.not. geom_nuc_charge) return
c
      total_charge = 0.0d0
      do i = 1, ncenter(geom)
         total_charge = total_charge + charge(i,geom)
      enddo
c
      end
c
      logical function geom_verify_coords(geom)
      implicit none
c
c::functions      
      logical geom_ncent
      logical geom_cent_get
      external geom_ncent
      external geom_cent_get
c::passed
      integer geom ! geometry handle
c::local
      integer nat, iat, jat, num2compare, atomi, atomj
      integer i,j
      parameter (num2compare = 2)
      character*16 name(num2compare)
      double precision xyz(3,num2compare)
      double precision chg(num2compare)
      double precision dist_min
      double precision dist_thresh
      parameter (dist_thresh = 1.0d-03)
      double precision dist2
c
      dist2(i,j) = (xyz(1,i)-xyz(1,j))**2 +
     &    (xyz(2,i)-xyz(2,j))**2 +
     &    (xyz(3,i)-xyz(3,j))**2
c
      if(.not.geom_ncent(geom, nat))
     &    call errquit('geom_verify_coords: geom_ncent failed',911)
c
      atomi = 1
      atomj = 2
      dist_min = 56565.89d00
      do 00100 iat = 1,nat
        do 00200 jat = 1,iat
          if (jat.lt.iat) then
            if(.not.geom_cent_get
     &          (geom,iat,name(atomi),xyz(1,atomi),chg(atomi)))
     &          call errquit
     &          ('geom_verify_coords: geom_cent_get<1> failed',911)
            if(.not. geom_cent_get
     &          (geom,jat,name(atomj),xyz(1,atomj),chg(atomj)))
     &          call errquit
     &          ('geom_verify_coords: geom_cent_get<2> failed',911)
            dist_min = min(dist_min, dist2(atomi,atomj))
          endif
00200   continue
00100 continue
      dist_min = sqrt(dist_min)
      geom_verify_coords = dist_min.gt.dist_thresh
*      write(6,*)' distance minimum =',
*     &    dist_min, geom_verify_coords
c
      if (geom_verify_coords) return
      write(6,*)' ************ WARNING ******************' 
      write(6,*)' at least two atoms are at the same physical location'
      end
c
c
c---> new functions added on incorporation on symmetry and solid state codes
c
c
      logical function geom_systype_get(geom, itype)
#include "geomP.fh"
      integer geom              ! [input]
      integer itype
      logical geom_check_handle
      external geom_check_handle

      geom_systype_get = geom_check_handle(geom, 'geom_systype_get')
      if (.not. geom_systype_get) return
c
c--> make the assignment
c
      itype=isystype(geom)
c
      geom_systype_get = .true.
      end
      logical function geom_latvec_get(geom,vectors)
      implicit none
#include "geomP.fh"
c     
      integer geom,i              ! [input]
      double precision vectors(3) ! [output]
      logical geom_check_handle
      external geom_check_handle
c
      geom_latvec_get = geom_check_handle(geom, 'geom_latvec_get')
      if (.not. geom_latvec_get) return

      do i=1,3
        vectors(i)=lattice_vectors(i,geom)
      enddo
      geom_latvec_get = .true.
      end
      logical function geom_latang_get(geom,angles)
      implicit none
#include "geomP.fh"
c     
      integer geom,i              ! [input]
      double precision angles(3) ! [output]
      logical geom_check_handle
      external geom_check_handle

      geom_latang_get = geom_check_handle(geom, 'geom_latang_get')
      if (.not. geom_latang_get) return
c
      do i=1,3
        angles(i)=lattice_angles(i,geom)
      enddo
      geom_latang_get = .true. 
      end
      logical function geom_recipvec_get(geom,rvectors)
      implicit none
#include "geomP.fh"
c     
      integer geom,i               ! [input]
      double precision rvectors(3) ! [output]
      logical geom_check_handle
      external geom_check_handle
c
      geom_recipvec_get = geom_check_handle(geom, 'geom_recipvec_get')
      if (.not. geom_recipvec_get) return

      do i=1,3
        rvectors(i)=recip_lat_vectors(i,geom)
      enddo
      write(*,*) 'in geom printing rcectors'
      write(*,*) (rvectors(i),i=1,3)
      geom_recipvec_get = .true.
      end
      logical function geom_recipang_get(geom,rangles)
      implicit none
#include "geomP.fh"
c     
      integer geom,i               ! [input]
      double precision rangles(3) ! [output]
      logical geom_check_handle
      external geom_check_handle
c
      geom_recipang_get = geom_check_handle(geom, 'geom_recipang_get')
      if (.not. geom_recipang_get) return

      do i=1,3
        rangles(i)=recip_lat_angles(i,geom)
      enddo
      geom_recipang_get = .true.
      end
      logical function geom_volume_get(geom,volume)
      implicit none
#include "geomP.fh"
c     
      integer geom              ! [input]
      double precision volume   ! [output]
      logical geom_check_handle
      external geom_check_handle
c
      geom_volume_get = geom_check_handle(geom, 'geom_volume_get')
      if (.not. geom_volume_get) return

      volume=volume_direct(geom)

      end
      logical function geom_amatrix_get(geom,amat)
      implicit none
#include "geomP.fh"
c     
      integer geom,i,j               ! [input]
      double precision amat(3,3)     ! [output]
      logical geom_check_handle
      external geom_check_handle

      geom_amatrix_get = geom_check_handle(geom, 'geom_amatrix_get')
      if (.not. geom_amatrix_get) return
c
      do i=1,3
        do j=1,3
          amat(i,j)=amatrix(i,j,geom)
        enddo
      enddo
      end
      logical function geom_amatinv_get(geom,amatinv)
      implicit none
#include "geomP.fh"
c     
      integer geom,i,j               ! [input]
      double precision amatinv(3,3)  ! [output]
      logical geom_check_handle
      external geom_check_handle

      geom_amatinv_get = geom_check_handle(geom, 'geom_amatinv_get')
      if (.not. geom_amatinv_get) return
c
      do i=1,3
        do j=1,3
          amatinv(i,j)=amatrix_inv(i,j,geom)
        enddo
      enddo
      end




