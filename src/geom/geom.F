      block data geom_data
C$Id: geom.F,v 1.57 1997-02-08 00:12:25 d3g681 Exp $
      implicit none
#include "geomP.fh"
c     
      integer i                 ! For implied do in data staements
      data ngeom_rtdb /0/
      data active /max_geom*.false./
c
c--> names of the 3-dimensional space groups
c
      data (sym_spgnames(i),i=1,95) /
     & 'P1','P-1','P2','P2_1','C2',
     & 'Pm','Pc','Cm','Cc','P2/m',
     & 'P2_1/m','C2/m','P2/c','P2_1/c','C2/c',
     & 'P222','P222_1','P2_12_12','P2_12_12_1','C222_1',
     & 'C222','F222','I222','I2_12_12_1','Pmm2',
     & 'Pmc2_1','Pcc2','Pma2','Pca2_1','Pnc2',
     & 'Pmn2_1','Pba2','Pna2_1','Pnn2','Cmm2',
     & 'Cmc2_1','Ccc2','Amm2','Abm2','Ama2',
     & 'Aba2','Fmm2','Fdd2','Imm2','Iba2',
     & 'Ima2','Pmmm','Pnnn','Pccm','Pban',
     & 'Pmma','Pnna','Pmna','Pcca','Pbam',
     & 'Pccn','Pbcm','Pnnm','Pmmn','Pbcn',
     & 'Pbca','Pnma','Cmcm','Cmca','Cmmm',
     & 'Cccm','Cmma','Ccca','Fmmm','Fddd',
     & 'Immm','Ibam','Ibca','Imma','P4',
     & 'P4_1','P4_2','P4_3','I4','I4_1',
     & 'P-4','I-4','P4/m','P4_2/m','P4/n',
     & 'P4_2/n','I4/m','I4_1/a','P422','P42_12',
     & 'P4_122','P4_12_12','P4_222','P4_22_12','P4_322'/
      data (sym_spgnames(i),i=96,190)/   
     & 'P4_32_12','I422','I4_122','P4mm','P4bm',
     & 'P4_2cm','P4_2nm','P4cc','P4nc','P4_2mc',
     & 'P4_2bc','I4mm','I4cm','I4_1md','I4_1cd',
     & 'P-42m','P-42c','P-42_1m','P-42_1c','P-4m2',
     & 'P-4c2','P-4b2','P-4n2','I-4m2','I-4c2',
     & 'I-42m','I-42d','P4/mmm','P4/mcc','P4/nbm',
     & 'P4/nnc','P4/mbm','P4/mnc','P4/nmm','P4/ncc',
     & 'P4_2/mmc','P4_2/mcm','P4_2/nbc','P4_2/nnm','P4_2/mbc',
     & 'P4_2/mnm','P4_2/nmc','P4_2/ncm','I4/mmm','I4/mcm',
     & 'I4_1/amd','I4_1/acd','P3','P3_1','P3_2',
     & 'R3','P-3','R-3','P312','P321',
     & 'P3_112','P3_121','P3_212','P3_221','R32',
     & 'P3m1','P31m','P3c1','P31c','R3m',
     & 'R3c','P-31m','P-31c','P-3m1','P-3c1',
     & 'R-3m','R-3c','P6','P6_1','P6_5',
     & 'P6_2','P6_4','P6_3','P-6','P6/m',
     & 'P6_3/m','P622','P6_122','P6_522','P6_222',
     & 'P6_422','P6_322','P6mm','P6cc','P6_3cm',
     & 'P6_3mc','P-6m2','P-6c2','P-62m','P-62c'/
      data (sym_spgnames(i),i=191,230)/
     & 'P6/mmm','P6/mcc','P6_3/mcm','P6_3/mmc','P23',
     & 'F23','I23','P2_13','I2_13','Pm-3',
     & 'Pn-3','Fm-3','Fd-3','Im-3','Pa-3',
     & 'Ia-3','P432','P4_232','F432','F4_132',
     & 'I432','P4_332','P4_132','I4_132','P-43m',
     & 'F-43m','I-43m','P-43n','F-43c','I-43d',
     & 'Pm-3m','Pn-3n','Pm-3n','Pn-3m','Fm-3m',
     & 'Fm-3c','Fd-3m','Fd-3c','Im-3m','Ia-3d'/
c
c--> names of the molecular point groups
c
      data sym_molgnames/
     & 'C1','Cs','Ci','C2','C3',
     & 'C4','C5','C6','C7','C8',
     & 'D2','D3','D4','D5','D6',
     & 'C2v','C3v','C4v','C5v','C6v',
     & 'C2h','C3h','C4h','C5h','C6h',
     & 'D2h','D3h','D4h','D5h','D6h',
     & 'D8h','D2d','D3d','D4d','D5d',
     & 'D6d','S4','S6','S8','T',
     & 'Th','Td','O','Oh','I',
     & 'Ih'/

c
*rak:oldest:      data angstrom_to_au /1.8897265d0/
*rak:older:      data angstrom_to_au /1.8897266d0/
*. match inverse of new standard. 0.529177249
      data angstrom_to_au /1.88972598858d0/
      data isystype / max_geom*0/
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
      end if
c     
      end
      logical function geom_check_cent(geom, msg, icent)
      implicit none
#include "geomP.fh"
c     
      integer geom              ! [input]
      character*(*) msg         ! [input]
      integer icent             ! [input]
      logical status, geom_print
      external geom_print
c     
      geom_check_cent = icent.gt.0 .and. icent.le.ncenter(geom)
      if (.not. geom_check_cent) then
         write(6,*) msg,': icent invalid ', icent,
     $        names(geom)(1:lenn(geom))
         call geom_err_info(msg)
         status = geom_print(geom)
      end if
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
            end do
            geom_rtdb_in = .true.
         end if
      end if
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
      end do
c     
c     Name is not present ... add and rewrite info
c     
      if (ngeom_rtdb .eq. max_geom_rtdb) then
         write(6,*) ' geom_rtdb_add: too many geometries on rtdb ', name
         geom_rtdb_add = .false.
         return
      end if
      ngeom_rtdb = ngeom_rtdb + 1
      names_rtdb(ngeom_rtdb) = name
      lenr(ngeom_rtdb) = ln
c     
      if (.not. geom_rtdb_out(rtdb)) then
         write(6,*) ' geom_rtdb_add: rtdb error adding ', name(1:ln)
         geom_rtdb_add = .false.
         return
      end if
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
      end do
      write(6,1) info, ngeom
 1    format(' ',a,': open geometies: ',i2)
      ngeom = 0
      do geom = 1, max_geom
         if (active(geom)) then
            write(6,2) geom, info, names(geom)(1:lenn(geom)),
     $           trans(geom)(1:lent(geom))
 2          format(' ',i2,' ',a,': "',a, '" -> "', a,'"')
         end if
      end do
      if (ngeom_rtdb .gt. 0) then
         write(6,3) info, ngeom_rtdb
 3       format(' ',a,': geometries in last accessed data base: ', i2)
         do geom = 1, ngeom_rtdb
            write(6,4) names_rtdb(geom)(1:lenr(geom))
 4          format(' ',a)
         end do
      end if
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

*#include "global.fh"
*      integer node
c     
      integer rtdb              ! [input]
      integer geom              ! [input]
      character*(*) name        ! [input]
c     
      character*256 tmp
      integer k, nelem, ma_type
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
      tmp(k:) = ':masses'
      status = status .and. 
     $     rtdb_get(rtdb, tmp, mt_dbl, max_cent, geom_mass(1,geom))
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
      tmp(k:) = ':group name'
      status = status .and.
     $     rtdb_cget(rtdb, tmp, 1, group_name(geom))
      tmp(k:) = ':user units'
      status = status .and.
     $     rtdb_cget(rtdb, tmp, 1, user_units(geom))
      tmp(k:) = ':angstrom_to_au'
      status = status .and. 
     $     rtdb_get(rtdb, tmp, mt_dbl, 1, angstrom_to_au)
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
      tmp(k:) = ':include_bqbq'
      status = status .and. 
     $     rtdb_get(rtdb, tmp, mt_log, 1, include_bqbq(geom))
*rak: probably only incore data.  built by basis and ecp functionality
*rak:      tmp(k:) = ':ecp centers'
*rak:      status = status .and. 
*rak:     $     rtdb_get(rtdb, tmp, mt_log, ncenter(geom), oecpcent(1,geom))
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
      if (sym_num_ops(geom) .gt. 0) then
c
c     If loading into an old geometry free this memory
c
         if (sym_center_map_handle(geom) .ne. -1) then
            if (.not. ma_free_heap(sym_center_map_handle(geom)))
     $           call errquit('geom_rtdb_load: free of atom map', 0)
         end if
c
         tmp(k:) = ':map_atoms'
         status = status .and. 
     $        rtdb_ma_get(rtdb, tmp, ma_type, nelem,
     $        sym_center_map_handle(geom))
         if (nelem .ne. sym_num_ops(geom)*ncenter(geom)) call errquit
     $        ('geom_rtdb_load: invalid no. of element in sym map',
     $        nelem)
         if (.not. ma_get_index(sym_center_map_handle(geom),
     $        sym_center_map_index(geom)))call errquit
     $        ('geom_rtdb_load: bad ma handle for sym map', 0)
      else
         sym_center_map_handle(geom) = -1
         sym_center_map_index(geom) = 1 ! Not used but address is created
      endif
c     
      if (.not. status) then
         write(6,*) ' geom_rtdb_load: not found or rtdb corrupt: ',
     $        names(geom)(1:lenn(geom)), ' -> ',
     $        trans(geom)(1:lent(geom))
         call geom_err_info('geom_rtdb_load')
         geom_rtdb_load = .false.
         return
      end if
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
c setup geometry related stuff particular to the dimension of the system
c
      if (isystype(geom) .eq. 3) then
         call geom_3d(geom)
      elseif(isystype(geom).eq.2) then
         call geom_2d(geom)
      elseif(isystype(geom).eq.1) then
         call geom_1d(geom)
      endif
c
*      do node = 0, ga_nnodes()-1
*         call ga_sync
*         if (ga_nodeid() .eq. node) then
*            write(6,*) ' node ', ga_nodeid()
*            call sym_print_all(geom, .true., .true., .true., .true., .true.)
*            call util_flush(6)
*         endif
*         call ga_sync
*      enddo
         
c     
      end
      subroutine geom_compute_values(geom)
      implicit none
#include "geomP.fh"
#include "util.fh"
#include "inp.fh"
      integer geom              ! [input]
c
c     compute effective nuclear repulsion energy, dipole and
c     interaction with external fields
c
c     eventually need to also make the symmetry info consistent
c     and make internals/cartesians consistent
c
      double precision e, e_nd_ef, r, rx, ry, rz
      integer i, j
      logical j_is_atom, i_is_atom
      logical geom_tag_to_element
      external geom_tag_to_element
      logical is_atom
      is_atom(i) = (.not. inp_compare(.false., 'bq', tags(i,geom)(1:2)))
c
      e = 0.0d0
      ndipole(1,geom) = 0.0d0
      ndipole(2,geom) = 0.0d0
      ndipole(3,geom) = 0.0d0
c
c     compute nuclear dipole moment and usual nuclear repulsion energy
c
      do i = 1,ncenter(geom)
         i_is_atom = is_atom(i)
         if (include_bqbq(geom) .or. i_is_atom) then
            do j = 1, 3
               ndipole(j,geom) = ndipole(j,geom) +
     $              charge(i,geom)*coords(j,i,geom)
            end do
         endif
         do j = i+1, ncenter(geom)
            j_is_atom = is_atom(j)
            if (include_bqbq(geom) .or. (i_is_atom.or.j_is_atom)) then

*               r = dsqrt(
*     $              (coords(1,i,geom)-coords(1,j,geom))**2 + 
*     $              (coords(2,i,geom)-coords(2,j,geom))**2 + 
*     $              (coords(3,i,geom)-coords(3,j,geom))**2)
              rx = coords(1,i,geom)-coords(1,j,geom)
              rx = rx*rx
              ry = coords(2,i,geom)-coords(2,j,geom)
              ry = ry*ry
              rz = coords(3,i,geom)-coords(3,j,geom)
              rz = rz*rz
              r  = sqrt(rx+ry+rz)
               e = e + charge(i,geom)*charge(j,geom)/r
            endif
         end do
      end do
c
c     add in interaction of nuclear dipole with external field
c
      e_nd_ef = ddot(3, ndipole(1,geom), 1, efield(1,geom), 1)
*:debug-s      
*debug:      write(6,*)' interaction of nuclear dipole ',
*debug:     &    'with external field is ',e_nd_ef
*:debug-e
      e = e + e_nd_ef
c
      erep(geom) = e
c
      if(isystype(geom).eq.0) then
         call sym_init_inv_op(geom)
      endif
c
      end
      logical function geom_include_bqbq(geom)
      implicit none
#include "geomP.fh"
      integer geom
      logical geom_check_handle
      external geom_check_handle
c
      if (.not. geom_check_handle(geom, 'geom_include_bqbq'))
     $     call errquit('geom_include_bqbq: bad handle',0)
      geom_include_bqbq = include_bqbq(geom)
c
      end
      logical function geom_set_bqbq(geom, value)
      implicit none
#include "geomP.fh"
      logical value
      integer geom
      logical geom_check_handle
      external geom_check_handle
c
      geom_set_bqbq = geom_check_handle(geom, 'geom_set_bqbq')
      if (.not. geom_set_bqbq) return
      include_bqbq(geom) = value
      call geom_compute_values(geom)
c
      end
      logical function geom_rtdb_store(rtdb, geom, name)
      implicit none
#include "geomP.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "context.fh"
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
      end if
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
      tmp(k:) = ':masses'
      status = status .and. 
     $     rtdb_put(rtdb, tmp, mt_dbl, ncenter(geom),
     $     geom_mass(1,geom))
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
      tmp(k:) = ':group name'
      status = status .and.
     $     rtdb_cput(rtdb, tmp, 1, group_name(geom))
      tmp(k:) = ':user units'
      status = status .and.
     $     rtdb_cput(rtdb, tmp, 1, user_units(geom))
      tmp(k:) = ':angstrom_to_au'
      status = status .and.
     &     rtdb_put(rtdb, tmp, mt_dbl, 1, angstrom_to_au)
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
      tmp(k:) = ':include_bqbq'
      status = status .and. 
     $     rtdb_put(rtdb, tmp, mt_log, 1, include_bqbq(geom))
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
      if (sym_num_ops(geom) .gt. 0) then
         tmp(k:) = ':map_atoms'
         status = status .and. 
     $        rtdb_put(rtdb, tmp, mt_int,
     $        ncenter(geom)*sym_num_ops(geom),
     $        int_mb(sym_center_map_index(geom)))
      endif
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
      end if
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
      end do
      write(6,*) ' geom_rtdb_delete: no such geometry ',
     $     name(1:inp_strlen(name)),' -> ', translation(1:lt)
 10   do geom2 = geom+1, ngeom_rtdb
         names_rtdb(geom2-1) = names_rtdb(geom2)
      end do
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
      tmp(k:) = ':masses'
      status = status .and. rtdb_delete(rtdb, tmp)
      tmp(k:) = ':efield'
      status = status .and. rtdb_delete(rtdb, tmp)
      tmp(k:) = ':tags'
      status = status .and. rtdb_delete(rtdb, tmp)
      tmp(k:) = ':lattice vectors'
      status = status .and. rtdb_delete(rtdb, tmp)
      tmp(k:) = ':lattice angles'
      status = status .and. rtdb_delete(rtdb, tmp)
      tmp(k:) = ':system type'
      status = status .and. rtdb_delete(rtdb, tmp)
      tmp(k:) = ':no. unique centers'
      status = status .and. rtdb_delete(rtdb, tmp)
      tmp(k:) = ':group number'
      status = status .and. rtdb_delete(rtdb, tmp)
      tmp(k:) = ':group name'
      status = status .and. rtdb_delete(rtdb, tmp)
      tmp(k:) = ':user units'
      status = status .and. rtdb_delete(rtdb, tmp)
      tmp(k:) = ':angstrom_to_au'
      status = status .and. rtdb_delete(rtdb, tmp)
      tmp(k:) = ':setting number'
      status = status .and. rtdb_delete(rtdb, tmp)
      tmp(k:) = ':recip vectors'
      status = status .and. rtdb_delete(rtdb, tmp)
      tmp(k:) = ':recip angles'
      status = status .and. rtdb_delete(rtdb, tmp)
      tmp(k:) = ':direct volume'
      status = status .and. rtdb_delete(rtdb, tmp)
      tmp(k:) = ':unique centers'
      status = status .and. rtdb_delete(rtdb, tmp)
      tmp(k:) = ':include_bqbq'
      status = status .and. rtdb_delete(rtdb, tmp)
      tmp(k:) = ':num_operators'
      status = status .and. rtdb_delete(rtdb, tmp)
      tmp(k:) = ':operators'
      if (sym_num_ops(geom) .gt. 0) then
         tmp(k:) = ':map_atoms'
      endif
c
c     check status of all rtdb operations
c
      if (.not. status) then
         write(6,*) ' geom_rtdb_delete: rtdb corrupt ',translation(1:lt)
         call geom_err_info('geom_rtdb_delete')
         geom_rtdb_delete = .false.
         return
      end if
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
      external geom_data  ! This for T3D linker
c
c     Assign the next free slot for a geometry
c
      do geom = 1, max_geom
         if (.not. active(geom)) goto 10
      end do
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
      sym_center_map_index(geom) = 1
      group_name(geom) = 'C1'
      sym_num_ops(geom) = 0
      user_units(geom) = 'angstrom'
      include_bqbq(geom) = .false.
c
      do i = 1, 3
        lattice_vectors(i,geom) = 0
        lattice_angles(i,geom) = 0
        do j = 1, 3
           amatrix(j,i,geom) = 0.0d0
           amatrix_inv(j,i,geom) = 0.0d0
           bmatrix(j,i,geom) = 0.0d0
        end do
        amatrix(i,i,geom) = 1.0d0
        amatrix_inv(i,i,geom) = 1.0d0
        bmatrix(i,i,geom) = 1.0d0
      end do
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
      end if
c
      end
      logical function geom_cart_set(geom, ncent, t, c, q)
      implicit none
#include "geomP.fh"
c     
      integer geom              ! [input]
      integer ncent             ! [input]
      character*16 t(ncent)     ! [input]
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
      end if
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
      end do
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
      character*16 t(ncent)     ! [output]
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
      end do
c
      end
      logical function geom_cent_get(geom, icent, t, c, q)
      implicit none
#include "geomP.fh"
c     
      integer geom              ! [input]
      integer icent             ! [input]
      character*16 t            ! [output]
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
      character*16 t            ! [input]
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
      logical function geom_ncent_unique(geom, ncent)
      implicit none
#include "geomP.fh"
c     
      integer geom              ! [input]
      integer ncent             ! [output]
      logical geom_check_handle, geom_check_cent
      external geom_check_handle, geom_check_cent
c
      geom_ncent_unique = geom_check_handle(geom, 'geom_ncent_unique')
      if (.not. geom_ncent_unique) return
      ncent = ncenter_unique(geom)
c      
      end
      logical function geom_cent_tag(geom, icent, tag)
      implicit none
#include "geomP.fh"
c     
      integer geom              ! [input]
      integer icent             ! [input]
      character*16 tag          ! [output]
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
      double precision ef(3)    ! [output]
      logical geom_check_handle
      external geom_check_handle
      integer i
c     
      if (.not. geom_check_handle(geom, 'geom_efield_get')) then
         geom_efield_get = .false.
         return
      end if
c
      if (oefield(geom)) then
         do i = 1, 3
            ef(i) = efield(i,geom)
         end do
      else
         do i = 1, 3
            ef(i) = 0.0d0
         end do
      endif
      geom_efield_get = .true.
      end
      logical function geom_print_xyz(geom, unit)
      implicit none
#include "geomP.fh"
#include "util.fh"
#include "inp.fh"
      integer geom, unit
      integer i, j, icent
      double precision scale, tmp(3)
      logical geom_check_handle
      external geom_check_handle
c     
      geom_print_xyz = .true.
      if (.not. geom_check_handle(geom, 'geom_print_xyz')) then
         geom_print_xyz = .false.
         return
      end if
c
      scale = 1.0d0 / angstrom_to_au
c
      write(unit,1) ncenter(geom)
 1    format(1x,i5)
      write(unit,2) names(geom)(1:inp_strlen(names(geom)))
 2    format(1x,a)
      do icent = 1, ncenter(geom)
c
c     Convert from cartesian to crystallographic coordinates
c
         do i = 1, 3
            tmp(i) = 0.0d0
            do j = 1, 3
               tmp(i) = tmp(i) + 
     $              amatrix_inv(i,j,geom)*coords(j,icent,geom)
            end do
            tmp(i) = tmp(i)*scale ! Scale to angstrom
         end do
         write(unit,3) tags(icent,geom), (tmp(j),j=1,3)
 3       format(1x,a16,1x,3f15.8)
      end do
c
      end
      logical function geom_print(geom)
      implicit none
#include "geomP.fh"
#include "util.fh"
#include "inp.fh"
c     
c     Basic printing of cartesian geometry 
c     needs support for internal coords, different formats, ...
c     
      integer geom              ! [input]
      integer icent, jcent
      logical geom_check_handle, geom_check_cent, geom_get_user_scale
      external geom_check_handle, geom_check_cent, geom_get_user_scale
      integer i
      double precision scale, tmp(3)
      character*80 buf
      logical oprint_uniq
c
      if (.not. geom_check_handle(geom, 'geom_print')) then
         geom_print = .false.
         return
      end if
c
c     All of the code seems to be commented out except for
c     molecules so just return if this is not a molecule (RJH)
c
c     ... it would be nice to have one routine that prints all
c     possible geometries but ...
c
      if (isystype(geom) .ne. 0) return
c     
      if (.not. geom_get_user_scale(geom, scale))
     $     call errquit('geom_print: user units?',0)
c     
      buf = ' '
      write(buf,1) 'Geometry',
     $        names(geom)(1:lenn(geom)), 
     $        trans(geom)(1:lent(geom))
 1    format(a,' "',a,'" -> "',a,'"')
      write(6,*)
      write(6,*)
      call util_print_centered(6,buf,40,.true.)
      write(6,*)
      write(6,2) user_units(geom)(1:inp_strlen(user_units(geom))),
     $     scale
 2    format(' Output coordinates in ', a,
     $     ' (scale by ',f12.9,' to convert to a.u.)')
      if (include_bqbq(geom))
     $     write(6,*) ' Include Bq-Bq interactions'
c
      write(6,*)
c     
      write(6,3) 
 3    format('  No.       Tag          Charge          X',
     $     '              Y              Z'/
     $     ' ---- ---------------- ---------- --------------',
     $     ' -------------- --------------')
      do icent = 1, ncenter(geom)
         do i = 1, 3
            tmp(i) = coords(i,icent,geom)/scale ! Scale units as necessary
         end do
         write(6,4) icent, tags(icent,geom), charge(icent,geom),
     $        (tmp(i),i=1,3)
 4       format(' ',i4,' ',a16,' ',f10.4,3f15.8)
      end do
c     
c     Only print out the masses for unique tags ... the structure
c     should actually only store the data for unique tags.
c     Also, keep all common output within 80 columns
c     
      write(6,*)
      write(6,*) '     Atomic Mass '
      write(6,*) '     ----------- '
      write(6,*)
      do icent = 1, ncenter(geom)
	 do jcent = 1, icent-1
	    if (tags(icent,geom) .eq. tags(jcent,geom)) goto 765
         enddo
	 write(6,44) tags(icent,geom), geom_mass(icent,geom)
 44      format('      ',a16,' ',f10.6)
 765  enddo
      write(6,*)
c     
      write(6,41) erep(geom)
 41   format(/' Effective nuclear repulsion energy (a.u.) ', f18.10/)
c     
      write(6,91)
 91   format('            Nuclear Dipole moment (a.u.) ')
      write(6,101)
 101  format('            ----------------------------')
      write(6,7)
 7    format('        X                 Y               Z'/
     $     ' ---------------- ---------------- ----------------')
      write(6,8) (ndipole(i,geom), i=1,3)
 8    format(3(1x,f16.10))
      write(6,*)
c     
      oprint_uniq = sym_num_ops(geom) .gt. 0
      if (oprint_uniq) then
         call sym_print_all(geom, .true., oprint_uniq, .false., .false.,
     $        .false.)
      endif
c     
      geom_print = .true.
c     
      end
      logical function geom_default_charge_with_ecp(atn, q)
      implicit none
c     
      integer atn               ! [input] atomic number
      double precision q        ! [output] charge
c
c     return a default for the effective nuclear charge
c     if an ecp is placed on a atom with atomic number atn
c
c     This is just a first guess at this routine
c
      geom_default_charge_with_ecp = .true.
      if (atn .le. 2) then
         q = atn
      else if (atn .le. 10) then
         q = atn - 2
      else if (atn .le. 18) then
         q = atn - 10
      else
         geom_default_charge_with_ecp = .false.
      endif
c
      end
      logical function geom_tag_to_element(tag, symbol, element, atn)
      implicit none
#include "inp.fh"
#include "geomP.fh"
      character*2 symbols(nelements)
      character*16 elements(nelements)
      character*16 tag          ! [input]
      character*(*) symbol      ! [output]
      character*(*) element     ! [output]
      integer atn               ! [output]
c
c     attempt to figure out which element a tag refers to
c     and return the symbol, name and atomic no.
c
      integer lbuf, ind
      character*16 buf
      character*1 sym1(14)      ! 1 character atomic symbols+atomic no.s
      integer atn1(14)
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
      end if
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
         end if
      end if
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
         end if
      end if
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
      end if
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
      end do
c
      if (qsum .eq. 0.0d0) return ! System is charge neutral
c
      do k = 1, 3
         shift(k) = ndipole(k,geom)/qsum
      end do
      do i = 1, ncenter(geom)
         do k = 1, 3
            coords(k,i,geom) = coords(k,i,geom) - shift(k)
         end do
      end do
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
      end do
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
          end if
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
      end do
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
      end do
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
      end do
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
      end do
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
        end do
      end do
      end
      logical function geom_bmatrix_get(geom,bmat)
      implicit none
#include "geomP.fh"
c     
      integer geom,i,j               ! [input]
      double precision bmat(3,3)     ! [output]
      logical geom_check_handle
      external geom_check_handle

      geom_bmatrix_get = geom_check_handle(geom, 'geom_bmatrix_get')
      if (.not. geom_bmatrix_get) return
c
      do i=1,3
        do j=1,3
          bmat(i,j)=bmatrix(i,j,geom)
        end do
      end do
      end
      logical function geom_amatinv_get(geom,amatinv)
      implicit none
#include "geomP.fh"
c     
      integer geom,i,j               ! [input]
      double precision amatinv(3,3)  ! [output]
      logical geom_check_handle
      external geom_check_handle
c
      geom_amatinv_get = geom_check_handle(geom, 'geom_amatinv_get')
      if (.not. geom_amatinv_get) return
c
      do i=1,3
        do j=1,3
          amatinv(i,j)=amatrix_inv(i,j,geom)
        end do
      end do
      end
      logical function geom_uniquecent_get(geom,ncent,uniquecent)
      implicit none
#include "geomP.fh"
c     
      integer geom,ncent             ! [input]
      integer uniquecent(ncent)      ! [output] indicies of unique centers
      integer i
      logical geom_check_handle
      external geom_check_handle
c
      geom_uniquecent_get=geom_check_handle(geom,'geom_uniquecent_get')
      if (.not. geom_uniquecent_get) return
c
      do i=1,ncent
         uniquecent(i)=unique_cent(i,geom)
      enddo
      end
      logical function geom_get_user_scale(geom, scale)
      implicit none
#include "geomP.fh"
      integer geom              ! [input]
      double precision scale    ! [output] Conversion from user units to A.U.
      logical geom_check_handle
      external geom_check_handle
c
      geom_get_user_scale = 
     $     geom_check_handle(geom, 'geom_get_user_scale')
c
      if (user_units(geom) .eq. 'a.u.') then
         scale = 1.0d0
      else if (user_units(geom) .eq. 'angstrom') then
         scale = angstrom_to_au
      else if (user_units(geom) .eq. 'nanometer') then
         scale = angstrom_to_au * 10.0d0
      else if (user_units(geom) .eq. 'picometer') then
         scale = angstrom_to_au * 0.01d0
      else
         call errquit('geom_get_user_scale: unknown units',0)
      endif
c
      end
      logical function geom_set_user_units(geom, units)
      implicit none
#include "geomP.fh"
      integer geom
      character*(*) units
      logical geom_check_handle
      external geom_check_handle
c
      geom_set_user_units =
     $     geom_check_handle(geom, 'geom_set_user_units')
      user_units(geom) = units
c
      end
      logical function geom_get_user_units(geom, units)
      implicit none
#include "geomP.fh"
      integer geom
      character*(*) units       ! [output] 
      logical geom_check_handle
      external geom_check_handle
c
      geom_get_user_units =
     $     geom_check_handle(geom, 'geom_get_user_units')
      units = user_units(geom)
c
      end
      logical function geom_tag_to_default_mass(tag,mass)
      implicit none
c
c this routine takes a tag matches it to the atomic number
c and returns the default atomic mass.
c
      character*16 tag          ! [input] geometry tag
      double precision mass  ! [output] corresponding elemental default mass
c
      logical geom_tag_to_element
      external geom_tag_to_element
      logical geom_atn_to_default_mass
      external geom_atn_to_default_mass
c
      character*2 tag_symbol
      character*16 tag_element
      integer tag_atomic_number
c
      geom_tag_to_default_mass = .false.
c
      if (.not. geom_tag_to_element(tag,tag_symbol, tag_element,
     &    tag_atomic_number)) call errquit
     &    ('geom_tag_to_default_mass: geom_tag_to_element failed ?',
     &    911)
      geom_tag_to_default_mass =
     &    geom_atn_to_default_mass(tag_atomic_number,mass)
      end
      logical function geom_atn_to_default_mass(atn,mass)
c
c This routine returns the default atomic mass from based on the atomic
c number.  The mass for each element comes from the book "The Elements" 
c by John Emsley, Oxford University Press, (C) 1989, ISBN 0-19-855237-8
c The specific mass chosen was the most abundant isotope with a known mass. 
c When the abundance was equal the isotope with the longest half life was 
c used.  
c
c RAK 11/95 PNNL/EMSL/HPCCG
c
      implicit none
#include "geomP.fh"
c
      integer atn             ! [input] atomic number of element
      double precision mass   ! [output] default elemental atomic mass.
c
      double precision def_masses(nelements)
c
      integer i
c
      data (def_masses(i),i=1,50) /
     &  1.007825d0, 4.0026d0,    7.016d0,    9.01218d0, 11.00931d0,
     & 12.0d0,     14.00307d0,  15.99491d0, 18.9984d0,  19.99244d0,
     & 22.9898d0,  23.98504d0,  26.98154d0, 27.97693d0, 30.97376d0,
     & 31.97207d0, 34.96885d0,  39.9624d0,  38.96371d0, 39.96259d0,
     & 44.95592d0, 45.948d0,    50.9440d0,  51.9405d0,  54.9381d0,
     & 55.9349d0,  58.9332d0,   57.9353d0,  62.9298d0,  63.9291d0,
     & 68.9257d0,  73.9219d0,   74.9216d0,  79.9165d0,  78.9183d0,
     & 83.912d0,   84.9117d0,   87.9056d0,  88.9054d0,  89.9043d0,
     & 92.9060d0,  97.9055d0,   97.9072d0, 101.9037d0, 102.9048d0,
     &105.9032d0, 106.90509d0, 113.9036d0, 114.9041d0, 117.9018d0/
      data (def_masses(i),i=51,103) /
     & 120.9038d0, 129.9067d0, 126.9004d0, 131.9042d0, 132.9051d0,
     & 137.9050d0, 138.9061d0, 139.9053d0, 140.9074d0, 143.9099d0,
     & 144.9128d0, 151.9195d0, 152.9209d0, 157.9241d0, 159.9250d0,
     & 163.9288d0, 164.9303d0, 165.9304d0, 168.9344d0, 173.9390d0,
     & 174.9409d0, 179.9468d0, 180.948d0,  183.9510d0, 186.9560d0,
     & 189.9586d0, 192.9633d0, 194.9648d0, 196.9666d0, 201.9706d0,
     & 204.9745d0, 207.9766d0, 208.9804d0, 209.9829d0, 210.9875d0,
     & 222.0175d0, 223.0198d0, 226.0254d0, 227.0278d0, 232.0382d0,
     & 231.0359d0, 238.0508d0, 237.048d0,  244.0642d0, 243.0614d0,
     & 247.0704d0, 247.0702d0, 249.0748d0, 254.0881d0, 253.0852d0,
     & 258.0d0,    259.0d0,    260.0/
c
      geom_atn_to_default_mass = .false.
c
      if (atn.lt.0) call errquit
     &    ('geom_atn_to_default_mass: negative atomic number',atn)
      if (atn.gt.nelements) call errquit
     &    ('geom_atn_to_default_mass: atomic number too large',atn)
c
      if (atn.eq.0) then
        mass = 0.0d00  ! Bq centers have no mass
      else
        mass = def_masses(atn)
      endif
      geom_atn_to_default_mass = .true.
c
      end
      logical function geom_masses_set(geom, ncent, masses)
      implicit none
#include "geomP.fh"
c
      integer geom                   ! [input] geometry handle
      integer ncent                  ! [input] number of centers
      double precision masses(ncent) ! [input] mass on each center
c
      integer i
c
      logical geom_check_handle
      external geom_check_handle
c
      geom_masses_set = geom_check_handle(geom, 'geom_masses_set')
      if (.not. geom_masses_set) return
c
      if (ncent.le.0 .or. ncent.gt.max_cent) then
         write(6,*) ' geom_masses_set: too many centers ',ncent,
     $        names(geom)(1:lenn(geom))
         return
      end if
c
      do i = 1, ncent
        geom_mass(i,geom) = masses(i)
      enddo
c
      end      
      logical function geom_masses_get(geom, ncent, masses)
      implicit none
#include "geomP.fh"
c
      integer geom                   ! [input] geometry handle
      integer ncent                  ! [input] number of centers
      double precision masses(ncent) ! [output] mass on each center
c
      integer i
c
      logical geom_check_handle
      external geom_check_handle
c
      geom_masses_get = geom_check_handle(geom, 'geom_masses_get')
      if (.not. geom_masses_get) return
c
      ncent = ncenter(geom)
      do i = 1, ncent
        masses(i) = geom_mass(i,geom)
      enddo
c
      end      
      logical function geom_mass_set(geom, icent, mass)
      implicit none
#include "geomP.fh"
c
      integer geom            ! [input] geometry handle
      integer icent           ! [input] number of center for mass
      double precision mass   ! [input] mass on center icent
c
      logical geom_check_handle
      external geom_check_handle
c
      geom_mass_set = geom_check_handle(geom, 'geom_mass_set')
      if (.not. geom_mass_set) return
c
      if (icent.le.0 .or. icent.gt.ncenter(geom)) then
         write(6,*) ' geom_mass_set: icent out of range',icent,
     &        ncenter(geom),
     $        names(geom)(1:lenn(geom))
         return
      end if
c
      geom_mass(icent,geom) = mass
c
      end      
      logical function geom_mass_get(geom, icent, mass)
      implicit none
#include "geomP.fh"
c
      integer geom            ! [input] geometry handle
      integer icent           ! [input] number of center for mass
      double precision mass   ! [output] mass on center icent
c
      logical geom_check_handle
      external geom_check_handle
c
      geom_mass_get = geom_check_handle(geom, 'geom_mass_get')
      if (.not. geom_mass_get) return
c
      if (icent.le.0 .or. icent.gt.ncenter(geom)) then
         write(6,*) ' geom_mass_get: icent out of range',icent,
     &        ncenter(geom),
     $        names(geom)(1:lenn(geom))
         return
      end if
c
      mass = geom_mass(icent,geom)
c
      end      
      logical function geom_set_ang2au(geom,value)
      implicit none
#include "geomP.fh"
c::functions
      logical geom_check_handle
      external geom_check_handle
c::passed
      integer geom               ! [input] geometry handle
      double precision value     ! [input] converts angstroms to au value ~1.8......
c
      geom_set_ang2au = geom_check_handle(geom,'geom_set_ang2au')
      if (.not. geom_set_ang2au) return
c
      angstrom_to_au = value
c
      end
      logical function geom_get_ang2au(geom,value)
      implicit none
#include "geomP.fh"
c::functions
      logical geom_check_handle
      external geom_check_handle
c::passed
      integer geom               ! [input] geometry handle
      double precision value     ! [output] converts angstroms to au value ~1.8......
c
      geom_get_ang2au = geom_check_handle(geom,'geom_get_ang2au')
      if (.not. geom_get_ang2au) return
c
      value = angstrom_to_au
c
      end
      logical function geom_set_au2ang(geom,value)
      implicit none
#include "geomP.fh"
c::functions
      logical geom_check_handle
      external geom_check_handle
c::passed
      integer geom               ! [input] geometry handle
      double precision value     ! [input] converts au to angstroms value ~0.52917
c
      geom_set_au2ang = geom_check_handle(geom,'geom_set_au2ang')
      if (.not. geom_set_au2ang) return
c
      angstrom_to_au = 1.0d00/value
c
      end
      logical function geom_get_au2ang(geom,value)
      implicit none
#include "geomP.fh"
c::functions
      logical geom_check_handle
      external geom_check_handle
c::passed
      integer geom               ! [input] geometry handle
      double precision value     ! [output] converts au to angstroms value ~0.52917
c
      geom_get_au2ang = geom_check_handle(geom,'geom_get_au2ang')
      if (.not. geom_get_au2ang) return
c
      value = 1.0d00/angstrom_to_au
c
      end
      logical function geom_ecp_allset(geom,ncenter_in,oecp)
      implicit none
c
#include "geomP.fh"
c
      integer geom                 ! [input] geometry handle
      integer ncenter_in           ! [input] number of centers
      logical oecp(ncenter_in)     ! [input] array of T/F for having ECPs
c
      logical geom_check_handle
      external geom_check_handle
c
      integer icenter
c
      geom_ecp_allset = geom_check_handle(geom, 'geom_ecp_allset')
c
      if (ncenter_in.ne.ncenter(geom)) call errquit
     &    (' too many or to few centers specified delta=',
     &    (ncenter(geom)-ncenter_in))
c
      do icenter = 1,ncenter_in
        oecpcent(icenter,geom) = oecp(icenter)
      enddo
c
      end
      logical function geom_ecp_allget(geom,ncenter_in,oecp)
      implicit none
c
#include "geomP.fh"
c
      integer geom                 ! [input] geometry handle
      integer ncenter_in           ! [input] number of centers
      logical oecp(ncenter_in)     ! [input] array of T/F for having ECPs
c
      logical geom_check_handle
      external geom_check_handle
c
      integer icenter
c
      geom_ecp_allget = geom_check_handle(geom, 'geom_ecp_allget')
c
      if (ncenter_in.ne.ncenter(geom)) call errquit
     &    (' too many or to few centers specified delta=',
     &    (ncenter(geom)-ncenter_in))
c
      do icenter = 1,ncenter_in
        oecp(icenter)= oecpcent(icenter,geom) 
      enddo
c
      end

      logical function geom_ecp_set(geom,icent,oecp)
      implicit none
c
#include "geomP.fh"
c
      integer geom     ! [input] geometry handle
      integer icent    ! [input] number of center to use
      logical oecp     ! [input] T/F for having ECPs
c
      logical geom_check_handle
      external geom_check_handle
c
      geom_ecp_set = geom_check_handle(geom, 'geom_ecp_set')
c
      if (.not.(icent.gt.0.and.icent.le.ncenter(geom)))
     &    call errquit
     &    (' icent out of range ncenter = ',
     &    ncenter(geom))
c
      oecpcent(icent,geom) = oecp
c
      end
      logical function geom_ecp_get(geom,icent)
      implicit none
c
#include "geomP.fh"
c
      integer geom     ! [input] geometry handle
      integer icent    ! [input] number of center to use
* return call is  [output] T/F for having ECPs
c
      logical geom_check_handle
      external geom_check_handle
c
      geom_ecp_get = geom_check_handle(geom, 'geom_ecp_get')
c
      if (.not.(icent.gt.0.and.icent.le.ncenter(geom)))
     &    call errquit
     &    (' icent out of range ncenter = ',
     &    ncenter(geom))
c
      geom_ecp_get = oecpcent(icent,geom) 
c
      end
      logical function geom_ncent_ecp(geom, ncent_ecp)
      implicit none
#include "geomP.fh"
c     
      integer geom              ! [input]
      integer ncent_ecp             ! [output]
      logical geom_check_handle
      external geom_check_handle
c
      integer icent
c
      geom_ncent_ecp = geom_check_handle(geom, 'geom_ncent_ecp')
      if (.not. geom_ncent_ecp) return
      ncent_ecp = 0
      do icent = 1,ncenter(geom)
        if (oecpcent(icent,geom)) ncent_ecp = ncent_ecp + 1
      enddo
c      
      end
      logical function geom_coords_ecp(geom, coords_ecp, ncent_in)
      implicit none
#include "geomP.fh"
c     
      integer geom                              ! [input]
      integer ncent_in                          ! [input]
      double precision coords_ecp(3,ncent_in)   ! [output]
c
      logical geom_check_handle
      external geom_check_handle
c
      integer icent, ncent_ecp
c
      geom_coords_ecp = geom_check_handle(geom, 'geom_coords_ecp')
      if (.not. geom_coords_ecp) return
      ncent_ecp = 0
      do icent = 1,ncenter(geom)
        if (oecpcent(icent,geom)) then
          ncent_ecp = ncent_ecp + 1
          if (ncent_ecp.gt.ncent_in) call errquit
     &        ('geom_coords_ecp: number of ecp centers is greater'//
     &         ' than the coord array dimension which is:',ncent_in)
*          write(6,*)' geom       = ',geom
*          write(6,*)' ncent_ecp  = ',ncent_ecp
*          write(6,*)' icent      = ',icent
*          write(6,*)' coords geom 1',coords(1,icent,geom)
*          write(6,*)' coords geom 2',coords(2,icent,geom)
*          write(6,*)' coords geom 3',coords(3,icent,geom)
          
          coords_ecp(1,ncent_ecp) = coords(1,icent,geom)
          coords_ecp(2,ncent_ecp) = coords(2,icent,geom)
          coords_ecp(3,ncent_ecp) = coords(3,icent,geom)

*          write(6,*)' coords ecp 1',coords_ecp(1,ncent_ecp)
*          write(6,*)' coords ecp 2',coords_ecp(2,ncent_ecp)
*          write(6,*)' coords ecp 3',coords_ecp(3,ncent_ecp)
        endif
      enddo
*      write(6,*)' coordinates inside geom_coords_ecp'
*      call output(coords_ecp,1,3,1,ncent_ecp,3,ncent_ecp,1)
c      
      end
      logical function geom_any_ecp(geom)
      implicit none
#include "geomP.fh"
c     
      integer geom              ! [input]
      logical geom_check_handle
      external geom_check_handle
c
      integer icent
c
      geom_any_ecp = geom_check_handle(geom, 'geom_any_ecp')
      if (.not. geom_any_ecp) return
      geom_any_ecp = .false.
      do icent = 1,ncenter(geom)
        if (oecpcent(icent,geom)) then
          geom_any_ecp = .true.
          return
        endif
      enddo
c      
      end
      logical function geom_ecp_center_list(geom, num_ecp_cent,
     &    ecp_cent)
      implicit none
#include "geomP.fh"
      logical geom_check_handle
      external geom_check_handle
c
      integer geom ! [input] geometry handle
      integer num_ecp_cent ! [input] dimension of ecp_cent 
*. . . . . . . . . . . . .           array from calling routine
      integer ecp_cent(num_ecp_cent)  ! [output] list of centers that
*. . . . . . . . . . . . . . . . . . .           have ECPs
*
      integer icent, num_ecp
*
      geom_ecp_center_list =
     &    geom_check_handle(geom,'geom_ecp_center_list')
      if (.not. geom_ecp_center_list) return
c
      num_ecp = 0
      do icent = 1, ncenter(geom)
        if (oecpcent(icent,geom)) then
          num_ecp = num_ecp + 1
          if (num_ecp.gt.num_ecp_cent) call errquit
     &        ('geom_ecp_center_list: number of ecp centers greater'//
     &         ' than array size passed in which is:',num_ecp_cent)
          ecp_cent(num_ecp) = icent
        endif
      enddo
      end
      logical function geom_nuc_dipole(geom,dip)
      implicit none
#include "geomP.fh"
      integer geom              ! [input]
      double precision dip(3)   ! [output] Returns the nuclear dipole in AU
c
      logical geom_check_handle
      external geom_check_handle
c
      geom_nuc_dipole = geom_check_handle(geom,'geom_nuc_dipole')
      if (.not. geom_nuc_dipole) return
c
      dip(1) = ndipole(1,geom)
      dip(2) = ndipole(2,geom)
      dip(3) = ndipole(3,geom)
c
      end
      logical function geom_calc_distance(a,b,ab)
      implicit none
c
* computes distance between two atoms
c
      double precision a(3) ! [input] coords of center a
      double precision b(3) ! [input] coords of center b
      double precision ab   ! [output] distance between centers a,b
c
      ab = (a(1)-b(1))*(a(1)-b(1))
      ab = (a(2)-b(2))*(a(2)-b(2)) + ab
      ab = (a(3)-b(3))*(a(3)-b(3)) + ab
      ab = sqrt(ab)
      geom_calc_distance = ab.ge.0.0d00
      end
c
      logical function geom_calc_angle(a,b,c,angle)
      implicit none
c
c computes the angle (in degrees) between 3 atoms in order given
c
c::-functions
      logical geom_calc_distance
      external geom_calc_distance
c::-passed
      double precision a(3) ! [input] coordinates of center a
      double precision b(3) ! [input] coordinates of center b
      double precision c(3) ! [input] coordinates of center c
      double precision angle ! [output] the angle (in degrees)
c::-local
      double precision ab, bc, ac, xcosine
      double precision pi
c
      pi = 2.0d00*acos(0.0d00)
      geom_calc_angle = geom_calc_distance(a,b,ab)
      geom_calc_angle = geom_calc_angle.and.geom_calc_distance(b,c,bc)
      geom_calc_angle = geom_calc_angle.and.geom_calc_distance(a,c,ac)
      if (.not.geom_calc_angle) call errquit
     &    ('geom_calc_angle:error computing a distance',911)
      
      xcosine = ab*ab + bc*bc - ac*ac
      xcosine = xcosine/(2.0d00*ab*bc)
   
      if( abs(xcosine) .gt. 1.00d00 ) xcosine = sign(1.0d00,xcosine)

      angle = (180.0d00/pi)*acos(xcosine)
      
      end
      logical function geom_calc_dihedral(a,b,c,d,dihedral)
      implicit none 
c
c computes the dihedral angle for the given 4 atom coordinates
c
c::-functions
      logical geom_calc_distance
      external geom_calc_distance
c::-passed      
      double precision a(3) ! [input] coordinates of center a
      double precision b(3) ! [input] coordinates of center b
      double precision c(3) ! [input] coordinates of center c
      double precision d(3) ! [input] coordinates of center d
      double precision dihedral ! [output] the dihedral angle (in degrees)
c::-local
      double precision ab,pab,qab,rab
      double precision pac,qac,rac
      double precision pbd,qbd,rbd
      double precision dot, aa, xcosine, xsine
      double precision pi
c
      pi = 2.0d00*acos(0.0d00)
      geom_calc_dihedral = geom_calc_distance(a,b,ab)
      if (.not.geom_calc_dihedral) call errquit
     &    ('geom_calc_dihedral: error computing distance',911)
      
      pab = (b(1) - a(1))/ab
      qab = (b(2) - a(2))/ab
      rab = (b(3) - a(3))/ab

      pac = (c(1) - a(1))
      qac = (c(2) - a(2))
      rac = (c(3) - a(3))
      dot = pac*pab + qac*qab + rac*rab
      pac = pac - dot*pab
      qac = qac - dot*qab
      rac = rac - dot*rab
      aa  = sqrt((pac*pac+qac*qac+rac*rac))
      pac = pac/aa
      qac = qac/aa
      rac = rac/aa

      pbd = d(1) - b(1)
      qbd = d(2) - b(2)
      rbd = d(3) - b(3)
      dot = pbd*pab + qbd*qab + rbd*rab
      pbd = pbd - dot*pab
      qbd = qbd - dot*qab
      rbd = rbd - dot*rab
      aa  = sqrt((pbd*pbd+qbd*qbd+rbd*rbd))
      pbd = pbd/aa
      qbd = qbd/aa
      rbd = rbd/aa
      
      xcosine = pac*pbd + qac*qbd + rac*rbd
      xsine   = pab*(qbd*rac - rbd*qac)
      xsine   = xsine + qab*(rbd*pac - pbd*rac)
      xsine   = xsine + rab*(pbd*qac - qbd*pac)
      
*      dihedral = 360.0d00 - ((180.0d00/pi)*atan2(xsine,xcosine))
      dihedral = ((180.0d00/pi)*atan2(xsine,xcosine))
*      if (dihedral.eq.360.0d00)dihedral = 0.0d00
      end
      logical function geom_print_distance(geom)
      implicit none
c
c prints arbitrary i>j atom distances
c
#include "geom.fh"
#include "stdio.fh"
#include "inp.fh"
c::-functions
      logical geom_calc_distance
      external geom_calc_distance
c::-passed      
      integer geom ! [input] geometry handle
c::-local
      double precision thresh
      integer nat ! number of atoms
      integer iat ! ith atom
      integer jat ! jth atom
      double precision chg   ! charge (ignored)
      double precision ci(3) ! coords of atom i
      character*16 tagi      ! tag of atom i
      double precision cj(3) ! coords of atom j
      character*16 tagj      ! tag of atom j
      integer lmtag
      double precision dij   ! distance between atoms i and j
      double precision dist_min
      double precision dist_ave
      double precision dist_max
      double precision dist_cnt 
      integer num_pos
      integer num_prt

c
      thresh = 3.00001d00
*rak:      thresh = thresh + 0.01d00
*rak:00001 continue
*rak:      thresh = thresh - 0.01d00
      if (.not.geom_ncent(geom,nat)) call errquit
     &    ('geom_print_distance: ',911)

      num_pos = nat*(nat-1)/2
      num_prt = 0
      dist_cnt = 0.0d00
      dist_min = 4040.0d00
      dist_max = 0.0d00
      dist_ave = 0.0d00
      write(luout,'(1x,a,f10.4)')
     &    ' geom_print_distance: threshold :',thresh
      do iat = 1,nat
        if (.not.geom_cent_get(geom,iat,tagi,ci,chg))call errquit
     &      ('geom_print_distance: ',911)
        lmtag = inp_strlen(tagi)
        do jat = 1,iat
          if (iat.ne.jat) then

            if (.not.geom_cent_get(geom,jat,tagj,cj,chg))call errquit
     &          ('geom_print_distance: ',911)

            if (.not.geom_calc_distance(ci,cj,dij)) call errquit
     &          ('geom_print_distance: ',911)

            dist_min = min(dist_min,dij)
            dist_max = max(dist_max,dij)
            dist_ave = dist_ave + dij
            dist_cnt = dist_cnt + 1.0d00
            if (dij.lt.thresh) then
              lmtag = max(lmtag,inp_strlen(tagj))
              write(luout,10000)
     &            tagi(1:lmtag),iat,
     &            tagj(1:lmtag),jat,dij
              num_prt = num_prt + 1
            endif
          endif
        enddo
      enddo
10000 format(1x,'distance(',a,'|',i4,',',a,'|',i4,') =',f12.6)
      dist_ave = dist_ave/dist_cnt
      write(luout,'(1x,a,f12.6)')
     &    ' minimum distance ',dist_min
      write(luout,'(1x,a,f12.6)')
     &    ' maximum distance ',dist_max
      write(luout,'(1x,a,f12.6)')
     &    ' average distance ',dist_ave
      write(luout,'(1x,a,i7)')
     &    'possible pairs to print:',num_pos
      write(luout,'(1x,a,i7)')
     &    '          pairs printed:',num_prt
      geom_print_distance = .true.
*      write(69,'(1x,f12.6,i10)')thresh,num_prt
*      if (num_prt.gt.0) goto 00001
      end
      logical function geom_print_angles(geom)
      implicit none
#include "geom.fh"
#include "inp.fh"
#include "stdio.fh"
c::-functions
      logical geom_calc_distance
      external geom_calc_distance
      logical geom_calc_angle
      external geom_calc_angle
c::-passed      
      integer geom ! [input] geometry handle
c::-local
      double precision thresh
      integer nat ! number of atoms
      integer iat ! ith atom
      integer jat ! jth atom
      integer kat ! kth atom
      double precision chg   ! charge (ignored)
      double precision ci(3) ! coords of atom i
      character*16 tagi      ! tag of atom i
      double precision cj(3) ! coords of atom j
      character*16 tagj      ! tag of atom j
      double precision ck(3) ! coords of atom k
      character*16 tagk      ! tag of atom k
      integer lmtag
      double precision dij   ! distance between atoms i and j
      double precision djk   ! distance between atoms j and k
      double precision dik   ! distance between atoms i and k
      double precision angle ! angle to be printed
      logical FF, FT         ! fortran true and false
      integer ngood          ! number of sides under threshold
      logical dij_okay       ! dij under threshold
      logical djk_okay       ! djk under threshold
      logical dik_okay       ! dik under threshold
      logical print_ijk      ! print angle i, j, k
      logical print_ikj      ! print angle i, k, j
      logical print_jik      ! print angle j, i, k
      logical should_print   ! should something be printed?
      integer num_pos
      integer num_prt
c
c initialize variables
      thresh = 3.00001d00
*      thresh = thresh + 0.01d00
*00001 continue
*      thresh = thresh - 0.01d00
      FF = .false.
      FT = .true.
      dij_okay = FF
      djk_okay = FF
      dik_okay = FF
      num_prt = 0
      
      if (.not.geom_ncent(geom,nat)) call errquit
     &    ('geom_print_angles: ',911)

      num_pos = nat*(nat-1)*(nat-2)/6

      if (nat.lt.3) return
      write(luout,'(1x,a,f10.4)')
     &    ' geom_print_angles: distance threshold :',thresh
      do iat = 1,nat
        if (.not.geom_cent_get(geom,iat,tagi,ci,chg))call errquit
     &      ('geom_print_angles: ',911)
        lmtag = inp_strlen(tagi)
        do jat = 1,iat
          if (iat.ne.jat) then

            if (.not.geom_cent_get(geom,jat,tagj,cj,chg))call errquit
     &          ('geom_print_angles: ',911)

            lmtag = max(lmtag,inp_strlen(tagj))
            if (.not.geom_calc_distance(ci,cj,dij)) call errquit
     &          ('geom_print_angles: ',911)

            dij_okay = dij.lt.thresh
            if (dij_okay) then
              do kat = 1,jat
                if (kat.ne.jat.and.kat.ne.iat) then
                  if (.not.geom_cent_get(geom,kat,tagk,ck,chg))
     &                call errquit
     &                ('geom_print_angles: ',911)
                  lmtag = max(lmtag,inp_strlen(tagk))
                  
                  if (.not.geom_calc_distance(ci,ck,dik)) call errquit
     &                ('geom_print_angles: ',911)
                  if (.not.geom_calc_distance(cj,ck,djk)) call errquit
     &                ('geom_print_angles: ',911)
                  dik_okay = dik.lt.thresh
                  djk_okay = djk.lt.thresh
                  ngood = 0
                  if (dij_okay) ngood = ngood + 1
                  if (dik_okay) ngood = ngood + 1
                  if (djk_okay) ngood = ngood + 1
*
* ngood is 0 or 1 then atoms too far apart to be interesting
*
                  print_ijk = FF ! a(ijk) = a(kji)
                  print_ikj = FF ! a(ikj) = a(jki)
                  print_jik = FF ! a(jik) = a(kji)
                  if (ngood.eq.2) then
* ngood = 2 then only one interesting angle
                    if     (dij_okay.and.dik_okay) then
                      print_jik = FT ! then angle should be j, i, k
                    elseif (dij_okay.and.djk_okay) then
                      print_ijk = FT ! then angle should be i, j, k
                    elseif (dik_okay.and.djk_okay) then
                      print_ikj = FT ! then angle should be i, k, j
                    else
                      call errquit(' should not get here 1',911)
                    endif
                  elseif (ngood.eq.3) then
                    
* if isocoles print angle between equal sides
                    if (dij.eq.djk) then
                      print_ijk = FT
                    else if (dij.eq.dik) then
                      print_jik = FT
                    else if (djk.eq.dik) then
                      print_ikj = FT
                      
* print angle with largest value.
                    else if (dij.gt.djk.and.dij.gt.dik) then
                      print_ikj = FT
                    else if (djk.gt.dij.and.djk.gt.dik) then
                      print_jik = FT
                    else if (dik.gt.dij.and.dik.gt.djk) then
                      print_ijk = FT
                    else
                      call errquit(' should not get here 2',911)
                    endif
                  endif
                  should_print = (ngood.eq.2.or.ngood.eq.3) .and.
     &                (print_ijk.or.print_ikj.or.print_jik)
                  if (print_ijk) then
                    if (.not.should_print) stop 'error'
                    if (.not.geom_calc_angle(ci,cj,ck,angle))
     &                  call errquit
     &                  (' error ',911)
                    num_prt =num_prt + 1
                    write(luout,10000)
     &                  tagi(1:lmtag),iat,
     &                  tagj(1:lmtag),jat,
     &                  tagk(1:lmtag),kat,angle
                  else if (print_ikj) then
                    if (.not.should_print) stop 'error'
                    if (.not.geom_calc_angle(ci,ck,cj,angle))
     &                  call errquit
     &                  (' error ',911)
                    num_prt =num_prt + 1
                    write(luout,10000)
     &                  tagi(1:lmtag),iat,
     &                  tagk(1:lmtag),kat,
     &                  tagj(1:lmtag),jat,angle
                  else if (print_jik) then
                    if (.not.should_print) stop 'error'
                    if (.not.geom_calc_angle(cj,ci,ck,angle))
     &                  call errquit
     &                  (' error ',911)
                    num_prt =num_prt + 1
                    write(luout,10000)
     &                  tagj(1:lmtag),jat,
     &                  tagi(1:lmtag),iat,
     &                  tagk(1:lmtag),kat,angle
                  endif
                endif
              enddo
            endif
          endif
        enddo
      enddo
10000 format(1x,'angle(',a,'|',i4,',',a,'|',i4,',',a,'|',i4,
     &    ') =',f9.3)
      write(luout,'(1x,a,i7)')
     &    'possible angles to print:',num_pos
      write(luout,'(1x,a,i7)')
     &    '         angles  printed:',num_prt
      geom_print_angles = FT
*      write(69,'(1x,f12.6,i10)')thresh,num_prt
*      if (num_prt.gt.0) goto 00001
      end
      logical function geom_print_dihedrals(geom)
      implicit none
#include "mafdecls.fh"
#include "geom.fh"
#include "stdio.fh"
#include "inp.fh"
c::-functions
      logical geom_calc_distance
      external geom_calc_distance
      logical geom_calc_dihedral
      external geom_calc_dihedral
c::-passed      
      integer geom ! [input] geometry handle
c::-local
      double precision thresh
      integer nat ! number of atoms
      integer iat ! ith atom
      integer jat ! jth atom
      integer kat ! kth atom
      integer lat ! lth atom
      double precision chg   ! charge (ignored)
      double precision ci(3) ! coords of atom i
      character*16 tagi      ! tag of atom i
      double precision cj(3) ! coords of atom j
      character*16 tagj      ! tag of atom j
      double precision ck(3) ! coords of atom k
      character*16 tagk      ! tag of atom k
      double precision cl(3) ! coords of atom k
      character*16 tagl      ! tag of atom k
      integer lmtag
      double precision dij   ! distance between atoms i and j
      double precision dik   ! distance between atoms i and k
      double precision dil   ! distance between atoms i and l
      double precision djk   ! distance between atoms j and k
      double precision djl   ! distance between atoms j and l
      double precision dkl   ! distance between atoms k and l
      double precision diangle ! dihedral angle to be printed
      logical FF, FT         ! fortran true and false
      logical dij_okay       ! dij under threshold
      logical dik_okay       ! dik under threshold
      logical dil_okay       ! dil under threshold
      logical djk_okay       ! djk under threshold
      logical djl_okay       ! djl under threshold
      logical dkl_okay       ! dkl under threshold
      logical should_print   ! should something be printed?
c
      integer k_printed      ! ma index to printed dihedral angles
      integer h_printed      ! ma handle to printed dihedral angles
      integer dh_printed
      integer nat_size
      integer indx
      integer ip
      logical will_print
      integer num_pos
      integer num_prt
      logical print4, print8, print12, print16
c::-statement functions
      integer si, sj, sk, sl
      integer isym2m, isym2, isym4mm, isym4m, isym4
      isym2m(si,sj)=max(si,sj)*((max(si,sj))-1)/2 + min(si,sj)
      isym2(si,sj) = si*(si-1)/2 + sj
      isym4mm(si,sj,sk,sl) =
     &    max(isym2m(si,sj),isym2m(sk,sl)) *
     &    (max(isym2m(si,sj),isym2m(sk,sl))-1)/2 +
     &    min(isym2m(si,sj),isym2m(sk,sl))
      isym4m(si,sj,sk,sl) = 
     &    max(isym2(si,sj),isym2(sk,sl)) *
     &    (max(isym2(si,sj),isym2(sk,sl))-1)/2 +
     &    min(isym2(si,sj),isym2(sk,sl))
      isym4(si,sj,sk,sl) = 
     &    isym2(si,sj)*(isym2(si,sj)-1)/2 +
     &    isym2(sk,sl)
c
      
      if (.not.geom_ncent(geom,nat)) call errquit
     &    ('geom_print_dihedrals: ',911)

      num_pos = nat*(nat-1)*(nat-2)*(nat-3)/24

      FF = .false.
      FT = .true.

      if (nat.lt.4) return

      thresh = 4.00001d00
*      thresh = thresh+0.01d00
*00001 continue
*      thresh = thresh-0.01d00
c initialize variables
      dij_okay = .false.      ! dij under threshold
      dik_okay = .false.      ! dik under threshold
      dil_okay = .false.      ! dil under threshold
      djk_okay = .false.      ! djk under threshold
      djl_okay = .false.      ! djl under threshold
      dkl_okay = .false.      ! dkl under threshold
      num_prt = 0
      write(luout,'(1x,a,f10.4)')
     &    ' geom_print_dihedrals: distance threshold :',thresh
c
*.. get ma array for printed dihedrals
      dh_printed = 0
*      nat_size = isym4mm(nat,nat,nat,1)
      nat_size = 1000
      if (.not.ma_push_get(mt_int,nat_size,'print list for dihedrals',
     &    h_printed, k_printed)) call errquit
     &    ('geom_print_dihedrals: not enough stack inc. needed',
     &    nat_size)
c
      lmtag = -1
      do iat = 1,nat
        if (.not.geom_cent_get(geom,iat,tagi,ci,chg))call errquit
     &      ('geom_print_dihedrals: ',911)
        lmtag = max(lmtag,inp_strlen(tagi))
      enddo
      print4  = FF
      print8  = FF
      print12 = FF
      print16 = FF
      if (lmtag.le.4) then
        print4 = FT
      elseif (lmtag.le.8) then
        print8 = FT
      elseif (lmtag.le.12) then
        print12 = FT
      else
        print16 = FT
      endif
      if (print4)  write(luout,10001)
      if (print8)  write(luout,10003)
      if (print12) write(luout,10005)
      if (print16) write(luout,10007)
      do iat = 1,nat
        if (.not.geom_cent_get(geom,iat,tagi,ci,chg))call errquit
     &      ('geom_print_dihedrals: ',911)
        do jat = 1,iat
          if (iat.ne.jat) then
            
            dh_printed = 0
            if (.not.geom_cent_get(geom,jat,tagj,cj,chg))call errquit
     &          ('geom_print_dihedrals: ',911)
            
            if (.not.geom_calc_distance(ci,cj,dij)) call errquit
     &          ('geom_print_dihedrals: ',911)
            
            dij_okay = dij.lt.thresh
            if (dij_okay) then
              do kat = 1,jat
                if (kat.ne.jat.and.kat.ne.iat) then
                  if (.not.geom_cent_get(geom,kat,tagk,ck,chg))
     &                call errquit
     &                ('geom_print_dihedrals: ',911)
                  
                  if (.not.geom_calc_distance(ci,ck,dik)) call errquit
     &                ('geom_print_dihedrals: ',911)
                  if (.not.geom_calc_distance(cj,ck,djk)) call errquit
     &                ('geom_print_dihedrals: ',911)
                  dik_okay = dik.lt.thresh
                  djk_okay = djk.lt.thresh
                  if (djk_okay) then
                    do lat = 1,kat
                      if(lat.ne.iat.and.lat.ne.jat.and.
     &                    lat.ne.kat) then
                        if (.not.geom_cent_get(geom,lat,tagl,cl,chg))
     &                      call errquit
     &                      ('geom_print_dihedrals: ',911)
                        
                        if (.not.geom_calc_distance(ci,cl,dil))
     &                      call errquit
     &                      ('geom_print_dihedrals: ',911)
                        if (.not.geom_calc_distance(cj,cl,djl))
     &                      call errquit
     &                      ('geom_print_dihedrals: ',911)
                        if (.not.geom_calc_distance(ck,cl,dkl))
     &                      call errquit
     &                      ('geom_print_dihedrals: ',911)
                        dil_okay = dil.lt.thresh
                        djl_okay = djl.lt.thresh
                        dkl_okay = dkl.lt.thresh
                        should_print =
     &                      dij_okay.and.djk_okay.and.dkl_okay
*                        should_print = should_print.and.dik_okay
*                        should_print = should_print.and.djl_okay
                        if (should_print) then
                          indx = isym4mm(iat,jat,kat,lat) 
                          will_print = FT
                          do ip = 0,(dh_printed-1)
                            if (indx.eq.int_mb(k_printed+ip))
     &                          will_print = FF
                          enddo
                          if (will_print) then
                            int_mb(k_printed+dh_printed) = indx
                            dh_printed = dh_printed + 1
                            if ((dh_printed+1).gt.nat_size)
     &                          stop ' error size'
                            if (.not.geom_calc_dihedral
     &                          (ci,cj,ck,cl,diangle))
     &                          call errquit('dih error',911)
                            num_prt = num_prt + 1
                            if (print4) write(luout,10002)
     &                          tagi(1:lmtag),iat,
     &                          tagj(1:lmtag),jat,
     &                          tagk(1:lmtag),kat,
     &                          tagl(1:lmtag),lat,diangle
                            if (print8) write(luout,10004)
     &                          tagi(1:lmtag),iat,
     &                          tagj(1:lmtag),jat,
     &                          tagk(1:lmtag),kat,
     &                          tagl(1:lmtag),lat,diangle
                            if (print12) write(luout,10006)
     &                          tagi(1:lmtag),iat,
     &                          tagj(1:lmtag),jat,
     &                          tagk(1:lmtag),kat,
     &                          tagl(1:lmtag),lat,diangle
                            if (print16) write(luout,10008)
     &                          tagi(1:lmtag),iat,
     &                          tagj(1:lmtag),jat,
     &                          tagk(1:lmtag),kat,
     &                          tagl(1:lmtag),lat,diangle
*                        write(luout,10000)
*     &                      tagi(1:lmtag),iat,
*     &                      tagj(1:lmtag),jat,
*     &                      tagk(1:lmtag),kat,
*     &                      tagl(1:lmtag),lat,diangle
                          endif
                        endif
                      endif
                    enddo
                  endif
                endif
              enddo
            endif
          endif
        enddo
      enddo
10000 format(1x,'dihedral angle(',
     &    a,'|',i5,',',a,'|',i5,',',a,'|',i5,',',a,'|',i5,
     &    ') =',f9.3)
10001 format(1x,33('-'),'dihedral angles',32('-'),/,
     &    2x,4(1x,'tag  number  '),3x,'dihedral angle',/,
     &    1x,80('-'),/)
10002 format(2x,4(a4,3x,i4,3x),5x,f9.3)
10003 format(1x,33('-'),'dihedral angles',32('-'),/,
*                 123456789012345 
     &    2x,4(1x,'  tag   number'),3x,' dihedral angle ',/,
     &    1x,80('-'),/)
10004 format(2x,4(a8,2x,i4,1x),6x,f9.3)
10005 format(1x,25('-'),'dihedral angles',25('-'),/,
     &    4(1x,'atom tag, number'),3x,' dihedral angle ',/,
     &    80('-'),/)
10007 format(1x,25('-'),'dihedral angles',25('-'),/,
     &    4(1x,'atom tag, number'),3x,' dihedral angle ',/,
     &    80('-'),/)
10006 format(4(2x,a12,i4),f9.3)
10008 format(4(2x,a16,i4),f9.3)
      write(luout,'(1x,a,i10)')
     &    'possible dihedral angles to print:',num_pos
      write(luout,'(1x,a,i10)')
     &    '          dihedral angles printed:',num_prt
*      write(69,'(1x,f12.6,i10)')thresh,num_prt
      geom_print_dihedrals = ma_pop_stack(h_printed)
*      if (num_prt.gt.0) goto 00001
      end
      

