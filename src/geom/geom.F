      block data geom_data
C$Id: geom.F,v 1.40 1995-11-19 19:32:19 d3g681 Exp $
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
      data angstrom_to_au /1.8897265d0/
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
      double precision e, r
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

               r = dsqrt(
     $              (coords(1,i,geom)-coords(1,j,geom))**2 + 
     $              (coords(2,i,geom)-coords(2,j,geom))**2 + 
     $              (coords(3,i,geom)-coords(3,j,geom))**2)
               e = e + charge(i,geom)*charge(j,geom)/r
            endif
         end do
      end do
c
c     add in interaction of nuclear dipole with external field
c
      e = e + ddot(3, ndipole(1,geom), 1, efield(1,geom), 1)
c
      erep(geom) = e
c
      if (isystype(geom) .eq. 3) then
         call geom_3d(geom)
      endif
c
      call sym_init_inv_op(geom)
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
      tmp(k:) = ':efield'
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
      user_units(geom) = 'a.u.'
      include_bqbq(geom) = .false.
c
      do i = 1, 3
        lattice_vectors(i,geom) = 0
        lattice_angles(i,geom) = 0
        do j = 1, 3
           amatrix(j,i,geom) = 0.0d0
           amatrix_inv(j,i,geom) = 0.0d0
        end do
        amatrix(i,i,geom) = 1.0d0
        amatrix_inv(i,i,geom) = 1.0d0
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
      end do
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
c     Basic printing of cartesian geometry ... needs support for
c     user defined units, internal coords, different formats, ...
c
      integer geom              ! [input]
      integer icent
      logical geom_check_handle, geom_check_cent
      external geom_check_handle, geom_check_cent
      integer i, j
      double precision scale, tmp(3)
      character*80 buf
      character*10 csys
      logical oprint_uniq
c
      if (.not. geom_check_handle(geom, 'geom_print')) then
         geom_print = .false.
         return
      end if
c
      if (user_units(geom) .eq. 'a.u.') then
         scale = 1.0d0
      else if (user_units(geom) .eq. 'fractional') then
         scale = 1.0d0
      else if (user_units(geom) .eq. 'angstrom') then
         scale = 1.0d0 / angstrom_to_au
      else
         call errquit('geom_print: user_units are unspecified', 0)
      end if
c
      if (isystype(geom) .eq. 0) then
         csys = 'Molecule'
      else if (isystype(geom) .eq. 1) then
         csys = 'Polymer'
      else if (isystype(geom) .eq. 2) then
         csys = 'Slab'
      else if (isystype(geom) .eq. 3) then
         csys = 'Crystal'
      else
         call errquit('geom_print: unknown system type', isystype(geom))
      end if
c
      buf = ' '
      write(buf,1) 
     $     csys(1:inp_strlen(csys)),
     $     names(geom)(1:lenn(geom)), 
     $     trans(geom)(1:lent(geom))
 1    format(a,' "',a,'" -> "',a,'"')
      write(6,*)
      write(6,*)
      call util_print_centered(6,buf,40,.true.)
      write(6,*)
      write(6,2) user_units(geom)(1:inp_strlen(user_units(geom)))
 2    format(' Units of output coordinates = ', a)
      if (include_bqbq(geom))
     $     write(6,*) ' Include Bq-Bq interactions'
      if (isystype(geom) .eq. 3) then
         write(6,21) 
     $        volume_direct(geom),
     $        (lattice_vectors(i,geom),i=1,3),
     $        (lattice_angles(i,geom),i=1,3),
     $        (recip_lat_vectors(i,geom),i=1,3),
     $        (recip_lat_angles(i,geom),i=1,3)
 21      format(
     $        ' Direct space cell volume (angstrom^3)   ', f12.6/
     $        ' Lattice vectors (angstrom)              ', 3f12.6/
     $        ' Lattice angles (degrees)                ', 3f12.6/
     $        ' Reciprocal lattice vectors (angs^-1)    ', 3f12.6/
     $        ' Reciprocal lattice angles (degrees)     ', 3f12.6/)
c
         write(6,*) ' Crystallographic to cartesian (angstrom) matrix'
         write(6,*)
         write(6,22) ((amatrix(i,j,geom)/angstrom_to_au,j=1,3),i=1,3)
 22      format(3(5x,3f12.6/))
c
         write(6,*) ' Cartesian (angstrom) to crystallographic matrix'
         write(6,*)
         write(6,22)((amatrix_inv(i,j,geom)*angstrom_to_au,j=1,3),i=1,3)
c
         write(6,*) ' Metric matrix (angstrom^2) for the direct lattice'
         write(6,*)
         write(6,22)((metric_matrix(i,j,geom),j=1,3),
     $        i=1,3)
c
      endif
      write(6,*)
c
      write(6,3) 
 3    format('  No.       Tag          Charge          X',
     $     '              Y              Z'/
     $       ' ---- ---------------- ---------- --------------',
     $     ' -------------- --------------')
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
            tmp(i) = tmp(i)*scale ! Scale units as necessary
         end do
         write(6,4) icent, tags(icent,geom), charge(icent,geom),
     $        (tmp(i),i=1,3)
 4       format(' ',i4,' ',a16,' ',f10.6,3f15.8)
      end do
      write(6,*)
c
      if (isystype(geom) .eq. 0) then
         write(6,41) erep(geom)
 41      format(/' Effective nuclear repulsion energy (au) ', f18.10/)
c
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
      end if
c
      if (oefield(geom)) then
         write(6,9)
 9       format('         Electric Field (au) ')
         write(6,10)
 10      format('         -------------------')
         write(6,7)
         write(6,8) (efield(i,geom), i=1,3)
         write(6,*)
      end if
c
      oprint_uniq = sym_num_ops(geom) .gt. 0
      call sym_print_all(geom, .true., oprint_uniq, .false., .false.,
     $     .false.)
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
      character*(*) tag         ! [input]
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
