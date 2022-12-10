c
C> \ingroup geom
C> @{
      block data geom_data
C$Id$
      implicit none
#include "nwc_const.fh"
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
c--> names of the extra 3-dimensional space groups
c
c     *** extra triclinic groups ***
      data (sym_spgnames(i),i=231,240)/
     & 'A1', 'B1', 'C1', 'F1', 'I1',
     & 'A-1','B-1','C-1','F-1','I-1'/
c
c     *** extra monoclinic groups - 45 groups - y-axis unique***
c      data (sym_spgnames(i),i=241,285)/
c     & 'P121', 'B121', 
c     & 'P12_11', 'B12_11',
c     & 'C121','A121',I121','F121',
c     & 'P1m1','B1m1',
c     & 'P1c1','P1a1','P1n1','B1a1','B1d1',
c     & 'C1m1','A1m1','I1m1','F1m1',
c     & 'C1c1','A1a1','I1a1','F1d1'
c     & 'P12/m1','B12/m1',
c     & 'P12_1/m1','B12_1/m1',
c     & 'C12/m1','A12/m1','I12/m1','F12/m1',
c     & 'P12/c1','P12/a1','P12/n1','B12/a1','B12/d1',
c     & 'P12_1/c1','P12_1/a1','P12_1/n1','B12_1/a1','B12_1/d1',
c     & 'C12/c1','A12/a1','I12/a1','F12/d1'/
c
c     *** extra monoclinic groups - 45 groups - z-axis unique ***
c      data (sym_spgnames(i),i=286,330)/
c     & 'P112','C112',
c     & 'P112_1','C112_1',
c     & 'A112','B112','I112','F112',
c     & 'P11m','C11m',
c     & 'P11a','P11b','P11n','C11a','C11d',
c     & 'A11m','B11m','I11m','F11m',
c     & 'A11a','B11b','I11a','F11d',
c     & 'P112/m','C112/m'
c     & 'P112_1/m','C112_1/m'
c     & 'A112/m','B112/m','I112/m','F112/m',
c     & 'P112/a','P112/b','P112/n','C112/a','C112/d'
c     & 'P112_1/a','P112_1/b','P112_1/n','C112_1/a','C112_1/d',
c     & 'A112/a','B112/b','I112/a','F112/d'/
c
c     *** extra monoclinic groups - 45 groups - x-axis unique ***
c      data (sym_spgnames(i),i=331,375)/
c     & 'P211','A211',
c     & 'P2_111','A2_111',
c     & 'B211','C211','I211','F211',
c     & 'Pm11','Am11',
c     & 'Pb11','Pc11','Pn11','Ab11','Ad11',
c     & 'Bm11','Cm11','Im11','Fm11',
c     & 'Bb11','Cc11','Ib11','Fd11',
c     & 'P2/m11','A2/m11',
c     & 'P2_1/m11','A2_1/m11',
c     & 'B2/m11','C2/m11','I2/m11','F2/m11',
c     & 'P2/b11','P2/c11','P2/n11','A2/b11','A2/d11',
c     & 'P2_1/b11','P2_1/c11','P2_1/n11','A2_1/b11','A2_1/d11',
c     & 'B2_1/b11','C2/c11','I2/b11','F2/d11'/
c
c    *** extra orthorhombic groups ***
c      data (sym_spgnames(i),i=???,???)/
c     & /
c
c    *** extra tetragonal groups ***
c      data (sym_spgnames(i),i=???,???)/
c     & /
c
c    *** extra trigonal groups ***
c      data (sym_spgnames(i),i=???,???)/
c     & /

c--> names of the 3-dimensional space groups, without _ for car files
c
      data (sym_carnames(i),i=1,95) /
     & 'P1','P-1','P2','P21','C2',
     & 'Pm','Pc','Cm','Cc','P2/m',
     & 'P21/m','C2/m','P2/c','P21/c','C2/c',
     & 'P222','P2221','P21212','P212121','C2221',
     & 'C222','F222','I222','I212121','Pmm2',
     & 'Pmc21','Pcc2','Pma2','Pca21','Pnc2',
     & 'Pmn21','Pba2','Pna21','Pnn2','Cmm2',
     & 'Cmc21','Ccc2','Amm2','Abm2','Ama2',
     & 'Aba2','Fmm2','Fdd2','Imm2','Iba2',
     & 'Ima2','Pmmm','Pnnn','Pccm','Pban',
     & 'Pmma','Pnna','Pmna','Pcca','Pbam',
     & 'Pccn','Pbcm','Pnnm','Pmmn','Pbcn',
     & 'Pbca','Pnma','Cmcm','Cmca','Cmmm',
     & 'Cccm','Cmma','Ccca','Fmmm','Fddd',
     & 'Immm','Ibam','Ibca','Imma','P4',
     & 'P41','P42','P43','I4','I41',
     & 'P-4','I-4','P4/m','P42/m','P4/n',
     & 'P42/n','I4/m','I41/a','P422','P4212',
     & 'P4122','P41212','P4222','P42212','P4322'/
      data (sym_carnames(i),i=96,190)/   
     & 'P43212','I422','I4122','P4mm','P4bm',
     & 'P42cm','P42nm','P4cc','P4nc','P42mc',
     & 'P42bc','I4mm','I4cm','I41md','I41cd',
     & 'P-42m','P-42c','P-421m','P-421c','P-4m2',
     & 'P-4c2','P-4b2','P-4n2','I-4m2','I-4c2',
     & 'I-42m','I-42d','P4/mmm','P4/mcc','P4/nbm',
     & 'P4/nnc','P4/mbm','P4/mnc','P4/nmm','P4/ncc',
     & 'P42/mmc','P42/mcm','P42/nbc','P42/nnm','P42/mbc',
     & 'P42/mnm','P42/nmc','P42/ncm','I4/mmm','I4/mcm',
     & 'I41/amd','I41/acd','P3','P31','P32',
     & 'R3','P-3','R-3','P312','P321',
     & 'P3112','P3121','P3212','P3221','R32',
     & 'P3m1','P31m','P3c1','P31c','R3m',
     & 'R3c','P-31m','P-31c','P-3m1','P-3c1',
     & 'R-3m','R-3c','P6','P61','P65',
     & 'P62','P64','P63','P-6','P6/m',
     & 'P63/m','P622','P6122','P6522','P6222',
     & 'P6422','P6322','P6mm','P6cc','P63cm',
     & 'P63mc','P-6m2','P-6c2','P-62m','P-62c'/
      data (sym_carnames(i),i=191,230)/
     & 'P6/mmm','P6/mcc','P63/mcm','P63/mmc','P23',
     & 'F23','I23','P213','I213','Pm-3',
     & 'Pn-3','Fm-3','Fd-3','Im-3','Pa-3',
     & 'Ia-3','P432','P4232','F432','F4132',
     & 'I432','P4332','P4132','I4132','P-43m',
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
c
C> \brief Check whether a given handle corresponds to a valid geometry instance
c
C> Checks whether a given handle corresponds to a valid and active geometry. If
C> not a message will be printed on the output.
c
C> \return Returns .true. is the handle corresponds to a valid geometry, 
C> and .false. otherwise.
c
      logical function geom_check_handle(geom, msg)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
#include "stdio.fh"
c     
      integer geom              !< [Input] the geometry handle
      character*(*) msg         !< [Input] the message to be included in the error
                                !< error message
c     
      geom_check_handle = geom.gt.0 .and. geom.le.max_geom
      if (geom_check_handle) geom_check_handle = geom_check_handle
     $     .and. active(geom)
c     
      if (.not. geom_check_handle) then
         write(LuOut,*) msg,': geometry handle invalid ', geom
         call geom_err_info(msg)
      end if
c     
      end
c
C> \brief Check whether center rank is valid
c
C> Tests whether a given center rank is a valid rank for a given
C> geometry instance. If it is not an error message is generated
C> on standard output.
c
C> \return Return .true. if the rank is valid, and .false.
C> otherwise.
      logical function geom_check_cent(geom, msg, icent)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
#include "stdio.fh"
c     
      integer geom      !< [Input] the geometry handle
      character*(*) msg !< [Input] the message string
      integer icent     !< [Input] the center rank
      logical status, geom_print
      external geom_print
c     
      geom_check_cent = icent.gt.0 .and. icent.le.ncenter(geom)
      if (.not. geom_check_cent) then
         write(LuOut,*) msg,': icent invalid ', icent,
     $        names(geom)(1:lenn(geom))
         call geom_err_info(msg)
         status = geom_print(geom)
      end if
c     
      end
c
C> \brief Prints summary information about every geometry in the RTDB
c
C> This routine extracts information about all geometries from the RTDB.
C> The information is summarized and printed on standard output.
c
      subroutine geom_print_known_geoms(rtdb)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "global.fh"
#include "stdio.fh"
c     
      integer rtdb              !< [Input] the RTDB handle
c     
      integer geom, ma_type, natom, nelem
      character*26 date
      character*32 name32
      logical geom_rtdb_in, ignore
      character*128 key
c
      ignore = geom_rtdb_in(rtdb)
      if (ga_nodeid() .eq. 0) then
         write(LuOut,*)
         call util_print_centered(LuOut,'Geometries in the database',
     $        23,.true.)
         write(LuOut,*)
         if (ngeom_rtdb .le. 0) then
            write(LuOut,*) ' There are no geometries in the database'
            write(LuOut,*)
         else
            if (ngeom_rtdb .gt. 0) write(LuOut,3)
 3          format(
     $           1x,4x,2x,'Name',28x,2x,'Natoms',2x,
     $           'Last Modified',/,
     $           1x,4x,2x,32('-'),2x,6('-'),2x,24('-'))
            do geom = 1, ngeom_rtdb
               key = ' '
               write(key,'(''geometry:'',a,'':ncenter'')')
     $              names_rtdb(geom)(1:lenr(geom))
               if (.not. rtdb_get(rtdb, key, mt_int, 1, natom)) then
                  write(LuOut,*) ' Warning: geometry ', geom, 
     $                 ' may be corrupt'
                  natom = -1
               endif
               if (.not. rtdb_get_info(rtdb, key, ma_type, 
     $              nelem, date)) then
                  write(LuOut,*) ' Warning: geometry ', geom, 
     $                 ' may be corrupt'
                  date = 'unknown'
               endif
               name32 = names_rtdb(geom)(1:lenr(geom))
               write(LuOut,4) geom, name32, natom, date
 4             format(1x,i4,2x,a32,2x,i6,2x,a26)
            end do
            if (ngeom_rtdb .gt. 0) then
               if (.not. rtdb_cget(rtdb,'geometry',1,key)) 
     $              key = 'geometry'
               write(LuOut,*)
               write(LuOut,5) key(1:inp_strlen(key))
 5             format(2x,'The geometry named "',a,
     $              '" is the default for restart')
            endif
            write(LuOut,*)
            write(LuOut,*)
         endif
         call util_flush(LuOut)
      endif
c     
      end
      logical function geom_rtdb_in(rtdb)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "stdio.fh"
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
         if (ngeom_rtdb . gt. 0) then
            if (.not. rtdb_cget(rtdb, 'geometry:names', ngeom_rtdb,
     $           names_rtdb)) then
               write(LuOut,*) 'geom_rtdb_in: rtdb corrupt'
            else
               do geom = 1, ngeom_rtdb
                  lenr(geom) = inp_strlen(names_rtdb(geom))
               end do
               geom_rtdb_in = .true.
            end if
         end if
      end if
c     
      end
      logical function geom_rtdb_out(rtdb)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "stdio.fh"
c     
      integer rtdb              ! [input]
c     
c     output to rtdb info about known geometries
c     
      geom_rtdb_out  =
     $     rtdb_put(rtdb, 'geometry:ngeom', mt_int, 1, ngeom_rtdb)
      if (ngeom_rtdb . gt. 0) then
         geom_rtdb_out  =  geom_rtdb_out  .and.
     $        rtdb_cput(rtdb, 'geometry:names', ngeom_rtdb, names_rtdb)
      endif
      if (.not. geom_rtdb_out) 
     $     write(LuOut,*) ' geom_rtdb_out: rtdb is corrupt '
c     
      end
      logical function geom_rtdb_add(rtdb, name)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "geomP.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "stdio.fh"
c     
      integer rtdb              ! [input]
      character*(*) name        ! [input]
      integer geom
      logical status
      integer ln
      logical geom_rtdb_in, geom_rtdb_out
      external geom_rtdb_in, geom_rtdb_out
c
      if (ngeom_rtdb.lt.0 .or. ngeom_rtdb.gt.max_geom_rtdb)
     $     call errquit('geom_rtdb_add: ngeom_rtdb?',ngeom_rtdb,
     &       RTDB_ERR)
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
         write(LuOut,*) ' geom_rtdb_add: too many geometries on rtdb ',
     &                    name
         geom_rtdb_add = .false.
         return
      end if
      ngeom_rtdb = ngeom_rtdb + 1
      names_rtdb(ngeom_rtdb) = name
      lenr(ngeom_rtdb) = ln
c     
      if (.not. geom_rtdb_out(rtdb)) then
         write(LuOut,*) ' geom_rtdb_add: rtdb error adding ', name(1:ln)
         geom_rtdb_add = .false.
         return
      end if
c     
      geom_rtdb_add = .true.
c     
      end
      subroutine geom_err_info(info)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
#include "stdio.fh"
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
      write(LuOut,1) info, ngeom
 1    format(' ',a,': open geometies: ',i2)
      ngeom = 0
      do geom = 1, max_geom
         if (active(geom)) then
            write(LuOut,2) geom, info, names(geom)(1:lenn(geom)),
     $           trans(geom)(1:lent(geom))
 2          format(' ',i2,' ',a,': "',a, '" -> "', a,'"')
         end if
      end do
      if (ngeom_rtdb .gt. 0) then
         write(LuOut,3) info, ngeom_rtdb
 3       format(' ',a,': geometries in last accessed data base: ', i2)
         do geom = 1, ngeom_rtdb
            write(LuOut,4) names_rtdb(geom)(1:lenr(geom))
 4          format(' ',a)
         end do
      end if
c     
      end
c
C> \brief Extract the number of centers directly from the RTDB
c
C> For a named geometry this function extracts the number of centers
C> directly from the RTDB. I.e. this routine bypasses all of the usual
C> geometry infra-structure and directly exploits the stored data 
C> format.
c
C> \return Return .true. if the number of centers was found
C> successfully, and .false. otherwise.
c
      logical function geom_rtdb_ncent(rtdb, name, ncent)
      implicit none
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
      integer rtdb              !< [Input] the RTDB handle
      character*(*) name        !< [Input] the geometry name
      integer ncent             !< [Output] the number of centers
c
c     Return the number of atoms in a geometry that is
c     stored on the database ... a convenience routine.
c
      character*128 trans, tmp
      integer lent
c
      if (.not. rtdb_cget(rtdb, name, 1, trans)) trans = name
      lent = inp_strlen(trans)
      tmp = 'geometry:'//trans(1:lent)
      lent = inp_strlen(tmp)
      tmp(lent+1:) = ':ncenter'
      geom_rtdb_ncent = rtdb_get(rtdb, tmp, mt_int, 1, ncent)
c
      end
c
C> \brief Load a geometry from the RTDB
C>
C> Load a geometry from the RTDB with a specified name.
C> The name is used also to define the name of the geometry.
C>
C> If no part of a geometry is found the code assumes that the geometry
C> simply is not stored on the RTDB. If some part of a geometry is 
C> present on the RTDB but not all parts then the code assumes the
C> RTDB has been corrupted and an additional message to that effect is
C> printed on standard output.
C>
C> If loading storing the geometry was successful return .true.,
C> return .false. otherwise.
C>
      logical function geom_rtdb_load(rtdb, geom, name)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "nwc_const.fh"
#include "geomP.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
#include "stdio.fh"
*     integer node
c     
      integer rtdb              !< [Input] the RTDB handle
      integer geom              !< [Input] the geometry handle
      character*(*) name        !< [Input] the name of the geometry
c     
      double precision scale
      character*256 tmp
      integer k, nelem, ma_type
      logical s
      logical geom_check_handle, geom_rtdb_in, geom_get_user_scale
      external geom_check_handle, geom_rtdb_in, geom_get_user_scale
      logical getsym
      double precision ddot
      external ddot
c     
      geom_rtdb_load = geom_check_handle(geom, 'geom_rtdb_load')
      if (.not. geom_rtdb_load) return
      s = geom_rtdb_in(rtdb)
c     
c     translate the provided name
c     
      names(geom) = name
      lenn(geom) = inp_strlen(name)
      trans(geom) = 'junk'
      if (.not. rtdb_cget(rtdb, name, 1, trans(geom)))
     $     trans(geom) = name
*     if (.not.context_rtdb_match(rtdb, name, trans(geom)))
*     $     trans(geom) = name
      lent(geom) = inp_strlen(trans(geom))
c     
c     now get the info from the data base
c     
      tmp = 'geometry:'//trans(geom)(1:lent(geom))
      k = inp_strlen(tmp)+1
      s = .true.
c     
      tmp(k:) = ' ' 
      tmp(k:) = ':ncenter'
      s = s .and. rtdb_get(rtdb, tmp, mt_int, 1, ncenter(geom))
      geom_rtdb_load = s
      if (.not.s) then
c
c       Even the first item is not present on the RTDB so this geometry
c       is simply not present. Return this information to the caller.
c
        return
      else
c
c       At least some information about this geometry is stored on the
c       RTDB. So a complete geometry specification must be found from
c       hereon, otherwise the RTDB is corrupt.
c
      endif
      tmp(k:) = ' ' 
      tmp(k:) = ':coords'
      s = s .and. rtdb_get(rtdb, tmp, mt_dbl, max_cent*3, 
     $     coords(1,1,geom))
      tmp(k:) = ' ' 
      tmp(k:) = ':vel'
      s = s .and. rtdb_get(rtdb, tmp, mt_dbl, max_cent*3, 
     $     velocities(1,1,geom))
      tmp(k:) = ' ' 
      tmp(k:) = ':charges'
      s = s .and. rtdb_get(rtdb, tmp, mt_dbl, max_cent, charge(1,geom))
      tmp(k:) = ' ' 
      tmp(k:) = ':masses'
      s = s .and. rtdb_get(rtdb, tmp, mt_dbl, max_cent, 
     $     geom_mass(1,geom))
      tmp(k:) = ' ' 
      tmp(k:) = ':atomct'
      s = s .and. rtdb_get(rtdb, tmp, mt_dbl, max_cent, 
     $     geom_atomct(1,geom))
C     new
      tmp(k:) = ' ' 
      tmp(k:) = ':inv nuc expon'
      s = s .and. rtdb_get(rtdb, tmp, mt_dbl, ncenter(geom),
     $     geom_invnucexp(1,geom))
C     end
      tmp(k:) = ' ' 
      tmp(k:) = ':efield'
      s = s .and. rtdb_get(rtdb, tmp, mt_dbl, 3, efield(1,geom))

c     tmp(k:) = ' ' 
c     tmp(k:) = ':lattice vectors'
c     s = s .and. rtdb_get(rtdb, tmp, mt_dbl, 3, 
c    $     lattice_vectors(1,geom))
c     tmp(k:) = ' ' 
c     tmp(k:) = ':lattice angles'
c     s = s .and. rtdb_get(rtdb, tmp, mt_dbl, 3, 
c    $     lattice_angles(1,geom))
      tmp(k:) = ' ' 
      tmp(k:) = ':amatrix'
      s = s .and. rtdb_get(rtdb, tmp, mt_dbl, 9, 
     $     amatrix(1,1,geom))

      tmp(k:) = ' ' 
      tmp(k:) = ':system type'
      s = s .and. rtdb_get(rtdb, tmp, mt_int, 1, isystype(geom))
      tmp(k:) = ' ' 
      tmp(k:) = ':no. unique centers'
      s = s .and. rtdb_get(rtdb, tmp, mt_int, 1, ncenter_unique(geom))
      tmp(k:) = ' ' 
      tmp(k:) = ':group number'
      s = s .and. rtdb_get(rtdb, tmp, mt_int, 1, group_number(geom))
      tmp(k:) = ' ' 
      tmp(k:) = ':group name'
      s = s .and. rtdb_cget(rtdb, tmp, 1, group_name(geom))
      tmp(k:) = ' '
      tmp(k:) = ':use_primitive'
      s = s .and. rtdb_get(rtdb, tmp, mt_log, 1, use_primitive(geom))
      tmp(k:) = ' '
      tmp(k:) = ':primitive_center'
      s = s .and. rtdb_cget(rtdb, tmp, 1, primitive_center(geom))
      tmp(k:) = ' ' 
      tmp(k:) = ':user units'
      s = s .and. rtdb_cget(rtdb, tmp, 1, user_units(geom))
      tmp(k:) = ' ' 
      tmp(k:) = ':angstrom_to_au'
      s = s .and. rtdb_get(rtdb, tmp, mt_dbl, 1, angstrom_to_au)
      tmp(k:) = ' ' 
      tmp(k:) = ':setting number'
      s = s .and. rtdb_get(rtdb, tmp, mt_int, 1, setting_number(geom))
      tmp(k:) = ' ' 
      tmp(k:) = ':recip vectors'
      s = s .and. rtdb_get(rtdb, tmp, mt_dbl, 3, 
     $     recip_lat_vectors(1,geom))
      tmp(k:) = ' ' 
      tmp(k:) = ':recip angles'
      s = s .and. rtdb_get(rtdb, tmp, mt_dbl, 3, 
     $     recip_lat_angles(1,geom))
      tmp(k:) = ' ' 
      tmp(k:) = ':direct volume'
      s = s .and. rtdb_get(rtdb, tmp, mt_dbl, 1, volume_direct(geom))
      tmp(k:) = ' ' 
      tmp(k:) = ':unique centers'
      s = s .and. rtdb_get(rtdb, tmp, mt_int, ncenter_unique(geom), 
     $     unique_cent(1,geom))
      tmp(k:) = ' ' 
      tmp(k:) = ':tags'
      s = s .and. rtdb_cget(rtdb, tmp, max_cent, tags(1,geom))
      tmp(k:) = ' ' 
      tmp(k:) = ':include_bqbq'
      s = s .and. rtdb_get(rtdb, tmp, mt_log, 1, include_bqbq(geom))
c     
c     Zmatrix info
c     
      tmp(k:) = ' ' 
      tmp(k:) = ':zmt_source'
      s = s .and. rtdb_cget(rtdb, tmp, 1, zmt_source(geom))
      if (zmt_source(geom) .ne. ' ') then
         tmp(k:) = ' ' 
         tmp(k:) = ':zmt_nizmat'
         s = s .and. rtdb_get(rtdb, tmp, mt_int, 1, zmt_nizmat(geom))
         tmp(k:) = ' ' 
         tmp(k:) = ':zmt_izmat'
         s = s .and. rtdb_get(rtdb, tmp, mt_int, zmt_nizmat(geom), 
     $        zmt_izmat(1,geom))
         tmp(k:) = ' ' 
         tmp(k:) = ':zmt_nzfrz'
         if (rtdb_get(rtdb, tmp, mt_int, 1, zmt_nzfrz(geom))) then
            tmp(k:) = ' ' 
            tmp(k:) = ':zmt_izfrz'
            s = s .and. rtdb_get(rtdb, tmp, mt_int, zmt_nzfrz(geom), 
     $           zmt_izfrz(1,geom))
            tmp(k:) = ' ' 
            tmp(k:) = ':zmt_izfrz_val'
            s = s .and. rtdb_get(rtdb, tmp, mt_dbl, zmt_nzfrz(geom), 
     $           zmt_izfrz_val(1,geom))
         endif
         tmp(k:) = ' ' 
         tmp(k:) = ':zmt_nzvar'
         s = s .and. rtdb_get(rtdb, tmp, mt_int, 1, zmt_nzvar(geom))
         tmp(k:) = ' ' 
         tmp(k:) = ':zmt_varsign'
         s = s .and. rtdb_get(rtdb, tmp, mt_dbl, zmt_nzvar(geom), 
     $        zmt_varsign(1,geom))
         tmp(k:) = ' ' 
         tmp(k:) = ':zmt_varname'
         s = s .and. rtdb_cget(rtdb, tmp, zmt_nzvar(geom), 
     $        zmt_varname(1,geom))
         tmp(k:) = ' ' 
         tmp(k:) = ':zmt_maxtor'
         s = s .and. rtdb_get(rtdb, tmp, mt_int, 1, zmt_maxtor(geom))
         tmp(k:) = ' ' 
         tmp(k:) = ':zmt_ijbond'
         s = s .and. rtdb_get(rtdb, tmp, mt_int, 2*max_zcoord, 
     $        zmt_ijbond(1,1,geom))
         tmp(k:) = ' ' 
         tmp(k:) = ':zmt_ijkang'
         s = s .and. rtdb_get(rtdb, tmp, mt_int, 3*max_zcoord, 
     $        zmt_ijkang(1,1,geom))
         tmp(k:) = ' ' 
         tmp(k:) = ':zmt_ijklto'
         s = s .and. rtdb_get(rtdb, tmp, mt_int, 4*max_zcoord, 
     $        zmt_ijklto(1,1,geom))
         tmp(k:) = ' ' 
         tmp(k:) = ':zmt_ijklop'
         s = s .and. rtdb_get(rtdb, tmp, mt_int, 4*max_zcoord, 
     $        zmt_ijklop(1,1,geom))
         tmp(k:) = ' ' 
         tmp(k:) = ':zmt_ijklnb'
         s = s .and. rtdb_get(rtdb, tmp, mt_int, 4*max_zcoord, 
     $        zmt_ijklnb(1,1,geom))
*     
         tmp(k:) = ' ' 
         tmp(k:) = ':zmt_ijbond_val'
         s = s .and. rtdb_get(rtdb, tmp, mt_dbl, max_zcoord, 
     $        zmt_ijbond_val(1,geom))
         tmp(k:) = ' ' 
         tmp(k:) = ':zmt_ijkang_val'
         s = s .and. rtdb_get(rtdb, tmp, mt_dbl, max_zcoord, 
     $        zmt_ijkang_val(1,geom))
         tmp(k:) = ' ' 
         tmp(k:) = ':zmt_ijklto_val'
         s = s .and. rtdb_get(rtdb, tmp, mt_dbl, max_zcoord, 
     $        zmt_ijklto_val(1,geom))
         tmp(k:) = ' ' 
         tmp(k:) = ':zmt_ijklop_val'
         s = s .and. rtdb_get(rtdb, tmp, mt_dbl, max_zcoord, 
     $        zmt_ijklop_val(1,geom))
         tmp(k:) = ' ' 
         tmp(k:) = ':zmt_ijklnb_val'
         s = s .and. rtdb_get(rtdb, tmp, mt_dbl, max_zcoord, 
     $        zmt_ijklnb_val(1,geom))
*     
         tmp(k:) = ' ' 
         tmp(k:) = ':zmt_ijbond_frz'
         s = s .and. rtdb_get(rtdb, tmp, mt_log, max_zcoord, 
     $        zmt_ijbond_frz(1,geom))
         tmp(k:) = ' ' 
         tmp(k:) = ':zmt_ijkang_frz'
         s = s .and. rtdb_get(rtdb, tmp, mt_log, max_zcoord, 
     $        zmt_ijkang_frz(1,geom))
         tmp(k:) = ' ' 
         tmp(k:) = ':zmt_ijklto_frz'
         s = s .and. rtdb_get(rtdb, tmp, mt_log, max_zcoord, 
     $        zmt_ijklto_frz(1,geom))
         tmp(k:) = ' ' 
         tmp(k:) = ':zmt_ijklop_frz'
         s = s .and. rtdb_get(rtdb, tmp, mt_log, max_zcoord, 
     $        zmt_ijklop_frz(1,geom))
         tmp(k:) = ' ' 
         tmp(k:) = ':zmt_ijklnb_frz'
         s = s .and. rtdb_get(rtdb, tmp, mt_log, max_zcoord, 
     $        zmt_ijklnb_val(1,geom))
*     
         tmp(k:) = ' ' 
         tmp(k:) = ':zmt_ijbond_nam'
         s = s .and. rtdb_cget(rtdb, tmp, max_zcoord, 
     $        zmt_ijbond_nam(1,geom))
         tmp(k:) = ' ' 
         tmp(k:) = ':zmt_ijkang_nam'
         s = s .and. rtdb_cget(rtdb, tmp, max_zcoord, 
     $        zmt_ijkang_nam(1,geom))
         tmp(k:) = ' ' 
         tmp(k:) = ':zmt_ijklto_nam'
         s = s .and. rtdb_cget(rtdb, tmp, max_zcoord, 
     $        zmt_ijklto_nam(1,geom))
         tmp(k:) = ' ' 
         tmp(k:) = ':zmt_ijklop_nam'
         s = s .and. rtdb_cget(rtdb, tmp, max_zcoord, 
     $        zmt_ijklop_nam(1,geom))
         tmp(k:) = ' ' 
         tmp(k:) = ':zmt_ijklnb_nam'
         s = s .and. rtdb_cget(rtdb, tmp, max_zcoord, 
     $        zmt_ijklnb_nam(1,geom))
*     
         tmp(k:) = ' ' 
         tmp(k:) = ':zmt_cvr_scaling'
         s = s .and. rtdb_get(rtdb, tmp, mt_dbl, 1,
     $        zmt_cvr_scaling(geom))
      endif
c     
c--   > get symmetry operators, number of operators and operator/atom
c     map from rtdb
c     

      tmp(k:) = ' ' 
      tmp(k:) = ' ' 
      tmp(k:) = ':num_operators'
      s = s .and. 
     $     rtdb_get(rtdb, tmp, mt_int, 1, sym_num_ops(geom))
      tmp(k:) = ' ' 
      tmp(k:) = ':operators'
      s = s .and. 
     $     rtdb_get(rtdb, tmp, mt_dbl, max_sym_ops*3*4,
     $     sym_ops(1,1,geom))
      if (sym_num_ops(geom) .gt. 0) then
c     
c     If loading into an old geometry free this memory
c     
         if (sym_center_map_handle(geom) .ne. -1) then
            if (.not. ma_free_heap(sym_center_map_handle(geom)))
     $           call errquit('geom_rtdb_load: free of atom map', 0,
     &       MA_ERR)
         end if
c     
         tmp(k:) = ' ' 
         tmp(k:) = ' ' 
         tmp(k:) = ':map_atoms'
         s = s .and. 
     $        rtdb_ma_get(rtdb, tmp, ma_type, nelem,
     $        sym_center_map_handle(geom))
         if (nelem .ne. sym_num_ops(geom)*ncenter(geom)) call errquit
     $        ('geom_rtdb_load: invalid no. of element in sym Tap',
     $        nelem, RTDB_ERR)
         if (.not. ma_get_index(sym_center_map_handle(geom),
     $        sym_center_map_index(geom)))call errquit
     $        ('geom_rtdb_load: bad ma handle for sym map', 0, MA_ERR)
      else
         sym_center_map_handle(geom) = -1
         sym_center_map_index(geom) = 1 ! Not used but address is created
      endif
c     
      if (.not. s) then
        if (ga_nodeid().eq.0) then
          write(LuOut,*)' geom_rtdb_load: rtdb corrupt: ',
     $                  names(geom)(1:lenn(geom)), ' -> ',
     $                  trans(geom)(1:lent(geom))
          call geom_err_info('geom_rtdb_load')
        endif
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
c     periodic systems: find conversion factor for geometrical parameters
c     
      if (isystype(geom) .gt. 0) then 
         if (.not. geom_get_user_scale(geom,scale))
     $        call errquit('geom_rtdb_load:failed get user scale',0,
     &       GEOM_ERR)
      endif
c     
c     setup geometry related stuff particular to the dimension of the system
c     
      if (isystype(geom) .eq. 3) then
         call geom_3d_amatrix(geom,scale)
      elseif(isystype(geom).eq.2) then
         call geom_2d_amatrix(geom,scale)
      elseif(isystype(geom).eq.1) then
         call geom_1d(geom,scale)
      endif
c
c     hack to fix numerical gradient issue when symmetry changes
c
      if(sym_num_ops(geom) .gt. 0) then
        if (rtdb_get(rtdb,'geom:getsym', mt_log, 1, getsym)) then
          if(getsym) then
            call geom_getsym(rtdb,geom,'geometry')
          endif
        endif
      endif
c     
*     do node = 0, ga_nnodes()-1
*     call ga_sync
*     if (ga_nodeid() .eq. node) then
*     write(LuOut,*) ' node ', ga_nodeid()
*     call sym_print_all(geom, .true., .true., .true., .true., .true.)
*     call util_flush(LuOut)
*     endif
*     call ga_sync
*     enddo
      
c     
      end
c
C> \brief Compute and store the nuclear repulsion and nuclear dipole -
C> external field energies
c
C> Compute and store the energies that depend on the atomic positions,
C> such as the nuclear - nuclear repulsion energy, and the nuclear
C> dipole - external field interaction energy.
c
      subroutine geom_compute_values(geom)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
#include "util.fh"
#include "inp.fh"
#include "stdio.fh"
      integer geom !< [Input] the geometry handle
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
      double precision ddot
      external ddot
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
*debug:      write(LuOut,*)' interaction of nuclear dipole ',
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
c
C> \brief Look up whether Bq - Bq interactions should be calculated
c
C> Point charges (Bq centers) have many uses in quantum chemistry
C> models. In some applications Bq - Bq interactions are an important
C> component of the energy expression, in other applications these
C> interactions should be omitted. This function returns what has been
C> specified for this interaction in this geometry instance.
c
C> \return Return .true. if Bq - Bq interactions should be evaluated,
C> and .false. otherwise.
c
      logical function geom_include_bqbq(geom)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "geomP.fh"
      integer geom !< [Input] the geometry handle
      logical geom_check_handle
      external geom_check_handle
c
      if (.not. geom_check_handle(geom, 'geom_include_bqbq'))
     $     call errquit('geom_include_bqbq: bad handle',0, GEOM_ERR)
      geom_include_bqbq = include_bqbq(geom)
c
      end
      logical function geom_set_bqbq(geom, value)
      implicit none
#include "nwc_const.fh"
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
c
C> \brief Store a geometry on the RTDB
c
C> Store a geometry onto the RTDB with a specified key, if no key is specified
C> a key will be constructed from the current geometry name.
C> If storing the geometry was successful return .true., return .false. otherwise.
      logical function geom_rtdb_store(rtdb, geom, name)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "stdio.fh"
***** #include "context.fh"
#include "inp.fh"
c     
      integer rtdb              !< [Input] the RTDB handle
      character*(*) name        !< [Input] the geometry RTDB key
      integer geom              !< [Input] the geometry handle
      logical geom_check_handle, geom_rtdb_add, geom_rtdb_delete
      external geom_check_handle, geom_rtdb_add, geom_rtdb_delete
      logical s
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
      s = geom_rtdb_delete(rtdb, name) ! Delete any old junk
c     
c     try to translate the name
c     
      trans(geom) = 'junk'
      if (.not. rtdb_cget(rtdb, name, 1, trans(geom)))
     $     trans(geom) = name
*     if (.not. context_rtdb_match(rtdb, name, trans(geom)))
*     $     trans(geom) = name
      lent(geom) = inp_strlen(trans(geom))
c     
c     now put the info into the data base
c     
      tmp = 'geometry:'//trans(geom)(1:lent(geom))
      k = inp_strlen(tmp)+1
      s = .true.
c     
      tmp(k:) = ' ' 
      tmp(k:) = ':ncenter'
      s = s .and. rtdb_put(rtdb, tmp, mt_int, 1, ncenter(geom))
      tmp(k:) = ' ' 
      tmp(k:) = ':coords'
      s = s .and. rtdb_put(rtdb, tmp, mt_dbl, ncenter(geom)*3,
     $     coords(1,1,geom))
      tmp(k:) = ' ' 
      tmp(k:) = ':vel'
      s = s .and. rtdb_put(rtdb, tmp, mt_dbl, ncenter(geom)*3,
     $     velocities(1,1,geom))
      tmp(k:) = ' ' 
      tmp(k:) = ':charges'
      s = s .and. rtdb_put(rtdb, tmp, mt_dbl, ncenter(geom),
     $     charge(1,geom))
      tmp(k:) = ' ' 
      tmp(k:) = ':masses'
      s = s .and. rtdb_put(rtdb, tmp, mt_dbl, ncenter(geom),
     $     geom_mass(1,geom))
      tmp(k:) = ' ' 
      tmp(k:) = ':atomct'
      s = s .and. rtdb_put(rtdb, tmp, mt_dbl, ncenter(geom),
     $     geom_atomct(1,geom))
C     new
      tmp(k:) = ' ' 
      tmp(k:) = ':inv nuc expon'
      s = s .and. rtdb_put(rtdb, tmp, mt_dbl, ncenter(geom),
     $     geom_invnucexp(1,geom))
C     end
      tmp(k:) = ' ' 
      tmp(k:) = ':efield'
      s = s .and. rtdb_put(rtdb, tmp, mt_dbl, 3, efield(1,geom))

c     tmp(k:) = ' ' 
c     tmp(k:) = ':lattice vectors'
c     s = s .and. rtdb_put(rtdb, tmp, mt_dbl, 3, 
c    $     lattice_vectors(1,geom))
c     tmp(k:) = ' ' 
c     tmp(k:) = ':lattice angles'
c     s = s .and. rtdb_put(rtdb, tmp, mt_dbl, 3, lattice_angles(1,geom))
      tmp(k:) = ' ' 
      tmp(k:) = ':amatrix'
      s = s .and. rtdb_put(rtdb, tmp, mt_dbl, 9, 
     $     amatrix(1,1,geom))

      tmp(k:) = ' ' 
      tmp(k:) = ':system type'
      s = s .and. rtdb_put(rtdb, tmp, mt_int, 1, isystype(geom))
      tmp(k:) = ' ' 
      tmp(k:) = ':no. unique centers'
      s = s .and. rtdb_put(rtdb, tmp, mt_int, 1, ncenter_unique(geom))
      tmp(k:) = ' ' 
      tmp(k:) = ':group number'
      s = s .and. rtdb_put(rtdb, tmp, mt_int, 1, group_number(geom))
      tmp(k:) = ' ' 
      tmp(k:) = ':group name'
      s = s .and. rtdb_cput(rtdb, tmp, 1, group_name(geom))
      tmp(k:) = ' '
      tmp(k:) = ':use_primitive'
      s = s .and. rtdb_put(rtdb, tmp, mt_log, 1, use_primitive(geom))
      tmp(k:) = ' '
      tmp(k:) = ':primitive_center'
      s = s .and. rtdb_cput(rtdb, tmp, 1, primitive_center(geom))
      tmp(k:) = ' ' 
      tmp(k:) = ':user units'
      s = s .and. rtdb_cput(rtdb, tmp, 1, user_units(geom))
      tmp(k:) = ' ' 
      tmp(k:) = ':angstrom_to_au'
      s = s .and. rtdb_put(rtdb, tmp, mt_dbl, 1, angstrom_to_au)
      tmp(k:) = ' ' 
      tmp(k:) = ':setting number'
      s = s .and. rtdb_put(rtdb, tmp, mt_int, 1, setting_number(geom))
      tmp(k:) = ' ' 
      tmp(k:) = ':recip vectors'
      s = s .and. rtdb_put(rtdb, tmp, mt_dbl, 3, 
     $     recip_lat_vectors(1,geom))
      tmp(k:) = ' ' 
      tmp(k:) = ':recip angles'
      s = s .and. rtdb_put(rtdb, tmp, mt_dbl, 3, 
     $     recip_lat_angles(1,geom))
      tmp(k:) = ' ' 
      tmp(k:) = ':direct volume'
      s = s .and. rtdb_put(rtdb, tmp, mt_dbl, 1, volume_direct(geom))
      tmp(k:) = ' ' 
      tmp(k:) = ':unique centers'
      s = s .and. rtdb_put(rtdb, tmp, mt_int, ncenter_unique(geom), 
     $     unique_cent(1,geom))
      tmp(k:) = ' ' 
      tmp(k:) = ':tags'
      s = s .and. rtdb_cput(rtdb, tmp, ncenter(geom), tags(1,geom))
      tmp(k:) = ' ' 
      tmp(k:) = ':include_bqbq'
      s = s .and. rtdb_put(rtdb, tmp, mt_log, 1, include_bqbq(geom))
c     
c     Zmatrix info
c     
      tmp(k:) = ' ' 
      tmp(k:) = ':zmt_source'
      s = s .and. rtdb_cput(rtdb, tmp, 1, zmt_source(geom))
      if (zmt_source(geom) .ne. ' ') then
         tmp(k:) = ' ' 
         tmp(k:) = ':zmt_nizmat'
         s = s .and. rtdb_put(rtdb, tmp, mt_int, 1, zmt_nizmat(geom))
         tmp(k:) = ' ' 
         tmp(k:) = ':zmt_izmat'
         s = s .and. rtdb_put(rtdb, tmp, mt_int, zmt_nizmat(geom), 
     $        zmt_izmat(1,geom))
         if (zmt_nzfrz(geom) .gt. 0) then
            tmp(k:) = ' ' 
            tmp(k:) = ':zmt_nzfrz'
            s = s .and. rtdb_put(rtdb, tmp, mt_int, 1, zmt_nzfrz(geom))
            tmp(k:) = ' ' 
            tmp(k:) = ':zmt_izfrz'
            s = s .and. rtdb_put(rtdb, tmp, mt_int, zmt_nzfrz(geom), 
     $           zmt_izfrz(1,geom))
            tmp(k:) = ' ' 
            tmp(k:) = ':zmt_izfrz_val'
            s = s .and. rtdb_put(rtdb, tmp, mt_dbl, zmt_nzfrz(geom), 
     $           zmt_izfrz_val(1,geom))
         endif
         tmp(k:) = ' ' 
         tmp(k:) = ':zmt_nzvar'
         s = s .and. rtdb_put(rtdb, tmp, mt_int, 1, zmt_nzvar(geom))
         tmp(k:) = ' ' 
         tmp(k:) = ':zmt_varsign'
         s = s .and. rtdb_put(rtdb, tmp, mt_dbl, zmt_nzvar(geom), 
     $        zmt_varsign(1,geom))
         tmp(k:) = ' ' 
         tmp(k:) = ':zmt_varname'
         s = s .and. rtdb_cput(rtdb, tmp, zmt_nzvar(geom), 
     $        zmt_varname(1,geom))
         tmp(k:) = ' ' 
         tmp(k:) = ':zmt_maxtor'
         s = s .and. rtdb_put(rtdb, tmp, mt_int, 1, zmt_maxtor(geom))
         tmp(k:) = ' ' 
         tmp(k:) = ':zmt_ijbond'
         s = s .and. rtdb_put(rtdb, tmp, mt_int, 2*max_zcoord, 
     $        zmt_ijbond(1,1,geom))
         tmp(k:) = ' ' 
         tmp(k:) = ':zmt_ijkang'
         s = s .and. rtdb_put(rtdb, tmp, mt_int, 3*max_zcoord, 
     $        zmt_ijkang(1,1,geom))
         tmp(k:) = ' ' 
         tmp(k:) = ':zmt_ijklto'
         s = s .and. rtdb_put(rtdb, tmp, mt_int, 4*max_zcoord, 
     $        zmt_ijklto(1,1,geom))
         tmp(k:) = ' ' 
         tmp(k:) = ':zmt_ijklop'
         s = s .and. rtdb_put(rtdb, tmp, mt_int, 4*max_zcoord, 
     $        zmt_ijklop(1,1,geom))
         tmp(k:) = ' ' 
         tmp(k:) = ':zmt_ijklnb'
         s = s .and. rtdb_put(rtdb, tmp, mt_int, 4*max_zcoord, 
     $        zmt_ijklnb(1,1,geom))
*     
         tmp(k:) = ' ' 
         tmp(k:) = ':zmt_ijbond_val'
         s = s .and. rtdb_put(rtdb, tmp, mt_dbl, max_zcoord, 
     $        zmt_ijbond_val(1,geom))
         tmp(k:) = ' ' 
         tmp(k:) = ':zmt_ijkang_val'
         s = s .and. rtdb_put(rtdb, tmp, mt_dbl, max_zcoord, 
     $        zmt_ijkang_val(1,geom))
         tmp(k:) = ' ' 
         tmp(k:) = ':zmt_ijklto_val'
         s = s .and. rtdb_put(rtdb, tmp, mt_dbl, max_zcoord, 
     $        zmt_ijklto_val(1,geom))
         tmp(k:) = ' ' 
         tmp(k:) = ':zmt_ijklop_val'
         s = s .and. rtdb_put(rtdb, tmp, mt_dbl, max_zcoord, 
     $        zmt_ijklop_val(1,geom))
         tmp(k:) = ' ' 
         tmp(k:) = ':zmt_ijklnb_val'
         s = s .and. rtdb_put(rtdb, tmp, mt_dbl, max_zcoord, 
     $        zmt_ijklnb_val(1,geom))
*     
         tmp(k:) = ' ' 
         tmp(k:) = ':zmt_ijbond_frz'
         s = s .and. rtdb_put(rtdb, tmp, mt_log, max_zcoord, 
     $        zmt_ijbond_frz(1,geom))
         tmp(k:) = ' ' 
         tmp(k:) = ':zmt_ijkang_frz'
         s = s .and. rtdb_put(rtdb, tmp, mt_log, max_zcoord, 
     $        zmt_ijkang_frz(1,geom))
         tmp(k:) = ' ' 
         tmp(k:) = ':zmt_ijklto_frz'
         s = s .and. rtdb_put(rtdb, tmp, mt_log, max_zcoord, 
     $        zmt_ijklto_frz(1,geom))
         tmp(k:) = ' ' 
         tmp(k:) = ':zmt_ijklop_frz'
         s = s .and. rtdb_put(rtdb, tmp, mt_log, max_zcoord, 
     $        zmt_ijklop_frz(1,geom))
         tmp(k:) = ' ' 
         tmp(k:) = ':zmt_ijklnb_frz'
         s = s .and. rtdb_put(rtdb, tmp, mt_log, max_zcoord, 
     $        zmt_ijklnb_val(1,geom))
*     
         tmp(k:) = ' ' 
         tmp(k:) = ':zmt_ijbond_nam'
         s = s .and. rtdb_cput(rtdb, tmp, max_zcoord, 
     $        zmt_ijbond_nam(1,geom))
         tmp(k:) = ' ' 
         tmp(k:) = ':zmt_ijkang_nam'
         s = s .and. rtdb_cput(rtdb, tmp, max_zcoord, 
     $        zmt_ijkang_nam(1,geom))
         tmp(k:) = ' ' 
         tmp(k:) = ':zmt_ijklto_nam'
         s = s .and. rtdb_cput(rtdb, tmp, max_zcoord, 
     $        zmt_ijklto_nam(1,geom))
         tmp(k:) = ' ' 
         tmp(k:) = ':zmt_ijklop_nam'
         s = s .and. rtdb_cput(rtdb, tmp, max_zcoord, 
     $        zmt_ijklop_nam(1,geom))
         tmp(k:) = ' ' 
         tmp(k:) = ':zmt_ijklnb_nam'
         s = s .and. rtdb_cput(rtdb, tmp, max_zcoord, 
     $        zmt_ijklnb_nam(1,geom))
*     
         tmp(k:) = ' ' 
         tmp(k:) = ':zmt_cvr_scaling'
         s = s .and. rtdb_put(rtdb, tmp, mt_dbl, 1, 
     $        zmt_cvr_scaling(geom))
      endif
c     
c--   > put symmetry operators, number of operators and operator/atom
c     map out to rtdb
c     
      tmp(k:) = ' ' 
      tmp(k:) = ':num_operators'
      s = s .and. 
     $     rtdb_put(rtdb, tmp, mt_int, 1, sym_num_ops(geom))
      tmp(k:) = ' ' 
      tmp(k:) = ':operators'
      s = s .and. 
     $     rtdb_put(rtdb, tmp, mt_dbl, max_sym_ops*3*4,
     $     sym_ops(1,1,geom))
      if (sym_num_ops(geom) .gt. 0) then
         tmp(k:) = ' ' 
         tmp(k:) = ':map_atoms'
         s = s .and. 
     $        rtdb_put(rtdb, tmp, mt_int,
     $        ncenter(geom)*sym_num_ops(geom),
     $        int_mb(sym_center_map_index(geom)))
      endif

c     
c     
c     insert translated name into list of known geometries
c     
      s = s .and. geom_rtdb_add(rtdb, trans(geom))
c     
c     check that all rtdb operations were successful
c     
      if (.not. s) then
         write(LuOut,*) ' geom_rtdb_store: write to rtdb failed',
     $        names(geom)(1:lenn(geom)), ' -> ',
     $        trans(geom)(1:lent(geom))
         call geom_err_info('geom_rtdb_store')
         geom_rtdb_store = .false.
         return
      end if
      geom_rtdb_store = .true.
c     
      end
c
C> \brief Delete a geometry from the RTDB
c
C> Delete a geometry with a given name from the RTDB.
c
C> \return Return .true. if the geometry was successfully deleted,
C> return .false. otherwise.
      logical function geom_rtdb_delete(rtdb, name)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "geomP.fh"
#include "rtdb.fh"
#include "inp.fh"
#include "global.fh"
c
      integer rtdb              !< [Input] the RTDB handle
      character*(*) name        !< [Input] the geometry name
      character*256 translation, tmp, test
      integer lt, geom, geom2, k
      logical status, mode
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
      goto 11
 10   do geom2 = geom+1, ngeom_rtdb ! Matched
         names_rtdb(geom2-1) = names_rtdb(geom2)
      end do
      ngeom_rtdb = ngeom_rtdb - 1
      status = geom_rtdb_out(rtdb)
c
c     Delete junk in rtdb even if did not find geometry in
c     the list just in case things are a little messed up
c
 11   if (ga_nodeid() .eq. 0) then
         mode = rtdb_parallel(.false.)
c     
c     delete each entry assoicated with a geometry in the database
c     
         tmp = 'geometry:'//translation(1:lt)
         k = inp_strlen(tmp)
         k = k + 1
         tmp(k:k) = ':'
c     
         status = rtdb_first(rtdb, test)
 20      if (status) then
            if (inp_compare(.true.,tmp(1:k),test(1:k))) then
               if (.not. rtdb_delete(rtdb,test)) call errquit
     $              ('geom_rtdb_delete:failed deleting known entry',0,
     &       RTDB_ERR)
            endif
            status = rtdb_next(rtdb, test)
            goto 20
         endif
         mode = rtdb_parallel(mode) ! Restore previous state
      endif
c
      geom_rtdb_delete = .true.
c
      end
c
C> \brief Reset the symmetry to C1 for a given geometry instance
c
C> Resets the symmetry information to C1 for the specified geometry
C> instance. This also frees any associated in-core data structures.
c
C> \return Return .true. if successfull, and .false. otherwise.
c
      logical function geom_strip_sym(geom)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "geomP.fh"
#include "mafdecls.fh"
c
c     Reset the given geometry to have just C1 symmetry, freeing
c     any associated in-core data structures.
c
      integer geom              !< [Input] the geometry handle
      integer i
      logical geom_check_handle
      external geom_check_handle
c     
      geom_strip_sym = geom_check_handle(geom, 'geom_strip_sym')
      if (.not. geom_strip_sym) return
c
c      isystype(geom) = 0
      group_number(geom) = 1
      setting_number(geom) = 0
      if (sym_center_map_handle(geom) .ne. -1) then
         if (.not. ma_free_heap(sym_center_map_handle(geom)))
     $        call errquit('geom_strip_sum: free of atom map', 0,
     &       MA_ERR)
      end if
      sym_center_map_handle(geom) = -1
      sym_center_map_index(geom) = 1
      group_name(geom) = 'C1'
      sym_num_ops(geom) = 0
      use_primitive(geom) = .false.
      primitive_center(geom) = 'x'
c
      ncenter_unique(geom) = ncenter(geom)
      do i = 1, ncenter_unique(geom)
         unique_cent(i,geom) = i
      end do
c
      end
c
C> \brief Destroy a geometry instance
c
C> Destroys a geometry instance. After this operation the geometry handle is
C> no longer valid.
C> Returns .true. if the instance was successfully destroyed,
C> returns .false. otherwise
      logical function geom_destroy(geom)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "geomP.fh"
#include "mafdecls.fh"
c
      integer geom              !< [Input] the geometry handle
      integer i
      logical geom_check_handle, geom_check_cent
      external geom_check_handle, geom_check_cent
c     
      geom_destroy = geom_check_handle(geom, 'geom_destroy')
      if (.not. geom_destroy) return
c
      active(geom) = .false.
* this is set for a geometry at every basis set load
* This info needs to be nullified when the geometry is gone
      do i = 1,ncenter(geom)
        oecpcent(i,geom) = .false.
      enddo
      geom_destroy = .true.
      if (sym_center_map_handle(geom) .ne. -1) then
         if (.not. ma_free_heap(sym_center_map_handle(geom)))
     $        call errquit('geom_destroy: free of atom map', 0, MA_ERR)
      end if
c
      end
c
C> \brief Set the point group symmetry for a geometry instance
c
C> Sets the point group or space group for a given geometry instance.
C> Whether the point group or the space group is used depends on the
C> type of geometry under consideration.
c
C> \return Return .true. if successfull, and .false. otherwise.
c
      logical function geom_group_set(geom, group)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
#include "inp.fh"
c     
      integer geom              !< [Input] the geometry handle
      character*(*) group       !< [Input] the point/space group
      logical geom_check_handle, geom_check_cent
      external geom_check_handle, geom_check_cent
c
      geom_group_set = geom_check_handle(geom, 'geom_group_set')
      if (.not. geom_group_set) return
c
      if (isystype(geom).eq.0) then
         geom_group_set = inp_match(46,.false.,group,sym_molgnames,
     $        group_number(geom))
      else
         geom_group_set = inp_match(240,.false.,group,sym_spgnames,
     $        group_number(geom))

c        try car file style names
         if (.not. geom_group_set ) then
            geom_group_set = inp_match(230,.false.,group,sym_carnames,
     $        group_number(geom))
         endif

      endif
c
      end
c
C> \brief Define velocities for the centers
c
C> In dynamics simulations the centers in the system move with a
C> velocity. This function allows these velocities to be stored in 
C> a geometry instance. It is assumed that the number of centers in 
C> the geometry has already been defined.
c
C> \return Returns .true. if the velocities were stored successfully,
C> and .false. otherwise.
c
      logical function geom_vel_set(geom, vel)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
c     
      integer geom               !< [Input] the geometry handle
      double precision vel(3, *) !< [Input] the velocities
      logical geom_check_handle, geom_check_cent
      external geom_check_handle, geom_check_cent
c
      geom_vel_set = geom_check_handle(geom, 'geom_vel_set')
      if (.not. geom_vel_set) return
c
      call dcopy(3*ncenter(geom), vel, 1, velocities(1,1,geom), 1)
c
      end
c
C> \brief Retrieve the velocities of the centers
c
C> Retrieves the velocities of the centers in the specified geometry
C> instance.
c
C> \return Return .true. if the velocities we found successfully,
C> and .false. otherwise.
c
      logical function geom_vel_get(geom, vel)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
c     
      integer geom               !< [Input] the geometry handle
      double precision vel(3, *) !< [Output] the center velocities
      logical geom_check_handle, geom_check_cent
      external geom_check_handle, geom_check_cent
c
      geom_vel_get = geom_check_handle(geom, 'geom_vel_get')
      if (.not. geom_vel_get) return
c
      call dcopy(3*ncenter(geom), velocities(1,1,geom), 1, vel, 1)
c
      end
      function geom_cart_set_gen(geom, i0,ncent,nt,ns, t, c, q)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
#include "stdio.fh"
c     
      logical geom_cart_set_gen
      integer geom              ! [input]
      integer i0                ! [input]
      integer ncent             ! [input]
      integer nt             ! [input]
      integer ns             ! [input]
      character*1 t(nt*ns)     ! [input]
      double precision c(nt,3) ! [input]
      double precision q(nt) ! [input]
      logical geom_check_handle, geom_check_cent
      external geom_check_handle, geom_check_cent
      integer i,j
      double precision scale
      integer k
      character*16 atag
c
      geom_cart_set_gen = geom_check_handle(geom, 'geom_cart_set_gen')
      if (.not. geom_cart_set_gen) return
c
      if (ncent.le.0) then
         write(LuOut,*) ' geom_cart_set_gen: too few centers ',ncent,
     $        names(geom)(1:lenn(geom))
         geom_cart_set_gen = .false.
         return
      else if (ncent.gt.max_cent) then
         write(LuOut,*) ' geom_cart_set_gen: too many centers ',ncent,
     $        names(geom)(1:lenn(geom))
         geom_cart_set_gen = .false.
         return
      end if
c
      if (ncenter(geom).ne.ncent) then
         ncenter_unique(geom) = ncent ! Assume symmetry is unchanged !!!!!
      endif
      ncenter(geom) = ncent
      scale = angstrom_to_au
      do i = 1, ncent
         j = i0+i-1
         atag = ""
         do k=1,16
           atag(k:k) = t((j-1)*ns+k)
         end do
         tags(i,geom) = atag
         charge(i,geom) = q(j)
         coords(1,i,geom) = scale*c(j,1)
         coords(2,i,geom) = scale*c(j,2)
         coords(3,i,geom) = scale*c(j,3)
         unique_cent(i,geom) = i
      end do
c
      end
      function geom_cart_set_gen1(geom, i0,ncent,nt,ns, t, c, q)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
#include "stdio.fh"
c     
      logical geom_cart_set_gen1
      integer geom              ! [input]
      integer i0                ! [input]
      integer ncent             ! [input]
      integer nt             ! [input]
      integer ns             ! [input]
      character*1 t(nt*ns)     ! [input]
      double precision c(3,nt) ! [input]
      double precision q(nt) ! [input]
      logical geom_check_handle, geom_check_cent
      external geom_check_handle, geom_check_cent
      integer i,j
      double precision scale
      integer k
      character*16 atag
c
      geom_cart_set_gen1 = geom_check_handle(geom, 'geom_cart_set')
      if (.not. geom_cart_set_gen1) return
c
      if (ncent.le.0) then
         write(LuOut,*) ' geom_cart_set: too few centers ',ncent,
     $        names(geom)(1:lenn(geom))
         geom_cart_set_gen1 = .false.
         return
      else if (ncent.gt.max_cent) then
         write(LuOut,*) ' geom_cart_set: too many centers ',ncent,
     $        names(geom)(1:lenn(geom))
         geom_cart_set_gen1 = .false.
         return
      end if
c
      if (ncenter(geom).ne.ncent) then
         ncenter_unique(geom) = ncent ! Assume symmetry is unchanged !!!!!
      endif
      ncenter(geom) = ncent
      scale = angstrom_to_au
      do i = 1, ncent
         j = i0+i-1
         atag = ""
         do k=1,16
           atag(k:k) = t((j-1)*ns+k)
         end do
         tags(i,geom) = atag
         charge(i,geom) = q(j)
         coords(1,i,geom) = scale*c(1,j)
         coords(2,i,geom) = scale*c(2,j)
         coords(3,i,geom) = scale*c(3,j)
         unique_cent(i,geom) = i
      end do
c
      end
      function geom_cart_set1(geom, i0,ncent,nt, t, c, q)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
#include "stdio.fh"
c     
      logical geom_cart_set1
      integer geom              ! [input]
      integer i0                ! [input]
      integer ncent             ! [input]
      integer nt             ! [input]
      character*16 t(nt)     ! [input]
      double precision c(nt,3) ! [input]
      double precision q(nt) ! [input]
      logical geom_check_handle, geom_check_cent
      external geom_check_handle, geom_check_cent
      integer i,j
      double precision scale
c
      geom_cart_set1 = geom_check_handle(geom, 'geom_cart_set')
      if (.not. geom_cart_set1) return
c
      if (ncent.le.0) then
         write(LuOut,*) ' geom_cart_set: too few centers ',ncent,
     $        names(geom)(1:lenn(geom))
         geom_cart_set1 = .false.
         return
      else if (ncent.gt.max_cent) then
         write(LuOut,*) ' geom_cart_set: too many centers ',ncent,
     $        names(geom)(1:lenn(geom))
         geom_cart_set1 = .false.
         return
      end if
c
      if (ncenter(geom).ne.ncent) then
         ncenter_unique(geom) = ncent ! Assume symmetry is unchanged !!!!!
      endif
      ncenter(geom) = ncent
      scale = angstrom_to_au
      do i = 1, ncent
         j = i0+i-1
         tags(i,geom) = t(j)
         charge(i,geom) = q(j)
         coords(1,i,geom) = scale*c(j,1)
         coords(2,i,geom) = scale*c(j,2)
         coords(3,i,geom) = scale*c(j,3)
         unique_cent(i,geom) = i
      end do
c
      end
      logical function geom_cart_set(geom, ncent, t, c, q)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
#include "stdio.fh"
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
      if (ncent.le.0) then
         write(LuOut,*) ' geom_cart_set: too few centers ',ncent,
     $        names(geom)(1:lenn(geom))
         geom_cart_set = .false.
         return
      else if (ncent.gt.max_cent) then
         write(LuOut,*) ' geom_cart_set: too many centers ',ncent,
     $        names(geom)(1:lenn(geom))
         geom_cart_set = .false.
         return
      end if
c
      if (ncenter(geom).ne.ncent) then
         ncenter_unique(geom) = ncent ! Assume symmetry is unchanged !!!!!
      endif
      ncenter(geom) = ncent
      do i = 1, ncent
         tags(i,geom) = t(i)
         charge(i,geom) = q(i)
         coords(1,i,geom) = c(1,i)
         coords(2,i,geom) = c(2,i)
         coords(3,i,geom) = c(3,i)
         unique_cent(i,geom) = i
      end do
c
      end
c
C> \brief Extract only the coordinates and charges from a geometry
C> instance
c
C> \return Return .true. if successfull, and .false. otherwise.
c
      logical function geom_efc_cart_get(geom, ncent, c, q)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
c
      integer geom                 ! [Input] the geometry handle
      integer ncent                ! [Output] the number of centers
      double precision c(3, ncent) ! [Output] the coordinates
      double precision q(ncent)    ! [Output] the charges
      logical geom_check_handle, geom_check_cent
      external geom_check_handle, geom_check_cent
      integer i
c
c-    geom_efc_cart_get = geom_check_handle(geom, 'geom_efc_cart_get')
c-    if (.not. geom_efc_cart_get) return
c
      ncent = ncenter(geom)
      do i = 1, ncent
         q(i) = charge(i,geom)
         c(1,i) = coords(1,i,geom)
         c(2,i) = coords(2,i,geom)
         c(3,i) = coords(3,i,geom)
      end do
      geom_efc_cart_get = .true.
c
      end
c
C> \brief Define only the coordinates and charges of a geometry
C> instance
c
C> \return Return .true. if successfull, and .false. otherwise.
c
      logical function geom_efc_cart_set(geom, ncent, c, q)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
c
      integer geom                 ! [Input] the geometry handle
      integer ncent                ! [Input] the number of centers
      double precision c(3, ncent) ! [Input] the coordinates
      double precision q(ncent)    ! [Input] the charges
      logical geom_check_handle, geom_check_cent
      external geom_check_handle, geom_check_cent
      integer i
c
c-    geom_efc_cart_set = geom_check_handle(geom, 'geom_efc_cart_set')
c-    if (.not. geom_efc_cart_set) return
c
      ncenter(geom) = ncent
      do i = 1, ncent
         charge(i,geom) = q(i)
         coords(1,i,geom) = c(1,i)
         coords(2,i,geom) = c(2,i)
         coords(3,i,geom) = c(3,i)
      end do
      geom_efc_cart_set = .true.
c
      end
c
C> \brief Extract the charges from a geometry instance
c
C> \return Return .true. if successfull, and .false. otherwise.
c
      function geom_cart_get_charges(geom, ncent,q)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
c     
      logical geom_cart_get_charges
      integer geom              !< [Input] the geometry handle
      integer ncent             !< [Output] the number of centers
      double precision q(ncent) !< [Output] the charges
      logical geom_check_handle, geom_check_cent
      external geom_check_handle, geom_check_cent
      integer i
c
      geom_cart_get_charges = geom_check_handle(geom, 'geom_cart_get')
      if (.not. geom_cart_get_charges) return
c
      ncent = ncenter(geom)
      do i = 1, ncent
         q(i) = charge(i,geom)
      end do
c
      end
c
C> \brief Extract the tags, coordinates and charges from a geometry
C> instance
c
C> \return Return .true. if successfull, and .false. otherwise.
c
      function geom_cart_get(geom, ncent, t, c, q)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
c     
      logical geom_cart_get
      integer geom                 !< [Input] the geometry handle
      integer ncent                !< [Output] the number of centers
      character*16 t(ncent)        !< [Output] the tags
      double precision c(3, ncent) !< [Output] the coordinates
      double precision q(ncent)    !< [Output] the charges
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

      function geom_cart_get1(geom, ncent, t, c)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
c     
      logical geom_cart_get1
      integer geom                 !< [Input] the geometry handle
      integer ncent                !< [Output] the number of centers
      character*16 t(ncent)        !< [Output] the tags
      double precision c(3, ncent) !< [Output] the coordinates
      logical geom_check_handle, geom_check_cent
      external geom_check_handle, geom_check_cent
      integer i
c
      geom_cart_get1 = geom_check_handle(geom, 'geom_cart_get')
      if (.not. geom_cart_get1) return
c
      ncent = ncenter(geom)
      do i = 1, ncent
         t(i) = tags(i,geom)
         c(1,i) = coords(1,i,geom)
         c(2,i) = coords(2,i,geom)
         c(3,i) = coords(3,i,geom)
      end do
c
      end
c
c
C> \brief Extract the tags, coordinates, charges and atomic numbers from
C> a geometry instance
c
C> Extracts the tags, coordinates, and charges from the geometry by
C> simply copying the data. The atomic numbers are based on a
C> translation of the tags. If the tag does not correspond to a chemical
C> element, e.g. 'X' or 'Bq', the atomic number is 0.
c
C> \return Return .true. if successfull, and .false. otherwise.
c
      function geom_cart_get2(geom, ncent, t, c, q, atnum)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
c
      logical geom_cart_get2
      logical status_tagi
      integer geom                 !< [Input] the geometry handle
      integer ncent                !< [Output] the number of centers
      character*16 t(ncent)        !< [Output] the tags
      double precision c(3, ncent) !< [Output] the coordinates
      double precision q(ncent)    !< [Output] the charges
      integer atnum(ncent)         !< [Output] the atomic numbers
      logical geom_check_handle, geom_check_cent
      external geom_check_handle, geom_check_cent
      integer i
      integer iatn
      character*2 symi
      character*16 elei
      logical geom_tag_to_element
      external geom_tag_to_element
c
      geom_cart_get2 = geom_check_handle(geom, 'geom_cart_get')
      if (.not. geom_cart_get2) return
c
      ncent = ncenter(geom)
      do i = 1, ncent
         t(i) = tags(i,geom)
         q(i) = charge(i,geom)
         c(1,i) = coords(1,i,geom)
         c(2,i) = coords(2,i,geom)
         c(3,i) = coords(3,i,geom)
         status_tagi = geom_tag_to_element(t(i),symi,elei,iatn)
         atnum(i) = iatn  ! iatn is 0 if status_tagi is false
      end do
c
      end
c
C> \brief Extracts the coordinates of all centers in a geometry
c
C> Extracts the coordinates of all centers in a geometry assuming
C> that the caller has made sure the buffer is big enough. If the
C> buffer is too small the results are undefined.
c
C> \return Returns .true. if the function was successfull,
C> and .false. otherwise.
      logical function geom_cart_coords_get(geom, c)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
c     
      integer geom             !< [Input] the geometry handle
      double precision c(3, *) !< [Output] the Cartesian coordinates
      logical geom_check_handle, geom_check_cent
      external geom_check_handle, geom_check_cent
      integer i, ncent
c
      geom_cart_coords_get = 
     $     geom_check_handle(geom, 'geom_cart_coords_get')
      if (.not. geom_cart_coords_get) return
c
      ncent = ncenter(geom)
      do i = 1, ncent
         c(1,i) = coords(1,i,geom)
         c(2,i) = coords(2,i,geom)
         c(3,i) = coords(3,i,geom)
      end do
c
      end
c
C> \brief Defines the coordinates of all centers in a geometry
c
C> Defines the coordinates of all centers in a geometry assuming
C> that the caller has previously defined how many centers there
C> are (see e.g. geom_ncent_set). 
c
C> \return Returns .true. if the function was successfull,
C> and .false. otherwise.
      logical function geom_cart_coords_set(geom, c)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
c     
      integer geom             !< [Input] the geometry handle
      double precision c(3, *) !< [Input] the coordinates
      logical geom_check_handle, geom_check_cent
      external geom_check_handle, geom_check_cent
      integer i, ncent
c
      geom_cart_coords_set = 
     $     geom_check_handle(geom, 'geom_cart_coords_set')
      if (.not. geom_cart_coords_set) return
c
      ncent = ncenter(geom)
      do i = 1, ncent
         coords(1,i,geom) = c(1,i)
         coords(2,i,geom) = c(2,i)
         coords(3,i,geom) = c(3,i)
      end do
c
      end
c
C> \brief Extracts the coordinates and velocities of all centers in a geometry
c
C> Extracts the coordinates and velocities of all centers in a geometry assuming
C> that the caller has made sure the buffer is big enough. If the
C> buffer is too small the results are undefined.
c
C> \return Returns .true. if the function was successfull,
C> and .false. otherwise.
      logical function geom_coords_vels_get(geom, c, v)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
c
      integer geom            !< [Input] the geometry handle
      double precision c(3,*) !< [Output] cartesian coordinates
      double precision v(3,*) !< [Output] velocities
      logical geom_check_handle, geom_check_cent
      external geom_check_handle, geom_check_cent
      integer i, ncent
c
      geom_coords_vels_get =
     $     geom_check_handle(geom, 'geom_coords_vels_get')
      if (.not. geom_coords_vels_get) return
c
      ncent = ncenter(geom)
      do i = 1, ncent
         c(1,i) = coords(1,i,geom)
         c(2,i) = coords(2,i,geom)
         c(3,i) = coords(3,i,geom)
         v(1,i) = velocities(1,i,geom)
         v(2,i) = velocities(2,i,geom)
         v(3,i) = velocities(3,i,geom)
      end do
c
      end
c
C> \brief Defines the coordinates and velocities of all centers in a geometry
c
C> Defines the coordinates and velocities of all centers in a geometry assuming
C> that the caller has previously defined how many centers there
C> are (see e.g. geom_ncent_set).
c
C> \return Returns .true. if the function was successfull,
C> and .false. otherwise.
      logical function geom_coords_vels_set(geom, c, v)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
c
      integer geom             !< [Input] the geometry handle
      double precision c(3, *) !< [Input] coordinates
      double precision v(3, *) !< [Input] velocities
      logical geom_check_handle, geom_check_cent
      external geom_check_handle, geom_check_cent
      integer i, ncent
c
      geom_coords_vels_set =
     $     geom_check_handle(geom, 'geom_coords_vels_set')
      if (.not. geom_coords_vels_set) return
c
      ncent = ncenter(geom)
      do i = 1, ncent
         coords(1,i,geom) = c(1,i)
         coords(2,i,geom) = c(2,i)
         coords(3,i,geom) = c(3,i)
         velocities(1,i,geom) = v(1,i)
         velocities(2,i,geom) = v(2,i)
         velocities(3,i,geom) = v(3,i)
      end do
c
      end
c
C> \brief Look up data of a specific atom
c
C> Extracts the data, such as the tag, coordinates and charge, of a
C> specific center of a geometry instance. The center of interest is
C> given by the rank of the center.
c
C> \return Returns .true. if the function was successful, .false. otherwise.
      logical function geom_cent_get(geom, icent, t, c, q)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
c     
      integer geom              !< [Input] the geometry handle
      integer icent             !< [Input] the center rank
      character*16 t            !< [Output] the center tag
      double precision c(3)     !< [Output] the center coordinates
      double precision q        !< [Output] the center charge
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
c
C> \brief Set the data of a specific atom
c
C> Defines the data, such as the tag, coordinates and charge, of a
C> specific center of a geometry instance. The center of interest is
C> given by the rank of the center.
c
C> \return Returns .true. if the function was successful, .false. otherwise.
      logical function geom_cent_set(geom, icent, t, c, q)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
c     
      integer geom              !< [Input] the geometry handle
      integer icent             !< [Input] the center rank
      character*16 t            !< [Input] the center tag
      double precision c(3)     !< [Input] the center coordinates
      double precision q        !< [Input] the center charge
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
c ***** commented  out by EJB, Please do not uncomment w/o talking to Eric ******
c      call geom_compute_values(geom)
c ***** commented  out by EJB, Please do not uncomment w/o talking to Eric ******
c
      end
c
C> \brief Look up data of a specific atom, including velocity
c
C> Extracts the data, such as the tag, coordinates, velocity and charge,
C> of a specific center of a geometry instance. The center of interest
C> is given by the rank of the center.
c
C> \return Returns .true. if the function was successful, .false. otherwise.
c
      logical function geom_centv_get(geom, icent, t, c, v, q)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
c     
      integer geom          !< [Input] the geometry handle
      integer icent         !< [Input] the center rank
      character*16 t        !< [Output] the center tag
      double precision c(3) !< [Output] the center coordinates
      double precision v(3) !< [Output] the center velocity
      double precision q    !< [Output] the center charge
      logical geom_check_handle, geom_check_cent
      external geom_check_handle, geom_check_cent
c
      geom_centv_get = geom_check_handle(geom, 'geom_centv_get')
      if (.not. geom_centv_get) return
      geom_centv_get = geom_check_cent(geom, 'geom_centv_get', icent)
      if (.not. geom_centv_get) return

c
      t = tags(icent,geom)
      c(1) = coords(1,icent,geom)
      c(2) = coords(2,icent,geom)
      c(3) = coords(3,icent,geom)
      v(1) = velocities(1,icent,geom)
      v(2) = velocities(2,icent,geom)
      v(3) = velocities(3,icent,geom)
      q = charge(icent,geom)
      geom_centv_get = .true.
c
      end
c
C> \brief Store the data of a specific atom, including velocity
c
C> Stores the data, such as the tag, coordinates, velocity and charge,
C> of a specific center of a geometry instance. The center of interest
C> is given by the rank of the center.
c
C> \return Returns .true. if the function was successful, .false. otherwise.
c
      logical function geom_centv_set(geom, icent, t, c, v, q)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
c     
      integer geom          !< [Input] the geometry handle
      integer icent         !< [Input] the center rank
      character*16 t        !< [Input] the center tag
      double precision c(3) !< [Input] the center coordinates
      double precision v(3) !< [Input] the center velocity
      double precision q    !< [Input] the center charge
      logical geom_check_handle, geom_check_cent
      external geom_check_handle, geom_check_cent
c
      geom_centv_set = geom_check_handle(geom, 'geom_centv_set')
      if (.not. geom_centv_set) return
      geom_centv_set = geom_check_cent(geom, 'geom_centv_set', icent)
      if (.not. geom_centv_set) return
c
      tags(icent,geom) = t
      coords(1,icent,geom) = c(1)
      coords(2,icent,geom) = c(2)
      coords(3,icent,geom) = c(3)
      velocities(1,icent,geom) = v(1)
      velocities(2,icent,geom) = v(2)
      velocities(3,icent,geom) = v(3)
      charge(icent,geom) = q
c
c     compute effective nuclear repulsion energy, dipole and
c     interaction with external fields
c
      call geom_compute_values(geom)
c
      end
c
C> \brief Get the number of centers of a geometry
c
C> Query a geometry for the number of centers in it.
C> The function returns .true. if geom holds a valid handle,
C> it returns .false. otherwise.
      logical function geom_ncent(geom, ncent)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
c     
      integer geom              !< [input] the geometry handle
      integer ncent             !< [output] the number of centers
      logical geom_check_handle, geom_check_cent
      external geom_check_handle, geom_check_cent
c
      geom_ncent = geom_check_handle(geom, 'geom_ncent')
      if (.not. geom_ncent) return
      ncent = ncenter(geom)
c      
      end
c
C> \brief Set the number of centers of a geometry
c
C> Define the number of centers in a geometry.
C> The function returns .true. if geom holds a valid handle,
C> it returns .false. otherwise.
      logical function geom_ncent_set(geom, ncent)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
c     
      integer geom              !< [input] the geometry handle
      integer ncent             !< [input] the number of centers
      logical geom_check_handle, geom_check_cent
      external geom_check_handle, geom_check_cent
c
      geom_ncent_set = geom_check_handle(geom, 'geom_ncent_set')
      if (.not. geom_ncent_set) return
      ncenter(geom) = ncent
c      
      end
c
C> \brief Extracts the number of symmetry unique centers
c
C> This function extracts the number of symmetry unique centers from
C> a geometry instance.
c
C> \return Returns .true. if successfull, and .false. otherwise.
c
      logical function geom_ncent_unique(geom, ncent)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
c     
      integer geom              !< [Input] the geometry handle
      integer ncent             !< [Output] the number of unique centers
      logical geom_check_handle, geom_check_cent
      external geom_check_handle, geom_check_cent
c
      geom_ncent_unique = geom_check_handle(geom, 'geom_ncent_unique')
      if (.not. geom_ncent_unique) return
      ncent = ncenter_unique(geom)
c      
      end
c
C> \brief Checks whether a given center is a point charge
c
C> \return Returns .true. if the center is a point charge, and .false.
C> otherwise.
c
      logical function geom_isbq(geom, icent)
      implicit none
#include "nwc_const.fh"
#include "inp.fh"
#include "geomP.fh"
c     
      integer geom              !< [Input] the geometry handle
      integer icent             !< [Input] the center rank
      logical status
      character*16 tag          
      logical geom_check_handle, geom_check_cent
      external geom_check_handle, geom_check_cent
c
      status = geom_check_handle(geom, 'geom_cent_tag')
      if (.not. status) then
      call errquit("no geometry handle",0,0)
      end if
      status = geom_check_cent(geom, 'geom_cent_tag', icent)
      if (.not. status) then
      call errquit("no geometry center",0,0)
      end if
c
      tag = tags(icent,geom)
      geom_isbq = inp_compare(0,tag,'bq')
c
      end
c
C> \brief Look the tag of a specific center up
c
C> Extracts the tag of a specified center from the geometry instance.
c
C> \return Return .true. if the tag was found, and .false. otherwise.
c
      logical function geom_cent_tag(geom, icent, tag)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
c     
      integer geom              !< [Input] the geometry handle
      integer icent             !< [Input] the center rank
      character*16 tag          !< [Output] the center tag
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
#include "errquit.fh"
#include "nwc_const.fh"
#include "geomP.fh"
c     
      integer geom              ! [input]
      double precision ef       ! [input]
c
      call errquit('geom_efield_set: not yet!', 0, GEOM_ERR)
c     call geom_set_values(geom)
      geom_efield_set = .false.
      end
      logical function geom_efield_get(geom, ef)
      implicit none
#include "nwc_const.fh"
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
c
C> \brief Print the geometry in XYZ + charge format
c
C> Write the specified geometry in ASCII to the specified file. The
C> geometry is written in the usual XYZ format but with explicit
C> charges added. This also implies converting
C> the units of the coordinates to Angstrom.
c
C> \return Returns .true. if the geometry was successfully printed,
C> and .false. otherwise.
c
      logical function geom_print_xyzq(geom, unit)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
#include "util.fh"
#include "inp.fh"
#include "stdio.fh"
      integer geom !< [Input] the geometry handle
      integer unit !< [Input] the file unit number
      integer j, icent
      double precision scale
      logical geom_check_handle
      external geom_check_handle
c
      geom_print_xyzq = .true.
      if (.not. geom_check_handle(geom, 'geom_print_xyzq')) then
         geom_print_xyzq = .false.
         return
      end if
c
      scale = 1.0d0 / angstrom_to_au
c
      do icent = 1, ncenter(geom)
       if(inp_compare(0,tags(icent,geom),'bq')) then
       write(unit,3) tags(icent,geom),
     $                (coords(j,icent,geom)*scale,j=1,3),
     $                charge(icent,geom)


 3       format(1x,a16,1x,3f15.8,3x,"charge",3x,f15.8)
      else
       write(unit,4) tags(icent,geom),
     $                (coords(j,icent,geom)*scale,j=1,3)


 4       format(1x,a16,1x,3f15.8)

      end if
      end do
c
      end
c
C> \brief Print a geometry in PDB format 
c
C> Write a geometry as ASCII and in PDB format to the specified
C> output unit.
c
C> \return Return .true. if successfull, and .false. otherwise.
c
      logical function geom_print_pdb(geom, unit)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
#include "util.fh"
#include "inp.fh"
#include "stdio.fh"
      integer geom !< [Input] the geometry handle
      integer unit !< [Input] the file unit number
      integer j, icent
      double precision scale
      logical geom_check_handle
      external geom_check_handle
c     
      geom_print_pdb = .true.
      if (.not. geom_check_handle(geom, 'geom_print_pdb')) then
         geom_print_pdb = .false.
         return
      end if
c
      scale = 1.0d0 / angstrom_to_au
c
      write(unit,1) 
 1    format("####",T11,"id",T13,"name",
     >        T38,"x",T46,"y",T54,"z",T57,"charge")

      do icent = 1, ncenter(geom)

        write(unit,3) icent,tags(icent,geom), 
     $                (coords(j,icent,geom)*scale,j=1,3),
     $                 charge(icent,geom)


 3      format("ATOM",T7,I5,T13,A4,T31,F8.3,T39,F8.3,T47,F8.3,T55,F6.2)
      end do
c
      end
c
C> \brief Print the geometry in XYZ format
c
C> Write the specified geometry in ASCII to the specified file. The
C> geometry is written in the usual XYZ format. This implies converting
C> the units of the coordinates to Angstrom.
c
C> \return Returns .true. if the geometry was successfully printed,
C> and .false. otherwise.
c
      logical function geom_print_xyz(geom, unit)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
#include "util.fh"
#include "inp.fh"
#include "stdio.fh"
      integer geom !< [Input] the geometry handle
      integer unit !< [Input] the file unit number
      integer j, icent
      double precision scale
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

cc EJB commented this out 
cc     Convert from cartesian to crystallographic coordinates
cc
c         do i = 1, 3
c            tmp(i) = 0.0d0
c            do j = 1, 3
c               tmp(i) = tmp(i) + 
c     $              amatrix_inv(i,j,geom)*coords(j,icent,geom)
c            end do
c            tmp(i) = tmp(i)*scale ! Scale to angstrom
c         end do
c         write(unit,3) tags(icent,geom), (tmp(j),j=1,3)
        write(unit,3) tags(icent,geom), 
     $                (coords(j,icent,geom)*scale,j=1,3)


 3       format(1x,a16,1x,3f15.8)
      end do
c
      end
      logical function mol_geom_print_xyz(geom, unit, energy)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
#include "util.fh"
#include "inp.fh"
#include "stdio.fh"
      integer geom, unit
      integer j, icent
      double precision scale, energy
      logical geom_check_handle
      external geom_check_handle
c     
      mol_geom_print_xyz = .true.
      if (.not. geom_check_handle(geom, 'mol_geom_print_xyz')) then
         mol_geom_print_xyz = .false.
         return
      end if
c
      scale = 1.0d0 / angstrom_to_au
c
      write(unit,1) ncenter(geom)
 1    format(1x,i5)
      write(unit,2) energy
 2    format(1x,f15.8)
      do icent = 1, ncenter(geom)
        write(unit,3) tags(icent,geom), 
     $                (coords(j,icent,geom)*scale,j=1,3)


 3       format(1x,a16,1x,3f15.8)
      end do
c
      end
c
C> \brief Print the contents of a geometry instance
c
C> Prints the contents of a geometry instance, irrespective
C> of whether it is a molecule or a crystal structure.
C> The output is always provided on standard output.
c
C> \return Return .true. if the structure could be printed,
C> and .false. otherwise.
      logical function geom_print(geom)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "geomP.fh"
#include "util.fh"
#include "inp.fh"
#include "stdio.fh"
c     
c     Basic printing of cartesian geometry 
c     needs support for internal coords, different formats, ...
c     
      integer geom              !< [Input] the geometry handle
      integer icent, jcent
      integer i
      double precision scale, tmp(3),twopi
      character*80 buf
      logical oprint_uniq,ofinite,oprint_crystal
c
c     external functions
      logical geom_check_handle, geom_check_cent, geom_get_user_scale,
     $     geom_print_zmatrix, geom_any_finuc
      external geom_check_handle, geom_check_cent, geom_get_user_scale,
     $     geom_print_zmatrix, geom_any_finuc
      double precision deter3
      external         deter3
      double precision ddot
      external ddot
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

       oprint_crystal = (isystype(geom).ne.0)
c      if (isystype(geom) .ne. 0) then 
c        geom_print = .true.
c        return
c      endif

c     
      if (.not. geom_get_user_scale(geom, scale))
     $     call errquit('geom_print: user units?',0, GEOM_ERR)
c     
      buf = ' '
      write(buf,1) 'Geometry',
     $        names(geom)(1:lenn(geom)), 
     $        trans(geom)(1:lent(geom))
 1    format(a,' "',a,'" -> "',a,'"')
      write(LuOut,*)
      write(LuOut,*)
      call util_print_centered(LuOut,buf,40,.true.)
      write(LuOut,*)
      write(LuOut,2) user_units(geom)(1:inp_strlen(user_units(geom))),
     $     scale
 2    format(' Output coordinates in ', a,
     $     ' (scale by ',f12.9,' to convert to a.u.)')
      if (include_bqbq(geom))
     $     write(LuOut,*) ' Include Bq-Bq interactions'
c
      write(LuOut,*)
c     
      write(LuOut,3) 
 3    format('  No.       Tag          Charge          X',
     $     '              Y              Z'/
     $     ' ---- ---------------- ---------- --------------',
     $     ' -------------- --------------')
      do icent = 1, ncenter(geom)
         do i = 1, 3
            tmp(i) = coords(i,icent,geom)/scale ! Scale units as necessary
         end do
         write(LuOut,4) icent, tags(icent,geom), charge(icent,geom),
     $        (tmp(i),i=1,3)
 4       format(' ',i4,' ',a16,' ',f10.4,3f15.8)
      end do
c
      if (ddot(3*ncenter(geom),velocities(1,1,geom),1,
     $     velocities(1,1,geom),1) .gt. 1d-10) then

         write(LuOut,*)
         write(LuOut,*)
         call util_print_centered(LuOut,'Velocities',40,.true.)
         write(LuOut,3) 
         do icent = 1, ncenter(geom)
            write(LuOut,4) icent, tags(icent,geom), charge(icent,geom),
     $           (velocities(i,icent,geom),i=1,3)
         end do
      endif
c
c     print out lattice parameters 
c
      if (oprint_crystal) then
        write(LuOut,*)
        write(LuOut,*) '     Lattice Parameters '
        write(LuOut,*) '     ------------------ '
        write(LuOut,*)
        write(LuOut,5) user_units(geom)(1:inp_strlen(user_units(geom))),
     >     scale
 5    format('      lattice vectors in ', a,
     $     ' (scale by ',f12.9,' to convert to a.u.)')
        write(LuOut,*)
        write(LuOut,1241) amatrix(1,1,geom)/scale,
     >                    amatrix(2,1,geom)/scale,
     >                    amatrix(3,1,geom)/scale
        write(LuOut,1242) amatrix(1,2,geom)/scale,
     >                    amatrix(2,2,geom)/scale,
     >                    amatrix(3,2,geom)/scale
        write(LuOut,1243) amatrix(1,3,geom)/scale,
     >                    amatrix(2,3,geom)/scale,
     >                    amatrix(3,3,geom)/scale

        write(LuOut,1232) lattice_vectors(1,geom),
     >                    lattice_vectors(2,geom),
     >                    lattice_vectors(3,geom),
     >                    lattice_angles(1,geom),
     >                    lattice_angles(2,geom),
     >                    lattice_angles(3,geom)
        write(LuOut,1231) deter3(amatrix(1,1,geom))/(scale**3)

      write(LuOut,*)
      write(LuOut,6) 
 6    format('      reciprocal lattice vectors in a.u.')
      write(LuOut,*)
        twopi = 8.0d0*datan(1.0d0)
        write(LuOut,1244) amatrix_inv(1,1,geom)*twopi,
     >                    amatrix_inv(1,2,geom)*twopi,
     >                    amatrix_inv(1,3,geom)*twopi
        write(LuOut,1245) amatrix_inv(2,1,geom)*twopi,
     >                    amatrix_inv(2,2,geom)*twopi,
     >                    amatrix_inv(2,3,geom)*twopi
        write(LuOut,1246) amatrix_inv(3,1,geom)*twopi,
     >                    amatrix_inv(3,2,geom)*twopi,
     >                    amatrix_inv(3,3,geom)*twopi


      end if

c     
c     Only print out the masses for unique tags ... the structure
c     should actually only store the data for unique tags.
c     Also, keep all common output within 80 columns
c
      ofinite = geom_any_finuc(geom)
      write(LuOut,*)
      if (ofinite) then
        write(LuOut,*) '     Atomic Mass and Nuclear Exponent '
        write(LuOut,*) '     -------------------------------- '
      else
        write(LuOut,*) '     Atomic Mass '
        write(LuOut,*) '     ----------- '
      end if
      write(LuOut,*)
      do icent = 1, ncenter(geom)
        if (abs(geom_mass(icent,geom)).lt.1.0d-07) goto 765
        do jcent = 1, icent-1
          if (tags(icent,geom) .eq. tags(jcent,geom)) goto 765
        enddo
        if (geom_invnucexp(icent,geom) .gt. 0.0d0) then
          write(LuOut,43) tags(icent,geom), geom_mass(icent,geom),
     &        1.0d0/geom_invnucexp(icent,geom)
   43     format('      ',a16,' ',f10.6,1pe20.6)
        else
          write(LuOut,44) tags(icent,geom), geom_mass(icent,geom)
   44     format('      ',a16,' ',f10.6)
        end if
765     continue
      enddo

      write(LuOut,*)
      if (.not.oprint_crystal) then
c     
      write(LuOut,41) erep(geom)
   41 format(/' Effective nuclear repulsion energy (a.u.) ', f18.10/)
c     
      write(LuOut,91)
 91   format('            Nuclear Dipole moment (a.u.) ')
      write(LuOut,101)
 101  format('            ----------------------------')
      write(LuOut,7)
 7    format('        X                 Y               Z'/
     $     ' ---------------- ---------------- ----------------')
      write(LuOut,8) (ndipole(i,geom), i=1,3)
 8    format(3(1x,f16.10))
      end if
      write(LuOut,*)
c     
      oprint_uniq = sym_num_ops(geom) .gt. 0
      if (oprint_uniq) then
         call sym_print_all(geom,.true.,oprint_uniq,.false.,
     >                      .false.,.false.)
      endif
c     
      if (zmt_source(geom) .ne. ' ' .and. 
     $     util_print('geomzmat',print_none)) then
c        JEM: Must pass an array, not a scalar
         geom_print = geom_print_zmatrix(geom,(/0.d0/),' ',.false.)
      else
         geom_print = .true.
      endif
c     
      return
 1231 FORMAT(5x,' omega=',f8.1)
 1232 FORMAT(5x,' a=    ',f8.3,' b=   ',f8.3,' c=    ',f8.3,
     >      /5x,' alpha=',f8.3,' beta=',f8.3,' gamma=',f8.3)
 1241 FORMAT(5x,' a1=<',3f8.3,' >')
 1242 FORMAT(5x,' a2=<',3f8.3,' >')
 1243 FORMAT(5x,' a3=<',3f8.3,' >')
 1244 FORMAT(5x,' b1=<',3f8.3,' >')
 1245 FORMAT(5x,' b2=<',3f8.3,' >')
 1246 FORMAT(5x,' b3=<',3f8.3,' >')

      end

*     **************************************************
*     *                                                *
*     *            geom_use_primitive                  *
*     *                                                *
*     **************************************************
      logical function geom_use_primitive(geom)
      implicit none
      integer geom  
#include "nwc_const.fh"
#include "geomP.fh"
      geom_use_primitive = use_primitive(geom)
      return
      end

*     **************************************************
*     *                                                *
*     *            geom_primitive_center               *
*     *                                                *
*     **************************************************
      character*1 function geom_primitive_center(geom)
      implicit none
      integer geom  
#include "nwc_const.fh"
#include "geomP.fh"
      geom_primitive_center = primitive_center(geom)
      return
      end


*     **************************************************
*     *                                                *
*     *            geom_is_conventional                *
*     *                                                *
*     **************************************************

      logical function geom_is_conventional(geom)
      implicit none
      integer geom

#include "errquit.fh"
#include "nwc_const.fh"
#include "geomP.fh"
#include "util.fh"
#include "inp.fh"
#include "stdio.fh"

      logical value,is_convention
      real*8 lat(6),radtodeg

      integer grp_num,crystal
      integer Triclinic,Monoclinic,Orthorhombic
      integer Tetragonal,Trigonal,Hexagonal,Cubic
      parameter (Triclinic=1,Monoclinic=2,Orthorhombic=3)
      parameter (Tetragonal=4,Trigonal=5,Hexagonal=6,Cubic=7)

      logical  geom_lattice_get
      external geom_lattice_get

      radtodeg = 180.0d0/(4.0d0*datan(1.0d0))

      value = geom_lattice_get(geom,lat)
      lat(4) = lat(4)*radtodeg
      lat(5) = lat(5)*radtodeg
      lat(6) = lat(6)*radtodeg
      grp_num = group_number(geom)

      if (grp_num.lt.3)                         crystal = Triclinic
      if ((grp_num.ge.3  ).and.(grp_num.lt.16 ))crystal = Monoclinic
      if ((grp_num.ge.16 ).and.(grp_num.lt.75 ))crystal = Orthorhombic
      if ((grp_num.ge.75 ).and.(grp_num.lt.143))crystal = Tetragonal
      if ((grp_num.ge.143).and.(grp_num.lt.168))crystal = Trigonal
      if ((grp_num.ge.168).and.(grp_num.lt.195))crystal = Hexagonal
      if ((grp_num.ge.195).and.(grp_num.lt.231))crystal = Cubic
      is_convention = .true.
      if (crystal.eq.Triclinic) then
         is_convention=is_convention.and.(dabs(lat(4)-90.0d0).gt.1.0d-3)
         is_convention=is_convention.and.(dabs(lat(5)-90.0d0).gt.1.0d-3)
         is_convention=is_convention.and.(dabs(lat(6)-90.0d0).gt.1.0d-3)
         is_convention=is_convention.and.(dabs(lat(1)-lat(2)).gt.1.0d-3)
         is_convention=is_convention.and.(dabs(lat(1)-lat(3)).gt.1.0d-3)
         is_convention=is_convention.and.(dabs(lat(2)-lat(3)).gt.1.0d-3)
      else if (crystal.eq.Monoclinic) then
         is_convention=is_convention.and.(dabs(lat(4)-90.0d0).lt.1.0d-3)
         is_convention=is_convention.and.(dabs(lat(5)-90.0d0).gt.1.0d-3)
         is_convention=is_convention.and.(dabs(lat(6)-90.0d0).lt.1.0d-3)
         is_convention=is_convention.and.(dabs(lat(1)-lat(2)).gt.1.0d-3)
         is_convention=is_convention.and.(dabs(lat(1)-lat(3)).gt.1.0d-3)
         is_convention=is_convention.and.(dabs(lat(2)-lat(3)).gt.1.0d-3)
      else if (crystal.eq.Orthorhombic) then
         is_convention=is_convention.and.(dabs(lat(4)-90.0d0).lt.1.0d-3)
         is_convention=is_convention.and.(dabs(lat(5)-90.0d0).lt.1.0d-3)
         is_convention=is_convention.and.(dabs(lat(6)-90.0d0).lt.1.0d-3)
         is_convention=is_convention.and.(dabs(lat(1)-lat(2)).gt.1.0d-3)
         is_convention=is_convention.and.(dabs(lat(1)-lat(3)).gt.1.0d-3)
         is_convention=is_convention.and.(dabs(lat(2)-lat(3)).gt.1.0d-3)
      else if (crystal.eq.Tetragonal) then
         is_convention=is_convention.and.(dabs(lat(4)-90.0d0).lt.1.0d-3)
         is_convention=is_convention.and.(dabs(lat(5)-90.0d0).lt.1.0d-3)
         is_convention=is_convention.and.(dabs(lat(6)-90.0d0).lt.1.0d-3)
         is_convention=is_convention.and.(dabs(lat(1)-lat(2)).lt.1.0d-3)
         is_convention=is_convention.and.(dabs(lat(1)-lat(3)).gt.1.0d-3)
      else if (crystal.eq.Trigonal) then
         is_convention=is_convention.and.(dabs(lat(4)-90.0d0).lt.1.0d-3)
         is_convention=is_convention.and.(dabs(lat(5)-90.0d0).lt.1.0d-3)
         is_convention=is_convention.and.(dabs(lat(6)-120.0d0).lt.1.d-3)
         is_convention=is_convention.and.(dabs(lat(1)-lat(2)).lt.1.0d-3)
      else if (crystal.eq.Hexagonal) then
         is_convention=is_convention.and.(dabs(lat(4)-90.0d0).lt.1.0d-3)
         is_convention=is_convention.and.(dabs(lat(5)-90.0d0).lt.1.0d-3)
         is_convention=is_convention.and.(dabs(lat(6)-120.0d0).lt.1.d-3)
         is_convention=is_convention.and.(dabs(lat(1)-lat(2)).lt.1.0d-3)
      else if (crystal.eq.Cubic) then
         is_convention=is_convention.and.(dabs(lat(4)-90.0d0).lt.1.0d-3)
         is_convention=is_convention.and.(dabs(lat(5)-90.0d0).lt.1.0d-3)
         is_convention=is_convention.and.(dabs(lat(6)-90.0d0).lt.1.0d-3)
         is_convention=is_convention.and.(dabs(lat(1)-lat(2)).lt.1.0d-3)
         is_convention=is_convention.and.(dabs(lat(1)-lat(3)).lt.1.0d-3)
         is_convention=is_convention.and.(dabs(lat(2)-lat(3)).lt.1.0d-3)
      end if

      geom_is_conventional = is_convention
      return
      end 

*     *****************************************************
*     *                                                   *
*     *            geom_convert_to_primitive              *
*     *                                                   *
*     *****************************************************
*
*   This routine converts conventional centered cells to primitive cells
*
*   The conventions used to relate the centered cells (abc) to the 
*  primitive cells (a'b'c') are as follows:
*
*    A-centered: The author's (DGC)
*                a'=a, b'=(b+c)/2, c'=(-b+c)/2
*    C-centered: Orthohexagonal cell convention (from Int. Tables, Fig. 5.8, p.70, cell #C2)
*                a'=(a+b)/2, b'=(-a+b)/2, c'=c
*    R-centered: Obverse Int. Tables, Table 5.1, p.78 - Cell #1
*                See also Figs. 5.7a & c, p.79)
*                   a' = 2a/3 + b/3 + c/3
*                   b' = -a/3 + b/3 + c/3
*                   c' = -a/3 -2b/3 + c/3
*    I-centered: Int. Tables, Table 5.1, p.76. See also Fig. 5.4, p.77.
*                   a' = -a/2 + b/2 + c/2
*                   b' =  a/2 - b/2 + c/2
*                   c' =  a/2 + b/2 - c/2
*    F-centered: Int. Tables, Table 5.1, p.77. See also Fig. 5.5, p.77.
*                a' = (b+c)/2, b' = (a+c)/2, c' = (a+b)/2
*
*
c---> The following is a list of the centered 3D groups:
c
c      A-centered (4 total): 38, 39, 40, 41 - orthorhombic
c
c      C-centered (16 total): 5, 8, 9, 12, 15, 20, 21, 35, 36, 37, 63,  - mono and ortho
c                              64, 65, 66, 67, 68
c
c      R-centered (7 total):146,148,155,160,161,166,167  - trig
c
c      I-centered (38 total):
c          23, 24, 44, 45, 46, 71, 72, 73, 74, 79, 80, 82, 87, 88, 97, 98, - ortho,tetra,cubic
c           107,108,109,110,119,120,121,122,139,140,141,142,197,199,204,206,
c           211,214,217,220,229,230
c
c      F-centered (16 total)
c            22, 42, 43, 69, 70,196,202,203,209,210,216,219,225,226,227,228 - ortho,cubic
c
      logical function geom_convert_to_primitive(geom,ctr_type,tf)
      implicit none
      integer geom
      character*1 ctr_type
      real*8 tf(6,3)

#include "errquit.fh"
#include "nwc_const.fh"
#include "geomP.fh"
#include "util.fh"
#include "inp.fh"
#include "stdio.fh"

*     **** local variables ****
      logical value,update_lattice,is_conventional
      integer i,ncent
      real*8 lat(6),degtorad,radtodeg
      real*8 amat(3,3),amat2(3,3),c(3)

      integer grp_num,crystal
      integer Triclinic,Monoclinic,Orthorhombic
      integer Tetragonal,Trigonal,Hexagonal,Cubic
      parameter (Triclinic=1,Monoclinic=2,Orthorhombic=3)
      parameter (Tetragonal=4,Trigonal=5,Hexagonal=6,Cubic=7)

*     **** external functions ****
      logical  geom_lattice_get,geom_lattice_set,geom_is_conventional
      external geom_lattice_get,geom_lattice_set,geom_is_conventional
      logical  geom_amatrix_get,geom_amatrix_set
      external geom_amatrix_get,geom_amatrix_set

      value = .true.
      if ((ctr_type.eq.'F').or.
     >   (ctr_type.eq.'A').or.
     >   (ctr_type.eq.'C').or.
     >   (ctr_type.eq.'R').or.
     >   (ctr_type.eq.'I')) then

      update_lattice = .False.
      degtorad = 4.0d0*datan(1.0d0)/180.0d0
      radtodeg = 180.0d0/(4.0d0*datan(1.0d0))
      value = geom_amatrix_get(geom,amat)
      value = value.and.geom_amatrix_get(geom,amat2)
      value = value.and.geom_lattice_get(geom,lat)
      lat(4) = lat(4)*radtodeg
      lat(5) = lat(5)*radtodeg
      lat(6) = lat(6)*radtodeg

      is_conventional = geom_is_conventional(geom)

c     *** F-centered: Int. Tables, Table 5.1, p.77. See also Fig. 5.5, p.77. ***
c     *** a' = (b+c)/2, b' = (a+c)/2, c' = (a+b)/2 ***
      if (ctr_type.eq.'F') then
         amat(1,1) = 0.50d0*(amat2(1,2)+amat2(1,3))
         amat(2,1) = 0.50d0*(amat2(2,2)+amat2(2,3))
         amat(3,1) = 0.50d0*(amat2(3,2)+amat2(3,3))

         amat(1,2) = 0.50d0*(amat2(1,1)+amat2(1,3))
         amat(2,2) = 0.50d0*(amat2(2,1)+amat2(2,3))
         amat(3,2) = 0.50d0*(amat2(3,1)+amat2(3,3))

         amat(1,3) = 0.50d0*(amat2(1,1)+amat2(1,2))
         amat(2,3) = 0.50d0*(amat2(2,1)+amat2(2,2))
         amat(3,3) = 0.50d0*(amat2(3,1)+amat2(3,2))
         update_lattice = is_conventional
         primitive_center(geom) = ctr_type

      end if

c     *** A-centered: The author's (DGC) ***
c     *** a'=a, b'=(b+c)/2, c'=(-b+c)/2 ***
      if (ctr_type.eq.'A') then
         amat(1,2) = 0.50d0*(amat2(1,2)+amat2(1,3))
         amat(2,2) = 0.50d0*(amat2(2,2)+amat2(2,3))
         amat(3,2) = 0.50d0*(amat2(3,2)+amat2(3,3))

         amat(1,3) = 0.50d0*(-amat2(1,2)+amat2(1,3))
         amat(2,3) = 0.50d0*(-amat2(2,2)+amat2(2,3))
         amat(3,3) = 0.50d0*(-amat2(3,2)+amat2(3,3))
         update_lattice = is_conventional
         primitive_center(geom) = ctr_type
      end if

c     *** C-centered: Orthohexagonal cell convention (from Int. Tables, Fig. 5.8, p.70, cell #C2) ***
c     *** a'=(a+b)/2, b'=(-a+b)/2, c'=c ***
      if (ctr_type.eq.'C') then
         amat(1,1) = 0.50d0*(amat2(1,1)+amat2(1,2))
         amat(2,1) = 0.50d0*(amat2(2,1)+amat2(2,2))
         amat(3,1) = 0.50d0*(amat2(3,1)+amat2(3,2))

         amat(1,2) = 0.50d0*(-amat2(1,1)+amat2(1,2))
         amat(2,2) = 0.50d0*(-amat2(2,1)+amat2(2,2))
         amat(3,2) = 0.50d0*(-amat2(3,1)+amat2(3,2))
         update_lattice = is_conventional
         primitive_center(geom) = ctr_type
      end if

c     *** R-centered: Obverse Int. Tables, Table 5.1, p.78 - Cell #1 ***
c     *** See also Figs. 5.7a & c, p.79) ***
c     *** a' = 2a/3 + b/3 + c/3 ***
c     *** b' = -a/3 + b/3 + c/3 ***
c     *** c' = -a/3 -2b/3 + c/3 ***
      if (ctr_type.eq.'R') then
         amat(1,1) = (2.0d0*amat2(1,1)+amat2(1,2)+amat2(1,3))/3.0d0
         amat(2,1) = (2.0d0*amat2(2,1)+amat2(2,2)+amat2(2,3))/3.0d0
         amat(3,1) = (2.0d0*amat2(3,1)+amat2(3,2)+amat2(3,3))/3.0d0

         amat(1,2) = (-amat2(1,1)+amat2(1,2)+amat2(1,3))/3.0d0
         amat(2,2) = (-amat2(2,1)+amat2(2,2)+amat2(2,3))/3.0d0
         amat(3,2) = (-amat2(3,1)+amat2(3,2)+amat2(3,3))/3.0d0

         amat(1,3) = (-amat2(1,1)-2.0d0*amat2(1,2)+amat2(1,3))/3.0d0
         amat(2,3) = (-amat2(2,1)-2.0d0*amat2(2,2)+amat2(2,3))/3.0d0
         amat(3,3) = (-amat2(3,1)-2.0d0*amat2(3,2)+amat2(3,3))/3.0d0
         update_lattice = is_conventional
         primitive_center(geom) = ctr_type
      end if

c     *** I-centered: Int. Tables, Table 5.1, p.76. See also Fig. 5.4, p.77. ***
c     *** a' = -a/2 + b/2 + c/2 ***
c     *** b' =  a/2 - b/2 + c/2 ***
c     *** c' =  a/2 + b/2 - c/2 ***
      if (ctr_type.eq.'I') then
         amat(1,1) = 0.5d0*(-amat2(1,1)+amat2(1,2)+amat2(1,3))
         amat(2,1) = 0.5d0*(-amat2(2,1)+amat2(2,2)+amat2(2,3))
         amat(3,1) = 0.5d0*(-amat2(3,1)+amat2(3,2)+amat2(3,3))

         amat(1,2) = 0.5d0*(amat2(1,1)-amat2(1,2)+amat2(1,3))
         amat(2,2) = 0.5d0*(amat2(2,1)-amat2(2,2)+amat2(2,3))
         amat(3,2) = 0.5d0*(amat2(3,1)-amat2(3,2)+amat2(3,3))

         amat(1,3) = 0.5d0*(amat2(1,1)+amat2(1,2)-amat2(1,3))
         amat(2,3) = 0.5d0*(amat2(2,1)+amat2(2,2)-amat2(2,3))
         amat(3,3) = 0.5d0*(amat2(3,1)+amat2(3,2)-amat2(3,3))
         update_lattice = is_conventional
         primitive_center(geom) = ctr_type
      end if

      if (update_lattice) then
c        *** convert the fractional coords ***
         ncent = ncenter(geom)
         do i = 1, ncent
            c(1) = coords(1,i,geom)
            c(2) = coords(2,i,geom)
            c(3) = coords(3,i,geom)
            coords(1,i,geom) = tf(1,1)*c(1)+tf(1,2)*c(2)+tf(1,3)*c(3)
            coords(2,i,geom) = tf(2,1)*c(1)+tf(2,2)*c(2)+tf(2,3)*c(3)
            coords(3,i,geom) = tf(3,1)*c(1)+tf(3,2)*c(2)+tf(3,3)*c(3)
         end do
         call xlattice_abc_abg(lat(1),lat(2),lat(3),
     >                         lat(4),lat(5),lat(6),amat)
       
          lat(4) = lat(4)*degtorad
          lat(5) = lat(5)*degtorad
          lat(6) = lat(6)*degtorad
         !value = value.and.geom_amatrix_set(geom,amat)
         value = value.and.geom_lattice_set(geom,lat)
      end if

      end if

      geom_convert_to_primitive = value
      return
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
c
C> \brief Convert an atom tag to a covalent radius
c
C> \return Returns .true. when successfull, and .false. otherwise.
c
      logical function geom_tag_to_covalent_radius(tag,radius)
      implicit none
#include "inp.fh"
c
c     Try to decode a tag and return the covalent radius (a.u.) for
c     the corresponding atom.
c
      character*16 tag          !< [Input] the atom tag
      double precision radius   !< [Output] the atom covalent radius
c
      character*2 symbol
      character*16 element, ttag
      integer atn
      logical geom_get_def_rcov, geom_tag_to_element
      external geom_get_def_rcov, geom_tag_to_element
c
      geom_tag_to_covalent_radius = .false.
c
      if (.not. geom_tag_to_element(tag, symbol, element, atn)) then
c
c        Is not an atom.  Try removing Bq or X.
c
         if (inp_compare(.false., tag(1:1), 'x')) then
            ttag = tag(2:)
         else if (inp_compare(.false., tag(1:2), 'bq')) then
            ttag = tag(3:)
         else
            return              ! Nothing recognizable
         endif
         if (.not. geom_tag_to_element(ttag, symbol, element, atn)) then
c
c           We found a "Bq" or "X" but it is not labeled with an element
c           (e.g. "XH" or "BqN") so we cannot associate any atomic
c           properties. 
c
            geom_tag_to_covalent_radius = .false.
            return
         else
c
c           We found a "Bq" or "X" with an atomic type indication and 
c           we now have the atomic number indicated, so we will go 
c           with that.
c
         endif
      else
c
c        We found an atom and now know its atomic number
c
      endif
c
c     atn should be set to something sensible
c
      geom_tag_to_covalent_radius = geom_get_def_rcov(atn, radius)
c
      end
c
C> \brief Converts an atom tag into a chemical element
c
C> Tags are names used to identify atoms in geometries. These tags allow
C> for considerable flexibility to specify centers, e.g. valid tags for
C> a Carbon atom are: Carbon, C, C14, C_alpha. However, in a chemistry
C> code at some point we just need to know that we are dealing with a
C> Carbon atom. This routine analyses a given tag and returns the
C> element it specifies.
c
C> \return Returns .true. if the chemical element could be established,
C> and .false. otherwise.
c
      logical function geom_tag_to_element(tag, symbol, element, atn)
      implicit none
#include "inp.fh"
#include "nwc_const.fh"
#include "geomP.fh"
      character*2 symbols(nelements)
      character*16 elements(nelements)
      character*16 tag          !< [Input] the tag, e.g. He232
      character*(*) symbol      !< [Output] the chemical symbol, e.g. He
      character*(*) element     !< [Output] the element, e.g. Helium
      integer atn               !< [Output] the nuclear charge, e.g. 2
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
     $     'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds',
     $     'Rg', 'Cn'/
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
     $     'Nobelium', 'Lawrencium','Rutherfordium','Dubnium',
     $     'Seaborgium','Bohrium','Hassium','Meitnerium',
     $     'Darmstadtium', 'Roentgenium', 'Copernicium'/
      data sym1/'H','B','C','N','O','F','P','S','K','V','Y','I','W','U'/
      data atn1/ 1 , 5 , 6 , 7 , 8 , 9 , 15, 16, 19, 23, 39, 53, 74, 92/

      geom_tag_to_element = .false.
c
c     eliminate conventions that refer to centers used for
c     computation purposes .. just bq and x for now
c
      buf = tag
      lbuf = inp_strlen(buf)
      if (lbuf .eq. 0) return
c
      call inp_lcase(buf)
      if (buf(1:2) .eq. 'bq' .or. 
     $     ((buf(1:1).eq.'x') .and. (buf(2:2).ne.'e'))) then ! X but not Xe
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
        do ind = 1,nelements
          if (inp_compare(.false.,buf(1:4),elements(ind)(1:4))) then
            symbol  = symbols(ind)
            element = elements(ind)
            atn     = ind
            geom_tag_to_element = .true.
            return
          endif
        enddo
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
      if (inp_match(14, .false., buf(2:2), sym1, ind)) then
         ind = atn1(ind)
         symbol  = symbols(ind)
         element = elements(ind)
         atn     = ind
         geom_tag_to_element = .true.
         return
      end if
cc
c     Nothing matched
c
      symbol = ' '
      element = ' '
      atn = 0
      return
c
      end
      function geom_tag_to_charge_gen(nt,ns,tag,q)
      implicit none
#include "inp.fh"
#include "nwc_const.fh"
#include "geomP.fh"
      logical geom_tag_to_charge_gen
      integer nt
      integer ns
      integer i
      logical match
      character*2 symbols(nelements)
      character*16 elements(nelements)
      character*1 tag(nt*ns)                  ! [input]
      double precision  q(nt)            ! [output]
c
c     attempt to figure out which element a tag refers to
c     and return the symbol, name and atomic no.
c
      integer j,offset
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
     $     'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds',
     $     'Rg', 'Cn'/
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
     $     'Nobelium', 'Lawrencium','Rutherfordium','Dubnium',
     $     'Seaborgium','Bohrium','Hassium','Meitnerium',
     $     'Darmstadtium', 'Roentgenium', 'Copernicium'/
      data sym1/'H','B','C','N','O','F','P','S','K','V','Y','I','W','U'/
      data atn1/ 1 , 5 , 6 , 7 , 8 , 9 , 15, 16, 19, 23, 39, 53, 74, 92/

      geom_tag_to_charge_gen = .false.
c
c     eliminate conventions that refer to centers used for
c     computation purposes .. just bq and x for now
c
      match = .false.
      do i=1,nt
        match = .false.
        offset = (i-1)*ns
        buf = " "
        do j=1,16
          buf(j:j) = tag(j+offset)
        end do
        lbuf = inp_strlen(buf)
        if (lbuf .eq. 0) goto 100 
c
        call inp_lcase(buf)
        if (buf(1:2) .eq. 'bq' .or. 
     $       ((buf(1:1).eq.'x') .and. (buf(2:2).ne.'e'))) then ! X but not Xe
           q(i)   = 0.0
           match = .true.
           goto 100
        end if
c
c       Attempt to match the first 4 characters of the
c       full names of the elements
c
        q(i) = 0.0
        if (lbuf .ge. 4) then
          do ind = 1,nelements
            if (inp_compare(.false.,buf(1:4),elements(ind)(1:4))) then
              q(i)    = dble(ind)
              match = .true.
              goto 100
            endif
          enddo
        end if
c
c       Failed ... attempt to match the first two characters
c       against two character element names
c
        if (buf(2:2) .ne. ' ') then
           if (inp_match(nelements,.false.,buf(1:2),symbols,ind)) then
              q(i)    = dble(ind)
              match = .true.
              goto 100
           end if
        end if
c
c      not Last ditch attempt ... match against 1 character symbols
c
        if (inp_match(14, .false., buf(1:1), sym1, ind)) then
           ind = atn1(ind)
           q(i)   = dble(ind)
           match = .true.
           goto 100
        end if
        if (inp_match(14, .false., buf(2:2), sym1, ind)) then
           ind = atn1(ind)
           q(i)   = dble(ind)
           match = .true.
           goto 100
        end if
100     continue
        if(.not.match) then
          write(*,*) "buffer",buf
          goto 101
        end if
      end do
101   continue
      
      geom_tag_to_charge_gen = match

      return
c
      end
      function geom_tag_to_charge(nt,tag,q)
      implicit none
#include "inp.fh"
#include "nwc_const.fh"
#include "geomP.fh"
      logical geom_tag_to_charge
      integer nt
      integer i
      logical match
      character*2 symbols(nelements)
      character*16 elements(nelements)
      character*16 tag(nt)                  ! [input]
      double precision  q(nt)            ! [output]
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
     $     'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds',
     $     'Rg', 'Cn'/
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
     $     'Nobelium', 'Lawrencium','Rutherfordium','Dubnium',
     $     'Seaborgium','Bohrium','Hassium','Meitnerium',
     $     'Darmstadtium', 'Roentgenium', 'Copernicium'/
      data sym1/'H','B','C','N','O','F','P','S','K','V','Y','I','W','U'/
      data atn1/ 1 , 5 , 6 , 7 , 8 , 9 , 15, 16, 19, 23, 39, 53, 74, 92/

      geom_tag_to_charge = .false.
c
c     eliminate conventions that refer to centers used for
c     computation purposes .. just bq and x for now
c
      match = .false.
      do i=1,nt
        match = .false.
        buf = tag(i)
        lbuf = inp_strlen(buf)
        if (lbuf .eq. 0) goto 100 
c
        call inp_lcase(buf)
        if (buf(1:2) .eq. 'bq' .or. 
     $       ((buf(1:1).eq.'x') .and. (buf(2:2).ne.'e'))) then ! X but not Xe
           q(i)   = 0.0
           match = .true.
           goto 100
        end if
c
c       Attempt to match the first 4 characters of the
c       full names of the elements
c
        q(i) = 0.0
        if (lbuf .ge. 4) then
          do ind = 1,nelements
            if (inp_compare(.false.,buf(1:4),elements(ind)(1:4))) then
              q(i)    = dble(ind)
              match = .true.
              goto 100
            endif
          enddo
        end if
c
c       Failed ... attempt to match the first two characters
c       against two character element names
c
        if (buf(2:2) .ne. ' ') then
           if (inp_match(nelements,.false.,buf(1:2),symbols,ind)) then
              q(i)    = dble(ind)
              match = .true.
              goto 100
           end if
        end if
c
c       Last ditch attempt ... match against 1 character symbols
c
        if (inp_match(14, .false., buf(1:1), sym1, ind)) then
           ind = atn1(ind)
           q(i)   = dble(ind)
           match = .true.
           goto 100
        end if
100     continue
        if(.not.match) goto 101
      end do
101   continue
      
      geom_tag_to_charge = match

      return
c
      end
      function geom_tag_to_atn(nt,tag,atn)
      implicit none
#include "inp.fh"
#include "nwc_const.fh"
#include "geomP.fh"
      logical geom_tag_to_atn
      integer nt
      integer i
      logical match
      character*2 symbols(nelements)
      character*16 elements(nelements)
c     TP: changed from character*(*) tag(nt) to character*(16) tag(nt)
      character*(16) tag(nt)         ! [input]
      integer atn(nt)            ! [output]
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
     $     'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds',
     $     'Rg', 'Cn'/
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
     $     'Nobelium', 'Lawrencium','Rutherfordium','Dubnium',
     $     'Seaborgium','Bohrium','Hassium','Meitnerium',
     $     'Darmstadtium', 'Roentgenium', 'Copernicium'/
      data sym1/'H','B','C','N','O','F','P','S','K','V','Y','I','W','U'/
      data atn1/ 1 , 5 , 6 , 7 , 8 , 9 , 15, 16, 19, 23, 39, 53, 74, 92/

      geom_tag_to_atn = .false.
c
c     eliminate conventions that refer to centers used for
c     computation purposes .. just bq and x for now
c
      match = .false.
      do i=1,nt
        match = .false.
c       TP: changed from buf = tag(nt) to buf = tag(i)
        buf = tag(i)
        lbuf = inp_strlen(buf)
        if (lbuf .eq. 0) goto 100 
c
        call inp_lcase(buf)
        if (buf(1:2) .eq. 'bq' .or. 
     $       ((buf(1:1).eq.'x') .and. (buf(2:2).ne.'e'))) then ! X but not Xe
           atn(i)   = 0
           match = .true.
           goto 100
        end if
c
c       Attempt to match the first 4 characters of the
c       full names of the elements
c
        atn(i) = 0
        if (lbuf .ge. 4) then
          do ind = 1,nelements
            if (inp_compare(.false.,buf(1:4),elements(ind)(1:4))) then
              atn(i)    = ind
              match = .true.
              goto 100
            endif
          enddo
        end if
c
c       Failed ... attempt to match the first two characters
c       against two character element names
c
        if (buf(2:2) .ne. ' ') then
           if (inp_match(nelements,.false.,buf(1:2),symbols,ind)) then
              atn(i)    = ind
              match = .true.
              goto 100
           end if
        end if
c
c       Last ditch attempt ... match against 1 character symbols
c
        if (inp_match(14, .false., buf(1:1), sym1, ind)) then
           ind = atn1(ind)
           atn(i)   = ind
           match = .true.
           goto 100
        end if
100     continue
        if(.not.match) goto 101
      end do
101   continue
      
      geom_tag_to_atn = match

      return
c
      end
c
      logical function geom_charge_center(geom)
      implicit none
#include "nwc_const.fh"
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
c
C> \brief Extracts the center of charge from a geometry instance
c
C> \return Returns .true. if the center of charge was successfully
C> established, and .false. otherwise.
c
      logical function geom_center_of_charge(geom, center)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
      integer geom               !< [Input] the geometry handle
      double precision center(3) !< [Output] the center of charge
      integer i, k
      double precision qsum
      logical geom_check_handle
      external geom_check_handle
c
c  Get the center of charge of the geometry
c  DOES not shift the center of charge to the origin
c
      geom_center_of_charge=geom_check_handle(geom,'geom_charge_center')
      if (.not. geom_center_of_charge) return
      call geom_compute_values(geom)
      qsum = 0.0d0
      do i = 1, ncenter(geom)
         qsum = qsum + charge(i,geom)
      end do
c
      if (qsum .eq. 0.0d0) then ! System is charge neutral
         do k = 1, 3
            center(k) = 0.0d0
         enddo
      else
         do k = 1, 3
            center(k) = ndipole(k,geom)/qsum
         end do
      endif
c
      geom_center_of_charge=.true.
      end
c
C> \brief Extracts the center of mass from a geometry instance
c
C> \return Returns .true. if the center of mass was successfully
C> established, and .false. otherwise.
c
      logical function geom_center_of_mass(geom, center)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "geomP.fh"
      integer geom               !< [Input] the geometry handle
      double precision center(3) !< [Output] the center of mass
      logical geom_ncent, geom_mass_get
      external geom_ncent, geom_mass_get
c
      integer i, iat, ncent
      double precision mass, amass
c
c  Get the center of mass of the geometry
c  DOES not shift the center of mass to the origin
c
      if (.not. geom_ncent(geom, ncent))
     &  call errquit('geom_mass_center: unable to get ncent',555,
     &       GEOM_ERR)
      do i = 1,3
        center(i) = 0.0d0
      enddo
      amass = 0.0d0
      do iat = 1, ncent
         if(.not.geom_mass_get(geom, iat, mass)) call
     &        errquit(' mass_get  failed ',iat, GEOM_ERR)
         amass = amass + mass
         do i=1,3
            center(i) = center(i) + mass*coords(i,iat,geom)
         enddo
      enddo
      do i = 1, 3
         center(i) = center(i)/amass
      enddo    
c
      geom_center_of_mass=.true.
      return
      end

C> \brief Extracts the centroid from a geometry instance
c
C> \return Returns .true. if the centroid was successfully
C> established, and .false. otherwise.
c
      logical function geom_centroid(geom, centroid)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "geomP.fh"
      integer geom               !< [Input] the geometry handle
      double precision centroid(3) !< [Output] the centroid
      logical geom_ncent
      external geom_ncent
c
      integer i, iat, ncent
c
      if (.not. geom_ncent(geom, ncent))
     &  call errquit('geom_centroid: unable to get ncent',555,
     &       GEOM_ERR)
c
      call dfill(3, 0.0d0, centroid, 1)
      do iat = 1, ncent
        do i=1,3
            centroid(i) = centroid(i) + coords(i,iat,geom)/ncent
        enddo
      enddo
c
      geom_centroid=.true.
      return
      end
c
C> \brief Extracts the nuclear repulsion energy from a geometry instance
c
C> \return Return .true. if successfull, and .false. otherwise.
c
      logical function geom_nuc_rep_energy(geom, energy)
#include "nwc_const.fh"
#include "geomP.fh"
#include "errquit.fh"
      integer geom            !< [Input] the geometry handle
      double precision energy !< [Output] the nuclear repulsion energy
      logical bq_add_nuc_rep_energy
      external bq_add_nuc_rep_energy
      logical geom_check_handle
      external geom_check_handle
      logical geom_extbq_on
      external geom_extbq_on
c
c     return the effective nuclear repulsion energy etc.
c
      geom_nuc_rep_energy = geom_check_handle(geom, 'geom_nuc_rep_e')
      if (.not. geom_nuc_rep_energy) return
      energy = erep(geom)
c
      if(geom_extbq_on()) then
      if(.not. bq_add_nuc_rep_energy(geom,energy))
     >   call errquit("failed bq_add_nuc_rep_energy",0,GEOM_ERR)
      end if
c
      end
c
C> \brief Extracts the total nuclear charge of a geometry instance
c
C> \return Return .true. if successfull, and .false. otherwise.
c
      logical function geom_nuc_charge(geom, total_charge)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
      integer geom                  !< [Input] the geometry handle
      double precision total_charge !< [Output] the accumulated nuclear charge
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
C> \brief Tests for super imposed atoms
c
C> A common error in geometries is to have two atoms essentially on
C> top of eachother. Usually this leads to problems in calculations.
C> This routine checks for such occurances.
c
C> \return Return .true. if no atoms are on top of eachother,
C> and .false. otherwise.
c
      logical function geom_verify_coords(geom)
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "util_params.fh"
c
c::functions      
      logical geom_ncent
      logical geom_cent_get
      external geom_ncent
      external geom_cent_get
c::passed
      integer geom !< [input] the geometry handle
c::local
      integer nat, iat, jat, num2compare, atomi, atomj
      integer i,j
      parameter (num2compare = 2)
      character*16 name(num2compare)
      double precision xyz(3,num2compare)
      double precision chg(num2compare)
      double precision dist_min
      double precision dist_thresh
      double precision dist_my
      parameter (dist_thresh = 0.5d0)
      double precision dist2
c
      dist2(i,j) =
     &    (xyz(1,i)-xyz(1,j))*(xyz(1,i)-xyz(1,j)) +
     &    (xyz(2,i)-xyz(2,j))*(xyz(2,i)-xyz(2,j)) +
     &    (xyz(3,i)-xyz(3,j))*(xyz(3,i)-xyz(3,j))
c
      if(.not.geom_ncent(geom, nat))
     &    call errquit('geom_verify_coords: geom_ncent failed',911,
     &       GEOM_ERR)
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
     &          ('geom_verify_coords: geom_cent_get<1> failed',911,
     &       GEOM_ERR)
            if(.not. geom_cent_get
     &          (geom,jat,name(atomj),xyz(1,atomj),chg(atomj)))
     &          call errquit
     &          ('geom_verify_coords: geom_cent_get<2> failed',911,
     &       GEOM_ERR)
            dist_my = dist2(atomi,atomj)
            if (dist_my.lt.dist_thresh) then
              write(luout,*)' atoms ',iat,' and ',jat,' are similar'
              write(luout,*)' atom ',iat,' coordinates',
     &            xyz(1,atomi),xyz(2,atomi),xyz(3,atomi)
              write(luout,*)' atom ',jat,' coordinates',
     &            xyz(1,atomj),xyz(2,atomj),xyz(3,atomj)
            endif
            dist_min = min(dist_min, dist_my)
          end if
00200   continue
00100 continue
      dist_min = sqrt(dist_min)
      geom_verify_coords = dist_min.gt.dist_thresh
*      write(LuOut,*)' distance minimum =',
*     &    dist_min, geom_verify_coords
c
      if (geom_verify_coords) return
      write(luout,*)'minimum distance ',dist_min
      write(luout,*)
     &    ' ************ WARNING ******************' 
      write(luout,'(A,F8.3,A)')
     &    ' at least two atoms are only ',dist_min/cang2au,
     A     ' angstrom apart '
      write(luout,*)
     &    'Please check your geometry input' 
      write(luout,*)
     &    ' ************ WARNING ******************' 
      write(luout,*)
     &    'If you like danger & want to skip this check' 
      write(luout,*)
     &    'add the following input line'
      write(luout,*)
     &    'set geom:dont_verify .true.'

      end
c
c
c---> new functions added on incorporation on symmetry and solid state codes
c
c
C> \brief Retrieve the system type for a given geometry
c
C> The system type indicates the number of dimensions in which the
C> system is periodic. 
c
C> \return Return .true. if the system type was found successfully, and
C> .false. otherwise.
c
      logical function geom_systype_get(geom, itype)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
      integer geom              !< [Input] the geometry handle
      integer itype             !< [Output] the system type
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
c
C> \brief Extracts the lattice vectors of a geometry instance
c
C> \return Return .true. if successfull, and .false. otherwise.
c
      logical function geom_latvec_get(geom,vectors)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
c     
      integer geom                !< [Input] the geometry handle
      double precision vectors(3) !< [Output] the lattice vectors
      integer i
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
c
C> \brief Extracts the lattice angles of a geometry instance
c
C> \return Return .true. if successfull, and .false. otherwise.
c
      logical function geom_latang_get(geom,angles)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
c     
      integer geom               !< [Input] the geometry handle
      double precision angles(3) !< [Output] the lattice angles
      integer i
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
c
C> \brief Extracts the reciprocal lattice vectors of a geometry instance
c
C> \return Return .true. if successfull, and .false. otherwise.
c
      logical function geom_recipvec_get(geom,rvectors)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
c     
      integer geom                 !< [Input] the geometry handle
      double precision rvectors(3) !< [Output] the reciprocal lattice vectors
      integer i
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
c
C> \brief Extracts the reciprocal lattice angles of a geometry instance
c
C> \return Return .true. if successfull, and .false. otherwise.
c
      logical function geom_recipang_get(geom,rangles)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
c
      integer geom                !< [Input] the geometry handle
      double precision rangles(3) !< [Output] the reciprocal lattice angles
      integer i
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
#include "nwc_const.fh"
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

      logical function geom_lattice_get(geom,lattice)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
#include "errquit.fh"
c     
      integer geom,i,j               ! [input]
      double precision lattice(6)    ! [output]
      double precision rad,scale
      logical  geom_check_handle,geom_get_user_scale
      external geom_check_handle,geom_get_user_scale

      geom_lattice_get = geom_check_handle(geom, 'geom_lattice_get')
      if (.not. geom_lattice_get) return
      if (.not. geom_get_user_scale(geom,scale))
     $     call errquit('geom_lattice_get: call eric!',0, GEOM_ERR)
c

      rad = 4.0d0*datan(1.0d0)/180.0d0
      lattice(1) = lattice_vectors(1,geom)*scale
      lattice(2) = lattice_vectors(2,geom)*scale
      lattice(3) = lattice_vectors(3,geom)*scale
      lattice(4) = lattice_angles(1,geom)*rad 
      lattice(5) = lattice_angles(2,geom)*rad 
      lattice(6) = lattice_angles(3,geom)*rad 
      end

      logical function geom_lattice_set(geom,lattice)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "geomP.fh"
c
      integer geom,i,j               ! [input]
      double precision lattice(6)    ! [output]
c
      integer iang
      double precision rad,dperm
      double precision c(3,3),vol,scale,amat(3,3),gmat(3,3)
      double precision c1,c2,c3,s3,cdist(3),cang(3)

*     *** external functions ***
      logical geom_check_handle, geom_get_user_scale
      external geom_check_handle,geom_get_user_scale
      double precision deter3
      external         deter3

c
      geom_lattice_set = geom_check_handle(geom, 'geom_lattice_set')
      if (.not. geom_lattice_set) return
      if (.not. geom_get_user_scale(geom,scale))
     $     call errquit('geom_lattice_set: call eric!',0, GEOM_ERR)

      rad = 180.0d0/(4.0d0*datan(1.0d0))
      lattice_vectors(1,geom) = lattice(1)/scale
      lattice_vectors(2,geom) = lattice(2)/scale
      lattice_vectors(3,geom) = lattice(3)/scale
      lattice_angles(1,geom)  = lattice(4)*rad
      lattice_angles(2,geom)  = lattice(5)*rad
      lattice_angles(3,geom)  = lattice(6)*rad
      cdist(1) = lattice(1)
      cdist(2) = lattice(2)
      cdist(3) = lattice(3)
      cang(1)  = lattice(4)
      cang(2)  = lattice(5)
      cang(3)  = lattice(6)
c
c--------> build the metrical matrix (atomic units)
c
      do 200 i=1,3
        gmat(i,i)=cdist(i)**2
  200 continue
      iang=3
      do 210 i=1,3
        do 220 j=i+1,3
          gmat(i,j)=cdist(i)*cdist(j)*dcos(cang(iang))
          gmat(j,i)=gmat(i,j)
          iang=iang-1
  220   continue
  210 continue
c
      do 230 i=1,3
        do 240 j=1,3
          metric_matrix(i,j,geom)=gmat(i,j)
  240   continue
  230 continue

      dperm = deter3(gmat)
*
      vol=dsqrt(dperm)
      volume_direct(geom)=vol
c

      c1=dcos(cang(1))
      c2=dcos(cang(2))
      c3=dcos(cang(3))
      s3=dsin(cang(3))
      amat(1,1) = cdist(1)*s3
      amat(1,2) = 0.0d+00
      amat(1,3) = (cdist(3)*(c2-c1*c3)/s3)
      amat(2,1) = cdist(1)*c3
      amat(2,2) = cdist(2)
      amat(2,3) = cdist(3)*c1
      amat(3,1) = 0.0d+00
      amat(3,2) = 0.0d+00
      amat(3,3) = (vol/(cdist(1)*cdist(2)*s3))
c
      do i=1,3
        do j=1,3
          amatrix(i,j,geom) = amat(i,j)
        end do
      end do
c
c     Mmmm ... the original code only set this stuff from the input
c     using the a,b,c,alpha,beta,gamma, but now we have changed
c     the amatrix ... need to update ainv and also recompute the
c     other crap ... for now just set the other crap to crap so that
c     we'll know if it is used
c
      do i = 1,3
         do j = 1,3
            metric_matrix(i,j,geom) = 1d300
            bmatrix(i,j,geom) = 1d300
         end do
         recip_lat_vectors(i,geom) = 1d300
         recip_lat_angles(i,geom) = 1d300
      end do
c
c     HERE SHOULD RECOMPUTE AMATRIX WITH STANDARD ORIENTATION
c     SINCE IF THE GEOMETRY IS STORED AND RELOADED THE
c     STANDARD ORIENTATION IS IMPOSED.
c
c     Update the amatrix inverse
c      - Since amat=[a1,a2,a3]
c              ainv=[b1,b2,b3]^t
c
      call dfill(9,0.0d0,c,1)
      c(1,1) = amat(2,2)*amat(3,3) - amat(3,2)*amat(2,3)  ! = b(1,1)
      c(1,2) = amat(3,2)*amat(1,3) - amat(1,2)*amat(3,3)  ! = b(2,1)
      c(1,3) = amat(1,2)*amat(2,3) - amat(2,2)*amat(1,3)  ! = b(3,1)
      c(2,1) = amat(2,3)*amat(3,1) - amat(3,3)*amat(2,1)  ! = b(1,2)
      c(2,2) = amat(3,3)*amat(1,1) - amat(1,3)*amat(3,1)  ! = b(2,2)
      c(2,3) = amat(1,3)*amat(2,1) - amat(2,3)*amat(1,1)  ! = b(3,2)
      c(3,1) = amat(2,1)*amat(3,2) - amat(3,1)*amat(2,2)  ! = b(1,3)
      c(3,2) = amat(3,1)*amat(1,2) - amat(1,1)*amat(3,2)  ! = b(2,3)
      c(3,3) = amat(1,1)*amat(2,2) - amat(2,1)*amat(1,2)  ! = b(3,3)
      vol = amat(1,1)*c(1,1)
     >    + amat(2,1)*c(1,2)

c
      call dscal(9,1.0d0/vol,c,1)
c
      call dcopy(9,c,1,amatrix_inv(1,1,geom),1)

      return
      end

c
C> \brief Retrieve the fractional to Cartesian transformation matrix
c
C> Look up the transformation matrix that transforms fractional to
C> Cartesian coordinates for a geometry instance. The matrix is a simple
C> 3x3 matrix.
c
C> \return Return .true. if the matrix was successfully found, and
C> .false. otherwise.
c
      logical function geom_amatrix_get(geom,amat)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
c     
      integer geom               !< [Input] the geometry handle
      double precision amat(3,3) !< [Output] the transformation matrix
      integer i,j
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
c
C> \brief Set the fractional to Cartesian transformation matrix
c
C> Stores the transformation matrix that transforms fractional to
C> Cartesian coordinates in a geometry instance. This matrix is a
C> simple 3x3 matrix. In addition the inverse transformation is
C> calculated and stored as well, and the lattice parameters are
C> updated.
c
C> \return Return .true. when the transformation was stored
C> successfully, and .false. otherwise
c
      logical function geom_amatrix_set(geom,amat)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "geomP.fh"
c     
      integer geom               !< [Input] the geometry handle
      double precision amat(3,3) !< [Input] the transformation matrix
      integer i,j
      logical geom_check_handle, geom_get_user_scale
      external geom_check_handle
c
      double precision c(3,3), vol, scale
c
      geom_amatrix_set = geom_check_handle(geom, 'geom_amatrix_set')
      if (.not. geom_amatrix_set) return
      if (.not. geom_get_user_scale(geom,scale)) 
     $     call errquit('geom_amtrix_set: call eric!',0, GEOM_ERR)
c
      do i=1,3
        do j=1,3
          amatrix(i,j,geom) = amat(i,j)
        end do
      end do
c
c     Mmmm ... the original code only set this stuff from the input
c     using the a,b,c,alpha,beta,gamma, but now we have changed
c     the amatrix ... need to update ainv and also recompute the
c     other crap ... for now just set the other crap to crap so that
c     we'll know if it is used
c
      do i = 1,3
         do j = 1,3
            metric_matrix(i,j,geom) = 1d300
            bmatrix(i,j,geom) = 1d300
         end do
         recip_lat_vectors(i,geom) = 1d300
         recip_lat_angles(i,geom) = 1d300
      end do
c
c     HERE SHOULD RECOMPUTE AMATRIX WITH STANDARD ORIENTATION
c     SINCE IF THE GEOMETRY IS STORED AND RELOADED THE 
c     STANDARD ORIENTATION IS IMPOSED.
c
c     Update the amatrix inverse
c      - Since amat=[a1,a2,a3]
c              ainv=[b1,b2,b3]^t
c
      call dfill(9,0.0d0,c,1)
      c(1,1) = amat(2,2)*amat(3,3) - amat(3,2)*amat(2,3)  ! = b(1,1)
      c(1,2) = amat(3,2)*amat(1,3) - amat(1,2)*amat(3,3)  ! = b(2,1)
      c(1,3) = amat(1,2)*amat(2,3) - amat(2,2)*amat(1,3)  ! = b(3,1)
      c(2,1) = amat(2,3)*amat(3,1) - amat(3,3)*amat(2,1)  ! = b(1,2)
      c(2,2) = amat(3,3)*amat(1,1) - amat(1,3)*amat(3,1)  ! = b(2,2)
      c(2,3) = amat(1,3)*amat(2,1) - amat(2,3)*amat(1,1)  ! = b(3,2)
      c(3,1) = amat(2,1)*amat(3,2) - amat(3,1)*amat(2,2)  ! = b(1,3)
      c(3,2) = amat(3,1)*amat(1,2) - amat(1,1)*amat(3,2)  ! = b(2,3)
      c(3,3) = amat(1,1)*amat(2,2) - amat(2,1)*amat(1,2)  ! = b(3,3)
      vol = amat(1,1)*c(1,1)
     >    + amat(2,1)*c(1,2)
     >    + amat(3,1)*c(1,3)
      volume_direct(geom) = vol
c
      call dscal(9,1.0d0/vol,c,1)
c
      call dcopy(9,c,1,amatrix_inv(1,1,geom),1)
c
c     Ooops ... must also update the pesky lattice parameters
c
      call xlattice_abc_abg(
     $     lattice_vectors(1,geom),
     $     lattice_vectors(2,geom),
     $     lattice_vectors(3,geom),
     $     lattice_angles(1,geom),
     $     lattice_angles(2,geom),
     $     lattice_angles(3,geom),amat)

      lattice_vectors(1,geom) = lattice_vectors(1,geom)/scale
      lattice_vectors(2,geom) = lattice_vectors(2,geom)/scale
      lattice_vectors(3,geom) = lattice_vectors(3,geom)/scale
c
      end
      subroutine xlattice_abc_abg(a,b,c,alpha,beta,gamma,lattice_unita)
      implicit none
      double precision a,b,c
      double precision alpha,beta,gamma,lattice_unita(3,3)

*     *** local variables ****
      double precision d2,pi

*     **** determine a,b,c,alpha,beta,gmma ***
      pi = 4.0d0*datan(1.0d0)
      a = dsqrt(lattice_unita(1,1)**2 
     >        + lattice_unita(2,1)**2 
     >        + lattice_unita(3,1)**2)
      b = dsqrt(lattice_unita(1,2)**2 
     >        + lattice_unita(2,2)**2 
     >        + lattice_unita(3,2)**2)
      c = dsqrt(lattice_unita(1,3)**2 
     >        + lattice_unita(2,3)**2 
     >        + lattice_unita(3,3)**2)
 
      d2 = (lattice_unita(1,2)-lattice_unita(1,3))**2 
     >   + (lattice_unita(2,2)-lattice_unita(2,3))**2 
     >   + (lattice_unita(3,2)-lattice_unita(3,3))**2
      alpha = (b*b + c*c - d2)/(2.0d0*b*c)
      alpha = dacos(alpha)*180.0d0/pi
 
      d2 = (lattice_unita(1,3)-lattice_unita(1,1))**2 
     >   + (lattice_unita(2,3)-lattice_unita(2,1))**2 
     >   + (lattice_unita(3,3)-lattice_unita(3,1))**2
      beta = (c*c + a*a - d2)/(2.0d0*c*a)
      beta = dacos(beta)*180.0d0/pi
 
      d2 = (lattice_unita(1,1)-lattice_unita(1,2))**2 
     >   + (lattice_unita(2,1)-lattice_unita(2,2))**2 
     >   + (lattice_unita(3,1)-lattice_unita(3,2))**2
      gamma = (a*a + b*b - d2)/(2.0d0*a*b)
      gamma = dacos(gamma)*180.0d0/pi

      return
      end




      logical function geom_bmatrix_get(geom,bmat)
      implicit none
#include "nwc_const.fh"
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
#include "nwc_const.fh"
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
c
C> \brief Extract the symmetry unique centers from a geometry instance
c
C> \return Return .true. if successfull, and .false. otherwise.
c
      logical function geom_uniquecent_get(geom,ncent,uniquecent)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
c     
      integer geom              !< [Input] the geometry handle
      integer ncent             !< [Input] the number of unique centers
      integer uniquecent(ncent) !< [Output] the indicies of unique centers
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
c
C> \brief Define the symmetry unique centers of a geometry instance
c
C> \return Return .true. if successfull, and .false. otherwise.
c
      logical function geom_uniquecent_set(geom,ncent,uniquecent)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
c     
      integer geom              !< [Input] the geometry handle
      integer ncent             !< [Input] the number of unique centers
      integer uniquecent(ncent) !< [Input] the indicies of unique centers
      integer i
      logical geom_check_handle
      external geom_check_handle
c
      geom_uniquecent_set=geom_check_handle(geom,'geom_uniquecent_set')
      if (.not. geom_uniquecent_set) return
c
      do i=1,ncent
         unique_cent(i,geom)=uniquecent(i)
      enddo
      end
c
C> \brief Retrieve the conversion factor from user units to atomic units
C> \return Return .true. if successfull, and .false. otherwise
      logical function geom_get_user_scale(geom, scale)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "geomP.fh"
      integer geom           !< [Input] the geometry handle
      double precision scale !< [Output] the unit conversion factor
      logical geom_check_handle
      external geom_check_handle
c
      geom_get_user_scale = 
     $     geom_check_handle(geom, 'geom_get_user_scale')
c
      if (user_units(geom) .eq. 'a.u.') then
         scale = 1.0d0
      else if (user_units(geom) .eq. 'angstroms') then
         scale = angstrom_to_au
      else if (user_units(geom) .eq. 'nanometer') then
         scale = angstrom_to_au * 10.0d0
      else if (user_units(geom) .eq. 'picometer') then
         scale = angstrom_to_au * 0.01d0
      else
         call errquit('geom_get_user_scale: unknown units',0, GEOM_ERR)
      endif
c
      end
c
C> \brief Defines the unit of the coordinates specified by the user
c
C> Internally the code always uses atomic units to store the
C> coordinates. However, in the input file the user may choose different
C> units, such as Angstrom, or nm. This function stores which units
C> the user used to specify the geometry. This allows, e.g. the 
C> geometry to be printed in the same units it was specified in.
c
C> \return Return .true. if the function was successfull, and .false.
C> otherwise.
c
      logical function geom_set_user_units(geom, units)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
      integer geom        !< [Input] the geometry handle
      character*(*) units !< [Input] the user units
      logical geom_check_handle
      external geom_check_handle
c
      geom_set_user_units =
     $     geom_check_handle(geom, 'geom_set_user_units')
      user_units(geom) = units
c
      end
c
C> \brief Retrieves the unit of the coordinates specified by the user
c
C> Internally the code always uses atomic units to store the
C> coordinates. However, in the input file the user may choose different
C> units, such as Angstrom, or nm. This function retrieves which units
C> the user used to specify the geometry. This allows, e.g. the 
C> geometry to be printed in the same units it was specified in.
c
C> \return Return .true. if the function was successfull, and .false.
C> otherwise.
c
      logical function geom_get_user_units(geom, units)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
      integer geom        !< [Input] the geometry handle
      character*(*) units !< [Output] the user units
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
#include "errquit.fh"
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
     &    911, GEOM_ERR)
      geom_tag_to_default_mass =
     &    geom_atn_to_default_mass(tag_atomic_number,mass)
      end
c
C> \brief Converts an atomic number to atomic mass
c
C> This routine returns the default atomic mass from based on the atomic
C> number.  The mass for each element comes from the book "The Elements" 
C> by John Emsley, Oxford University Press, (C) 1989, ISBN 0-19-855237-8
C> The specific mass chosen was the most abundant isotope with a known mass. 
C> When the abundance was equal the isotope with the longest half life was 
C> used.  
c
C> \return Return .true. if the conversion was successfull, and .false. 
C> otherwise.
c
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
c Updated 09/99 KG Dyall, correcting some transuranics and adding new
c values from the WebElements website www.webelements.com
c
c
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "geomP.fh"
c
      integer atn           !< [Input] the atomic number of element
      double precision mass !< [Output] the default elemental atomic mass.
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
     & 68.9257d0,  73.9219d0,   74.9216d0,  78.9183d0,  79.9165d0, 
     & 83.912d0,   84.9117d0,   87.9056d0,  88.9054d0,  89.9043d0,
     & 92.9060d0,  97.9055d0,   97.9072d0, 101.9037d0, 102.9048d0,
     &105.9032d0, 106.90509d0, 113.9036d0, 114.9041d0, 117.9018d0/
      data (def_masses(i),i=51,109) /
     & 120.9038d0, 129.9067d0, 126.9004d0, 131.9042d0, 132.9051d0,
     & 137.9050d0, 138.9061d0, 139.9053d0, 140.9074d0, 143.9099d0,
     & 144.9128d0, 151.9195d0, 152.9209d0, 157.9241d0, 159.9250d0,
     & 163.9288d0, 164.9303d0, 165.9304d0, 168.9344d0, 173.9390d0,
     & 174.9409d0, 179.9468d0, 180.948d0,  183.9510d0, 186.9560d0,
     & 189.9586d0, 192.9633d0, 194.9648d0, 196.9666d0, 201.9706d0,
     & 204.9745d0, 207.9766d0, 208.9804d0, 209.9829d0, 210.9875d0,
     & 222.0175d0, 223.0198d0, 226.0254d0, 227.0278d0, 232.0382d0,
     & 231.0359d0, 238.0508d0, 237.0482d0, 244.0642d0, 243.0614d0,
     & 247.0704d0, 247.0703d0, 251.0796d0, 252.0829d0, 257.0950d0,
     & 258.0986d0, 259.1009d0, 262.1100d0, 261.1087d0, 262.1138d0,
     & 266.1219d0, 262.1229d0, 267.1318d0, 268.1388d0 /
c
      geom_atn_to_default_mass = .false.
c
      if (atn.lt.0) call errquit
     &    ('geom_atn_to_default_mass: negative atomic number',atn,
     &       GEOM_ERR)
      if (atn.gt.nelements) call errquit
     &    ('geom_atn_to_default_mass: atomic number too large',atn,
     &       GEOM_ERR)
c
      if (atn.eq.0) then
        mass = 0.0d00  ! Bq centers have no mass
      else
        mass = def_masses(atn)
      endif
      geom_atn_to_default_mass = .true.
c
      end
C>
C> \brief Define the atomic masses of the centers in a geometry instance
C>
C> \return Return .true. if successfull, and .false. otherwise
C>
      logical function geom_masses_set(geom, ncent, masses)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
#include "stdio.fh"
c
      integer geom                   !< [Input] the geometry handle
      integer ncent                  !< [Input] the number of centers
      double precision masses(ncent) !< [Input] the mass on each center
c
      integer i
c
      logical geom_check_handle
      external geom_check_handle
c
      geom_masses_set = geom_check_handle(geom, 'geom_masses_set')
      if (.not. geom_masses_set) return
c
      if (ncent.le.0) then
         write(LuOut,*) ' geom_masses_set: too few centers ',ncent,
     $        names(geom)(1:lenn(geom))
         geom_masses_set = .false.
         return
      else if (ncent.gt.max_cent) then
         write(LuOut,*) ' geom_masses_set: too many centers ',ncent,
     $        names(geom)(1:lenn(geom))
         geom_masses_set = .false.
         return
      end if
c
      do i = 1, ncent
        geom_mass(i,geom) = masses(i)
      enddo
c
      end      
c
C> \brief Retrieve the masses of the centers
c
C> Retrieves the masses associated with the centers in a geometry
C> instance.
c
C> \return Return .true. if the function was successfull, and .false.
C> otherwise.
c
      logical function geom_masses_get(geom, ncent, masses)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
c
      integer geom                   ! [Input] the geometry handle
      integer ncent                  ! [Output] the number of centers
      double precision masses(ncent) ! [Output] the mass on each center
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
c
C> \brief Define the atomic mass of the specified center in a geometry instance
c
C> \return Return .true. if successfull, and .false. otherwise
c
      logical function geom_mass_set(geom, icent, mass)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
#include "stdio.fh"
c
      integer geom            !< [Input] the geometry handle
      integer icent           !< [Input] the center rank
      double precision mass   !< [Input] the mass on center icent
c
      logical geom_check_handle
      external geom_check_handle
c
      geom_mass_set = geom_check_handle(geom, 'geom_mass_set')
      if (.not. geom_mass_set) return
c
      if (icent.le.0 .or. icent.gt.ncenter(geom)) then
         write(LuOut,*) ' geom_mass_set: icent out of range',icent,
     &        ncenter(geom),
     $        names(geom)(1:lenn(geom))
         return
      end if
c
      geom_mass(icent,geom) = mass
c
      end      
c
C> \brief Retrieve the atomic mass of the specified center in a geometry instance
c
C> \return Return .true. if successfull, and .false. otherwise
c
      logical function geom_mass_get(geom, icent, mass)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
#include "stdio.fh"
c
      integer geom            !< [Input] the geometry handle
      integer icent           !< [Input] the center rank
      double precision mass   !< [Output] the mass on center icent
c
      logical geom_check_handle
      external geom_check_handle
c
      geom_mass_get = geom_check_handle(geom, 'geom_mass_get')
      if (.not. geom_mass_get) return
c
      if (icent.le.0 .or. icent.gt.ncenter(geom)) then
         write(LuOut,*) ' geom_mass_get: icent out of range',icent,
     &        ncenter(geom),
     $        names(geom)(1:lenn(geom))
         return
      end if
c
      mass = geom_mass(icent,geom)
c
      end      
C>
C> \brief Define the atom constraint type of the centers in a geometry
C> instance
C>
C> \return Return .true. if successfull, and .false. otherwise
C>
      logical function geom_atomct_set(geom, ncent, atomct)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
#include "stdio.fh"
c
      integer geom                   !< [Input] the geometry handle
      integer ncent                  !< [Input] the number of centers
      double precision atomct(ncent) !< [Input] the atom constraint type
                                     !< on each center
c
      integer i
c
      logical geom_check_handle
      external geom_check_handle
c
      geom_atomct_set = geom_check_handle(geom, 'geom_atomct_set')
      if (.not. geom_atomct_set)  return
c
      if (ncent.le.0) then
         write(LuOut,*) ' geom_atomct_set: too few centers ',ncent,
     $        names(geom)(1:lenn(geom))
         geom_atomct_set = .false.
         return
      else if (ncent.gt.max_cent) then
         write(LuOut,*) ' geom_atomct_set: too many centers ',ncent,
     $        names(geom)(1:lenn(geom))
         geom_atomct_set = .false.
         return
      end if
c
      do i = 1, ncent
        geom_atomct(i,geom) = atomct(i)
      enddo
c
      end      
C>
C> \brief Retrieve the atom constraint type of the centers in a geometry
C> instance
C>
C> \return Return .true. if successfull, and .false. otherwise
C>
      logical function geom_atomct_get(geom, ncent, atomct)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
#include "stdio.fh"
c
      integer geom                   !< [Input] the geometry handle
      integer ncent                  !< [Output] the number of centers
      double precision atomct(ncent) !< [Output] the atom constraint
                                     !< type on each center
c
      integer i
c
      logical geom_check_handle
      external geom_check_handle
c
      geom_atomct_get = geom_check_handle(geom, 'geom_atomct_set')
      if (.not. geom_atomct_get) return
c
      ncent = ncenter(geom)
c
      do i = 1, ncent
        atomct(i) = geom_atomct(i,geom)
      enddo
c
      end      
c
C> \brief Set the Angstrom to Bohr conversion factor for a geometry
C> instance
c
C> \return Return .true. if successfull, and .false. otherwise
c
      logical function geom_set_ang2au(geom,value)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
c::functions
      logical geom_check_handle
      external geom_check_handle
c::passed
      integer geom               !< [Input] the geometry handle
      double precision value     !< [Input] the conversion factor from
c                                !< angstroms to au value ~1.8...
c
      geom_set_ang2au = geom_check_handle(geom,'geom_set_ang2au')
      if (.not. geom_set_ang2au) return
c
      angstrom_to_au = value
c
      end
c
C> \brief Retrieve the Angstrom to Bohr conversion factor for a geometry
C> instance
c
C> \return Return .true. if successfull, and .false. otherwise
c
      logical function geom_get_ang2au(geom,value)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
c::functions
      logical geom_check_handle
      external geom_check_handle
c::passed
      integer geom           !< [Input] the geometry handle
      double precision value !< [Output] the conversion factor from
c                            !< angstroms to au value ~1.8......
c
      geom_get_ang2au = geom_check_handle(geom,'geom_get_ang2au')
      if (.not. geom_get_ang2au) return
c
      value = angstrom_to_au
c
      end
      logical function geom_set_au2ang(geom,value)
      implicit none
#include "nwc_const.fh"
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
#include "nwc_const.fh"
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
c
C> \brief Define the centers in a geometry instance that have an ECP
c
C> \return Return .true. if successfull, and .false. otherwise
c
      logical function geom_ecp_allset(geom,ncenter_in,oecp)
      implicit none
#include "errquit.fh"
c
#include "nwc_const.fh"
#include "geomP.fh"
c
      integer geom             ! [Input] the geometry handle
      integer ncenter_in       ! [Input] the number of centers
      logical oecp(ncenter_in) ! [Input] array of T/F for having ECPs
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
     &    (ncenter(geom)-ncenter_in), GEOM_ERR)
c
      do icenter = 1,ncenter_in
        oecpcent(icenter,geom) = oecp(icenter)
      enddo
c
      end
c
C> \brief Retrieve the centers in a geometry instance that have an ECP
c
C> \return Return .true. if successfull, and .false. otherwise
c
      logical function geom_ecp_allget(geom,ncenter_in,oecp)
      implicit none
#include "errquit.fh"
c
#include "nwc_const.fh"
#include "geomP.fh"
c
      integer geom             !< [Input] the geometry handle
      integer ncenter_in       !< [Input] the number of centers
      logical oecp(ncenter_in) !< [Output] array of T/F for having ECPs
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
     &    (ncenter(geom)-ncenter_in), GEOM_ERR)
c
      do icenter = 1,ncenter_in
        oecp(icenter)= oecpcent(icenter,geom) 
      enddo
c
      end

      logical function geom_ecp_set(geom,icent,oecp)
      implicit none
#include "errquit.fh"
c
#include "stdio.fh"
#include "nwc_const.fh"
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
      if (.not.(icent.gt.0.and.icent.le.ncenter(geom))) then
        write(luout,*)' icent   = ',icent
        write(luout,*)' ncenter = ',ncenter(geom)
        call errquit('geom_ecp_set: icent out of range ncenter = ',911,
     &       GEOM_ERR)
      endif
c
      oecpcent(icent,geom) = oecp
c
      end
      logical function geom_ecp_get(geom,icent)
      implicit none
#include "errquit.fh"
c
#include "stdio.fh"
#include "nwc_const.fh"
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
      if (.not.(icent.gt.0.and.icent.le.ncenter(geom))) then
        write(luout,*)' icent   = ',icent
        write(luout,*)' ncenter = ',ncenter(geom)
        call errquit('geom_ecp_get: icent out of range ncenter = ',911,
     &       GEOM_ERR)
      endif
c
      geom_ecp_get = oecpcent(icent,geom) 
c
      end
      logical function geom_ncent_ecp(geom, ncent_ecp)
      implicit none
#include "nwc_const.fh"
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
#include "errquit.fh"
#include "nwc_const.fh"
#include "geomP.fh"
#include "stdio.fh"
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
     &         ' than the coord array dimension which is:',ncent_in,
     &       GEOM_ERR)
*          write(LuOut,*)' geom       = ',geom
*          write(LuOut,*)' ncent_ecp  = ',ncent_ecp
*          write(LuOut,*)' icent      = ',icent
*          write(LuOut,*)' coords geom 1',coords(1,icent,geom)
*          write(LuOut,*)' coords geom 2',coords(2,icent,geom)
*          write(LuOut,*)' coords geom 3',coords(3,icent,geom)
          
          coords_ecp(1,ncent_ecp) = coords(1,icent,geom)
          coords_ecp(2,ncent_ecp) = coords(2,icent,geom)
          coords_ecp(3,ncent_ecp) = coords(3,icent,geom)

*          write(LuOut,*)' coords ecp 1',coords_ecp(1,ncent_ecp)
*          write(LuOut,*)' coords ecp 2',coords_ecp(2,ncent_ecp)
*          write(LuOut,*)' coords ecp 3',coords_ecp(3,ncent_ecp)
        endif
      enddo
*      write(LuOut,*)' coordinates inside geom_coords_ecp'
*      call output(coords_ecp,1,3,1,ncent_ecp,3,ncent_ecp,1)
c      
      end
      logical function geom_any_ecp(geom)
      implicit none
#include "nwc_const.fh"
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
#include "errquit.fh"
#include "nwc_const.fh"
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
     &         ' than array size passed in which is:',num_ecp_cent,
     &       GEOM_ERR)
          ecp_cent(num_ecp) = icent
        endif
      enddo
      end
      logical function geom_nuc_dipole(geom,dip)
      implicit none
#include "nwc_const.fh"
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
#include "errquit.fh"
c
c computes the angle (in degrees) between 3 atoms in order given
c
#include "stdio.fh"
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
      double precision thresh
      parameter (thresh=1.0d-6)
c::-statement function
      logical is_it_close_to  
      double precision value,test
*---          is value close to test?
      is_it_close_to(value,test) = (abs(value-test).lt.thresh)
c
      pi = 2.0d00*acos(0.0d00)
      geom_calc_angle = geom_calc_distance(a,b,ab)
      geom_calc_angle = geom_calc_angle.and.geom_calc_distance(b,c,bc)
      geom_calc_angle = geom_calc_angle.and.geom_calc_distance(a,c,ac)
      if (.not.geom_calc_angle) call errquit
     &    ('geom_calc_angle:error computing a distance',911, GEOM_ERR)
      
      xcosine = ab*ab + bc*bc - ac*ac
      if (is_it_close_to(ab,0.0d00).or.
     &    is_it_close_to(bc,0.0d00)) then
        write(luout,*)' fatal error in geom_calc_angle '
        write(luout,*)' distance ab ',ab
        write(luout,*)' distance ac ',ac
        write(luout,*)' distance bc ',bc
        write(luout,*)' please report this data to:'
        write(luout,*)'         nwchem-users@emsl.pnl.gov'
        geom_calc_angle = .false.
        angle = -565.6589d00
        return
      endif
      xcosine = xcosine/(2.0d00*ab*bc)
   
      if( abs(xcosine) .gt. 1.00d00 ) xcosine = sign(1.0d00,xcosine)

      angle = (180.0d00/pi)*acos(xcosine)
      
      end
      logical function geom_calc_dihedral(ain,bin,cin,din,dihedral)
      implicit none 
#include "errquit.fh"
c
c computes the dihedral angle for the given 4 atom coordinates
c
c::-includes
#include "stdio.fh"
c::-functions
      logical geom_calc_angle
      external geom_calc_angle
c::-passed      
      double precision ain(3) ! [input] coordinates of center a
      double precision bin(3) ! [input] coordinates of center b
      double precision cin(3) ! [input] coordinates of center c
      double precision din(3) ! [input] coordinates of center d
      double precision dihedral ! [output] the dihedral angle (in degrees)
c::-local
      double precision abc, bcd, abd, acd
      double precision a(3),b(3),c(3),d(3)
      double precision pi
      double precision BA(3), BC(3), CB(3), CD(3)
      double precision BAxBC(3), CBxCD(3) 
      double precision mbaxbc, mcbxcd
      double precision cosangle
      double precision threshcos
      parameter (threshcos=1.0d-6)
      double precision thresh
      parameter (thresh = 1.0d-3)
      logical linear1
      logical linear2
c
c::-statement function
      logical is_it_close_to  
      double precision value,test
*---          is value close to test?
      is_it_close_to(value,test) = (abs(value-test).lt.thresh)
c
      pi = 2.0d00*acos(0.0d00)
      geom_calc_dihedral = .true.
      dihedral = -565.6589d00
* compute appropriate angles 
      geom_calc_dihedral = geom_calc_angle(ain,bin,cin,abc)
      geom_calc_dihedral = geom_calc_dihedral.and.
     &    geom_calc_angle(bin,cin,din,bcd) 
      geom_calc_dihedral = geom_calc_dihedral.and.
     &    geom_calc_angle(ain,bin,din,abd)
      geom_calc_dihedral = geom_calc_dihedral.and.
     &    geom_calc_angle(ain,cin,din,acd)
      if (.not.geom_calc_dihedral) then
        write(luout,*)' angle   abc ',abc
        write(luout,*)' angle   bcd ',bcd
        write(luout,*)' angle   abd ',abd
        write(luout,*)' angle   acd ',acd
        write(luout,*)' please report this data to:'
        write(luout,*)'         nwchem-users@emsl.pnl.gov'
        call util_flush(luout)
        call errquit
     &    ('geom_calc_dihedral: fatal angle error',1, GEOM_ERR)
      endif
*
* check special cases  a,b,c or b,c,d are linear
      linear1 =            is_it_close_to(abc,0.0d00)
      linear1 = linear1.or.is_it_close_to(abc,180.0d00)
      linear1 = linear1.or.is_it_close_to(bcd,0.0d00)
      linear1 = linear1.or.is_it_close_to(bcd,180.0d00)
      if (linear1) then
        dihedral = 0.0d00 
        return
      endif
* a,b,d or a,c,d are linear
      linear2 =            is_it_close_to(abd,0.0d00)
      linear2 = linear2.or.is_it_close_to(acd,0.0d00)
      if (linear2) then
        dihedral = 180.0d00 
        return
      endif
c
*... abc (b center)
      call dcopy(3,ain,1,a,1)
      call dcopy(3,bin,1,b,1)
      call dcopy(3,cin,1,c,1)
* form vectors BA and BC (make B the origin)
      BA(1) = a(1)-b(1)
      BA(2) = a(2)-b(2)
      BA(3) = a(3)-b(3)
      BC(1) = c(1)-b(1)
      BC(2) = c(2)-b(2)
      BC(3) = c(3)-b(3)
* form cross product of BA and BC
      BAxBC(1) = BA(2)*BC(3)-BA(3)*BC(2)
      BAxBC(2) = BA(3)*BC(1)-BA(1)*BC(3)
      BAxBC(3) = BA(1)*BC(2)-BA(2)*BC(1)
* find magnitude of BAxBC
      mbaxbc = BAxBC(1)*BAxBC(1) + BAxBC(2)*BAxBC(2) + BAxBC(3)*BAxBC(3)
      mbaxbc = sqrt(mbaxbc)
c
*... bcd (c center)  ! right hand screw!!
      call dcopy(3,bin,1,b,1)
      call dcopy(3,cin,1,c,1)
      call dcopy(3,din,1,d,1)
* form vectors CB and CD (make C the origin)
      CB(1) = b(1) - c(1)
      CB(2) = b(2) - c(2)
      CB(3) = b(3) - c(3)
      CD(1) = d(1) - c(1)
      CD(2) = d(2) - c(2)
      CD(3) = d(3) - c(3)
* form cross product of CB and CD
      CBxCD(1) = CB(2)*CD(3)-CB(3)*CD(2)
      CBxCD(2) = CB(3)*CD(1)-CB(1)*CD(3)
      CBxCD(3) = CB(1)*CD(2)-CB(2)*CD(1)
* now find the angle between two vectors BAxBC and CBxCD      
* find magnitude of CBxCD      
      mcbxcd = CBxCD(1)*CBxCD(1) + CBxCD(2)*CBxCD(2) + CBxCD(3)*CBxCD(3)
      mcbxcd = sqrt(mcbxcd)
*
      cosangle = BAxBC(1)*CBxCD(1) + BAxBC(2)*CBxCD(2) +
     &    BAxBC(3)*CBxCD(3) 
      if (is_it_close_to(mbaxbc,0.0d00).or.
     &    is_it_close_to(mcbxcd,0.0d00)) then
        write(luout,*)' fatal error in geom_calc_dihedral '
        write(luout,*)' mbaxbc      ',mbaxbc
        write(luout,*)' mcbxcd      ',mcbxcd
        write(luout,*)'a coordinates',ain
        write(luout,*)'b coordinates',bin
        write(luout,*)'c coordinates',cin
        write(luout,*)'d coordinates',din
        write(luout,*)' angle   abc ',abc
        write(luout,*)' angle   bcd ',bcd
        write(luout,*)' angle   abd ',abd
        write(luout,*)' angle   acd ',acd
        write(luout,*)' please report this data to:'
        write(luout,*)'         nwchem-users@emsl.pnl.gov'
        call util_flush(luout)
        geom_calc_dihedral = .false.
        return
      endif
      cosangle = cosangle/mbaxbc/mcbxcd
      if (cosangle.gt.1.0d00) then
        abc = cosangle - 1.0d00
        if (abs(abc).lt.threshcos) cosangle = cosangle - abc
      endif
      if (cosangle.lt.-1.0d00) then
        abc = -1.0d00 - cosangle
        if (abs(abc).lt.threshcos) cosangle = cosangle + abc
      endif
      dihedral = acos(cosangle)
      dihedral = dihedral*180.0d00/pi
      end
      logical function geom_print_distances(geom)
      implicit none
#include "errquit.fh"
c
c prints arbitrary i>j atom distances
c
#include "stdio.fh"
#include "inp.fh"
c::-functions
      logical geom_get_user_units
      logical geom_get_user_scale
      logical geom_ncent
      logical geom_cent_get
      logical geom_tag_to_element
      logical geom_calc_distance
      logical geom_get_def_rcov
      external geom_get_user_units
      external geom_get_user_scale
      external geom_ncent
      external geom_cent_get
      external geom_tag_to_element
      external geom_calc_distance
      external geom_get_def_rcov
c::-passed      
      integer geom ! [input] geometry handle
c::-local
      integer nat ! number of atoms
      integer iat ! ith atom
      integer jat ! jth atom
      double precision chg   ! charge (ignored)
      double precision ci(3) ! coords of atom i
      character*16 tagi      ! tag of atom i
      double precision cj(3) ! coords of atom j
      character*16 tagj      ! tag of atom j
      logical status_tagi, status_tagj ! return status of call to geom-2-element
      integer iatn, jatn ! atomic numbers for atom i and j
      character*2 symi, symj ! atomic symbols for atom i and j
      character*16 elei, elej ! atomic names for atom i and j
      double precision i_rcov, j_rcov ! covalent radii for atom i and j
      double precision rcov ! combined covalent radii
      double precision rscale ! scale factor
      integer lmtag
      double precision dij   ! distance between atoms i and j
      character*10 usr_units ! units user used as input
      double precision usr_scale ! unit scale factor
      character*128 emsg
      integer num_prt
      logical header
      logical debug
      integer ludbg
c
c
      geom_print_distances = .false.
      ludbg = 69
      debug = .false.
      header = .false.
      num_prt = 0
      rscale = 1.1d00
      if (.not.geom_get_user_units(geom,usr_units)) call errquit
     &    ('geom_print_distances: geom_get_user_units failed',911,
     &       GEOM_ERR)
      if (.not.geom_get_user_scale(geom,usr_scale)) call errquit
     &    ('geom_print_distances: geom_get_user_scale failed',911,
     &       GEOM_ERR)
      if (.not.geom_ncent(geom,nat)) call errquit
     &    ('geom_print_distances: geom_ncent failed',911,
     &       GEOM_ERR)
      if (nat.eq.1) then
        geom_print_distances = .true.
        return
      endif
      do iat = 1,nat
        if (.not.geom_cent_get(geom,iat,tagi,ci,chg)) call errquit
     &      ('geom_print_distances: geom_cent_get failed:i',911,
     &       GEOM_ERR)
        status_tagi = geom_tag_to_element(tagi,symi,elei,iatn)
        if ((symi.eq.'bq').and.
     &      (.not.status_tagi))status_tagi = .true.
        if (.not.status_tagi)call errquit
     &      ('geom_print_distances:geom_tag_to_element failed:i',911,
     &       GEOM_ERR)
        if (.not.geom_get_def_rcov(iatn,i_rcov)) call errquit
     &      ('geom_print_distances: geom_get_def_rcov failed atom i', 
     &      911, GEOM_ERR)
        lmtag = inp_strlen(tagi)
        do jat = 1,iat
          if (iat.ne.jat) then

            if (.not.geom_cent_get(geom,jat,tagj,cj,chg)) call errquit
     &          ('geom_print_distances: geom_cent_get failed:j',911,
     &       GEOM_ERR)
            status_tagj = geom_tag_to_element(tagj,symj,elej,jatn)
            if ((symj.eq.'bq').and.
     &          (.not.status_tagj))status_tagj = .true.
            if (.not.status_tagj) call errquit
     &      ('geom_print_distances:geom_tag_to_element failed:j',911,
     &       GEOM_ERR)
            if (.not.geom_get_def_rcov(jatn,j_rcov)) then
              emsg = 'geom_print_distances: '//
     &            'geom_get_def_rcov failed atom j'
              call errquit(emsg,911, GEOM_ERR)
            endif
            if (.not.geom_calc_distance(ci,cj,dij)) call errquit
     &          ('geom_print_distances: ',911, GEOM_ERR)

            rcov = rscale*(j_rcov+i_rcov)
            if (debug) then
              write(ludbg,*)'**************** iat,jat',iat,jat
              write(ludbg,*)' rcov ',rcov
              write(ludbg,*)' rscale ',rscale
              write(ludbg,*)' i_rcov ',i_rcov
              write(ludbg,*)' j_rcov ',j_rcov
              write(ludbg,10002)
     &            tagi(1:lmtag),symi,iat,
     &            tagj(1:lmtag),symj,jat,dij
            endif
            if ((dij.lt.rcov).or.debug) then
              lmtag = max(lmtag,inp_strlen(tagj))
              if (.not.header) then
                write(luout,10000)usr_units(1:inp_strlen(usr_units))
                header = .true.
              endif
              num_prt = num_prt + 1
              write(luout,10001)
     &            iat,tagi,
     &            jat,tagj,
     &            dij,(dij/usr_scale)
            endif
          endif
        enddo
      enddo
      if (header) then
        write(luout,10003)
        write(luout,10004) num_prt
        write(luout,10005)
        write(luout,10006)
      endif
10000 format(1x,78('='),/,
     &    32x,'internuclear distances',/,1x,78('-'),/,
     &    7x,'center one',6x,'|',
     &    6x,'center two',6x,'|',
     &    ' atomic units |',1x,a10,
     &    /,1x,78('-'))
10001 format(1x,
     &    i4,1x,a16,1x,'|',
     &    i4,1x,a16,1x,'|',
     &    1x,f11.5,2x,'|',1x,f11.5)
10002 format(1x,'debug:distance(',
     &    a,'|',a2,'|',i4,',',
     &    a,'|',a2,'|',i4,') =',f12.6)
10003 format(1x,78('-'))
10004 format(25x,'number of included internuclear distances: ',i10)
10005 format(1x,78('='))
10006 format(/,/)
      geom_print_distances = .true.
      end
      logical function geom_print_angles(geom)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      logical geom_prt_angles
      logical geom_ncent
      external geom_prt_angles
      external geom_ncent
      integer geom
      integer nat
*     integer max_netp
*     parameter (max_netp=24)
      integer max_net
      integer h_xnet, k_xnet, h_xlist, k_xlist
*
      if (.not.geom_ncent(geom,nat)) call errquit
     &    ('geom_print_angles: geom_ncent',911, GEOM_ERR)

*24 seems to break      max_net = min(max_netp,nat)
      max_net = nat
      if (.not.ma_push_get(mt_int,(max_net*nat),'p_xnet',
     &    h_xnet,k_xnet)) call errquit(
     &    'geom_print_angles: ma get xnet failed',911, MA_ERR)

      if (.not.ma_push_get(mt_int,(nat),'p_xlist',
     &    h_xlist,k_xlist)) call errquit(
     &    'geom_print_angles: ma get xlist failed',911, MA_ERR)

      geom_print_angles =
     &    geom_prt_angles(geom,nat,max_net,
     &    int_mb(k_xnet),int_mb(k_xlist))
      geom_print_angles = geom_print_angles .and.
     &    ma_pop_stack(h_xlist)
      geom_print_angles = geom_print_angles .and.
     &    ma_pop_stack(h_xnet)
      end
      logical function geom_prt_angles(geom,nat,max_net,xnet,xlist)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "mafdecls.fh"
c::-functions
      logical geom_cent_get
      logical geom_tag_to_element
      logical geom_calc_distance
      logical geom_calc_angle
      logical geom_get_def_rcov
      external geom_cent_get
      external geom_tag_to_element
      external geom_calc_distance
      external geom_calc_angle
      external geom_get_def_rcov
c::-passed      
      integer geom ! [input] geometry handle
      integer nat ! number of atoms
      integer max_net ! maximum number of "connected" atoms for a given atom
      integer xlist(nat)
      integer xnet(max_net,nat)
c::-local
      double precision rscale
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
*. . . . . . . . . . . . . . ! return status of call to geom-2-element
      logical status_tagi, status_tagj, status_tagk 
      integer iatn, jatn, katn ! atomic numbers for atom i, j and k
      character*2 symi, symj, symk ! atomic symbols for atom i, j and k
      character*16 elei, elej, elek ! atomic names for atom i, j and k
*. . . . . . . . . . . . . . . . . ! covalent radii for atom i, j and k
      character*128 emsg
      double precision i_rcov, j_rcov, k_rcov 
      integer num_prt
      integer itmp, jtmp, ktmp
      logical header
      integer ludbg
      logical debug
c
c initialize variables
      ludbg = 69
      debug = .false.
      header = .false.
      rscale = 1.1d00
      FF = .false.
      FT = .true.
      dij_okay = FF
      djk_okay = FF
      dik_okay = FF
      num_prt = 0
      
      geom_prt_angles = FF
      if (nat.lt.3) then
        geom_prt_angles = FT
        return
      endif
      call ifill((max_net*nat),0,xnet,1)
      call ifill(nat,0,xlist,1)
      do iat = 1,nat
        if (.not.geom_cent_get(geom,iat,tagi,ci,chg)) call errquit
     &      ('geom_prt_angles: geom_cent_get:i',911, GEOM_ERR)
        status_tagi = geom_tag_to_element(tagi,symi,elei,iatn)
        if ((symi.eq.'bq').and.
     &      (.not.status_tagi))status_tagi = .true.
        if (.not.status_tagi) call errquit
     &      ('geom_prt_angles:geom_tag_to_element failed:i',911,
     &       GEOM_ERR)
        if (.not.geom_get_def_rcov(iatn,i_rcov)) call errquit
     &      ('geom_prt_angles: geom_get_def_rcov failed atom i',911,
     &       GEOM_ERR)
        do jat = 1,nat
          
          if (iat.ne.jat) then
            
            if (.not.geom_cent_get(geom,jat,tagj,cj,chg)) call errquit
     &          ('geom_prt_angles:geom_cent_get:j ',911, GEOM_ERR)
            
            status_tagj = geom_tag_to_element(tagj,symj,elej,jatn)
            if ((symj.eq.'bq').and.
     &          (.not.status_tagj))status_tagj = .true.
            if (.not.status_tagj) call errquit
     &          ('geom_prt_angles:geom_tag_to_element failed:j',911,
     &       GEOM_ERR)
            if (.not.geom_get_def_rcov(jatn,j_rcov)) call errquit
     &          ('geom_prt_angles: geom_get_def_rcov failed atom j',
     &          911, GEOM_ERR)
            if (.not.geom_calc_distance(ci,cj,dij)) call errquit
     &          ('geom_prt_angles:geom_calc_distance:ij ',911, GEOM_ERR)
            
            if (dij.lt.(rscale*(i_rcov+j_rcov))) then
              itmp = xlist(iat) + 1
              if(itmp.gt.max_net) call errquit(
     &            'geom_prt_angles:max_net is too small ',max_net,
     &       GEOM_ERR)
              xlist(iat) = itmp
              xnet(itmp,iat) = jat
            endif
          endif
        enddo
      enddo
*rak:      write(LuOut,*)' xlist: ', xlist
*rak:      do iat = 1,nat
*rak:        write(LuOut,*)' xnet: ',iat,':',(xnet(jat,iat),jat=1,max_net)
*rak:      enddo
*
      lmtag = 0
*
      do iat = 1,nat
        if (.not.geom_cent_get(geom,iat,tagi,ci,chg)) call errquit
     &      ('geom_prt_angles: geom_cent_get:i',911, GEOM_ERR)
        status_tagi = geom_tag_to_element(tagi,symi,elei,iatn)
        if ((symi.eq.'bq').and.
     &      (.not.status_tagi))status_tagi = .true.
        if (.not.status_tagi) call errquit
     &      ('geom_prt_angles:geom_tag_to_element failed:i',911,
     &       GEOM_ERR)
        if (.not.geom_get_def_rcov(iatn,i_rcov)) call errquit
     &      ('geom_prt_angles: geom_get_def_rcov failed atom i',911,
     &       GEOM_ERR)
        if (xlist(iat).gt.1) then
          do jtmp = 1,xlist(iat)
            jat = xnet(jtmp,iat)
            if (iat.ne.jat) then
              
              if (.not.geom_cent_get(geom,jat,tagj,cj,chg))
     &            call errquit
     &            ('geom_prt_angles:geom_cent_get:j ',911, GEOM_ERR)
              
              status_tagj = geom_tag_to_element(tagj,symj,elej,jatn)
              if ((symj.eq.'bq').and.
     &            (.not.status_tagj))status_tagj = .true.
              if (.not.status_tagj) call errquit
     &            ('geom_prt_angles:geom_tag_to_element failed:j',
     &            911, GEOM_ERR)
              if (.not.geom_get_def_rcov(jatn,j_rcov)) call errquit
     &            ('geom_prt_angles:geom_get_def_rcov fail atom j',
     &            911, GEOM_ERR)
              if (.not.geom_calc_distance(ci,cj,dij)) call errquit
     &            ('geom_prt_angles:geom_calc_distance:ij ',911,
     &       GEOM_ERR)
              
              dij_okay = dij.lt.(rscale*(i_rcov+j_rcov))
              if (dij_okay.or.debug) then
                do ktmp = jtmp+1,xlist(iat)
                  kat = xnet(ktmp,iat)
                  if (kat.ne.jat.and.kat.ne.iat) then
                    if (.not.geom_cent_get(geom,kat,tagk,ck,chg))
     &                  call errquit
     &                  ('geom_prt_angles:geom_cent_get:k ',911,
     &       GEOM_ERR)
                    status_tagk =
     &                  geom_tag_to_element(tagk,symk,elek,katn)
                    if ((symk.eq.'bq').and.
     &                  (.not.status_tagk))status_tagk = .true.
                    if (.not.status_tagk) then 
                      emsg = 'geom_prt_angles: '//
     &                    'geom_tag_to_element failed:k'
                      call errquit(emsg,911, GEOM_ERR)
                    endif
                    if (.not.geom_get_def_rcov(katn,k_rcov)) then
                      emsg = 'geom_prt_angles: '//
     &                    'geom_egt_def_rcov failed atom k'
                      call errquit(emsg,911, GEOM_ERR)
                    endif
                    lmtag = max(lmtag,inp_strlen(tagk))
                    
                    if (.not.geom_calc_distance(ci,ck,dik))
     &                  call errquit
     &                  ('geom_prt_angles:geom_calc_distance:ik ',
     &                  911, GEOM_ERR)
                    if (.not.geom_calc_distance(cj,ck,djk))
     &                  call errquit
     &                  ('geom_prt_angles:geom_calc_distance:jk ',
     &                  911, GEOM_ERR)
                    dik_okay = dik.lt.(rscale*(i_rcov+k_rcov))
                    djk_okay = djk.lt.(rscale*(j_rcov+k_rcov))
                    ngood = 0
                    if (dij_okay) ngood = ngood + 1
                    if (dik_okay) ngood = ngood + 1
                    if (djk_okay) ngood = ngood + 1
                    if (debug) then
                      write(ludbg,*)'**************** iat,jat,kat',
     &                    iat,jat,kat
                      write(ludbg,*)' ngood   : ',ngood
                      write(ludbg,*)' dij_okay: ',dij_okay
                      write(ludbg,*)' dik_okay: ',dik_okay
                      write(ludbg,*)' djk_okay: ',djk_okay
                      write(ludbg,*)' dij     : ',dij
                      write(ludbg,*)' dik     : ',dik
                      write(ludbg,*)' djk     : ',djk
                      write(ludbg,*)' rij     : ',
     &                    rscale*(i_rcov+j_rcov)
                      write(ludbg,*)' rik     : ',
     &                    rscale*(i_rcov+k_rcov)
                      write(ludbg,*)' rjk     : ',
     &                    rscale*(j_rcov+k_rcov)
                    endif
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
                        emsg = 'geom_prt_angles: '//
     &                      'should not get here 1'
                        call errquit(emsg,911, GEOM_ERR)
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
                        emsg = 'geom_prt_angles: '//
     &                      'should not get here 2'
                        call errquit(emsg,911, GEOM_ERR)
                      endif
                    endif
                    should_print = (ngood.eq.2.or.ngood.eq.3) .and.
     &                  (print_ijk.or.print_ikj.or.print_jik)
                    if (should_print.and.(.not.header)) then
                      write(luout,10000)
                      header = .true.
                    endif
                    if (print_ijk) then
                      if (.not.should_print) call errquit(
     &                    'geom_prt_angles "should_print" error',
     &                    911, GEOM_ERR)
                      if (.not.geom_calc_angle(ci,cj,ck,angle))
     &                    call errquit
     &                    ('geom_prt_angles:geom_calc_angle failed',
     &                    911, GEOM_ERR)
                      num_prt =num_prt + 1
                      write(luout,10001)
     &                    iat, tagi,
     &                    jat, tagj,
     &                    kat, tagk,angle
                    else if (print_ikj) then
                      if (.not.should_print) call errquit(
     &                    'geom_prt_angles "should_print" error',
     &                    911, GEOM_ERR)
                      if (.not.geom_calc_angle(ci,ck,cj,angle))
     &                    call errquit
     &                    ('geom_prt_angles:geom_calc_angle failed',
     &                    911, GEOM_ERR)
                      num_prt =num_prt + 1
                      write(luout,10001)
     &                    iat, tagi,
     &                    kat, tagk,
     &                    jat, tagj,angle
                    else if (print_jik) then
                      if (.not.should_print) call errquit(
     &                    'geom_prt_angles "should_print" error',
     &                    911, GEOM_ERR)
                      if (.not.geom_calc_angle(cj,ci,ck,angle))
     &                    call errquit
     &                    ('geom_prt_angles:geom_calc_angle failed',
     &                    911, GEOM_ERR)
                      num_prt =num_prt + 1
                      write(luout,10001)
     &                    jat, tagj,
     &                    iat, tagi,
     &                    kat, tagk,angle
                    endif
                  endif
                enddo
              endif
            endif
          enddo
        endif
      enddo
      if (header) then
        write(luout,10002)
        write(luout,10003) num_prt
        write(luout,10004)
        write(luout,10005)
      endif
10000 format(1x,78('='),/,
     &    33x,'internuclear angles',/,1x,78('-'),/,
     &    8x,'center 1',7x,'|',
     &    7x,'center 2',7x,'|',
     &    7x,'center 3',7x,'|',
     &    '  degrees',
     &    /,1x,78('-'))
10001 format(1x,
     &    i4,1x,a16,1x,'|',
     &    i4,1x,a16,1x,'|',
     &    i4,1x,a16,1x,'|',
     &    1x,f8.2)
10002 format(1x,78('-'))
10003 format(28x,'number of included internuclear angles: ',i10)
10004 format(1x,78('='))
10005 format(/,/)
      geom_prt_angles = FT
      end
*B4-xnet:      logical function geom_print_angles(geom)
*B4-xnet:      implicit none
*B4-xnet:#include "errquit.fh"
*B4-xnet:#include "inp.fh"
*B4-xnet:#include "stdio.fh"
*B4-xnet:c::-functions
*B4-xnet:      logical geom_calc_distance
*B4-xnet:      external geom_calc_distance
*B4-xnet:      logical geom_calc_angle
*B4-xnet:      external geom_calc_angle
*B4-xnet:      logical geom_get_def_rcov
*B4-xnet:      external geom_get_def_rcov
*B4-xnet:c::-passed      
*B4-xnet:      integer geom ! [input] geometry handle
*B4-xnet:c::-local
*B4-xnet:      double precision rscale
*B4-xnet:      integer nat ! number of atoms
*B4-xnet:      integer iat ! ith atom
*B4-xnet:      integer jat ! jth atom
*B4-xnet:      integer kat ! kth atom
*B4-xnet:      double precision chg   ! charge (ignored)
*B4-xnet:      double precision ci(3) ! coords of atom i
*B4-xnet:      character*16 tagi      ! tag of atom i
*B4-xnet:      double precision cj(3) ! coords of atom j
*B4-xnet:      character*16 tagj      ! tag of atom j
*B4-xnet:      double precision ck(3) ! coords of atom k
*B4-xnet:      character*16 tagk      ! tag of atom k
*B4-xnet:      integer lmtag
*B4-xnet:      double precision dij   ! distance between atoms i and j
*B4-xnet:      double precision djk   ! distance between atoms j and k
*B4-xnet:      double precision dik   ! distance between atoms i and k
*B4-xnet:      double precision angle ! angle to be printed
*B4-xnet:      logical FF, FT         ! fortran true and false
*B4-xnet:      integer ngood          ! number of sides under threshold
*B4-xnet:      logical dij_okay       ! dij under threshold
*B4-xnet:      logical djk_okay       ! djk under threshold
*B4-xnet:      logical dik_okay       ! dik under threshold
*B4-xnet:      logical print_ijk      ! print angle i, j, k
*B4-xnet:      logical print_ikj      ! print angle i, k, j
*B4-xnet:      logical print_jik      ! print angle j, i, k
*B4-xnet:      logical should_print   ! should something be printed?
*B4-xnet:*. . . . . . . . . . . . . . ! return status of call to geom-2-element
*B4-xnet:      logical status_tagi, status_tagj, status_tagk 
*B4-xnet:      integer iatn, jatn, katn ! atomic numbers for atom i, j and k
*B4-xnet:      character*2 symi, symj, symk ! atomic symbols for atom i, j and k
*B4-xnet:      character*16 elei, elej, elek ! atomic names for atom i, j and k
*B4-xnet:*. . . . . . . . . . . . . . . . . ! covalent radii for atom i, j and k
*B4-xnet:      character*128 emsg
*B4-xnet:      double precision i_rcov, j_rcov, k_rcov 
*B4-xnet:      integer num_prt
*B4-xnet:      logical header
*B4-xnet:      integer ludbg
*B4-xnet:      logical debug
*B4-xnet:c
*B4-xnet:c initialize variables
*B4-xnet:      ludbg = 69
*B4-xnet:      debug = .false.
*B4-xnet:      header = .false.
*B4-xnet:      rscale = 1.1d00
*B4-xnet:      FF = .false.
*B4-xnet:      FT = .true.
*B4-xnet:      dij_okay = FF
*B4-xnet:      djk_okay = FF
*B4-xnet:      dik_okay = FF
*B4-xnet:      num_prt = 0
*B4-xnet:      
*B4-xnet:      if (.not.geom_ncent(geom,nat)) call errquit
*B4-xnet:     &    ('geom_print_angles: geom_ncent',911, GEOM_ERR)
*B4-xnet:
*B4-xnet:      if (nat.lt.3) then
*B4-xnet:        geom_print_angles = FT
*B4-xnet:        return
*B4-xnet:      endif
*B4-xnet:      do iat = 1,nat
*B4-xnet:        if (.not.geom_cent_get(geom,iat,tagi,ci,chg)) call errquit
*B4-xnet:     &      ('geom_print_angles: geom_cent_get:i',911, GEOM_ERR)
*B4-xnet:        status_tagi = geom_tag_to_element(tagi,symi,elei,iatn)
*B4-xnet:        if ((symi.eq.'bq').and.
*B4-xnet:     &      (.not.status_tagi))status_tagi = .true.
*B4-xnet:        if (.not.status_tagi) call errquit
*B4-xnet:     &      ('geom_print_angles:geom_tag_to_element failed:i',911,
*     &       GEOM_ERR)
*B4-xnet:        if (.not.geom_get_def_rcov(iatn,i_rcov)) call errquit
*B4-xnet:     &      ('geom_print_angles: geom_get_def_rcov failed atom i',911,
*     &       GEOM_ERR)
*B4-xnet:
*B4-xnet:        lmtag = inp_strlen(tagi)
*B4-xnet:        do jat = 1,nat
*B4-xnet:          if (iat.ne.jat) then
*B4-xnet:
*B4-xnet:            if (.not.geom_cent_get(geom,jat,tagj,cj,chg)) call errquit
*B4-xnet:     &          ('geom_print_angles:geom_cent_get:j ',911)
*B4-xnet:
*B4-xnet:            status_tagj = geom_tag_to_element(tagj,symj,elej,jatn)
*B4-xnet:            if ((symj.eq.'bq').and.
*B4-xnet:     &          (.not.status_tagj))status_tagj = .true.
*B4-xnet:            if (.not.status_tagj) call errquit
*B4-xnet:     &          ('geom_print_angles:geom_tag_to_element failed:j',911)
*B4-xnet:            if (.not.geom_get_def_rcov(jatn,j_rcov)) call errquit
*B4-xnet:     &          ('geom_print_angles: geom_get_def_rcov failed atom j',
*B4-xnet:     &          911)
*B4-xnet:            lmtag = max(lmtag,inp_strlen(tagj))
*B4-xnet:            if (.not.geom_calc_distance(ci,cj,dij)) call errquit
*B4-xnet:     &          ('geom_print_angles:geom_calc_distance:ij ',911)
*B4-xnet:
*B4-xnet:            dij_okay = dij.lt.(rscale*(i_rcov+j_rcov))
*B4-xnet:            if (dij_okay.or.debug) then
*B4-xnet:              do kat = 1,min(iat,jat)
*B4-xnet:                if (kat.ne.jat.and.kat.ne.iat) then
*B4-xnet:                  if (.not.geom_cent_get(geom,kat,tagk,ck,chg))
*B4-xnet:     &                call errquit
*B4-xnet:     &                ('geom_print_angles:geom_cent_get:k ',911)
*B4-xnet:                  status_tagk =
*B4-xnet:     &                geom_tag_to_element(tagk,symk,elek,katn)
*B4-xnet:                  if ((symk.eq.'bq').and.
*B4-xnet:     &                (.not.status_tagk))status_tagk = .true.
*B4-xnet:                  if (.not.status_tagk) then 
*B4-xnet:                    emsg = 'geom_print_angles: '//
*B4-xnet:     &                  'geom_tag_to_element failed:k'
*B4-xnet:                    call errquit(emsg,911)
*B4-xnet:                  endif
*B4-xnet:                  if (.not.geom_get_def_rcov(katn,k_rcov)) then
*B4-xnet:                    emsg = 'geom_print_angles: '//
*B4-xnet:     &                  'geom_egt_def_rcov failed atom k'
*B4-xnet:                    call errquit(emsg,911)
*B4-xnet:                  endif
*B4-xnet:                  lmtag = max(lmtag,inp_strlen(tagk))
*B4-xnet:                  
*B4-xnet:                  if (.not.geom_calc_distance(ci,ck,dik)) call errquit
*B4-xnet:     &                ('geom_print_angles:geom_calc_distance:ik ',911)
*B4-xnet:                  if (.not.geom_calc_distance(cj,ck,djk)) call errquit
*B4-xnet:     &                ('geom_print_angles:geom_calc_distance:jk ',911)
*B4-xnet:                  dik_okay = dik.lt.(rscale*(i_rcov+k_rcov))
*B4-xnet:                  djk_okay = djk.lt.(rscale*(j_rcov+k_rcov))
*B4-xnet:                  ngood = 0
*B4-xnet:                  if (dij_okay) ngood = ngood + 1
*B4-xnet:                  if (dik_okay) ngood = ngood + 1
*B4-xnet:                  if (djk_okay) ngood = ngood + 1
*B4-xnet:                  if (debug) then
*B4-xnet:                    write(ludbg,*)'**************** iat,jat,kat',
*B4-xnet:     &                  iat,jat,kat
*B4-xnet:                    write(ludbg,*)' ngood   : ',ngood
*B4-xnet:                    write(ludbg,*)' dij_okay: ',dij_okay
*B4-xnet:                    write(ludbg,*)' dik_okay: ',dik_okay
*B4-xnet:                    write(ludbg,*)' djk_okay: ',djk_okay
*B4-xnet:                    write(ludbg,*)' dij     : ',dij
*B4-xnet:                    write(ludbg,*)' dik     : ',dik
*B4-xnet:                    write(ludbg,*)' djk     : ',djk
*B4-xnet:                    write(ludbg,*)' rij     : ',rscale*(i_rcov+j_rcov)
*B4-xnet:                    write(ludbg,*)' rik     : ',rscale*(i_rcov+k_rcov)
*B4-xnet:                    write(ludbg,*)' rjk     : ',rscale*(j_rcov+k_rcov)
*B4-xnet:                  endif
*B4-xnet:*
*B4-xnet:* ngood is 0 or 1 then atoms too far apart to be interesting
*B4-xnet:*
*B4-xnet:                  print_ijk = FF ! a(ijk) = a(kji)
*B4-xnet:                  print_ikj = FF ! a(ikj) = a(jki)
*B4-xnet:                  print_jik = FF ! a(jik) = a(kji)
*B4-xnet:                  if (ngood.eq.2) then
*B4-xnet:* ngood = 2 then only one interesting angle
*B4-xnet:                    if     (dij_okay.and.dik_okay) then
*B4-xnet:                      print_jik = FT ! then angle should be j, i, k
*B4-xnet:                    elseif (dij_okay.and.djk_okay) then
*B4-xnet:                      print_ijk = FT ! then angle should be i, j, k
*B4-xnet:                    elseif (dik_okay.and.djk_okay) then
*B4-xnet:                      print_ikj = FT ! then angle should be i, k, j
*B4-xnet:                    else
*B4-xnet:                      emsg = 'geom_print_angles: '//
*B4-xnet:     &                    'should not get here 1'
*B4-xnet:                      call errquit(emsg,911)
*B4-xnet:                    endif
*B4-xnet:                  elseif (ngood.eq.3) then
*B4-xnet:                    
*B4-xnet:* if isocoles print angle between equal sides
*B4-xnet:                    if (dij.eq.djk) then
*B4-xnet:                      print_ijk = FT
*B4-xnet:                    else if (dij.eq.dik) then
*B4-xnet:                      print_jik = FT
*B4-xnet:                    else if (djk.eq.dik) then
*B4-xnet:                      print_ikj = FT
*B4-xnet:                      
*B4-xnet:* print angle with largest value.
*B4-xnet:                    else if (dij.gt.djk.and.dij.gt.dik) then
*B4-xnet:                      print_ikj = FT
*B4-xnet:                    else if (djk.gt.dij.and.djk.gt.dik) then
*B4-xnet:                      print_jik = FT
*B4-xnet:                    else if (dik.gt.dij.and.dik.gt.djk) then
*B4-xnet:                      print_ijk = FT
*B4-xnet:                    else
*B4-xnet:                      emsg = 'geom_print_angles: '//
*B4-xnet:     &                    'should not get here 2'
*B4-xnet:                      call errquit(emsg,911)
*B4-xnet:                    endif
*B4-xnet:                  endif
*B4-xnet:                  should_print = (ngood.eq.2.or.ngood.eq.3) .and.
*B4-xnet:     &                (print_ijk.or.print_ikj.or.print_jik)
*B4-xnet:                  if (should_print.and.(.not.header)) then
*B4-xnet:                    write(luout,10000)
*B4-xnet:                    header = .true.
*B4-xnet:                  endif
*B4-xnet:                  if (print_ijk) then
*B4-xnet:                    if (.not.should_print) call errquit(
*B4-xnet:     &                  'geom_print_angles "should_print" error',911)
*B4-xnet:                    if (.not.geom_calc_angle(ci,cj,ck,angle))
*B4-xnet:     &                  call errquit
*B4-xnet:     &                  ('geom_print_angles:geom_calc_angle failed',
*B4-xnet:     &                  911)
*B4-xnet:                    num_prt =num_prt + 1
*B4-xnet:                    write(luout,10001)num_prt,
*B4-xnet:     &                  iat, tagi,
*B4-xnet:     &                  jat, tagj,
*B4-xnet:     &                  kat, tagk,angle
*B4-xnet:                  else if (print_ikj) then
*B4-xnet:                    if (.not.should_print) call errquit(
*B4-xnet:     &                  'geom_print_angles "should_print" error',911)
*B4-xnet:                    if (.not.geom_calc_angle(ci,ck,cj,angle))
*B4-xnet:     &                  call errquit
*B4-xnet:     &                  ('geom_print_angles:geom_calc_angle failed',
*B4-xnet:     &                  911)
*B4-xnet:                    num_prt =num_prt + 1
*B4-xnet:                    write(luout,10001)num_prt,
*B4-xnet:     &                  iat, tagi,
*B4-xnet:     &                  kat, tagk,
*B4-xnet:     &                  jat, tagj,angle
*B4-xnet:                  else if (print_jik) then
*B4-xnet:                    if (.not.should_print) call errquit(
*B4-xnet:     &                  'geom_print_angles "should_print" error',911)
*B4-xnet:                    if (.not.geom_calc_angle(cj,ci,ck,angle))
*B4-xnet:     &                  call errquit
*B4-xnet:     &                  ('geom_print_angles:geom_calc_angle failed',
*B4-xnet:     &                  911)
*B4-xnet:                    num_prt =num_prt + 1
*B4-xnet:                    write(luout,10001)num_prt,
*B4-xnet:     &                  jat, tagj,
*B4-xnet:     &                  iat, tagi,
*B4-xnet:     &                  kat, tagk,angle
*B4-xnet:                  endif
*B4-xnet:                endif
*B4-xnet:              enddo
*B4-xnet:            endif
*B4-xnet:          endif
*B4-xnet:        enddo
*B4-xnet:      enddo
*B4-xnet:      if (header) write(luout,10002)
*B4-xnet:10000 format(1x,86('='),/,
*B4-xnet:     &    33x,'internuclear angles',/,1x,86('-'),/,
*B4-xnet:     &    1x,'count |',
*B4-xnet:     &    7x,'center 1',7x,'|',
*B4-xnet:     &    7x,'center 2',7x,'|',
*B4-xnet:     &    7x,'center 3',7x,'|',
*B4-xnet:     &    '  degrees',
*B4-xnet:     &    /,1x,86('-'))
*B4-xnet:10001 format(1x,i5,1x,'|',
*B4-xnet:     &    i4,1x,a16,1x,'|',
*B4-xnet:     &    i4,1x,a16,1x,'|',
*B4-xnet:     &    i4,1x,a16,1x,'|',
*B4-xnet:     &    1x,f8.2)
*B4-xnet:10002 format(1x,86('='),/,/)
*B4-xnet:      geom_print_angles = FT
*B4-xnet:      end
      logical function geom_print_dihedrals(geom)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      logical geom_ncent
      logical geom_prt_dihedrals
      external geom_ncent
      external geom_prt_dihedrals
      integer geom
      integer nat
**    integer max_netp
**    parameter (max_netp=24)
      integer max_net
      integer h_xnet, k_xnet, h_xlist, k_xlist
*
      if (.not.geom_ncent(geom,nat)) call errquit
     &    ('geom_print_dihedrals: geom_ncent',911, GEOM_ERR)

* 24 seems to break      max_net = min(max_netp,nat)
      max_net = nat
      if (.not.ma_push_get(mt_int,(max_net*nat),'p_xnet',
     &    h_xnet,k_xnet)) call errquit(
     &    'geom_print_dihedrals: ma get xnet failed',911, MA_ERR)

      if (.not.ma_push_get(mt_int,(nat),'p_xlist',
     &    h_xlist,k_xlist)) call errquit(
     &    'geom_print_dihedrals: ma get xlist failed',911, MA_ERR)

      geom_print_dihedrals =
     &    geom_prt_dihedrals(geom,nat,max_net,
     &    int_mb(k_xnet),int_mb(k_xlist))
      geom_print_dihedrals = geom_print_dihedrals .and.
     &    ma_pop_stack(h_xlist)
      geom_print_dihedrals = geom_print_dihedrals .and.
     &    ma_pop_stack(h_xnet)
      end
      logical function geom_prt_dihedrals(geom,nat,max_net,xnet,xlist)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "inp.fh"
c::-functions
      logical geom_calc_distance
      logical geom_calc_dihedral
      logical geom_get_def_rcov
      logical geom_cent_get
      logical geom_tag_to_element
      external geom_calc_distance
      external geom_calc_dihedral
      external geom_get_def_rcov
      external geom_cent_get
      external geom_tag_to_element
c::-passed      
      integer geom ! [input] geometry handle
      integer nat ! number of atoms
      integer max_net
      integer xlist(nat), xnet(max_net,nat)
c::-local
      double precision rscale, tscale
      integer iat ! ith atom
      integer jat ! jth atom
      integer kat ! kth atom
      integer lat ! lth atom
      integer ipat,jpat,kpat,lpat
      double precision chg   ! charge (ignored)
      double precision ci(3),pci(3) ! coords of atom i
      character*16 tagi      ! tag of atom i
      character*8  ptagi     ! tag of atom i
      double precision cj(3),pcj(3) ! coords of atom j
      character*16 tagj      ! tag of atom j
      character*8  ptagj     ! tag of atom j
      double precision ck(3),pck(3) ! coords of atom k
      character*16 tagk      ! tag of atom k
      character*8  ptagk     ! tag of atom k
      double precision cl(3),pcl(3) ! coords of atom k
      character*16 tagl      ! tag of atom k
      character*8  ptagl     ! tag of atom k
*      double precision c_all(3,4) ! all coords
*      double precision dall(6) ! all distances
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
*rak:      logical all_okay
      logical switch_jk
c
      logical status_tagi, status_tagj, status_tagk, status_tagl
      character*2 symi, symj, symk, syml
      character*16 elei, elej, elek, elel
      integer iatn, jatn, katn, latn
      integer itmp, jtmp, ktmp, ltmp
      double precision i_rcov, j_rcov, k_rcov, l_rcov
c
*      integer ngood
      integer num_pos
      integer num_prt
      logical header
*      
      FF = .false.
      FT = .true.
      num_pos = nat*(nat-1)*(nat-2)*(nat-3)/24
      geom_prt_dihedrals = FF
      if (nat.lt.4) then
        geom_prt_dihedrals = FT        
        return
      endif
c initialize variables
      rscale = 1.1d00
      tscale = 1.1d00
      dij_okay = FF      ! dij under threshold
      dik_okay = FF      ! dik under threshold
      dil_okay = FF      ! dil under threshold
      djk_okay = FF      ! djk under threshold
      djl_okay = FF      ! djl under threshold
      dkl_okay = FF      ! dkl under threshold
      header = FF
      num_prt = 0
c
      call ifill((max_net*nat),0,xnet,1)
      call ifill(nat,0,xlist,1)
      do iat = 1,nat
        if (.not.geom_cent_get(geom,iat,tagi,ci,chg)) call errquit
     &      ('geom_prt_angles: geom_cent_get:i',911, GEOM_ERR)
        status_tagi = geom_tag_to_element(tagi,symi,elei,iatn)
        if ((symi.eq.'bq').and.
     &      (.not.status_tagi))status_tagi = .true.
        if (.not.status_tagi) call errquit
     &      ('geom_prt_angles:geom_tag_to_element failed:i',911,
     &       GEOM_ERR)
        if (.not.geom_get_def_rcov(iatn,i_rcov)) call errquit
     &      ('geom_prt_angles: geom_get_def_rcov failed atom i',911,
     &       GEOM_ERR)
        do jat = 1,nat
          
          if (iat.ne.jat) then
            
            if (.not.geom_cent_get(geom,jat,tagj,cj,chg)) call errquit
     &          ('geom_prt_angles:geom_cent_get:j ',911, GEOM_ERR)
            
            status_tagj = geom_tag_to_element(tagj,symj,elej,jatn)
            if ((symj.eq.'bq').and.
     &          (.not.status_tagj))status_tagj = .true.
            if (.not.status_tagj) call errquit
     &          ('geom_prt_angles:geom_tag_to_element failed:j',911,
     &       GEOM_ERR)
            if (.not.geom_get_def_rcov(jatn,j_rcov)) call errquit
     &          ('geom_prt_angles: geom_get_def_rcov failed atom j',
     &          911, GEOM_ERR)
            if (.not.geom_calc_distance(ci,cj,dij)) call errquit
     &          ('geom_prt_angles:geom_calc_distance:ij ',911, GEOM_ERR)
            
            if (dij.lt.(rscale*(i_rcov+j_rcov))) then
              itmp = xlist(iat) + 1
              if(itmp.gt.max_net) call errquit(
     &            'geom_prt_angles:max_net is too small ',max_net,
     &       GEOM_ERR)
              xlist(iat) = itmp
              xnet(itmp,iat) = jat
            endif
          endif
        enddo
      enddo
*rak:      write(LuOut,*)' xlist: ', xlist
*rak:      do iat = 1,nat
*rak:        write(LuOut,*)' xnet: ',iat,':',(xnet(jat,iat),jat=1,max_net)
*rak:      enddo
*rak:      write(LuOut,*)'b4 dih loop'
*rak:      itmp = 0
*rak:      do iat = 1,nat
*rak:        do jtmp = 1,xlist(iat)
*rak:          jat = xnet(jtmp,iat)
*rak:          if (iat.ne.jat) then
*rak:            do ktmp = jtmp+1,xlist(iat)
*rak:              kat = xnet(ktmp,iat)
*rak:              if (kat.ne.jat.and.kat.ne.iat) then
*rak:                do ltmp = ktmp + 1,xlist(iat)
*rak:                  lat = xnet(ltmp,iat)
*rak:                  if (lat.ne.kat.and.lat.ne.jat.and.lat.ne.iat) then
*rak:                    itmp = itmp + 1
*rak:                    write(LuOut,*)'dihang:i:  ',itmp,':',iat,jat,kat,lat
*rak:                  endif
*rak:                enddo
*rak:*rak:                do ltmp = 1,xlist(jat)
*rak:*rak:                  lat = xnet(ltmp,jat)
*rak:*rak:                  if (lat.ne.kat.and.lat.ne.jat.and.lat.ne.iat) then
*rak:*rak:                    itmp = itmp + 1
*rak:*rak:                    write(LuOut,*)'dihang:j:  ',itmp,':',iat,jat,kat,lat
*rak:*rak:                  endif
*rak:*rak:                enddo
*rak:                do ltmp = 1,xlist(kat)
*rak:                  lat = xnet(ltmp,kat)
*rak:                  if (lat.ne.kat.and.lat.ne.jat.and.lat.ne.iat) then
*rak:                    itmp = itmp + 1
*rak:                    write(LuOut,*)'dihang:k:  ',itmp,':',iat,jat,kat,lat
*rak:                  endif
*rak:                enddo
*rak:              endif
*rak:            enddo
*rak:          endif
*rak:        enddo
*rak:      enddo
*rak:      write(LuOut,*)'after dih loop'
      do iat = 1,nat
        if (.not.geom_cent_get(geom,iat,tagi,ci,chg)) call errquit
     &      ('geom_prt_dihedrals:geom_cent_get:i ',911, GEOM_ERR)
        status_tagi = geom_tag_to_element(tagi,symi,elei,iatn)
        if ((symi.eq.'bq').and.(.not.status_tagi))
     &      status_tagi = FT
        if (.not.status_tagi) call errquit
     &      ('geom_prt_dihedrals:tag2element failed:i',911, GEOM_ERR)
        if (.not.geom_get_def_rcov(iatn,i_rcov)) call errquit
     &      ('geom_prt_dihedrals:defrcov failed:i',911, GEOM_ERR)
        do jtmp = 1,xlist(iat)
          jat = xnet(jtmp,iat)
          if (iat.ne.jat) then
            if (.not.geom_cent_get(geom,jat,tagj,cj,chg)) call errquit
     &          ('geom_prt_dihedrals:geom_cent_get:j ',911, GEOM_ERR)
            status_tagj = geom_tag_to_element(tagj,symj,elej,jatn)
            if ((symj.eq.'bq').and.(.not.status_tagj))
     &          status_tagj = FT
            if (.not.status_tagj) call errquit
     &          ('geom_prt_dihedrals:tag2element failed:j',911,
     &       GEOM_ERR)
            if (.not.geom_get_def_rcov(jatn,j_rcov)) call errquit
     &          ('geom_prt_dihedrals:defrcov failed:j',911, GEOM_ERR)
            
            if (.not.geom_calc_distance(ci,cj,dij)) call errquit
     &          ('geom_prt_dihedrals:geom_calc_distance:ij ',911,
     &       GEOM_ERR)
            
            dij_okay = dij.lt.(rscale*(i_rcov+j_rcov))
            if (dij_okay) then
              do ktmp = jtmp+1,xlist(iat)
                kat = xnet(ktmp,iat)
                if (kat.ne.jat.and.kat.ne.iat) then
                  if (.not.geom_cent_get(geom,kat,tagk,ck,chg))
     &                call errquit
     &                ('geom_prt_dihedrals:geom_cent_get:k ',911,
     &       GEOM_ERR)
                  status_tagk =
     &                geom_tag_to_element(tagk,symk,elek,katn)
                  if ((symk.eq.'bq').and.(.not.status_tagk))
     &                status_tagk = FT
                  if (.not.status_tagk) call errquit
     &                ('geom_prt_dihedrals:tag2element failed:k',
     &                911, GEOM_ERR)
                  if (.not.geom_get_def_rcov(katn,k_rcov))
     &                call errquit
     &                ('geom_prt_dihedrals:defrcov failed:k',911,
     &       GEOM_ERR)
                  
                  if (.not.geom_calc_distance(ci,ck,dik)) call errquit
     &                ('geom_prt_dihedrals:geom_calc_distance:ik ',
     &                911, GEOM_ERR)
                  if (.not.geom_calc_distance(cj,ck,djk)) call errquit
     &                ('geom_prt_dihedrals:geom_calc_distance:jk ',
     &                911, GEOM_ERR)
                  
                  dik_okay = dik.lt.(rscale*(i_rcov+k_rcov))
                  djk_okay = djk.lt.(rscale*(j_rcov+k_rcov))
                  switch_jk = dik.lt.dij.and.dik_okay
                  do ltmp = ktmp + 1,xlist(iat)
                    lat = xnet(ltmp,iat)
                    if (lat.ne.kat.and.
     &                  lat.ne.jat.and.lat.ne.iat) then
                      if (.not.geom_cent_get(geom,lat,tagl,cl,chg))
     &                    call errquit
     &                    ('geom_prt_dihedrals:geom_cent_get:l ',
     &                    911, GEOM_ERR)
                      status_tagl =
     &                    geom_tag_to_element(tagl,syml,elel,latn)
                      if ((syml.eq.'bq').and.(.not.status_tagl))
     &                    status_tagl = FT
                      if (.not.status_tagl) call errquit
     &                    ('geom_prt_dihedrals:tag2elmnt fail:l',
     &                    911, GEOM_ERR)
                      if (.not.geom_get_def_rcov(latn,l_rcov))
     &                    call errquit
     &                    ('geom_prt_dihedrals:defrcov fail:l',
     &                    911, GEOM_ERR)
                      
                      if (.not.geom_calc_distance(ci,cl,dil))
     &                    call errquit
     &                    ('geom_prt_dihedrals:calc_distance:il',
     &                    911, GEOM_ERR)
                      if (.not.geom_calc_distance(cj,cl,djl))
     &                    call errquit
     &                    ('geom_prt_dihedrals:calc_distance:jl',
     &                    911, GEOM_ERR)
                      if (.not.geom_calc_distance(ck,cl,dkl))
     &                    call errquit
     &                    ('geom_prt_dihedrals:calc_distance:kl',
     &                    911, GEOM_ERR)
                      dil_okay = dil.lt.(rscale*(i_rcov+l_rcov))
                      djl_okay = djl.lt.
     &                    (tscale*rscale*(j_rcov+l_rcov))
                      dkl_okay = dkl.lt.
     &                    (tscale*rscale*(k_rcov+l_rcov))
                      num_prt = num_prt + 1
                      ipat = lat
                      jpat = iat
                      call dcopy(3,cl,1,pci,1)
                      call dcopy(3,ci,1,pcj,1)
                      ptagi = tagl
                      ptagj = tagi
                      if (switch_jk) then
                        kpat = kat
                        lpat = jat
                        call dcopy(3,ck,1,pck,1)
                        call dcopy(3,cj,1,pcl,1)
                        ptagk = tagk
                        ptagl = tagj
                      else
                        kpat = jat
                        lpat = kat
                        call dcopy(3,cj,1,pck,1)
                        call dcopy(3,ck,1,pcl,1)
                        ptagk = tagj
                        ptagl = tagk
                      endif
                      if (.not.geom_calc_dihedral
     &                    (pci,pcj,pck,pcl,diangle)) call errquit
     &                    ('geom_print_dih:geom_calc_dih death',
     &                    911, GEOM_ERR)
                      if (.not.header) then
                        write(luout,10000)
                        header = FT
                      endif ! .not.header
                      write(luout,10001)
     &                    ipat,ptagi,jpat,ptagj,
     &                    kpat,ptagk,lpat,ptagl,
     &                    diangle
*rak:                        write(LuOut,*)'i',pci
*rak:                        write(LuOut,*)'j',pcj
*rak:                        write(LuOut,*)'k',pck
*rak:                        write(LuOut,*)'l',pcl
*rak:                        write(LuOut,*)'dihang::i::',num_prt,':',
*rak:     &                      ipat,jpat,kpat,lpat,diangle
                    endif
                  enddo
*rak:                  do ltmp = 1,xlist(jat)
*rak:                    lat = xnet(ltmp,jat)
*rak:                    if (lat.ne.kat.and.
*rak:     &                  lat.ne.jat.and.lat.ne.iat) then
*rak:                      if (.not.geom_cent_get(geom,lat,tagl,cl,chg))
*rak:     &                    call errquit
*rak:     &                    ('geom_prt_dihedrals:geom_cent_get:l ',
*rak:     &                    911)
*rak:                      status_tagl =
*rak:     &                    geom_tag_to_element(tagl,syml,elel,latn)
*rak:                      if ((syml.eq.'bq').and.(.not.status_tagl))
*rak:     &                    status_tagl = FT
*rak:                      if (.not.status_tagl) call errquit
*rak:     &                    ('geom_prt_dihedrals:tag2elmnt fail:l',
*rak:     &                    911)
*rak:                      if (.not.geom_get_def_rcov(latn,l_rcov))
*rak:     &                    call errquit
*rak:     &                    ('geom_prt_dihedrals:defrcov fail:l',
*rak:     &                    911)
*rak:                      
*rak:                      if (.not.geom_calc_distance(ci,cl,dil))
*rak:     &                    call errquit
*rak:     &                    ('geom_prt_dihedrals:calc_distance:il',
*rak:     &                    911)
*rak:                      if (.not.geom_calc_distance(cj,cl,djl))
*rak:     &                    call errquit
*rak:     &                    ('geom_prt_dihedrals:calc_distance:jl',
*rak:     &                    911)
*rak:                      if (.not.geom_calc_distance(ck,cl,dkl))
*rak:     &                    call errquit
*rak:     &                    ('geom_prt_dihedrals:calc_distance:kl',
*rak:     &                    911)
*rak:                      dil_okay = dil.lt.(rscale*(i_rcov+l_rcov))
*rak:                      djl_okay = djl.lt.
*rak:     &                    (tscale*rscale*(j_rcov+l_rcov))
*rak:                      dkl_okay = dkl.lt.
*rak:     &                    (tscale*rscale*(k_rcov+l_rcov))
*rak:                      num_prt = num_prt + 1
*rak:                      ipat = iat
*rak:                      call dcopy(3,ci,1,pci,1)
*rak:                      ptagi = tagi
*rak:                      if (switch_jk) then
*rak:                        jpat = kat
*rak:                        kpat = jat
*rak:                        lpat = lat
*rak:                        call dcopy(3,ck,1,pcj,1)
*rak:                        call dcopy(3,cj,1,pck,1)
*rak:                        call dcopy(3,cl,1,pcl,1)
*rak:                        ptagj = tagk
*rak:                        ptagk = tagj
*rak:                        ptagl = tagl
*rak:                      else
*rak:                        jpat = jat
*rak:                        call dcopy(3,cj,1,pcj,1)
*rak:                        ptagj = tagj
*rak:                        if (djk.gt.djl) then
*rak:                          kpat = kat
*rak:                          lpat = lat
*rak:                          call dcopy(3,ck,1,pck,1)
*rak:                          call dcopy(3,cl,1,pcl,1)
*rak:                          ptagk = tagk
*rak:                          ptagl = tagl
*rak:                        else
*rak:                          kpat = lat
*rak:                          lpat = kat
*rak:                          call dcopy(3,cl,1,pck,1)
*rak:                          call dcopy(3,ck,1,pcl,1)
*rak:                          ptagk = tagl
*rak:                          ptagl = tagk
*rak:                        endif
*rak:                      endif
*rak:                      if (.not.geom_calc_dihedral
*rak:     &                    (pci,pcj,pck,pcl,diangle)) call errquit
*rak:     &                    ('geom_print_dih:geom_calc_dih death',
*rak:     &                    911)
*rak:                      if (.not.header) then
*rak:                        write(luout,10000)
*rak:                        header = FT
*rak:                      endif ! .not.header
*rak:                      write(luout,10001)num_prt,
*rak:     &                    ipat,ptagi,jpat,ptagj,
*rak:     &                    kpat,ptagk,lpat,ptagl,
*rak:     &                    diangle
*rak:*rak:                      write(LuOut,*)'i',pci
*rak:*rak:                      write(LuOut,*)'j',pcj
*rak:*rak:                      write(LuOut,*)'k',pck
*rak:*rak:                      write(LuOut,*)'l',pcl
*rak:*rak:                      write(LuOut,*)'dihang::j::',num_prt,':',
*rak:*rak:     &                    ipat,jpat,kpat,lpat,diangle
*rak:                    endif
*rak:                  enddo
                  do ltmp = 1,xlist(kat)
                    lat = xnet(ltmp,kat)
                    if (lat.ne.kat.and.
     &                  lat.ne.jat.and.lat.ne.iat) then
                      if (.not.geom_cent_get(geom,lat,tagl,cl,chg))
     &                    call errquit
     &                    ('geom_prt_dihedrals:geom_cent_get:l ',
     &                    911, GEOM_ERR)
                      status_tagl =
     &                    geom_tag_to_element(tagl,syml,elel,latn)
                      if ((syml.eq.'bq').and.(.not.status_tagl))
     &                    status_tagl = FT
                      if (.not.status_tagl) call errquit
     &                    ('geom_prt_dihedrals:tag2elmnt fail:l',
     &                    911, GEOM_ERR)
                      if (.not.geom_get_def_rcov(latn,l_rcov))
     &                    call errquit
     &                    ('geom_prt_dihedrals:defrcov fail:l',
     &                    911, GEOM_ERR)
                      
                      if (.not.geom_calc_distance(ci,cl,dil))
     &                    call errquit
     &                    ('geom_prt_dihedrals:calc_distance:il',
     &                    911, GEOM_ERR)
                      if (.not.geom_calc_distance(cj,cl,djl))
     &                    call errquit
     &                    ('geom_prt_dihedrals:calc_distance:jl',
     &                    911, GEOM_ERR)
                      if (.not.geom_calc_distance(ck,cl,dkl))
     &                    call errquit
     &                    ('geom_prt_dihedrals:calc_distance:kl',
     &                    911, GEOM_ERR)
                      dil_okay = dil.lt.(rscale*(i_rcov+l_rcov))
                      djl_okay = djl.lt.
     &                    (tscale*rscale*(j_rcov+l_rcov))
                      dkl_okay = dkl.lt.
     &                    (tscale*rscale*(k_rcov+l_rcov))
                      num_prt = num_prt + 1
                      ipat = iat
                      call dcopy(3,ci,1,pci,1)
                      ptagi = tagi
                      if (switch_jk) then
                        jpat = kat
                        call dcopy(3,ck,1,pcj,1)
                        ptagj = tagk
                        if (djk.gt.djl) then
                          kpat = jat
                          lpat = lat
                          call dcopy(3,cj,1,pck,1)
                          call dcopy(3,cl,1,pcl,1)
                          ptagk = tagj
                          ptagl = tagl
                        else
                          kpat = lat
                          lpat = jat
                          call dcopy(3,cl,1,pck,1)
                          call dcopy(3,cj,1,pcl,1)
                          ptagk = tagl
                          ptagl = tagj
                        endif
                      else
                        jpat = jat
                        kpat = kat
                        lpat = lat
                        call dcopy(3,cj,1,pcj,1)
                        call dcopy(3,ck,1,pck,1)
                        call dcopy(3,cl,1,pcl,1)
                        ptagj = tagj
                        ptagk = tagk
                        ptagl = tagl
                      endif
                      if (.not.geom_calc_dihedral
     &                    (pci,pcj,pck,pcl,diangle)) call errquit
     &                    ('geom_print_dih:geom_calc_dih death',
     &                    911, GEOM_ERR)
                      if (.not.header) then
                        write(luout,10000)
                        header = FT
                      endif ! .not.header
                      write(luout,10001)
     &                    ipat,ptagi,jpat,ptagj,
     &                    kpat,ptagk,lpat,ptagl,
     &                    diangle
*rak:                      write(LuOut,*)'i',pci
*rak:                      write(LuOut,*)'j',pcj
*rak:                      write(LuOut,*)'k',pck
*rak:                      write(LuOut,*)'l',pcl
*rak:                      write(LuOut,*)'dihang::k::',num_prt,':',
*rak:     &                    ipat,jpat,kpat,lpat,diangle
                    endif
                  enddo
                endif
              enddo
            endif
          endif
        enddo
      enddo
      if (header) then
        write(luout,10002)
        write(luout,10003)num_prt
        write(luout,10004)
        write(luout,10005)
      endif
10000 format(1x,78('='),/,
     &    29x,'internuclear dihedral angles',/,1x,78('-'),/,
     &    4x,'center 1',3x,'|',
     &    3x,'center 2',3x,'|',
     &    3x,'center 3',3x,'|',
     &    3x,'center 4',3x,'|',
     &    '  degrees',
     &    /,1x,78('-'))
10001 format(1x,
     &    i4,1x,a8,1x,'|',
     &    i4,1x,a8,1x,'|',
     &    i4,1x,a8,1x,'|',
     &    i4,1x,a8,1x,'|',
     &    1x,f8.2)
10002 format(1x,78('-'))
10003 format(28x,'number of included dihedral angles: ',i10)
10004 format(1x,78('='))
10005 format(/,/)
      geom_prt_dihedrals = .true.
      end
*B4-xnet:      logical function geom_print_dihedrals(geom)
*B4-xnet:      implicit none
*B4-xnet:#include "errquit.fh"
*B4-xnet:#include "mafdecls.fh"
*B4-xnet:#include "stdio.fh"
*B4-xnet:#include "inp.fh"
*B4-xnet:c::-functions
*B4-xnet:      logical geom_calc_distance
*B4-xnet:      external geom_calc_distance
*B4-xnet:      logical geom_calc_dihedral
*B4-xnet:      external geom_calc_dihedral
*B4-xnet:      logical geom_get_def_rcov
*B4-xnet:      external geom_get_def_rcov
*B4-xnet:c::-passed      
*B4-xnet:      integer geom ! [input] geometry handle
*B4-xnet:c::-local
*B4-xnet:      double precision rscale, tscale
*B4-xnet:      integer nat ! number of atoms
*B4-xnet:      integer iat ! ith atom
*B4-xnet:      integer jat ! jth atom
*B4-xnet:      integer kat ! kth atom
*B4-xnet:      integer lat ! lth atom
*B4-xnet:      integer ipat,jpat,kpat,lpat
*B4-xnet:      double precision chg   ! charge (ignored)
*B4-xnet:      double precision ci(3) ! coords of atom i
*B4-xnet:      character*16 tagi      ! tag of atom i
*B4-xnet:      character*8  ptagi     ! tag of atom i
*B4-xnet:      double precision cj(3) ! coords of atom j
*B4-xnet:      character*16 tagj      ! tag of atom j
*B4-xnet:      character*8  ptagj     ! tag of atom j
*B4-xnet:      double precision ck(3) ! coords of atom k
*B4-xnet:      character*16 tagk      ! tag of atom k
*B4-xnet:      character*8  ptagk     ! tag of atom k
*B4-xnet:      double precision cl(3) ! coords of atom k
*B4-xnet:      character*16 tagl      ! tag of atom k
*B4-xnet:      character*8  ptagl     ! tag of atom k
*B4-xnet:*      double precision c_all(3,4) ! all coords
*B4-xnet:*      double precision dall(6) ! all distances
*B4-xnet:      double precision dij   ! distance between atoms i and j
*B4-xnet:      double precision dik   ! distance between atoms i and k
*B4-xnet:      double precision dil   ! distance between atoms i and l
*B4-xnet:      double precision djk   ! distance between atoms j and k
*B4-xnet:      double precision djl   ! distance between atoms j and l
*B4-xnet:      double precision dkl   ! distance between atoms k and l
*B4-xnet:      double precision diangle ! dihedral angle to be printed
*B4-xnet:      logical FF, FT         ! fortran true and false
*B4-xnet:      logical dij_okay       ! dij under threshold
*B4-xnet:      logical dik_okay       ! dik under threshold
*B4-xnet:      logical dil_okay       ! dil under threshold
*B4-xnet:      logical djk_okay       ! djk under threshold
*B4-xnet:      logical djl_okay       ! djl under threshold
*B4-xnet:      logical dkl_okay       ! dkl under threshold
*B4-xnet:      logical all_okay
*B4-xnet:      logical switch_jk
*B4-xnet:c
*B4-xnet:      logical status_tagi, status_tagj, status_tagk, status_tagl
*B4-xnet:      character*2 symi, symj, symk, syml
*B4-xnet:      character*16 elei, elej, elek, elel
*B4-xnet:      integer iatn, jatn, katn, latn
*B4-xnet:      double precision i_rcov, j_rcov, k_rcov, l_rcov
*B4-xnet:c
*B4-xnet:*      integer ngood
*B4-xnet:      integer num_pos
*B4-xnet:      integer num_prt
*B4-xnet:      logical header
*B4-xnet:*      
*B4-xnet:      if (.not.geom_ncent(geom,nat)) call errquit
*B4-xnet:     &    ('geom_print_dihedrals: geom_ncent failed',911)
*B4-xnet:
*B4-xnet:      num_pos = nat*(nat-1)*(nat-2)*(nat-3)/24
*B4-xnet:
*B4-xnet:      FF = .false.
*B4-xnet:      FT = .true.
*B4-xnet:
*B4-xnet:      geom_print_dihedrals = FF
*B4-xnet:      if (nat.lt.4) then
*B4-xnet:        geom_print_dihedrals = FT        
*B4-xnet:        return
*B4-xnet:      endif
*B4-xnet:c initialize variables
*B4-xnet:      rscale = 1.1d00
*B4-xnet:      tscale = 1.1d00
*B4-xnet:      header = FF
*B4-xnet:      dij_okay = FF      ! dij under threshold
*B4-xnet:      dik_okay = FF      ! dik under threshold
*B4-xnet:      dil_okay = FF      ! dil under threshold
*B4-xnet:      djk_okay = FF      ! djk under threshold
*B4-xnet:      djl_okay = FF      ! djl under threshold
*B4-xnet:      dkl_okay = FF      ! dkl under threshold
*B4-xnet:      num_prt = 0
*B4-xnet:c
*B4-xnet:      do iat = 1,nat
*B4-xnet:        if (.not.geom_cent_get(geom,iat,tagi,ci,chg)) call errquit
*B4-xnet:     &      ('geom_print_dihedrals:geom_cent_get:i ',911)
*B4-xnet:        status_tagi = geom_tag_to_element(tagi,symi,elei,iatn)
*B4-xnet:        if ((symi.eq.'bq').and.(.not.status_tagi))
*B4-xnet:     &      status_tagi = FT
*B4-xnet:        if (.not.status_tagi) call errquit
*B4-xnet:     &      ('geom_print_dihedrals:tag2element failed:i',911)
*B4-xnet:        if (.not.geom_get_def_rcov(iatn,i_rcov)) call errquit
*B4-xnet:     &      ('geom_print_dihedrals:defrcov failed:i',911)
*B4-xnet:        do jat = 1,nat
*B4-xnet:          if (iat.ne.jat) then
*B4-xnet:            
*B4-xnet:            if (.not.geom_cent_get(geom,jat,tagj,cj,chg)) call errquit
*B4-xnet:     &          ('geom_print_dihedrals:geom_cent_get:j ',911)
*B4-xnet:            status_tagj = geom_tag_to_element(tagj,symj,elej,jatn)
*B4-xnet:            if ((symj.eq.'bq').and.(.not.status_tagj))
*B4-xnet:     &          status_tagj = FT
*B4-xnet:            if (.not.status_tagj) call errquit
*B4-xnet:     &          ('geom_print_dihedrals:tag2element failed:j',911)
*B4-xnet:            if (.not.geom_get_def_rcov(jatn,j_rcov)) call errquit
*B4-xnet:     &          ('geom_print_dihedrals:defrcov failed:j',911)
*B4-xnet:            
*B4-xnet:            if (.not.geom_calc_distance(ci,cj,dij)) call errquit
*B4-xnet:     &          ('geom_print_dihedrals:geom_calc_distance:ij ',911)
*B4-xnet:            
*B4-xnet:            dij_okay = dij.lt.(rscale*(i_rcov+j_rcov))
*B4-xnet:            if (dij_okay) then
*B4-xnet:              do kat = 1,nat
*B4-xnet:                if (kat.ne.jat.and.kat.ne.iat) then
*B4-xnet:                  if (.not.geom_cent_get(geom,kat,tagk,ck,chg))
*B4-xnet:     &                call errquit
*B4-xnet:     &                ('geom_print_dihedrals:geom_cent_get:k ',911)
*B4-xnet:                  status_tagk =
*B4-xnet:     &                geom_tag_to_element(tagk,symk,elek,katn)
*B4-xnet:                  if ((symk.eq.'bq').and.(.not.status_tagk))
*B4-xnet:     &                status_tagk = FT
*B4-xnet:                  if (.not.status_tagk) call errquit
*B4-xnet:     &                ('geom_print_dihedrals:tag2element failed:k',
*B4-xnet:     &                911)
*B4-xnet:                  if (.not.geom_get_def_rcov(katn,k_rcov))
*B4-xnet:     &                call errquit
*B4-xnet:     &                ('geom_print_dihedrals:defrcov failed:k',911)
*B4-xnet:                  
*B4-xnet:                  if (.not.geom_calc_distance(ci,ck,dik)) call errquit
*B4-xnet:     &                ('geom_print_dihedrals:geom_calc_distance:ik ',
*B4-xnet:     &                911)
*B4-xnet:                  if (.not.geom_calc_distance(cj,ck,djk)) call errquit
*B4-xnet:     &                ('geom_print_dihedrals:geom_calc_distance:jk ',
*B4-xnet:     &                911)
*B4-xnet:
*B4-xnet:                  dik_okay = dik.lt.(rscale*(i_rcov+k_rcov))
*B4-xnet:                  djk_okay = djk.lt.(rscale*(j_rcov+k_rcov))
*B4-xnet:                  switch_jk = dik.lt.dij.and.dik_okay
*B4-xnet:                  if (djk_okay)then
*B4-xnet:                    do lat = 1,nat
*B4-xnet:                      if(lat.ne.iat.and.lat.ne.jat.and.
*B4-xnet:     &                    lat.ne.kat) then
*B4-xnet:                        if (.not.geom_cent_get(geom,lat,tagl,cl,chg))
*B4-xnet:     &                      call errquit
*B4-xnet:     &                      ('geom_print_dihedrals:geom_cent_get:l ',
*B4-xnet:     &                      911)
*B4-xnet:                        status_tagl =
*B4-xnet:     &                      geom_tag_to_element(tagl,syml,elel,latn)
*B4-xnet:                        if ((syml.eq.'bq').and.(.not.status_tagl))
*B4-xnet:     &                      status_tagl = FT
*B4-xnet:                        if (.not.status_tagl) call errquit
*B4-xnet:     &                      ('geom_print_dihedrals:tag2elmnt fail:l',
*B4-xnet:     &                      911)
*B4-xnet:                        if (.not.geom_get_def_rcov(latn,l_rcov))
*B4-xnet:     &                      call errquit
*B4-xnet:     &                      ('geom_print_dihedrals:defrcov fail:l',
*B4-xnet:     &                      911)
*B4-xnet:
*B4-xnet:                        if (.not.geom_calc_distance(ci,cl,dil))
*B4-xnet:     &                      call errquit
*B4-xnet:     &                      ('geom_print_dihedrals:calc_distance:il',
*B4-xnet:     &                      911)
*B4-xnet:                        if (.not.geom_calc_distance(cj,cl,djl))
*B4-xnet:     &                      call errquit
*B4-xnet:     &                      ('geom_print_dihedrals:calc_distance:jl',
*B4-xnet:     &                      911)
*B4-xnet:                        if (.not.geom_calc_distance(ck,cl,dkl))
*B4-xnet:     &                      call errquit
*B4-xnet:     &                      ('geom_print_dihedrals:calc_distance:kl',
*B4-xnet:     &                      911)
*B4-xnet:                        dil_okay = dil.lt.(rscale*(i_rcov+l_rcov))
*B4-xnet:                        djl_okay = djl.lt.
*B4-xnet:     &                      (tscale*rscale*(j_rcov+l_rcov))
*B4-xnet:                        dkl_okay = dkl.lt.
*B4-xnet:     &                      (tscale*rscale*(k_rcov+l_rcov))
*B4-xnet:* collect info calculate dihedral angle
*B4-xnet:                        ipat = iat
*B4-xnet:                        ptagi = tagi
*B4-xnet:                        lpat = lat
*B4-xnet:                        ptagl = tagl
*B4-xnet:                        if (switch_jk) then
*B4-xnet:                          jpat = kat
*B4-xnet:                          ptagj = tagk
*B4-xnet:                          kpat = jat
*B4-xnet:                          ptagk = tagk
*B4-xnet:                          all_okay = dij_okay.and.djk_okay.and.
*B4-xnet:     &                        djl_okay
*B4-xnet:                          if (all_okay) then
*B4-xnet:                            if (.not.geom_calc_dihedral
*B4-xnet:     &                          (ci,ck,cj,cl,diangle)) call errquit
*B4-xnet:     &                          ('geom_print_dih:geom_calc_dih death',
*B4-xnet:     &                          911)
*B4-xnet:                          endif
*B4-xnet:                        else
*B4-xnet:                          jpat = jat
*B4-xnet:                          ptagj = tagj
*B4-xnet:                          kpat = kat
*B4-xnet:                          ptagk = tagk
*B4-xnet:                          all_okay = dij_okay.and.djk_okay.and.
*B4-xnet:     &                        dkl_okay
*B4-xnet:                          if (all_okay) then
*B4-xnet:                            if (.not.geom_calc_dihedral
*B4-xnet:     &                          (ci,cj,ck,cl,diangle)) call errquit
*B4-xnet:     &                          ('geom_print_dih:geom_calc_dih death',
*B4-xnet:     &                          911)
*B4-xnet:                          endif
*B4-xnet:                        endif ! switch_jk
*B4-xnet:                        if (all_okay) then
*B4-xnet:                          num_prt = num_prt + 1
*B4-xnet:                          if (.not.header) then
*B4-xnet:                            write(luout,10000)
*B4-xnet:                            header = FT
*B4-xnet:                          endif ! .not.header
*B4-xnet:                          write(luout,10001)num_prt,
*B4-xnet:     &                        ipat,ptagi,jpat,ptagj,
*B4-xnet:     &                        kpat,ptagk,lpat,ptagl,
*B4-xnet:     &                        diangle
*B4-xnet:                        endif ! all_okay
*B4-xnet:                      endif ! lat != iat,jat,kat
*B4-xnet:                    enddo   ! lat loop
*B4-xnet:                  endif     ! djk_okay
*B4-xnet:                endif       ! kat != iat,jat
*B4-xnet:              enddo         ! kat loop
*B4-xnet:            endif           ! dij_okay
*B4-xnet:          endif             ! jat != iat
*B4-xnet:        enddo               ! jat loop
*B4-xnet:      enddo                 ! iat loop
*B4-xnet:      if (header) write(luout,10002)
*B4-xnet:10000 format(1x,86('='),/,
*B4-xnet:     &    29x,'internuclear dihedral angles',/,1x,86('-'),/,
*B4-xnet:     &    1x,'count |',
*B4-xnet:     &    3x,'center 1',3x,'|',
*B4-xnet:     &    3x,'center 2',3x,'|',
*B4-xnet:     &    3x,'center 3',3x,'|',
*B4-xnet:     &    3x,'center 4',3x,'|',
*B4-xnet:     &    '  degrees',
*B4-xnet:     &    /,1x,86('-'))
*B4-xnet:10001 format(1x,i5,1x,'|',
*B4-xnet:     &    i4,1x,a8,1x,'|',
*B4-xnet:     &    i4,1x,a8,1x,'|',
*B4-xnet:     &    i4,1x,a8,1x,'|',
*B4-xnet:     &    i4,1x,a8,1x,'|',
*B4-xnet:     &    1x,f8.2)
*B4-xnet:10002 format(1x,86('='),/,/)
*B4-xnet:      geom_print_dihedrals = .true.
*B4-xnet:      end
      logical function geom_get_def_rcov(atn,rcoval)
      implicit none
#include "errquit.fh"
c
c routine to return the default covalent radii (in a.u.) for the given
c   atomic number.   
c
c Written by: R. A. Kendall, PNNL, December 1996
c
#include "stdio.fh"
#include "nwc_const.fh"
#include "geomP.fh"
      integer atn ! [input] atomic number of element
      double precision rcoval ! [output] estimate of covalent 
                              !          radii for atom
c
      integer i
      double precision def_rcov(nelements)
C
C Data for 1-96 From "Covalent radii revisited", Cordero et al, Dalton Trans. 2832 (2008)
C data for 97-103 RA Kendall 
* Guess = 1.2*atomic: Fr<87>, Ra<88>, Ac<89>, Th<90>, Pa<91>, 
* 
* Guess = U<92> 3.000
*
* Guess = 1.2*atomic:  Np<93>, Pu<94>, Am<95>
*
* Guess = 1.3*largest cation radii: Bk<97>, Cf<98>, Es<99>, 
*                                   Fm<100>, Md<101>, No<102>, Lr<103>
*
* Added elements 104-109 with dummy values of 1.4 - KG Dyall.
*
*  Note: values in data structure are in Angstroms.
*
      data (def_rcov(i), i=1,2)
     &    /0.31D+00,0.28D+00/
      data (def_rcov(i),i=3,10)
     &    /1.28D+00,0.96D+00,0.84D+00,0.76D+00,
     &     0.71D+00,0.66D+00,0.57D+00,0.58D+00/
      data (def_rcov(i),i=11,18)
     &    /1.66D+00,1.41D+00,1.21D+00,1.11D+00,
     &     1.07D+00,1.05D+00,1.02D+00,1.06D+00/
      data (def_rcov(i),i=19,36)
     &    /2.03D+00,1.76D+00,
     &     1.70D+00,1.60D+00,1.53D+00,1.39D+00,1.39D+00,
     &     1.32D+00,1.26D+00,1.24D+00,1.32D+00,1.22D+00,
     &     1.22D+00,1.20D+00,1.19D+00,1.20D+00,1.20D+00,1.16D+00/
      data (def_rcov(i),i=37,54)
     &    /2.20D+00,1.95D+00,
     &     1.90D+00,1.75D+00,1.64D+00,1.54D+00,1.47D+00,
     &     1.46D+00,1.42D+00,1.39D+00,1.45D+00,1.44D+00,
     &     1.42D+00,1.39D+00,1.39D+00,1.38D+00,1.39D+00,1.40D+00/
      data (def_rcov(i),i=55,86)
     &    /2.44D+00,2.15D+00,
     &     2.07D+00,2.04D+00,2.03D+00,2.01D+00,1.99D+00,
     &     1.98D+00,1.98D+00,
     &     1.96D+00,1.94D+00,1.92D+00,1.92D+00,1.89D+00,
     &     1.90D+00,1.87D+00,
     &     1.87D+00,1.75D+00,1.70D+00,1.62D+00,1.51D+00,
     &     1.44D+00,1.41D+00,1.36D+00,1.36D+00,1.32D+00,
     &     1.45D+00,1.46D+00,1.48D+00,1.40D+00,1.50D+00,1.50D+00/
      data (def_rcov(i),i=87,109) /
     &     2.60d00, 2.21d00, 2.15d00, 2.06d00, 2.00d00,
     &     1.96d00, 1.90d00, 1.87d00, 1.80d00, 1.69d00,
     &     1.42d00, 1.40d00, 1.39d00, 1.38d00, 1.37d00,
     &     1.36d00, 1.34d00, 1.40d00, 1.40d00, 1.40d00,
     &     1.40d00, 1.40d00, 1.40d00/
      geom_get_def_rcov = .false.
      if (atn.eq.0) then
        rcoval = 2.0d00  ! dummy center sees lots of things?
      elseif (atn.gt.0.and.atn.le.nelements) then
        rcoval = def_rcov(atn)
      else
        write(luout,*)' geom_get_def_rcov: atomic number:',atn
        write(luout,*)' out of range 0 -> ',nelements
        call errquit('geom_get_def_rcov: fatal error',911, GEOM_ERR)
      endif
      rcoval = rcoval*angstrom_to_au
      geom_get_def_rcov = .true.
      end
c
C> \brief Create an new geometry instance
c
C> Create a new geometry instance with the specified name.
C> Return .true. if the instance was successfully created, 
C> return .false. otherwise.
      logical function geom_create(geom, name)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
#include "inp.fh"
#include "stdio.fh"
c     
      integer geom              !< [Output] the handle of the new geometry
      character*(*) name        !< [Input] the geometry name
c
      integer i,j
      external geom_data  ! This for T3D linker
c
c     Assign the next free slot for a geometry
c
      do geom = 1, max_geom
         if (.not. active(geom)) goto 10
      end do
      write(LuOut,1) name
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
      user_units(geom) = 'angstroms'
      include_bqbq(geom) = .false.
      use_primitive(geom) = .true.
      primitive_center(geom) = 'x'
c
      zmt_nizmat(geom) = 0
      zmt_nzvar(geom) = 0
      zmt_nzfrz(geom) = 0
      zmt_source(geom) = ' '
      zmt_maxtor(geom) = 100
      zmt_cvr_scaling(geom) = 0d0   ! Indicates no user zcoord input
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
*      call dfill((3*max_cent),0.0d00,coord(1,1,geom),1)
*      call dfill(max_cent,0.0d00,charge(1,geom),1)
*      call dfill(3,0.0d00,efield(1,geom),1)
*      erep(geom) = 0.0d00
*      call dfill(3,0.0d00,ndipole(1,geom),1)
*      do i = 1,max_cent
*        oecpcent(i,geom) = .false.
*      enddo
c
      end
************************************************************************
      logical function geom_disable_zmatrix(geom)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
      integer geom
      logical geom_check_handle
      external geom_check_handle
c
      geom_disable_zmatrix = geom_check_handle(geom, 'disable_zmat')
      if (geom_disable_zmatrix) then
         zmt_nizmat(geom) = 0
         zmt_nzvar(geom) = 0
         zmt_nzfrz(geom) = 0
         zmt_source(geom) = ' '
         zmt_maxtor(geom) = 100
         zmt_cvr_scaling(geom) = 0d0 ! Indicates no user zcoord input
      end if
c
      end
************************************************************************
      logical function geom_nucexps_set(geom, ncent, invnucexp)
      implicit none
#include "nwc_const.fh"
#include "stdio.fh"
#include "geomP.fh"
c
      integer geom                      ! [input] geometry handle
      integer ncent                     ! [input] number of centers
      double precision invnucexp(ncent) ! [input] inverse nuclear exponent on each center
c
      integer i
c
      logical geom_check_handle
      external geom_check_handle
c
      geom_nucexps_set = geom_check_handle(geom, 'geom_nucexps_set')
      if (.not. geom_nucexps_set) return
c
      if (ncent.le.0) then
         write(luout,*) ' geom_nucexps_set: too few centers ',ncent,
     $        names(geom)(1:lenn(geom))
         geom_nucexps_set = .false.
         return
      else if (ncent.gt.max_cent) then
         write(luout,*) ' geom_nucexps_set: too many centers ',ncent,
     $        names(geom)(1:lenn(geom))
         geom_nucexps_set = .false.
         return
      end if
c
      do i = 1, ncent
        geom_invnucexp(i,geom) = invnucexp(i)
      enddo
c
      end      
************************************************************************
      logical function geom_nucexps_get(geom, ncent, invnucexp)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
c
      integer geom                      ! [input] geometry handle
      integer ncent                     ! [input] number of centers
      double precision invnucexp(ncent) ! [output] inverse nuclear exponent on each center
c
      integer i
c
      logical geom_check_handle
      external geom_check_handle
c
      geom_nucexps_get = geom_check_handle(geom, 'geom_nucexps_get')
      if (.not. geom_nucexps_get) return
c
      ncent = ncenter(geom)
      do i = 1, ncent
        invnucexp(i) = geom_invnucexp(i,geom)
      enddo
c
      end      
************************************************************************
      logical function geom_nucexp_set(geom, icent, invnucexp)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
#include "stdio.fh"
c
      integer geom                 ! [input] geometry handle
      integer icent                ! [input] index of center for invnucexp
      double precision invnucexp   ! [input] inverse nuclear exponent on center icent
c
      logical geom_check_handle
      external geom_check_handle
c
      geom_nucexp_set = geom_check_handle(geom, 'geom_nucexp_set')
      if (.not. geom_nucexp_set) return
c
      if (icent.le.0 .or. icent.gt.ncenter(geom)) then
         write(luout,*) ' geom_nucexp_set: icent out of range',icent,
     &        ncenter(geom),names(geom)(1:lenn(geom))
         geom_nucexp_set = .false.
      else
         geom_invnucexp(icent,geom) = invnucexp
      end if
c
      return
      end      
************************************************************************
      logical function geom_nucexp_get(geom, icent, invnucexp)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
#include "stdio.fh"
c
      integer geom            ! [input] geometry handle
      integer icent           ! [input] index of center for invnucexp
      double precision invnucexp   ! [output] inverse nuclear exponent on center icent
c
      logical geom_check_handle
      external geom_check_handle
c
      geom_nucexp_get = geom_check_handle(geom, 'geom_nucexp_get')
      if (.not. geom_nucexp_get) return
c
      if (icent.le.0 .or. icent.gt.ncenter(geom)) then
         write(luout,*) ' geom_nucexp_get: icent out of range',icent,
     &        ncenter(geom),names(geom)(1:lenn(geom))
         geom_nucexp_get = .false.
      else
         invnucexp = geom_invnucexp(icent,geom)
      end if
c
      return
      end      
************************************************************************
      logical function geom_mass_to_invnucexp (mass, invnucexp)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "geomP.fh"
c
      double precision mass      ! [input] nuclear mass
      double precision invnucexp ! [output] inverse nuclear exponent
c--local
      double precision athird
c
      geom_mass_to_invnucexp = mass .gt. 0.0d0
c
      if (mass .gt. 0.0d0) then
        athird = anint(mass)**(1.0d0/3.0d0)
        if (angstrom_to_au .eq. 0.0d0) call errquit(
     &      'geom_mass_to_invnucexp:zero conversion factor',911,
     &       GEOM_ERR)
        invnucexp = ((0.836d0*athird+0.570d0)*angstrom_to_au)**2/1.5d10
      end if
c
      end      
************************************************************************
c
C> \brief Are there any finite sized nucleii
c
C> Assesses whether there are any finite sized nucleii in the specified
C> geometry instance.
c
C> \return Return .true. if there are finite size nucleii, and .false.
C> otherwise.
      logical function geom_any_finuc (geom)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
c
      integer geom !< [Input] the geometry handle
      integer i
      double precision sum
c
      sum = 0.0d0
      do i = 1,ncenter(geom)
        sum = sum+geom_invnucexp(i,geom)
      end do
      geom_any_finuc = sum .gt. 1.0d-20
      return
      end
      subroutine geom_momint0(geom,coord,natoms,ci,AI,oprint,
     ,     considerbq,lautosym)
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "geom.fh"
#include "inp.fh"
C
C     ----- CENTER AND MOMENTS OF INERTIA ----- 
C
      integer geom ! [in]
      integer natoms ! [in]
      double precision coord(3,*) ! [in]
      double precision ci(3),ai(3,3) ! [out] ctr of mass and inertua tensor
      logical oprint, considerbq,lautosym
      character*16 element
      character*16 tag
      character*2 symbol
c
      integer iat,i,j
      double precision mass,x,y,z
      integer ibq,maxbqtype,mybq,lll
      parameter(maxbqtype=20)
      character*6 tagbq(maxbqtype)
      logical lisbq
c
c
      do j=1,3
         do i=1,3
            ai(j,i)=0d0
         enddo
      enddo
      do j=1,maxbqtype
         tagbq(j)=' '
      enddo
      ibq=0
c
      if (.not.geom_center_of_mass(geom,ci)) call errquit
     &   ('geom_momint0: could not get center of mass',555, GEOM_ERR)
c
      do iat=1,natoms
         if (.not. geom_cent_tag(geom,iat,tag)) call 
     &        errquit(' momint0 hosed ',0, GEOM_ERR)
         lisbq=inp_compare(.false.,tag(1:2),'bq')
         if (considerbq.and.lisbq) then 
c
c     ahah bq
c
            if(tag(3:3).ne.' ') then
            if(.not.geom_tag_to_default_mass(tag(3:),mass))
     .              call errquit(' momint fails ',2, GEOM_ERR)
         else
            mass=0d0
         endif
         else  
            if(.not.geom_mass_get(geom, iat, mass)) call
     &           errquit(' mass_get  failed ',iat, GEOM_ERR)
c
c     assign some mass to bqs
c
            if(mass.eq.0d0.and.lisbq.and.lautosym) then
c
c     check if we alreayd have this bq
c
               lll=inp_strlen(tag)
               do j=1,ibq
                  if(tagbq(j).eq.tag(3:lll)) then
                     mybq=j
                     goto 123
                  endif
               enddo
               ibq=ibq+1
               if(ibq.gt.maxbqtype) call errquit(
     *              ' momint0: maxbqtype too small ',ibq,0)
               tagbq(ibq)=tag(3:lll)
               mybq=ibq
 123           mass=mybq*1d0
            endif
         endif
         x =coord(1,iat) - ci(1)
         y =coord(2,iat) - ci(2)
         z =coord(3,iat) - ci(3)
         ai(1,1)=ai(1,1)+mass*(y*y+z*z)
         ai(2,1)=ai(2,1)-mass* x*y
         ai(1,2)=ai(2,1)
         ai(3,1)=ai(3,1)-mass* x*z
         ai(1,3)=ai(3,1)
         ai(2,2)=ai(2,2)+mass*(x*x+z*z)
         ai(3,2)=ai(3,2)-mass* y*z
         ai(2,3)=ai(3,2)
         ai(3,3)=ai(3,3)+mass*(x*x+y*y)
      enddo
      if(oprint) then
         write(luout,9999)
         write(luout,9998) (ci(i),i=1,3)
         write(luout,9997)
         do i=1,3
         write(luout,'(3f25.12)') (ai(i,j),j=1,3)
         enddo
      endif
c
      return
 9999 format(/,1x,'center of mass',/,1x,14(1h-))
 9998 format(' x = ',f12.8,' y = ',f12.8,' z = ',f12.8)
 9997 format(/,1x,'moments of inertia (a.u.)',/,1x,18(1h-))
      end
      subroutine geom_momint(geom)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "geom.fh"
      integer geom
c
      integer natoms,l_coord,k_coord,
     ,     k_charge,l_charge,k_tag,l_tag
      logical oprint
      double precision ci(3),ai(3,3)
      oprint = .true.
c
c     print moment of inertia
c
      if ( .NOT. geom_ncent(geom, natoms) ) call errquit(
     $     'rohf: problem with call to geom_ncent', geom , GEOM_ERR)
      if (.not. ma_push_get(mt_dbl,3*natoms,'tcoords',l_coord,k_coord))
     $     call errquit('uhf_analyze: ma failed on tmp',  natoms,
     &       MA_ERR)
      if (.not. ma_push_get(mt_dbl,natoms,'coords',l_charge,k_charge))
     $     call errquit('uhf_analyze: ma failed on tmp',  natoms,
     &       MA_ERR)
      if (.not. ma_push_get(mt_byte,natoms*16,'coords',l_tag,k_tag))
     $     call errquit('uhf_analyze: ma failed on tmp',  natoms,
     &       MA_ERR)
      if (.not. geom_cart_get(geom, natoms, byte_mb(k_tag), 
     .     dbl_mb(k_coord), dbl_mb(k_charge)))
     $     call errquit('uhf_anal: geom_cent_tag failed',0,
     &       MA_ERR)
      if (.not. ma_chop_stack(l_charge))
     $     call errquit('uhf_analyze: pop failed', 0,
     &       MA_ERR)
      call geom_momint0(geom,dbl_mb(k_coord),natoms,ci,AI,oprint,
     &  .false.,.false.)
      if (.not. ma_chop_stack(l_coord))
     $     call errquit('uhf_analyze: pop failed', 0,
     &       MA_ERR)
      return
      end
c
C> \brief Converts center coordinates from Cartesian to fractional
C> coordinates
c
C> In finite systems Cartesian coordinates are in common use whereas
C> in crystal structures fractional coordinates are used. This routine
C> converts a set of Cartesian coordinates into the corresponding 
C> fractional coordinates based on a transformation that is stored
C> within the geometry instance.
c
C> \return Return .true. if the conversion was successful, and .false.
C> otherwise.
c
      logical function geom_cart_to_frac(geom, c)
      implicit none
#include "errquit.fh"
      integer geom            !< [Input] the geometry handle
      double precision c(3,*) !< [Input|Output] the center coordinates
c
      integer iat, nat, i, j
      logical geom_check_handle, geom_amatinv_get, geom_ncent
      external geom_check_handle, geom_amatinv_get
      double precision ainv(3,3), t(3)
c
      geom_cart_to_frac = geom_check_handle(geom, 'geom_cart_to_frac')
      if (.not. geom_cart_to_frac) return
      if (.not. geom_ncent(geom,nat)) 
     $     call errquit('geom_cart_to_frac: nat', 0, GEOM_ERR)
      if (.not. geom_amatinv_get(geom, ainv))
     $     call errquit('geom_cart_to_frac: ainv', 0, GEOM_ERR)
*      write(6,*) ' The amatrix inverse'
*      call output(ainv, 1, 3, 1, 3, 3, 3, 1)
c 
      do iat = 1, nat
*         write(6,*) 'c2f before ', iat, (c(i,iat),i=1,3)
         do i = 1, 3
            t(i) = 0.0d0
            do j = 1, 3
               t(i) = t(i) + ainv(i,j)*c(j,iat)
            end do
         end do
         do i = 1, 3
            c(i,iat) = t(i)
         end do
*         write(6,*) 'c2f after  ', iat, (c(i,iat),i=1,3)
      end do
c
      end
c
C> \brief Converts center coordinates from fractional to Cartesian
C> coordinates
c
C> In finite systems Cartesian coordinates are in common use whereas
C> in crystal structures fractional coordinates are used. This routine
C> converts a set of fractional coordinates into the corresponding 
C> Cartesian coordinates based on a transformation that is stored
C> within the geometry instance.
c
C> \return Return .true. if the conversion was successful, and .false.
C> otherwise.
c
      logical function geom_frac_to_cart(geom, c)
      implicit none
#include "errquit.fh"
      integer geom            !< [Input] the geometry handle
      double precision c(3,*) !< [Input|Output] the center coordinates
c
      integer iat, nat, i, j
      logical geom_check_handle, geom_amatrix_get, geom_ncent
      external geom_check_handle, geom_amatrix_get
      double precision a(3,3), t(3)
c
      geom_frac_to_cart = geom_check_handle(geom, 'geom_frac_to_cart')
      if (.not. geom_frac_to_cart) return
      if (.not. geom_ncent(geom,nat)) 
     $     call errquit('geom_frac_to_cart: nat', 0, GEOM_ERR)
      if (.not. geom_amatrix_get(geom, a))
     $     call errquit('geom_frac_to_cart: a', 0, GEOM_ERR)
c
      do iat = 1, nat
         do i = 1, 3
            t(i) = 0.0d0
            do j = 1, 3
               t(i) = t(i) + a(i,j)*c(j,iat)
            end do
         end do
         do i = 1, 3
            c(i,iat) = t(i)
         end do
      end do
c
      end
      logical function geom_grad_cart_to_frac(geom, c)
      implicit none
#include "errquit.fh"
      integer geom
      double precision c(3,*)
c
      integer iat, nat, i, j
      logical geom_check_handle, geom_amatrix_get, geom_ncent
      external geom_check_handle, geom_amatrix_get
      double precision a(3,3), t(3)
c
      geom_grad_cart_to_frac = 
     $     geom_check_handle(geom, 'geom_grad_cart_to_frac')
      if (.not. geom_grad_cart_to_frac) return
      if (.not. geom_ncent(geom,nat)) 
     $     call errquit('geom_grad_cart_to_frac: nat', 0, GEOM_ERR)
      if (.not. geom_amatrix_get(geom, a))
     $     call errquit('geom_grad_cart_to_frac: a', 0, GEOM_ERR)
c
      do iat = 1, nat
         do i = 1, 3
            t(i) = 0.0d0
            do j = 1, 3
               t(i) = t(i) + a(j,i)*c(j,iat)
            end do
         end do
         do i = 1, 3
            c(i,iat) = t(i)
         end do
      end do
c
      end
      logical function geom_makec1(geom1, geom2)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
c
c  Creates a new geometry which is like the old one, but has C1 symmetry
c
      integer geom1   ! [in] Geometry potentially with symmetry
      integer geom2   ! [out] New geometry without symmetry
      integer ncenter          ! no. of centers
      character*16 tags(nw_max_atom)
      double precision coords(3,nw_max_atom)
      double precision charge(nw_max_atom), mass(nw_max_atom)
c
      logical geom_create, geom_set_user_units 
      logical geom_cart_get, geom_cart_set
      logical geom_masses_get, geom_masses_set
      external geom_create, geom_set_user_units 
      external geom_cart_get, geom_cart_set
      external geom_masses_get, geom_masses_set
c
      geom_makec1 = .false.
c
      if (.not.geom_create(geom2,'geometrytemp'))
     &  call errquit('geom_makec1: geom_create failed',555, GEOM_ERR)
      if (.not.geom_set_user_units(geom2,'a.u.'))
     &      call errquit('geom_makec1: geom_set_user_units failed',555,
     &       GEOM_ERR)
      if (.not.geom_cart_get(geom1,ncenter,tags,coords,charge))
     &  call errquit('geom_makec1: failed to get geom1',555, GEOM_ERR)
      if (.not.geom_cart_set(geom2,ncenter,tags,coords,charge))
     &      call errquit('geom_makec1: geom_cart_set failed',555,
     &       GEOM_ERR)
      if (.not.geom_masses_get(geom1,ncenter,mass))
     &      call errquit('geom_makec1:geom_masses_get failed',555,
     &       GEOM_ERR)
      if (.not.geom_masses_set(geom2,ncenter,mass))
     &      call errquit('geom_makec1:geom_masses_set failed',555,
     &       GEOM_ERR)
c
      geom_makec1 = .true.
      return
      end  
c
C> \brief Returns whether there is an active external Bq instance.
c
C> \return Returns .true. if there is an active external Bq instance,
C> and .false. otherwise.
      function geom_extbq_on()
      implicit none
#include "bq.fh"
      logical geom_extbq_on
      geom_extbq_on = bq_on() 
      return
      end
c
C> \brief Look up the number of centers in the external Bq instance
c
C> A problem with geometries is that the maximum number of centers is
C> fixed. In particular for QM/MM calculations this is problematic 
C> as the embedding requires up to thousands of point charges. To 
C> address this issue the point charges can be stored in an instance
C> outside of the geometry. These external Bq instances require an
C> interface of their own to interact with them. This particular
C> function extracts the number of Bq centers in the active
C> external Bq instance.
C
C> \returns The number of Bq centers in the active external Bq 
C> instance.
      function geom_extbq_ncenter()
      implicit none
#include "bq.fh"
#include "errquit.fh"
      integer  geom_extbq_ncenter
c
      integer bq_handle
      integer bq_ncent
      character*32 pname

      pname = "geom_extbq_ncenter"

      if(.not.bq_get_active(bq_handle))
     >   call errquit(pname//'no active bq handle',0,0) 
      if(.not.bq_ncenter(bq_handle,bq_ncent))
     >   call errquit(pname//':no bq centers',0,0)
       
       geom_extbq_ncenter = bq_ncent

      return
      end

c
C> \brief Look up the index of the array holding the Bq charges
c
C> The charges of the Bq centers are stored in an array that can
C> be accessed through an offset in a common block (array dbl_mb in 
C> mafdecls.fh). This function returns that offset for active 
C> external Bq instance.
c
C> \returns The offset of the Bq charges for the currently
C> active Bq instance.
      function geom_extbq_charge()
      implicit none
#include "bq.fh"
#include "errquit.fh"
      integer  geom_extbq_charge
c
      integer bq_handle
      integer i_qbq
      character*32 pname

      pname = "geom_extbq_charge"

      if(.not.bq_get_active(bq_handle))
     >   call errquit(pname//':no active bq handle',0,0) 
      if(.not.bq_index_charge(bq_handle,i_qbq))
     >   call errquit(pname//':no bq coords',0,0)

      geom_extbq_charge = i_qbq
      return
      end
c
C> \brief Look up the index of the array holding the Bq coordinates
c
C> The coordinates of the Bq centers are stored in an array that can
C> be accessed through an offset in a common block (array dbl_mb in 
C> mafdecls.fh). This function returns that offset for active 
C> external Bq instance.
c
C> \returns The offset of the Bq coordinates for the currently
C> active Bq instance.
      function geom_extbq_coord()
      implicit none
#include "bq.fh"
#include "errquit.fh"
      integer  geom_extbq_coord
c
      integer bq_handle
      integer i_cbq
      character*32 pname

      pname = "geom_extbq_coord"

      if(.not.bq_get_active(bq_handle))
     >   call errquit(pname//':no active bq handle',0,0) 
      if(.not.bq_index_coord(bq_handle,i_cbq))
     >   call errquit(pname//':no bq coords',0,0)

      geom_extbq_coord = i_cbq
      return
      end

      function geom_create_from_file(in_xyz,irtdb)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "msgids.fh"
#include "global.fh"
#include "inp.fh"
#include "stdio.fh"
#include "util.fh"
      character*(*) in_xyz
      integer irtdb
      logical geom_create_from_file
c     local variables
      integer ns
      integer i,j
      integer k
      logical otitle
      integer i_t,h_t
      integer i_m,h_m
      integer i_q,h_q
      integer i_ctmp,h_ctmp
      integer atn
      character*32 pname
      character*72 title
      character*16 tag
      character*16 buf
      character*255 filename
      character*255 xyzfile
      character*255 trjfile
      character*255 message

      integer fn_xyz,fn_trj
      logical end_of_file
      logical master
      integer geom              ! handle for geometry
      character*255 geomname    ! for name of geometry

      logical geom_create,geom_print
      external geom_create,geom_print
      logical geom_tag_to_element
      external geom_tag_to_element
      logical geom_cart_set,geom_masses_set
      external geom_cart_set,geom_masses_set
      logical geom_tag_to_default_mass
      external geom_tag_to_default_mass
      logical geom_rtdb_store,geom_destroy
      external geom_rtdb_store,geom_destroy
 
      master = ga_nodeid().eq.0
      pname = "geom_create_from_file"
c
c     we assume that xyz file has a title
c     -----------------------------------
      otitle = .true.
      geom_create_from_file = .false.
c
      xyzfile = in_xyz(1:inp_strlen(in_xyz))
      call util_file_name_resolve(xyzfile, .false.)
c
      filename = in_xyz(1:inp_strlen(xyzfile))
      if(master) 
     + call util_print_centered(luout,
     + "reading external xyz file "//
     +   filename,
     +   40,.true.)
c
c     prepare files for reading/writing
c     ---------------------------------
      if(.not.util_get_io_unit(fn_xyz)) 
     >       call errquit("cannot get file number",0,0)
      filename = xyzfile
      open(fn_xyz,file=filename,form='formatted',status='old',
     $          err=133)
c
c     get number of atoms
c     ------------------
      message = " number of atoms "
      read(fn_xyz,*,err=134) ns
c
c     temporary stack memory
c     ----------------------
      if(.not.ma_push_get(mt_byte,16*ns,'t',h_t,i_t))
     + call errquit(pname//'Failed to allocate memory for t',ns,
     &       MA_ERR)

      if(.not.ma_push_get(mt_dbl,3*ns,'ctmp',h_ctmp,i_ctmp))
     + call errquit( pname//'Failed to allocate memory for ctmp',
     + 3*ns, MA_ERR)

      if(.not.ma_push_get(mt_dbl,ns,'q',h_q,i_q))
     + call errquit(pname//'Failed to allocate memory for q',ns,
     &       MA_ERR)

      if(.not.ma_push_get(mt_dbl,ns,'m',h_m,i_m))
     + call errquit('qmmm: Failed to allocate memory for m',ns,
     &       MA_ERR)

c
c     read the coords 
c     --------------------------------
      message = " title field"
      if(otitle) 
     +  read(fn_xyz,*,err=134,end=135) title

      do i=1,ns
        tag = " "
        read(fn_xyz,*,err=134,end=135) tag,
     +        (dbl_mb(i_ctmp+3*(i-1)+k-1),k=1,3)
        do j=1,16
            byte_mb(i_t+16*(i-1)+j-1)=tag(j:j)
        end do

        if (.not.
     &      geom_tag_to_default_mass(tag,dbl_mb(i_m+i-1)))
     &      call errquit(pname//'default mass failed',
     &      911, INPUT_ERR)

        if (.not.
     &      geom_tag_to_element(tag,buf,buf,atn))
     &      call errquit(pname//'default atn failed',
     &      911, INPUT_ERR)

        dbl_mb(i_q+i-1)=atn
      end do
c     call dscal(3*ns,1/cau2ang,dbl_mb(i_ctmp),1)
c
c
      geomname = "geometry"
      if (.not. geom_create(geom, geomname)) call errquit
     $     (pname//'geom_create failed !', 0, GEOM_ERR)
c
      if(.not.geom_cart_set(geom,ns,byte_mb(i_t),
     + dbl_mb(i_ctmp),dbl_mb(i_q)))
     + call errquit('qmmm: Failed to initialize geometry',0, GEOM_ERR)
c
      if(.not.geom_masses_set(geom,ns,dbl_mb(i_m)))
     + call errquit('qmmm: Failed to initialize masses',0, GEOM_ERR)
      call geom_compute_values(geom)
c
      if(.not.geom_print(geom))
     +   call errquit('qmmm: Failed to print geom',0, RTDB_ERR)
c
      if(.not.geom_rtdb_store(irtdb,geom,geomname))
     + call errquit('qmmm: Failed to store geom to rtdb',0, RTDB_ERR)

       if(.not.geom_destroy(geom))
     + call errquit('qmmm: Failed to destroy geometry',0, GEOM_ERR)

c
      if(.not.ma_pop_stack(h_m))
     & call errquit(pname//' 
     >              Failed to deallocate stack c_tmp',ns,
     &       MA_ERR)

      if(.not.ma_pop_stack(h_q))
     & call errquit(pname//' 
     >              Failed to deallocate stack c_tmp',ns,
     &       MA_ERR)

      if(.not.ma_pop_stack(h_ctmp))
     & call errquit(pname//' 
     >              Failed to deallocate stack c_tmp',ns,
     &       MA_ERR)

      if(.not.ma_pop_stack(h_t))
     & call errquit(pname//'
     >              Failed to deallocate stack i_itmp',ns,
     &       MA_ERR)

      close(fn_xyz)
      geom_create_from_file = .true.
      return

 133  call errquit(pname//'error opening/closing '//filename,0, 0)
 134  call errquit(pname//'error reading xyz file'//message,0, 0)
 135  call errquit(pname//'error end of file at'//message,0, 0)

      end

      function geom_create_from_trj(in_xyz,nf,irtdb)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "msgids.fh"
#include "global.fh"
#include "inp.fh"
#include "stdio.fh"
#include "util.fh"
      character*(*) in_xyz
      integer nf
      integer irtdb
      logical geom_create_from_trj
c     local variables
      integer ns
      integer i,j
      integer k
      logical otitle
      integer i_t,h_t
      integer i_m,h_m
      integer i_q,h_q
      integer i_ctmp,h_ctmp
      integer atn
      character*32 pname
      character*72 title
      character*16 tag
      character*16 buf
      character*255 filename
      character*255 xyzfile
      character*255 trjfile
      character*255 message

      integer fn_xyz,fn_trj
      logical end_of_file
      logical master
      integer geom              ! handle for geometry
      character*255 geomname    ! for name of geometry

      logical geom_create,geom_print
      external geom_create,geom_print
      logical geom_tag_to_element
      external geom_tag_to_element
      logical geom_cart_set,geom_masses_set
      external geom_cart_set,geom_masses_set
      logical geom_tag_to_default_mass
      external geom_tag_to_default_mass
      logical geom_rtdb_store,geom_destroy
      external geom_rtdb_store,geom_destroy
 

      master = ga_nodeid().eq.0
      pname = "geom_create_from_trj"
      geom_create_from_trj = .false.
c
c     we assume that xyz file has a title
c     -----------------------------------
      otitle = .true.
c
c      if(.not.util_xyz_nframes(in_xyz,fn_xyz))
c     >       call errquit("cannot get number of frames",0,0)
c      write(*,*) "number of frames ", fn_xyz

      xyzfile = in_xyz(1:inp_strlen(in_xyz))
      call util_file_name_resolve(xyzfile, .false.)
c
      filename = in_xyz(1:inp_strlen(xyzfile))
      if(master) 
     + call util_print_centered(luout,
     + "reading external xyz file "//
     +   filename,
     +   40,.true.)
c
c     prepare files for reading/writing
c     ---------------------------------
      if(.not.util_get_io_unit(fn_xyz)) 
     >       call errquit("cannot get file number",0,0)
      filename = xyzfile
      open(fn_xyz,file=filename,form='formatted',status='old',
     $          err=133)
c
c
c     seek frame
c     ----------
      if(.not.util_xyz_seek(fn_xyz,nf))
     >       call errquit("cannot get frame",0,0)

c
c     get number of atoms
c     ------------------
      message = " number of atoms "
      read(fn_xyz,*,err=134) ns
c
c     temporary stack memory
c     ----------------------
      if(.not.ma_push_get(mt_byte,16*ns,'t',h_t,i_t))
     + call errquit(pname//'Failed to allocate memory for t',ns,
     &       MA_ERR)

      if(.not.ma_push_get(mt_dbl,3*ns,'ctmp',h_ctmp,i_ctmp))
     + call errquit( pname//'Failed to allocate memory for ctmp',
     + 3*ns, MA_ERR)

      if(.not.ma_push_get(mt_dbl,ns,'q',h_q,i_q))
     + call errquit(pname//'Failed to allocate memory for q',ns,
     &       MA_ERR)

      if(.not.ma_push_get(mt_dbl,ns,'m',h_m,i_m))
     + call errquit('qmmm: Failed to allocate memory for m',ns,
     &       MA_ERR)


c     read the coords 
c     --------------------------------
      message = " title field"
      if(otitle) 
     +  read(fn_xyz,*,err=134,end=135) title

      do i=1,ns
        tag = " "
        read(fn_xyz,*,err=134,end=135) tag,
     +        (dbl_mb(i_ctmp+3*(i-1)+k-1),k=1,3)
        do j=1,16
            byte_mb(i_t+16*(i-1)+j-1)=tag(j:j)
        end do

        if (.not.
     &      geom_tag_to_default_mass(tag,dbl_mb(i_m+i-1)))
     &      call errquit(pname//'default mass failed',
     &      911, INPUT_ERR)

        if (.not.
     &      geom_tag_to_element(tag,buf,buf,atn))
     &      call errquit(pname//'default atn failed',
     &      911, INPUT_ERR)

        dbl_mb(i_q+i-1)=atn
      end do
c     call dscal(3*ns,1/cau2ang,dbl_mb(i_ctmp),1)
c
c
      geomname = "geometry"
      if (.not. geom_create(geom, geomname)) call errquit
     $     (pname//'geom_create failed !', 0, GEOM_ERR)
c
      if(.not.geom_cart_set(geom,ns,byte_mb(i_t),
     + dbl_mb(i_ctmp),dbl_mb(i_q)))
     + call errquit('qmmm: Failed to initialize geometry',0, GEOM_ERR)
c
      if(.not.geom_masses_set(geom,ns,dbl_mb(i_m)))
     + call errquit('qmmm: Failed to initialize masses',0, GEOM_ERR)
      call geom_compute_values(geom)
c
      if(.not.geom_print(geom))
     +   call errquit('qmmm: Failed to print geom',0, RTDB_ERR)
c
      if(.not.geom_rtdb_store(irtdb,geom,geomname))
     + call errquit('qmmm: Failed to store geom to rtdb',0, RTDB_ERR)

       if(.not.geom_destroy(geom))
     + call errquit('qmmm: Failed to destroy geometry',0, GEOM_ERR)

c
      if(.not.ma_pop_stack(h_m))
     & call errquit(pname//' 
     >              Failed to deallocate stack c_tmp',ns,
     &       MA_ERR)

      if(.not.ma_pop_stack(h_q))
     & call errquit(pname//' 
     >              Failed to deallocate stack c_tmp',ns,
     &       MA_ERR)

      if(.not.ma_pop_stack(h_ctmp))
     & call errquit(pname//' 
     >              Failed to deallocate stack c_tmp',ns,
     &       MA_ERR)

      if(.not.ma_pop_stack(h_t))
     & call errquit(pname//'
     >              Failed to deallocate stack i_itmp',ns,
     &       MA_ERR)

      close(fn_xyz)
      geom_create_from_trj = .true.
      return

 133  call errquit(pname//'error opening/closing '//filename,0, 0)
 134  call errquit(pname//'error reading xyz file'//message,0, 0)
 135  call errquit(pname//'error end of file at'//message,0, 0)

      end

C**********************************************************************

      integer function geom_get_group_number(geom)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
      integer geom
      geom_get_group_number = group_number(geom)
      return
      end
C> @}

