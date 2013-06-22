C
C     rt_tddft_input_excite.F
C
C     Input deck parser to rules about applying fields to geometries.
C
      subroutine rt_tddft_input_excite (rtdb, nexcites, nfields)
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "util.fh"
#include "rtdb.fh"
#include "cdft.fh"
#include "inp.fh"
#include "geomP.fh"
#include "geom.fh"
#include "bas.fh"
#include "sym.fh"
#include "rt_tddft.fh"


C     == Inputs ==
      integer, intent(in) :: rtdb
      integer, intent(in) :: nexcites
      integer, intent(in) :: nfields

C     == Parameters ==
      character(len=*), parameter :: pname = "rt_tddft_input_excite: "


C     == Variables ==
      character*255 opt
      type(rt_field_t) field
      type(rt_excite_t) excite
      integer igeom, ifield
      logical found_geom, found_field
      integer ao_indx


      if (nexcites .gt. rt_max_excites)
     $     call errquit (pname//"cannot exceed max excites", 0, 0)


C
C     Get geom info from rtdb.  Note: geomP.fh supplies ngeom_rtdb,
C     names_rtdb, max_geom_rtdb, ndipole
C
      if (.not.rtdb_get (rtdb, "geometry:ngeom",
     $     mt_int, 1, ngeom_rtdb))
     $     call errquit (pname//"failed to get ngeom from rtdb",0,0)
      
      if (ngeom_rtdb .lt. 1)
     $     call errquit (pname//"invalid ngeom_rtdb < 1",0,0)

c$$$      if (ngeom_rtdb .gt. nw_max_geoms)
c$$$     $     call errquit (pname//"invalid ngeom_rtdb > nw_max_geoms",0,0)

      if (ngeom_rtdb .gt. max_geom_rtdb)
     $     call errquit(pname//"invalid ngeom_rtdb > max_geom_rtdb",0,0)
 
      if (.not. rtdb_cget (rtdb,'geometry:names',
     $     ngeom_rtdb,names_rtdb))
     $     call errquit (pname//"failed to read names from rtdb",0,0)



C
C     Read in geom/atomic orbital to excite
C
c$$$      if (.not. inp_a (opt)) call errquit (pname//
c$$$     $     "incorrect syntax--expecting: "//
c$$$     $     "excite [geometry,ao] <target> with <field>",0,0)
c$$$
c$$$      if (inp_compare(.false., opt, "geometry")) then  !user want to excite a geometry
c$$$         excite%ao_indx = -1           !excite all AO's 
c$$$         
c$$$         if (.not. inp_a (opt)) call errquit (pname//
c$$$     $        "expecting geometry tag",0,0)
c$$$
c$$$         found_geom = .false.
c$$$         do igeom = 1, ngeom_rtdb
c$$$            if (names_rtdb(igeom) .eq. opt) then
c$$$               if (found_geom) then 
c$$$                  call errquit (pname//
c$$$     $                 "multiple matches to geom name: "//opt, 0, 0)
c$$$               else
c$$$                  found_geom = .true.
c$$$                  excite%geom_indx = igeom
c$$$               endif
c$$$            endif
c$$$         enddo
c$$$         
c$$$         if (.not. found_geom)
c$$$     $        call errquit (pname//
c$$$     $        "Specified geom to excite does not exist: "//opt, 0, 0)
c$$$
c$$$      elseif (inp_compare(.false., opt, "ao")) then  !user want to excite a particular atomic orbital
c$$$         call errquit (pname//"ao excite not implemented yet",0,0)
c$$$         excite%geom_indx = -1         !excite all geoms
c$$$         
c$$$         if (.not. inp_i (ao_indx)) call errquit (pname//
c$$$     $        "ao takes a positive integer",0,0)
c$$$         if (ao_indx < 1) call errquit (pname//
c$$$     $        "ao takes a positive integer",0,0)
c$$$         excite%ao_indx = ao_indx
c$$$
c$$$      else
c$$$         call errquit (pname//"invalid excitation target type: "//
c$$$     $        trim(opt),0,0)
c$$$      endif

      
      if (.not. inp_a (opt)) call errquit (pname//
     $     "expecting geometry tag",0,0)

      found_geom = .false.
      do igeom = 1, ngeom_rtdb
         if (names_rtdb(igeom) .eq. opt) then
            if (found_geom) then 
               call errquit (pname//
     $              "multiple matches to geom name: "//opt, 0, 0)
            else
               found_geom = .true.
               excite%geom_indx = igeom
            endif
         endif
      enddo
         
      if (.not. found_geom)
     $     call errquit (pname//
     $     "Specified geom to excite does not exist: "//opt, 0, 0)

      excite%ao_indx = -1       !excite all AO's 

C
C     Then read in field to use.
C
      if (.not. inp_a (opt)) call errquit (pname//  !get "with"
     $     "incorrect syntax--expecting: "//
     $     "excite <geom> with <field>",0,0)
      
      if (opt .ne. "with") call errquit (pname//
     $     "incorrect syntax--expecting: "//
     $     "excite <geom> with <field>",0,0)


      if (.not. inp_a (opt)) call errquit (pname// !get field
     $     "incorrect syntax--expecting: "//
     $     "excite <geom> with <field>",0,0)


      found_field = .false.
      do ifield = 1, nfields
         call rt_tddft_field_rtdb_get (rtdb, ifield, field)
         if (field%name .eq. opt) then
            if (found_field) then 
               call errquit (pname//
     $              "multiple matches to field name: "//opt, 0, 0)
            else
               found_field = .true.
               excite%field_indx = ifield
            endif
         endif
      enddo
      
      if (.not. found_field)
     $     call errquit (pname//
     $     "Specified field does not exist: "//opt, 0, 0)


C     XXX FIX GEOM TO EXCITE IF -1


c$$$C
c$$$C     Read in geometry to excite, check that it exists, and if so XXX
c$$$C
c$$$      if (.not. inp_a (opt)) call errquit (pname// !get geom
c$$$     $     "incorrect syntax--expecting: excite <geom> with <field>")
c$$$      
c$$$      found_geom = .false.
c$$$      do igeom = 1, ngeom_rtdb
c$$$         if (names_rtdb(igeom) .eq. opt) then
c$$$            if (found_geom) then 
c$$$               call errquit (pname//
c$$$     $              "multiple matches to geom name: "//opt, 0, 0)
c$$$            else
c$$$               found_geom = .true.
c$$$               excite%geom_indx = igeom
c$$$            endif
c$$$         endif
c$$$      enddo
c$$$
c$$$      if (.not. found_geom)
c$$$     $     call errquit (pname//
c$$$     $     "Specified geom to excite does not exist: "//opt, 0, 0)
c$$$
c$$$      
c$$$
c$$$C


C
C     Check for additional options
C
c$$$      ao_indx = -99
c$$$      if (inp_a (opt)) then
c$$$         if (inp_compare (.false., opt, "ao")) then
c$$$            if (.not. inp_i (ao_indx)) call errquit (pname//
c$$$     $           "ao takes a positive integer",0,0)
c$$$            if (ao_indx < 1) call errquit (pname//
c$$$     $           "ao takes a positive integer",0,0)
c$$$            excite%ao_indx = ao_indx
c$$$         endif
c$$$      endif


C
C     Save excite rule to rtdb
C
      call rt_tddft_excite_rtdb_put (rtdb, nexcites, excite)

      end subroutine


C====================================================================
C
C     Generate entry name for field rtdb stuff (hack)
C
      subroutine rt_tddft_excite_rtdb_entry_name (i, name)
      implicit none

#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "rt_tddft.fh"

      
C     == Inputs ==
      integer, intent(in) :: i

      
C     == Outputs ==
      character(len=*), intent(out) :: name   !was 17

      
C     == Parameters ==
      character(len=*), parameter :: pname =
     $     "rt_tddft_excite_rtdb_entry_name"

      
C     == Variables ==
      character*5 istring       !note length 5 limit size of int

      
      if ( (i .gt. rt_max_excites).or.(i .lt. 1) )
     $     call errquit(pname//"i must be between 1,rt_max_excites",0,0)

      if (rt_max_fields .gt. 999) call errquit(pname//
     $     "rt_max_excites too large; fix formatting", 0, 0)

      write (istring, "(i0.5)") i
      
      name = "rt_tddft:excite_"//trim(istring)//"_"

      end subroutine


C====================================================================      
C
C     Load excite into rtbd.  This is an ugly hack, but it's easier than
C     adding a custom struct to the rtdb routines.
C     
      subroutine rt_tddft_excite_rtdb_put (rtdb, i, excite)
      implicit none

#include "rt_tddft.fh"
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "stdio.fh"

      
C     == Inputs ==
      integer, intent(in)           :: rtdb
      integer, intent(in)           :: i         !index for the excite
      type(rt_excite_t), intent(in) :: excite


C     == Parameters ==
      character(len=*), parameter :: pname="rt_tddft_excite_rtdb_put: "


C     == Variables ==
      character*32 basename
      character*32 entry_name

      if ( (i .gt. rt_max_excites).or.(i .lt. 1) )
     $     call errquit(pname//"i must be between 1,rt_max_excites",0,0)

      call rt_tddft_excite_rtdb_entry_name (i, basename)

      entry_name = trim(basename) // "geom_indx"
      if (.not.rtdb_put(rtdb,entry_name,mt_int,1,excite%geom_indx))
     $     call errquit(pname//'Write failed to geom_indx rtdb',
     $     0,RTDB_ERR)

      entry_name = trim(basename) // "field_indx"
      if (.not.rtdb_put(rtdb,entry_name,mt_int,1,excite%field_indx))
     $     call errquit(pname//'Write failed to field_indx rtdb',
     $     0,RTDB_ERR)

      entry_name = trim(basename) // "ao_indx"
      if (.not.rtdb_put(rtdb,entry_name,mt_int,1,excite%ao_indx))
     $     call errquit(pname//'Write failed to ao_indx rtdb',
     $     0,RTDB_ERR)

      end subroutine


C====================================================================
C
C     Get excite from rtdb and put into struct
C
      subroutine rt_tddft_excite_rtdb_get (rtdb, i, excite)
      implicit none

#include "rt_tddft.fh"
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "stdio.fh"

      
C     == Inputs ==
      integer, intent(in) :: rtdb
      integer, intent(in) :: i                 !index for the excite

      
C     == Outputs ==
      type(rt_excite_t), intent(out) :: excite



C     == Parameters ==
      character(len=*), parameter :: pname="rt_tddft_excite_rtdb_get: "


C     == Variables ==
      character*32 basename
      character*32 entry_name


      if ( (i .gt. rt_max_excites).or.(i .lt. 1) )
     $     call errquit(pname//"i must be between 1,rt_max_excites",0,0)

      call rt_tddft_excite_rtdb_entry_name (i, basename)

      entry_name = trim(basename) // "geom_indx"
      if (.not.rtdb_get(rtdb,entry_name,mt_int,1,excite%geom_indx))
     $     call errquit(pname//'Read failed for geom_indx rtdb',
     $     0,RTDB_ERR)

      entry_name = trim(basename) // "field_indx"
      if (.not.rtdb_get(rtdb,entry_name,mt_int,1,excite%field_indx))
     $     call errquit(pname//'Read failed for field_indx rtdb',
     $     0,RTDB_ERR)

      entry_name = trim(basename) // "ao_indx"
      if (.not.rtdb_get(rtdb,entry_name,mt_int,1,excite%ao_indx))
     $     call errquit(pname//'Read failed for ao_indx rtdb',
     $     0,RTDB_ERR)
      
      end subroutine
c $Id$
