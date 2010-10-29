      logical function argos_prepare_protein(sgmnam)
c
c $Id$
c
      implicit none
c
#include "argos_prepare_common.fh"
c
      character*5 sgmnam
c
      argos_prepare_protein=.true.
c
      if(sgmnam(1:3).eq.'ALA'.or.sgmnam(1:3).eq.'ARG'.or.
     + sgmnam(1:3).eq.'ASN'.or.sgmnam(1:3).eq.'ASP'.or.
     + sgmnam(1:3).eq.'CYS'.or.sgmnam(1:3).eq.'GLN'.or.
     + sgmnam(1:3).eq.'GLU'.or.sgmnam(1:3).eq.'GLY'.or.
     + sgmnam(1:3).eq.'HIS'.or.sgmnam(1:3).eq.'ILE'.or.
     + sgmnam(1:3).eq.'LEU'.or.sgmnam(1:3).eq.'LYS'.or.
     + sgmnam(1:3).eq.'MET'.or.sgmnam(1:3).eq.'PHE'.or.
     + sgmnam(1:3).eq.'PRO'.or.sgmnam(1:3).eq.'SER'.or.
     + sgmnam(1:3).eq.'THR'.or.sgmnam(1:3).eq.'TRP'.or.
     + sgmnam(1:3).eq.'TYR'.or.sgmnam(1:3).eq.'VAL') return
      if(ffield(1:5).eq.'amber') then
      if(sgmnam(1:3).eq.'ASH'.or.sgmnam(1:3).eq.'CYX'.or.
     + sgmnam(1:3).eq.'GLH'.or.sgmnam(1:3).eq.'HID'.or.
     + sgmnam(1:3).eq.'HIE'.or.sgmnam(1:3).eq.'HIP') return
      endif
      if(ffield(1:6).eq.'charmm') then
      if(sgmnam(1:3).eq.'CYX'.or.sgmnam(1:3).eq.'HSD'.or.
     + sgmnam(1:3).eq.'HSE'.or.sgmnam(1:3).eq.'HSP') return
      endif
c
      argos_prepare_protein=.false.
      return
      end
