C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     geom_num_core -- return the number of orbitals considered to be
C     in the "core" in a system
C
C SYNOPSIS
      Logical Function Geom_Num_Core(rtdb, Geom, module, NCore)
*
* $Id: geom_numcore.F,v 1.7 2003-02-20 17:22:00 edo Exp $
*
      Implicit NONE
C
C DESCRIPTION
C     Determines the number of core orbitals in a system based on the
C     simply the constituent atoms and the standard general chemistry
C     ideas of core and valance.
C
C BUGS
C     Without any knowledge of the basis we have no choice but to 
C     assume the physical 5d, 7f shells, but if the calculation is
c     done with cartesian representations of these shells, this approach
c     will underestimate the actual number of core orbitals in the
c     calculation.
C
C INCLUDE FILES
#include "nwc_const.fh"
#include "geomP.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
C 
C EXTERNAL ROUTINES
      Logical Geom_Check_Handle, Geom_Tag_to_Element
      External Geom_Check_Handle, Geom_Tag_to_Element

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
      Integer I, AtNr ! Atomic number
      Character*2 Sy  ! Atomic symbol
      Character*16 El ! Element name
      Logical Is_At   ! Is this tag an atom or something else?
      
      Integer num_elem
      Integer tagfrz(60)
      Integer numfrz(60)
      Integer Geom     ! [in]  Geometry handle
      Integer NCore    ! [out] Number of orbitals considered core
      Integer rtdb     ! [in]  Runtime database handle
      Integer tag
      character*(*) module
      character*80 key 
      character*20 space
      Integer mlen
C
C PRIVATE DATA
      Integer PerTbl(NElements)
C
C     NOTE: This version of the table assumes 5d & 7f representations in
c     the basis.
c
c     rjh: To avoid problems with ECPs and atoms that aren't really
c     .    atoms 'cos the user is messing with their charge, simply
c     .    ignore atoms that don't have a charge = atomic number.
C
      Data PerTbl /2*0,     ! H-He, <none> 
     $    8* 1,             ! Li-Ne, 1s
     $    8* 5,             ! Na-Ar,  +2s2p
     $   18* 9,             ! K-Kr,       +3s3p
     $   18*18,             ! Rb-Xe,          +4s3d4p 
     $   32*27,             ! Cs-Rn,                +5s4d5p
     $   23*43/             ! Fr-Mt,                      +6s4f5d6p
C
      
      
      space = 'number tag elements'
      key = ' '
      mlen = inp_strlen(module) 
      write (key, '(a,'':'',a,'':'',a)') module(1:mlen),
     $    space(1:inp_strlen(space)), 'freeze by atoms'

      if(rtdb_get(rtdb, key, mt_int, 1, num_elem)) then
          
          space = 'tag freeze'
          key = ' '
          write (key, '(a,'':'',a,'':'',a)') module(1:mlen),
     $        space(1:inp_strlen(space)), 'freeze by atoms' 
      
          if(rtdb_get(rtdb, key, mt_int, num_elem, tagfrz)) then
            space = 'num freeze'
            write (key, '(a,'':'',a,'':'',a)') module(1:mlen),
     $        space(1:inp_strlen(space)), 'freeze by atoms'
            
             if(rtdb_get(rtdb, key, mt_int, num_elem, numfrz)) then
              do i = 1, num_elem
                tag = tagfrz(i)
                PerTbl(tag) = numfrz(i)
              enddo
            endif
          endif
       endif  

      Geom_Num_Core = Geom_Check_Handle(Geom, 'Geom_Num_Core') 
      If ( .NOT. Geom_Num_Core) Return
C
      NCore = 0
      Do I = 1, NCenter(Geom)
         Is_At = Geom_Tag_To_Element( Tags(i, Geom), Sy, El, AtNr)
         If ( Is_At ) then
            if (abs(charge(i,geom)-dble(atnr)) .lt. 1d-6) then
               NCore = NCore + PerTbl( AtNr )
            endif
         endif
      EndDo
C
      Return
      End
