C$Id:
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     freeze_orbitals -- freeze number of orbitals in the molecule 
C     atom, return true if successfull, and false if not
C
C SYNOPSIS
C> \ingroup geom
C> @{
c
C> \brief Retrieve the number of frozen core orbitals
C
C> Retrieves the number of frozen core orbitals for the specified
C> geometry instance to be used in a particular code module. I.e. in
C> an SCF calculation frozen orbitals are never used, but non-zero 
C> numbers of frozen orbitals are routinely used in many-body methods,
C> such as MP2 or CCSD(T).
c
C> \return The function returns .true. if this number is found, and
C> .false. otherwise.
c
      Logical Function geom_freeze (rtdb, geom, module, ncore)

      Implicit NONE
      Integer rtdb                       !< [Input] the RTDB handle
      Integer geom                       !< [Input] the geometry handle
      Character*(*) module               !< [Input] the name of the code module
      Integer ncore                      !< [Output] the number of frozen core orbitals
      Character*255 a_name               ! freeze by atoms
      Character*255 c_name               ! freeze by core
      Logical status                     ! status variable for rtdb
C
C
C
C INCLUDE FILES
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
C 
C EXTERNAL ROUTINES
      Logical Geom_Num_Core
      External Geom_Num_Core
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      a_name = module // ':frozen core:freeze by atoms'
      c_name = module // ':frozen core'

      Geom_Freeze = .true.

      if (rtdb_get(rtdb, a_name , MT_LOG, 1, status)) then
           if (.not. geom_num_core(rtdb, geom, module, ncore)) then
                Geom_Freeze = .false.
           endif
      else if (rtdb_get(rtdb, c_name, MT_INT, 1, ncore)) then
           continue
      else
           ncore = 0
      endif

      Return
      End      
C> @}
