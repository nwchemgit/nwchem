      Subroutine Module_Census
      include "errquit.fh"
* $Id$
C
C     Eventually this routine should print out a list of modules or
C     tasks or something which are incorporated into the present
C     build of the program.
C
C     For now, let's not worry too much
C
      Write (6, *) 'I dont''t know which modules are in this build.'
      Return
      End
      Subroutine debtest(RTDB)
      include "errquit.fh"
* $Id$
      Integer RTDB
      Call ErrQuit(
     $   'debtest is not in this build of NWChem', 0, CAPMIS_ERR)
      Return
      End
      Subroutine gap_parse(RTDB)
      include "errquit.fh"
* $Id$
      Integer RTDB
      Call ErrQuit(
     $   'gap_parse is not in this build of NWChem', 0, CAPMIS_ERR)
      Return
      End
      Subroutine gapss(RTDB)
      include "errquit.fh"
* $Id$
      Integer RTDB
      Call ErrQuit(
     $   'gapss is not in this build of NWChem', 0, CAPMIS_ERR)
      Return
      End
      Subroutine ideaz_input(RTDB)
      include "errquit.fh"
* $Id$
      Integer RTDB
      Call ErrQuit(
     $   'ideaz_input is not in this build of NWChem', 0, CAPMIS_ERR)
      Return
      End
      Subroutine jantest(RTDB)
      include "errquit.fh"
* $Id$
      Integer RTDB
      Call ErrQuit(
     $   'jantest is not in this build of NWChem', 0, CAPMIS_ERR)
      Return
      End
      Subroutine jvltest(RTDB)
      include "errquit.fh"
* $Id$
      Integer RTDB
      Call ErrQuit(
     $   'jvltest is not in this build of NWChem', 0, CAPMIS_ERR)
      Return
      End
      Subroutine jvltest2(RTDB)
      include "errquit.fh"
* $Id$
      Integer RTDB
      Call ErrQuit(
     $   'jvltest2 is not in this build of NWChem', 0, CAPMIS_ERR)
      Return
      End
      Subroutine kgdtest(RTDB)
      include "errquit.fh"
* $Id$
      Integer RTDB
      Call ErrQuit(
     $   'kgdtest is not in this build of NWChem', 0, CAPMIS_ERR)
      Return
      End
      Subroutine mdmain(RTDB)
      include "errquit.fh"
* $Id$
      Integer RTDB
      Call ErrQuit(
     $   'mdmain is not in this build of NWChem', 0, CAPMIS_ERR)
      Return
      End
      Subroutine oimp2(RTDB)
      include "errquit.fh"
* $Id$
      Integer RTDB
      Call ErrQuit(
     $   'oimp2 is not in this build of NWChem', 0, CAPMIS_ERR)
      Return
      End
      Subroutine plnwv(RTDB)
      include "errquit.fh"
* $Id$
      Integer RTDB
      Call ErrQuit(
     $   'plnwv is not in this build of NWChem', 0, CAPMIS_ERR)
      Return
      End
      Subroutine plnwv_input(RTDB)
      include "errquit.fh"
* $Id$
      Integer RTDB
      Call ErrQuit(
     $   'plnwv_input is not in this build of NWChem', 0, CAPMIS_ERR)
      Return
      End
      Subroutine raktest(RTDB)
      include "errquit.fh"
* $Id$
      Integer RTDB
      Call ErrQuit(
     $   'raktest is not in this build of NWChem', 0, CAPMIS_ERR)
      Return
      End
      Subroutine rimp2g(RTDB)
      include "errquit.fh"
* $Id$
      Integer RTDB
      Call ErrQuit(
     $   'rimp2g is not in this build of NWChem', 0, CAPMIS_ERR)
      Return
      End
      Subroutine rjhtest(RTDB)
      include "errquit.fh"
* $Id$
      Integer RTDB
      Call ErrQuit(
     $   'rjhtest is not in this build of NWChem', 0, CAPMIS_ERR)
      Return
      End
#ifndef WITHPYTHON
      Subroutine python_input(RTDB)
      include "errquit.fh"
* $Id$
      Integer RTDB
      Call ErrQuit(
     $   'python_input is not in this build of NWChem', 0, CAPMIS_ERR)
      Return
      End
      Subroutine task_python(RTDB)
      include "errquit.fh"
* $Id$
      Integer RTDB
      Call ErrQuit(
     $   'task_python is not in this build of NWChem', 0, CAPMIS_ERR)
      Return
      End
#endif
      Subroutine task_vscf(RTDB)
      include "errquit.fh"
* $Id$
      Integer RTDB
      Call ErrQuit(
     $   'task_vscf is not in this build of NWChem', 0, CAPMIS_ERR)
      Return
      End
      Subroutine uccsdt(RTDB)
      include "errquit.fh"
* $Id$
      Integer RTDB
      Call ErrQuit(
     $   'uccsdt is not in this build of NWChem', 0, CAPMIS_ERR)
      Return
      End
      Subroutine uccsdtest(RTDB)
      include "errquit.fh"
* $Id$
      Integer RTDB
      Call ErrQuit(
     $   'uccsdtest is not in this build of NWChem', 0, CAPMIS_ERR)
      Return
      End
      Subroutine vscf_input(RTDB)
      include "errquit.fh"
* $Id$
      Integer RTDB
      Call ErrQuit(
     $   'vscf_input is not in this build of NWChem', 0, CAPMIS_ERR)
      Return
      End
* $Id$
      Subroutine mcstep(RTDB)
      include "errquit.fh"
      Integer RTDB
      Call ErrQuit(
     $   'mcstep is not in this build of NWChem', 0, CAPMIS_ERR)
      Return
      End
* $Id$
      Subroutine mcstep_main(RTDB)
      include "errquit.fh"
      Integer RTDB
      Call ErrQuit(
     $   'mcstep_main is not in this build of NWChem', 0, CAPMIS_ERR)
      Return
      End
