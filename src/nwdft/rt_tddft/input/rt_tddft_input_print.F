C
C     rt_tddft_input_print.F
C     
C     Parses input deck for rt-tddft print options.
C     
C
      subroutine rt_tddft_input_print (rtdb)
      implicit none

#include "rt_tddft.fh"
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "stdio.fh"

      
C     == Inputs ==
      integer, intent(in) :: rtdb


C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_input_print: "


C     == Variables ==
      logical got_opt
      character*255 curr_popt
      logical print_dipole
      logical print_quad
      logical print_field
      logical print_moocc
      logical print_dipcontribs
      logical print_energy
      logical print_cputime
      logical print_charge
      logical print_convergence
      logical print_s2
      

C     (zero print options if print directive found)
      print_dipole = .false.
      print_quad = .false.
      print_field = .false. 
      print_moocc = .false.
      print_dipcontribs = .false.
      print_energy = .false. 
      print_cputime = .false.
      print_charge = .false.
      print_convergence = .false.
      print_s2 = .false.
            

C     (scan over all words found after "print" directive)
      got_opt = .true.
      
      do while (got_opt)
         got_opt = inp_a (curr_popt)
         
         if ( trim(curr_popt) .eq. "") exit
         
C     Now that have an option, toggle print flag if it matches.  If no
C     match, print an error.
         if (curr_popt .eq. "dipole") then 
            print_dipole = .true.
         elseif (curr_popt .eq. "quadrupole") then 
            print_quad = .true.
         elseif (curr_popt .eq. "field") then 
            print_field = .true.
         elseif (curr_popt .eq. "moocc") then 
            print_moocc = .true.
         elseif (curr_popt .eq. "dipcontribs") then 
            print_dipcontribs = .true.
         elseif (curr_popt .eq. "energy") then 
            print_energy = .true.
         elseif (curr_popt .eq. "cputime") then 
            print_cputime = .true.
         elseif (curr_popt .eq. "charge") then 
            print_charge = .true.
         elseif (curr_popt .eq. "convergence") then 
            print_convergence = .true.
         elseif (curr_popt .eq. "s2") then 
            print_s2 = .true.
         elseif (curr_popt .eq. "*") then 
            print_dipole = .true.
            print_quad = .true.
            print_field = .true. 
            print_moocc = .true.
            print_dipcontribs = .true.
            print_energy = .true. 
            print_cputime = .true.
            print_charge = .true.
            print_convergence = .true.
            print_s2 = .true.
         else
            call errquit (pname//
     $           "unknown print option: "//trim(curr_popt)//
     $           ".  Allowed options are: dipole, field, moocc,"//
     $           " energy, cputime, charge, convergence, "//
     $           "s2, dipcontribs, *.",
     $           0, 0)
         endif
      enddo
      if (print_dipole) then
      if (.not.rtdb_put(rtdb,'rt_tddft:print_dipole',
     $     mt_log,1,print_dipole))
     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)
      endif
      if (print_quad) then
      if (.not.rtdb_put(rtdb,'rt_tddft:print_quad',
     $     mt_log,1,print_quad))
     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)
      endif
      if (print_field) then
      if (.not.rtdb_put(rtdb,'rt_tddft:print_field',
     $     mt_log,1,print_field))
     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)
      endif
      if(print_moocc) then
      if (.not.rtdb_put(rtdb,'rt_tddft:print_moocc',
     $     mt_log,1,print_moocc))
     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)
      endif
      if(print_dipcontribs) then
      if (.not.rtdb_put(rtdb,'rt_tddft:print_dipcontribs',
     $     mt_log,1,print_dipcontribs))
     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)
      endif
      if(print_energy) then
      if (.not.rtdb_put(rtdb,'rt_tddft:print_energy',
     $     mt_log,1,print_energy))
     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)
      endif
      if(print_cputime) then
      if (.not.rtdb_put(rtdb,'rt_tddft:print_cputime',
     $     mt_log,1,print_cputime))
     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)
      endif
      if(print_charge) then
      if (.not.rtdb_put(rtdb,'rt_tddft:print_charge',
     $     mt_log,1,print_charge))
     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)
      endif
      if(print_convergence) then
      if (.not.rtdb_put(rtdb,'rt_tddft:print_convergence',
     $     mt_log,1,print_convergence))
     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)
      endif
      if(print_s2) then
      if (.not.rtdb_put(rtdb,'rt_tddft:print_s2',
     $     mt_log,1,print_s2))
     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)
      endif
     
      end subroutine
      
c $Id$
