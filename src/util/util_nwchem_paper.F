c
      subroutine util_nwchem_paper()
c
      implicit none
#include "stdio.fh"
c
      call util_print_centered(LuOut, 'CITATION', 40, .true.)
      call util_print_centered(LuOut, 
     &     'Please cite the following reference when publishing',
     &     40, .false.)
      call util_print_centered(LuOut, 
     &     'results obtained with NWChem:',
     &     40, .false.)
      write(LuOut,*)
      call util_print_centered(LuOut, 
     &     'M. Valiev, E.J. Bylaska, N. Govind, K. Kowalski,',
     &     40, .false.)
      call util_print_centered(LuOut, 
     &     'T.P. Straatsma, H.J.J. van Dam, D. Wang, J. Nieplocha,',
     &     40, .false.)
      call util_print_centered(LuOut, 
     &     'E. Apra, T.L. Windus, W.A. de Jong',
     &     40, .false.)
      call util_print_centered(LuOut, 
     &     '"NWChem: a comprehensive and scalable open-source',
     &     40, .false.)
      call util_print_centered(LuOut, 
     &     'solution for large scale molecular simulations"',
     &     40, .false.)
      call util_print_centered(LuOut, 
     &     'Comput. Phys. Commun. 181, 1477 (2010)',
     &     40, .false.)
      call util_print_centered(LuOut, 
     &     'doi:10.1016/j.cpc.2010.04.018',
     &     40, .false.)
c
      end 
