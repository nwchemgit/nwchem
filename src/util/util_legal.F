      subroutine util_legal()
* $Id: util_legal.F,v 1.2 1996-03-07 02:25:41 d3e129 Exp $
*------------------------------------------------*
* COPYRIGHT (C) 1994, 1995, 1996,                *
*         Pacific Northwest National Laboratory, * 
*         Battelle Memorial Institute.           *
*------------------------------------------------*
*------------> All Rights Reserved <-------------*
*------------------------------------------------*
c
c routine to print legal disclaimer to stdout
c
      implicit none
#include "stdio.fh"
c
c
      integer center_pt
      logical FF, FT
c
      FF = .false.
      FT = .true.
      center_pt = 40
c
      call util_print_centered(luout,
     &    ' COPYRIGHT (C) 1994, 1995, 1996,          ',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    '    Pacific Northwest National Laboratory, ',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    '    Battelle Memorial Institute.           ',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    ' ',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    ' >>> All Rights Reserved <<< ',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    ' ',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    ' '      ,
     &    center_pt, FF)
      call util_print_centered(luout,
     &    'DISCLAIMER',
     &    center_pt, FT)
      call util_print_centered(luout,
     &    ' ',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    ' This material was prepared as an account of work',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    ' sponsored by an agency of the United States ',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    ' Government.  Neither the United States Govern-',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    ' ment nor the United States Department of Energy,',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    ' nor Battelle, nor any of their employees, MAKES',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    ' ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    ' LEGAL LIABILITY OR RESPONSIBILITY FOR THE ',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    ' ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY ',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    ' INFORMATION, APPARATUS, PRODUCT, SOFTWARE, OR ',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    ' PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE ',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    ' WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    ' ',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    ' ',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    'LIMITED USE',
     &    center_pt, FT)
      call util_print_centered(luout,
     &    ' ',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    ' This software (including any documentation) is ',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    ' being made available to you for your internal ',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    ' use only, solely for use in performance of work ',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    ' directly for the U.S. Federal Government or work',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    ' under contracts with the U.S. Department of ',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    ' Energy or other U.S. Federal Government agencies.',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    ' This software is a version which has not yet ',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    ' been evaluated and cleared for commercialization.',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    ' Adherence to this notice may be necessary for the',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    ' author, Battelle Memorial Institute, to success-',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    ' fully assert copyright in and commercialize this ',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    ' software.  This software is not intended for ',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    ' duplication or distribution to third parties ',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    ' without the permission of the Manager of Software',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    ' Products at Pacific Northwest National Labora-',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    ' tory, Richland, Washington, 99352.',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    ' ',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    ' ',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    'ACKNOWLEDGMENT',
     &    center_pt, FT)
      call util_print_centered(luout,
     &    ' This software and its documentation were produced',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    ' with Government support under Contract Number ',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    ' DE-AC06-76RLO-1830 awarded by the United States ',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    ' Department of Energy.  The Government retains a ',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    ' paid-up non-exclusive, irrevocable worldwide ',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    ' license to reproduce, prepare derivative works, ',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    ' perform publicly and display publicly by or for ',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    ' the Government, including the right to distribute',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    ' to other Government contractors.',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    ' ',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    ' ',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    ' ',
     &    center_pt, FF)
      end
      subroutine util_cite()
      implicit none
c
#include "stdio.fh"
c
      integer center_pt
      logical FF, FT
c
      FF = .false.
      FT = .true.
      center_pt = 40
c
      write(luout,*)' '
      write(luout,*)' '
      write(luout,*)' '
      call util_print_centered(luout,
     &    'If you wish to cite NWChem in the references section of ',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    'a publication, please use the following citation:',
     &    center_pt, FF)
      write(luout,*)' '
      call util_print_centered(luout,
     &    'High Performance Computational Chemistry Group, ',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    'NWChem, A Computational Chemistry Package ',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    'for Parallel Computers, ',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    'Version 1.2 (1995), ',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    'Pacific Northwest National Laboratory, ',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    'Richland, Washington 99352, ',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    'USA',
     &    center_pt, FF)
      write(luout,*)' '
      write(luout,*)' '
      write(luout,*)' '
c
      end
