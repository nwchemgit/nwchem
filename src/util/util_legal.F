      subroutine util_legal()
* $Id: util_legal.F,v 1.3 1996-03-12 22:43:18 d3g681 Exp $
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
     &    center_pt+1, FF)
      call util_print_centered(luout,
     &    '    Pacific Northwest National Laboratory, ',
     &    center_pt, FF)
      call util_print_centered(luout,
     &    '    Battelle Memorial Institute.           ',
     &    center_pt, FF)
      write(luout,*)
      call util_print_centered(luout,
     &    ' >>> All Rights Reserved <<< ',
     &    center_pt, FF)
      write(luout,*)
      write(luout,*)
c
      call util_print_centered(luout,
     &    'DISCLAIMER',
     &    center_pt, FT)
      write(luout,1)
 1    format(/
     &     12x, 
     &     'This material was prepared as an account of work sponsored'/
     &     12x, 
     &     'by an agency of the United States Government.  Neither the'/
     &     12x, 
     &     'United States Government nor the United States Department'/
     &     12x, 
     &     'of Energy, nor Battelle, nor any of their employees, MAKES'/
     &     12x, 
     &     'ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL'/
     &     12x, 
     &    'LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS,'/
     &     12x, 
     &     'OR USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT,'/
     &     12x, 
     &     'SOFTWARE, OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE'/
     &     12x, 
     &     'WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.'//)
c
      call util_print_centered(luout, 'LIMITED USE', center_pt, FT)
      write(luout,2)
 2    format(/
     &     12x,
     &     'This software (including any documentation) is being made'/
     &     12x, 
     &    'available to you for your internal use only, solely for use'/
     &     12x, 
     &     'in performance of work directly for the U.S. Federal'/
     &     12x, 
     &    'Government or work under contracts with the U.S. Department'/
     &     12x, 
     &     'of Energy or other U.S. Federal Government agencies.  This'/
     &     12x, 
     &     'software is a version which has not yet been evaluated and'/
     &     12x, 
     &     'cleared for commercialization.  Adherence to this notice'/
     &     12x, 
     &     'may be necessary for the author, Battelle Memorial'/
     &     12x, 
     &     'Institute, to successfully assert copyright in and'/
     &     12x, 
     &     'commercialize this software. This software is not intended'/
     &     12x, 
     &     'for duplication or distribution to third parties without'/
     &     12x, 
     &     'the permission of the Manager of Software Products at'/
     &     12x, 
     &     'Pacific Northwest National Laboratory, Richland,'/
     &     12x, 
     &     'Washington, 99352.'//)
c
      call util_print_centered(luout, 'ACKNOWLEDGMENT', center_pt, FT)
      write(luout,3)
 3    format(/
     &     12x, 
     &     'This software and its documentation were produced with'/
     &     12x, 
     &    'Government support under Contract Number DE-AC06-76RLO-1830'/
     &     12x, 
     &     'awarded by the United States Department of Energy.  The'/
     &     12x, 
     &     'Government retains a paid-up non-exclusive, irrevocable'/
     &     12x, 
     &     'worldwide license to reproduce, prepare derivative works,'/
     &     12x, 
     &     'perform publicly and display publicly by or for the'/
     &     12x, 
     &     'Government, including the right to distribute to other'/
     &     12x, 
     &     'Government contractors.'//)
c
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
      write(luout,1)
 1    format(//
     $     1x,
     $     'If you wish to cite NWChem in the references section of a'/
     &     1x, 
     &     'publication, please use the following citation:'//
     &     1x, 
     &    '   High Performance Computational Chemistry Group, "NWChem,'/
     &     1x, 
     &     '   A Computational Chemistry Package for Parallel'/
     &     1x, 
     &     '   Computers, Version 1.2," (1995), Pacific Northwest'/
     &     1x, 
     &     '   National Laboratory, Richland, Washington 99352, USA.'//)
c
      end
