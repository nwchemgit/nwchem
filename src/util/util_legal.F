      subroutine util_legal()
* $Id: util_legal.F,v 1.1 1996-03-07 01:23:22 d3e129 Exp $
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
      write(luout,*)' COPYRIGHT (C) 1994, 1995, 1996, '
      write(luout,*)' Pacific Northwest National Laboratory, '
      write(luout,*)' Battelle Memorial Institute.  '
      write(luout,*)' '
      write(luout,*)' >>> All Rights Reserved <<< '
      write(luout,*)' '
      write(luout,*)' '      
      write(luout,*)'                   DISCLAIMER'
      write(luout,*)' '
      write(luout,*)' This material was prepared as an account of work'
      write(luout,*)' sponsored by an agency of the United States '
      write(luout,*)' Government.  Neither the United States Govern-'
      write(luout,*)' ment nor the United States Department of Energy,'
      write(luout,*)' nor Battelle, nor any of their employees, MAKES'
      write(luout,*)' ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY'
      write(luout,*)' LEGAL LIABILITY OR RESPONSIBILITY FOR THE '
      write(luout,*)' ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY '
      write(luout,*)' INFORMATION, APPARATUS, PRODUCT, SOFTWARE, OR '
      write(luout,*)' PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE '
      write(luout,*)' WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.'
      write(luout,*)' '
      write(luout,*)' '
      write(luout,*)'                   LIMITED USE'
      write(luout,*)' '
      write(luout,*)' This software (including any documentation) is '
      write(luout,*)' being made available to you for your internal '
      write(luout,*)' use only, solely for use in performance of work '
      write(luout,*)' directly for the U.S. Federal Government or work'
      write(luout,*)' under contracts with the U.S. Department of '
      write(luout,*)' Energy or other U.S. Federal Government agencies.'
      write(luout,*)' This software is a version which has not yet '
      write(luout,*)' been evaluated and cleared for commercialization.'
      write(luout,*)' Adherence to this notice may be necessary for the'
      write(luout,*)' author, Battelle Memorial Institute, to success-'
      write(luout,*)' fully assert copyright in and commercialize this '
      write(luout,*)' software.  This software is not intended for '
      write(luout,*)' duplication or distribution to third parties '
      write(luout,*)' without the permission of the Manager of Software'
      write(luout,*)' Products at Pacific Northwest National Labora-'
      write(luout,*)' tory, Richland, Washington, 99352.'
      write(luout,*)' '
      write(luout,*)' '
      write(luout,*)'                   ACKNOWLEDGMENT'
      write(luout,*)' This software and its documentation were produced'
      write(luout,*)' with Government support under Contract Number '
      write(luout,*)' DE-AC06-76RLO-1830 awarded by the United States '
      write(luout,*)' Department of Energy.  The Government retains a '
      write(luout,*)' paid-up non-exclusive, irrevocable worldwide '
      write(luout,*)' license to reproduce, prepare derivative works, '
      write(luout,*)' perform publicly and display publicly by or for '
      write(luout,*)' the Government, including the right to distribute'
      write(luout,*)' to other Government contractors.'
      write(luout,*)' '
      write(luout,*)' '
      write(luout,*)' '
      end
