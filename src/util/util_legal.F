      subroutine util_legal()
* $Id: util_legal.F,v 1.23 2001-08-09 20:33:37 edo Exp $
*---------------------------------------------------------------*
* COPYRIGHT (C) 1994, 1995, 1996, 1997, 1998, 1999, 2000, 2001  *
*         Pacific Northwest National Laboratory,                *
*         Battelle Memorial Institute.                          *
*---------------------------------------------------------------*
*------------> All Rights Reserved <----------------------------*
*---------------------------------------------------------------*
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
     & 'COPYRIGHT (C) 1994, 1995, 1996, 1997, 1998, 1999, 2000, 2001',
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
#include "util.fh"
c
      integer center_pt
      logical FF, FT
c
      FF = .false.
      FT = .true.
      center_pt = 40
c
      write(luout,'(/,/)')
      call util_print_centered(luout, 'ACKNOWLEDGEMENT', center_pt, FT)
      write(luout,10000)
      call util_print_centered(luout, 'CITATION', center_pt, FT)
      write(luout,10001)
      call util_print_centered(luout, 'REFERENCES', center_pt, FT)
      write(luout,10002)
      write(luout,10101) 1
      if (util_print('references',print_high)) then
        write(luout,10102) 2
        write(luout,10103) 3
        write(luout,10104) 4
        write(luout,10105) 5
        write(luout,10106) 6
        write(luout,10107) 7
        write(luout,10108) 8
      endif
10000 format(
     &    /,12x,
     &    'Please use the following acknowledgement where appropriate ',
     &    /,12x,
     &    'for results obtained with NWChem:',
     &    /,/,12x,
     &    'High Performance Computational Chemistry Group, "NWChem, A',
     &    /,12x,
     &    'Computational Chemistry Package for Parallel Computers, ',
     &    /,12x,
     &    'Version 4.0.1" (2001), ',
     &    'Pacific Northwest National Laboratory,',
     &    /,12x,
     &    'Richland, Washington 99352-0999, USA.'
     &    /,/)
10001 format(
     &    /,12x,
     &    'Please use the following citation when publishing results',
     &    /,12x,
     &    'obtained with NWChem:'
     &    /,/,
     + 12x,'R. J. Harrison, J. A. Nichols, T. P. Straatsma, M. Dupuis,',
     + ' E. J. Bylaska,',/,
     + 12x,'G. I. Fann, T. L. Windus, E. Apra, J. Anchell, ',
     + 'D. Bernholdt, P. Borowski,',/,
     + 12x,'T. Clark, D. Clerc, H. Dachsel, B. de Jong, M. Deegan,',
     + ' K. Dyall, D. Elwood,',/,
     + 12x,'H. Fruchtl, E. Glendenning, M. Gutowski, A. Hess,',
     + ' J. Jaffe, B. Johnson, J. Ju,',/,
     + 12x,'R. Kendall, R. Kobayashi, R, Kutteh, Z. Lin,',
     + ' R. Littlefield, X. Long, B. Meng,',/,
     + 12x,'J. Nieplocha, S. Niu, M. Rosing, G. Sandrone, M. Stave,',
     + ' H. Taylor, G. Thomas,',/,
     + 12x,'J. van Lenthe, K. Wolinski, A. Wong, and Z. Zhang,',/,
     + 12x,'NWChem, A Computational Chemistry Package for Parallel ',
     &    'Computers, ',
     &    /,12x,
     &    'Version 4.0.1" (2001),',
     &    /,22x,'Pacific Northwest National Laboratory,',
     &    /,22x,'Richland, Washington 99352-0999, USA.'
     &    /,/)
10002 format(
     &    /,12x,
     &    'Where appropriate, please cite the following reference(s) ',
     &    /,12x,
     &    'when publishing results obtained with NWChem:'
     &    /,/,12x)
10101 format(1x,i5,'.',1x,
     &    'D. E. Bernholdt, E. Apra, H. A. Fruchtl, M.F. Guest, ',
     &    'R. J. Harrison, ',/,8x,'R. A. Kendall, R. A. Kutteh, ',
     &    'X. Long, J. B. Nicholas, J. A. Nichols, ',/,8x,
     &    'H. L. Taylor, A. T. Wong, G. I. Fann, R. J. Littlefi',
     &    'eld and J. Nieplocha,',/,8x,'"Parallel Computational ',
     &    'Chemistry Made Easier: The Development of NWChem,"',/,8x,
     &    'Int. J. Quantum Chem. Symposium  29, 475-483 ',
     &    '(1995).'
     &    /,/)
10102 format(1x,i5,'.',1x,
     &    'M. F. Guest, E. Apra, D. E. Bernholdt, H. A. Fruchtl, ',
     &    'R. J. Harrison, ',/,8x,'R. A. Kendall, R. A. Kutteh, ',
     &    'X. Long, J. B. Nicholas, J. A. Nichols, ',/,8x,
     &    'H. L. Taylor, A. T. Wong, G. I. Fann, R. J. Littlefield',
     &    ' and J. Nieplocha,',/,8x,'"High Performance Computatio',
     &    'nal Chemistry: NWChem and Fully Distributed',/,8x,
     &    'Parallel Applications", in Advances in Parallel Comput',
     &    'ing, 10, High',/,8x,'Performance Computing: Technology,',
     &    ' Methods, and Applications, Eds. ',/,8x,'J. Dongarra, ',
     &    'L. Gradinetti, G. Joubert, and J. Kowalik, ',/,8x,
     &    '(Elsevier Science B. V. pp. 395-427, 1995).'
     &    /,/)
10103 format(1x,i5,'.',1x,
     &    'M. F. Guest, E. Apra, D. E. Bernholdt, H. A. Fruchtl, ',
     &    'R. J. Harrison, ',/,8x,'R. A. Kendall, R. A. Kutteh, ',
     &    'J. B. Nicholas, J. A. Nichols, M. S. Stave,',
     &    /,8x,'A. T. Wong, R. J. Littlefield and J. Nieplocha, ',
     &    '"High Performance ',/,8x,'Computational Chemistry: ',
     &    'Towards Fully Distributed Parallel Algorithms",',
     &    /,8x,'in High Performance Computing: Symposium 1995, ',
     &    'Grand Challenges in ',/,8x,'Computer Simulation, ',
     &    'Adrian M. Tentner, Editor, Proceedings of the 1995 ',/,8x,
     &    'Simulation Multiconference, April 9-13, 1995, Phoenix, ',
     &    'Arizona, ',/,8x,'Simulation Councils, Inc., The Society',
     &    ' for Computer Simulation, ',/,8x,
     &    'San Diego, CA, 1995. pp. 511-521.',
     &    /,/)
10104 format(1x,i5,'.',1x,
     &    'M. F. Guest, E. Apra, D. E. Bernholdt, H. A. Fruchtl, ',
     &    'R. J. Harrison,',/,8x,'R. A. Kendall, R. A. Kutteh, ',
     &    'X. Long, J. B. Nicholas, J. A. Nichols, ',/,8x,'H. L. ',
     &    'Taylor, A. T. Wong, G. I. Fann, R. J. Littlefield and ',
     &    'J. Nieplocha,',/,8x,'"Advances in Parallel Distributed',
     &    ' Data Software; Computational Chemistry ',/,8x,
     &    'and NWChem", in Applied Parallel Computing.  Computati',
     &    'ons in Physics,',/,8x,'Chemistry, and Engineering ',
     &    'Science, Eds. J. Wasnieski, J. Dongarra, and ',/,8x,
     &    'K. Madsen, Lecture Notes in Computer Science, 1041,',/,8x,
     &    '(Springer-Verlag, Berlin, 278-294, 1996).'
     &    /,/)
10105 format(1x,i5,'.',1x,
     &    'M. F. Guest, E. Apra, D. E. Bernholdt, H. A. Fruchtl, ',
     &    'R. J. Harrison,',/,8x,'R. A. Kendall, R. A. Kutteh, ',
     &    'X. Long, J. B. Nicholas, J. A. Nichols, ',/,8x,
     &    'H. L. Taylor, A. T. Wong, G. I. Fann, R. J. Littlefield',
     &    ' and J. Nieplocha, ',/,8x,
     &    '"High-Performance Computing in Chemistry; NWChem", ',
     &    'Future ',/,8x,
     &    'Generations Computer Systems  12(4), 273, (1996).'
     &    /,/)
10106 format(1x,i5,'.',1x,
     &    'R. A. Kendall, E. Apra, D. E. Bernholdt, E. J. Bylaska, ',
     &    'M. Dupuis,',/,8x,'G. I. Fann, R. J. Harrison, J. Ju, ',
     &    'J. A. Nichols, J. Nieplocha,',/,8x,'T. P. Straatsma, ',
     &    'T. L. Windus, A. T. Wong, Computer Phys. Comm. 128, ',
     &    '260 (2000).'
     &    /,/)
10107 format(1x,i5,'.',1x,
     + 'T. P. Straatsma and V. Helms, in Molecular Dynamics on',
     + ' Parallel Computers',/,8x,
     + '2000, R. Esser et al. editors, p 70.',//)
10108 format(1x,i5,'.',1x,
     + ' T. P. Straatsma, M. Philippopoulos and J. A. McCammon,',
     + /,8x,' Computer Phys. Comm. 128, ',
     &    '377 (2000).'
     &    /,/)
c
      end
