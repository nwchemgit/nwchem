c
      subroutine util_legal()
c
* $Id: util_legal.F,v 1.54 2009-03-09 22:30:58 niri Exp $
c
c routine to print legal disclaimer to stdout
c
      implicit none
#include "stdio.fh"
c
      call util_print_centered(luout, 'ACKNOWLEDGMENT', 40, .true.)
      write(luout,1)
 1    format(/
     &12x,
     &'This software and its documentation were developed at the'/
     &12x, 
     &'EMSL at Pacific Northwest National Laboratory, a multiprogram'/
     &12x, 
     &'national laboratory, operated for the U.S. Department of Energy'/
     &12x, 
     &'by Battelle under Contract Number DE-AC05-76RL01830. Support'/
     &12x, 
     &'for this work was provided by the Department of Energy Office'/
     &12x, 
     &'of Biological and Environmental Research, Office of Basic'/
     &12x, 
     &'Energy Sciences, and the Office of Advanced Scientific Computing.
     &'//)
c
      end
      subroutine util_cite()
      implicit none
c
#include "stdio.fh"
#include "util.fh"
c
      character*80 chr01,chr02,chr03,chr04,chr05,chr06,chr07,chr08
      character*80 chr09,chr10,chr11,chr12,chr13
      data chr01/
     +'E. J. Bylaska, W. A. de Jong, N. Govind, K. Kowalski, T. P. Straa
     +tsma,'/
      data chr02/
     +'M. Valiev, H. J. J. van Dam, D. Wang, E. Apra, T. L. Windus, J. H
     +ammond,'/
      data chr03/
     +'J. Autschbach, P. Nichols, S. Hirata, M. T. Hackler, Y. Zhao, P.-
     +D. Fan,'/
      data chr04/
     +'R. J. Harrison, M. Dupuis, D. M. A. Smith, K. Glaesemann, J. Niep
     +locha,'/
      data chr05/
     +'V. Tipparaju, M. Krishnan, A. Vazquez-Mayagoitia, L. Jensen, M. S
     +wart,'/
      data chr06/
     +'Q. Wu, T. Van Voorhis, A. A. Auer, M. Nooijen, L. D. Crosby, E. B
     +rown,'/
      data chr07/
     +'G. Cisneros, G. I. Fann, H. Fruchtl, J. Garza, K. Hirao,'/
      data chr08/
     +'R. Kendall, J. A. Nichols, K. Tsemekhman, K. Wolinski, J. Anchell
     +,'/
      data chr09/
     +'D. Bernholdt, P. Borowski, T. Clark, D. Clerc, H. Dachsel, M. Dee
     +gan,'/
      data chr10/
     +'K. Dyall, D. Elwood, E. Glendening, M. Gutowski, A. Hess, J. Jaff
     +e,'/
      data chr11/
     +'B. Johnson, J. Ju, R. Kobayashi, R. Kutteh, Z. Lin, R. Littlefiel
     +d,'/
      data chr12/
     +'X. Long, B. Meng, T. Nakajima, S. Niu, L. Pollack, M. Rosing, G. 
     +Sandrone,'/
      data chr13/
     +'M. Stave, H. Taylor, G. Thomas, J. H. van Lenthe, A. Wong, Z. Zha
     +ng.'/

c
c     == default ==
      write(luout,*)
      write(luout,*)
      call util_nwchem_paper()
      write(luout,*)
      call util_print_centered(luout, 'AUTHORS & CONTRIBUTORS', 40, 
     & .true.)
      call util_print_centered(luout, chr01, 40, .false.)
      call util_print_centered(luout, chr02, 40, .false.)
      call util_print_centered(luout, chr03, 40, .false.)
      call util_print_centered(luout, chr04, 40, .false.)
      call util_print_centered(luout, chr05, 40, .false.)
      call util_print_centered(luout, chr06, 40, .false.)
      call util_print_centered(luout, chr07, 40, .false.)
      call util_print_centered(luout, chr08, 40, .false.)
      call util_print_centered(luout, chr09, 40, .false.)
      call util_print_centered(luout, chr10, 40, .false.)
      call util_print_centered(luout, chr11, 40, .false.)
      call util_print_centered(luout, chr12, 40, .false.)
      call util_print_centered(luout, chr13, 40, .false.)
c
c     == print high ==
      if (util_print('references',print_high)) then
        write(luout,*)
        call util_print_centered(luout, 'REFERENCES', 40, .true.)
        write(luout,10002)
        write(luout,10101) 1
        write(luout,10102) 2
        write(luout,10103) 3
        write(luout,10104) 4
        write(luout,10105) 5
        write(luout,10106) 6
        write(luout,10107) 7
        write(luout,10108) 8
        write(luout,10109) 9
        write(luout,10110) 10
        write(luout,10111) 11
      endif
c
10002 format(
     &    /,10x,
     &    'Where appropriate, please cite the following reference(s) ',
     &    /,10x,
     &    'when publishing results obtained with NWChem:'
     &    /,/,10x)
c
10101 format(1x,i5,'.',1x,
     &    'D. E. Bernholdt, E. Apra, H. A. Fruchtl, M.F. Guest, ',
     &    'R. J. Harrison, ',/,8x,'R. A. Kendall, R. A. Kutteh, ',
     &    'X. Long, J. B. Nicholas, J. A. Nichols, ',/,8x,
     &    'H. L. Taylor, A. T. Wong, G. I. Fann, R. J. Littlefi',
     &    'eld and',/,8x,'J. Nieplocha,',/,8x,
     &    '"Parallel Computational Chemistry Made Easier: ',/,8x,
     &    'The Development of NWChem,"',/,8x,
     &    'Int. J. Quantum Chem. Symposium  29, 475-483 ',
     &    '(1995).'
     &    /,/)
c
10102 format(1x,i5,'.',1x,
     &    'M. F. Guest, E. Apra, D. E. Bernholdt, H. A. Fruchtl, ',
     &    'R. J. Harrison, ',/,8x,'R. A. Kendall, R. A. Kutteh, ',
     &    'X. Long, J. B. Nicholas, J. A. Nichols, ',/,8x,
     &    'H. L. Taylor, A. T. Wong, G. I. Fann, R. J. Littlefield',
     &    ' and',/,8x,'J. Nieplocha,',/,8x,
     &    '"High Performance Computatio',
     &    'nal Chemistry: NWChem and Fully Distributed',/,8x,
     &    'Parallel Applications", in Advances in Parallel Comput',
     &    'ing, 10, High',/,8x,'Performance Computing: Technology,',
     &    ' Methods, and Applications, Eds. ',/,8x,'J. Dongarra, ',
     &    'L. Gradinetti, G. Joubert, and J. Kowalik, ',/,8x,
     &    '(Elsevier Science B. V. pp. 395-427, 1995).'
     &    /,/)
c
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
c
10104 format(1x,i5,'.',1x,
     &    'M. F. Guest, E. Apra, D. E. Bernholdt, H. A. Fruchtl, ',
     &    'R. J. Harrison,',/,8x,'R. A. Kendall, R. A. Kutteh, ',
     &    'X. Long, J. B. Nicholas, J. A. Nichols, ',/,8x,'H. L. ',
     &    'Taylor, A. T. Wong, G. I. Fann, R. J. Littlefield and',/,8x,
     &    'J. Nieplocha,',/,8x,'"Advances in Parallel Distributed',
     &    ' Data Software; Computational Chemistry',/,8x,
     &    'and NWChem", in Applied Parallel Computing.  Computati',
     &    'ons in Physics,',/,8x,'Chemistry, and Engineering ',
     &    'Science, Eds. J. Wasnieski, J. Dongarra, and ',/,8x,
     &    'K. Madsen, Lecture Notes in Computer Science, 1041,',/,8x,
     &    '(Springer-Verlag, Berlin, 278-294, 1996).'
     &    /,/)
c
10105 format(1x,i5,'.',1x,
     &    'M. F. Guest, E. Apra, D. E. Bernholdt, H. A. Fruchtl, ',
     &    'R. J. Harrison,',/,8x,'R. A. Kendall, R. A. Kutteh, ',
     &    'X. Long, J. B. Nicholas, J. A. Nichols, ',/,8x,
     &    'H. L. Taylor, A. T. Wong, G. I. Fann, R. J. Littlefield',
     &    ' and',/,8x,'J. Nieplocha, ',/,8x,
     &    '"High-Performance Computing in Chemistry; NWChem", ',
     &    'Future ',/,8x,
     &    'Generations Computer Systems  12(4), 273, (1996).'
     &    /,/)
c
10106 format(1x,i5,'.',1x,
     &    'R. A. Kendall, E. Apra, D. E. Bernholdt, E. J. Bylaska, ',
     &    'M. Dupuis,',/,8x,'G. I. Fann, R. J. Harrison, J. Ju, ',
     &    'J. A. Nichols, J. Nieplocha,',/,8x,'T. P. Straatsma, ',
     &    'T. L. Windus, A. T. Wong,',/,8x,'Computer Phys. Comm. 128, ',
     &    '260 (2000).'
     &    /,/)
c
10107 format(1x,i5,'.',1x,
     + 'T. P. Straatsma and V. Helms, in Molecular Dynamics on',
     + ' Parallel',/,8x,'Computers',
     + '2000, R. Esser et al. editors, p 70.',//)
c
10108 format(1x,i5,'.',1x,
     + ' T. P. Straatsma, M. Philippopoulos and J. A. McCammon,',
     + /,8x,' Computer Phys. Comm. 128, ',
     &    '377 (2000).'
     &    /,/)
c
c    ***PAW references***
10109 format(1x,i5,'.',1x,
     + ' M. Valiev, E. J. Bylaska, A. Gramada, and J. H. Weare,'
     + /,8x,' Reviews in Modern Quantum Chemistry, pages 1684-1734',
     + /,8x,' (World Scientific, Singapore, 2002).'
     &    /,/)
c
10110 format(1x,i5,'.',1x,
     + ' E. J. Bylaska, M. Valiev, R. Kawai, and J. H. Weare,',
     + /,8x,' Computer Phys. Comm. 143, ',
     &    '11 (2002).'
     &    /,/)
c
10111 format(1x,i5,'.',1x,
     + ' M. Valiev and J. H. Weare,',
     + /,8x,' J. Phys. Chem. A, 103, ',
     &    '10588 (1999).'
     &    /,/)
c
      end
