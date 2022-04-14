c
      subroutine util_legal()
c
* $Id$
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
#include "inp.fh"
#include "util.fh"
c
      logical  next_line
      external next_line
      integer  ibegin,iend,ilen,i
      character*4096 newchr
      integer maxno_authors,len_newchr
      parameter(maxno_authors=106)
      character*64 name_author(maxno_authors)
      data name_author/
     +     'E. Apra', 'E. J. Bylaska', 
     +     'N. Govind', 'K. Kowalski',
     +     'M. Valiev',
     +     'D. Mejia-Rodriguez', 'A. Kunitsa',
     +     'N. P. Bauman', 'A. Panyala', 
     +     'W. A. de Jong', 'T. P. Straatsma',  'H. J. J. van Dam',
     +     'D. Wang', 'T. L. Windus',
     +     'J. Hammond', 'J. Autschbach', 'K. Bhaskaran-Nair', 
     +     'J. Brabec','K. Lopata','S. A. Fischer','S. Krishnamoorthy',
     +     'M. Jacquelin',
     +     'W. Ma', 'M. Klemm', 'O. Villa', 'Y. Chen', 'V. Anisimov',
     +     'F. Aquino', 'S. Hirata', 'M. T. Hackler',
     +     'Eric Hermes', 'L. Jensen', 'J. E. Moore', 'J. C. Becca',
     +     'V. Konjkov', 'T. Risthaus', 
     +     'M. Malagoli', 'A. Marenich', 'A. Otero-de-la-Roza', 
     +     'J. Mullin', 'P. Nichols', 'R. Peverati', 'J. Pittner', 
     +     'Y. Zhao', 'P.-D. Fan', 'A. Fonari', 'M. J. Williamson',
     +     'R. J. Harrison', 'J. R. Rehr',
     +     'M. Dupuis', 'D. Silverstein', 'D. M. A. Smith', 
     +     'J. Nieplocha', 'V. Tipparaju', 'M. Krishnan', 
     +     'B. E. Van Kuiken', 'A. Vazquez-Mayagoitia', 
     +     'M. Swart', 'Q. Wu', 'T. Van Voorhis', 'A. A. Auer', 
     +     'M. Nooijen', 'L. D. Crosby', 'E. Brown', 
     +     'G. Cisneros', 'G. I. Fann', 'H. Fruchtl', 'J. Garza', 
     +     'K. Hirao', 'R. A. Kendall', 'J. A. Nichols', 
     +     'K. Tsemekhman', 'K. Wolinski', 'J. Anchell', 
     +     'D. E. Bernholdt', 'P. Borowski', 'T. Clark', 
     +     'D. Clerc', 'H. Dachsel', 'M. J. O. Deegan', 
     +     'K. Dyall', 'D. Elwood', 'E. Glendening', 
     +     'M. Gutowski', 'A. C. Hess', 'J. Jaffe', 'B. G. Johnson', 
     +     'J. Ju', 'R. Kobayashi', 'R. Kutteh', 'Z. Lin', 
     +     'R. Littlefield', 'X. Long', 'B. Meng', 
     +     'T. Nakajima', 'S. Niu', 'L. Pollack',
     +     'M. Rosing', 'K. Glaesemann', 'G. Sandrone', 'M. Stave', 
     +     'H. Taylor', 'G. Thomas', 'J. H. van Lenthe', 
     +     'A. T. Wong', 'Z. Zhang.'/
c
c     == default ==
      write(luout,*)
      write(luout,*)
      call util_nwchem_paper()
      write(luout,*)
      call util_print_centered(luout, 'AUTHORS', 40, 
     & .true.)
c
      call util_print_authors(maxno_authors,name_author)

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
c
C> \brief Return the next line from a string of comma separated items
C>
C> The list of authors and contributors is getting rather long and it
C> has become too cumbersome to manually move names around to fit them
C> into the list. This function partitions the comma separated list of
C> names into individual lines. It is indented to be used in the code
C> as
C> \code
C>       ibegin    = 0
C>       mlinelen  = 80
C>       separator = ','
C>       do while (next_line(ibegin,iend,mlinelin,separator,str_names))
C>         write(*,*)str_names(ibegin:iend)
C>       enddo
C> \endcode
C> I.e. the string is partitioned into chunks of at most `mlinelen`
C> characters.
C>
C> \return Returns .true. if another line of names was found, and
C> .false. otherwise.
C>
      logical function next_line(ibegin,iend,mlinelen,separators,
     +                           str_names)
      implicit none
c
#include "inp.fh"
c
      character*(*) str_names  !< [Input] The string of names
      character*(*) separators !< [Input] The string of separator
                               !< characters
      integer mlinelen         !< [Input] The maximum line length
      integer iend             !< [In/Output] On input the last
                               !< character of the previous line,
                               !< on output the end character of the
                               !< next line
      integer ibegin           !< [Output] The beginning of the next
                               !< line
c
c     Local variables
c
      integer inambeg !< The beginning of a name
      integer inamend !< The end of a name
      integer ii, jj  !< Temp variable
      integer ilen
      logical found
c
      ilen = inp_strlen(str_names)
      if (ibegin.eq.0) then
        iend    = ibegin
        inambeg = ibegin
        inamend = ibegin
      else
        ibegin  = iend
        inambeg = ibegin
        inamend = ibegin
      endif
      do while (inp_strtok(str_names,separators,inambeg,inamend))
        if (min(inamend,ilen)-ibegin+2.lt.mlinelen) then
          iend = min(inamend+1,ilen)  ! We need to include the
                                      ! separator char
        else
          exit
        endif
      enddo
c
c     Increase ibegin to skip any separators or spaces
c
      ii = ibegin+1
      do ibegin = ii, iend+1
        found = " ".eq.str_names(ibegin:ibegin)
        do jj = 1, len(separators)
          found = found
     +    .or.(str_names(ibegin:ibegin).eq.separators(jj:jj))
        enddo
        if (.not.found) exit
      enddo
c
      next_line = ibegin.le.iend
c
      return
      end
c
      subroutine util_print_authors(maxno_authors,name_author)
      implicit none
      integer maxno_authors     ! [in]
      character*64 name_author(maxno_authors)
c
      integer ibegin,iend,ilen,i
      integer len_newchr
      character*4096 newchr
#include "stdio.fh"
#include "inp.fh"
#include "util.fh"
      logical  next_line
      external next_line
c      
      ibegin = 1
      do i=1,maxno_authors
         ilen=inp_strlen(name_author(i))
         newchr(ibegin:ibegin+ilen-1)=name_author(i)
         len_newchr=ibegin+ilen-1
         newchr(ibegin+ilen:ibegin+ilen+1)=', '
         ibegin=ibegin+ilen+2
      enddo
      ibegin = 0
      iend   = 0 
      ilen   = 80
      do while (next_line(ibegin,iend,ilen,',',newchr(1:len_newchr)))
        call util_print_centered(luout,newchr(ibegin:iend),40,.false.)
      enddo
      return
      end
