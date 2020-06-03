c
      subroutine util_nwchem_paper()
c
      implicit none
      integer maxno_authors
      parameter(maxno_authors=114)
      character*64 name_author(maxno_authors)
      data name_author/
     +     'E. Aprà', 'E. J. Bylaska', 'W. A. de Jong', 'N. Govind',
     +     'K. Kowalski', 'T. P. Straatsma', 'M. Valiev',
     +     'H. J. J. van Dam', 'Y. Alexeev', 'J. Anchell',
     +     'V. Anisimov', 'F. W. Aquino', 'R. Atta-Fynn',
     +     'J. Autschbach', 'N. P. Bauman', 'J. C. Becca',
     +     'D. E. Bernholdt', 'K. Bhaskaran-Nair', 'S. Bogatko',
     +     'P. Borowski', 'J. Boschen', 'J. Brabec', 'A. Bruner',
     +     'E. Cauët', 'Y. Chen', 'G. N. Chuev', 'C. J. Cramer',
     +     'J. Daily', 'M. J. O. Deegan', 'T. H. Dunning Jr.',
     +     'M. Dupuis', 'K. G. Dyall', 'G. I. Fann', 'S. A. Fischer',
     +     'A. Fonari', 'H. Früchtl', 'L. Gagliardi', 'J. Garza',
     +     'N. Gawande', 'S. Ghosh', 'K. Glaesemann', 'A. W. Götz',
     +     'J. Hammond', 'V. Helms', 'E. D. Hermes', 'K. Hirao',
     +     'S. Hirata', 'M. Jacquelin', 'L. Jensen', 'B. G. Johnson',
     +     'H. Jónsson', 'R. A. Kendall', 'M. Klemm', 'R. Kobayashi',
     +     'V. Konkov', 'S. Krishnamoorthy', 'M. Krishnan', 'Z. Lin',
     +     'R. D. Lins', 'R. J. Littlefield', 'A. J. Logsdail',
     +     'K. Lopata', 'W. Ma', 'A. V. Marenich',
     +     'J. Martin del Campo', 'D. Mejia-Rodriguez', 'J. E. Moore',
     +     'J. M. Mullin', 'T. Nakajima', 'D. R. Nascimento',
     +     'J. A. Nichols', 'P. J. Nichols', 'J. Nieplocha',
     +     'A. Otero-de-la-Roza', 'B. Palmer', 'A. Panyala',
     +     'T. Pirojsirikul', 'B. Peng', 'R. Peverati', 'J. Pittner',
     +     'L. Pollack', 'R. M. Richard', 'P. Sadayappan',
     +     'G. C. Schatz', 'W. A. Shelton', 'D. W. Silverstein',
     +     'D. M. A. Smith', 'T. A. Soares', 'D. Song', 'M. Swart',
     +     'H. L. Taylor', 'G. S. Thomas', 'V. Tipparaju',
     +     'D. G. Truhlar', 'K. Tsemekhman', 'T. Van Voorhis',
     +     'Á. Vázquez-Mayagoitia', 'P. Verma', 'O. Villa',
     +     'A. Vishnu', 'K. D. Vogiatzis', 'D. Wang', 'J. H. Weare',
     +     'M. J. Williamson', 'T. L. Windus', 'K. Woliński',
     +     'A. T. Wong', 'Q. Wu', 'C. Yang', 'Q. Yu', 'M. Zacharias',
     +     'Z. Zhang', 'Y. Zhao', 'and R. J. Harrison'/
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
      call util_print_authors(maxno_authors,name_author)
      
      call util_print_centered(LuOut, 
     &     '"NWChem: Past, present, and future',
     &     40, .false.)
      call util_print_centered(LuOut, 
     &     'J. Chem. Phys. 152, 184102 (2020)',
     &     40, .false.)
      call util_print_centered(LuOut, 
     &     'doi:10.1063/5.0004997',
     &     40, .false.)
c
      end 
c $Id$
