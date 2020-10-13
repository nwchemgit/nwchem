      logical function pre_namiup(lfnout,
     + lseq,cseq,mseq,nseq,latm,catm,matm,natm)
c
c $Id$
c
c     Function to comply with IUPAC atom naming [1-3].
c
c     in  : lfnout    = output file logical file number
c           kse
c           mseq      = dimension of the sequence list
c           nseq      = length of the sequence list
c
c     out : lseq(1,*) = segment numbers
c           lseq(2,*) = number of atoms in segment
c           lseq(3,*) = index to list of unique segments
c           cseq      = segment names
c           nseq      = actual length of the sequence list
c
c     References
c
c     [1] "IUPAC-IUB Commission on Biochemical Nomenclature. 
c         Abbreviations and Symbols for the Description of the
c         Conformation of Polypeptide Chains. Tentative Rules (1969)".
c         Biochemistry (1970) 9, 3471-3479. DOI: 10.1021/bi00820a001
c
c     [2] Charles Hoogstraten, "Correlation of hydrogen atom naming
c         systems, including diastereotopic protons." Web-page:
c         http://www.bmrb.wisc.edu/ref_info/atom_nom.tbl 
c         [accessed Jan 6, 2018]
c
c     [3] IUPAC-IUBMB Joint Commission on Biochemical Nomenclature
c         and Nomenclature Commission of IUBMB, "Biochemical
c         Nomenclature and Related Documents (White Book)." 2nd Edition,
c         Portland Press, 1992, pp. 39-69. ISBN: 1-85578-005-4.
c         Part: "Amino acids, peptides and proteins", Section 3AA-1
c         "Names of Common alpha-Amino Acids". 
c         DOI: 10.1351/pac198456050595
c
      implicit none
c
#include "util.fh"
#include "pre_common.fh"
c
      integer pre_atnum
      external pre_atnum
      logical pre_protein
      external pre_protein
c
      integer lfnout
      integer mseq,matm
      integer nseq,natm
      integer lseq(6,mseq),latm(5,matm)
      character*6 catm(3,matm)
      character*5 cseq(2,mseq)
      character*2 elemnt
      character*3 cres   ! residue name
c
      integer i
      logical flag1
c
      if(util_print('sequence',print_debug)) then
      write(lfnout,2000)
 2000 format(/,'pre_namiup ')
      endif
c
      flag1=.false.
      do 1 i=1,natm
      if(pre_protein(cseq(1,latm(1,i)))) then
      if(catm(1,i)(1:2).eq.'1H') flag1=.true.
      endif
      catm(2,i)=catm(1,i)
    1 continue
c
c     IUPAC hydrogen naming
c
      if(flag1) then
      do 2 i=1,natm
      if(pre_protein(cseq(1,latm(1,i)))) then
      if(catm(2,i)(1:2).eq.'3H') catm(2,i)(1:2)='4H'
      if(catm(2,i)(1:2).eq.'2H') catm(2,i)(1:2)='3H'
      if(catm(2,i)(1:2).eq.'1H') catm(2,i)(1:2)='2H'
      endif
    2 continue
      if(util_print('sequence',print_high)) then
      write(lfnout,2001)
 2001 format(/,'Hydrogen atom prefix corrected ')
      endif
      endif
c
 2002 format('Renamed ',a3,i3.3,x,a)
c
c     ALA  HN  ->  H 
c          HB1 -> 2HB
c          HB2 -> 3HB
c          HB3 -> 4HB
c
      cres='ALA'
      do 3 i=1,natm
      if(cseq(1,latm(1,i))(1:3).eq.'ALA') then
      if(catm(2,i)(1:4).eq.' HN ') then
      catm(2,i)(1:4)=' H  '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HN  to H'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HB1') then
      catm(2,i)(1:4)='2HB '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HB1 to 2HB'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HB2') then
      catm(2,i)(1:4)='3HB '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HB2 to 3HB'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HB3') then
      catm(2,i)(1:4)='4HB '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HB3 to 4HB'
      endif
      endif
      endif
    3 continue
c
c     ARG  HN  ->  H 
c          HB1 -> 2HB
c          HB2 -> 3HB
c          HB3 -> 2HB
c          HG1 -> 2HG
c          HG2 -> 3HG
c          HG3 -> 2HG
c          HD1 -> 2HD
c          HD2 -> 3HD
c          HD3 -> 2HD
c         HH11 -> 2HH1
c         HH12 -> 3HH1
c         HH21 -> 2HH2
c         HH22 -> 3HH2
c
c     Some PDB files list HB1, HB2, some list HB2, HB3, ultimately both
c     cases need to map to 2HB, 3HB. Hence the funky mapping of HB1,
c     HB2, and HB3. Similary for HG1, HG2 and HD1, HD2.
c
      cres='ARG'
      do 4 i=1,natm
      if(cseq(1,latm(1,i))(1:3).eq.'ARG') then
      if(catm(2,i)(1:4).eq.' HN ') then
      catm(2,i)(1:4)=' H  '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HN  to H'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HB1') then
      catm(2,i)(1:4)='2HB '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HB1 to 2HB'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HB2') then
      catm(2,i)(1:4)='3HB '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HB2 to 3HB'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HB3') then
      catm(2,i)(1:4)='2HB '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HB3 to 2HB'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HG1') then
      catm(2,i)(1:4)='2HG '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HG1 to 2HG'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HG2') then
      catm(2,i)(1:4)='3HG '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HG2 to 3HG'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HG3') then
      catm(2,i)(1:4)='2HG '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HG3 to 2HG'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HD1') then
      catm(2,i)(1:4)='2HD '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HD1 to 2HD'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HD2') then
      catm(2,i)(1:4)='3HD '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HD2 to 3HD'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HD3') then
      catm(2,i)(1:4)='2HD '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HD3 to 2HD'
      endif
      endif
      if(catm(2,i)(1:4).eq.'HH11') then
      catm(2,i)(1:4)='2HH1'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HH11 to 2HH1'
      endif
      endif
      if(catm(2,i)(1:4).eq.'HH12') then
      catm(2,i)(1:4)='3HH1'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HH12 to 3HH1'
      endif
      endif
      if(catm(2,i)(1:4).eq.'HH21') then
      catm(2,i)(1:4)='2HH2'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HH21 to 2HH2'
      endif
      endif
      if(catm(2,i)(1:4).eq.'HH22') then
      catm(2,i)(1:4)='3HH2'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HH22 to 3HH2'
      endif
      endif
      endif
    4 continue
c
c     ASN  HN  ->  H 
c          HB1 -> 2HB
c          HB2 -> 3HB
c          HB3 -> 2HB
c         HD21 -> 2HD2
c         HD22 -> 3HD2
c
c     Some PDB files list HB1, HB2, some list HB2, HB3, ultimately both
c     cases need to map to 2HB, 3HB. Hence the funky mapping of HB1,
c     HB2, and HB3.
c
      cres='ASN'
      do 5 i=1,natm
      if(cseq(1,latm(1,i))(1:3).eq.'ASN') then
      if(catm(2,i)(1:4).eq.' HN ') then
      catm(2,i)(1:4)=' H  '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HN  to H'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HB1') then
      catm(2,i)(1:4)='2HB '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HB1 to 2HB'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HB2') then
      catm(2,i)(1:4)='3HB '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HB2 to 3HB'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HB3') then
      catm(2,i)(1:4)='2HB '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HB3 to 2HB'
      endif
      endif
      if(catm(2,i)(1:4).eq.'HD21') then
      catm(2,i)(1:4)='2HD2'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HD21 to 2HD2'
      endif
      endif
      if(catm(2,i)(1:4).eq.'HD22') then
      catm(2,i)(1:4)='3HD2'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HD22 to 3HD2'
      endif
      endif
      endif
    5 continue
c
c     ASP  HN  ->  H 
c          HB1 -> 2HB
c          HB2 -> 3HB
c          HB3 -> 2HB
c
c     Some PDB files list HB1, HB2, some list HB2, HB3, ultimately both
c     cases need to map to 2HB, 3HB. Hence the funky mapping of HB1,
c     HB2, and HB3.
c
      cres='ASP'
      do 6 i=1,natm
      if(cseq(1,latm(1,i))(1:3).eq.'ASP') then
      if(catm(2,i)(1:4).eq.' HN ') then
      catm(2,i)(1:4)=' H  '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HN  to H'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HB1') then
      catm(2,i)(1:4)='2HB '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HB1 to 2HB'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HB2') then
      catm(2,i)(1:4)='3HB '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HB2 to 3HB'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HB3') then
      catm(2,i)(1:4)='2HB '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HB3 to 2HB'
      endif
      endif
      endif
    6 continue
c
c     CYS  HN  ->  H 
c          HB1 -> 2HB
c          HB2 -> 3HB
c          HB3 -> 2HB
c
c     Some PDB files list HB1, HB2, some list HB2, HB3, ultimately both
c     cases need to map to 2HB, 3HB. Hence the funky mapping of HB1,
c     HB2, and HB3.
c
      cres='CYS'
      do 7 i=1,natm
      if(cseq(1,latm(1,i))(1:3).eq.'CYS'.or.
     +   cseq(1,latm(1,i))(1:3).eq.'CYX') then
      if(catm(2,i)(1:4).eq.' HN ') then
      catm(2,i)(1:4)=' H  '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HN  to H'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HB1') then
      catm(2,i)(1:4)='2HB '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HB1 to 2HB'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HB2') then
      catm(2,i)(1:4)='3HB '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HB2 to 3HB'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HB3') then
      catm(2,i)(1:4)='2HB '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HB3 to 2HB'
      endif
      endif
      endif
    7 continue
c
c     GLN  HN  ->  H 
c          HB1 -> 2HB
c          HB2 -> 3HB
c          HB3 -> 2HB
c          HG1 -> 2HG
c          HG2 -> 3HG
c          HG3 -> 2HG
c         HE21 -> 2HE2
c         HE22 -> 3HE2
c
c     Some PDB files list HB1, HB2, some list HB2, HB3, ultimately both
c     cases need to map to 2HB, 3HB. Hence the funky mapping of HB1,
c     HB2, and HB3. Likewise for the HG atoms.
c
      cres='GLN'
      do 8 i=1,natm
      if(cseq(1,latm(1,i))(1:3).eq.'GLN') then
      if(catm(2,i)(1:4).eq.' HN ') then
      catm(2,i)(1:4)=' H  '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HN  to H'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HB1') then
      catm(2,i)(1:4)='2HB '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HB1 to 2HB'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HB2') then
      catm(2,i)(1:4)='3HB '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HB2 to 3HB'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HB3') then
      catm(2,i)(1:4)='2HB '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HB3 to 2HB'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HG1') then
      catm(2,i)(1:4)='2HG '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HG1 to 2HG'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HG2') then
      catm(2,i)(1:4)='3HG '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HG2 to 3HG'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HG3') then
      catm(2,i)(1:4)='2HG '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HG3 to 2HG'
      endif
      endif
      if(catm(2,i)(1:4).eq.'HE21') then
      catm(2,i)(1:4)='2HE2'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HE21 to 2HE2'
      endif
      endif
      if(catm(2,i)(1:4).eq.'HE22') then
      catm(2,i)(1:4)='3HE2'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HE22 to 3HE2'
      endif
      endif
      endif
    8 continue
c
c     GLU  HN  ->  H 
c          HB1 -> 2HB
c          HB2 -> 3HB
c          HB3 -> 2HB
c          HG1 -> 2HG
c          HG2 -> 3HG
c          HG3 -> 2HG
c
c     Some PDB files list HB1, HB2, some list HB2, HB3, ultimately both
c     cases need to map to 2HB, 3HB. Hence the funky mapping of HB1,
c     HB2, and HB3. Likewise for the HG atoms.
c
      cres='GLU'
      do 9 i=1,natm
      if(cseq(1,latm(1,i))(1:3).eq.'GLU') then
      if(catm(2,i)(1:4).eq.' HN ') then
      catm(2,i)(1:4)=' H  '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HN  to H'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HB1') then
      catm(2,i)(1:4)='2HB '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HB1 to 2HB'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HB2') then
      catm(2,i)(1:4)='3HB '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HB2 to 3HB'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HB3') then
      catm(2,i)(1:4)='2HB '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HB3 to 2HB'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HG1') then
      catm(2,i)(1:4)='2HG '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HG1 to 2HG'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HG2') then
      catm(2,i)(1:4)='3HG '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HG2 to 3HG'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HG3') then
      catm(2,i)(1:4)='2HG '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HG3 to 2HG'
      endif
      endif
      endif
    9 continue
c
c     GLY  HN  ->  H 
c          HA1 -> 2HA
c          HA2 -> 3HA
c          HA3 -> 2HA
c
c     Some PDB files list HA1, HA2, some list HA2, HA3, ultimately both
c     cases need to map to 2HA, 3HA. Hence the funky mapping of HA1,
c     HA2, and HA3.
c
      cres='GLY'
      do 10 i=1,natm
      if(cseq(1,latm(1,i))(1:3).eq.'GLY') then
      if(catm(2,i)(1:4).eq.' HN ') then
      catm(2,i)(1:4)=' H  '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HN  to H'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HA1') then
      catm(2,i)(1:4)='2HA '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HA1 to 2HA'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HA2') then
      catm(2,i)(1:4)='3HA '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HA2 to 3HA'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HA3') then
      catm(2,i)(1:4)='2HA '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HA3 to 2HA'
      endif
      endif
      endif
   10 continue
c
c     HIS  HN  ->  H 
c          HB1 -> 2HB
c          HB2 -> 3HB
c          HB3 -> 2HB
c
c     Some PDB files list HB1, HB2, some list HB2, HB3, ultimately both
c     cases need to map to 2HB, 3HB. Hence the funky mapping of HB1,
c     HB2, and HB3.
c
      cres='HIS'
      do 11 i=1,natm
      if(cseq(1,latm(1,i))(1:3).eq.'HIS'.or.
     +   cseq(1,latm(1,i))(1:3).eq.'HID'.or. ! amber
     +   cseq(1,latm(1,i))(1:3).eq.'HIE'.or.
     +   cseq(1,latm(1,i))(1:3).eq.'HIP'.or.
     +   cseq(1,latm(1,i))(1:3).eq.'HSD'.or. ! charmm
     +   cseq(1,latm(1,i))(1:3).eq.'HSE'.or.
     +   cseq(1,latm(1,i))(1:3).eq.'HSP') then
      if(catm(2,i)(1:4).eq.' HN ') then
      catm(2,i)(1:4)=' H  '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HN  to H'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HB1') then
      catm(2,i)(1:4)='2HB '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HB1 to 2HB'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HB2') then
      catm(2,i)(1:4)='3HB '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HB2 to 3HB'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HB3') then
      catm(2,i)(1:4)='2HB '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HB3 to 2HB'
      endif
      endif
      endif
   11 continue
c
c     ILE  HN  ->  H 
c         HG21 -> 2HG2
c         HG22 -> 3HG2
c         HG23 -> 4HG2
c         HG11 -> 2HG1
c         HG12 -> 3HG1
c         HG13 -> 2HG1
c          CD1 ->  CD
c          HD1 -> 2HD
c          HD2 -> 3HD
c          HD3 -> 4HD
c         HD11 -> 2HD
c         HD12 -> 3HD
c         HD13 -> 4HD
c
c     Some PDB files list HG11, HG12, some list HG12, HG13, ultimately
c     both cases need to map to 2HG1, 3HG1. Hence the funky mapping of 
c     HG11, HG12, and HG13.
c     Also some PDB files list HD1, HD2, and HD3, some instead list
c     HD11, HD12, and HD13.
c
      cres='ILE'
      do 12 i=1,natm
      if(cseq(1,latm(1,i))(1:3).eq.'ILE') then
      if(catm(2,i)(1:4).eq.' HN ') then
      catm(2,i)(1:4)=' H  '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HN  to H'
      endif
      endif
      if(catm(2,i)(1:4).eq.'HG21') then
      catm(2,i)(1:4)='2HG2'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HG21 to 2HG2'
      endif
      endif
      if(catm(2,i)(1:4).eq.'HG22') then
      catm(2,i)(1:4)='3HG2'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HG22 to 3HG2'
      endif
      endif
      if(catm(2,i)(1:4).eq.'HG23') then
      catm(2,i)(1:4)='4HG2'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HG22 to 4HG2'
      endif
      endif
      if(catm(2,i)(1:4).eq.'HG11') then
      catm(2,i)(1:4)='2HG1'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HG11 to 2HG1'
      endif
      endif
      if(catm(2,i)(1:4).eq.'HG12') then
      catm(2,i)(1:4)='3HG1'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HG12 to 3HG1'
      endif
      endif
      if(catm(2,i)(1:4).eq.'HG13') then
      catm(2,i)(1:4)='2HG1'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HG13 to 2HG1'
      endif
      endif
      if(catm(2,i)(1:4).eq.' CD1') then
      catm(2,i)(1:4)=' CD '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'CD1 to CD'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HD1') then
      catm(2,i)(1:4)='2HD '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HD1 to 2HD'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HD2') then
      catm(2,i)(1:4)='3HD '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HD2 to 3HD'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HD3') then
      catm(2,i)(1:4)='4HD '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HD3 to 4HD'
      endif
      endif
      if(catm(2,i)(1:4).eq.'HD11') then
      catm(2,i)(1:4)='2HD '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HD11 to 2HD'
      endif
      endif
      if(catm(2,i)(1:4).eq.'HD12') then
      catm(2,i)(1:4)='3HD '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HD12 to 3HD'
      endif
      endif
      if(catm(2,i)(1:4).eq.'HD13') then
      catm(2,i)(1:4)='4HD '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HD13 to 4HD'
      endif
      endif
      endif
   12 continue
c
c     LEU  HN  ->  H 
c          HB1 -> 2HB
c          HB2 -> 3HB
c          HB3 -> 2HB
c         HD11 -> 2HD1
c         HD12 -> 3HD1
c         HD13 -> 4HD1
c         HD21 -> 2HD2
c         HD22 -> 3HD2
c         HD23 -> 4HD2
c
c     Some PDB files list HB1, HB2, some list HB2, HB3, ultimately both
c     cases need to map to 2HB, 3HB. Hence the funky mapping of HB1,
c     HB2, and HB3.
c
      cres='LEU'
      do 13 i=1,natm
      if(cseq(1,latm(1,i))(1:3).eq.'LEU') then
      if(catm(2,i)(1:4).eq.' HN ') then
      catm(2,i)(1:4)=' H  '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HN  to H'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HB1') then
      catm(2,i)(1:4)='2HB '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HB1 to 2HB'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HB2') then
      catm(2,i)(1:4)='3HB '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HB2 to 3HB'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HB3') then
      catm(2,i)(1:4)='2HB '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HB3 to 2HB'
      endif
      endif
      if(catm(2,i)(1:4).eq.'HD11') then
      catm(2,i)(1:4)='2HD1'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HD11 to 2HD1'
      endif
      endif
      if(catm(2,i)(1:4).eq.'HD12') then
      catm(2,i)(1:4)='3HD1'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HD12 to 3HD1'
      endif
      endif
      if(catm(2,i)(1:4).eq.'HD13') then
      catm(2,i)(1:4)='4HD1'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HD13 to 4HD1'
      endif
      endif
      if(catm(2,i)(1:4).eq.'HD21') then
      catm(2,i)(1:4)='2HD2'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HD21 to 2HD2'
      endif
      endif
      if(catm(2,i)(1:4).eq.'HD22') then
      catm(2,i)(1:4)='3HD2'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HD22 to 3HD2'
      endif
      endif
      if(catm(2,i)(1:4).eq.'HD23') then
      catm(2,i)(1:4)='4HD2'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HD23 to 4HD2'
      endif
      endif
      endif
   13 continue
c
c     LYS  HN  ->  H 
c          HB1 -> 2HB
c          HB2 -> 3HB
c          HB3 -> 2HB
c          HG1 -> 2HG
c          HG2 -> 3HG
c          HG3 -> 2HG
c          HD1 -> 2HD
c          HD2 -> 3HD
c          HD3 -> 2HD
c          HE1 -> 2HE
c          HE2 -> 3HE
c          HE3 -> 2HE
c          HZ1 -> 2HZ
c          HZ2 -> 3HZ
c          HZ3 -> 4HZ
c
c     Some PDB files list HB1, HB2, some list HB2, HB3, ultimately both
c     cases need to map to 2HB, 3HB. Hence the funky mapping of HB1,
c     HB2, and HB3. Similarly for HG, HD, and HE atoms.
c
      cres='LYS'
      do 14 i=1,natm
      if(cseq(1,latm(1,i))(1:3).eq.'LYS') then
      if(catm(2,i)(1:4).eq.' HN ') then
      catm(2,i)(1:4)=' H  '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HN  to H'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HB1') then
      catm(2,i)(1:4)='2HB '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HB1 to 2HB'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HB2') then
      catm(2,i)(1:4)='3HB '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HB2 to 3HB'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HB3') then
      catm(2,i)(1:4)='2HB '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HB3 to 2HB'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HG1') then
      catm(2,i)(1:4)='2HG '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HG1 to 2HG'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HG2') then
      catm(2,i)(1:4)='3HG '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HG2 to 3HG'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HG3') then
      catm(2,i)(1:4)='2HG '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HG3 to 2HG'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HD1') then
      catm(2,i)(1:4)='2HD '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HD1 to 2HD'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HD2') then
      catm(2,i)(1:4)='3HD '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HD2 to 3HD'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HD3') then
      catm(2,i)(1:4)='2HD '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HD3 to 2HD'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HE1') then
      catm(2,i)(1:4)='2HE '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HE1 to 2HE'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HE2') then
      catm(2,i)(1:4)='3HE '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HE2 to 3HE'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HE3') then
      catm(2,i)(1:4)='2HE '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HE3 to 2HE'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HZ1') then
      catm(2,i)(1:4)='2HZ '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HZ1 to 2HZ'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HZ2') then
      catm(2,i)(1:4)='3HZ '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HZ2 to 3HZ'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HZ3') then
      catm(2,i)(1:4)='4HZ '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HZ3 to 4HZ'
      endif
      endif
      endif
   14 continue
c
c     MET  HN  ->  H 
c          HB1 -> 2HB
c          HB2 -> 3HB
c          HB3 -> 2HB
c          HG1 -> 2HG
c          HG2 -> 3HG
c          HG3 -> 2HG
c          HE1 -> 2HE
c          HE2 -> 3HE
c          HE3 -> 4HE
c
c     Some PDB files list HB1, HB2, some list HB2, HB3, ultimately both
c     cases need to map to 2HB, 3HB. Hence the funky mapping of HB1,
c     HB2, and HB3. Similarly HG1, HG2.
c
      cres='MET'
      do 15 i=1,natm
      if(cseq(1,latm(1,i))(1:3).eq.'MET') then
      if(catm(2,i)(1:4).eq.' HN ') then
      catm(2,i)(1:4)=' H  '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HN  to H'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HB1') then
      catm(2,i)(1:4)='2HB '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HB1 to 2HB'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HB2') then
      catm(2,i)(1:4)='3HB '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HB2 to 3HB'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HB3') then
      catm(2,i)(1:4)='2HB '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HB3 to 2HB'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HG1') then
      catm(2,i)(1:4)='2HG '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HG1 to 2HG'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HG2') then
      catm(2,i)(1:4)='3HG '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HG2 to 3HG'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HG3') then
      catm(2,i)(1:4)='2HG '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HG3 to 2HG'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HE1') then
      catm(2,i)(1:4)='2HE '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HE1 to 2HE'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HE2') then
      catm(2,i)(1:4)='3HE '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HE2 to 3HE'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HE3') then
      catm(2,i)(1:4)='4HE '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HE3 to 4HE'
      endif
      endif
      endif
   15 continue
c
c     PHE  HN  ->  H 
c          HB1 -> 2HB
c          HB2 -> 3HB
c          HB3 -> 2HB
c
c     Some PDB files list HB1, HB2, some list HB2, HB3, ultimately both
c     cases need to map to 2HB, 3HB. Hence the funky mapping of HB1,
c     HB2, and HB3.
c
      cres='PHE'
      do 16 i=1,natm
      if(cseq(1,latm(1,i))(1:3).eq.'PHE') then
      if(catm(2,i)(1:4).eq.' HN ') then
      catm(2,i)(1:4)=' H  '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HN  to H'
      endif
      endif
      if(catm(2,i)(1:4).eq.' H1 ') then
      catm(2,i)(1:4)='2H  '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HN  to H'
      endif
      endif
      if(catm(2,i)(1:4).eq.' H2 ') then
      catm(2,i)(1:4)='3H  '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HN  to H'
      endif
      endif
      if(catm(2,i)(1:4).eq.' H3 ') then
      catm(2,i)(1:4)='4H  '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HN  to H'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HB1') then
      catm(2,i)(1:4)='2HB '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HB1 to 2HB'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HB2') then
      catm(2,i)(1:4)='3HB '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HB2 to 3HB'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HB3') then
      catm(2,i)(1:4)='2HB '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HB3 to 2HB'
      endif
      endif
      endif
   16 continue
c
c     PRO  HB1 -> 2HB
c          HB2 -> 3HB
c          HB3 -> 2HB
c          HG1 -> 2HG
c          HG2 -> 3HG
c          HG3 -> 2HG
c          HD1 -> 2HD
c          HD2 -> 3HD
c          HD3 -> 2HD
c
c     Some PDB files list HB1, HB2, some list HB2, HB3, ultimately both
c     cases need to map to 2HB, 3HB. Hence the funky mapping of HB1,
c     HB2, and HB3.
c
      cres='PRO'
      do 17 i=1,natm
      if(cseq(1,latm(1,i))(1:3).eq.'PRO') then
      if(catm(2,i)(1:4).eq.' HB1') then
      catm(2,i)(1:4)='2HB '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HB1 to 2HB'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HB2') then
      catm(2,i)(1:4)='3HB '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HB2 to 3HB'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HB3') then
      catm(2,i)(1:4)='2HB '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HB3 to 2HB'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HG1') then
      catm(2,i)(1:4)='2HG '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HG1 to 2HG'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HG2') then
      catm(2,i)(1:4)='3HG '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HG2 to 3HG'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HG3') then
      catm(2,i)(1:4)='2HG '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HG3 to 2HG'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HD1') then
      catm(2,i)(1:4)='2HD '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HD1 to 2HD'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HD2') then
      catm(2,i)(1:4)='3HD '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HD2 to 3HD'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HD3') then
      catm(2,i)(1:4)='2HD '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HD3 to 2HD'
      endif
      endif
      endif
   17 continue
c
c     SER  HN  ->  H 
c          H1  -> 2H
c          H2  -> 3H
c          H3  -> 4H
c          HB1 -> 2HB
c          HB2 -> 3HB
c          HB3 -> 2HB
c          HG1 ->  HG
c
c     Some PDB files list HB1, HB2, some list HB2, HB3, ultimately both
c     cases need to map to 2HB, 3HB. Hence the funky mapping of HB1,
c     HB2, and HB3.
c
      cres='SER'
      do 18 i=1,natm
      if(cseq(1,latm(1,i))(1:3).eq.'SER') then
      if(catm(2,i)(1:4).eq.' HN ') then
      catm(2,i)(1:4)=' H  '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HN  to H'
      endif
      endif
      if(catm(2,i)(1:4).eq.' H1 ') then
      catm(2,i)(1:4)='2H  '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'H1  to 2H'
      endif
      endif
      if(catm(2,i)(1:4).eq.' H2 ') then
      catm(2,i)(1:4)='3H  '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'H2  to 3H'
      endif
      endif
      if(catm(2,i)(1:4).eq.' H3 ') then
      catm(2,i)(1:4)='4H  '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'H3  to 4H'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HB1') then
      catm(2,i)(1:4)='2HB '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HB1 to 2HB'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HB2') then
      catm(2,i)(1:4)='3HB '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HB2 to 3HB'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HB3') then
      catm(2,i)(1:4)='2HB '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HB2 to 3HB'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HG1') then
      catm(2,i)(1:4)=' HG '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HG1 to HG'
      endif
      endif
      endif
   18 continue
c
c     THR  HN  ->  H 
c         HG21 -> 2HG2
c         HG22 -> 3HG2
c         HG23 -> 4HG2
c
      cres='THR'
      do 19 i=1,natm
      if(cseq(1,latm(1,i))(1:3).eq.'THR') then
      if(catm(2,i)(1:4).eq.' HN ') then
      catm(2,i)(1:4)=' H  '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HN  to H'
      endif
      endif
      if(catm(2,i)(1:4).eq.'HG21') then
      catm(2,i)(1:4)='2HG2'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HG21 to 2HG2'
      endif
      endif
      if(catm(2,i)(1:4).eq.'HG22') then
      catm(2,i)(1:4)='3HG2'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HG22 to 3HG2'
      endif
      endif
      if(catm(2,i)(1:4).eq.'HG23') then
      catm(2,i)(1:4)='4HG2'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HG22 to 4HG2'
      endif
      endif
      endif
   19 continue
c
c     TRP  HN  ->  H 
c          HB1 -> 2HB
c          HB2 -> 3HB
c          HB3 -> 2HB
c
c     Some PDB files list HB1, HB2, some list HB2, HB3, ultimately both
c     cases need to map to 2HB, 3HB. Hence the funky mapping of HB1,
c     HB2, and HB3.
c
      cres='TRP'
      do 20 i=1,natm
      if(cseq(1,latm(1,i))(1:3).eq.'TRP') then
      if(catm(2,i)(1:4).eq.' HN ') then
      catm(2,i)(1:4)=' H  '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HN  to H'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HB1') then
      catm(2,i)(1:4)='2HB '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HB1 to 2HB'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HB2') then
      catm(2,i)(1:4)='3HB '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HB2 to 3HB'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HB3') then
      catm(2,i)(1:4)='2HB '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HB3 to 2HB'
      endif
      endif
      endif
   20 continue
c
c     TYR  HN  ->  H 
c          HB1 -> 2HB
c          HB2 -> 3HB
c          HB3 -> 2HB
c
c     Some PDB files list HB1, HB2, some list HB2, HB3, ultimately both
c     cases need to map to 2HB, 3HB. Hence the funky mapping of HB1,
c     HB2, and HB3.
c
      cres='TYR'
      do 21 i=1,natm
      if(cseq(1,latm(1,i))(1:3).eq.'TYR') then
      if(catm(2,i)(1:4).eq.' HN ') then
      catm(2,i)(1:4)=' H  '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HN  to H'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HB1') then
      catm(2,i)(1:4)='2HB '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HB1 to 2HB'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HB2') then
      catm(2,i)(1:4)='3HB '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HB2 to 3HB'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HB3') then
      catm(2,i)(1:4)='2HB '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HB3 to 2HB'
      endif
      endif
      endif
   21 continue
c
c     VAL  HN  ->  H 
c          H1  -> 2H
c          H2  -> 3H
c          H3  -> 4H
c         HG11 -> 2HG1
c         HG12 -> 3HG1
c         HG13 -> 4HG1
c         HG21 -> 2HG2
c         HG22 -> 3HG2
c         HG23 -> 4HG2
c
      cres='VAL'
      do 22 i=1,natm
      if(cseq(1,latm(1,i))(1:3).eq.'VAL') then
      if(catm(2,i)(1:4).eq.' HN ') then
      catm(2,i)(1:4)=' H  '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HN  to H'
      endif
      endif
      if(catm(2,i)(1:4).eq.' H1 ') then
      catm(2,i)(1:4)='2H  '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'H1  to 2H'
      endif
      endif
      if(catm(2,i)(1:4).eq.' H2 ') then
      catm(2,i)(1:4)='3H  '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'H2  to 3H'
      endif
      endif
      if(catm(2,i)(1:4).eq.' H3 ') then
      catm(2,i)(1:4)='4H  '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'H3  to 4H'
      endif
      endif
      if(catm(2,i)(1:4).eq.'HG11') then
      catm(2,i)(1:4)='2HG1'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HG11 to 2HG1'
      endif
      endif
      if(catm(2,i)(1:4).eq.'HG12') then
      catm(2,i)(1:4)='3HG1'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HG12 to 3HG1'
      endif
      endif
      if(catm(2,i)(1:4).eq.'HG13') then
      catm(2,i)(1:4)='4HG1'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HG13 to 4HG1'
      endif
      endif
      if(catm(2,i)(1:4).eq.'HG21') then
      catm(2,i)(1:4)='2HG2'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HG21 to 2HG2'
      endif
      endif
      if(catm(2,i)(1:4).eq.'HG22') then
      catm(2,i)(1:4)='3HG2'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HG22 to 3HG2'
      endif
      endif
      if(catm(2,i)(1:4).eq.'HG23') then
      catm(2,i)(1:4)='4HG2'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HG22 to 4HG2'
      endif
      endif
      endif
   22 continue
c
c     DA   O3' ->  O3*
c          O5' ->  O5*
c          C5' ->  C5*
c          H5' -> 2H5*
c         H5'' -> 3H5*
c          C4' ->  C4*
c          H4' ->  H4*
c          O4' ->  O4*
c          C1' ->  C1*
c          H1' ->  H1*
c          C3' ->  C3*
c          H3' ->  H3*
c          C2' ->  C2*
c          H2' -> 2H2*
c         H2'' -> 3H2*
c          OP1 ->  O1P
c          OP2 ->  O2P
c          H61 -> 2H6
c          H62 -> 3H6
c         HO3' ->  H3T
c         HO5' ->  H5T
c
      cres='DA'
      do 23 i=1,natm
      if(cseq(1,latm(1,i))(1:3).eq.'DA '.or.
     +   cseq(1,latm(1,i))(1:3).eq.'DA3'.or.
     +   cseq(1,latm(1,i))(1:3).eq.'DA5') then
      if(catm(2,i)(1:4).eq." O3'") then
      catm(2,i)(1:4)=' O3*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"O3' to O3*"
      endif
      endif
      if(catm(2,i)(1:4).eq." O5'") then
      catm(2,i)(1:4)=' O5*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"O5' to O5*"
      endif
      endif
      if(catm(2,i)(1:4).eq." C5'") then
      catm(2,i)(1:4)=' C5*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"C5' to C5*"
      endif
      endif
      if(catm(2,i)(1:4).eq." H5'") then
      catm(2,i)(1:4)='2H5*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"H5' to 2H5*"
      endif
      endif
      if(catm(2,i)(1:4).eq."H5''") then
      catm(2,i)(1:4)='3H5*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"H5'' to 3H5*"
      endif
      endif
      if(catm(2,i)(1:4).eq." C4'") then
      catm(2,i)(1:4)=' C4*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"C4' to C4*"
      endif
      endif
      if(catm(2,i)(1:4).eq." H4'") then
      catm(2,i)(1:4)=' H4*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"H4' to H4*"
      endif
      endif
      if(catm(2,i)(1:4).eq." O4'") then
      catm(2,i)(1:4)=' O4*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"O4' to O4*"
      endif
      endif
      if(catm(2,i)(1:4).eq." C1'") then
      catm(2,i)(1:4)=' C1*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"C1' to C1*"
      endif
      endif
      if(catm(2,i)(1:4).eq." H1'") then
      catm(2,i)(1:4)=' H1*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"H1' to H1*"
      endif
      endif
      if(catm(2,i)(1:4).eq." C3'") then
      catm(2,i)(1:4)=' C3*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"C3' to C3*"
      endif
      endif
      if(catm(2,i)(1:4).eq." H3'") then
      catm(2,i)(1:4)=' H3*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"H3' to H3*"
      endif
      endif
      if(catm(2,i)(1:4).eq." C2'") then
      catm(2,i)(1:4)=' C2*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"C2' to C2*"
      endif
      endif
      if(catm(2,i)(1:4).eq." H2'") then
      catm(2,i)(1:4)='2H2*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"H2' to 2H2*"
      endif
      endif
      if(catm(2,i)(1:4).eq."H2''") then
      catm(2,i)(1:4)='3H2*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"H2'' to 3H2*"
      endif
      endif
      if(catm(2,i)(1:4).eq." OP1") then
      catm(2,i)(1:4)=' O1P'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"OP1 to O1P"
      endif
      endif
      if(catm(2,i)(1:4).eq." OP2") then
      catm(2,i)(1:4)=' O2P'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"OP2 to O2P"
      endif
      endif
      if(catm(2,i)(1:4).eq." H61") then
      catm(2,i)(1:4)='2H6 '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"H61 to 2H6"
      endif
      endif
      if(catm(2,i)(1:4).eq." H62") then
      catm(2,i)(1:4)='3H6 '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"H61 to 2H6"
      endif
      endif
      if(catm(2,i)(1:4).eq."HO3'") then
      catm(2,i)(1:4)=' H3T'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"HO3' to  H3T"
      endif
      endif
      if(catm(2,i)(1:4).eq."HO5'") then
      catm(2,i)(1:4)=' H5T'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"HO5' to  H5T"
      endif
      endif
      endif
   23 continue
c
c     DT   O3' ->  O3*
c          O5' ->  O5*
c          C5' ->  C5*
c          H5' -> 2H5*
c         H5'' -> 3H5*
c          C7  ->  C5M
c          H71 -> 2H5M
c          H72 -> 3H5M
c          H73 -> 4H5M
c          C4' ->  C4*
c          H4' ->  H4*
c          O4' ->  O4*
c          C1' ->  C1*
c          H1' ->  H1*
c          C3' ->  C3*
c          H3' ->  H3*
c          C2' ->  C2*
c          H2' -> 2H2*
c         H2'' -> 3H2*
c          OP1 ->  O1P
c          OP2 ->  O2P
c         HO3' ->  H3T
c         HO5' ->  H5T
c
      cres='DT'
      do 24 i=1,natm
      if(cseq(1,latm(1,i))(1:3).eq.'DT '.or.
     +   cseq(1,latm(1,i))(1:3).eq.'DT3'.or.
     +   cseq(1,latm(1,i))(1:3).eq.'DT5') then
      if(catm(2,i)(1:4).eq." O3'") then
      catm(2,i)(1:4)=' O3*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"O3' to O3*"
      endif
      endif
      if(catm(2,i)(1:4).eq." O5'") then
      catm(2,i)(1:4)=' O5*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"O5' to O5*"
      endif
      endif
      if(catm(2,i)(1:4).eq." C5'") then
      catm(2,i)(1:4)=' C5*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"C5' to C5*"
      endif
      endif
      if(catm(2,i)(1:4).eq." H5'") then
      catm(2,i)(1:4)='2H5*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"H5' to 2H5*"
      endif
      endif
      if(catm(2,i)(1:4).eq."H5''") then
      catm(2,i)(1:4)='3H5*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"H5'' to 3H5*"
      endif
      endif
      if(catm(2,i)(1:4).eq." C7 ") then
      catm(2,i)(1:4)=' C5M'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"C7  to C5M"
      endif
      endif
      if(catm(2,i)(1:4).eq." H71") then
      catm(2,i)(1:4)='2H5M'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"H71 to 2C5M"
      endif
      endif
      if(catm(2,i)(1:4).eq." H72") then
      catm(2,i)(1:4)='3H5M'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"H72 to 3C5M"
      endif
      endif
      if(catm(2,i)(1:4).eq." H73") then
      catm(2,i)(1:4)='4H5M'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"H73 to 4C5M"
      endif
      endif
      if(catm(2,i)(1:4).eq." C4'") then
      catm(2,i)(1:4)=' C4*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"C4' to C4*"
      endif
      endif
      if(catm(2,i)(1:4).eq." H4'") then
      catm(2,i)(1:4)=' H4*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"H4' to H4*"
      endif
      endif
      if(catm(2,i)(1:4).eq." O4'") then
      catm(2,i)(1:4)=' O4*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"O4' to O4*"
      endif
      endif
      if(catm(2,i)(1:4).eq." C1'") then
      catm(2,i)(1:4)=' C1*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"C1' to C1*"
      endif
      endif
      if(catm(2,i)(1:4).eq." H1'") then
      catm(2,i)(1:4)=' H1*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"H1' to H1*"
      endif
      endif
      if(catm(2,i)(1:4).eq." C3'") then
      catm(2,i)(1:4)=' C3*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"C3' to C3*"
      endif
      endif
      if(catm(2,i)(1:4).eq." H3'") then
      catm(2,i)(1:4)=' H3*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"H3' to H3*"
      endif
      endif
      if(catm(2,i)(1:4).eq." C2'") then
      catm(2,i)(1:4)=' C2*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"C2' to C2*"
      endif
      endif
      if(catm(2,i)(1:4).eq." H2'") then
      catm(2,i)(1:4)='2H2*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"H2' to 2H2*"
      endif
      endif
      if(catm(2,i)(1:4).eq."H2''") then
      catm(2,i)(1:4)='3H2*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"H2'' to 3H2*"
      endif
      endif
      if(catm(2,i)(1:4).eq." OP1") then
      catm(2,i)(1:4)=' O1P'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"OP1 to O1P"
      endif
      endif
      if(catm(2,i)(1:4).eq." OP2") then
      catm(2,i)(1:4)=' O2P'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"OP2 to O2P"
      endif
      endif
      if(catm(2,i)(1:4).eq."HO3'") then
      catm(2,i)(1:4)=' H3T'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"HO3' to  H3T"
      endif
      endif
      if(catm(2,i)(1:4).eq."HO5'") then
      catm(2,i)(1:4)=' H5T'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"HO5' to  H5T"
      endif
      endif
      endif
   24 continue
c
c     DC   O3' ->  O3*
c          O5' ->  O5*
c          C5' ->  C5*
c          H5' -> 2H5*
c         H5'' -> 3H5*
c          C4' ->  C4*
c          H4' ->  H4*
c          O4' ->  O4*
c          C1' ->  C1*
c          H1' ->  H1*
c          C3' ->  C3*
c          H3' ->  H3*
c          C2' ->  C2*
c          H2' -> 2H2*
c         H2'' -> 3H2*
c         HO3' ->  H3T
c         HO5' ->  H5T
c          OP1 ->  O1P
c          OP2 ->  O2P
c          H41 ->  2H4
c          H42 ->  3H4
c
      cres='DC'
      do 25 i=1,natm
      if(cseq(1,latm(1,i))(1:3).eq.'DC '.or.
     +   cseq(1,latm(1,i))(1:3).eq.'DC3'.or.
     +   cseq(1,latm(1,i))(1:3).eq.'DC5') then
      if(catm(2,i)(1:4).eq." O3'") then
      catm(2,i)(1:4)=' O3*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"O3' to O3*"
      endif
      endif
      if(catm(2,i)(1:4).eq." O5'") then
      catm(2,i)(1:4)=' O5*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"O5' to O5*"
      endif
      endif
      if(catm(2,i)(1:4).eq." C5'") then
      catm(2,i)(1:4)=' C5*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"C5' to C5*"
      endif
      endif
      if(catm(2,i)(1:4).eq." H5'") then
      catm(2,i)(1:4)='2H5*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"H5' to 2H5*"
      endif
      endif
      if(catm(2,i)(1:4).eq."H5''") then
      catm(2,i)(1:4)='3H5*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"H5'' to 3H5*"
      endif
      endif
      if(catm(2,i)(1:4).eq." C4'") then
      catm(2,i)(1:4)=' C4*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"C4' to C4*"
      endif
      endif
      if(catm(2,i)(1:4).eq." H4'") then
      catm(2,i)(1:4)=' H4*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"H4' to H4*"
      endif
      endif
      if(catm(2,i)(1:4).eq." O4'") then
      catm(2,i)(1:4)=' O4*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"O4' to O4*"
      endif
      endif
      if(catm(2,i)(1:4).eq." C1'") then
      catm(2,i)(1:4)=' C1*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"C1' to C1*"
      endif
      endif
      if(catm(2,i)(1:4).eq." H1'") then
      catm(2,i)(1:4)=' H1*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"H1' to H1*"
      endif
      endif
      if(catm(2,i)(1:4).eq." C3'") then
      catm(2,i)(1:4)=' C3*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"C3' to C3*"
      endif
      endif
      if(catm(2,i)(1:4).eq." H3'") then
      catm(2,i)(1:4)=' H3*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"H3' to H3*"
      endif
      endif
      if(catm(2,i)(1:4).eq." C2'") then
      catm(2,i)(1:4)=' C2*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"C2' to C2*"
      endif
      endif
      if(catm(2,i)(1:4).eq." H2'") then
      catm(2,i)(1:4)='2H2*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"H2' to 2H2*"
      endif
      endif
      if(catm(2,i)(1:4).eq."H2''") then
      catm(2,i)(1:4)='3H2*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"H2'' to 3H2*"
      endif
      endif
      if(catm(2,i)(1:4).eq."HO3'") then
      catm(2,i)(1:4)=' H3T'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"HO3' to  H3T"
      endif
      endif
      if(catm(2,i)(1:4).eq."HO5'") then
      catm(2,i)(1:4)=' H5T'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"HO5' to  H5T"
      endif
      endif
      if(catm(2,i)(1:4).eq." OP1") then
      catm(2,i)(1:4)=' O1P'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i)," OP1 to  O1P"
      endif
      endif
      if(catm(2,i)(1:4).eq." OP2") then
      catm(2,i)(1:4)=' O2P'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i)," OP2 to  O2P"
      endif
      endif
      if(catm(2,i)(1:4).eq." H41") then
      catm(2,i)(1:4)='2H4 '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i)," H41 to 2H4 "
      endif
      endif
      if(catm(2,i)(1:4).eq." H42") then
      catm(2,i)(1:4)='3H4 '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i)," H42 to 3H4 "
      endif
      endif
      endif
   25 continue
c
c     DG   O3' ->  O3*
c          O5' ->  O5*
c          C5' ->  C5*
c          H5' -> 2H5*
c         H5'' -> 3H5*
c          C4' ->  C4*
c          H4' ->  H4*
c          O4' ->  O4*
c          C1' ->  C1*
c          H1' ->  H1*
c          C3' ->  C3*
c          H3' ->  H3*
c          C2' ->  C2*
c          H2' -> 2H2*
c         H2'' -> 3H2*
c          OP1 ->  O1P
c          OP2 ->  O2P
c          H21 -> 2H2
c          H22 -> 3H2
c         HO3' ->  H3T
c         HO5' ->  H5T
c
      cres='DG'
      do 26 i=1,natm
      if(cseq(1,latm(1,i))(1:3).eq.'DG '.or.
     +   cseq(1,latm(1,i))(1:3).eq.'DG3'.or.
     +   cseq(1,latm(1,i))(1:3).eq.'DG5') then
      if(catm(2,i)(1:4).eq." O3'") then
      catm(2,i)(1:4)=' O3*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"O3' to O3*"
      endif
      endif
      if(catm(2,i)(1:4).eq." O5'") then
      catm(2,i)(1:4)=' O5*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"O5' to O5*"
      endif
      endif
      if(catm(2,i)(1:4).eq." C5'") then
      catm(2,i)(1:4)=' C5*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"C5' to C5*"
      endif
      endif
      if(catm(2,i)(1:4).eq." H5'") then
      catm(2,i)(1:4)='2H5*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"H5' to 2H5*"
      endif
      endif
      if(catm(2,i)(1:4).eq."H5''") then
      catm(2,i)(1:4)='3H5*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"H5'' to 3H5*"
      endif
      endif
      if(catm(2,i)(1:4).eq." C4'") then
      catm(2,i)(1:4)=' C4*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"C4' to C4*"
      endif
      endif
      if(catm(2,i)(1:4).eq." H4'") then
      catm(2,i)(1:4)=' H4*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"H4' to H4*"
      endif
      endif
      if(catm(2,i)(1:4).eq." O4'") then
      catm(2,i)(1:4)=' O4*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"O4' to O4*"
      endif
      endif
      if(catm(2,i)(1:4).eq." C1'") then
      catm(2,i)(1:4)=' C1*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"C1' to C1*"
      endif
      endif
      if(catm(2,i)(1:4).eq." H1'") then
      catm(2,i)(1:4)=' H1*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"H1' to H1*"
      endif
      endif
      if(catm(2,i)(1:4).eq." C3'") then
      catm(2,i)(1:4)=' C3*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"C3' to C3*"
      endif
      endif
      if(catm(2,i)(1:4).eq." H3'") then
      catm(2,i)(1:4)=' H3*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"H3' to H3*"
      endif
      endif
      if(catm(2,i)(1:4).eq." C2'") then
      catm(2,i)(1:4)=' C2*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"C2' to C2*"
      endif
      endif
      if(catm(2,i)(1:4).eq." H2'") then
      catm(2,i)(1:4)='2H2*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"H2' to 2H2*"
      endif
      endif
      if(catm(2,i)(1:4).eq."H2''") then
      catm(2,i)(1:4)='3H2*'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"H2'' to 3H2*"
      endif
      endif
      if(catm(2,i)(1:4).eq." OP1") then
      catm(2,i)(1:4)=' O1P'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i)," OP1 to  O1P"
      endif
      endif
      if(catm(2,i)(1:4).eq." OP2") then
      catm(2,i)(1:4)=' O2P'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i)," OP2 to  O2P"
      endif
      endif
      if(catm(2,i)(1:4).eq." H21") then
      catm(2,i)(1:4)='2H2 '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i)," H21 to 2H2"
      endif
      endif
      if(catm(2,i)(1:4).eq." H22") then
      catm(2,i)(1:4)='3H2 '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i)," H22 to 3H2"
      endif
      endif
      if(catm(2,i)(1:4).eq."HO3'") then
      catm(2,i)(1:4)=' H3T'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"HO3' to  H3T"
      endif
      endif
      if(catm(2,i)(1:4).eq."HO5'") then
      catm(2,i)(1:4)=' H5T'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),"HO5' to  H5T"
      endif
      endif
      endif
   26 continue
c
c     Fix N-terminal hydrogens
c
c     XAA  HT1 -> 2H
c          HT2 -> 3H
c          HT3 -> 4H
c
      do 27 i=1,natm
      cres = cseq(1,latm(1,i))(1:3)
      if(pre_protein(cseq(1,latm(1,i)))) then
      if(catm(2,i)(1:4).eq.' HT1') then
      catm(2,i)(1:4)='2H  '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HT1 to 2H'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HT2') then
      catm(2,i)(1:4)='3H  '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HT2 to 3H'
      endif
      endif
      if(catm(2,i)(1:4).eq.' HT3') then
      catm(2,i)(1:4)='4H  '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'HT3 to 4H'
      endif
      endif
      endif
   27 continue
c
c     Fix C-terminal oxygens
c
c     XAA  OT1 ->  O
c          OT2 ->  OXT
c
      do 28 i=1,natm
      cres = cseq(1,latm(1,i))(1:3)
      if(pre_protein(cseq(1,latm(1,i)))) then
      if(catm(2,i)(1:4).eq.' OT1') then
      catm(2,i)(1:4)=' O  '
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'OT1 to O'
      endif
      endif
      if(catm(2,i)(1:4).eq.' OT2') then
      catm(2,i)(1:4)=' OXT'
      if(util_print('sequence',print_high)) then
      write(lfnout,2002) cres,latm(1,i),'OT2 to OXT'
      endif
      endif
      endif
   28 continue
c
      do 29 i=1,natm
      if(latm(2,i).le.0) latm(2,i)=pre_atnum(catm(2,i)(1:2))
      if(latm(2,i).eq.0) then
      elemnt(1:1)=' '
      elemnt(2:2)=catm(2,i)(2:2)     
      latm(2,i)=pre_atnum(elemnt)
      if(util_print('sequence',print_high)) then
      write(lfnout,2004) latm(1,i)
 2004 format(/,'Unrecognized element in segment ',i4)
      endif
      endif
   29 continue
c
      pre_namiup=.true.
c
      return
      end
