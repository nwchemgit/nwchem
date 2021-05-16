      program pdb_amber2nwchem
      implicit none
#include "pre_common.fh"
      character(len=1024) :: inputfile
      character(len=1024) :: outputfile
      character(len=100)  :: line
      character(len=12)   :: l1        ! characters  1-12
      character(len=6)    :: catm(3,1) ! characters 13-16
      character(len=1)    :: l2        ! characters 17
      character(len=5)    :: cseq(2,1) ! characters 18-20
      character(len=100)  :: l3        ! characters 21-72
      integer, external   :: inp_strlen
      logical, external   :: pre_namiup
      integer :: lfninp
      integer :: lfnout
      integer :: iline
      integer :: lseq(6,1)
      integer :: latm(5,1)
      logical :: ll
      ffield = "amber"
      call pdb_readcommandlinearguments(inputfile,outputfile,ffield)
      lfninp = 5
      lfnout = 6
      lseq   = 1
      latm   = 1
      if (inputfile.ne."-") then
        close(lfninp)
        open(unit=lfninp,file=inputfile,err=100,form="formatted",
     +       status="old")
      endif
 100  continue
      if (outputfile.ne."-") then
        close(lfnout)
        open(unit=lfnout,file=outputfile,err=200,form="formatted")
      endif
 200  continue
      iline = 0
      line  = ""
      do while (line(1:3).ne."END")
        catm = "      "
        cseq = "     "
        iline = iline + 1
        read(lfninp,'(a)')line
        select case (line(1:4))
          case ("ATOM  ", "HETA")
            read(line,1002,end=300,err=300)l1,catm(1,1)(1:4),l2,
     +                                     cseq(1,1)(1:3),l3
            ll = pre_namiup(0,lseq,cseq,1,1,latm,catm,1,1)
            call namseq(catm,cseq)
            write(lfnout,1002)l1,catm(2,1)(1:4),l2,cseq(2,1)(1:3),
     +                        l3(1:inp_strlen(l3))
          case default
            write(lfnout,'(a)')line(1:inp_strlen(line))
        end select
      enddo
 300  continue
      close(lfninp)
      close(lfnout)
 1002 format(a12,a4,a1,a3,a)
      end
!
!-----------------------------------------------------------------------
!
      subroutine namseq(catm,cseq)
      implicit none
      character(len=6), intent(inout) :: catm(3,1)
      character(len=5), intent(inout) :: cseq(2,1)
!
      cseq(2,1)(1:3) = cseq(1,1)(1:3)
      if (cseq(2,1)(1:3).eq."WAT") then
        cseq(2,1)(1:3) = "HOH"
        if (catm(2,1)(1:4).eq." H1 ") then
          catm(2,1)(1:4) = "2H  "
        endif
        if (catm(2,1)(1:4).eq." H2 ") then
          catm(2,1)(1:4) = "3H  "
        endif
        if (catm(2,1)(1:4).eq." H3 ") then
          catm(2,1)(1:4) = "2H  "
        endif
      endif
      if (cseq(2,1)(1:3).eq."Cl-") then
        cseq(2,1)(1:3) = " Cl"
        if (catm(2,1)(1:4).eq."Cl- ") then
          catm(2,1)(1:4) = "Cl  "
        endif
      endif
      if (cseq(2,1)(1:3).eq."Na+") then
        cseq(2,1)(1:3) = " Na"
        if (catm(2,1)(1:4).eq."Na+ ") then
          catm(2,1)(1:4) = "Na  "
        endif
      endif
      if (cseq(2,1)(1:3).eq." K+") then
        cseq(2,1)(1:3) = " K "
        if (catm(2,1)(1:4).eq." K+ ") then
          catm(2,1)(1:4) = " K  "
        endif
      endif
      return
      end
