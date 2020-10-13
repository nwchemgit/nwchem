      program pdb_amber2nwchem
      implicit none
      character(len=1024) :: inputfile
      character(len=1024) :: outputfile
      character(len=100)  :: line
      character(len=20)   :: l1        ! characters  1-20
      character(len=1)    :: chain     ! characters 22
      character(len=6)    :: resnum    ! characters 23-26
      character(len=100)  :: l2        ! characters 27-72
      character(len=6)    :: old_resnum
      integer, external   :: inp_strlen
      logical, external   :: pre_namiup
      integer :: lfninp
      integer :: lfnout
      integer :: iline
      integer :: resid
      logical :: largepdb
      logical :: writtenlarge
      call pdb_readcommandlinearguments(inputfile,outputfile)
      lfninp = 5
      lfnout = 6
      largepdb = .false.
      writtenlarge = .false.
      old_resnum = "-1"
      resid = 0
      chain = "A"
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
        iline = iline + 1
        read(lfninp,'(a)')line
        select case (line(1:4))
          case ("LRGP") ! LRGPDB
            largepdb = .true.
          case ("ATOM", "HETA")
            if (largepdb) then
              read(line,2002,end=300,err=300)l1,resid,l2
              write(lfnout,3002)l1,chain,mod(resid-1,9999)+1,
     +                          l2(1:inp_strlen(l2))
            else
              if (.not.writtenlarge) then
                write(lfnout,'("LRGPDB")')
                writtenlarge = .true.
              endif
              read(line,1002,end=300,err=300)l1,chain,resnum,l2
              if (resnum.ne.old_resnum) then
                resid = resid + 1
                old_resnum = resnum
              endif
              write(lfnout,2002)l1,resid,l2(1:inp_strlen(l2))
            endif
          case ("TER ")
            if (largepdb) then
              read(line,2002,end=300,err=300)l1,resid,l2
              write(lfnout,3002)l1,chain,mod(resid-1,9999)+1,
     +                          l2(1:inp_strlen(l2))
            else
              read(line,1002,end=300,err=300)l1,chain,resnum,l2
              write(lfnout,2002)l1,resid,l2(1:inp_strlen(l2))
            endif
            if (chain.eq."Z") then
              chain = "A"
            else
              chain = char(ichar(chain)+1)
            endif
          case default
            write(lfnout,'(a)')line(1:inp_strlen(line))
        end select
      enddo
 300  continue
      close(lfninp)
      close(lfnout)
 1002 format(a20,1x,a1,a4,a)
 2002 format(a20,i6,a)
 3002 format(a20,1x,a1,i4,a)
      end
