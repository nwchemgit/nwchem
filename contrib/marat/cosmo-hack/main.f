       program cosmo_hack
       implicit none
       character*255 message
       character*255 solutefile
       character*255 cosmofile
       character*255 outfile
       character*16  tar
       character*4 a4
       character*255 buf
       logical ofile
       integer fns,fnc,fno      
       integer iat,ires
       integer nat
       integer i
       character*5, dimension(:), allocatable :: atag
       double precision, dimension(:,:), allocatable :: c
       double precision, dimension(:), allocatable :: q
c
       external util_get_io_unit
       logical util_get_io_unit
c      ---------------------
c      read input parameters
c      ---------------------
       inquire(file="cosmo_hack.dat",exist=ofile)
       if(ofile) then
          open(10,file="cosmo_hack.dat",
     $            form='formatted',status='unknown')
          read(10,*,ERR=911) solutefile
          read(10,*,ERR=911) cosmofile
          read(10,*,ERR=911) outfile
          close(10)
       else
          write(6,*) "Enter solute PDB file"
          read(5,*) solutefile
          write(6,*) "Enter cosmo file"
          read(5,*) cosmofile
          write(6,*) "Enter output file"
          read(5,*) outfile
       end if 
     
       write(*,*) solutefile
       write(*,*) cosmofile
       write(*,*) outfile

c      -------------
c      open io units
c      -------------
       if(.not.util_get_io_unit(fns)) then
           message = "no free file units"
           goto 911
       else
           message = "opening "//trim(solutefile)
           open(fns,file=trim(solutefile),
     $                form='formatted',status='old',err=911)
       end if 

       if(.not.util_get_io_unit(fnc)) then
           message = "no free file units"
           goto 911
       else
           message = "opening "//trim(cosmofile)
           open(fnc,file=trim(cosmofile),
     $                form='formatted',status='unknown',err=911)
       end if 
       if(.not.util_get_io_unit(fno)) then
           message = "no free file units"
           goto 911
       else
           message = "opening "//trim(outfile)
           open(fno,file=trim(outfile),
     $                form='formatted',status='unknown',err=911)
       end if 
c      -------------------------------------------
c      extract solute portion from solute pdb file
c      -------------------------------------------
       rewind(fns)
       do 
         message = "reading "//trim(solutefile)
         read(fns,200,ERR=911,END=911) buf
         write(*,*) trim(buf)
         IF(buf(1:3).eq."TER") exit
         IF(buf(1:3).eq."END") exit
         if(buf(1:4).eq."ATOM") then
           read(buf,*) a4,iat,a4,a4,ires
           write(*,*) "iat,ires",iat,ires
         end if
         write(fno,'(A)') trim(buf)
       end do
c      --------------------------------------------
c      insert cosmo charges as individual fragments
c      --------------------------------------------
       call xyz_read_natoms(nat,fnc)
       write(*,*) "nat",nat
       allocate(c(3,nat))
       allocate(q(nat))
       allocate(atag(nat))
       rewind(fnc)
       call xyzq_read(nat,c,q,atag,fnc)
       write(*,*) "nat",nat
       do i=1,nat
         iat = iat + 1
         ires = ires +1 
         write(a4,'(I3.3)') i
         write(fno,FMT=9000)iat," Q  ",a4,ires,c(1,i),c(2,i),c(3,i)
         write(*,*) "q=",q(i)
c        create fragment
         call create_simple_fragment(1,a4," Q    ","Q     ",q(i))
      end do
c      ------------------
c      add solvent if any
c      ------------------
      do 
         message = "reading "//trim(solutefile)
         read(fns,200,ERR=20,END=20) buf
         write(*,*) trim(buf)
         if(buf(1:4).ne."ATOM") exit
         write(fno,'(A)') trim(buf)
       end do
20     continue
       write(fno,'(A)') "END"
       stop
c      if you reach this you are in trouble
911    continue
       write(*,*) "Emergency STOP"
       write(*,*) message
       stop
c      format statements
c      -----------------
200    FORMAT(A255)
9000   FORMAT("ATOM",T7,I5,T13,A4,T18,A3,T23,I4,T31,
     >        F8.3,T39,F8.3,T47,F8.3,T55,F6.2)
        end

      function util_get_io_unit(fn)

      implicit none
      integer fn
      logical util_get_io_unit
c 
      integer k
      logical ostatus
c
      do k=80,90
        INQUIRE(UNIT=k,OPENED=ostatus)
        ostatus = .not.ostatus
        if(ostatus) 
     >    INQUIRE(UNIT=k,EXIST=ostatus)
        if(ostatus) then
          fn = k
          util_get_io_unit = .true.
          return
        end if 
      end do
      util_get_io_unit = .false.
      return
      end

       subroutine create_simple_fragment(n,ar,an,at,q)
       implicit none
       integer n             !> number of atoms
       character*(*) ar      !> name of the residue
       character(*) an(n)   !> name of atoms in the residue
       character*(*) at(n)   !> name of atom type in the residue
       double precision q(n) !> charges
c
       integer i
       integer fnf
       character*30 pname 
       character*255 message
       external util_get_io_unit
       logical util_get_io_unit
c
       pname = "create_simple_fragment"
       if(.not.util_get_io_unit(fnf)) then
           message = "no free file units"
           goto 911
       else
           message = "opening "//trim(ar)//".frg"
           open(fnf,file=trim(ar)//".frg",
     $                 form='formatted',status='unknown',err=911)
       end if         
       write(fnf,'(A)') "$"//trim(ar)
       write(fnf,'(4I5)') n,1,1,0
       write(fnf,'(A)') trim(ar)
       do i=1,n
         write(fnf,'(I5,A6,A6,5I5,2f12.6)') i,an(i),
     >         at(i),
     >         0,0,0,1,1,q(i),0.0d0
       end do

       close(fnf)
       write(*,*) "out "//pname
       return
c      -------------
c      error section
c      -------------
911    continue
       write(*,*) "ERROR STOP"
       write(*,*) "subroutine:" //pname(1:len_trim(pname))
       write(*,*) "message" //message(1:len_trim(message))
       stop
c
c      format statements
c      -----------------
1030     FORMAT(180A1)
1020     FORMAT(A180)
       end
       subroutine xyz_read_natoms(n,fn_in)
       implicit none
       integer n
       integer fn_in
c       
       character*30 pname
       character*180 bigbuf
       logical is_integer
       external is_integer
c
       pname = "xyz_read_natoms"
       n = 0
c      ------------------------------------------
c      get number of atoms (skipping empty lines)
c      ------------------------------------------
       bigbuf = " "
       do 
         read(fn_in,20,end=10) bigbuf
         if(is_integer(bigbuf)) exit
       end do
       read(bigbuf,*,err=136) n
10     continue
       return
c      -------------
c      error section
c      -------------
136    continue
       write(*,*) "error in ",pname(1:len_trim(pname))
       write(*,*) "current buffer: ",bigbuf(1:len_trim(bigbuf))
       stop
c
c      format statements
c      -----------------
20     FORMAT(A180)
       end

      function is_integer(string)
      implicit none
      character*(*) string
      logical is_integer
c
      character(len=10) :: anumber = "0123456789"
      integer i,l
      is_integer = .false.
      if(string.eq." ") return
      l  = len_trim(string)
      do i=1,l
       if(string(i:i).ne." ") exit
      end do
      is_integer = verify(string(i:l), anumber).eq.0
      end function


       subroutine xyzq_read(n0,c,q,atag,fn_in)
       implicit none
       integer n0
       double precision c(3,n0)
       double precision q(n0)
       integer fn_in
       character*(*) atag(n0)
c
       integer i,k,n
       character*180 bigbuf
       character*180 message
       logical is_integer
       external is_integer
       character*30 pname 
c
       pname = "xyz_read"
       c = 0.0d0
       atag = " "
       n=0
c
       call xyz_read_natoms(n,fn_in)
       if(n.eq.0) then
          n0=0
          goto 20
       end if
       if(n.gt.n0) then
          message = "too small array size"
          goto 911
       end if
c
c      -------------------
c      get title field(if any) 
c      -------------------
       read(fn_in,1020) bigbuf
c
c       ------------------------------------
c       read coordinates(skipping empty lines)
c       ------------------------------------
       do i=1,n
        do
          read(fn_in,1020,end=20) bigbuf
          if(bigbuf.eq."") then
             cycle
          else if(bigbuf(1:1).eq."#") then
             cycle
          else
             exit
          end if
        end do
        message = "reading "//bigbuf
        read(bigbuf,*,err=911) atag(i),(c(k,i),k=1,3),q(i)
        n0=i
       end do
c      check if all atoms were read
       if(n0.lt.n) then
         message = "could not find all the atoms"
         goto 911
       end if
20     continue
       write(*,*) "out "//pname
       return
c      -------------
c      error section
c      -------------
911    continue
       write(*,*) "ERROR STOP"
       write(*,*) "subroutine:" //pname(1:len_trim(pname))
       write(*,*) "message" //message(1:len_trim(message))
       stop
c
c      format statements
c      -----------------
1030     FORMAT(180A1)
1020     FORMAT(A180)
       end
c $Id$
