       program main
       implicit none
#include "amber_params.fh"       
       logical ofile
       integer nmax
       parameter(nmax=20)
       double precision lmax(nmax),sigma(nmax),sum(nmax)
       double precision f(nmax),e(nmax)
       character*(180) buffer
       integer itype
       integer i,n
       double precision l1,l2,tmp
       integer u_angle
       integer u_bond
       integer u_dihedral
       integer u_improper
       integer u_mass
       integer u_nonbonded
       integer u_param
c
       buffer = "amber.in"
       inquire(file=buffer,exist=ofile)
       if(ofile) then
          open(10,file=buffer,
     $            form='formatted',status='old')
       else
         write(6,*) "no file found"//buffer
         return
       end if 

       u_bond = 11
       open(u_bond,file="bond.dat",
     $            form='formatted',status='unknown')
       write(u_bond,'(A5)') 
     $ "Bonds"
       u_angle = 12
       open(u_angle,file="angle.dat",
     $            form='formatted',status='unknown')
       write(u_angle,'(A6)') 
     $ "Angles"
       u_dihedral = 13
       open(u_dihedral,file="dihedral.dat",
     $            form='formatted',status='unknown')
       write(u_dihedral,'(A16)') 
     $ "Proper dihedrals"
       u_improper = 14
       open(u_improper,file="improper.dat",
     $            form='formatted',status='unknown')
       write(u_improper,'(A18)') 
     $ "Improper dihedrals"
       u_mass     = 15
       open(u_mass,file="mass.dat",
     $            form='formatted',status='unknown')
       u_nonbonded     = 16
       open(u_nonbonded,file="nonbonded.dat",
     $            form='formatted',status='unknown')
       write(u_nonbonded,'(A5)') 
     $ "Atoms"
       u_param     = 17
       open(u_param,file="amber.par",
     $            form='formatted',status='unknown')
       buffer = " "
c
c      parse all the sections except nonbonding
c      ----------------------------------------
10     continue 
       read(10,'(A180)',ERR=30,END=30) buffer
       call test_for_type(buffer,itype)
       if(buffer.eq." ") goto 10
       if(itype.eq.tbond) then
         call parse_bond(u_bond,buffer)
       else if(itype.eq.tmass) then
         call parse_mass(u_mass,buffer)
       else if(itype.eq.tangle) then
         call parse_angle(u_angle,buffer)
       else if(itype.eq.tdihedral) then
         call parse_dihedral(u_dihedral,buffer)
       else if(itype.eq.timproper) then
         call parse_improper(u_improper,buffer)
       end if
       go to 10
30     continue
c
c      now parse nonbonding section
c      ----------------------------
       write(*,*) "parsing nonbonding section"
       rewind(10)
11     continue
       read(10,'(A180)',ERR=31,END=31) buffer
       call test_for_type(buffer,itype)
       if(buffer.eq." ") goto 11
       if(itype.eq.tnonbond) then
         call parse_nonbonded(u_nonbonded,u_mass,buffer)
       end if
       goto 11
31     continue
       close(10)
       call append_section(u_param,u_nonbonded)
       call append_section(u_param,u_bond)
       call append_section(u_param,u_angle)
       call append_section(u_param,u_dihedral)
       call append_section(u_param,u_improper)
       end
      subroutine test_for_type(buffer,itype)
      implicit none
#include "amber_params.fh"
      character*180 buffer
      integer itype
      if(INDEX(buffer(1:4),"BOND").ne.0) then
        itype = tbond 
        buffer=" "
        write(*,*) "found BOND"
      else if(INDEX(buffer(1:4),"MASS").ne.0) then
        itype = tmass
        write(*,*) "found MASS"
        buffer=" "
      else if(INDEX(buffer(1:5),"ANGLE").ne.0) then
        itype = tangle
        write(*,*) "found ANGLE"
        buffer=" "
      else if(INDEX(buffer(1:8),"NONBOND").ne.0) then
        itype = tnonbond
        write(*,*) "found NONBOND"
        buffer=" "
      else if(INDEX(buffer(1:8),"IMPROPER").ne.0) then
        itype = timproper
        write(*,*) "found IMPROPER"
        buffer=" "
      else if(INDEX(buffer(1:8),"DIHEDRAL").ne.0) then
        itype = tdihedral
        buffer=" "
        write(*,*) "found DIHEDRAL"
      end if
      end subroutine
      subroutine parse_mass(un,buffer)
      implicit none
#include "amber_params.fh"
      integer un
      character*180 buffer
      integer itype
      character*5 atyp1
      double precision m
      
      read(buffer,*) atyp1,m
      write(un,*) atyp1,m
      end subroutine
      subroutine locate_mass(un,atyp,m)
      implicit none
#include "amber_params.fh"
      integer un
      character*5 atyp   ![IN] atom name
      double precision m ![OUT] mass
c      
      character*5 atyp1
      m = -1.0
      rewind(un)
10    continue
      read(un,*,ERR=30,END=30) atyp1,m
      if(atyp.eq.atyp1) return
      go to 10
30    continue
      return
      end subroutine
      subroutine parse_bond(un,buffer)
      implicit none
#include "amber_params.fh"
      integer un
      character*180 buffer
      integer itype
      character*5 atyp1,atyp2
      double precision p1,p2
      double precision conv1,conv2
      
      read(buffer,1000) atyp1,atyp2,p2,p1
 1000 FORMAT(A2,1X,A2,2F10.2)       
      conv1 = ang_to_nm
      conv2 = 2*kcal_to_kj/(ang_to_nm)**2
      write(*,*) "conv1,conv2",conv1,conv2
      p1 = p1*conv1
      p2 = p2*conv2
      write(un,1001) atyp1,atyp2,p1,p2
 1001 format(a5,"-",a5,f10.5,e12.5)
      end subroutine
      subroutine parse_angle(un,buffer)
      implicit none
#include "amber_params.fh"
      integer un
      character*180 buffer
      integer itype
      character*5 atyp1,atyp2,atyp3
      double precision p1,p2
      double precision conv1,conv2
      
      read(buffer,1000) atyp1,atyp2,atyp3,p2,p1
 1000 FORMAT(A2,1X,A2,1X,A2,2F10.2)       
      conv1 = deg_to_rad
      conv2 = 2*kcal_to_kj
      write(*,*) "conv1,conv2",conv1,conv2
      p1 = p1*conv1
      p2 = p2*conv2
      write(un,1001) atyp1,atyp2,atyp3,p1,p2
 1001 format(a5,"-",a5,"-",a5,f10.5,e12.5)
      end subroutine
      subroutine parse_dihedral(un,buffer)
      implicit none
#include "amber_params.fh"
      integer un
      character*180 buffer
      integer itype
      character*5 atyp1,atyp2,atyp3,atyp4
      integer idivf
      double precision pk,phase,pn
      double precision conv1,conv2
      
      read(buffer,1000) atyp1,atyp2,atyp3,atyp4,idivf,pk,phase,pn
 1000 FORMAT(A2,1X,A2,1X,A2,1X,A2,I4,3F15.2)
      conv1 = deg_to_rad
      conv2 = kcal_to_kj
      write(*,*) "conv1,conv2",conv1,conv2
      phase = phase*conv1
      pk = pk*conv2/idivf
      if(atyp1.eq."X") atyp1=" "
      if(atyp4.eq."X") atyp4=" "
      write(un,1001) atyp1,atyp2,atyp3,atyp4,phase,pk,INT(pn)
 1001 format(a5,"-",a5,"-",a5,"-",a5,f10.5,e12.5,i5)
      end subroutine
      subroutine parse_improper(un,buffer)
      implicit none
#include "amber_params.fh"
      integer un
      character*180 buffer
      integer itype
      character*5 atyp1,atyp2,atyp3,atyp4
      integer idivf
      double precision pk,phase,pn
      double precision conv1,conv2
      
      read(buffer,1000) atyp1,atyp2,atyp3,atyp4,pk,phase,pn
 1000 FORMAT(A2,1X,A2,1X,A2,1X,A2,3F15.2)
      conv1 = deg_to_rad
      conv2 = kcal_to_kj
      write(*,*) "conv1,conv2",conv1,conv2
      phase = phase*conv1
      pk = pk*conv2
      if(atyp1.eq."X") atyp1=" "
      if(atyp2.eq."X") atyp2=" "
      write(un,1001) atyp1,atyp2,atyp3,atyp4,phase,pk,INT(pn)
 1001 format(a5,"-",a5,"-",a5,"-",a5,f10.5,e12.5,i5)
      end subroutine
      subroutine parse_nonbonded(un,un_mass,buffer)
      implicit none
#include "amber_params.fh"
      integer un,un_mass
      character*180 buffer
      character*180 buffer1
      integer itype
      character*5 atyp1
      integer idivf
      double precision p1,p2
      double precision rmass
      double precision conv1,conv2
      integer atn
      
      read(buffer,*) atyp1,p1,p2
      conv1 = ang_to_nm
      conv2 = kcal_to_kj
      write(*,*) "conv1,conv2",conv1,conv2
      p1 = p1*conv1
      p2 = p2*conv2
      call locate_mass(un_mass,atyp1,rmass)
      call get_atn_from_mass(atn,rmass)
      write(un,1001) atyp1,rmass,p2,p1
 1001 format(a5,f10.5,2e12.5,24X,4x,"1",1x,"1111111111")
      write(un,1002) atn,0.5d0*p2,p1
 1002 format(10x,i5,2e12.5)
      end subroutine
      subroutine  get_atn_from_mass(atn,mass)
c
      implicit none
c
      integer atn            
      double precision mass  
c
      double precision def_masses(109)
      double precision tolerance
c
      integer i
c
      data (def_masses(i),i=1,50) /
     &  1.007825d0, 4.0026d0,    7.016d0,    9.01218d0, 11.00931d0,
     & 12.0d0,     14.00307d0,  15.99491d0, 18.9984d0,  19.99244d0,
     & 22.9898d0,  23.98504d0,  26.98154d0, 27.97693d0, 30.97376d0,
     & 31.97207d0, 34.96885d0,  39.9624d0,  38.96371d0, 39.96259d0,
     & 44.95592d0, 45.948d0,    50.9440d0,  51.9405d0,  54.9381d0,
     & 55.9349d0,  58.9332d0,   57.9353d0,  62.9298d0,  63.9291d0,
     & 68.9257d0,  73.9219d0,   74.9216d0,  79.9165d0,  78.9183d0,
     & 83.912d0,   84.9117d0,   87.9056d0,  88.9054d0,  89.9043d0,
     & 92.9060d0,  97.9055d0,   97.9072d0, 101.9037d0, 102.9048d0,
     &105.9032d0, 106.90509d0, 113.9036d0, 114.9041d0, 117.9018d0/
      data (def_masses(i),i=51,109) /
     & 120.9038d0, 129.9067d0, 126.9004d0, 131.9042d0, 132.9051d0,
     & 137.9050d0, 138.9061d0, 139.9053d0, 140.9074d0, 143.9099d0,
     & 144.9128d0, 151.9195d0, 152.9209d0, 157.9241d0, 159.9250d0,
     & 163.9288d0, 164.9303d0, 165.9304d0, 168.9344d0, 173.9390d0,
     & 174.9409d0, 179.9468d0, 180.948d0,  183.9510d0, 186.9560d0,
     & 189.9586d0, 192.9633d0, 194.9648d0, 196.9666d0, 201.9706d0,
     & 204.9745d0, 207.9766d0, 208.9804d0, 209.9829d0, 210.9875d0,
     & 222.0175d0, 223.0198d0, 226.0254d0, 227.0278d0, 232.0382d0,
     & 231.0359d0, 238.0508d0, 237.0482d0, 244.0642d0, 243.0614d0,
     & 247.0704d0, 247.0703d0, 251.0796d0, 252.0829d0, 257.0950d0,
     & 258.0986d0, 259.1009d0, 262.1100d0, 261.1087d0, 262.1138d0,
     & 266.1219d0, 262.1229d0, 267.1318d0, 268.1388d0 /
c
      tolerance = 0.2
      do i=1,109
       if(abs(mass-def_masses(i)).lt.tolerance) then
         atn = i
         return
       end if
      end do
c
      end

      subroutine append_section(un,un1)
      implicit none
#include "amber_params.fh"
      integer un
      integer un1
c      
      character*80 buffer
      rewind(un1)
10    continue
      read(un1,'(A80)',ERR=30,END=30) buffer
      write(un,*) buffer
      go to 10
30    continue
      return
      end subroutine
c      function geom_tag_to_element(tag, symbol, element, atn)
c      implicit none
c      logical geom_tag_to_element
c      character*2 symbols(nelements)
c      character*16 elements(nelements)
c      character*16 tag          ! [input]
c      character*(*) symbol      ! [output]
c      character*(*) element     ! [output]
c      integer atn               ! [output]
cc
cc     attempt to figure out which element a tag refers to
cc     and return the symbol, name and atomic no.
cc
c      integer lbuf, ind
c      character*16 buf
c      character*1 sym1(14)      ! 1 character atomic symbols+atomic no.s
c      integer atn1(14)
c      character*2 amber_types_name(40)
c      integer amber_types_atn(40)
c
c      data amber_types_name/ 
c     $ 'CT','C ','CA','CM','CC','CV','CW','CR','CB','C*','CN','CK','CQ',
c     $ 'N ' ,'NA' ,'NB' ,'NC' ,'N*' ,'N2' ,'N3',
c     $ 'OW' ,'OH' ,'OS' ,'O ' ,'O2', 
c     $ 'S ' ,'SH' ,
c     $ 'P ',
c     $ 'H ','HW','HO','HS','HA','HC','H1','H2','H3','HP','H4','H5'/ 
c      data amber_types_atn/ 
c     $ 6,6,6,6,6,6,6,6,6,6,6,6,6,
c     $ 7,7,7,7,7,7,7,
c     $ 8,8,8,8,8,
c     $ 16,16,
c     $ 15,
c     $ 1,1,1,1,1,1,1,1,1,1,1,1/
c
c      data sym1/'H','B','C','N','O','F','P','S','K','V','Y','I','W','U'/
c      data atn1/ 1 , 5 , 6 , 7 , 8 , 9 , 15, 16, 19, 23, 39, 53, 74, 92/
c
c      geom_tag_to_element = .false.
cc
cc     eliminate conventions that refer to centers used for
cc     computation purposes .. just bq and x for now
cc
c      buf = tag
c      lbuf = inp_strlen(buf)
c      if (lbuf .eq. 0) return
cc
c      call inp_lcase(buf)
c      if (buf(1:2) .eq. 'bq' .or. 
c     $     ((buf(1:1).eq.'x') .and. (buf(2:2).ne.'e'))) then ! X but not Xe
c         element = 'point charge' ! Note that false is returned
c         symbol  = 'bq'
c         atn     = 0
c         return
c      end if
cc
cc     Attempt to match the first 4 characters of the
cc     full names of the elements
cc
c      atn = 0
c      if (lbuf .ge. 4) then
c        do ind = 1,nelements
c          if (inp_compare(.false.,buf(1:4),elements(ind)(1:4))) then
c            symbol  = symbols(ind)
c            element = elements(ind)
c            atn     = ind
c            geom_tag_to_element = .true.
c            return
c          endif
c        enddo
c      end if
cc
cc     Failed ... attempt to match the first two characters
cc     against two character element names
cc
c      if (buf(2:2) .ne. ' ') then
c         if (inp_match(nelements,.false.,buf(1:2),symbols,ind)) then
c            symbol  = symbols(ind)
c            element = elements(ind)
c            atn     = ind
c            geom_tag_to_element = .true.
c            return
c         end if
c      end if
cc
cc     Last ditch attempt ... match against 1 character symbols
cc
c      if (inp_match(14, .false., buf(1:1), sym1, ind)) then
c         ind = atn1(ind)
c         symbol  = symbols(ind)
c         element = elements(ind)
c         atn     = ind
c         geom_tag_to_element = .true.
c         return
c      end if
cc
c      if (inp_match(14, .false., buf(2:2), sym1, ind)) then
c         ind = atn1(ind)
c         symbol  = symbols(ind)
c         element = elements(ind)
c         atn     = ind
c         geom_tag_to_element = .true.
c         return
c      end if
ccc
cc     Nothing matched
cc
c      symbol = ' '
c      element = ' '
c      atn = 0
c      return
cc
c      end

