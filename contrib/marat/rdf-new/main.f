       program pair_correlation
       implicit none
       integer i,l
       double precision lat(3,3),latv(3)
       double precision rmax
       integer i1,nb
       character*255 infile
       character*5 atag
       logical ofile
       CHARACTER(LEN=11) :: anumber = "0123456789."
       logical is_integer
       external is_integer
       logical is_number
       external is_number
       character*(180) buffer
       character*(180) message
       integer istatus
       logical overb,ohelp
       character*255 file_lattice,file_in,file_out
       character*16 atom1_tag
       character*16 atom2_tag
       character*5 aformat
       integer atom1_id,atom2_id
       integer k 
       double precision f
       integer fn_in
       integer fn_out
       integer n,n1,n2
c      allocatable arrays
       double precision, dimension(:,:), allocatable :: c
       double precision, dimension(:,:), allocatable :: c1
       double precision, dimension(:,:), allocatable :: c2
       double precision, dimension(:), allocatable :: gr
c
c      --------------------------------------------      
c      beging parsing command line arguments if any
c      --------------------------------------------      
       aformat = " "
       file_in  = " "
       file_out = " "
       overb = .false.
       file_lattice = " "
       latv = -1.0
       lat = 0.0
       atom1_id = -1
       atom2_id = -1
       atom1_tag = " "
       atom2_tag = " "
       file_in = " "
       file_out = " "
       rmax = -1
       nb = -1
       i = 0
16     continue
       i = i+1
       call get_command_argument(i,buffer,l,istatus)
       if(istatus.ne.0) goto 18
c       write(*,*) "argument ",i,buffer
       if(buffer.eq."-lattice") then
c          write(*,*) "reading lattice"
          do k=1,3
            i = i+1
            call get_command_argument(i,buffer,l,istatus)
            if(istatus.ne.0) goto 18
c            write(*,*) "argument ",i,buffer,is_number(buffer)
            if(is_number(buffer)) then
              read(buffer,*) f
              latv(k:3) = f
c              write(*,*) "lat",k,latv(k)
            else
              i = i-1
              cycle
            end if
          end do
c          write(*,*) "done reading lattice"
          go to 16
       else if(buffer.eq."-v") then
          overb=.true.
          go to 16
c       else if(buffer.eq."-help") then
c          ohelp=.true.
c          write(*,1000)
c          go to 14
       else if(buffer.eq."-atom1") then
          i = i+1
          call get_command_argument(i,buffer,l,istatus)
          if(istatus.ne.0) goto 18
          if(is_number(buffer)) then
            message = "Only atom tags are supported now"
            goto 911
            read(buffer,*) atom1_id
          else
            atom1_tag = buffer
          end if
          go to 16
       else if(buffer.eq."-atom2") then
          i = i+1
          call get_command_argument(i,buffer,l,istatus)
          if(istatus.ne.0) goto 18
          if(is_number(buffer)) then
            message = "Only atom tags are supported now"
            goto 911
            read(buffer,*) atom2_id
          else
            atom2_tag = buffer
          end if
          go to 16
       else if(buffer.eq."-rmax") then
          i = i+1
          call get_command_argument(i,buffer,l,istatus)
          if(istatus.ne.0) goto 18
          if(is_number(buffer)) then
            read(buffer,*) rmax
          else
            message = "Maximum radius gas to be a number"
            goto 911
          end if
          go to 16
        else if(buffer.eq."-nbins") then
          i = i+1
          call get_command_argument(i,buffer,l,istatus)
          if(istatus.ne.0) goto 18
          if(is_integer(buffer)) then
            read(buffer,*) nb
          else
            message = "Number of bins have to be integer"
            goto 911
          end if
          go to 16
       else 
          if(file_in .eq. ' ') then
            file_in=buffer
            go to 16
          else if(file_out.eq." ") then
            file_out=buffer
            go to 16
          end if
       end if 
c      ---------------------------      
c      end of command line parsing
c      ---------------------------      
18     continue
       rmax = 4.0
       nb = 100
       atom1_tag = "O1"
       atom2_tag = "O"
       latv = 10.0014453
       file_in = "test.xyz"
       file_out = "test.out"
c      ---------------------------      
c      start checks/balances
c      ---------------------------      
c      atom tags
       if(atom1_tag.eq." ") then
         message = "please provide central atom tag"
         goto 911
       end if
       if(atom2_tag.eq." ") then
         atom2_tag = "*"
       end if
c      the lattice
       if(all(latv.lt.0)) then
          message = "problems with lattice input"
          goto 911
       end if
       lat = 0.0
       do k=1,3
         lat(k,k) = latv(k)
       end do
c      maximum radius
       if(rmax.lt.0) then
         rmax = 0.5*MINVAL(latv)
       end if
c      number of bins
       if(nb.lt.0) then
         nb = INT(rmax/0.1)
       end if
c      input file
       if(file_in.eq." ") then
         file_in = "traj.xyz"
       end if
c      output file
       if(file_out.eq." ") then
         file_out = "gr.dat"
       end if
c      open io channels
       inquire(file=file_in,exist=ofile)
       fn_in = 10
       if(ofile) then
          open(fn_in,file=file_in,
     $            form='formatted',status='old',err=911)
       else
           message = "no file found: "//file_in
           goto 911
       end if 

       fn_out = 11
       if(ofile) then
               message = "opening output file "//file_out
          open(fn_out,file=file_out,
     $            form='formatted',status='unknown',err=911)
       else
           message = "cannot open output file: "//file_out
           goto 911
       end if 

       if(overb) then
         write(*,*) "lattice vectors"
         write(*,*) (lat(i,1),i=1,3)
         write(*,*) (lat(i,2),i=1,3)
         write(*,*) (lat(i,3),i=1,3)
         write(*,*) "1st atom tag: ",atom1_tag
         write(*,*) "2nd atom tag: ",atom2_tag
         write(*,*) "Number of bins:",nb
         write(*,*) "Maximum radius:",rmax
         write(*,*) "Input file: " ,file_in
         write(*,*) "Output file: ",file_out
       end if

c      figure out format for trajectory file
       if(aformat.eq." ") then
        i=INDEX(file_in,".",.true.) 
        aformat=file_in(i+1:)
        write(*,*) "format is ",aformat
       end if

       call xyz_read_natoms(n,fn_in)
       rewind(fn_in)
       allocate(c(n,3))
       allocate(c1(n,3))
       allocate(c2(n,3))
       allocate(gr(nb))
       write(*,*) "number of atoms", n
c       call xyz_read_coords_byname(n,c,atom1_tag,fn_in)
30     continue
       n1=n
       n2=n
       call xyz_read_coords_byname2(n1,c1,atom1_tag,
     +                              n2,c2,atom2_tag,
     +                              fn_in)
       if(n1.eq.0.or.n2.eq.0) goto 31
       do i=1,n1
         write(34,*) (c1(i,k),k=1,3)
       end do
       do i=1,n2
         write(35,*) (c2(i,k),k=1,3)
       end do
c
       call rdf_compute(n, n1,c1,atom1_tag,
     +                  n2,c2,atom2_tag,
     +                  nb,rmax,lat,
     +                  gr,
     +                  fn_in)
       goto 30
31     continue
       write(*,*) "came to the end of the file"
       stop

       call pair_compute(rmax,lat,i1,atag,nb,infile)
       return
911    continue       
c      if you reach this you are in trouble
       write(*,*) "Emergency STOP"
       write(*,*) message
       stop
       end program

       subroutine xyz_read_natoms(n,fn_in)
       implicit none
       integer n
       integer fn_in
c       
       character*180 bigbuf
       character*180 message
       integer l
       logical is_integer
       external is_integer
c      ------------------------------------------
c      get number of atoms (skipping empty lines)
c      ------------------------------------------
11     continue
       message = "looking for number of atoms"
       read(fn_in,20,end=138) bigbuf
c      increment line number
       l=l+1
       if(.not.is_integer(bigbuf)) goto 11
       read(bigbuf,*,err=136) n
       return
c      -------------
c      error section
c      -------------
136    continue
       write(*,*) "error number of atoms field"
       write(*,*) "line number: ", l
       write(*,*) "buffer: ",bigbuf
       stop
       return
138    continue
       write(*,*) message
       stop
c
c      format statements
c      -----------------
30     FORMAT(180A1)
20     FORMAT(A180)
       end

       subroutine xyz_read_coords_byname(n0,c,atag,fn_in)
       implicit none
       integer n0
       double precision c(n0,3)
       integer fn_in
       character*(*) atag
c       
       double precision rmax
       integer i1,nb
c
       integer j
       integer i,n,a,nd,k
       double precision rlat(3,3),r1(3)
       character*30 buf
       character*30 tag
       character*180 bigbuf
       character*180 message
       integer l,sl
       logical is_integer
       external is_integer
       character*30 pname 
c
       pname = "xyz_read_coords_byname"
c      ------------------------------------------
c      get number of atoms (skipping empty lines)
c      ------------------------------------------
11     continue
       message = "looking for number of atoms "//pname
       read(fn_in,20,end=138) bigbuf
       write(*,*) bigbuf(1:len_trim(bigbuf)),is_integer(bigbuf)
c      increment line number
       l=l+1
       write(*,*) "is it integer",is_integer(bigbuf)
       if(.not.is_integer(bigbuf)) goto 11
       write(*,*) bigbuf(1:len_trim(bigbuf))
       read(bigbuf,*,err=136) n
c
c      -------------------
c      get title field(if any) 
c      -------------------
       read(fn_in,20) bigbuf
       l=l+1
c      processing coodinates
       j = 0
c
c       ------------------------------------
c       read coordinates(skipping empty lines)
c       ------------------------------------
       do i=1,n
        read(fn_in,20) bigbuf
        l=l+1
        j = j+1
        message = "verifying array size"
        if(j.gt.n0) goto 138
        if(bigbuf.eq."") cycle
        read(bigbuf,*,err=137) tag,(c(j,k),k=1,3)
c        
        sl = len_trim(atag)
        if(index(tag,atag(1:sl)).eq.0) then
          j = j-1
        end if
       end do
       n0=j
       return
c      -------------
c      error section
c      -------------
136    continue
       write(*,*) "subroutine:" //pname(1:len_trim(pname))
       write(*,*) "error number of atoms field"
       write(*,*) "line number: ", l
       write(*,*) "buffer: ",bigbuf
       stop
       return
137    continue
       write(*,*) "error reading coords field"
       write(*,*) "line number: ", l
       write(*,*) "buffer: ", bigbuf
       stop
       return
138    continue
       write(*,*) message
       stop
c
c      format statements
c      -----------------
30     FORMAT(180A1)
20     FORMAT(A180)
       end

       subroutine xyz_read_coords_byname2(n1,c1,atag1,
     +                                    n2,c2,atag2,
     +                                    fn_in)
       implicit none
       integer n1,n2
       double precision c1(n1,3),c2(n2,3)
       integer fn_in
       character*(*) atag1,atag2
c       
       double precision rmax
       integer i1,nb
c
       integer j1,j2
       integer i,n,a,nd,k
       double precision c3(3)
       character*30 buf
       character*30 tag
       character*180 bigbuf
       character*180 message
       integer l
       logical is_integer
       external is_integer
       character*30 pname 
c
       pname = "xyz_read_coords_byname"
       write(*,*) "in "//pname
       j1 = 0
       j2 = 0
c      ------------------------------------------
c      get number of atoms (skipping empty lines)
c      ------------------------------------------
11     continue
       message = "looking for number of atoms "//pname
       read(fn_in,20,end=135) bigbuf
       write(*,*) bigbuf(1:len_trim(bigbuf)),is_integer(bigbuf)
c      increment line number
       l=l+1
       if(.not.is_integer(bigbuf)) goto 11
       read(bigbuf,*,err=136) n
c
c      -------------------
c      get title field(if any) 
c      -------------------
       read(fn_in,20) bigbuf
       l=l+1
c      processing coodinates
c
c       ------------------------------------
c       read coordinates(skipping empty lines)
c       ------------------------------------
       do i=1,n
        read(fn_in,20) bigbuf
        l=l+1
        message = "verifying array size"
        if(bigbuf.eq."") cycle
        read(bigbuf,*,err=137) tag,(c3(k),k=1,3)
c        
        if(index(tag,atag1(1:len_trim(atag1))).ne.0) then
          j1 = j1+1
          if(j1.gt.n1) goto 138
          c1(j1,:)=c3
        else if(index(tag,atag2(1:len_trim(atag2))).ne.0) then
          j2 = j2+1
          if(j2.gt.n2) goto 138
          c2(j2,:)=c3
        end if
       end do
135    continue
       n1=j1
       n2=j2
       write(*,*) "out "//pname
       return
c      -------------
c      error section
c      -------------
136    continue
       write(*,*) "subroutine:" //pname(1:len_trim(pname))
       write(*,*) "error number of atoms field"
       write(*,*) "line number: ", l
       write(*,*) "buffer: ",bigbuf
       stop
       return
137    continue
       write(*,*) "error reading coords field"
       write(*,*) "line number: ", l
       write(*,*) "buffer: ", bigbuf
       stop
       return
138    continue
       write(*,*) message
       stop
c
c      format statements
c      -----------------
30     FORMAT(180A1)
20     FORMAT(A180)
       end

       subroutine rdf_compute(n0,n1,c1,atag1,
     +                        n2,c2,atag2,
     +                        nb,rmax,lat,
     +                        gr,
     +                        fn_in)
       implicit none
       integer n0,n1,n2
       double precision c1(n0,3),c2(n0,3)
       integer fn_in
       character*(*) atag1,atag2
       integer nb
       double precision rmax
       double precision gr(nb)
       double precision lat(3,3)
c       
       integer i1,i2
       double precision rd(n2+n1,3)
c
       integer j
       integer i,n,a,nd,k
       double precision c3(3)
       double precision rlat(3,3)
       character*30 buf
       character*30 tag
       character*180 bigbuf
       character*180 message
       integer l
       logical is_integer
       external is_integer
       character*30 pname 
       double precision d,dr,const
       double precision fourpi,x
       double precision ru,rl,vol,norm
       pname = "rdf_compute"
c
c      ----------------------
c      compute lattice params
c      ----------------------
       call smd_lat_invrt(lat,rlat)
       call smd_latt_vol(lat,vol)
c      ----------------------------
c      construct relative distances
c      ----------------------------
c       do i=1,n1
c         write(64,*) (c1(i,k),k=1,3)
c       end do
c       do i=1,n2
c         write(65,*) (c2(i,k),k=1,3)
c       end do
c
       i=0
       do i1=1,n1
         do i2=1,n2
           i = i+1
           rd(i,:) = c2(i2,:)-c1(i1,:)
            write(57,*) rd(i,:)
         end do
       end do
       n = i
       call smd_util_rebox(n,n,lat,rlat,rd)
c
c      compute histogram
c      -----------------
c      size of the bin
       dr = rmax/nb
       do j=1,n
         d=sqrt(rd(j,1)*rd(j,1)+rd(j,2)*rd(j,2)+rd(j,3)*rd(j,3))
         k = int(d/dr)+1
         if(k.le.nb ) then
          gr(k) = gr(k) + 1
         end if
       end do
c
       do j=1,nb
         write(55,*) gr(j)
       end do
c
       fourpi = 16.0*atan(1.0)
       const = fourpi*real(n1)*real(n2)/(vol*3.0d0)
       do k=1,nb
         rl = real(k-1)*dr
         ru = rl+dr
         norm = const*(ru**3-rl**3)
         gr(k) = gr(k)/norm
       end do
       do j=1,nb
         write(56,*) gr(j)
       end do
       return
c      -------------
c      error section
c      -------------
136    continue
       write(*,*) "subroutine:" //pname(1:len_trim(pname))
       write(*,*) "error number of atoms field"
       write(*,*) "line number: ", l
       write(*,*) "buffer: ",bigbuf
       stop
       return
137    continue
       write(*,*) "error reading coords field"
       write(*,*) "line number: ", l
       write(*,*) "buffer: ", bigbuf
       stop
       return
138    continue
       write(*,*) message
       stop
c
c      format statements
c      -----------------
30     FORMAT(180A1)
20     FORMAT(A180)
       end

       subroutine pair_compute(rmax,lat,i1,atag,nb,infile)
       implicit none
       double precision lat(3,3)
       double precision rmax
       integer i1,nb
       character*(*) atag
       character*(*) infile
c
       integer j
       integer i,n,a,nd,k,n0
       double precision rlat(3,3),r1(3)
       double precision, dimension(:,:), allocatable :: c
       double precision, dimension(:,:), allocatable :: rd
       character*5, dimension(:), allocatable :: tag
       integer, dimension(:), allocatable :: gr
       double precision, dimension(:), allocatable :: grs
       character*30 buf
       character*180 bigbuf
       double precision d,dr
       double precision fourpi,x
       integer  nf 
       double precision vol,sum,norm,const,rl,ru,rho
       integer l,sl
 
c
c
       nf = 0
c
c      size of the bin
       dr = rmax/nb
c      g(r) array
       allocate(gr(nb))
       allocate(grs(nb))
       gr = 0
       grs = 0.0d0
c      get reciprocal lattice vectors 
       call smd_lat_invrt(lat,rlat)
       write(*,*) "Inverse Lattice vectors"
       write(*,*) (rlat(i,1),i=1,3)
       write(*,*) (rlat(i,2),i=1,3)
       write(*,*) (rlat(i,3),i=1,3)
c      get lattice volume
       call smd_latt_vol(lat,vol)
c
c      allocate other memory
       open(10,file=infile,
     $            form='formatted',status='old')
       read(10,*) n
       rewind(10)
       allocate(c(n,3))
       allocate(rd(n,3))
       allocate(tag(n))
c      start calculating distances
       l = 0
       n0 = 0
10     continue
c
c      get number of atoms (skipping empty lines)
c      ------------------------------------------
11     continue
       read(10,20,end=135) bigbuf
c      increment line number
       l=l+1
       if(bigbuf.eq."") goto 11
       read(bigbuf,*,err=136) n
c      check if number of atoms consistent with prior
       if(n0.eq.0) then
         n0=n
       else if(n0.ne.n) then
         goto 136
       end if
c
c      get title field(if any) 
c      -------------------
       read(10,20) bigbuf
       l=l+1
c      processing coodinates
       j = 0
       do i=1,n
c
c       read coordinates(skipping empty lines)
c       ------------------------------------
12      continue
        read(10,20) bigbuf
        l=l+1
        if(bigbuf.eq."") goto 12
        read(bigbuf,*,err=137) tag(i),(c(i,k),k=1,3)
c
c       extract coords of secondary atoms
c       with "atag" name (i1 is a central atom)
c       ---------------------------------------
c        if((tag(i).eq.atag).and.i.ne.i1) then
        sl = len_trim(atag)
        if(index(tag(i),atag(1:sl)).ne.0) then
        if(i.ne.i1) then
          j = j+1
          rd(j,:) = c(i,:)
        end if
        end if
       end do
       
c
c      construct relative distances
c      ----------------------------
       nd = j 
       r1 = c(i1,:)
       do j=1,nd
         rd(j,:) = rd(j,:)-r1
       end do

c
c      fold the distances into the box
c      -------------------------------
       call smd_util_rebox(nd,n,lat,rlat,rd)

c
c      compute histogram
c      -----------------
       do j=1,nd
         d=sqrt(rd(j,1)*rd(j,1)+rd(j,2)*rd(j,2)+rd(j,3)*rd(j,3))
         k = int(d/dr)+1
         if(k.le.nb ) then
          gr(k) = gr(k) + 1
         end if
       end do
       nf = nf + 1

       goto 10
      
135    continue
        write(*,*) "Number of frames processed",nf
c
c      smoothing g(r) per Allan/Tild. p 204 (6.48)
c      ------------------------------------------
       grs(1)=(69.0*gr(1)+4.0*gr(2)-6.0*gr(3)+
     >         4.0*gr(4)-gr(5))/70.0d0
       grs(2)=(2.0*gr(1)+27.0*gr(2)+12.0*gr(3)
     >  -8.0*gr(4)+2.0*gr(5))/35.0d0
       do i=3,nb-2
         grs(i) =(-3.0*gr(i-2)+12.0*gr(i-1)+17.0*gr(i)+
     >            12.0*gr(i+1)-3.0*gr(i+2))/35.0d0 
       end do
       grs(nb-1)=(2.0*gr(nb)+27*gr(nb-1)+12.0*gr(nb-2)-
     >            8.0*gr(nb-3)+2.0*gr(nb-4))/35.0d0
       grs(nb)=(69.0*gr(nb)+4.0*gr(nb-1)-6.0*gr(nb-2)+
     >          4.0*gr(nb-3)-gr(nb-4))/70.0d0
       open(12,file="gr.dat",
     $            form='formatted',status='unknown')

       fourpi = 16.0*atan(1.0)
       rho = real(nd)/vol
       x = 1.0/real(nf)
       const = fourpi*rho/3.0d0
       sum = 0.0d0
       do k=1,nb
         rl = real(k-1)*dr
         ru = rl+dr
         norm = const*(ru**3-rl**3)
         sum = sum + grs(k)*x
c         write(12,'(i5,F12.6,I5,F12.6,F12.6)') 
c     >             k,ru,gr(k),sum,real(gr(k))*x/norm
         write(12,'(3F12.6)') 
     >             ru,real(grs(k))*x/norm,sum
       end do
       return
c      error section
c      -------------
136    continue
       write(*,*) "error number of atoms field"
       write(*,*) "line number: ", l
       write(*,*) "buffer: ",bigbuf
       return
137    continue
       write(*,*) "error reading coords field"
       write(*,*) "line number: ", l
       write(*,*) "buffer: ", bigbuf
       return
c
c      format statements
c      -----------------
30     FORMAT(180A1)
20     FORMAT(A180)
       end

      subroutine smd_util_rebox(n,nmax,latt,rlatt,aaa)

      implicit none

      integer n,nmax
      double precision rlatt(3,3),latt(3,3)
      double precision  aaa(nmax,3)
c
      integer i
      double precision  ssx,ssy,ssz,xss,yss,zss
      logical oprint


      if(n.eq.1) then
       oprint =.true.
      else
       oprint = .false.
      end if
      oprint = .false.
      do i=1,n

       if(oprint)
     >          write(*,*) "rebox",aaa(i,1),aaa(i,2),aaa(i,3)
       ssx=(rlatt(1,1)*aaa(i,1)+rlatt(1,2)*aaa(i,2)+rlatt(1,3)*aaa(i,3))
       ssy=(rlatt(2,1)*aaa(i,1)+rlatt(2,2)*aaa(i,2)+rlatt(2,3)*aaa(i,3))
       ssz=(rlatt(3,1)*aaa(i,1)+rlatt(3,2)*aaa(i,2)+rlatt(3,3)*aaa(i,3))

       xss=ssx-nint(ssx)
       yss=ssy-nint(ssy)
       zss=ssz-nint(ssz)

       aaa(i,1)=(latt(1,1)*xss+latt(1,2)*yss+latt(1,3)*zss)
       aaa(i,2)=(latt(2,1)*xss+latt(2,2)*yss+latt(2,3)*zss)
       aaa(i,3)=(latt(3,1)*xss+latt(3,2)*yss+latt(3,3)*zss)

      enddo

      return

      END

      subroutine smd_latt_vol(latt,vol)
      implicit none
      real*8 x,y,z,latt,vol

      dimension latt(3,3)

      x=latt(2,2)*latt(3,3)-latt(2,3)*latt(2,3)
      y=latt(3,2)*latt(1,3)-latt(1,2)*latt(3,3)
      z=latt(1,2)*latt(2,3)-latt(2,2)*latt(1,3)

      vol=abs(latt(1,1)*x+latt(2,1)*y+latt(3,1)*z)

      return

      END

      subroutine smd_lat_invrt(latt,rlatt)
      implicit none
      double precision  latt(3,3),rlatt(3,3)
c
      double precision  det

      rlatt(1,1)=latt(2,2)*latt(3,3)-latt(3,2)*latt(2,3)
      rlatt(2,1)=latt(3,1)*latt(2,3)-latt(2,1)*latt(3,3)
      rlatt(3,1)=latt(2,1)*latt(3,2)-latt(3,1)*latt(2,2)
      rlatt(1,2)=latt(3,2)*latt(1,3)-latt(1,2)*latt(3,3)
      rlatt(2,2)=latt(1,1)*latt(3,3)-latt(3,1)*latt(1,3)
      rlatt(3,2)=latt(3,1)*latt(1,2)-latt(1,1)*latt(3,2)
      rlatt(1,3)=latt(1,2)*latt(2,3)-latt(2,2)*latt(1,3)
      rlatt(2,3)=latt(2,1)*latt(1,3)-latt(1,1)*latt(2,3)
      rlatt(3,3)=latt(1,1)*latt(2,2)-latt(2,1)*latt(1,2)

      det=latt(1,1)*rlatt(1,1)+latt(1,2)*rlatt(2,1)+latt(1,3)*rlatt(3,1)
      if(abs(det).gt.0.d0)det=1.d0/det

      rlatt(1,1)=det*rlatt(1,1)
      rlatt(2,1)=det*rlatt(2,1)
      rlatt(3,1)=det*rlatt(3,1)
      rlatt(1,2)=det*rlatt(1,2)
      rlatt(2,2)=det*rlatt(2,2)
      rlatt(3,2)=det*rlatt(3,2)
      rlatt(1,3)=det*rlatt(1,3)
      rlatt(2,3)=det*rlatt(2,3)
      rlatt(3,3)=det*rlatt(3,3)

      return
      end

      function is_number(string)
      implicit none
      character*(*) string
      logical is_number
c
      character(len=11) :: anumber = "0123456789."
      integer i,l
      is_number = .false.
      if(string.eq." ") return
      l = len_trim(string)
      do i=1,l
       if(string(i:i).ne." ") exit
      end do
      is_number = verify(string(i:l), anumber).eq.0
      end function

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
       write(13,*) i,string(i:i)
       if(string(i:i).ne." ") exit
      end do
      is_integer = verify(string(i:l), anumber).eq.0
      end function
