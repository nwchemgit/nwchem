       program rotate
       implicit none
       character*255 message
       character*255 solutefile
       character*255 cosmofile
       character*255 outfile
       character*255 infile
       character*255 buffer
       character*16  tar
       character*4 a4
       character*255 buf
       logical ofile
       logical overb
       integer fns,fnc,fno      
       integer iat,ires,ires0
       integer nrs
       integer nat
       integer i
       character*5, dimension(:), allocatable :: atag
       double precision, dimension(:,:), allocatable :: c
       double precision, dimension(:), allocatable :: q
       integer id(100),nid,ida(2)
       double precision angle,pi,x(3),v(3),v1(3)
       integer istatus,k,l
c
       external util_get_io_unit
       logical util_get_io_unit
       external is_integer
       logical is_integer
       logical is_number
       external is_number
c      --------------------------------------------      
c      beging parsing command line arguments if any
c      --------------------------------------------      
       overb = .false.
       infile = " "
       outfile = " "
       nid = 0
       i = 0
       ida = 0
       angle = 1000
16     continue
       i = i+1
       call my_get_command_argument(i,buffer,l,istatus)
       if(istatus.ne.0) goto 18
       if(buffer.eq."-atoms") then
          do 
            i = i+1
            call my_get_command_argument(i,buffer,l,istatus)
            if(istatus.ne.0) goto 18
            if(is_integer(buffer)) then
              nid = nid+1
              read(buffer,*) id(nid)
            else
              i = i-1 
              exit
            end if
          end do
          go to 16
       else if(buffer.eq."-v") then
          overb = .true.
          goto 16
       else if(buffer.eq."-help") then
          write(*,1000)
          stop
       else if(buffer.eq."-angle") then
          i = i+1
          call my_get_command_argument(i,buffer,l,istatus)
          if(is_number(buffer)) then
            read(buffer,*) angle
            goto 16
          else
            goto 911
          end if
       else if(buffer.eq."-axis") then
          do k=1,2
            i = i+1
            call my_get_command_argument(i,buffer,l,istatus)
            if(istatus.ne.0) goto 18
            if(is_integer(buffer)) then
              read(buffer,*) ida(k)
            else
              message = "two atom ids are required to define the axis"
              goto 911
            end if
           end do
          go to 16
       else 
          if(infile.eq." ") then
              infile = buffer
              goto 16
          else if (outfile.eq." ") then
              outfile = buffer
              goto 18
          end if
       end if 
c       ---------------------------      
c      end of command line parsing
c      ---------------------------      
18     continue
       if(nid.eq.0) then
          message = "no atoms to rotate"
          goto 911
       end if
       if(ANY(ida.eq.0)) then
          message = "missing complete axis specification"
          goto 911
       end if
       if(angle.ge.1000) then
          message = "mising angle specification"
          goto 911
       end if
       if(infile.eq." ") then
          message = "missing input file"
          goto 911
       end if
       if(outfile.eq." ") then
          message = "missing output file"
          goto 911
       end if
       if(ABS(angle).gt.360) then
          message = "angle should be between -360 and 360 degrees"
          goto 911
       end if
       if(overb) then
          write(*,*) "atom id ",(id(k),k=1,nid)
          write(*,*) "axis id ",(ida(k),k=1,2)
          write(*,*) "angle",angle
          write(*,*) "infile ",trim(infile)
          write(*,*) "outfile ",trim(outfile)
       end if
c      
c      ---------------------
c      read input parameters
c      ---------------------
       pi = 4.0*atan(1.0)
c       inquire(file="rotate.dat",exist=ofile)
c       if(ofile) then
c          open(10,file="rotate.dat",
c     $            form='formatted',status='unknown')
c          read(10,*,ERR=911) infile
c          read(10,*,ERR=911) outfile
c          read(10,*,ERR=911) angle
c          do i=1,1000
c            read(10,*,ERR=911,END=30) id(i)
c            nid = i
c          end do
c          close(10)
c       else
c          goto 911
c       end if 
30     continue
c      -------------
c      open io units
c      -------------
       if(.not.util_get_io_unit(fns)) then
           message = "no free file units"
           goto 911
       else
           message = "opening "//trim(infile)
           open(fns,file=trim(infile),
     $                form='formatted',status='old',err=911)
       end if 
       if(.not.util_get_io_unit(fno)) then
           message = "no free file units"
           goto 911
       else
           message = "opening "//trim(outfile)
           open(fno,file=trim(outfile),
     $                form='formatted',status='unknown',err=911)
       end if 
c      -------------------
c      read xyz input file
c      -------------------
       call xyz_read_natoms(nat,fns)
       if(overb) write(*,*) "total number of atoms",nat
       allocate(c(3,nat))
       allocate(atag(nat))
       rewind(fns)
       call xyz_read(nat,c,atag,fns)
       angle = angle*pi/180.0
       v  = c(:,ida(1)) 
       v1 = c(:,ida(2))-v
       if(overb) write(*,*) "rotating by",angle
       do i=1,nid
         x = c(:,id(i))
         call rot_around_axis(v(1),v(2),v(3),v1(1),v1(2),v1(3),
     >                        x(1),x(2),x(3),angle)
         c(:,id(i))=x
       end do
       call xyz_write(nat,c,atag,fno)
       if(overb) write(*,*) "succesfull completion"
       stop
c      if you reach this you are in trouble
1000   format(
     > "NAME:",/
     > " rotate rotate atoms ",
     > " around specified axis               ",//,
     >  "SYNOPSIS",/,
     > " rdf [-help] ",
     > " [-atoms ] [-axis ] [-angle ] ",
     > " input file(s) outpiut file ",//,
     >  "DESCRIPTION",/,
     > " This utility rotates selected atoms around the axis",
     > " from the xyz files",//,
     > " The required option are as follows:",/,
     > " -atoms   atom ids to be rotated",/,
     > " -axis   pair of atom ids used as the axis",/,
     > " -angle rotation angle",/,
     > " other options are as follows:",/,
     > " -help  prints out this message",/,
     > " -v  generate more verbose output"
     > )
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

       subroutine xyz_read(n0,c,atag,fn_in)
       implicit none
       integer n0
       double precision c(3,n0)
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
        read(bigbuf,*,err=911) atag(i),(c(k,i),k=1,3)
        n0=i
       end do
c      check if all atoms were read
       if(n0.lt.n) then
         message = "could not find all the atoms"
         goto 911
       end if
20     continue
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

       subroutine xyz_write(n0,c,atag,fn_in)
       implicit none
       integer n0
       double precision c(3,n0)
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
       pname = "xyz_write"
       write(fn_in,*) n0
       write(fn_in,*) "rotated"
c      
       do i=1,n0
        write(fn_in,*) atag(i),(c(k,i),k=1,3)
       end do
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

       subroutine RotatePointAroundVector(x,y,z,p1,p2,angle)
       implicit none
       double precision p1(3),p2(3)
       double precision x,y,z,u,v,w,angle
       double precision ux,uy,uz
       double precision vy,vx,vz
       double precision wx,wy,wz
       double precision sa,ca
       double precision d,d2
       double precision u2,v2,w2,f1,bv,cw,au
       double precision a,b,c

       u = p1(1)-p2(1)
       v = p1(2)-p2(2)
       w = p1(3)-p2(3)
       a = p2(1)
       b = p2(2)
       c = p2(3)
       ux=u*x
       uy=u*y
       uz=u*z
       vx=v*x
       vy=v*y
       vz=v*z
       wx=w*x
       wy=w*y
       wz=w*z
       u2 = u*u
       v2 = v*v
       w2 = w*w
       sa=sin(angle)
       ca=cos(angle)
       d2 = u*u+v*v+w*w
       d = sqrt(d2)
       write(*,*) "before rotation",x,y,z
       f1 = ux+vy+wz
       bv = b*v
       cw=c*w
       au=a*u
       x = (a*(v2+w2)-u*(bv+cw-f1))*(1-ca)+d2*x*ca+
     >      d*(-c*v+b*w-w*y+v*z)*sa
       x = x/d2
       z = (b*(u2+w2)-v*(au+cw-f1))*(1-ca)+d2*y*ca+
     >       d*(c*u -a*w+w*x-u*z)*sa
       y = y/d2
       z = (c*(u2+v2)-w*(au+bv-f1))*(1-ca)+d2*z*ca+
     >        d*(-b*u+a*v-v*x+u*y)*sa
       z = z/d2
       write(*,*) "after rotation",x,y,z
       end subroutine

       subroutine RotatePointAroundVector1(x,y,z,u,v,w,a)
       implicit none
       double precision x,y,z,u,v,w,a
       double precision ux,uy,uz
       double precision vy,vx,vz
       double precision wx,wy,wz
       double precision sa,ca
       ux=u*x
       uy=u*y
       uz=u*z
       vx=v*x
       vy=v*y
       vz=v*z
       wx=w*x
       wy=w*y
       wz=w*z
       sa=sin(a)
       ca=cos(a)
       x=u*(ux+vy+wz)+(x*(v*v+w*w)-u*(vy+wz))*ca+(-wy+vz)*sa
       y=v*(ux+vy+wz)+(y*(u*u+w*w)-v*(ux+wz))*ca+(wx-uz)*sa
       z=w*(ux+vy+wz)+(z*(u*u+v*v)-w*(ux+vy))*ca+(-vx+uy)*sa
       end subroutine

       subroutine rot_around_axis(a,b,c,u,v,w,x,y,z,theta)
        
        double precision a,b,c,u,v,w,x,y,z,theta
        double precision u2,v2,w2
        double precision cost, onemcost,sint
        double precision nrm
        double precision x1,y1,z1
        nrm = u*u+v*v+w*w
        nrm = sqrt(nrm)
        u = u/nrm
        v = v/nrm
        w = w/nrm
        u2 = u*u
        v2 = v*v
        w2 = w*w
        cost = cos(theta)
        onemcost = 1 - cost
        sint = sin(theta)

        x1 = (a*(v2 + w2) - u*(b*v + c*w - u*x - v*y - w*z)) * onemcost
     +      + x*cost
     +      + (-c*v + b*w - w*y + v*z)*sint;

        y1 = (b*(u2 + w2) - v*(a*u + c*w - u*x - v*y - w*z)) * onemcost
     +      + y*cost
     +      + (c*u - a*w + w*x - u*z)*sint;

        z1= (c*(u2 + v2) - w*(a*u + b*v - u*x - v*y - w*z)) * onemcost
     +      + z*cost
     +      + (-b*u + a*v - v*x + u*y)*sint;

         x = x1
         y = y1
         z = z1
       end subroutine

      subroutine my_get_command_argument(i,buffer,l,istatus)
      implicit none
      integer i
      integer l
      character*(*) buffer
      integer istatus
c
      istatus = 0
      l = 0
      call getarg(i,buffer)
      if(buffer.eq." ") istatus = 1
c      call get_command_argument(i,buffer,l,istatus)
      end subroutine

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


c $Id$
