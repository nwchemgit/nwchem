       program rotate
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
       integer iat,ires,ires0
       integer nrs
       integer nat
       integer i
       character*5, dimension(:), allocatable :: atag
       double precision, dimension(:,:), allocatable :: c
       double precision, dimension(:), allocatable :: q
       integer id(100),nid
       double precision angle,pi,x(3),v(3),v1(3)
c
       external util_get_io_unit
       logical util_get_io_unit
c      ---------------------
c      read input parameters
c      ---------------------
       pi = 4.0*atan(1.0)
       inquire(file="rotate.dat",exist=ofile)
       if(ofile) then
          open(10,file="rotate.dat",
     $            form='formatted',status='unknown')
          read(10,*,ERR=911) solutefile
          read(10,*,ERR=911) outfile
          read(10,*,ERR=911) angle
          do i=1,1000
            read(10,*,ERR=911,END=30) id(i)
            nid = i
          end do
          close(10)
       else
          goto 911
       end if 
30     continue
       write(*,*) solutefile
       write(*,*) outfile
       write(*,*) "index ",(id(i),i=1,nid)
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
       write(*,*) "nat",nat
       allocate(c(3,nat))
       allocate(atag(nat))
       rewind(fns)
       call xyz_read(nat,c,atag,fns)
       write(*,*) "nat",nat
       do i=1,nid-2
         v  = c(:,id(nid-1)) 
         v1 = c(:,id(nid))-v
         x = c(:,id(i))
         write(*,*) "rotating by",angle
         angle = angle*pi/180.0
         write(*,*) "rotating by",angle
         call rot_around_axis(v(1),v(2),v(3),v1(1),v1(2),v1(3),
     >                        x(1),x(2),x(3),angle)
c         call
c     >   RotatePointAroundVector(x(1),x(2),x(3),
c     >                           v,v1,angle)
         c(:,id(i))=x
       end do
       call xyz_write(nat,c,atag,fno)
       angle = 62
       angle = angle*pi/180.0
       v(1) = -4
       v(2) = 0
       v(3) = -5
       v1(1) = 3
       v1(2) = 2
       v1(3) = -5
       x(1) = 4
       x(2) = 0
       x(3) = -5
         call rot_around_axis(v(1),v(2),v(3),v1(1),v1(2),v1(3),
     >                        x(1),x(2),x(3),angle)
         write(*,*) "x=",x(1),x(2),x(3)
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
        write(*,*) "nrm=",nrm
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
