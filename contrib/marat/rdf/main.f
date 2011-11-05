       program pair_correlation
       implicit none
       integer i,l
       double precision lat(3,3)
       double precision rmax
       integer i1,nb
       character*255 infile
       character*5 atag
       logical ofile
       
c
       inquire(file="pair.dat",exist=ofile)
       if(ofile) then
          open(10,file="pair.dat",
     $            form='formatted',status='old')
c         max distance
          read(10,*) rmax
c         lattice vectors a,b,c
          read(10,*) (lat(i,1),i=1,3)
          read(10,*) (lat(i,2),i=1,3)
          read(10,*) (lat(i,3),i=1,3)
c         index of central ion
          read(10,*) i1
c         tag of secondary ions
          read(10,*) atag
c         number of bins
          read(10,*) nb
c         xyz trajectory file name
          read(10,*) infile
          close(10)
       else
       write(6,*) "Please provide pair.dat file in the following format"
       write(*,*) "Maximum Distance:"
       write(*,*) "Lattice vectors"
       write(*,*) "Central ion index:"
       write(*,*) "Secondary ion tag:"
       write(*,*) "Number of bins:"
       write(*,*) "XYZ input file: "

       stop
       end if 
c      check if xyzfile does exist
       inquire(file=infile,exist=ofile)     
       if(.not.ofile) then
         write(6,*) "Cannot find xyz input file: ", infile
         stop
       end if
       write(*,*) "Maximum Distance:",rmax 
       write(*,*) "Lattice vectors"
       write(*,*) (lat(i,1),i=1,3)
       write(*,*) (lat(i,2),i=1,3)
       write(*,*) (lat(i,3),i=1,3)
       write(*,*) "Index of central ion:",i1
       write(*,*) "Secondary ion tag:",atag
       write(*,*) "Number of bins:",nb
       write(*,*) "XYZ input file: ",infile
c
       call pair_compute(rmax,lat,i1,atag,nb,infile)
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
       open(22,file="tmp.xyz",
     $            form='formatted',status='unknown')
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
c $Id$
