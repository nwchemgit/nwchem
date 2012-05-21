       program transform_pdb
       implicit none
       character*255 infile
       character*255 outfile
       character*16  tar
       logical ofile
       integer i,l,istatus 
       character*72 buffer
       integer ntot
c
       i=1
       call my_get_command_argument(i,buffer,l,istatus)
       if(istatus.ne.0) goto 18
       read(buffer,*) infile
18     continue
       write(*,*) "input file is",trim(infile)
       call cluster(infile)
c
       end

       subroutine cluster(infile)
       implicit none
       character*(*) infile
c
       integer ntot,nres
       double precision, dimension(:,:), allocatable :: c
       double precision, dimension(:,:), allocatable :: cg
       double precision, dimension(:), allocatable :: cd
       integer, dimension(:), allocatable :: ir
       integer, dimension(:), allocatable :: is
       integer, dimension(:), allocatable :: im
       integer, dimension(:,:), allocatable :: p2
       integer, dimension(:,:), allocatable :: p3
       integer, dimension(:,:), allocatable :: p4
       integer, dimension(:,:), allocatable :: p5
       integer, dimension(:,:), allocatable :: p6
       integer, dimension(:,:), allocatable :: p7
       integer, dimension(:,:), allocatable :: p8
       integer, dimension(:,:), allocatable :: p9
       character*16 , dimension(:), allocatable :: ta,tr
       integer i,il
       integer n2
       integer n3
       integer n4
       integer n5
       integer n6
       integer n7
       integer n8
       integer n9
 
       integer nmax
       nmax = 500
c
       call smd_pdb_natoms(infile,ntot)
       allocate(c(3,ntot))
       allocate(ta(ntot))
       allocate(tr(ntot))
       allocate(ir(ntot))
c
       call smd_pdb_read_coords(infile,ntot,c)
       call smd_pdb_read_atomres(infile,ntot,ta,tr,ir)
       call smd_pdb_sort_byres(ntot,ta,tr,ir,c)
       call smd_pdb_nres0(nres,ntot,ir)
       allocate(cg(3,nres))
       call smd_pdb_cog(ntot,nres,ir,c,cg)
       n2 = nres*(nres-1)/2
       allocate(p2(2,n2))
       call parse_dimers(nres,cg,n2,p2)
       allocate(p3(3,nmax))
       allocate(p4(4,nmax))
       allocate(p5(5,nmax))
       allocate(p6(6,nmax))
       allocate(p7(7,nmax))
       allocate(p8(8,nmax))
       allocate(p9(9,nmax))
       n3 = nmax
       n4 = nmax
       n5 = nmax
       n6 = nmax
       n7 = nmax
       n8 = nmax
       n9 = nmax

       call parse_trimers(n2,p2,n2,2,p2,n3,p3)
c       do i=1,n3
c         write(*,*) (p3(il,i),il=1,3)
c       end do
       call parse_trimers(n2,p2,n3,3,p3,n4,p4)
       do i=1,n4
         write(*,'(100I4)') (p4(il,i),il=1,4)
       end do
       call parse_trimers(n2,p2,n4,4,p4,n5,p5)
       do i=1,n5
         write(*,'(100I4)') (p5(il,i),il=1,5)
       end do
       call parse_trimers(n2,p2,n5,5,p5,n6,p6)
       do i=1,n6
         write(*,'(100I4)') (p6(il,i),il=1,6)
       end do
!!!      == deallocate arrays ==
       deallocate(p2)
       deallocate(cg)
       deallocate(ir)
       deallocate(tr)
       deallocate(ta)
       deallocate(c)
       end

      subroutine parse_dimers(nres,cg,nb,p2)
      implicit none
      integer nres
      double precision cg(3,nres)
      integer nb
      integer p2(2,nb)
c 
      integer i,j
      integer k
      double precision rd
      double precision dist
      double precision r1(3),r2(3)
      external dist


      k = 0
      do i=1,nres
      r1 = cg(:,i)
      do j=i+1,nres
        r2 = cg(:,j)
        rd  = dist(r1,r2)
        if(rd.lt.3.0) then
          k = k + 1
          p2(1,k) = i
          p2(2,k) = j
          write(12,*) i,j,rd
        end if 
      end do
      end do
      nb = k
      end subroutine

      subroutine parse_trimers(nb2,p2,nb,l,p,nx,px)
      implicit none
      integer nb2,nb,l,nx
      integer p2(2,nb2)
      integer p(l,nb)
      integer px(l+1,nx)
c 
      integer ptmp(l+1)
      integer ptmp0(l+1)
      integer i,j
      integer k
      integer tail, head
      logical ohead, otail
      integer ic

      integer i2,i3,j3,ik
      integer fn_out
      logical ofound,omatch

      fn_out = 78
c       open(fn_out,file="trimer.dat",
c     $            form='formatted',status='unknown')

      open(fn_out, form='formatted', status='scratch')
      k = 0
      do i2=1,nb2
      head = p2(1,i2)
      tail = p2(2,i2)
      do i=1,nb
        ic = 0
        do j=1,l
           ptmp(j) = p(j,i)
           if(p(j,i).eq.head) then
             ic = ic + 1
             ptmp(l+1) = tail
           end if
           if(p(j,i).eq.tail) then
             ic = ic + 1
             ptmp(l+1) = head
           end if
        end do
        if(ic.eq.1) then
          call sort(l+1,ptmp)
          if(k.ne.0) rewind(fn_out)
          ofound = .false.
          do ik=1,k
            read(fn_out,*) ptmp0
            if(ALL(ptmp0.eq.ptmp)) ofound = .true.
          end do
          if(.not.ofound) then
             write(fn_out,*) ptmp
             k = k +1
             px(:,k) = ptmp
          end if
        end if
      end do
      end do
      nx = k
      close(fn_out)
      end subroutine

      subroutine parse_fourmers(nb2,p2,nb,l,p,nx,px)
      implicit none
      integer nb2,nb,l,nx
      integer p2(2,nb2)
      integer p(l,nb)
      integer px(l+1,nx)
c 
      integer ptmp(l+1)
      integer ptmp0(l+1)
      integer i,j
      integer k
      integer tail, head
      logical ohead, otail
      integer ic

      integer i2,i3,j3,ik
      integer fn_out
      logical ofound,omatch
      character*73 message

      message = "opening file"
      fn_out = 33
       open(unit=fn_out,file="fourimer.dat",
     $            form='formatted',status='unknown',err=911)

      k = 0
      do i2=1,nb2
      head = p2(1,i2)
      tail = p2(2,i2)
      do i=1,nb
        ic = 0
        do j=1,l
           ptmp(j) = p(j,i)
           if(p(j,i).eq.head) then
             ic = ic + 1
             ptmp(l+1) = tail
           end if
           if(p(j,i).eq.tail) then
             ic = ic + 1
             ptmp(l+1) = head
           end if
        end do
        if(ic.eq.1) then
          call sort(l+1,ptmp)
          if(k.ne.0) rewind(fn_out)
          ofound = .false.
          do ik=1,k
            write(*,*) "ik=",ik
            read(fn_out,*) ptmp0
            write(*,*) "ptmp0=",ptmp0
            if(ALL(ptmp0.eq.ptmp)) ofound = .true.
          end do
          if(.not.ofound) then
             k = k +1
             write(*,*) k,ptmp
             write(fn_out,*) ptmp
             px(:,k) = ptmp
          end if
        end if
      end do
      end do
      close(fn_out)
      nx = k
      return
911    continue       
c      if you reach this you are in trouble
       write(*,*) "Emergency STOP "
       write(*,*) trim(message)
       stop
       end subroutine
      subroutine sort(n,a)
      implicit none
      integer n
      integer a(n)
c
c     local variables:
      integer i
      integer pass
      integer sorted
      integer temp
      character*32 pname

      pass = 1
      sorted = 0
      do while(sorted .eq. 0)
        sorted = 1
        do 2 i = 1,n-pass
          if(a(i) .gt. a(i+1)) then
            temp = a(i)
            a(i) = a(i+1)
            a(i+1) = temp
            sorted = 0
          endif
 2      continue
        pass = pass +1
      end do
      do i=1,n-1
       if(a(i).eq.a(i+1)) a(i)=-1
      end do

      return

      end

      function dist(r1,r2)
      implicit none
      double precision r1(3),r2(3)
      double precision dist
      dist = (r1(1)-r2(1))**2+
     +       (r1(2)-r2(2))**2+
     +       (r1(3)-r2(3))**2

      dist = sqrt(dist)
      return
      end function

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

      subroutine parse_trimers0(nb2,p2,nb,l,p,px)
      implicit none
      integer nb2,nb,l
      integer p2(2,nb2)
      integer px(l+1,nb)
      integer p(l,nb)
c 
      integer ptmp(l+1)
      integer i,j
      integer k
      integer tail, head
      logical ohead, otail
      integer ic

      integer i2,i3,j3

      k = 0
      do i2=1,nb2
      head = p2(1,i2)
      tail = p2(2,i2)
      do i=1,nb
        ic = 0
        do j=1,l
           ptmp(j) = p(j,i)
           if(p(j,i).eq.head) then
             ic = ic + 1
             ptmp(l+1) = tail
           end if
           if(p(j,i).eq.tail) then
             ic = ic + 1
             ptmp(l+1) = head
           end if
        end do
        if(ic.eq.1) then
           k = k+1
c          px(1:l,k) = p(:,i)
c          px(l+1,k) = tail
c         call sort(l+1,ptmp)
          write(14,*) (p(j,i),j=1,l),head,tail
          write(15,*) (ptmp(j),j=1,l+1)
        end if
        do j=1,l+1
          write(*,*) j,k,nb
          px(j,k) = ptmp(j)
        end do
      end do
      end do
      end subroutine

      subroutine parse_trimers1(nb2,p2,nb,l,p)
      implicit none
      integer nb2,nb,l
      integer p2(2,nb2)
      integer p(l,nb)
c 
      integer ptmp0(l+1)
      integer ptmp(l+1)
      integer i,j
      integer tail, head
      logical ohead, otail
      integer ic

      integer ik,in,il,jn,k,h0,h1,kk
      character*256 buf0,buf1
      logical ofound,omatch
      integer fn_out

      fn_out = 78
       open(fn_out,file="trimer.dat",
     $            form='formatted',status='unknown')
      kk = 0
      buf0 = " "
      ptmp0 = 0
      do in=1,nb
      do il=1,l
         h0 = p(il,in)
         h1 = p(min(il+1,l),in)
         do jn=1,nb2
           omatch = .false.
           if(p2(1,jn).eq.h0) then
             if(p2(2,jn).ne.h1) then
                     omatch = .true.
             end if
           else if(p2(2,jn).eq.h0) then
            omatch = .true.
           end if
           if(omatch) then
                write(14,*) (p(k,in),k=1,l),p2(1,jn),p2(2,jn)
                ptmp(2:l+1)= p(:,in)
                ptmp(1) = p2(2,jn)
                call sort(l+1,ptmp)
                write(buf1,*) (ptmp(k),k=1,l+1)
                ofound = .false.
                if(kk.ne.0) rewind(fn_out)
                do ik=1,kk
                  read(fn_out,*) ptmp0
                  if(ALL(ptmp0.eq.ptmp)) ofound = .true.
                end do
                if(.not.ofound) then
                   write(fn_out,*) ptmp
                   kk = kk +1
                end if
           end if
         end do
      end do
      end do
      rewind(fn_out)
      read(fn_out,'(a)') buf0
      write(*,*) "buf0",buf0
      end subroutine
c $Id$
