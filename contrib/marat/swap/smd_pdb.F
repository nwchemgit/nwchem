      subroutine smd_pdb_natoms(filename,nt)
      implicit none
      character*(*) filename
      integer nt
c
      character*(4) buffer
      integer un
c
      un = 70
      open(unit=un,status="old",form="formatted",file=filename)

      nt = 0
100   continue
      read(un,'(A4)',end=200) buffer
      if(buffer(1:4).eq."ATOM") then
        nt = nt + 1
      end if
      goto 100
200   continue
      close(un)

      end

      subroutine smd_pdb_read_coords(filename,nt,c)
      implicit none
      character*(*) filename
      integer nt
      double precision c(3,nt)
c
      character*(180) buffer
      character*(4)  tag
      integer i
      integer un
      character*(30) pname

      pname = "smd_pdb_read_coords"
c
      un = 70
c
      open(unit=un,status="old",form="formatted",file=filename)
      i = 0
100   continue
      read(un,'(A180)',end=200) buffer
      if(buffer(1:4).eq."ATOM") then
        i = i +1
        read(buffer,*) tag,tag,tag,tag,tag,
     >                 c(1,i),c(2,i),c(3,i)
       end if
      goto 100
200   continue
      close(un)

      return

      end

      subroutine smd_pdb_read_atomres(filename,nt,ta,tr,ir)
      implicit none
      character*(*) filename
      integer nt
      character*(*) ta(nt)
      character*(*) tr(nt)
      integer ir(nt)
c
      character*(180) buffer
      integer i
      integer un
      character*(30) pname

      pname = "smd_pdb_read_atomres"
c
      un = 70
c
      open(unit=un,status="old",form="formatted",file=filename)
      i = 0
100   continue
      read(un,'(A180)',end=200) buffer
      if(buffer(1:4).eq."ATOM") then
        i = i +1
        ta(i) = buffer(13:16)
        tr(i) = buffer(18:20)
        read(buffer(23:26),*) ir(i)
       end if
      goto 100
200   continue
      close(un)

      return

      end

      subroutine smd_pdb_sort_byres(nt,ta,tr,ir,c)
      implicit none
      integer nt
      character*(*) ta(nt)
      character*(*) tr(nt)
      integer ir(nt)
      double precision c(3,nt)
c
      character*(180) buffer
      integer i
      integer un
      character*(30) pname
      integer pass
      integer sorted
      integer itemp
      double precision ftemp
      character*16 stemp

      pass = 1
      sorted = 0
      do while(sorted .eq. 0) 
        sorted = 1
        do i = 1,nt-pass
          if(ir(i) .gt. ir(i+1)) then
            itemp = ir(i)
            ir(i) = ir(i+1)
            ir(i+1) = itemp
c
            stemp = ta(i)
            ta(i) = ta(i+1)
            ta(i+1) = stemp
c
            stemp = tr(i)
            tr(i) = tr(i+1)
            tr(i+1) = stemp
c
            ftemp = c(1,i)
            c(1,i) = c(1,i+1)
            c(1,i+1) = ftemp
c
            ftemp = c(2,i)
            c(2,i) = c(2,i+1)
            c(2,i+1) = ftemp
c
            ftemp = c(3,i)
            c(3,i) = c(3,i+1)
            c(3,i+1) = ftemp
c
            sorted = 0
          endif
        enddo
        pass = pass +1
      end do

      return

      end

      subroutine smd_pdb_sort_seq_distance(nr,is,im,d)
      implicit none
      integer nr
      integer is(nr)
      integer im(nr)
      double precision d(nr)
c
      integer i
      integer pass
      integer sorted
      integer itemp
      double precision ftemp

      pass = 1
      sorted = 0
      do while(sorted .eq. 0) 
        sorted = 1
        do i = 1,nr-pass
          if(d(i) .gt. d(i+1)) then
            ftemp = d(i)
            d(i) = d(i+1)
            d(i+1) = ftemp
c
            itemp = is(i)
            is(i) = is(i+1)
            is(i+1) = itemp
c
c            itemp = im(i)
c            im(i) = im(i+1)
c            im(i+1) = itemp

            sorted = 0
          endif
        enddo
        pass = pass +1
      end do

      return

      end

      subroutine smd_pdb_read(filename,nt,ta,tr,ir,c)
      implicit none
      character*(*) filename
      integer nt
      character*(16) ta(nt)
      character*(16) tr(nt)
      integer ir(nt)
      double precision c(3,nt)
c
      character*(180) buffer
      character*(4)  tag
      integer i
      integer un
      character*(30) pname

      pname = "smd_pdb_read"
c
      un = 70
c
      open(unit=un,status="old",form="formatted",file=filename)
      i = 0
100   continue
      read(un,'(A180)',end=200) buffer
      if(buffer(1:4).eq."ATOM") then
        i = i +1
        read(buffer,*) tag,tag,ta(i),tr(i),ir(i),
     >                 c(1,i),c(2,i),c(3,i)
       end if
      goto 100
200   continue
      close(un)

      return

      end

      subroutine smd_pdb_write_byseq(filename,nr,nt,im,is,ta,tr,c)
      implicit none
      character*(*) filename
      integer nr
      integer im(nr+1)
      integer is(nr)
      integer nt
      character*(16) ta(nt)
      character*(16) tr(nt)
      double precision c(3,nt)
c
      character*(180) buffer
      character*(4)  tag
      integer i
      integer un
      integer i0,j0,j,jb,je
c
      un = 70
c
      open(unit=un,status="unknown",form="formatted",file=filename)
      j0=0
      do i=1,nr
       i0=is(i)
       jb = im(i0)
       je = im(i0+1)-1
       do j=jb,je
         j0=j0+1
         write(un,FMT=9000)j0,ta(j),tr(j),i,c(1,j),c(2,j),c(3,j)
       end do
      end do
      write(un,'(A3)') "END"
      close(un)

      return
9000   FORMAT("ATOM",T7,I5,T13,A4,T18,A3,T23,I4,T31,
     >        F8.3,T39,F8.3,T47,F8.3,T55,F6.2)

      end

      subroutine smd_pdb_nres(filename,nr)
      implicit none
      character*(*) filename
      integer nr
c
      character*(180) buffer
      character*(4)  tag
      integer ir0,ir
      integer un
c
      un = 70
c
      open(unit=un,status="old",form="formatted",file=filename)

c      reset residue arrays to be the size of number of residues only
      nr = 0
      ir0 = 0
100   continue
      read(un,'(A180)',end=200) buffer
      if(buffer(1:4).eq."ATOM") then
        read(buffer,*) tag,tag,tag,tag,ir
        if(ir0.ne.ir) then
          nr = nr + 1
          ir0=ir
        end if
      end if
      goto 100
200   continue
      close(un)

      end

      subroutine smd_pdb_nres0(nr,nt,ir)
      implicit none
      integer nr
      integer nt
      integer ir(nt)
c
      integer i,ir0
      nr = 0
      ir0 = 0
      do i=1,nt
      if(ir0.ne.ir(i)) then
        nr = nr + 1
        ir0=ir(i)
      end if
      end do

      end

      subroutine smd_pdb_cog(nt,nr,ir,c,cg)
      implicit none
      integer nt
      integer nr
      integer ir(nt)
      double precision c(3,nt)
      double precision  cg(3,nr)
c
      integer i,j,ir0,nm
      j=1
      ir0=ir(1)
      nm=0
      cg=0.0d0
      do i=1,nt
       if(ir(i).ne.ir0) then
         cg(:,j)=cg(:,j)/nm
         j=j+1
         ir0=ir(i)
         nm=0
       end if 
       cg(:,j)=cg(:,j)+c(:,i)
       nm=nm+1
      end do
      cg(:,j)=cg(:,j)/nm
      return
      end

      subroutine smd_pdb_cog_byname(nt,nr,tar,ta,ir,c,cg)
C     nt total number of atoms
C     nr total number of residues
C     tar atom name
C     ta  atom name array
C     ir  residue index array
C     c   coordinate array
C     cg  center of mass array
      implicit none
      integer nt
      integer nr
      character*(*) tar
      character*(*) ta(nt)
      integer ir(nt)
      double precision c(3,nt)
      double precision  cg(3,nr)

c
      integer i,j,ir0,nm
      integer s0
      integer length
      external length

      s0 = length(tar)
      j=1
C     first residue
      ir0=ir(1)
      nm=0
      cg=0.0d0
      do i=1,nt
       if(ir(i).ne.ir0) then
         if(nm.ne.0) then
           cg(:,j)=cg(:,j)/nm
         else
           write(*,*) "no atoms maching pattern ",tar(1:s0)
           stop
         end if
         j=j+1
         ir0=ir(i)
         nm=0
       end if 
       write(*,*) "atom name",tar
C      check for name match
       if(INDEX(ta(i),tar(1:s0)).gt.0) then
       cg(:,j)=cg(:,j)+c(:,i)
       nm=nm+1
       end if
      end do
      cg(:,j)=cg(:,j)/nm
      return
      end

      subroutine smd_pdb_sequence_bounds(nt,nr,ir,is,im)
      implicit none
      integer nt
      integer nr
      integer ir(nt)
      integer is(nr)
      integer im(nr+1)
c
      integer i,j,ir0,nm
      j = 1
      is(1) = ir(1)
      im(1) = 1
      do i=1,nt
       if(ir(i).ne.is(j)) then
         j=j+1
         is(j) = ir(i)
         im(j) = i
       end if 
      end do
      im(j+1)=i


      return
      end

      subroutine smd_pdb_distance(nt,c,cor,cd)
      implicit none
      integer nt
      double precision c(3,nt)
      double precision  cor(3)
      double precision cd(nt)
c
      double precision c1(3)
      integer i

      do i=1,nt
        c1(:)=c(:,i)-cor(:)
        c1(:)=c1(:)*c1(:)
        cd(i)=SUM(c1)
      end do
 
      cd = sqrt(cd)
      return
      end

      subroutine smd_pdb_read_res(filename,nt,nr,tr,ir,nm)
      implicit none
      character*(*) filename
      integer nt,nr,nc
      character*(16) tr(nr)
      integer ir(nt)
      integer nm(nr)
c
      character*(30) pname
      character*(180) buffer
      character*(4)  tag
      character*(16)  rtag,rtag0
      integer ir0,nr0
      integer ncenter
      integer un
c
      pname = "sg_read_res"
c
      un = 70
c
      open(unit=un,status="old",form="formatted",file=filename)

      ncenter = 0
      nr0 = 0
      rtag0 = " "
      ir0 = 0
100   continue
      read(un,'(A180)',end=200) buffer
      if(buffer(1:4).eq."ATOM") then
        ncenter = ncenter + 1
        read(buffer,*) tag,tag,tag,rtag,ir(ncenter)
        if(ir0.ne.ir(ncenter)) then
          ir0=ir(ncenter)
          nr0 = nr0 + 1
          tr(nr0) = rtag
          rtag0=rtag 
        end if
        ir(ncenter) = nr0
        nm(nr0) = nm(nr0) + 1
      end if
      goto 100
200   continue

      close(un)

      return
      end

      function length(string) 
      implicit none
      integer length
*returns length of string ignoring trailing blanks 
      character*(*) string 
      integer i
      do i = len(string), 1, -1 
         if(string(i:i) .ne. ' ') go to 20 
      end do
20    length = i 
      end 

c $Id$
