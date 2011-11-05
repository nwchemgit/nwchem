       program transform_pdb
       implicit none
       character*255 infile
       character*255 outfile
       character*16  tar
       integer i,n
       logical ofile
       
       integer, dimension(:), allocatable :: ir1,ir2

       inquire(file="swap_pdb.dat",exist=ofile)
       if(ofile) then
          open(10,file="swap_pdb.dat",
     $            form='formatted',status='unknown')
          read(10,*) infile
          read(10,*) outfile
          read(10,*) n
          allocate(ir1(n))
          allocate(ir2(n))
          do i=1,n
            read(10,*) ir1(i), ir2(i)
          end do
          close(10)
       else
          write(*,*) "please provide swap_pdb.dat file"
          stop
c       write(6,*) "Enter input PDB file"
c       read(5,*) infile
c       write(6,*) "Enter output PDB file"
c       read(5,*) outfile
c       write(6,*) "Enter ir1" 
c       read(5,*) ir1
c       write(6,*) "Enter ir2" 
c       read(5,*) ir2
c       write(6,*) "If sorting by atom name enter it now"
c       read(5,*) tar
       end if 
     
       call swap(infile,outfile,n,ir1,ir2)

       deallocate(ir2)
       deallocate(ir1)

       end

       subroutine swap(infile,outfile,in,ir1,ir2)
       implicit none
       character*(*) infile
       character*(*) outfile
       integer in
       integer ir1(in),ir2(in)
c
       integer ntot,nres
       double precision, dimension(:,:), allocatable :: c
       integer, dimension(:), allocatable :: ir
       character*16 , dimension(:), allocatable :: ta,tr
       integer, dimension(:), allocatable :: is
       integer, dimension(:), allocatable :: im
       integer i,j
       double precision ctmp(3)
 
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
       allocate(is(nres))
       allocate(im(nres+1))
       call smd_pdb_sequence_bounds(ntot,nres,ir,is,im)
       do i=1,ntot
        do j=1,in
        if(ir(i).eq.ir1(j)) then
          ir(i)=ir2(j)
        else if(ir(i).eq.ir2(j)) then
          ir(i)=ir1(j)
        end if
        end do
       end do
       call smd_pdb_sort_byres(ntot,ta,tr,ir,c)
c
       call smd_pdb_write_byseq(outfile,nres,ntot,im,is,ta,tr,c)
!
!      == deallocate arrays ==
       deallocate(im)
       deallocate(is)
       deallocate(ir)
       deallocate(tr)
       deallocate(ta)
       deallocate(c)
       end


c $Id$
