       program transform_pdb
       implicit none
       character*255 infile
       character*255 outfile
       character*16  tar
       logical ofile
       
c
       inquire(file="transform_pdb.dat",exist=ofile)
       if(ofile) then
          open(10,file="transform_pdb.dat",
     $            form='formatted',status='unknown')
          read(10,*) infile
          read(10,*) outfile
          read(10,*) tar
          close(10)
       else
       write(6,*) "Enter input PDB file"
       read(5,*) infile
       write(6,*) "Enter output PDB file"
       read(5,*) outfile
       write(6,*) "If sorting by atom name enter it now"
       read(5,*) tar
       end if 
     
       call reorder(tar,infile,outfile)

       end

       subroutine reorder(aname,infile,outfile)
       implicit none
       character*(*) infile
       character*(*) outfile
       character*(*) aname
c
       integer ntot,nres
       double precision, dimension(:,:), allocatable :: c
       double precision, dimension(:,:), allocatable :: cg
       double precision, dimension(:), allocatable :: cd
       integer, dimension(:), allocatable :: ir
       integer, dimension(:), allocatable :: is
       integer, dimension(:), allocatable :: im
       character*16 , dimension(:), allocatable :: ta,tr
       integer i
 
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
       allocate(cd(nres))
       allocate(is(nres))
       allocate(im(nres+1))
c
       if(aname.eq."*") then
         call smd_pdb_cog(ntot,nres,ir,c,cg)
       else
         call smd_pdb_cog_byname(ntot,nres,aname,ta,ir,c,cg)
       end if
       call smd_pdb_sequence_bounds(ntot,nres,ir,is,im)
       call smd_pdb_distance(nres,cg,cg(:,1),cd)
       call smd_pdb_sort_seq_distance(nres,is,im,cd)
       call smd_pdb_write_byseq(outfile,nres,ntot,im,is,ta,tr,c)
!
!      == deallocate arrays ==
       deallocate(ir)
       deallocate(tr)
       deallocate(ta)
       deallocate(c)
       deallocate(cg)
       deallocate(cd)
9000   FORMAT("ATOM",T7,I5,T13,A4,T18,A3,T23,I4,T31,
     >        F8.3,T39,F8.3,T47,F8.3,T55,F6.2)
       end


c $Id$
