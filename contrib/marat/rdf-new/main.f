c $Id$

       program pair_correlation
       implicit none
       integer i,l
       double precision lat(3,3),latv(3),lat1(3,3)
       double precision rmax
       integer i1,j1,nb
       character*5 tvr
c       character*5 atag
       logical ofile
       CHARACTER(LEN=11) :: anumber = "0123456789."
       logical is_integer
       external is_integer
       logical is_number
       external is_number
       integer  my_command_argument_count
       external my_command_argument_count
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
       logical ok
       logical osmooth
       double precision f
       integer fn_in
       integer fn_out
       integer n,n1,n2
       double precision dr
       integer nf
       integer nfrm
       integer naver
       integer nu,nvm,nva,nv,nprec
       integer ifil
c      allocatable arrays
       integer nc
       integer nfil
       logical, dimension(:), allocatable :: oc1,oc2
       double precision, dimension(:,:), allocatable :: c
       double precision, dimension(:,:), allocatable :: c1
       double precision, dimension(:,:), allocatable :: c2
       double precision, dimension(:), allocatable :: gr0
       double precision, dimension(:), allocatable :: gr
       double precision, dimension(:), allocatable :: occ,occ0
       character*5, dimension(:), allocatable :: atag
       character*5 , dimension(:), allocatable :: tva,tua,tur
       integer, dimension(:), allocatable :: iur
       character*180 , dimension(:), allocatable :: infile
       logical util_get_io_unit
       external util_get_io_unit
c
c      --------------------------------------------      
c      set some default values
c      --------------------------------------------      
       fn_in = 10
       fn_out = 21
       nfrm = 0
       aformat = " "
       overb = .false.
       file_lattice = " "
       latv = -1.0
       lat = 0.0
       atom1_id = 0
       atom2_id = 0
       atom1_tag = " "
       atom2_tag = " "
       rmax = -1
       nb = -1
       osmooth = .false.
c
c      ------------------------
c      allocate io file array
c      ------------------------
       nc = my_command_argument_count()
       call get_carg_nfiles(nc)
       allocate(infile(nc))
       infile = " "
c      --------------------------------------------      
c      beging parsing command line arguments if any
c      --------------------------------------------      
       i = 0
       nfil = 0
16     continue
       i = i+1
       call my_get_command_argument(i,buffer,l,istatus)
       if(istatus.ne.0) goto 18
       if(buffer.eq."-lattice") then
          do k=1,3
            i = i+1
            call my_get_command_argument(i,buffer,l,istatus)
            if(istatus.ne.0) goto 18
            if(is_number(buffer)) then
              read(buffer,*) f
              latv(k:3) = f
            else
              if(k.gt.1) then
                i = i-1 
                goto 16
              else if(buffer(1:1).eq."-") then
                i = i-1
                goto 16
              else
                file_lattice = buffer
                exit
              end if
            end if
          end do
          go to 16
       else if(buffer.eq."-v") then
          overb=.true.
          go to 16
       else if(buffer.eq."-help") then
          ohelp=.true.
          write(*,1000)
          stop
       else if(buffer.eq."-atom1") then
          i = i+1
          call my_get_command_argument(i,buffer,l,istatus)
          if(istatus.ne.0) goto 18
          if(is_integer(buffer)) then
            read(buffer,*) atom1_id
          else
            atom1_tag = buffer
          end if
          go to 16
       else if(buffer.eq."-atom2") then
          i = i+1
          call my_get_command_argument(i,buffer,l,istatus)
          if(istatus.ne.0) goto 18
          if(is_number(buffer)) then
            read(buffer,*) atom2_id
          else
            atom2_tag = buffer
          end if
          go to 16
       else if(buffer.eq."-smooth") then
          osmooth = .true.
          goto 16
       else if(buffer.eq."-rmax") then
          i = i+1
          call my_get_command_argument(i,buffer,l,istatus)
          if(istatus.ne.0) goto 18
          if(is_number(buffer)) then
            read(buffer,*) rmax
          else
            message = "Maximum radius gas to be a number"
            goto 911
          end if
          go to 16
        else if(buffer.eq."-nframes") then
          i = i+1
          call my_get_command_argument(i,buffer,l,istatus)
          if(istatus.ne.0) goto 18
          if(is_integer(buffer)) then
            read(buffer,*) nfrm
          else
            message = "Number of frames have to be integer"
            goto 911
          end if
          go to 16
        else if(buffer.eq."-nbins") then
          i = i+1
          call my_get_command_argument(i,buffer,l,istatus)
          if(istatus.ne.0) goto 18
          if(is_integer(buffer)) then
            read(buffer,*) nb
          else
            message = "Number of bins have to be integer"
            goto 911
          end if
          go to 16
       else if(buffer(1:1).eq."-") then
          message = "Unknown keyword: "//trim(buffer)
          goto 911
       else 
          nfil = nfil+1
          infile(nfil) = buffer
          go to 16
       end if 
       
c      ---------------------------      
c      end of command line parsing
c      ---------------------------      
18     continue
c      ---------------------------      
c      start checks/balances
c      ---------------------------      
c      files
       if(nfil.lt.2) then
         message = "please provide input and output files"
         goto 911
       end if
c      atom tags
       if(atom1_id.eq.0.and.atom1_tag.eq." ") then
         message = "please provide central atom id or tag"
         goto 911
       end if
c      the lattice
       if(file_lattice.ne." ") then
         if(overb) write(*,*) "reading lattice from
     +       file"//trim(file_lattice)
         call lattice_read_file(file_lattice,lat)
       else
         do k=1,3
           lat(k,k) = latv(k)
         end do
       end if
c      maximum radius
       if(rmax.lt.0) then
         message = "please provide rmax"
         goto 911
       end if
c      number of bins
       if(nb.lt.0) then
         nb = INT(rmax/0.1)
       end if
c      output file
       file_out = infile(nfil)
       message = "opening output file "//file_out
       if(.not.util_get_io_unit(fn_out)) then
           message = "no free file units"
           goto 911
       end if 
       open(fn_out,file=file_out,
     $            form='formatted',status='unknown',err=911)

       if(overb) then
         write(*,*) "lattice vectors"
         write(*,*) (lat(i,1),i=1,3)
         write(*,*) (lat(i,2),i=1,3)
         write(*,*) (lat(i,3),i=1,3)
         if(atom1_id.ne.0) then
           write(*,*) "1st atom ID: ",atom1_id
         else
           write(*,*) "1st atom tag: ",atom1_tag
         end if
         if(atom2_id.ne.0) then
           write(*,*) "2nd atom ID: ",atom2_id
         else
           write(*,*) "2nd atom tag: ",atom2_tag
         end if
         write(*,*) "Number of bins:",nb
         write(*,*) "Maximum radius:",rmax
         write(*,*) "Input file(s)"
         do i=1,nfil-1
           write(*,*) trim(infile(i))
         end do
         write(*,*) "Output file:"
         write(*,*) trim(infile(nfil))
       end if

       allocate(gr(nb))
       allocate(gr0(nb))
       allocate(occ(nb))
       allocate(occ0(nb))
       naver = 0
       gr = 0
       occ = 0
c      figure out format for trajectory file
       do ifil=1,nfil-1
       file_in = infile(ifil)
c      open input channel
       inquire(file=file_in,exist=ofile)
       if(.not.util_get_io_unit(fn_in)) then
           message = "no free file units"
           goto 911
       end if 
       if(ofile) then
          open(fn_in,file=file_in,
     $            form='formatted',status='old',err=911)
       else
           message = "no file found: "//file_in
           goto 911
       end if 

       i=INDEX(file_in,".",.true.) 
       aformat=file_in(i+1:)

       if(aformat.eq."trj") then
         call trj_read_header(fn_in,nvm,nva,nu,nprec)
         allocate(tva(nva))
         allocate(tua(nu))
         allocate(tur(nu))
         allocate(iur(nu))
         call trj_read_solvent_specs(fn_in,nva,tvr,tva)
         call trj_read_solute_specs(fn_in,nu,tur,tua,iur)
         n = nva*nvm + nu
       else if(aformat.eq."xyz") then
         call xyz_read_natoms(n,fn_in)
       end if
       rewind(fn_in)
       allocate(c(3,n))
       allocate(atag(n))
       allocate(oc1(n))
       allocate(oc2(n))
c       
       if(aformat.eq."trj") then
         i=0
         do i1=1,nvm
         do j1=1,nva
           i=i+1
           atag(i) = tva(j1)
         end do
         end do
         do i1=1,nu
           i=i+1
           atag(i) = tua(i1)
         end do
       end if
c      loop over frames
       nf = 0
       do
         if(aformat.eq."trj") then
           call trj_read(fn_in,n,c,lat1,ok)
         else if(aformat.eq."xyz")  then
           call xyz_read(n,c,atag,fn_in)
           ok=n.ne.0
           lat1=lat
         end if
c        check lattice
         if(any(lat1.lt.0)) then
            message = "problems with lattice input"
            goto 911
         end if
         if(.not.ok) exit
         oc1=.false.
         call mask_all(n,c,atag,oc1,atom1_id,atom1_tag)
         oc2=.false.
         call mask_all(n,c,atag,oc2,atom2_id,atom2_tag)
c        if possible make atom1 and atom2 distinct
         oc2=oc2.and.(.not.oc1)
         if(count(oc2).eq.0) oc2=oc1
c
         if(any(oc2.and.oc1)) oc2=oc2.and.(.not.oc1)
         gr0 = 0
         call rdf_compute2(n,c,oc1,oc2,
     +                    nb,rmax,lat1,
     +                    gr0,
     +                    occ0,
     +                    fn_in)
         gr = gr + gr0
         occ = occ + occ0
         nf = nf+1
         if(nfrm.gt.0.and.nf.ge.nfrm) then
            if(overb) 
     +      write(*,*) "exiting because exceeded number of frames",nf,nfrm
            exit
         end if
       end do
       close(fn_in)
       naver = naver + nf
c      ---------------------------
c      end looping over the frames
c      ---------------------------
31     continue
       deallocate(oc2)
       deallocate(oc1)
       deallocate(atag)
       deallocate(c)
       
       if(aformat.eq."trj") then
         deallocate(iur)
         deallocate(tur)
         deallocate(tua)
         deallocate(tva)
       end if
c      end of loop of input files
       end do
c       
c      average RDF over the frames
       if(overb) 
     +   write(*,*) "averaging over ",naver," smaples"
       gr = gr/real(naver)
       occ = occ/real(naver)
c
c      smoothing g(r) per Allan/Tild. p 204 (6.48)
c      ------------------------------------------
       if(osmooth) then
       if(overb) 
     +   write(*,*) "smoothing rdf"
       gr0=0.0d0
       gr0(1)=(69.0*gr(1)+4.0*gr(2)-6.0*gr(3)+
     >         4.0*gr(4)-gr(5))/70.0d0
       gr0(2)=(2.0*gr(1)+27.0*gr(2)+12.0*gr(3)
     >  -8.0*gr(4)+2.0*gr(5))/35.0d0
       do i=3,nb-2
         gr0(i) =(-3.0*gr(i-2)+12.0*gr(i-1)+17.0*gr(i)+
     >            12.0*gr(i+1)-3.0*gr(i+2))/35.0d0 
       end do
       gr0(nb-1)=(2.0*gr(nb)+27*gr(nb-1)+12.0*gr(nb-2)-
     >            8.0*gr(nb-3)+2.0*gr(nb-4))/35.0d0
       gr0(nb)=(69.0*gr(nb)+4.0*gr(nb-1)-6.0*gr(nb-2)+
     >          4.0*gr(nb-3)-gr(nb-4))/70.0d0
       gr = gr0
       end if


       dr = rmax/nb
       do k=1,nb
         write(fn_out,'(3F12.6)') 
     >             real(k)*dr,gr(k),occ(k)
       end do
       stop
411    continue
1000   format(
     > "NAME:",/
     > " rdf constructs radial density function ",
     > " from NWChem trajectory or xyz files",//,
     >  "SYNOPSIS",/,
     > " rdf [-help] [-v ]",
     > " [-frames ] [-lattice ] [-nbins ] [-rmax] [-atom1] [-atom2]",
     > " input file(s) outpiut file ",//,
     >  "DESCRIPTION",/,
     > " This utility builds rdf from NWChem trajectory",
     > " or xyz files",//,
     > " The required option are as follows:",/,
     > " -atom1   first atom id or name",/,
     > " -atom2   second atom id or name",/,
     > " -lattice lattice file or vector(3)",/,
     > " -rmax    maximum distance ",//,
     > " other options are as follows:",/,
     > " -help  prints out this message",/,
     > " -v  generate more verbose output",/,
     > " -frames  maximum number of frames default is all ",/,
     > " -nbins number of bins with default of rmax/0.1"
     > )
 911    continue       
c      if you reach this you are in trouble
       write(*,*) "Emergency STOP"
       write(*,*) message
       stop
       end program

      subroutine trj_read(fn,n,c,lat,ok)
      implicit none
      integer, intent(in)    :: fn  
      integer, intent(in)   :: n
      logical, intent(out)  :: ok
      double precision, intent(out)   :: c(3,n)
      double precision, intent(out)   :: lat(3,3)
c      
      character*80 card
      character*80 message
      character*30 token
      character*1 sep
      character*1 a1
      integer i0
      integer i,k,lb,lt
      logical ostatus
c
      ok = .false.
      sep = " "
c      
c     look for the start of the fame record
c     -------------------------------------
  100 continue
      read(fn,1000,end=11,err=911) card
 1000 format(a)
      if(card(1:5).ne.'frame') goto 100
c
c     skip one lines here
c     --------------------
      read(fn,1000,err=911,end=911) card
c     read lattice
      message = "reading lattice"
      do i=1,3
        read(fn,1000,err=911,end=911) card
        i0=1
        do k=1,3
          message = message(1:len_trim(message)+1)//token
          call get_next_token(i0,card,sep,token,ostatus)
          read(token,*,ERR=911,END=911) lat(i,k)
        end do
      end do
      lat = lat*10.0d0
      read(fn,1000,err=911,end=911) card
      do i=1,n
        read(fn,1000,err=911,end=911) card
        i0=1
        do k=1,3
          call get_next_token(i0,card,sep,token,ostatus)
          read(token,*,ERR=911,END=911) c(k,i)
        end do
      end do
      c = c*10.0d0
      ok = .true.
11    continue
      return
911    continue       
c      if you reach this you are in trouble
       write(*,*) "found error"
       write(*,*) message
       stop

      end subroutine

      subroutine trj_read_header(fn,nvm,nva,nu,nprec)
      implicit none
      integer, intent(in)    :: fn  
      integer, intent(out)   :: nvm
      integer, intent(out)   :: nva
      integer, intent(out)   :: nu
      integer, intent(out)   :: nprec
c      
      character*80 card
      character*80 message
      character*30 token
      character*1 sep
      character*1 a1
      integer i0
      integer i,k,lb,lt
      logical ostatus
c
      rewind(fn)
c      
c     look for the "header"
c     ---------------------
  100 continue
      read(fn,1000,end=911) card
c      write(*,*) card
 1000 format(a)
      if(card(1:6).ne.'header') goto 100
c     parse the line      
      read(fn,1000,err=911,end=911) card
      i0 = 1
      sep = " "
      call get_next_token(i0,card,sep,token,ostatus)
      read(token,*) nva
      call get_next_token(i0,card,sep,token,ostatus)
      read(token,*) nu
      call get_next_token(i0,card,sep,token,ostatus)
      call get_next_token(i0,card,sep,token,ostatus)
      call get_next_token(i0,card,sep,token,ostatus)
      read(token,*) nvm
      call get_next_token(i0,card,sep,token,ostatus)
      read(token,*) nprec
      return
911    continue       
c      if you reach this you are in trouble
       write(*,*) "found error"
       write(*,*) message
       stop

      end subroutine

      subroutine get_next_token(i0,buffer,sep,token,ostatus)
      implicit none
      integer ,         intent(inout) :: i0
      character(LEN=*), intent(in)    :: buffer
      character(LEN=*), intent(out)   :: token
      character*1,      intent(in)    :: sep
      logical,          intent(out)   :: ostatus
c      
      character*1 a1
      integer i,k,lb,lt
c      
      token = " "
      ostatus = .false.
      lb=len(buffer)
      lt=len(token)
      k=0
      do i=i0,lb
       a1=buffer(i:i)
       if(a1.eq." ") then
          if(k.gt.0.and.sep.eq." ") goto 30
       else if(a1.eq.sep.and.sep.ne." ") then
          goto 30
       else
         k = k+1
         if(k.gt.lt) goto 30
         token(k:k) = a1
       end if
      end do
30    continue
      i0=i+1
      ostatus = k.gt.0
      end subroutine

       subroutine lattice_read_file(fname,lat)
       implicit none
       character*(*) fname
       double precision lat(3,3)
c       
       integer ifn
       logical ofile
       character*180 bigbuf
c       
       character*180 message
       character*30 pname
       integer i,k
       logical is_integer
       external is_integer
       logical util_get_io_unit
       external util_get_io_unit
c
       pname = "lattice_read_file"
       lat = 0.0d0
       inquire(file=fname,exist=ofile)
       if(.not.util_get_io_unit(ifn)) then
           message = "no free file units"
           goto 911
       end if 
       if(ofile) then
               open(ifn,file=fname,
     $            form='formatted',status='old',err=911)
       else
           message = "no file found: "//fname
           goto 911
       end if 
c 
       i = 0
       do 
         message = "reading lattice file"
         read(ifn,20,end=10,err=911) bigbuf
c        increment line number
         if(bigbuf.eq." ") cycle
         i = i+1
         read(bigbuf,*,err=10) (lat(k,i),k=1,3)
         if(i.eq.3) exit
       end do
       close(ifn)
10     continue
       if(i.eq.0.or.i.eq.2) then
         message = "something wrong with format"
         goto 911
       end if
       if(i.eq.1) then
          lat(2,2) = lat(2,1)
          lat(3,3) = lat(3,1)
          lat(2,1) = 0.00d0
          lat(3,1) = 0.00d0
       end if
       return
c
911    continue       
c      if you reach this you are in trouble
       write(*,*) "Emergency STOP in"
       write(*,*) pname(1:len_trim(pname))
       write(*,*) message(1:len_trim(message))
       stop
c
c      format statements
c      -----------------
30     FORMAT(180A1)
20     FORMAT(A180)
       end

       subroutine compute_default_rmax(rmax,lat)
       implicit none
       double precision rmax,lat(3,3)
c       
       double precision d(3)

       d(1) = sqrt(SUM(lat(:,1)*lat(:,1)))
       d(2) = sqrt(SUM(lat(:,2)*lat(:,2)))
       d(3) = sqrt(SUM(lat(:,3)*lat(:,3)))
       rmax = 0.5*MINVAL(d)
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

       subroutine mask_all(n,c,atag,oc,mid,mtag)
       implicit none
       integer n
       double precision c(3,n)
       character*(*) atag(n)
       logical oc(n)
       integer mid
       character*(*) mtag
c
       integer i
       if(mid.ne.0) then
         oc(mid)=.true.
        else if (mtag.ne. " ") then
c         oc=atag(i).eq.mtag(1:len_trim(mtag))
         oc=index(atag,mtag(1:len_trim(mtag))).ne.0
       end if
       end

       subroutine rdf_compute2(n,c,oc1,oc2,
     +                        nb,rmax,lat,
     +                        gr,
     +                        occ,
     +                        fn_in)
       implicit none
       integer n
       double precision c(3,n)
       logical oc1(n),oc2(n)
       integer fn_in
       integer nb
       double precision rmax
       double precision gr(nb)
       double precision occ(nb)
       double precision lat(3,3)
c       
       integer i1,i2
       double precision rd(3)
       double precision rd1(3)
c
       integer j
       integer i,a,nd,k
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
       double precision  sx,sy,sz,xs,ys,zs
       double precision  ss(3),rr(3)
       integer ic
       double precision ang,d1,ang0,ang1
       logical ocons,ofile
       integer fnc
       logical util_get_io_unit
       external util_get_io_unit
       double precision pi
       double precision asum
c
       pname = "rdf_compute"
       pi = acos(-1.0)
c
c      ----------------------
c      compute lattice params
c      ----------------------
       call smd_lat_invrt(lat,rlat)
       call smd_latt_vol(lat,vol)
c
c      any constraints?
       ic = 0
       ang0=0.0d0
       ang1=0.0d0
       inquire(file="cons.dat",exist=ofile)
       if(ofile) then
c         write(*,*) "found constraint file"
         if(.not.util_get_io_unit(fnc)) then
             message = "no free file units"
             goto 911
         end if 

          open(fnc,file="cons.dat",
     $            form='formatted',status='old',err=911)
          read(fnc,*,err=911) ic,ang0,ang1
          close(fnc)
          ang0 = pi*ang0/180.0
          ang1= pi*ang1/180.0
          ang0=cos(ang0)
          ang1=cos(ang1)
       end if 
c      ----------------------------
c      construct relative distances
c      and bin rdf
c      ----------------------------
       dr = rmax/nb
       i=0
       do i1=1,n
         if(oc1(i1)) then
         do i2=1,n
           if(oc2(i2)) then
           i = i+1
           rd =  c(:,i2)-c(:,i1)

           ss = matmul(rlat,rd)
           rr = ss-nint(ss)
           rd = matmul(lat,rr)

           d = sqrt(sum(rd*rd))
           ocons =.true.
           if(ic.gt.0) then
             rd1 = c(:,ic)-c(:,i1)
             d1= sqrt(sum(rd1*rd1))
             ang = SUM(rd*rd1)/(d1*d)
             ocons = ang.ge.ang0
             ocons = ocons.and.(ang.le.ang1)
           end if
           if(d.gt.0.01) k = int(d/dr)+1
           if(k.le.nb.and.ocons ) then
                   gr(k) = gr(k) + 1
           end if
           end if
         end do
         end if
       end do
c       do j=1,nb
c         write(55,*) gr(j)
c       end do
       i1 = count(oc1)
       i2 = count(oc2)
c
       fourpi = 16.0*atan(1.0)
       const = fourpi*real(i1*i2)/(vol*3.0d0)
       asum = 0.0d0
       do k=1,nb
         asum   = asum+gr(k)
         occ(k) = asum
       end do
       occ = occ/i1
       do k=1,nb
         rl = real(k-1)*dr
         ru = rl+dr
         norm = const*(ru**3-rl**3)
         gr(k) = gr(k)/norm
       end do
       return
c      -------------
c      error section
c      -------------
911    continue
       write(*,*) "subroutine:" //pname(1:len_trim(pname))
       stop
       return
       end

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

      function my_command_argument_count()
      implicit none
      integer my_command_argument_count
c
      my_command_argument_count = iargc()
c      my_command_argument_count = command_argument_count() 
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
       if(string(i:i).ne." ") exit
      end do
      is_integer = verify(string(i:l), anumber).eq.0
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

      subroutine get_carg_files(n,buffer,nb,istatus)
      implicit none
      integer n
      integer nb
      character*(*) buffer(n)
      integer istatus
c
      integer i,nc
      character*1 buf1
      nc = iargc()
      n = 0
      do i=nc,1,-1
        n = n+1
        call getarg(i,buffer(n))
        buf1 = adjustl(buffer(n))
        if(buf1.eq."-") then
          buffer(n) = " "
          n = n-1
          exit
        end if
      end do
      end subroutine

      subroutine get_carg_nfiles(n)
      implicit none
      integer n
      integer istatus
c
      integer i,nc
      character*1 a1
      nc = iargc()
      n = 0
      do i=nc,1,-1
        n = n+1
        call getarg(i,a1)
        if(a1.eq."-") then
          n = n-1
          exit
        end if
      end do
      end subroutine

      subroutine trj_read_solvent_specs(fn,nva,tvr,tva)
      implicit none
      integer, intent(in)    :: fn  
      integer, intent(in)   :: nva
      character*5, intent(out) :: tvr
      character*5, intent(out) :: tva(nva)
c      
      character*80 card
      character*80 message
      character*30 token
      character*1 sep
      character*1 a1
      integer i0
      integer i,k,lb,lt
      logical ostatus
c
      rewind(fn)
c      
c     look for the "header"
c     ---------------------
  100 continue
      read(fn,1000,end=911) card
 1000 format(a)
      if(card(1:6).ne.'header') goto 100
c     skip one line      
      read(fn,1000,err=911,end=911) card
c      
      sep = " "
      do i=1,nva
        read(fn,1000,err=911,end=911) card
        i0 = 1
        call get_next_token(i0,card,sep,token,ostatus)
        read(token,*) tvr
        call get_next_token(i0,card,sep,token,ostatus)
        read(token,*) tva(i)
      end do
c      do i=1,nva
c        write(*,*) tvr,tva(i)
c      end do
      return
911    continue       
c      if you reach this you are in trouble
       write(*,*) "found error"
       write(*,*) message
       stop

      end subroutine

      subroutine trj_read_solute_specs(fn,n,tur,tua,iur)
      implicit none
      integer, intent(in)    :: fn  
      integer, intent(in)   :: n
      character*5, intent(out) :: tur(n)
      character*5, intent(out) :: tua(n)
      integer, intent(out) :: iur(n)
c      
      integer nu,nva
      character*80 card
      character*80 message
      character*30 token
      character*1 sep
      character*1 a1
      integer i0
      integer i,k,lb,lt
      logical ostatus
c
      rewind(fn)
      sep = " "
c      
c     look for the "header"
c     ---------------------
  100 continue
      read(fn,1000,end=911) card
 1000 format(a)
      if(card(1:6).ne.'header') goto 100
c     find how many lines to skip 
      read(fn,1000,err=911,end=911) card
      i0 = 1
      call get_next_token(i0,card,sep,token,ostatus)
      read(token,*) nva
      call get_next_token(i0,card,sep,token,ostatus)
      read(token,*) nu
      if(nu.gt.n) then
        message = "insufficient solute size "
        goto 911
      end if
c     skip solvent block      
      do i=1,nva
        read(fn,1000,err=911,end=911) card
      end do
      do i=1,nu
        read(fn,1000,err=911,end=911) card
        i0 = 1
        call get_next_token(i0,card,sep,token,ostatus)
        read(token,*) tur(i)
        call get_next_token(i0,card,sep,token,ostatus)
        read(token,*) tua(i)
        call get_next_token(i0,card,sep,token,ostatus)
        read(token,*) iur(i)
      end do
c      do i=1,nu
c        write(*,*) tur(i),tua(i),iur(i)
c      end do
      return
911    continue       
c      if you reach this you are in trouble
       write(*,*) "found error"
       write(*,*) message
       stop

      end subroutine

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
