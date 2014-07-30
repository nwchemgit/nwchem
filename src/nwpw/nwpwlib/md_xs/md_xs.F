c $Id: md_xs.F 25353 2014-03-17 14:12:25Z jhammond $


c     **********************************************
c     *                                            *
c     *                md_xs_analysis              *
c     *                                            *
c     **********************************************
      subroutine md_xs_analysis(rtdb,filename)
      implicit none
      integer rtdb
      character*(*) filename

#include "mafdecls.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "stdio.fh"

*     **** local variables ****
      integer taskid_i,taskid_j,npi,npj
      integer taskid,np,MASTER
      parameter (MASTER=0)

      logical oprint,value
      character*20 codeversion
      character*2 tag2
      character*4 tag4
      integer nframes,natoms,stride,natoms1,center,nkf
      integer iii,i,k,icount,icpu,pto,pfrom
      real*8  rmax,x,y,z,vx,vy,vz,time1,volume,unita(3,3)

      integer nkf0
      parameter (nkf0 = 601)
      integer kf(2),chif(2),rion(2),tag(2)

      character*30 spectroscopy,absorber,edge
      character*80 chi_filename
      character*255 header,full_chi_filename

*     **** external functions ****
      integer  Parallel2d_convert_taskid_i
      external Parallel2d_convert_taskid_i
       

      call Parallel_np(np)
      call Parallel_taskid(taskid)
      oprint = (taskid.eq.MASTER)

      call Parallel2d_np_i(npi)
      call Parallel2d_np_j(npj)
      call Parallel2d_taskid_i(taskid_i)
      call Parallel2d_taskid_j(taskid_j)

#ifdef FEFF
      call feff_codeversion(codeversion)

*     **** get feff options from rtdb ****
      if (.not.rtdb_cget(rtdb,'md_xs:spectroscopy',1,spectroscopy)) 
     >   spectroscopy = 'exafs'
      if (.not.rtdb_cget(rtdb,'md_xs:absorber',1,absorber))
     >   absorber = ''
      if (.not.rtdb_cget(rtdb,'md_xs:edge',1,edge))
     >   edge = 'k'
      if (.not.rtdb_get(rtdb,'md_xs:center',mt_int,1,center))
     >   center = 1
      if (.not.rtdb_get(rtdb,'md_xs:rmax',mt_dbl,1,rmax))
     >   rmax = 10.0d0
      call util_file_name_noprefix('md_xs',.true.,.true.,header)

      if (.not.rtdb_get(rtdb,'md_xs:stride',mt_int,1,stride))
     >   stride = 1


*     **** define the chi data filename ****
      if (.not.rtdb_cget(rtdb,'md_xs:chi_filename',1,chi_filename))
     >    call util_file_prefix('chi',chi_filename)
      call util_file_name_noprefix(chi_filename,.false.,
     >                             .false.,
     >                             full_chi_filename)


      if (oprint) then
         write(luout,*) 
     >   "==== ","MD XS analysis: ",codeversion," ===="
         write(luout,*)
         write(luout,*) 'MOTION filename:  ',trim(filename)
      end if

      if (taskid.eq.MASTER) then
         call MOTION_nframes(filename,nframes,natoms,volume,unita)
      end if
      call Parallel_Brdcst_ivalue(0,nframes)
      call Parallel_Brdcst_ivalue(0,natoms)
      call Parallel_Brdcst_value(0,volume)
      call Parallel_Brdcst_values(0,9,unita)

      if (oprint) then
         write(luout,*) 'number of frames: ',nframes
         write(luout,*) 'stride:           ',stride
         write(luout,*) 'natoms:           ',natoms
         write(luout,*) 'volume:           ',volume
         write(luout,*)
         write(luout,*) 'spectroscopy:     ',spectroscopy
         write(luout,*) 'absorber:         ',absorber
         write(luout,*) 'edge:             ',edge
         write(luout,*) 'center:           ',center
         write(luout,*) 'rmax:             ',rmax
         write(luout,*)
      end if

*     **** allocate memory from stack ****
      value =  MA_push_get(mt_dbl,nkf0,'kf',kf(2),kf(1))
      value = value.and.
     >        MA_push_get(mt_dbl,nkf0,'chif',chif(2),chif(1))
      value = value.and.
     >        MA_push_get(mt_dbl,3*natoms,'rion',rion(2),rion(1)) 
      value = value.and.
     >        MA_push_get(mt_byte,2*natoms,'tag',tag(2),tag(1)) 
      if (.not.value) 
     >   call errquit("md_xs_analysis: out of stack",1,MA_ERR)


      if (taskid_i.eq.MASTER)
     >   open(unit=19,file=filename,form='formatted')

      call dcopy(nkf0,0.0d0,0,dbl_mb(kf(1)),1)
      call dcopy(nkf0,0.0d0,0,dbl_mb(chif(1)),1)
      icount = 0
      icpu   = 0
      do k=1,nframes
         if (taskid_i.eq.MASTER) then
            read(19,*) time1,natoms1,volume,unita
            do i=1,natoms1
               tag2 = '  '
               read(19,*) iii,tag2,tag4,x,y,z,vx,vy,vz
               byte_mb(tag(1)+2*(i-1))   = tag2(1:1)
               byte_mb(tag(1)+2*(i-1)+1) = tag2(2:2)
               dbl_mb(rion(1)+3*(i-1))   = x
               dbl_mb(rion(1)+3*(i-1)+1) = y
               dbl_mb(rion(1)+3*(i-1)+2) = z
            end do
         end if

         if (mod(k,stride).eq.0) then
            if (mod(icount,npj).eq.taskid_j) then
                 if (taskid_i.eq.MASTER) then
                    if (icpu.eq.MASTER) then
                       !*** compute md_xs with natoms1,tag,rion +--> kf,chif ***
                       call md_xs(header,spectroscopy,absorber,edge,
     >                            center,rmax,
     >                         natoms1,byte_mb(tag(1)),dbl_mb(rion(1)),
     >                         nkf,nkf0,dbl_mb(kf(1)),dbl_mb(chif(1)))
                    else
                       !*** send natoms1,tag,rion to icpu ***
                       pto = Parallel2d_convert_taskid_i(icpu)
                       call Parallel_send_ivalues(pto,3*k,1,natoms1)
                       call Parallel_send_characters(pto,3*k+1,
     >                                        2*natoms1,byte_mb(tag(1)))
                       call Parallel_send_values(pto,3*k+2,
     >                                        3*natoms1,dbl_mb(rion(1)))
                    end if
                 else
                    if (icpu.eq.taskid_i) then
                       !*** receive natoms1,tag,rion from MASTER ***
                       pfrom = Parallel2d_convert_taskid_i(MASTER)
                       call Parallel_recv_ivalues(pfrom,3*k,1,natoms1)
                       call Parallel_recv_characters(pfrom,3*k+1,
     >                                        2*natoms1,byte_mb(tag(1)))
                       call Parallel_recv_values(pfrom,3*k+2,
     >                                        3*natoms1,dbl_mb(rion(1)))

                       !*** compute md_xs with natoms1,tag,rion +--> kf,chif ***
                       call md_xs(header,spectroscopy,absorber,edge,
     >                            center,rmax,
     >                          natoms1,byte_mb(tag(1)),dbl_mb(rion(1)),
     >                          nkf,nkf0,dbl_mb(kf(1)),dbl_mb(chif(1)))
                    end if
                 end if
                 icpu = mod(icpu+1,npi)
            end if
            icount = icount + 1
         end if
      end do
      if (taskid_i.eq.MASTER) close(19)

      x = 1.0d0/dble(icount)
      call Parallel_Vector_SumAll(nkf0,dbl_mb(chif(1)))
      call dscal(nkf0,x,dbl_mb(chif(1)),1)

      
*     **** write out chi data file ****
      if (oprint) write(luout,*)
     >            'Chi data filename:',trim(full_chi_filename)

      if (taskid.eq.MASTER) then
         open(unit=19,file=full_chi_filename,form='formatted')
         write(19,'(A,A)') '# Chi data filename:',
     >                     trim(full_chi_filename)
         write(19,'(A,I4)') "# nkf =",nkf
         write(19,'(A,A)') "#spectroscopy:     ",trim(spectroscopy)
         write(19,'(A,A)') "#absorber:         ",trim(absorber)
         write(19,'(A,A)') "#edge:             ",trim(edge)
         write(19,'(A,I6)') "#center=",center
         write(19,'(A,F16.8)') "#rmax=",rmax
         write(19,'(A,A)') '# MOTION filename:  ',trim(filename)
         write(19,'(A,I8)') '# number of frames=',nframes
         write(19,'(A,I8)') '# number of frames used, icount=',icount
         write(19,'(A,I8)') '# stride=',stride
         write(19,'(A,I8)') '# natoms=',natoms
         write(19,'(A,F16.8)') '#volume=',volume
         do i=1,nkf
            write(19,'(4F16.8)') dbl_mb(kf(1)+i-1),dbl_mb(chif(1)+i-1)
         end do
         close(19)
      end if

c     **** deallocate stack ****
      value =           MA_pop_stack(tag(2))
      value = value.and.MA_pop_stack(rion(2))
      value = value.and.MA_pop_stack(chif(2))
      value = value.and.MA_pop_stack(kf(2))
      if (.not.value)
     >   call errquit("md_xs_analysis: popping stack",3,MA_ERR)

#else
      if (oprint) then
         write(luout,*) "==== ",
     >     "FEFF library not inlucded in compile: ",
     >     "MD XS analysis not performed"," ===="
      end if
#endif
      return
      end



#ifdef FEFF
c     **********************************************
c     *                                            *
c     *                md_xs                       *
c     *                                            *
c     **********************************************
      subroutine md_xs(header,spectroscopy,absorber,edge,
     >                  center,rmax,nion,symb,rion,nkf,nkf0,kf,chi)
      implicit none
      character*(*) header,spectroscopy,absorber,edge
      integer center
      real*8  rmax
      integer nion
      character*2 symb(*)
      real*8      rion(3,*)
      integer nkf,nkf0
      real*8 kf(*)
      real*8 chi(*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "stdio.fh"

*     **** local variables ****
      integer taskid,MASTER,np
      parameter (MASTER=0)

      logical value,found
      integer nkf4,j,ii,ia,zi
      integer zion(2),katm(2),zkatm(2),nkatm
      integer zion4(2),katm4(2),zkatm4(2),tchi(2)
      character*2 item

      character*2 symbols(112)
      data symbols/
     $     'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne',
     $     'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', 'K ', 'Ca',
     $     'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
     $     'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y ', 'Zr',
     $     'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
     $     'Sb', 'Te', 'I ', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
     $     'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
     $     'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
     $     'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
     $     'Pa', 'U ', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',
     $     'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds',
     $     'Rg', 'Cn'/

      call Parallel_taskid(taskid)
      call Parallel_np(np)

*     ***** allocate stack ****
      value = MA_push_get(mt_int,nion,'md_xs:zion',zion(2),zion(1))
      value = value.and.
     >        MA_push_get(mt_int,nion,'md_xs:zkatm',zkatm(2),zkatm(1))
      value = value.and.
     >        MA_push_get(mt_int,nion,'md_xs:katm',katm(2),katm(1))
      value = value.and.
     >        MA_push_get(mt_int,nion,'md_xs:zion4',zion4(2),zion4(1))
      value = value.and.
     >       MA_push_get(mt_int,nion,'md_xs:zkatm4',zkatm4(2),zkatm4(1))
      value = value.and.
     >        MA_push_get(mt_int,nion,'md_xs:katm4',katm4(2),katm4(1))
      value = value.and.
     >        MA_push_get(mt_dbl,nkf0,'md_xs:tchi',tchi(2),tchi(1))
      if (.not.value)
     >   call errquit("md_xs: pushing stack",0,MA_ERR)
      call dcopy(nkf0,0.0d0,0,dbl_mb(tchi(1)),1)

      do ii=1,nion
         item = '  '
         item = symb(ii)
         zi = -1
         do j=1,112
            if (item.eq.symbols(j)) zi = j
         end do
         int_mb(zion(1)+ii-1) = zi
      end do

      nkatm = 0
      do ii=1,nion
         found = .false.
         do j=1,nkatm
            if (int_mb(zion(1)+ii-1).eq.int_mb(zkatm(1)+j-1)) then
               found = .true.
               ia = j
            end if
         end do
         if (found) then
            int_mb(katm(1)+ii-1) = ia
         else
            nkatm = nkatm + 1
            int_mb(zkatm(1)+nkatm-1) = int_mb(zion(1)+ii-1)
            int_mb(katm(1)+ii-1) = nkatm
         end if
      end do



#ifdef FEFF8
      call icopy(nion,int_mb(zion(1)),1,int_mb(zion4(1)),1)
      call icopy(nion,int_mb(katm(1)),1,int_mb(katm4(1)),1)
      call icopy(nkatm,int_mb(zkatm(1)),1,int_mb(zkatm4(1)),1)
#else
      call ftupid_icopy8to4(nion,int_mb(zion(1)),int_mb(zion4(1)))
      call ftupid_icopy8to4(nion,int_mb(katm(1)),int_mb(katm4(1)))
      call ftupid_icopy8to4(nkatm,int_mb(zkatm(1)),int_mb(zkatm4(1)))
#endif
      call feff_fortran(header,spectroscopy,absorber,edge,
     >                  center,rmax,
     >                  nkatm,int_mb(katm4(1)),int_mb(zkatm4(1)),
     >                  nion,int_mb(zion4(1)),rion,
     >                  nkf4,kf,dbl_mb(tchi(1)))
#ifdef FEFF8
      nkf = nkf4
#else
      call ftupid_icopy4to8(1,nkf4,nkf)
#endif

*     **** chi += tchi  ****
      call daxpy(nkf0,1.0d0,dbl_mb(tchi(1)),1,chi,1)


*     **** popping stack ****
      value =           MA_pop_stack(tchi(2))
      value = value.and.MA_pop_stack(katm4(2))
      value = value.and.MA_pop_stack(zkatm4(2))
      value = value.and.MA_pop_stack(zion4(2))
      value = value.and.MA_pop_stack(katm(2))
      value = value.and.MA_pop_stack(zkatm(2))
      value = value.and.MA_pop_stack(zion(2))
      if (.not.value)
     >   call errquit("md_xs: popping stack",1,MA_ERR)
     
      return
      end


cc     **********************************************
cc     *                                            *
cc     *                md_xs_create_json           *
cc     *                                            *
cc     **********************************************
c      subroutine md_xs_create_json(rtdb,header,nion,symbols,rion,json)
c      implicit none
c      integer rtdb
c      character*(*) header
c      integer nion
c      character*2 symbols(*)
c      real*8 rion(3,*)
c      character*(*) json
c#include "mafdecls.fh"
c#include "rtdb.fh"
c#include "errquit.fh"
c#include "stdio.fh"
c
c      integer ii,center
c      real*8 rmax
c      character*30 tstr1,tstr2,tstr3
c
c*     **** open json ****
c      json = "{"
c
c*     **** "spectroscopy": ****
c      if (rtdb_cget(rtdb,'md_xs:spectroscopy',1,tstr1)) then
c         json = trim(json)//"""spectroscopy"":"
c         json = trim(json)//" "//""""//trim(tstr1)//""","
c      end if
c
c*     **** "absorber": ****
c      if (rtdb_cget(rtdb,'md_xs:absorber',1,tstr1)) then
c         json = trim(json)//" ""absorber"": ["
c         json = trim(json)//""""//trim(tstr1)//"""],"
c      end if
c
c*     **** "edge": ****
c      if (rtdb_cget(rtdb,'md_xs:edge',1,tstr1)) then
c         json = trim(json)//" ""edge"":"
c         json = trim(json)//" "//""""//trim(tstr1)//""","
c      end if
c
c*     **** "center": ****
c      if (rtdb_get(rtdb,'md_xs:center',mt_int,1,center)) then
c         write(tstr1,'(I6)') center
c         json = trim(json)//" ""center"": ["
c         json = trim(json)//trim(adjustl(tstr1))//"],"
c      end if
c
c*     **** "rmax": ****
c      if (rtdb_get(rtdb,'md_xs:rmax',mt_dbl,1,rmax)) then
c         write(tstr1,'(F20.6)') rmax
c         json = trim(json)//" ""rmax"":"
c         json = trim(json)//" "//trim(adjustl(tstr1))//","
c      end if
c
c*     **** "geometry": ****
c      json = trim(json)//" ""geometry"": ["
c      do ii=1,nion
c         write(tstr1,'(E22.15)') rion(1,ii)
c         write(tstr2,'(E22.15)') rion(2,ii)
c         write(tstr3,'(E22.15)') rion(3,ii)
c         if (ii.eq.1) then
c            json = trim(json)//"["""//trim(symbols(ii))//""","
c         else
c            json = trim(json)//" ["""//trim(symbols(ii))//""","
c         end if
c         json = trim(json)//" "//trim(adjustl(tstr1))//","
c         json = trim(json)//" "//trim(adjustl(tstr2))//","
c         json = trim(json)//" "//trim(adjustl(tstr3))//"]"
c         if (ii.lt.nion) json = trim(json)//","
c      end do
c      json = trim(json)//"], "
c
c*     **** "scratch_dir": ****
c      json = trim(json)//" ""scratch_dir"":"
c      json = trim(json)//" "//""""//trim(header)//""""
c
c*     **** close json ****
c      json = trim(json)//"}"
c
c      return
c      end


      subroutine ftupid_icopy8to4(n,array8,array4)
      implicit none
      integer n
      integer*8 array8(*)
      integer*4 array4(*)
      integer i
      do i=1,n
          array4(i) = array8(i)
      end do
      return
      end

      subroutine ftupid_icopy4to8(n,array4,array8)
      implicit none
      integer n
      integer*4 array4(*)
      integer*8 array8(*)
      integer i
      do i=1,n
          array8(i) = array4(i)
      end do
      return
      end

#endif

