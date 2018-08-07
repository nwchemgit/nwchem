C$Id$

cccccceccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    Tee matrix SYMOPS contains the matrix reps. of all group operators
c   except the identity. The variable NOPS holds the number of operators
c   in SYMOPS.
c
c***********************************************************************
      subroutine gensym_extra(itype,numgrp,numset,symops,nops,oprint,
     >                        group_name)
      implicit none
      integer maxops
      parameter(maxops=192)

      integer itype,numgrp,numset
      real*8 symops(maxops*3,4)
      integer nops
      character*(*) group_name
      logical oprint

#include "inp.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "util.fh"
#include "mafdecls.fh"

*     **** local varialbles ****
      logical mprint,hprint,debug,does_it_exist
      logical from_environment,from_compile,from_nwchemrc
      integer iop,lgth,unitf,print_level,i,j
      character*255 extra_group_libname,basis_library,test

*     **** external functions ****
      logical  util_find_dir,util_io_unit
      external util_find_dir,util_io_unit

*
* order of precedence for choosing name
* 1) value of NWCHEM_EXTRA_GROUPS environment variable
* 2) value of NWCHEM_EXTRA_GROUPS set in $HOME/.nwchemrc file
* 3) value of the compiled in library name
*
      call util_print_get_level(print_level)
      mprint = print_medium.le.print_level
      hprint = print_high  .le.print_level
      debug  = print_debug .le.print_level
      from_environment = .false.
      from_nwchemrc    = .false.
      from_compile     = .false.


      nops = 0


*     **** Try to get from NWCHEM_EXTRA_GROUPS environment variable ****
      call util_getenv('NWCHEM_EXTRA_GROUPS',extra_group_libname)
      if (lgth.gt.0) then
         if (util_find_dir(extra_group_libname)) then
            from_environment = .true.
            goto 99
         else
            write(luout,*)' warning:::::::::::::: from_environment'
            write(luout,*)' NWCHEM_EXTRA_GROUPS set to: <',
     &       extra_group_libname(1:inp_strlen(extra_group_libname)),'>'
            write(luout,*)' but file does not exist !'
            write(luout,*)' using compiled library'
         end if
      end if


*     **** Try to get from NWCHEM_EXTRA_GROUPS nwchemrc ****
*2: check for NWCHEM_EXTRA_GROUPS defined in users .nwchemrc file
*   assumed structure in .nwchemrc file is variable [whitespace] value
*   one setting per line
*
      basis_library='nwchem_extra_group'
      call inp_save_state() ! save state of any inp unit
      if(.not.util_nwchemrc_get(basis_library,extra_group_libname)) then
        if (debug) then
          write(luout,*)'util_nwchemrc_get failed'
        endif
      else
        does_it_exist=util_find_dir(extra_group_libname)
        if (does_it_exist)then
          from_nwchemrc = .true.
          call inp_restore_state() ! restore state of any inp unit
          goto 99
        else
          write(luout,*)' warning:::::::::::::: from_nwchemrc'
          write(luout,*)' NWCHEM_EXTRA_GROUPS set to: <',
     &     extra_group_libname(1:inp_strlen(extra_group_libname)),'>'
          write(luout,*)' but file does not exist !'
          write(luout,*)' using compiled in library'
        endif
      endif
      call inp_restore_state() ! restore state of any inp unit



*     **** Try to get from compile ****
      from_compile = .true.
      call util_nwchem_srcdir(extra_group_libname)
      extra_group_libname
     >     =extra_group_libname(1:inp_strlen(extra_group_libname))
     >     //"/symmetry/extragroups.dat"
      if (util_find_dir(extra_group_libname)) then
         goto 99
      else
         write(luout,*)' warning:::::::::::::: from_compile'
         write(luout,*)' NWCHEM_EXTRA_GROPUS is: <',
     &    extra_group_libname(1:inp_strlen(extra_group_libname)),'>'
         write(luout,*)' but file does ','not exist or you ', 
     &        'do not have ','access to it !'
      endif
      call errquit('gensym_extra: extragroups.dat not found',0,0)

 99   continue

      if (from_environment) then
          if (mprint)
     >     write(luout,*)
     >     ' extra_group_library name resolved from: environment'
      else if (from_nwchemrc) then
          if (mprint)
     >     write(luout,*)
     >     ' extra_group_library name resolved from: .nwchemrc'
      else
          if (mprint)
     >     write(luout,*)
     >     ' extra_group_library name resolved from: compiled reference'
      endif
      if (mprint) then
         write(luout,*) ' NWCHEM_EXTRA_GROUPS set to: <',
     >    extra_group_libname(1:inp_strlen(extra_group_libname)),'>'
      end if

      if(.not.util_io_unit(80,90,unitf)) 
     >  call errquit("gensym_extra cannot get io unit",0,DISK_ERR)

      open(unit=unitf,file=extra_group_libname,status='old',
     >     form='formatted',ERR=999)
      call inp_init(unitf,luout)


*     **** look for symmetry group ****
 100  if (.not. inp_read()) goto 200
      if (.not. inp_a(test))
     > call errquit('gensym_extra: failed string', 0,INPUT_ERR)
      if (.not.inp_compare(.false.,test,group_name)) goto 100

*     **** found the symmetry group - define the symmetry ops ****
      if (inp_read()) then
         if (inp_i(numgrp).and.(inp_read())) then
            if (inp_i(nops).and.(inp_read())) then
               do iop=0,nops-1
                  do i=1,3
                     do j=1,4
                        if (.not.inp_f(symops(3*iop+i,j)))
     >                  call errquit('gensym_extra:reading symops', 
     >                            0,INPUT_ERR)
                     end do
                     if (.not.inp_read())
     >               call errquit('gensym_extra:failed', 
     >                            1,INPUT_ERR)
                  end do
               end do
            end if
         end if
      end if

 200  continue
      call inp_init(luin,luout)
      close(unitf)



******* Extra Triclinic groups ******
c 
c    A1 
c 
      if (numgrp.eq.231) then
      group_name = 'A1'
      nops = 1
      iop = 0
      symops(3*iop+1,1) = 1.0d0
      symops(3*iop+2,1) = 0.0d0
      symops(3*iop+3,1) = 0.0d0
      symops(3*iop+1,2) = 0.0d0
      symops(3*iop+2,2) = 1.0d0
      symops(3*iop+3,2) = 0.0d0
      symops(3*iop+1,3) = 0.0d0
      symops(3*iop+2,3) = 0.0d0
      symops(3*iop+3,3) = 1.0d0
      symops(3*iop+1,4) = 0.0d0
      symops(3*iop+2,4) = 0.5d0
      symops(3*iop+3,4) = 0.5d0
      
      end if
c 
c    B1 
c 
      if (numgrp.eq.232) then
      group_name = 'B1'
      nops = 1
      iop = 0
      symops(3*iop+1,1) = 1.0d0
      symops(3*iop+2,1) = 0.0d0
      symops(3*iop+3,1) = 0.0d0
      symops(3*iop+1,2) = 0.0d0
      symops(3*iop+2,2) = 1.0d0
      symops(3*iop+3,2) = 0.0d0
      symops(3*iop+1,3) = 0.0d0
      symops(3*iop+2,3) = 0.0d0
      symops(3*iop+3,3) = 1.0d0
      symops(3*iop+1,4) = 0.5d0
      symops(3*iop+2,4) = 0.0d0
      symops(3*iop+3,4) = 0.5d0
      end if
c 
c    C1 
c 
      if (numgrp.eq.233) then
      group_name = 'C1'
      nops = 1
      iop = 0
      symops(3*iop+1,1) = 1.0d0
      symops(3*iop+2,1) = 0.0d0
      symops(3*iop+3,1) = 0.0d0
      symops(3*iop+1,2) = 0.0d0
      symops(3*iop+2,2) = 1.0d0
      symops(3*iop+3,2) = 0.0d0
      symops(3*iop+1,3) = 0.0d0
      symops(3*iop+2,3) = 0.0d0
      symops(3*iop+3,3) = 1.0d0
      symops(3*iop+1,4) = 0.5d0
      symops(3*iop+2,4) = 0.5d0
      symops(3*iop+3,4) = 0.0d0
      end if
c 
c    F1 
c 
      if (numgrp.eq.234) then
      group_name = 'F1'
      nops = 3
      iop = 0
      symops(3*iop+1,1) = 1.0d0
      symops(3*iop+2,1) = 0.0d0
      symops(3*iop+3,1) = 0.0d0
      symops(3*iop+1,2) = 0.0d0
      symops(3*iop+2,2) = 1.0d0
      symops(3*iop+3,2) = 0.0d0
      symops(3*iop+1,3) = 0.0d0
      symops(3*iop+2,3) = 0.0d0
      symops(3*iop+3,3) = 1.0d0
      symops(3*iop+1,4) = 0.0d0
      symops(3*iop+2,4) = 0.5d0
      symops(3*iop+3,4) = 0.5d0
      iop = 1
      symops(3*iop+1,1) = 1.0d0
      symops(3*iop+2,1) = 0.0d0
      symops(3*iop+3,1) = 0.0d0
      symops(3*iop+1,2) = 0.0d0
      symops(3*iop+2,2) = 1.0d0
      symops(3*iop+3,2) = 0.0d0
      symops(3*iop+1,3) = 0.0d0
      symops(3*iop+2,3) = 0.0d0
      symops(3*iop+3,3) = 1.0d0
      symops(3*iop+1,4) = 0.5d0
      symops(3*iop+2,4) = 0.0d0
      symops(3*iop+3,4) = 0.5d0
      iop = 2
      symops(3*iop+1,1) = 1.0d0
      symops(3*iop+2,1) = 0.0d0
      symops(3*iop+3,1) = 0.0d0
      symops(3*iop+1,2) = 0.0d0
      symops(3*iop+2,2) = 1.0d0
      symops(3*iop+3,2) = 0.0d0
      symops(3*iop+1,3) = 0.0d0
      symops(3*iop+2,3) = 0.0d0
      symops(3*iop+3,3) = 1.0d0
      symops(3*iop+1,4) = 0.5d0
      symops(3*iop+2,4) = 0.5d0
      symops(3*iop+3,4) = 0.0d0
      end if
c 
c    I1 
c 
      if (numgrp.eq.235) then
      group_name = 'I1'
      nops = 1
      iop = 0
      symops(3*iop+1,1) = 1.0d0
      symops(3*iop+2,1) = 0.0d0
      symops(3*iop+3,1) = 0.0d0
      symops(3*iop+1,2) = 0.0d0
      symops(3*iop+2,2) = 1.0d0
      symops(3*iop+3,2) = 0.0d0
      symops(3*iop+1,3) = 0.0d0
      symops(3*iop+2,3) = 0.0d0
      symops(3*iop+3,3) = 1.0d0
      symops(3*iop+1,4) = 0.5d0
      symops(3*iop+2,4) = 0.5d0
      symops(3*iop+3,4) = 0.5d0
      end if

c 
c    A-1
c 
      if (numgrp.eq.236) then
      group_name = 'A-1'
      nops = 3
      iop = 0
      symops(3*iop+1,1) = -1.0d0
      symops(3*iop+2,1) = 0.0d0
      symops(3*iop+3,1) = 0.0d0
      symops(3*iop+1,2) = 0.0d0
      symops(3*iop+2,2) = -1.0d0
      symops(3*iop+3,2) = 0.0d0
      symops(3*iop+1,3) = 0.0d0
      symops(3*iop+2,3) = 0.0d0
      symops(3*iop+3,3) = -1.0d0
      symops(3*iop+1,4) = 0.0d0
      symops(3*iop+2,4) = 0.0d0
      symops(3*iop+3,4) = 0.0d0
      iop = 1
      symops(3*iop+1,1) = 1.0d0
      symops(3*iop+2,1) = 0.0d0
      symops(3*iop+3,1) = 0.0d0
      symops(3*iop+1,2) = 0.0d0
      symops(3*iop+2,2) = 1.0d0
      symops(3*iop+3,2) = 0.0d0
      symops(3*iop+1,3) = 0.0d0
      symops(3*iop+2,3) = 0.0d0
      symops(3*iop+3,3) = 1.0d0
      symops(3*iop+1,4) = 0.0d0
      symops(3*iop+2,4) = 0.5d0
      symops(3*iop+3,4) = 0.5d0
      iop = 2
      symops(3*iop+1,1) = -1.0d0
      symops(3*iop+2,1) = 0.0d0
      symops(3*iop+3,1) = 0.0d0
      symops(3*iop+1,2) = 0.0d0
      symops(3*iop+2,2) = -1.0d0
      symops(3*iop+3,2) = 0.0d0
      symops(3*iop+1,3) = 0.0d0
      symops(3*iop+2,3) = 0.0d0
      symops(3*iop+3,3) = -1.0d0
      symops(3*iop+1,4) = 0.0d0
      symops(3*iop+2,4) = 0.5d0
      symops(3*iop+3,4) = 0.5d0
      end if
c 
c    B-1
c 
      if (numgrp.eq.237) then
      group_name = 'B-1'
      nops = 3
      iop = 0
      symops(3*iop+1,1) = -1.0d0
      symops(3*iop+2,1) = 0.0d0
      symops(3*iop+3,1) = 0.0d0
      symops(3*iop+1,2) = 0.0d0
      symops(3*iop+2,2) = -1.0d0
      symops(3*iop+3,2) = 0.0d0
      symops(3*iop+1,3) = 0.0d0
      symops(3*iop+2,3) = 0.0d0
      symops(3*iop+3,3) = -1.0d0
      symops(3*iop+1,4) = 0.0d0
      symops(3*iop+2,4) = 0.0d0
      symops(3*iop+3,4) = 0.0d0
      iop = 1
      symops(3*iop+1,1) = 1.0d0
      symops(3*iop+2,1) = 0.0d0
      symops(3*iop+3,1) = 0.0d0
      symops(3*iop+1,2) = 0.0d0
      symops(3*iop+2,2) = 1.0d0
      symops(3*iop+3,2) = 0.0d0
      symops(3*iop+1,3) = 0.0d0
      symops(3*iop+2,3) = 0.0d0
      symops(3*iop+3,3) = 1.0d0
      symops(3*iop+1,4) = 0.5d0
      symops(3*iop+2,4) = 0.0d0
      symops(3*iop+3,4) = 0.5d0
      iop = 2
      symops(3*iop+1,1) = -1.0d0
      symops(3*iop+2,1) = 0.0d0
      symops(3*iop+3,1) = 0.0d0
      symops(3*iop+1,2) = 0.0d0
      symops(3*iop+2,2) = -1.0d0
      symops(3*iop+3,2) = 0.0d0
      symops(3*iop+1,3) = 0.0d0
      symops(3*iop+2,3) = 0.0d0
      symops(3*iop+3,3) = -1.0d0
      symops(3*iop+1,4) = 0.5d0
      symops(3*iop+2,4) = 0.0d0
      symops(3*iop+3,4) = 0.5d0
      end if
c 
c    C-1 
c 
      if (numgrp.eq.238) then
      group_name = 'C-1'
      nops = 3
      iop = 0
      symops(3*iop+1,1) = -1.0d0
      symops(3*iop+2,1) = 0.0d0
      symops(3*iop+3,1) = 0.0d0
      symops(3*iop+1,2) = 0.0d0
      symops(3*iop+2,2) = -1.0d0
      symops(3*iop+3,2) = 0.0d0
      symops(3*iop+1,3) = 0.0d0
      symops(3*iop+2,3) = 0.0d0
      symops(3*iop+3,3) = -1.0d0
      symops(3*iop+1,4) = 0.0d0
      symops(3*iop+2,4) = 0.0d0
      symops(3*iop+3,4) = 0.0d0
      iop = 1
      symops(3*iop+1,1) = 1.0d0
      symops(3*iop+2,1) = 0.0d0
      symops(3*iop+3,1) = 0.0d0
      symops(3*iop+1,2) = 0.0d0
      symops(3*iop+2,2) = 1.0d0
      symops(3*iop+3,2) = 0.0d0
      symops(3*iop+1,3) = 0.0d0
      symops(3*iop+2,3) = 0.0d0
      symops(3*iop+3,3) = 1.0d0
      symops(3*iop+1,4) = 0.5d0
      symops(3*iop+2,4) = 0.5d0
      symops(3*iop+3,4) = 0.0d0
      iop = 2
      symops(3*iop+1,1) = -1.0d0
      symops(3*iop+2,1) = 0.0d0
      symops(3*iop+3,1) = 0.0d0
      symops(3*iop+1,2) = 0.0d0
      symops(3*iop+2,2) = -1.0d0
      symops(3*iop+3,2) = 0.0d0
      symops(3*iop+1,3) = 0.0d0
      symops(3*iop+2,3) = 0.0d0
      symops(3*iop+3,3) = -1.0d0
      symops(3*iop+1,4) = 0.5d0
      symops(3*iop+2,4) = 0.5d0
      symops(3*iop+3,4) = 0.0d0
      end if
c 
c    F-1 
c 
      if (numgrp.eq.239) then
      group_name = 'F-1'
      nops = 7
      iop = 0
      symops(3*iop+1,1) = -1.0d0
      symops(3*iop+2,1) = 0.0d0
      symops(3*iop+3,1) = 0.0d0
      symops(3*iop+1,2) = 0.0d0
      symops(3*iop+2,2) = -1.0d0
      symops(3*iop+3,2) = 0.0d0
      symops(3*iop+1,3) = 0.0d0
      symops(3*iop+2,3) = 0.0d0
      symops(3*iop+3,3) = -1.0d0
      symops(3*iop+1,4) = 0.0d0
      symops(3*iop+2,4) = 0.0d0
      symops(3*iop+3,4) = 0.0d0
      iop = 1
      symops(3*iop+1,1) = 1.0d0
      symops(3*iop+2,1) = 0.0d0
      symops(3*iop+3,1) = 0.0d0
      symops(3*iop+1,2) = 0.0d0
      symops(3*iop+2,2) = 1.0d0
      symops(3*iop+3,2) = 0.0d0
      symops(3*iop+1,3) = 0.0d0
      symops(3*iop+2,3) = 0.0d0
      symops(3*iop+3,3) = 1.0d0
      symops(3*iop+1,4) = 0.0d0
      symops(3*iop+2,4) = 0.5d0
      symops(3*iop+3,4) = 0.5d0
      iop = 2
      symops(3*iop+1,1) = -1.0d0
      symops(3*iop+2,1) = 0.0d0
      symops(3*iop+3,1) = 0.0d0
      symops(3*iop+1,2) = 0.0d0
      symops(3*iop+2,2) = -1.0d0
      symops(3*iop+3,2) = 0.0d0
      symops(3*iop+1,3) = 0.0d0
      symops(3*iop+2,3) = 0.0d0
      symops(3*iop+3,3) = -1.0d0
      symops(3*iop+1,4) = 0.0d0
      symops(3*iop+2,4) = 0.5d0
      symops(3*iop+3,4) = 0.5d0
      iop = 3
      symops(3*iop+1,1) = 1.0d0
      symops(3*iop+2,1) = 0.0d0
      symops(3*iop+3,1) = 0.0d0
      symops(3*iop+1,2) = 0.0d0
      symops(3*iop+2,2) = 1.0d0
      symops(3*iop+3,2) = 0.0d0
      symops(3*iop+1,3) = 0.0d0
      symops(3*iop+2,3) = 0.0d0
      symops(3*iop+3,3) = 1.0d0
      symops(3*iop+1,4) = 0.5d0
      symops(3*iop+2,4) = 0.0d0
      symops(3*iop+3,4) = 0.5d0
      iop = 4
      symops(3*iop+1,1) = -1.0d0
      symops(3*iop+2,1) = 0.0d0
      symops(3*iop+3,1) = 0.0d0
      symops(3*iop+1,2) = 0.0d0
      symops(3*iop+2,2) = -1.0d0
      symops(3*iop+3,2) = 0.0d0
      symops(3*iop+1,3) = 0.0d0
      symops(3*iop+2,3) = 0.0d0
      symops(3*iop+3,3) = -1.0d0
      symops(3*iop+1,4) = 0.5d0
      symops(3*iop+2,4) = 0.0d0
      symops(3*iop+3,4) = 0.5d0
      iop = 5
      symops(3*iop+1,1) = 1.0d0
      symops(3*iop+2,1) = 0.0d0
      symops(3*iop+3,1) = 0.0d0
      symops(3*iop+1,2) = 0.0d0
      symops(3*iop+2,2) = 1.0d0
      symops(3*iop+3,2) = 0.0d0
      symops(3*iop+1,3) = 0.0d0
      symops(3*iop+2,3) = 0.0d0
      symops(3*iop+3,3) = 1.0d0
      symops(3*iop+1,4) = 0.5d0
      symops(3*iop+2,4) = 0.5d0
      symops(3*iop+3,4) = 0.0d0
      iop = 6
      symops(3*iop+1,1) = -1.0d0
      symops(3*iop+2,1) = 0.0d0
      symops(3*iop+3,1) = 0.0d0
      symops(3*iop+1,2) = 0.0d0
      symops(3*iop+2,2) = -1.0d0
      symops(3*iop+3,2) = 0.0d0
      symops(3*iop+1,3) = 0.0d0
      symops(3*iop+2,3) = 0.0d0
      symops(3*iop+3,3) = -1.0d0
      symops(3*iop+1,4) = 0.5d0
      symops(3*iop+2,4) = 0.5d0
      symops(3*iop+3,4) = 0.0d0
      end if
c 
c    I-1 
c 
      if (numgrp.eq.240) then
      group_name = 'I-1'
      nops = 3
      iop = 0
      symops(3*iop+1,1) = -1.0d0
      symops(3*iop+2,1) = 0.0d0
      symops(3*iop+3,1) = 0.0d0
      symops(3*iop+1,2) = 0.0d0
      symops(3*iop+2,2) = -1.0d0
      symops(3*iop+3,2) = 0.0d0
      symops(3*iop+1,3) = 0.0d0
      symops(3*iop+2,3) = 0.0d0
      symops(3*iop+3,3) = -1.0d0
      symops(3*iop+1,4) = 0.0d0
      symops(3*iop+2,4) = 0.0d0
      symops(3*iop+3,4) = 0.0d0
      iop = 1
      symops(3*iop+1,1) = 1.0d0
      symops(3*iop+2,1) = 0.0d0
      symops(3*iop+3,1) = 0.0d0
      symops(3*iop+1,2) = 0.0d0
      symops(3*iop+2,2) = 1.0d0
      symops(3*iop+3,2) = 0.0d0
      symops(3*iop+1,3) = 0.0d0
      symops(3*iop+2,3) = 0.0d0
      symops(3*iop+3,3) = 1.0d0
      symops(3*iop+1,4) = 0.5d0
      symops(3*iop+2,4) = 0.5d0
      symops(3*iop+3,4) = 0.5d0
      iop = 2
      symops(3*iop+1,1) = -1.0d0
      symops(3*iop+2,1) = 0.0d0
      symops(3*iop+3,1) = 0.0d0
      symops(3*iop+1,2) = 0.0d0
      symops(3*iop+2,2) = -1.0d0
      symops(3*iop+3,2) = 0.0d0
      symops(3*iop+1,3) = 0.0d0
      symops(3*iop+2,3) = 0.0d0
      symops(3*iop+3,3) = -1.0d0
      symops(3*iop+1,4) = 0.5d0
      symops(3*iop+2,4) = 0.5d0
      symops(3*iop+3,4) = 0.5d0
      end if


      return

 999  continue

      return
      end

