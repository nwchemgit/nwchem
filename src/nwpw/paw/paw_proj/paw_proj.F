*
* $Id: paw_proj.F,v 1.4 2002-08-01 21:27:27 bylaska Exp $
*

!**************************************************
!
!	Name: paw_proj_init
!
!	Purpose: initializes the paw projectors
!
!	Created:	7/30/2002
!**************************************************
      subroutine paw_proj_init()
      implicit none

#include "mafdecls.fh"
#include "paw_proj.fh" 
  
*     **** local variables ****
      logical value,found
      integer ia,ntotal,npack1
      character*4 element
      character*20 fname

*     **** external functions ****
      logical     nwpw_filefind
      integer     ion_nkatm
      character*2 ion_atom
      external    nwpw_filefind
      external    ion_nkatm
      external    ion_atom

      nkatm = ion_nkatm()
      value = MA_alloc_get(mt_int,(nkatm),
     >                     'nbasis',nbasis(2),nbasis(1))
      if (.not.value) call errquit('paw_proj_init: alloc heap',0)


*     *** make sure prj names are formatted correctly ***
      call paw_proj_check_format()

*     *** get the size of prj common block and allocate it ***
      call paw_proj_get_nbasis(int_mb(nbasis(1)))
      ntotal = 0
      do ia=1,nkatm
        ntotal = ntotal + int_mb(nbasis(1)+ia-1)
      end do
      call Pack_npack(1,npack1)
      value = MA_alloc_get(mt_dcpl,(ntotal*npack1),
     >                     'prj',prj(2),prj(1))
      if (.not.value) call errquit('paw_proj_init: alloc heap',1)

*     *** read in formatted prj's to prj common block ***

      return
      end

!**************************************************
!
!	Name: paw_proj_end
!
!	Purpose: removes space used by the paw projectors
!
!	Created:	7/30/2002
!**************************************************
      subroutine paw_proj_end()
      implicit none
      
#include "mafdecls.fh"
#include "paw_proj.fh"

      logical value

      value =           MA_free_heap(nbasis(2))
      value = value.and.MA_free_heap(prj(2))
      if (.not.value) call errquit('paw_proj_end: dealloc heap',0)
      return
      end


!**************************************************
!
!	Name: paw_proj_nbasis
!
!	Purpose: returns the number of the paw projectors
!                for this kind of atom.
!
!	Created:	7/30/2002
!**************************************************
      integer function paw_proj_nbasis(ia)
      implicit none
      integer ia
      
#include "mafdecls.fh"
#include "paw_proj.fh"

      paw_proj_nbasis = int_mb(nbasis(1)+ia-1)
      return
      end

!**************************************************
!
!	Name: paw_proj_check_format
!
!	Purpose: 
!
!	Created:	7/30/2002
!**************************************************
      subroutine paw_proj_check_format()
      implicit none

*     **** local variables ****
      logical value,found
      integer ia,ntotal,nkatm
      character*4 element
      character*20 fname

*     **** external functions ****
      logical     nwpw_filefind,paw_proj_format_ok
      integer     ion_nkatm
      external    nwpw_filefind,paw_proj_format_ok
      external    ion_nkatm

      nkatm = ion_nkatm()
      do ia=1,nkatm

*       **** define formatted prj name ****
        call ion_atom_plus_suffix(ia,'.jpp',fname)
      
*       **** make sure prj names are formatted correctly ****
        found = .false.
        do while (.not.found)
          if (nwpw_filefind(fname)) then
             if (paw_proj_format_ok(fname)) found = .true.
          end if

*         **** generate formatted projectors atom.jpp ****
          if (.not.found) then
            call paw_proj_formatter_auto(ia)
          end if
        end do

      end do

      return
      end 

!**************************************************
!	Name: paw_proj_get_nbasis
!
!	Purpose: returns the number of basis functions
!                for each of the formatted projector 
!                files.
!
!	Created:	7/30/2002
!**************************************************
      subroutine paw_proj_get_nbasis(nbasis)
      implicit none
      integer nbasis(*)

*     **** local variables ****
      integer ia,nkatm,ngrid(3)
      real*8  unita(3,3)
      character*20 fname

*     **** external functions ****
      integer  ion_nkatm
      external ion_nkatm

      nkatm = ion_nkatm()
      do ia=1,nkatm
        call ion_atom_plus_suffix(ia,'.jpp',fname)
        call paw_proj_read_header(fname,ngrid,nbasis(ia),unita)
      end do

      return
      end


!**************************************************
!
!	Name: paw_proj_format_ok
!
!	Purpose: returns true if header of the formatted 
!                projector file agrees with control.
!
!	Created:	7/30/2002
!**************************************************
      logical function paw_proj_format_ok(fname)
      implicit none
      character*(*) fname


*     **** local variables ****
      logical correct_box
      integer ngrid(3),nbasis
      real*8  unita(3,3)

*     **** external functions ****
      integer  control_ngrid
      real*8   control_unita
      external control_ngrid
      external control_unita

      correct_box = .true.
      call paw_proj_read_header(fname,ngrid,nbasis,unita)
      if ( (ngrid(1).ne.control_ngrid(1)) .or.
     >     (ngrid(2).ne.control_ngrid(2)) .or.
     >     (ngrid(3).ne.control_ngrid(3)) .or.
     >     (unita(1,1).ne.control_unita(1,1)) .or.
     >     (unita(2,1).ne.control_unita(2,1)) .or.
     >     (unita(3,1).ne.control_unita(3,1)) .or.
     >     (unita(1,2).ne.control_unita(1,2)) .or.
     >     (unita(2,2).ne.control_unita(2,2)) .or.
     >     (unita(3,2).ne.control_unita(3,2)) .or.
     >     (unita(1,3).ne.control_unita(1,3)) .or.
     >     (unita(2,3).ne.control_unita(2,3)) .or.
     >     (unita(3,3).ne.control_unita(3,3))) then
              correct_box = .false.
           end if

      paw_proj_format_ok = correct_box
      return
      end



!**************************************************
!
!	Name: paw_proj_read_header
!
!	Purpose: read in the header of the formatted 
!                projector file.
!
!	Created:	7/30/2002
!**************************************************
      subroutine paw_proj_read_header(fname,ngrid,nbasis,unita)
      implicit none
      character*(*) fname
      integer ngrid(3),nbasis
      real*8 unita(3,3)

#include "tcgmsg.fh"
#include "msgtypesf.h"

*    *** local variables ***
      integer MASTER,taskid
      parameter(MASTER=0)
      integer msglen,l
      character*255 full_filename

      call Parallel_taskid(taskid)

*     **** open fname binary file ****
      if (taskid.eq.MASTER) then
         call util_file_name_noprefix(fname,.false.,
     >                             .false.,
     >                       full_filename)
         l = index(full_filename,' ') - 1
         call openfile(5,full_filename,l,'r',l)
         call iread(5,ngrid,3)
         call dread(5,unita,9)
         call iread(5,nbasis,1)
         call closefile(5)
      end if

*     **** send header data to all processors ****
      msglen = 3
      call BRDCST(9+MSGINT,ngrid,mitob(msglen),MASTER)
      msglen = 9
      call BRDCST(9+MSGDBL,unita,mdtob(msglen),MASTER)
      msglen = 1
      call BRDCST(9+MSGINT,nbasis,mitob(msglen),MASTER)

      return
      end



!**************************************************
!
!	Name: paw_proj_formatter_auto
!
!	Purpose: read in the header of the formatted 
!                projector file.
!
!	Created:	7/30/2002
!**************************************************
      subroutine paw_proj_formatter_auto(ia)
      implicit none
      integer ia    

#include "mafdecls.fh"

      !*** local variables ***
      integer MASTER,taskid
      parameter (MASTER=0)
      integer nbasis,l,ngrid(3)
      integer i_rgrid,i_prj
      real*8  unita(3,3)
      character*20 jppname,atomname
      character*255 full_filename
    
      !*** external functions ***
      integer  get_i_rgrid,get_i_prj_ps,get_i_nbasis
      real*8   lattice_unita
      external get_i_rgrid,get_i_prj_ps,get_i_nbasis
      external lattice_unita

      call ion_atom_plus_suffix(ia,'_basis',atomname)
      call ion_atom_plus_suffix(ia,'.jpp',jppname)

      write(*,*) "atomname,jppname:",atomname,jppname

      !*** read in projectors from _basis file ***
      i_rgrid  = get_i_rgrid(ia)
      i_prj    = get_i_prj_ps(ia)
      nbasis = int_mb(get_i_nbasis(ia))
      write(*,*) "nbasis:",nbasis
      
      call Parallel_taskid(taskid)
      unita(1,1) = lattice_unita(1,1)
      unita(2,1) = lattice_unita(2,1)
      unita(3,1) = lattice_unita(3,1)
      unita(1,2) = lattice_unita(1,2)
      unita(2,2) = lattice_unita(2,2)
      unita(3,2) = lattice_unita(3,2)
      unita(1,3) = lattice_unita(1,3)
      unita(2,3) = lattice_unita(2,3)
      unita(3,3) = lattice_unita(3,3)
      call D3dB_nx(1,ngrid(1))
      call D3dB_ny(1,ngrid(2))
      call D3dB_nz(1,ngrid(3))


*     **** open jppname binary file ****
      if (taskid.eq.MASTER) then
         call util_file_name_noprefix(jppname,.false.,
     >                             .false.,
     >                       full_filename)
         l = index(full_filename,' ') - 1
         call openfile(6,full_filename,l,'w',l)
         call iwrite(6,ngrid,3)
         call dwrite(6,unita,9)
         call iwrite(6,nbasis,1)
         call closefile(6)
      end if

      !stop

      return
      end
