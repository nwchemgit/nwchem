*
* $Id: paw_proj.F,v 1.26 2003-03-05 23:16:32 bylaska Exp $
*
* $Log: not supported by cvs2svn $
* Revision 1.25  2003/02/13 01:58:55  bylaska
* ...EJB
*
* Revision 1.24  2003/02/11 01:41:24  edo
* eliminated f90-isms
*
* Revision 1.23  2003/02/09 21:27:07  marat
* commented out some of the write statements
* MV
*
* Revision 1.22  2003/02/06 06:13:09  marat
* moved legendre and spher. harm. stuff to paw_special_functions dir
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
#include "paw_proj_data.fh" 
#include "paw_geom.fh" 
 
      integer paw_proj_nbasis
      external paw_proj_nbasis

*     **** local variables ****
      logical value,found
      integer ia,npack1
      character*4 element
      character*20 fname


*     *** get number of diff atom types ***
      prj_nkatm = ion_nkatm()

*     *** number of basis functions for each atom kind ***
      value = MA_alloc_get(mt_int,(prj_nkatm),
     >                     'prj_nbasis',prj_nbasis(2),prj_nbasis(1))
      value = value.and.
     >        MA_alloc_get(mt_int,(prj_nkatm),
     >                     'prj_indx',prj_indx(2),prj_indx(1))
      value = value.and.
     >        MA_alloc_get(mt_int,(2*prj_nkatm),
     >                     'i_prj_l',i_prj_l(2),i_prj_l(1))
      value = value.and.
     >        MA_alloc_get(mt_int,(2*prj_nkatm),
     >                     'i_prj_m',i_prj_m(2),i_prj_m(1))
      if (.not.value) call errquit('paw_proj_init: alloc heap',0)


*     *** generate formatted projector file if needed ***
      call paw_proj_check_format()

*     *** read number of (nlm) projectors for each atom
*     from headers off the formatted files
      call paw_proj_nbasis_read(int_mb(prj_nbasis(1)))

*     *** total number of projectors for diff types of atoms ***
      prj_ntot = 0
      do ia=1,prj_nkatm
        prj_ntot = prj_ntot + int_mb(prj_nbasis(1)+ia-1)
      end do
c      write(*,*) "paw_proj:prj_ntot",prj_ntot

      do ia=1,prj_nkatm
        value = MA_alloc_get(mt_int,(int_mb(prj_nbasis(1)+ia-1)),
     >                       'sub_i_prj_l',
     >                       int_mb(i_prj_l(1)+2*(ia-1)+1),
     >                       int_mb(i_prj_l(1)+2*(ia-1)))
        value = value.and.
     >          MA_alloc_get(mt_int,(int_mb(prj_nbasis(1)+ia-1)),
     >                       'sub_i_prj_m',
     >                       int_mb(i_prj_m(1)+2*(ia-1)+1),
     >                       int_mb(i_prj_m(1)+2*(ia-1)))

      end do

c     *** allocate kspace arrays***
      call Pack_npack(1,npack1)
      value = MA_alloc_get(mt_dcpl,(prj_ntot*npack1),
     >                     'prj',prj(2),prj(1))
      if (.not.value) call errquit('paw_proj_init: alloc heap',1)

*     *** set prj_indx ***
      int_mb(prj_indx(1)) = 0
      do ia=2,prj_nkatm
         int_mb(prj_indx(1)+ia-1) = int_mb(prj_indx(1)+ia-2) 
     >                            + int_mb(prj_nbasis(1)  +ia-2)*npack1
      end do

*     *** read in formatted prj's to prj common block ***
      do ia=1,prj_nkatm
         call paw_proj_read(ia,
     >            int_mb(int_mb(i_prj_l(1)+2*(ia-1))),
     >            int_mb(int_mb(i_prj_m(1)+2*(ia-1))),
     >            npack1,
     >            dcpl_mb(prj(1)+int_mb(prj_indx(1)+ia-1)))
      end do

     
      return
      end

!**************************************************
!
!	Name: paw_proj_i_prj
!
!	Purpose: returns the dcpl_mb ma index of
!                the paw projectors.
!
!	Created:	7/30/2002
!**************************************************
      integer function paw_proj_i_prj()
      implicit none
      
#include "paw_proj_data.fh"

      paw_proj_i_prj = prj(1)
      return
      end

!**************************************************
!
!	Name: paw_proj_i_prj_atom
!
!	Purpose: returns the dcpl_mb ma index of
!                the paw projectors.
!
!	Created:	7/30/2002
!**************************************************
      integer function paw_proj_i_prj_atom(ia)
      implicit none
      integer ia
      
#include "mafdecls.fh"
#include "paw_proj_data.fh"

      paw_proj_i_prj_atom = prj(1) + int_mb(prj_indx(1)+ia-1)
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
#include "paw_proj_data.fh"

      paw_proj_nbasis = int_mb(prj_nbasis(1)+ia-1)
      return
      end


!**************************************************
!
!	Name: paw_proj_total_nbasis
!
!	Purpose: returns the number of the paw projectors
!                for this kind of atom.
!
!	Created:	7/30/2002
!**************************************************
c      integer function paw_proj_total_nbasis()
c      implicit none
c      integer ia
      
c#include "mafdecls.fh"
c#include "paw_proj_data.fh"
c
c      paw_proj_total_nbasis = prj_ntot
c      return
c      end




!**************************************************
!
!	Name: paw_proj_l
!
!	Purpose: returns the orbital quantum number
!                for the nth projector of the iath 
!                kind of atom.
!
!	Created:	7/30/2002
!**************************************************
      integer function paw_proj_l(n,ia)
      implicit none
      integer n,ia
      
#include "mafdecls.fh"
#include "paw_proj_data.fh"

      paw_proj_l = int_mb(int_mb(i_prj_l(1)+2*(ia-1))+n-1)
      return
      end

!**************************************************
!
!	Name: paw_proj_m
!
!	Purpose: returns the magnetic quantum number
!                for the nth projector of the iath 
!                kind of atom.
!
!	Created:	7/30/2002
!**************************************************
      integer function paw_proj_m(n,ia)
      implicit none
      integer n,ia
      
#include "mafdecls.fh"
#include "paw_proj_data.fh"


      paw_proj_m = int_mb(int_mb(i_prj_m(1)+2*(ia-1))+n-1)
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

#include "paw_proj_data.fh"

*     **** local variables ****
      logical value,found
      integer ia,nkatm
      character*4 element
      character*20 fname

*     **** external functions ****
      logical     nwpw_filefind,paw_proj_format_ok
      external    nwpw_filefind,paw_proj_format_ok

      do ia=1,prj_nkatm

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
!	Name: paw_proj_nbasis_read
!
!	Purpose: returns the number of basis functions
!                for each of the formatted projector 
!                files.
!
!	Created:	7/30/2002
!**************************************************
      subroutine paw_proj_nbasis_read(nbasis)
      implicit none
      integer nbasis(*)

#include "paw_proj_data.fh"

*     **** local variables ****
      integer ia,ngrid(3)
      real*8  unita(3,3)
      character*20 fname


      do ia=1,prj_nkatm
        call ion_atom_plus_suffix(ia,'.jpp',fname)
        call paw_proj_read_header(fname,ngrid,nbasis(ia),unita)
      end do

      return
      end



!**************************************************
!	Name: paw_proj_paw_basis_tot_nbasis
!
!	Purpose: returns the number of basis functions
!                for each of the formatted projector 
!                files.
!
!	Created:	7/30/2002
!**************************************************
      integer function paw_proj_tot_nbasis()
      implicit none

#include "paw_proj_data.fh"

      external paw_proj_nbasis
      integer paw_proj_nbasis

*     **** local variables ****
      integer ia
      integer tot_nbasis

*     **** external functions ****
      integer  ion_natm
      external ion_natm

c     !*** calculate total number of (n,l,m) projectors  ***
      tot_nbasis = 0
c     !-- loop over diff kinds of atoms --
      do ia=1,prj_nkatm
c        !-- ion_natm(ia) is number of atoms of kind ia --
         tot_nbasis = tot_nbasis
     >                + paw_proj_nbasis(ia)*ion_natm(ia)
      end do

      paw_proj_tot_nbasis = tot_nbasis
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
!	Name: paw_proj_read
!
!	Purpose: read in the formatted 
!                projector file.
!
!	Created:	8/06/2002
!**************************************************
      subroutine paw_proj_read(ia,proj_l,proj_m,npack1,prj)
      implicit none
      integer ia
      integer proj_l(*),proj_m(*)
      integer    npack1
      complex*16 prj(npack1,*)

#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"

*    *** local variables ***
      integer MASTER,taskid
      parameter(MASTER=0)

      logical value
      integer ii,msglen,l,nbasis,ngrid(3),nfft3d
      real*8  unita(3,3)
      character*20  fname
      character*255 full_filename
      integer tmp1(2),tmp2(2)

      call Parallel_taskid(taskid)

*     **** open fname binary file ****
      call ion_atom_plus_suffix(ia,'.jpp',fname)
      if (taskid.eq.MASTER) then

         call util_file_name_noprefix(fname,.false.,
     >                             .false.,
     >                       full_filename)
         l = index(full_filename,' ') - 1
         call openfile(5,full_filename,l,'r',l)
         call iread(5,ngrid,3)
         call dread(5,unita,9)
         call iread(5,nbasis,1)
      end if

*     **** send header data to all processors ****
      msglen = 3
      call BRDCST(9+MSGINT,ngrid,mitob(msglen),MASTER)
      msglen = 9
      call BRDCST(9+MSGDBL,unita,mdtob(msglen),MASTER)
      msglen = 1
      call BRDCST(9+MSGINT,nbasis,mitob(msglen),MASTER)
      
      call D3dB_nfft3d(1,nfft3d)
      value = MA_push_get(mt_dcpl,nfft3d,'tmp1',tmp1(2),tmp1(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,nfft3d,'tmp2',tmp2(2),tmp2(1))
      if (.not.value) call errquit('paw_proj_read: push stack',0)

      do ii=1,nbasis
         if (taskid.eq.MASTER) then
            call iread(5,proj_l(ii),1)
            call iread(5,proj_m(ii),1)
         end if
         call BRDCST(9+MSGINT,proj_l(ii),mitob(msglen),MASTER) 
         call BRDCST(9+MSGINT,proj_m(ii),mitob(msglen),MASTER)
         call D3dB_c_read(1,5,dcpl_mb(tmp1(1)),dcpl_mb(tmp2(1)))
         call Pack_c_pack(1,dcpl_mb(tmp1(1)))
         call Pack_c_Copy(1,dcpl_mb(tmp1(1)),prj(1,ii))
      end do
      value =           MA_pop_stack(tmp2(2))
      value = value.and.MA_pop_stack(tmp1(2))
      if (.not.value) call errquit('paw_proj_read: pop stack',1)

      if (taskid.eq.MASTER) then
        call closefile(5)
      end if

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
#include "bessel_transform.fh"
#include "paw_params.fh"
#include "paw_basis.fh"


      !*** local variables ***
      integer MASTER,taskid
      parameter (MASTER=0)

      real*8  small
      parameter (small=1.0d-9)

      logical value
      integer ii,i,npack1,nfft3d,nr,basis_nbasis,l,m,nbasis,ngrid(3),jj
      integer i_rgrid,i_prj
      real*8  unita(3,3),gg,log_amesh
      integer tmp(2),prj(2),rayleigh(2),Gx(2),Gy(2),Gz(2)

      character*20 jppname,atomname
      character*255 full_filename

      !*** external functions ***
      integer  G_indx
      real*8   lattice_unita,lattice_omega
      external G_indx
      external lattice_unita,lattice_omega


      call ion_atom_plus_suffix(ia,'_basis',atomname)
      call ion_atom_plus_suffix(ia,'.jpp',jppname)


      !*** read in projectors from _basis file ***
      i_rgrid   = paw_basis_i_rgrid(ia)
      i_prj     = paw_basis_i_prj_ps(ia)
      basis_nbasis    = int_mb(paw_basis_i_nbasis(ia))
      nr        = int_mb(paw_basis_i_ngrid(ia))
      log_amesh = dbl_mb(paw_basis_i_log_amesh(ia))

      
      nbasis = 0
      do ii=1,basis_nbasis
         l =  int_mb(paw_basis_i_orb_l(ia)+ii-1)
         nbasis = nbasis + 2*l+1
      end do


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


*     **** open jppname binary file and write header ****
      if (taskid.eq.MASTER) then
         call util_file_name_noprefix(jppname,.false.,
     >                             .false.,
     >                       full_filename)
         l = index(full_filename,' ') - 1
         call openfile(6,full_filename,l,'w',l)
         call iwrite(6,ngrid,3)
         call dwrite(6,unita,9)
         call iwrite(6,nbasis,1)
      end if

*     **** compute bessel transforms ****
      call D3dB_nfft3d(1,nfft3d)
      call Pack_npack(1,npack1)
      value = MA_push_get(mt_dbl,npack1,
     >                    'rayleigh',rayleigh(2),rayleigh(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,nfft3d,
     >                    'prj',prj(2),prj(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,nfft3d,
     >                    'tmp',tmp(2),tmp(1))
      value = value.and.
     >        MA_push_get(mt_dbl,nfft3d,'Gx',Gx(2),Gx(1))
      value = value.and.
     >        MA_push_get(mt_dbl,nfft3d,'Gy',Gy(2),Gy(1))
      value = value.and.
     >        MA_push_get(mt_dbl,nfft3d,'Gz',Gz(2),Gz(1))
      if (.not.value) 
     > call errquit('paw_proj_formatter_auto, push stack',0)

      call dcopy(nfft3d,dbl_mb(G_indx(1)),1,dbl_mb(Gx(1)),1)
      call dcopy(nfft3d,dbl_mb(G_indx(2)),1,dbl_mb(Gy(1)),1)
      call dcopy(nfft3d,dbl_mb(G_indx(3)),1,dbl_mb(Gz(1)),1)
      call Pack_t_pack(1,dbl_mb(Gx(1)))
      call Pack_t_pack(1,dbl_mb(Gy(1)))
      call Pack_t_pack(1,dbl_mb(Gz(1)))

      do ii=1,basis_nbasis
        l =  int_mb(paw_basis_i_orb_l(ia)+ii-1)
        !write(*,*) "formatting prj ",ii," l:",l
        !write(*,*) "nr=",nr
        !write(*,*) "log_amesh=",log_amesh

        do i=1,npack1
          gg  = dsqrt( dbl_mb(Gx(1)+i-1)*dbl_mb(Gx(1)+i-1)
     >               + dbl_mb(Gy(1)+i-1)*dbl_mb(Gy(1)+i-1)
     >               + dbl_mb(Gz(1)+i-1)*dbl_mb(Gz(1)+i-1) )             
  
          if (gg.gt.small) then
            dbl_mb(rayleigh(1)+i-1) 
     >      = spher_bessel_transform(l,nr,log_amesh,
     >                               dbl_mb(i_rgrid),
     >                               dbl_mb(i_prj+(ii-1)*nr),
     >                               gg)
          else
            dbl_mb(rayleigh(1)+i-1) 
     >      = spher_bessel0_transform(l,nr,log_amesh,
     >                                dbl_mb(i_rgrid),
     >                                dbl_mb(i_prj+(ii-1)*nr))
          end if
c        write(*,*) "gg=",gg
c        write(*,*) "bess=",dbl_mb(rayleigh(1)+i-1)

        end do
        do m=-l,l
*          *** generate Ylm ***
           call spher_harmonics_generate(l,m,npack1,
     >                      dbl_mb(Gx(1)),
     >                      dbl_mb(Gy(1)),
     >                      dbl_mb(Gz(1)),
     >                      dcpl_mb(prj(1)))

c           write(*,*) "ylm(1)=",dcpl_mb(prj(1))
c           write(*,*) "bess(1)=",dbl_mb(rayleigh(1))

*          *** multiply Ylm and spherical besse
           call Pack_tc_Mul(1,dbl_mb(rayleigh(1)),
     >                     dcpl_mb(prj(1)),
     >                     dcpl_mb(prj(1)))
c           write(*,*) "ylm(1)*bess(1)=",dcpl_mb(prj(1))
           call rayleigh_itol_scaling(l,npack1,
     >                                dcpl_mb(prj(1)),
     >                                dcpl_mb(prj(1)))
c           write(*,*) "ylm(1)*bess(1)*coeff=",dcpl_mb(prj(1))
c           pause 
           if (taskid.eq.MASTER) then
             call iwrite(6,l,1)
             call iwrite(6,m,1)
           end if
           call Pack_c_unpack(1,dcpl_mb(prj(1)))
           call D3dB_c_write(1,6,dcpl_mb(prj(1)),dcpl_mb(tmp(1)))
        end do ! m loop
      end do ! basis_nbasis loop
      if (taskid.eq.MASTER) then 
         call closefile(6)
         write(6,*)
         write(6,*) " Generated formatted filename:"
         write(6,*) full_filename
      end if
      value =           MA_pop_stack(Gz(2))
      value = value.and.MA_pop_stack(Gy(2))
      value = value.and.MA_pop_stack(Gx(2))
      value = value.and.MA_pop_stack(tmp(2))
      value = value.and.MA_pop_stack(prj(2))
      value = value.and.MA_pop_stack(rayleigh(2))
      if (.not.value) 
     >  call errquit('paw_proj_formatter_auto, pop stack',1)
      
      return
      end

!**************************************************
!
!	Name: rayleigh_itol_scaling
!
!	Purpose: This routine should be rewritten
!                to improve performance.
!  This routine scales Y <-- (-i)**l * Y
!
!	Created:	8/06/2002
!**************************************************
      subroutine rayleigh_itol_scaling(l,nfft3d,Y)
      implicit none
      integer l,nfft3d
      complex*16 Y(*)

      real*8 lattice_omega
      external lattice_omega
#include "paw_params.fh"

      !*** local variables ***
      integer k
      complex*16 coef

      if      (mod(l,4).eq.0) then
         coef = dcmplx(1.0d0,0.0d0)
      else if (mod(l,4).eq.1) then
         coef = dcmplx(0.0d0,-1.0d0)
      else if (mod(l,4).eq.2) then
         coef = dcmplx(-1.0d0,0.0d0)
      else if (mod(l,4).eq.3) then
         coef = dcmplx(0.0d0,1.0d0)
      end if

      coef = coef*fourpi/sqrt(lattice_omega())
      do k=1,nfft3d
        Y(k) = coef*Y(k)
      end do
      return
      end

      subroutine paw_proj_weghts_set()

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
#include "paw_proj_data.fh"

      !*** local variables ***
      logical value
      integer ia
      
      value = .true.
      do ia=1,prj_nkatm
         value = value.and.MA_free_heap(int_mb(i_prj_l(1)+2*(ia-1)+1))
         value = value.and.MA_free_heap(int_mb(i_prj_m(1)+2*(ia-1)+1))
      end do
      value = value.and.MA_free_heap(i_prj_l(2))
      value = value.and.MA_free_heap(i_prj_m(2))
      value = value.and.MA_free_heap(prj_indx(2))
      value = value.and.MA_free_heap(prj_nbasis(2))
      value = value.and.MA_free_heap(prj(2))
      if (.not.value) call errquit('paw_proj_end: dealloc heap',0)
      return
      end



