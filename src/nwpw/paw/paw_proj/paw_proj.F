*
* $Id: paw_proj.F,v 1.11 2002-09-12 01:01:44 marat Exp $
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
      integer ia,npack1
      character*4 element
      character*20 fname

*     **** external functions ****
      logical     nwpw_filefind
      integer     ion_nkatm
      character*2 ion_atom
      external    nwpw_filefind
      external    ion_nkatm
      external    ion_atom

*     *** get number of diff atom types ***
      nkatm = ion_nkatm()

*     *** number of basis functions for each atom kind ***
      value = MA_alloc_get(mt_int,(nkatm),
     >                     'nbasis',nbasis(2),nbasis(1))
      value = value.and.
     >        MA_alloc_get(mt_int,(nkatm),
     >                     'prj_indx',prj_indx(2),prj_indx(1))
      value = value.and.
     >        MA_alloc_get(mt_int,(2*nkatm),
     >                     'i_prj_l',i_prj_l(2),i_prj_l(1))
      value = value.and.
     >        MA_alloc_get(mt_int,(2*nkatm),
     >                     'i_prj_m',i_prj_m(2),i_prj_m(1))
      if (.not.value) call errquit('paw_proj_init: alloc heap',0)



*     *** make sure prj files are formatted correctly ***
*     *** generate them if needed ***
      call paw_proj_check_format()

*     *** get the size of prj common block and allocate it ***
      call paw_proj_get_nbasis(int_mb(nbasis(1)))
      ntotal = 0
      do ia=1,nkatm
        ntotal = ntotal + int_mb(nbasis(1)+ia-1)
        value = MA_alloc_get(mt_int,(int_mb(nbasis(1)+ia-1)),
     >                       'sub_i_prj_l',
     >                       int_mb(i_prj_l(1)+2*(ia-1)+1),
     >                       int_mb(i_prj_l(1)+2*(ia-1)))
        value = value.and.
     >          MA_alloc_get(mt_int,(int_mb(nbasis(1)+ia-1)),
     >                       'sub_i_prj_m',
     >                       int_mb(i_prj_m(1)+2*(ia-1)+1),
     >                       int_mb(i_prj_m(1)+2*(ia-1)))

      end do
      call Pack_npack(1,npack1)
      value = MA_alloc_get(mt_dcpl,(ntotal*npack1),
     >                     'prj',prj(2),prj(1))
      if (.not.value) call errquit('paw_proj_init: alloc heap',1)

*     *** set prj_indx ***
      int_mb(prj_indx(1)) = 0
      do ia=2,nkatm
         int_mb(prj_indx(1)+ia-1) = int_mb(prj_indx(1)+ia-2) 
     >                            + int_mb(nbasis(1)  +ia-2)*npack1
      end do

*     *** read in formatted prj's to prj common block ***
      do ia=1,nkatm
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

      !*** local variables ***
      logical value
      integer ia
      
      value = .true.
      do ia=1,nkatm
         value = value.and.MA_free_heap(int_mb(i_prj_l(1)+2*(ia-1)+1))
         value = value.and.MA_free_heap(int_mb(i_prj_m(1)+2*(ia-1)+1))
      end do
      value = value.and.MA_free_heap(i_prj_l(2))
      value = value.and.MA_free_heap(i_prj_m(2))
      value = value.and.MA_free_heap(prj_indx(2))
      value = value.and.MA_free_heap(nbasis(2))
      value = value.and.MA_free_heap(prj(2))
      if (.not.value) call errquit('paw_proj_end: dealloc heap',0)
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
      
#include "paw_proj.fh"

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
#include "paw_proj.fh"

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
#include "paw_proj.fh"

      paw_proj_nbasis = int_mb(nbasis(1)+ia-1)
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
      integer function paw_proj_total_nbasis()
      implicit none
      integer ia
      
#include "mafdecls.fh"
#include "paw_proj.fh"

      paw_proj_total_nbasis = ntotal
      return
      end




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
#include "paw_proj.fh"

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
#include "paw_proj.fh"


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

*     **** local variables ****
      logical value,found
      integer ia,nkatm
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
      value =           MA_push_get(mt_dcpl,nfft3d,'tmp1',tmp1(2),tmp1(1))
      value = value.and.MA_push_get(mt_dcpl,nfft3d,'tmp2',tmp2(2),tmp2(1))
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
#include "paw_basis.fh"

      !*** external functions ***
      external G_indx
      external lattice_unita
      integer  G_indx
      real*8   lattice_unita

      !*** local variables ***
      integer MASTER,taskid
      parameter (MASTER=0)

      real*8  small
      parameter (small=1.0d-9)

      logical value
      integer ii,i,npack1,nfft3d,nr,basis_nbasis,l,m,nbasis,ngrid(3)
      integer i_rgrid,i_prj
      real*8  unita(3,3),gg,log_amesh
      integer tmp(2),prj(2),rayleigh(2),Gx(2),Gy(2),Gz(2)

      character*20 jppname,atomname
      character*255 full_filename


      call ion_atom_plus_suffix(ia,'_basis',atomname)
      call ion_atom_plus_suffix(ia,'.jpp',jppname)

      write(*,*) "atomname,jppname:",atomname,jppname

      !*** read in projectors from _basis file ***
      i_rgrid   = get_i_rgrid(ia)
      i_prj     = get_i_prj_ps(ia)
      basis_nbasis    = int_mb(get_i_nbasis(ia))
      nr        = int_mb(get_i_ngrid(ia))
      log_amesh = dbl_mb(get_i_log_amesh(ia))
      write(*,*) "basis_nbasis:",basis_nbasis," ngrid:",nr
      
      nbasis = 0
      do ii=1,basis_nbasis
         l =  int_mb(get_i_orb_l(ia)+ii-1)
         nbasis = nbasis + 2*l+1
      end do
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
        l =  int_mb(get_i_orb_l(ia)+ii-1)
        write(*,*) "formatting prj ",ii," l:",l
     
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
        end do
        do m=-l,l
           call spher_harmonics_generate(l,m,npack1,
     >                      dbl_mb(Gx(1)),
     >                      dbl_mb(Gy(1)),
     >                      dbl_mb(Gz(1)),
     >                      dcpl_mb(prj(1)))
           call Pack_tc_Mul(1,dbl_mb(rayleigh(1)),
     >                     dcpl_mb(prj(1)),
     >                     dcpl_mb(prj(1)))
           call rayleigh_itol_scaling(l,npack1,
     >                                dcpl_mb(prj(1)),
     >                                dcpl_mb(prj(1)))
        
           if (taskid.eq.MASTER) then
             write(*,*) "writing formatted prj:",ii,"  l,m",l,m
             call iwrite(6,l,1)
             call iwrite(6,m,1)
           end if
           call Pack_c_unpack(1,dcpl_mb(prj(1)))
           call D3dB_c_write(1,6,dcpl_mb(prj(1)),dcpl_mb(tmp(1)))
        end do
      end do
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

      do k=1,nfft3d
        Y(k) = coef*Y(k)
      end do
      return
      end

!**************************************************
!
!	Name: spher_harmonics_generate
!
!	Purpose: 
!
!	Created:	8/05/2002
!**************************************************
      subroutine spher_harmonics_generate(l,m,nfft3d,
     >                                    Gx,Gy,Gz,
     >                                    Ylm)
      implicit none
      integer    l,m,nfft3d
      real*8     Gx(*),Gy(*),Gz(*)
      complex*16 Ylm(*)

*     **** local variables ****
      real*8 small
      parameter (small=1.0d-9)

      integer i,k,mabs
      real*8     g,gxy,cs_theta,sgn,coeff,fourpi
      complex*16 exp_phi,coef2

*     **** external functions ****
      real*8   legendre_lm
      external legendre_lm


      fourpi = 16.0d0*datan(1.0d0)
      sgn    = 1.0d0
      if (m.ge.0) then
        mabs = m
      else
        mabs = -m
        if (mod(mabs,2).eq.1) sgn = -1.0d0
      end if
     

      !*** find coefficient ***
      if(mabs.eq.0) then
         coeff = 1.0d0
      else
         coeff = 1.0d0
         do i=1,2*mabs
           coeff = coeff/(l-mabs+i)
         end do
      end if
      coeff = coeff*(2*l+1.0d0)/fourpi
      coeff = sgn*dsqrt(coeff)

      if (m.ge.0) then
        do k=1,nfft3d
          gxy = Gx(k)**2 + Gy(k)**2
          g   = gxy + Gz(k)**2
          gxy = dsqrt(gxy)
          g   = dsqrt(g) 

          if (gxy.gt.small)then
            exp_phi = (dcmplx(Gx(k),Gy(k))/gxy)**mabs 
          else
            exp_phi = dcmplx(1.0d0,0.0d0)
          end if 

          if (g.gt.small) then
            cs_theta = Gz(k)/g
          else
            cs_theta = 0.0d0
          end if 

          Ylm(k) = coeff*legendre_lm(l,mabs,cs_theta)
     >                  *exp_phi
        end do
      else
        do k=1,nfft3d
          gxy = Gx(k)**2 + Gy(k)**2
          g   = gxy + Gz(k)**2
          gxy = dsqrt(gxy)
          g   = dsqrt(g) 

          if (gxy.gt.small)then
            exp_phi = (dcmplx(Gx(k),Gy(k))/gxy)**mabs 
          else
            exp_phi = dcmplx(1.0d0,0.0d0)
          end if 

          if (g.gt.small) then
            cs_theta = Gz(k)/g
          else
            cs_theta = 0.0d0
          end if 
       
          Ylm(k) = coeff*legendre_lm(l,mabs,cs_theta)
     >                  *dconjg(exp_phi)
        end do
      end if
     
      return
      end

!**************************************************
!
!	Name: legendre_lm
!
!	Purpose: 
!
!	Created:	8/05/2002
!**************************************************
    
      real*8 function legendre_lm(l,m,x)
      implicit none
      integer l,m
      real*8 x

*     *** local variables ****
      integer i,ll
      real*8  fact,pll,pmm,pmmp1,somx2,Plm

      if ((m.lt.0).or.(m.gt.l).or.(dabs(x).gt.1.0d0))
     >  call errquit('legendre_lm: bad argments',0)

      pmm = 1.0d0
      if (m.gt.0) then
        somx2 = dsqrt((1.0d0-x)*(1.0d0+x))
        fact  = 1.0d0
        do i=1,m
          pmm  = -pmm*fact*somx2
          fact = fact + 2.0d0
        end do
      end if

      if (l.eq.m) then
        Plm = pmm
      else
         pmmp1 = x*(2*m+1.0d0)*pmm
         if (l.eq.(m+1)) then
           Plm = pmmp1
         else
            do ll=m+2,l
              pll   = (x*(2*ll-1.0d0)*pmmp1-(ll+m-1.0d0)*pmm)/dble(ll-m)
              pmm   = pmmp1
              pmmp1 = pll
            end do
            Plm = pll
         end if
      end if

      legendre_lm = Plm
      return
      end
