*
* $Id$
*

*     *******************************************
*     *                                         *
*     *            lcao_write_psi               *
*     *                                         *
*     *******************************************
      subroutine lcao_write_psi(filename,
     >                          version,
     >                          ngrid,
     >                          unita,
     >                          ispin,ne,
     >                          psimatrix,phi)
      implicit none
      character*(*) filename
      integer      version
      integer      ngrid(3)
      real*8       unita(3,3)
      integer      ispin,ne(2)
      integer      psimatrix(*)
      complex*16    phi(*)

#include "bafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "errquit.fh"

*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)
      logical value,oprint
      integer n,l,nbasis,ms,nfft3d,occupation
      integer psi(2),tmp(2),coef(2)
      character*255 full_filename


*     ***** external functions ****
      logical  control_print
      external control_print
      integer  aorbs_nbasis
      external aorbs_nbasis

      call Parallel_taskid(taskid)
      oprint = (taskid.eq.MASTER).and.control_print(print_medium)
      call D3dB_nfft3d(1,nfft3d)
      nbasis = aorbs_nbasis()

      value = BA_push_get(mt_dcpl,nfft3d,'psi',psi(2),psi(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,(nfft3d),'tmp',tmp(2),tmp(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(nbasis),'coef',coef(2),coef(1))
      if (.not. value)
     >  call errquit('lcao_write_psi:out of stack memory',0, MA_ERR)


*     **** open ELCIN binary file ****
      if (oprint) then
         write(*,*)
         write(*,*) "lcao guess, initial psi:",filename
         write(*,*) "-   spin, nalpha, nbeta:",ispin,ne
      end if

      if (taskid.eq.MASTER) then
         call util_file_name_noprefix(filename,.false.,
     >                                .false.,
     >                       full_filename)
         l = index(full_filename,' ') -1
         call openfile(6,full_filename,l,'w',l)
         call iwrite(6,version,1)
         call iwrite(6,ngrid,3)
         call dwrite(6,unita,9)
         call iwrite(6,ispin,1)
         call iwrite(6,ne,2)
         occupation = -1
         call iwrite(6,occupation,1)
      end if

*     *** write out 3d blocks ***
      do ms=1,ispin
      do n=1,ne(ms)

*        **** calculate electronic wavefunction ****
         call ga_get(psimatrix(ms),1,nbasis,n,n,
     >               dbl_mb(coef(1)),nbasis)
         call ga_sync()
         call lcao_make_psi(dbl_mb(coef(1)),
     >                      dcpl_mb(psi(1)),
     >                      phi)

         call Pack_c_unpack(1,dcpl_mb(psi(1)))
         call D3dB_c_write(1,6,dcpl_mb(psi(1)),
     >                         dcpl_mb(tmp(1)),0)
      end do
      end do

*     *** close ELCIN binary file ***
      if (taskid.eq.MASTER) then
        call closefile(6)
      end if

      value = BA_pop_stack(coef(2))
      value = value.and.
     >        BA_pop_stack(tmp(2))
      value = value.and.
     >        BA_pop_stack(psi(2))
      if (.not. value)
     >  call errquit('lcao_write_psi:popping stack memory',0, MA_ERR)
      return
      end



*     *******************************************
*     *                                         *
*     *            lcao_write_epsi              *
*     *                                         *
*     *******************************************

      subroutine lcao_write_epsi(filename,
     >                          version,
     >                          ngrid,
     >                          unita,
     >                          ispin,ne,ne_excited,
     >                          psimatrix,phi)
      implicit none
      character*(*) filename
      integer      version
      integer      ngrid(3)
      real*8       unita(3,3)
      integer      ispin,ne(2),ne_excited(2)
      integer      psimatrix(*)
      complex*16    phi(*)

#include "bafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "errquit.fh"

*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)
      logical value,oprint
      integer n,l,nbasis,ms,nfft3d,occupation
      integer psi(2),tmp(2),coef(2)
      character*255 full_filename


*     ***** external functions ****
      logical  control_print
      external control_print
      integer  aorbs_nbasis
      external aorbs_nbasis

      call Parallel_taskid(taskid)
      oprint = (taskid.eq.MASTER).and.control_print(print_medium)
      call D3dB_nfft3d(1,nfft3d)
      nbasis = aorbs_nbasis()

      value = BA_push_get(mt_dcpl,nfft3d,'psi',psi(2),psi(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,(nfft3d),'tmp',tmp(2),tmp(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(nbasis),'coef',coef(2),coef(1))
      if (.not. value)
     >  call errquit('lcao_write_epsi:out of stack memory',0,0)


*     **** open ELCIN binary file ****
      if (oprint) then
         write(*,*)
         write(*,*) "lcao guess, initial epsi:",filename
         write(*,*) "-   spin, nalpha, nbeta:",ispin,ne_excited
      end if

      if (taskid.eq.MASTER) then
         call util_file_name_noprefix(filename,.false.,
     >                                .false.,
     >                       full_filename)
         l = index(full_filename,' ') -1
         call openfile(6,full_filename,l,'w',l)
         call iwrite(6,version,1)
         call iwrite(6,ngrid,3)
         call dwrite(6,unita,9)
         call iwrite(6,ispin,1)
         call iwrite(6,ne_excited,2)
         occupation = -1
         call iwrite(6,occupation,1)
      end if

*     *** write out 3d blocks ***
      do ms=1,ispin
      do n=ne(ms)+1,ne(ms)+ne_excited(ms)

*        **** calculate electronic wavefunction ****
         call ga_get(psimatrix(ms),1,nbasis,n,n,
     >               dbl_mb(coef(1)),nbasis)
         call lcao_make_psi(dbl_mb(coef(1)),
     >                      dcpl_mb(psi(1)),
     >                      phi)

         call Pack_c_unpack(1,dcpl_mb(psi(1)))
         call D3dB_c_write(1,6,dcpl_mb(psi(1)),
     >                         dcpl_mb(tmp(1)),0)
      end do
      end do

*     *** close ELCIN binary file ***
      if (taskid.eq.MASTER) then
        call closefile(6)
      end if

      value =           BA_pop_stack(coef(2))
      value = value.and.BA_pop_stack(tmp(2))
      value = value.and.BA_pop_stack(psi(2))
      if (.not. value)
     >  call errquit('lcao_write_epsi:popping stack memory',0,0)
      return
      end

