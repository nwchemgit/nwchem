
*
* $Id$
*

#define TCGMSG

*     ***********************************
*     *					*
*     *	 	     aorbs_init		*
*     *					*
*     ***********************************
      logical function aorbs_init()      
      implicit none

#include "bafdecls.fh"

***** AORBS common block ****
#include "aorbs.fh"

*     **** local variables ****
      integer npack1,nion
      logical value

*     **** external functions *****
      integer  ion_nkatm_qm,ion_nion_qm
      external ion_nkatm_qm,ion_nion_qm
 

      call Pack_npack(1,npack1)
      nion   = ion_nion_qm()
      nkatmx = ion_nkatm_qm()

      value = BA_alloc_get(mt_dbl,(norbs_max*nkatmx*npack1),
     >                    'aorbs',aorbs(2),aorbs(1))
      value = value.and.
     >        BA_alloc_get(mt_int,(nkatmx),'lmmax',lmmax(2),lmmax(1))
      value = value.and.
     >        BA_alloc_get(mt_int,(nkatmx),'lmax',lmax(2),lmax(1))
      value = value.and.
     >        BA_alloc_get(mt_int,(nkatmx),'locp',locp(2),locp(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,(nkatmx),'rcut',rcut(2),rcut(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,(nkatmx),'lmbda',lmbda(2),lmbda(1))

      value = value.and.
     >        BA_alloc_get(mt_int,(nion*norbs_max),
     >                    'lmaorb',lmaorb(2),lmaorb(1))
      value = value.and.
     >        BA_alloc_get(mt_int,(nion*norbs_max),
     >                    'iaaorb',iaaorb(2),iaaorb(1))
      value = value.and.
     >        BA_alloc_get(mt_int,(nion*norbs_max),
     >                    'iiaorb',iiaorb(2),iiaorb(1))
      value = value.and.
     >        BA_alloc_get(mt_int,(nion*norbs_max),
     >                    'basisaorb',basisaorb(2),basisaorb(1))

      if (value) 
     >  call dcopy(norbs_max*nkatmx*npack1,0.0d0,0,dbl_mb(aorbs(1)),1)
      aorbs_init = value
      return
      end


*     ***********************************
*     *					*
*     *	 	     aorbs_end		*
*     *					*
*     ***********************************
      subroutine aorbs_end()      
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"
#include "aorbs.fh"

*     **** local variables ****
      logical value

      value = BA_free_heap(aorbs(2))
      value = value.and.BA_free_heap(lmmax(2))
      value = value.and.BA_free_heap(lmax(2))
      value = value.and.BA_free_heap(locp(2))
      value = value.and.BA_free_heap(rcut(2))
      value = value.and.BA_free_heap(lmbda(2))
      value = value.and.BA_free_heap(lmaorb(2))
      value = value.and.BA_free_heap(iaaorb(2))
      value = value.and.BA_free_heap(iiaorb(2))
      value = value.and.BA_free_heap(basisaorb(2))
      if (.not. value) call errquit('aorbs_end:freeing heap memory',0,
     &       MA_ERR)

      return
      end


*     ***********************************
*     *					*
*     *	 	     aorbs_norbs	*
*     *					*
*     ***********************************

      integer function aorbs_norbs(ia)
      implicit none
      integer  ia

#include "bafdecls.fh"

***** AORBS common block ****
#include "aorbs.fh"

      aorbs_norbs = int_mb(lmmax(1)+ia-1)
      return
      end


*     ***********************************
*     *					*
*     *	 	   aorbs_nbasis		*
*     *					*
*     ***********************************
      integer function aorbs_nbasis()
      implicit none

***** AORBS common block ****
#include "aorbs.fh"

      aorbs_nbasis = ibasis
      return
      end

*     ***********************************
*     *                                 *
*     *              aorbs_lmax         *
*     *                                 *
*     ***********************************

      integer function aorbs_lmax(ia)
      implicit none
      integer  ia

#include "bafdecls.fh"

***** AORBS common block ****
#include "aorbs.fh"

      aorbs_lmax = int_mb(lmax(1)+ia-1)-1
      return
      end


*     ***********************************
*     *                                 *
*     *              aorbs_rcut         *
*     *                                 *
*     ***********************************
      real*8 function aorbs_rcut(ia)
      implicit none
      integer  ia

#include "bafdecls.fh"
#include "aorbs.fh"

      aorbs_rcut = dbl_mb(rcut(1)+ia-1)
      return
      end

*     ***********************************
*     *                                 *
*     *           aorbs_lmbda           *
*     *                                 *
*     ***********************************
      real*8 function aorbs_lmbda(ia)
      implicit none
      integer  ia

#include "bafdecls.fh"
#include "aorbs.fh"

      aorbs_lmbda = dbl_mb(lmbda(1)+ia-1)
      return
      end


*     ***********************************
*     *                                 *
*     *              aorbs_l            *
*     *                                 *
*     ***********************************

      integer function aorbs_l(ia,n)
      implicit none
      integer  ia
      integer n          ! basis number

#include "bafdecls.fh"

***** AORBS common block ****
#include "aorbs.fh"

*     *** local variables ***
      integer l,m,lm
      
      lm = int_mb(lmaorb(1)+n-1)
      l = 0
      if (lm.eq.1) l = 0 

      if (lm.eq.2) l = 1
      if (lm.eq.3) l = 1
      if (lm.eq.4) l = 1

      if (lm.eq.5) l = 2
      if (lm.eq.6) l = 2
      if (lm.eq.7) l = 2
      if (lm.eq.8) l = 2
      if (lm.eq.9) l = 2

      if (lm.eq.10) l = 3
      if (lm.eq.11) l = 3
      if (lm.eq.12) l = 3
      if (lm.eq.13) l = 3
      if (lm.eq.14) l = 3
      if (lm.eq.15) l = 3
      if (lm.eq.16) l = 3

      aorbs_l = l
      return
      end


*     ***********************************
*     *                                 *
*     *              aorbs_m            *
*     *                                 *
*     ***********************************

      integer function aorbs_m(ia,n)
      implicit none
      integer  ia
      integer n          ! basis number

#include "bafdecls.fh"

***** AORBS common block ****
#include "aorbs.fh"

*     *** local variables ***
      integer l,m,lm

      lm = int_mb(lmaorb(1)+n-1)
      m = 0
      if (lm.eq.1) m = 0

      if (lm.eq.2) m = -1
      if (lm.eq.3) m =  0
      if (lm.eq.4) m =  1

      if (lm.eq.5) m = -2
      if (lm.eq.6) m = -1
      if (lm.eq.7) m =  0
      if (lm.eq.8) m =  1
      if (lm.eq.9) m =  2

      if (lm.eq.10) m = -3
      if (lm.eq.11) m = -2
      if (lm.eq.12) m = -1
      if (lm.eq.13) m =  0
      if (lm.eq.14) m =  1
      if (lm.eq.15) m =  2
      if (lm.eq.16) m =  3

      aorbs_m = m
      return
      end




*     ***********************************
*     *                                 *
*     *           aorbs_katm            *
*     *                                 *
*     ***********************************
      integer function aorbs_katm(n)
      implicit none
      integer n

#include "bafdecls.fh"
#include "aorbs.fh"

      aorbs_katm = int_mb(iaaorb(1)+n-1)
      return
      end


*     ***********************************
*     *                                 *
*     *           aorbs_ii              *
*     *                                 *
*     ***********************************
      integer function aorbs_ii(n)
      implicit none
      integer n

#include "bafdecls.fh"
#include "aorbs.fh"

      aorbs_ii = int_mb(iiaorb(1)+n-1)
      return
      end



*     ***********************************
*     *                                 *
*     *     aorbs_get_basis_number      *
*     *                                 *
*     ***********************************

      integer function aorbs_get_basis_number(ii,lm)
      implicit none
      integer ii,lm

#include "bafdecls.fh"
#include "aorbs.fh"

      aorbs_get_basis_number=int_mb(basisaorb(1)+lm-1+(ii-1)*norbs_max)

      return
      end



*     ***********************************
*     *					*
*     *	 	   aorbs_normalize  	*
*     *					*
*     ***********************************
      subroutine aorbs_normalize()
      implicit none

#include "bafdecls.fh"

***** AORBS common block ****
#include "aorbs.fh"

*     **** local variables ****
      integer n,lm,ia,npack1
      real*8  sum

      call Pack_npack(1,npack1)

*     **** Normalize atomic orbitals in k space *****
      do n=1,ibasis
         lm = int_mb(lmaorb(1)+n-1)
         ia = int_mb(iaaorb(1)+n-1)
         call Pack_tt_dot(1,
     >                    dbl_mb(aorbs(1) 
     >                          + (lm-1)*npack1
     >                          + (ia-1)*npack1*norbs_max),
     >                    dbl_mb(aorbs(1) 
     >                          + (lm-1)*npack1
     >                          + (ia-1)*npack1*norbs_max),
     >                    sum)
         sum = 1.0d0/dsqrt(sum)
c         call Pack_t_SMul(1,sum,
c     >                    dbl_mb(aorbs(1) 
c     >                          + (lm-1)*npack1
c     >                          + (ia-1)*npack1*norbs_max),
c     >                    dbl_mb(aorbs(1) 
c     >                          + (lm-1)*npack1
c     >                          + (ia-1)*npack1*norbs_max))
         call Pack_t_SMul1(1,sum,
     >                    dbl_mb(aorbs(1) 
     >                          + (lm-1)*npack1
     >                          + (ia-1)*npack1*norbs_max))
      end do

      return
      end

*     ***********************************
*     *					*
*     *	 	   aorbs_weight		*
*     *					*
*     ***********************************
      real*8 function aorbs_weight(n)
      implicit none
      integer n          ! basis number

#include "bafdecls.fh"

***** AORBS common block ****
#include "aorbs.fh"

*     **** local variables ****
      integer ia
      real*8 zv,zcount

*     **** external functions ****
      real*8   psp_zv
      external psp_zv

      ia     = int_mb(iaaorb(1)+n-1)
      zcount = int_mb(lmmax(1)+ia-1)
      zv     = psp_zv(ia)

      aorbs_weight = (zv/zcount)
      return
      end

*     ***********************************
*     *					*
*     *	 	   aorbs_aorb  		*
*     *					*
*     ***********************************

      subroutine aorbs_aorb(n,aorb)
      implicit none
#include "errquit.fh"
      integer n          ! basis number
      complex*16 aorb(*) ! return orbital

#include "bafdecls.fh"

***** AORBS common block ****
#include "aorbs.fh"

*     **** local variables ****
      logical value
      integer lm,ia,ii
      integer nfft3d,npack1
      integer exi(2)

*     **** external functions ****
      logical  is_sORd
      external is_sORd


      call D3dB_nfft3d(1,nfft3d)
      call Pack_npack(1,npack1)

      value = BA_push_get(mt_dcpl,nfft3d,'exi', exi(2), exi(1))
      if (.not. value) call errquit('aorbs_aorb:out of heap memory',0,
     &       MA_ERR)

*     **** structure factor ****
      lm = int_mb(lmaorb(1)+n-1)
      ia = int_mb(iaaorb(1)+n-1)
      ii = int_mb(iiaorb(1)+n-1)
      call strfac(ii,dcpl_mb(exi(1)))
      call Pack_c_pack(1,dcpl_mb(exi(1)))

*     **** phase factor does not matter therefore ****
*     **** (-i)^l is the same as (i)^l in the     ****
*     **** Rayleigh scattering formula            ****

*     *** current function is s or d ****
      if (is_sORd(lm,int_mb(lmax(1)+ia-1),
     >                   int_mb(locp(1)+ia-1))
     >        ) then
         call Pack_tc_Mul(1,dbl_mb(aorbs(1) 
     >                          + (lm-1)*npack1
     >                          + (ia-1)*npack1*norbs_max),
     >                   dcpl_mb(exi(1)),
     >                   aorb)

*     *** current function is p or f ****
      else
         call Pack_itc_Mul(1,dbl_mb(aorbs(1)
     >                           +(lm-1)*npack1
     >                           +(ia-1)*npack1*norbs_max),
     >                   dcpl_mb(exi(1)),
     >                   aorb)
      end if


      value = BA_pop_stack(exi(2))
      if (.not. value) call errquit('aorbs_aorb:freeing heap memory',0,
     &       MA_ERR)

      return 
      end


*     ***********************************
*     *					*
*     *	 	   aorbs_read 		*
*     *					*
*     ***********************************
      subroutine aorbs_read(fname,
     >                      version,
     >                       nfft,unita,
     >						 atom,
     >                       lmmax,lmax,locp,rcut,lmbda,
     >                       npack1,aorbs,
     >                       tmp,tmp2,
     >                       ierr)
      implicit none 
      character*50 fname
      integer version
      integer nfft(3)
      real*8  unita(3,3)
      character*2 atom
      integer lmmax,lmax,locp
      real*8 rcut,lmbda
      integer npack1
      real*8 aorbs(npack1,*)
      complex*16 tmp(*)
      real*8     tmp2(*)
      integer ierr

#ifdef TCGMSG
#include "tcgmsg.fh"
#include "msgtypesf.h"
#endif

*    *** local variables ***
      logical pio
      integer MASTER,taskid,taskid_i,taskid_p,com_p
      parameter(MASTER=0)
      integer n,l
      integer msglen
      integer iatom(2)
      character*255 full_filename

      logical  control_parallel_io
      external control_parallel_io


      call Parallel_taskid(taskid)
      call Parallel2d_taskid_i(taskid_i)

      pio = control_parallel_io()
      if (pio) then
         taskid_p = taskid_i
         com_p = 1
      else
         taskid_p = taskid
         com_p = 0
      end if

*     **** open fname binary file ****
      if (taskid_p.eq.MASTER) then
         call util_file_name_noprefix(fname,.false.,
     >                             .false.,
     >                       full_filename)
         l = index(full_filename,' ') - 1
         call openfile(5,full_filename,l,'r',l)
         call iread(5,version,1)
         call iread(5,nfft,3)
         call dread(5,unita,9)
         call cread(5,atom,2)
         call iread(5,lmax,1)
         call iread(5,locp,1)
         call dread(5,rcut,1)
         call dread(5,lmbda,1)
         lmmax=(lmax+1)**2 - (2*locp+1)
      end if

*     **** send header data to all processors ****
      msglen = 1
      call Parallela_Brdcst_ivalues(com_p,MASTER,msglen,version)
      msglen = 3
      call Parallela_Brdcst_ivalues(com_p,MASTER,msglen,nfft)
      msglen = 9
      call Parallela_Brdcst_values(com_p,MASTER,msglen,unita)

      iatom(1) = ichar(atom(1:1))
      iatom(2) = ichar(atom(2:2))
      msglen = 2
      call Parallela_Brdcst_ivalues(com_p,MASTER,msglen,iatom)
      atom(1:1) = char(iatom(1))
      atom(2:2) = char(iatom(2))

      msglen = 1
      call Parallela_Brdcst_ivalues(com_p,MASTER,msglen,lmax)
      call Parallela_Brdcst_ivalues(com_p,MASTER,msglen,locp)
      call Parallela_Brdcst_values(com_p,MASTER,msglen,rcut)
      call Parallela_Brdcst_values(com_p,MASTER,msglen,lmbda)
      lmmax=(lmax+1)**2 - (2*locp+1)


*     **** read in aorb 3d blocks ****
      do n=1,lmmax
         if (pio) then
            call D3dB_t_read_pio(1,5,tmp2,tmp,-1)
         else
            call D3dB_t_read(1,5,tmp2,tmp,-1)
         end if
         call Pack_t_pack(1,tmp2)
         call Pack_t_Copy(1,tmp2,aorbs(1,n))
      end do


*     *** close fname binary file ***
      if (taskid_p.eq.MASTER) then
c       close(11)
         call closefile(5)
      end if

      ierr = 0
      return
      end

*     ***********************************
*     *					*
*     *	 	  aorbs_readall		*
*     *					*
*     ***********************************

      logical function aorbs_readall()
      implicit none
      
#include "bafdecls.fh"

***** AORBS common block ****
#include "aorbs.fh"

*     **** local variables ****
      integer ngp(3),version,nfft3d,npack1
      integer ia,l,lm,ii,icount
      real*8 unita(3,3)
      integer tmp(2),tmp2(2),ierr
      logical value,found,correct_box,value2
      character*2  atom
      character*4  element
      character*50 fname

*     **** parallel i/o variable ****
      integer MASTER,taskid
      parameter(MASTER=0)

*     **** external functions ****
      logical      nwpw_filefind
      integer      control_ngrid
      real*8       control_unita
      character*12 control_boundry
      character*4  ion_atom_qm
      external     nwpw_filefind
      external     control_ngrid
      external     control_unita
      external     control_boundry
      external     ion_atom_qm
      integer      ion_nion_qm,ion_katm_qm
      external     ion_nion_qm,ion_katm_qm
      

      call D3dB_nfft3d(1,nfft3d)
      call Pack_npack(1,npack1)
      call Parallel_taskid(taskid)

      value = BA_push_get(mt_dbl,(2*nfft3d),'tmp',tmp(2),tmp(1))
      if (.not. value) go to 9000

      value = BA_push_get(mt_dbl,(nfft3d),'tmp2',tmp2(2),tmp2(1))
      if (.not. value) go to 9000

*     **** read pseudopotentials ****
      do ia=1,nkatmx

*      **** define formatted aorb name ****
       element = '    '
       element = ion_atom_qm(ia)
       l = index(element,' ') - 1
       fname = element(1:l)//'.aorb'
        

       found = .false.
       do while (.not.found)

         if (nwpw_filefind(fname)) then
            call aorbs_read(fname,
     >                  version,
     >                  ngp,unita,
     >                  atom,
     >                  int_mb(lmmax(1)+ia-1),
     >                  int_mb(lmax(1)+ia-1),
     >                  int_mb(locp(1)+ia-1),
     >                  dbl_mb(rcut(1)+ia-1),
     >                  dbl_mb(lmbda(1)+ia-1),
     >                  npack1,
     >                  dbl_mb(aorbs(1)+ (ia-1)*npack1*norbs_max),
     >                  dbl_mb(tmp(1)),dbl_mb(tmp2(1)),
     >                  ierr)

           if (ierr.gt.0) then
              value = .false.
              go to 9000
           end if

*          **************************************************************
*          ***** logic for finding out if psp is correctly formatted ****
*          **************************************************************
           correct_box = .true.
           if ( (ngp(1).ne.control_ngrid(1)) .or.
     >       (ngp(2).ne.control_ngrid(2)) .or.
     >       (ngp(3).ne.control_ngrid(3)) .or. 
     >       (unita(1,1).ne.control_unita(1,1)) .or.
     >       (unita(2,1).ne.control_unita(2,1)) .or.
     >       (unita(3,1).ne.control_unita(3,1)) .or.
     >       (unita(1,2).ne.control_unita(1,2)) .or.
     >       (unita(2,2).ne.control_unita(2,2)) .or.
     >       (unita(3,2).ne.control_unita(3,2)) .or.
     >       (unita(1,3).ne.control_unita(1,3)) .or.
     >       (unita(2,3).ne.control_unita(2,3)) .or.
     >       (unita(3,3).ne.control_unita(3,3))) then
              correct_box = .false.
              if (taskid.eq.MASTER) then
              write(6,*) "atomic orbitals are not correctly formatted:",
     >                    fname
              end if
           end if
           if (correct_box) found = .true.

         end if

*        **** generate formatted pseudopotential atom.aorb *****
         if (.not.found) then
             call aorbs_formatter_auto(ion_atom_qm(ia),0.0d0,0.0d0)
         end if

       end do !***do while ****


      end do

*     ***********************************************
*     **** set up the index for the atomic basis ****
*     ***********************************************
      icount = 0
      do ii=1,ion_nion_qm()
         ia = ion_katm_qm(ii)

         do lm=1,int_mb(lmmax(1)+ia-1)
            icount = icount + 1

            int_mb(lmaorb(1)+icount-1)  = lm
            int_mb(iaaorb(1)+icount-1)  = ia
            int_mb(iiaorb(1)+icount-1)  = ii
            int_mb(basisaorb(1)+lm-1+(ii-1)*norbs_max) = icount
         end do
      end do
      ibasis = icount
      call aorbs_normalize()

 9000 value2 = BA_pop_stack(tmp2(2))
      value2 = BA_pop_stack(tmp(2))
      

      aorbs_readall = value
      return
      end

