
*
* $Id$
*

*     ***********************************
*     *					*
*     *	 	     borbs_init		*
*     *					*
*     ***********************************
      logical function borbs_init()      
      implicit none

#include "bafdecls.fh"
#include "borbs.fh"

*     **** local variables ****
      integer npack1,nion
      logical value

*     **** external functions *****
      integer  ion_nkatm_qm,ion_nion_qm
      external ion_nkatm_qm,ion_nion_qm
 

      call Pack_npack(1,npack1)
      nion   = ion_nion_qm()
      nkatmx = ion_nkatm_qm()

*     **** allocate borb datastructure  and borb tag lists****
      call borb_projector_init(2*nkatmx)
      value = BA_alloc_get(mt_int,(nkatmx),'borbs',borbs(2),borbs(1))

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
     >                    'lmborb',lmborb(2),lmborb(1))
      value = value.and.
     >        BA_alloc_get(mt_int,(nion*norbs_max),
     >                    'iaborb',iaborb(2),iaborb(1))
      value = value.and.
     >        BA_alloc_get(mt_int,(nion*norbs_max),
     >                    'iiborb',iiborb(2),iiborb(1))
      value = value.and.
     >        BA_alloc_get(mt_int,(nion*norbs_max),
     >                    'basisborb',basisborb(2),basisborb(1))

      borbs_init = value
      return
      end


*     ***********************************
*     *					*
*     *	 	     borbs_end		*
*     *					*
*     ***********************************
      subroutine borbs_end()      
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"
#include "borbs.fh"

*     **** local variables ****
      logical value

*     **** deallocate projector data ****
      call borb_projector_end()

      value = BA_free_heap(borbs(2))
      value = value.and.BA_free_heap(lmmax(2))
      value = value.and.BA_free_heap(lmax(2))
      value = value.and.BA_free_heap(locp(2))
      value = value.and.BA_free_heap(rcut(2))
      value = value.and.BA_free_heap(lmbda(2))
      value = value.and.BA_free_heap(lmborb(2))
      value = value.and.BA_free_heap(iaborb(2))
      value = value.and.BA_free_heap(iiborb(2))
      value = value.and.BA_free_heap(basisborb(2))
      if (.not. value) call errquit('borbs_end:freeing heap memory',0,
     &       MA_ERR)

      return
      end


*     ***********************************
*     *					*
*     *	 	     borbs_norbs	*
*     *					*
*     ***********************************

      integer function borbs_norbs(ia)
      implicit none
      integer  ia

#include "bafdecls.fh"
#include "borbs.fh"

      borbs_norbs = int_mb(lmmax(1)+ia-1)
      return
      end


*     ***********************************
*     *					*
*     *	 	   borbs_nbasis		*
*     *					*
*     ***********************************
      integer function borbs_nbasis()
      implicit none

#include "borbs.fh"

      borbs_nbasis = ibasis
      return
      end

*     ***********************************
*     *                                 *
*     *              borbs_lmax         *
*     *                                 *
*     ***********************************

      integer function borbs_lmax(ia)
      implicit none
      integer  ia

#include "bafdecls.fh"
#include "borbs.fh"

      borbs_lmax = int_mb(lmax(1)+ia-1)-1
      return
      end


*     ***********************************
*     *                                 *
*     *              borbs_rcut         *
*     *                                 *
*     ***********************************
      real*8 function borbs_rcut(ia)
      implicit none
      integer  ia

#include "bafdecls.fh"
#include "borbs.fh"

      borbs_rcut = dbl_mb(rcut(1)+ia-1)
      return
      end

*     ***********************************
*     *                                 *
*     *           borbs_lmbda           *
*     *                                 *
*     ***********************************
      real*8 function borbs_lmbda(ia)
      implicit none
      integer  ia

#include "bafdecls.fh"
#include "borbs.fh"

      borbs_lmbda = dbl_mb(lmbda(1)+ia-1)
      return
      end


*     ***********************************
*     *                                 *
*     *              borbs_l            *
*     *                                 *
*     ***********************************

      integer function borbs_l(ia,n)
      implicit none
      integer  ia
      integer n          ! basis number

#include "bafdecls.fh"
#include "borbs.fh"

*     *** local variables ***
      integer l,m,lm
      
      lm = int_mb(lmborb(1)+n-1)
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

      borbs_l = l
      return
      end


*     ***********************************
*     *                                 *
*     *     borbs_get_basis_number      *
*     *                                 *
*     ***********************************

      integer function borbs_get_basis_number(ii,lm)
      implicit none
      integer ii,lm

#include "bafdecls.fh"
#include "borbs.fh"

      borbs_get_basis_number=int_mb(basisborb(1)+lm-1+(ii-1)*norbs_max)

      return
      end



*     ***********************************
*     *					*
*     *	 	   borbs_normalize  	*
*     *					*
*     ***********************************
      subroutine borbs_normalize()
      implicit none

#include "bafdecls.fh"
#include "borbs.fh"

*     **** local variables ****
      integer lm,ia,npack1,nbrillq,nbq,shift
      real*8  sum

*     **** external functions ****
      integer  Pneb_nbrillq,borb_projector_get_ptr
      external Pneb_nbrillq,borb_projector_get_ptr


      nbrillq = Pneb_nbrillq()

*     **** Normalize atomic orbitals in k space *****
      do nbq=1,nbrillq
         do ia=1,nkatmx
         do lm=1,int_mb(lmmax(1)+ia-1)
         
            shift = borb_projector_get_ptr(int_mb(borbs(1)+ia-1),nbq,lm)
            call Cram_rr_dot(nbq,dbl_mb(shift),dbl_mb(shift),sum)
            sum = 1.0d0/dsqrt(sum)
            call Cram_r_SMul1(nbq,sum,dbl_mb(shift))
         end do
         end do
      end do

      return
      end

*     ***********************************
*     *					*
*     *	 	   borbs_weight		*
*     *					*
*     ***********************************
      real*8 function borbs_weight(n)
      implicit none
      integer n          ! basis number

#include "bafdecls.fh"
#include "borbs.fh"

*     **** local variables ****
      integer ia
      real*8 zv,zcount

*     **** external functions ****
      real*8   psp_zv
      external psp_zv

      ia     = int_mb(iaborb(1)+n-1)
      zcount = int_mb(lmmax(1)+ia-1)
      zv     = psp_zv(ia)

      borbs_weight = (zv/zcount)
      return
      end

*     ***********************************
*     *					*
*     *	 	   borbs_borb  		*
*     *					*
*     ***********************************

      subroutine borbs_borb(nbq,n,borb)
      implicit none
      integer nbq
      integer n          ! basis number
      complex*16 borb(*) ! return orbital

#include "bafdecls.fh"
#include "errquit.fh"
#include "borbs.fh"

*     **** local variables ****
      logical value
      integer lm,ia,ii
      integer nfft3d,npack1,shift
      integer exi(2)
      complex*16 cxr

*     **** external functions ****
      logical  is_sORd
      external is_sORd
      integer  borb_projector_get_ptr
      external borb_projector_get_ptr

      call C3dB_nfft3d(1,nfft3d)
      call Cram_npack(nbq,npack1)

      value = BA_push_get(mt_dcpl,npack1,'exi', exi(2), exi(1))
      if (.not. value) call errquit('borbs_borb:out of heap memory',0,
     &       MA_ERR)

*     **** structure factor ****
      lm = int_mb(lmborb(1)+n-1)
      ia = int_mb(iaborb(1)+n-1)
      ii = int_mb(iiborb(1)+n-1)
      call cstrfac_pack(nbq,ii,dcpl_mb(exi(1)))
      call cstrfac_k(ii,nbq,cxr)
      call zscal(npack1,cxr,dcpl_mb(exi(1)),1)

      shift = borb_projector_get_ptr(int_mb(borbs(1)+ia-1),nbq,lm)


*     **** phase factor does not matter therefore ****
*     **** (-i)^l is the same as (i)^l in the     ****
*     **** Rayleigh scattering formula            ****

*     *** current function is s or d ****
      if (is_sORd(lm,int_mb(lmax(1)+ia-1),
     >                   int_mb(locp(1)+ia-1))
     >        ) then
         call Cram_rc_Mul(nbq,dbl_mb(shift),dcpl_mb(exi(1)),borb)

*     *** current function is p or f ****
      else
         call Cram_irc_Mul(nbq,dbl_mb(shift),dcpl_mb(exi(1)),borb)
      end if


      value = BA_pop_stack(exi(2))
      if (.not. value) call errquit('borbs_borb:freeing heap memory',0,
     &       MA_ERR)

      return 
      end


*     ***********************************
*     *					*
*     *	 	   borbs_read 		*
*     *					*
*     ***********************************
      subroutine borbs_read(fname,
     >                      version,
     >                       nfft,unita,
     >                       atom,
     >                       lmmax,lmax,locp,rcut,lmbda,
     >                       npack1,borbs_tag,
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
      integer borbs_tag
      complex*16 tmp(*)
      real*8     tmp2(*)
      integer ierr

#include "bafdecls.fh"
#include "btdb.fh"
#include "util.fh"

*    *** local variables ***
      logical mprint,value
      integer MASTER,taskid,taskid_k
      parameter(MASTER=0)
      integer n,l,nbrillioun,nb,nbq,pk
      integer msglen
      integer iatom(2)
      character*255 full_filename
      real*8 kv(3)
      

*     ***** external functions ****
      integer  brillioun_nbrillioun,borb_projector_alloc
      external brillioun_nbrillioun,borb_projector_alloc
      integer  brillioun_nbrillq
      external brillioun_nbrillq
      real*8   brillioun_all_k
      external brillioun_all_k
      logical  control_print
      external control_print

      call Parallel_taskid(taskid)
      call Parallel3d_taskid_k(taskid_k)
      mprint = (taskid.eq.MASTER).and.control_print(print_medium)


*     **** open fname binary file ****
      if (taskid.eq.MASTER) then
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

         call iread(5,nbrillioun,1)
         ierr = 0
         if (nbrillioun.eq.brillioun_nbrillioun()) then
           do nb=1,nbrillioun
               call dread(5,kv,3)
               if ((brillioun_all_k(1,nb).ne.kv(1)).or.
     >             (brillioun_all_k(2,nb).ne.kv(2)).or.
     >             (brillioun_all_k(3,nb).ne.kv(3)))
     >          ierr = 1
           end do
           if (ierr.eq.1) then
              if (mprint) then
                 write(*,*)"Brillioun Zone Vectors do not match!"
                 call flush(6)
              end if
           end if
         else
           if (mprint) then
           write(*,*)"Brillioun Zone Points do not match!"
           write(*,*)"NB = ",nbrillioun," not equal ",
     >       brillioun_nbrillioun()
           call flush(6)
           end if
           ierr = 1
         end if
      end if

      msglen = 1
      call Parallel_Brdcst_ivalues(MASTER,msglen,ierr)
      if (ierr.ne.0) then
         if (taskid.eq.MASTER) call closefile(5)
         return
      end if


*     **** send header data to all processors ****
      msglen = 1
      call Parallel_Brdcst_ivalues(MASTER,msglen,version)
      msglen = 3
      call Parallel_Brdcst_ivalues(MASTER,msglen,nfft)
      msglen = 9
      call Parallel_Brdcst_values(MASTER,msglen,unita)

      iatom(1) = ichar(atom(1:1))
      iatom(2) = ichar(atom(2:2))
      msglen = 2
      call Parallel_Brdcst_ivalues(MASTER,msglen,iatom)
      atom(1:1) = char(iatom(1))
      atom(2:2) = char(iatom(2))

      msglen = 1
      call Parallel_Brdcst_ivalues(MASTER,msglen,lmax)
      call Parallel_Brdcst_ivalues(MASTER,msglen,locp)
      call Parallel_Brdcst_ivalues(MASTER,msglen,nbrillioun)
      call Parallel_Brdcst_values(MASTER,msglen,rcut)
      call Parallel_Brdcst_values(MASTER,msglen,lmbda)
      lmmax=(lmax+1)**2 - (2*locp+1)



*     **** read in borb 3d blocks ****
      borbs_tag = borb_projector_alloc(brillioun_nbrillq(),lmmax,npack1)
      do nb=1,brillioun_nbrillioun()
         call K1dB_ktoqp(nb,nbq,pk)
         do n=1,lmmax
            call C3dB_r_read(1,5,tmp2,tmp,-1,pk,.true.)
            if (pk.eq.taskid_k) then
               call Cram_r_pack(nbq,tmp2)
               call borb_projector_add(borbs_tag,nbq,n,tmp2)
            end if
         end do
      end do

*     *** close fname binary file ***
      if (taskid.eq.MASTER) then
         call closefile(5)
      end if

      ierr = 0
      return
      end

*     ***********************************
*     *					*
*     *	 	  borbs_readall		*
*     *					*
*     ***********************************
      logical function borbs_readall()
      implicit none
      
#include "bafdecls.fh"
#include "borbs.fh"

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
      
      call C3dB_nfft3d(1,nfft3d)
      call Cram_max_npack(npack1)
      call Parallel_taskid(taskid)

      value = BA_push_get(mt_dbl,(4*nfft3d),'tmp',tmp(2),tmp(1))
      if (.not. value) go to 9000

      value = BA_push_get(mt_dbl,(2*nfft3d),'tmp2',tmp2(2),tmp2(1))
      if (.not. value) go to 9000

*     **** read pseudopotentials ****
      do ia=1,nkatmx

*      **** define formatted borb name ****
       element = '    '
       element = ion_atom_qm(ia)
       l = index(element,' ') - 1
       fname = element(1:l)//'.borb'
        

       found = .false.
       do while (.not.found)

         if (nwpw_filefind(fname)) then
            call borbs_read(fname,
     >                  version,
     >                  ngp,unita,
     >                  atom,
     >                  int_mb(lmmax(1)+ia-1),
     >                  int_mb(lmax(1)+ia-1),
     >                  int_mb(locp(1)+ia-1),
     >                  dbl_mb(rcut(1)+ia-1),
     >                  dbl_mb(lmbda(1)+ia-1),
     >                  npack1,
     >                  int_mb(borbs(1)+ ia-1),
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

*        **** generate formatted pseudopotential atom.borb *****
         if (.not.found) then
             call borbs_formatter_auto(ion_atom_qm(ia),0.0d0,0.0d0)
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

            int_mb(lmborb(1)+icount-1)  = lm
            int_mb(iaborb(1)+icount-1)  = ia
            int_mb(iiborb(1)+icount-1)  = ii
            int_mb(basisborb(1)+lm-1+(ii-1)*norbs_max) = icount
         end do
      end do
      ibasis = icount
      call borbs_normalize()

 9000 value2 = BA_pop_stack(tmp2(2))
      value2 = BA_pop_stack(tmp(2))
      

      borbs_readall = value
      return
      end

