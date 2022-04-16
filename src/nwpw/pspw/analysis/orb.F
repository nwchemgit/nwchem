c
c $Id$
c

*     ***********************************
*     *                 		*
*     *           Orb_Init          	*
*     *                 		*
*     ***********************************

      subroutine Orb_Init(rtdb,ispin_in,ne_in)
      implicit none
      integer rtdb
      integer ispin_in
      integer ne_in(2)

#include "bafdecls.fh"
#include "msgtypesf.h"
#include "btdb.fh"
#include "errquit.fh"


*     **** orb common block ****
      integer lorb(2),b0(2),a(2),sum(2)
      integer ne(2),ispin,npack1,nemax
      common / ORB_BLOCK / sum,a,b0,lorb,ne,ispin,npack1,nemax

*     **** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)

      logical value,found,mulliken_kawai
      integer nion,npsp,k
      character*50 filename
      character*255 id,full_filename

*     **** external functions ****
      character*4 ion_atom_qm
      external    ion_atom_qm
      integer  ion_nion_qm,ion_nkatm_qm
      external ion_nion_qm,ion_nkatm_qm
      logical  control_mulliken_kawai,aorbs_init,aorbs_readall
      external control_mulliken_kawai,aorbs_init,aorbs_readall

      call Parallel_taskid(taskid)
      ispin = ispin_in
      ne(1) = ne_in(1)
      ne(2) = ne_in(2)
      nemax = ne(1)+ne(2)
      nion = ion_nion_qm()
      npsp = ion_nkatm_qm()
      call Pack_npack(1,npack1)
      mulliken_kawai = control_mulliken_kawai()

      value = BA_alloc_get(mt_int,npsp,'lorb',lorb(2),lorb(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,6*npsp,'b0',b0(2),b0(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,36*nemax*nion,'a',a(2),a(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,nemax,'sum',sum(2),sum(1))


*     ****************************************
*     **** read in expansion coefficients ****
*     ****************************************
      do k=1,npsp
          id = 'analysis:lorb'//ion_atom_qm(k)
          if (.not. btdb_get(rtdb,id,mt_int,1,int_mb(lorb(1)+k-1)))
     >        call errquit(
     >        'analysis: btdb_get lorb failed', 0, RTDB_ERR)

          id = 'analysis:expansion'//ion_atom_qm(k)
          if (.not. btdb_get(rtdb,id,mt_dbl,(int_mb(lorb(1)+k-1)+1),
     >                                       dbl_mb(b0(1)+(k-1)*6)))
     >        call errquit(
     >        'analysis_input: btdb_put failed', 0, RTDB_ERR)
      end do


      if (.not.mulliken_kawai) then
         value = aorbs_init()
         value = value.and.aorbs_readall()
         if (.not.value) 
     >    call errquit('analysis_input: aorbs_init failed',0,0)
      end if



      if (.not.btdb_cget(rtdb,'cpmd:omotion_filename',1,filename))
     > call util_file_prefix('omotion',filename)

      call util_file_name_noprefix(filename,.false.,.false.,
     >                    full_filename)
      if (taskid.eq.MASTER) then
        inquire(file=full_filename,exist=found)
        if (found) then
          open(unit=33, file=full_filename,form='formatted',
     >         status='old')
           do while (found)
             read(18,*,end=100) 
           end do
  100      continue
#if defined(FUJITSU)||defined(GCC46)
           backspace 18
#endif
        else
          open(unit=33,file=full_filename,form='formatted',
     >         status='new')
        end if
      end if

      return
      end

*     ***********************************
*     *                 		*
*     *           Orb_End 	       	*
*     *                 		*
*     ***********************************

      subroutine Orb_End()
      implicit none

#include "bafdecls.fh"
#include "msgtypesf.h"

*     **** orb common block ****
      integer lorb(2),b0(2),a(2),sum(2)
      integer ne(2),ispin,npack1,nemax
      common / ORB_BLOCK / sum,a,b0,lorb,ne,ispin,npack1,nemax

      integer taskid,MASTER
      parameter (MASTER=0)

      logical value,mulliken_kawai
      logical  control_mulliken_kawai
      external control_mulliken_kawai


      call Parallel_taskid(taskid)

      if (taskid.eq.MASTER) close(unit=33)

      mulliken_kawai = control_mulliken_kawai()
      if (.not.mulliken_kawai) call aorbs_end()

*     **** free heap space ****
      value = BA_free_heap(lorb(2))
      value = BA_free_heap(b0(2))
      value = BA_free_heap(a(2))
      value = BA_free_heap(sum(2))

      return
      end

*     ***********************************
*     *                 		*
*     *           Orb_Write	       	*
*     *                 		*
*     ***********************************

      subroutine Orb_Write(psi)
      implicit none
      complex*16 psi(*)

#include "bafdecls.fh"
#include "msgtypesf.h"

*     **** orb common block ****
      integer lorb(2),b0(2),a(2),sum(2)
      integer ne(2),ispin,npack1,nemax
      common / ORB_BLOCK / sum,a,b0,lorb,ne,ispin,npack1,nemax

      call Orb_Analysis(33,0,ispin,ne,npack1,nemax,psi,
     >                        int_mb(lorb(1)),
     >                        dbl_mb(b0(1)),
     >                        dbl_mb(a(1)),
     >                        dbl_mb(sum(1)))

      return
      end


