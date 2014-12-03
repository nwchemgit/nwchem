*
* $Id$
*

*     *******************************
*     *                             *
*     *     psi_history_read        *
*     *                             *
*     *******************************

      subroutine psi_history_read(ispin,ne,psi,psitmp)
      implicit none 
      integer ispin,ne(2)
      complex*16 psi(*)
      complex*16 psitmp(*)

#include "bafdecls.fh"

*     **** local variables ****
      integer MASTER,taskid
      parameter(MASTER=0)

      logical value
      character*50 tag1
      integer      ptr,i,size,nfft3d,neall,mm(2),ms,psi2(2)
      real*8       coef(15),scal

*     **** external functions ****
      integer  psi_history,psi_history_size,psi_history_pointer_start
      external psi_history,psi_history_size,psi_history_pointer_start


      call Pack_npack(1,nfft3d)
      neall = ne(1)+ne(2)

*     **** write out psi and ion histories ****
      if (psi_history().gt.0) then
      if (psi_history_size().gt.1) then

        size = psi_history_size()
        call psi_history_generate_coef(coef)

c        call Parallel_taskid(taskid)
c        if (taskid.eq.MASTER) then
c        write(*,*)
c        write(*,*) "Number of psi histories used:",size
c        do i=1,size-1
c          write(*,*) "coef(i):",i,coef(i)
c        end do
c        write(*,*)
c        end if

        value = BA_push_get(mt_dbl,(ne(1)*ne(1)),'mm',mm(2),mm(1))
        value = BA_push_get(mt_dcpl,(neall*nfft3d),
     >                      'psi2',psi2(2),psi2(1))


        ptr = psi_history_pointer_start()
        ptr = ptr-1
        if (ptr.le.0) ptr = size

        call Grsm_gg_Copy(nfft3d,neall,psi,dcpl_mb(psi2(1)))
c       scal = coef(1)
c       call Grsm_gg_daxpy(nfft3d,neall,scal,psi,psi)
        do ms=1,ispin
           call Grsm_ggm_dot(nfft3d,ne(ms),
     >                       dcpl_mb(psi2(1)+(ms-1)*nfft3d*ne(1)),
     >                       psitmp(1+(ms-1)*nfft3d*ne(1)),
     >                       dbl_mb(mm(1)))
           scal = coef(1)
           call dscal(ne(ms)*ne(ms),scal,dbl_mb(mm(1)),1)
           call Grsm_gg_Copy(nfft3d,ne(ms),
     >                       dcpl_mb(psi2(1)+(ms-1)*nfft3d*ne(1)),
     >                       psitmp(1+(ms-1)*nfft3d*ne(1)))
           call Grsm_gmg_daxpy(nfft3d,ne(ms),
     >                       psitmp(1+(ms-1)*nfft3d*ne(1)),
     >                       dbl_mb(mm(1)),
     >                       psi(1+(ms-1)*nfft3d*ne(1)))
        end do

        do i=2,(size-1)
          ptr = ptr-1
          if (ptr.le.0) ptr = size
          tag1 = 'nwpw_history_movecs'//CHAR(ICHAR('a')+ptr-1)
          call psi_read_filename(tag1,ispin,ne,psitmp)

c         scal = coef(i)-coef(i-1)
c         call Grsm_gg_daxpy(nfft3d,neall,scal,psitmp,psi)
          do ms=1,ispin
             call Grsm_ggm_dot(nfft3d,ne(ms),
     >                         dcpl_mb(psi2(1)+(ms-1)*nfft3d*ne(1)),
     >                         psitmp(1+(ms-1)*nfft3d*ne(1)),
     >                         dbl_mb(mm(1)))
             scal = coef(i)-coef(i-1)
             call dscal(ne(ms)*ne(ms),scal,dbl_mb(mm(1)),1)
             call Grsm_gg_Copy(nfft3d,ne(ms),
     >                         dcpl_mb(psi2(1)+(ms-1)*nfft3d*ne(1)),
     >                         psitmp(1+(ms-1)*nfft3d*ne(1)))
             call Grsm_gmg_daxpy(nfft3d,ne(ms),
     >                         psitmp(1+(ms-1)*nfft3d*ne(1)),
     >                         dbl_mb(mm(1)),
     >                         psi(1+(ms-1)*nfft3d*ne(1)))
          end do

        end do

        ptr = ptr-1
        if (ptr.le.0) ptr = size
        tag1 = 'nwpw_history_movecs'//CHAR(ICHAR('a')+ptr-1)
        call psi_read_filename(tag1,ispin,ne,psitmp)
c       scal = -coef(size-1)
c       call Grsm_gg_daxpy(nfft3d,neall,scal,psitmp,psi)
        do ms=1,ispin
           call Grsm_ggm_dot(nfft3d,ne(ms),
     >                       dcpl_mb(psi2(1)+(ms-1)*nfft3d*ne(1)),
     >                       psitmp(1+(ms-1)*nfft3d*ne(1)),
     >                       dbl_mb(mm(1)))
           scal = -coef(size-1)
           call dscal(ne(ms)*ne(ms),scal,dbl_mb(mm(1)),1)
           call Grsm_gg_Copy(nfft3d,ne(ms),
     >                       dcpl_mb(psi2(1)+(ms-1)*nfft3d*ne(1)),
     >                       psitmp(1+(ms-1)*nfft3d*ne(1)))
           call Grsm_gmg_daxpy(nfft3d,ne(ms),
     >                       psitmp(1+(ms-1)*nfft3d*ne(1)),
     >                       dbl_mb(mm(1)),
     >                       psi(1+(ms-1)*nfft3d*ne(1)))
        end do

        value = BA_pop_stack(psi2(2))
        value = BA_pop_stack(mm(2))
      end if
      end if

      return
      end



*     *******************************
*     *                             *
*     *  psi_history_generate_coef  *
*     *                             *
*     *******************************

      subroutine psi_history_generate_coef(coef)
      implicit none 
#include "errquit.fh"
      real*8 coef(*)

#include "bafdecls.fh"

*     **** local variables ****
      integer MASTER,taskid
      parameter(MASTER=0)

      logical      value
      character*50 tag2
      integer      ptr,i,ii,j,size,nion,start
      integer M(2),ipiv(2),info
      integer rion1(2),rion2(2),gion(2)
      real*8 gx,gy,gz
      real*8 g0x,g0y,g0z
      real*8 grx,gry,grz
      real*8 gsx,gsy,gsz

*     **** external functions ****
      integer  psi_history_size,psi_history_pointer_start
      integer  ion_nion
      real*8   ion_rion
       
      external psi_history_size,psi_history_pointer_start
      external ion_nion
      external ion_rion

      call Parallel_taskid(taskid)
      size = psi_history_size()
      ptr  = psi_history_pointer_start()
      start = ptr
      call dcopy((size-1),0.0d0,0,coef,1)

      if (taskid.eq.MASTER) then

      nion = ion_nion()
      value = BA_push_get(mt_dbl,3*nion,'rion1',rion1(2),rion1(1))
      value = value.and.
     >        BA_push_get(mt_dbl,3*nion,'rion2',rion2(2),rion2(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(3*nion*size),'gion',gion(2),gion(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(size*size),'M',M(2),M(1))
      value = value.and.
     >        BA_push_get(mt_int,(size),'ipiv',ipiv(2),ipiv(1))
      if (.not. value) call errquit('psi_history_generate_coef',1,
     &       MA_ERR)

    
      do ii=1,nion
         dbl_mb(rion1(1)+3*(ii-1))   = ion_rion(1,ii)
         dbl_mb(rion1(1)+3*(ii-1)+1) = ion_rion(2,ii)
         dbl_mb(rion1(1)+3*(ii-1)+2) = ion_rion(3,ii)
      end do

      do i=1,size
         ptr = ptr - 1
         if (ptr.le.0) ptr = size
         tag2 = 'nwpw_history_ion'//CHAR(ICHAR('a')+ptr-1)
         call ion_read_filename(tag2,dbl_mb(rion2(1)))

*        **** generate gion ****
         do ii=1,nion 
            gx = dbl_mb(rion1(1)+3*(ii-1)) 
     >         - dbl_mb(rion2(1)+3*(ii-1))
            gy = dbl_mb(rion1(1)+3*(ii-1)+1) 
     >         - dbl_mb(rion2(1)+3*(ii-1)+1)
            gz = dbl_mb(rion1(1)+3*(ii-1)+2) 
     >         - dbl_mb(rion2(1)+3*(ii-1)+2)
            dbl_mb(gion(1)+3*nion*(i-1)+ 3*(ii-1))   = gx
            dbl_mb(gion(1)+3*nion*(i-1)+ 3*(ii-1)+1) = gy
            dbl_mb(gion(1)+3*nion*(i-1)+ 3*(ii-1)+2) = gz
         end do

         call dcopy(3*nion,dbl_mb(rion2(1)),1,dbl_mb(rion1(1)),1)

      end do

*     **** generate b ****
      do i=1,size-1
         coef(i) = 0.0d0
         do ii=1,nion
            g0x = dbl_mb(gion(1)+ 3*(ii-1))   
            g0y = dbl_mb(gion(1)+ 3*(ii-1)+1) 
            g0z = dbl_mb(gion(1)+ 3*(ii-1)+2) 
            grx = dbl_mb(gion(1)+3*nion*(i)+ 3*(ii-1))   
            gry = dbl_mb(gion(1)+3*nion*(i)+ 3*(ii-1)+1) 
            grz = dbl_mb(gion(1)+3*nion*(i)+ 3*(ii-1)+2) 
            coef(i) = coef(i) + g0x*grx + g0y*gry + g0z*grz
         end do
      end do

*     **** generate M ****
      do j=1,size-1
      do i=1,size-1
         dbl_mb(M(1)+(i-1)+(j-1)*(size-1))  =0.0d0
         do ii=1,nion
            grx = dbl_mb(gion(1)+3*nion*(i)+ 3*(ii-1))   
            gry = dbl_mb(gion(1)+3*nion*(i)+ 3*(ii-1)+1) 
            grz = dbl_mb(gion(1)+3*nion*(i)+ 3*(ii-1)+2) 
            gsx = dbl_mb(gion(1)+3*nion*(j)+ 3*(ii-1))   
            gsy = dbl_mb(gion(1)+3*nion*(j)+ 3*(ii-1)+1) 
            gsz = dbl_mb(gion(1)+3*nion*(j)+ 3*(ii-1)+2) 

            dbl_mb(M(1)+(i-1)+(j-1)*(size-1)) 
     >         = dbl_mb(M(1)+(i-1)+(j-1)*(size-1)) 
     >         + grx*gsx + gry*gsy + grz*gsz
         end do
      end do
      end do

*     **** solve Mx=b ****
      call dgesv((size-1),1,dbl_mb(M(1)),(size-1),
     >           int_mb(ipiv(1)),coef,(size-1),info)
      if (info.ne.0) then
        call dcopy((size-1),0.0d0,0,coef,1)
      end if
      do i=1,size-1
        if (dabs(coef(i)) .gt.100.0d0) 
     >    call dcopy((size-1),0.0d0,0,coef,1)
      end do


*     **** write out estimated ion positions ****
c      ptr = start - 1
c      if (ptr.le.0) ptr = size
c      tag2 = 'nwpw_history_ion'//CHAR(ICHAR('a')+ptr-1)
c      call ion_read_filename(tag2,dbl_mb(rion1(1)))

c      gx = 0.0d0
c      write(*,*) 'old and new ion positions'
c      do ii=1,nion
c        write(*,1500) ii, 
c     >                dbl_mb(rion1(1)+3*(ii-1)),
c     >                dbl_mb(rion1(1)+3*(ii-1)+1),
c     >                dbl_mb(rion1(1)+3*(ii-1)+2),
c     >                ion_rion(1,ii),
c     >                ion_rion(2,ii),
c     >                ion_rion(3,ii) 
c       gx = gx
c     >    + (dbl_mb(rion1(1)+3*(ii-1))  -ion_rion(1,ii))**2
c     >    + (dbl_mb(rion1(1)+3*(ii-1)+1)-ion_rion(2,ii))**2
c     >    + (dbl_mb(rion1(1)+3*(ii-1)+2)-ion_rion(3,ii))**2
c      end do
c      write(*,*) "diff:",dsqrt(gx)

c      do i=1,(size-1)
c      do ii=1,nion
c        dbl_mb(rion1(1)+3*(ii-1))   = dbl_mb(rion1(1)+3*(ii-1))
c     >   + coef(i)*dbl_mb(gion(1) + 3*nion*(i) + 3*(ii-1))
c        dbl_mb(rion1(1)+3*(ii-1)+1) = dbl_mb(rion1(1)+3*(ii-1)+1)
c     >   + coef(i)*dbl_mb(gion(1) + 3*nion*(i) + 3*(ii-1)+1)
c        dbl_mb(rion1(1)+3*(ii-1)+2) = dbl_mb(rion1(1)+3*(ii-1)+2)
c     >   + coef(i)*dbl_mb(gion(1) + 3*nion*(i) + 3*(ii-1)+2)
c      end do
c      end do

c      gx = 0.0d0
c      write(*,*) 'estimated and new ion positions'
c      do ii=1,nion
c        write(*,1500) ii, 
c     >                dbl_mb(rion1(1)+3*(ii-1)),
c     >                dbl_mb(rion1(1)+3*(ii-1)+1),
c     >                dbl_mb(rion1(1)+3*(ii-1)+2),
c     >                ion_rion(1,ii),
c     >                ion_rion(2,ii),
c     >                ion_rion(3,ii) 
c       gx = gx
c     >    + (dbl_mb(rion1(1)+3*(ii-1))  -ion_rion(1,ii))**2
c     >    + (dbl_mb(rion1(1)+3*(ii-1)+1)-ion_rion(2,ii))**2
c     >    + (dbl_mb(rion1(1)+3*(ii-1)+2)-ion_rion(3,ii))**2
c      end do
c      write(*,*) "diff:",dsqrt(gx)

 1500 format(5x,I3, '(', 3F8.3, ') (' ,3F8.3, ')' )

      value = BA_pop_stack(ipiv(2))
      value = value.and.BA_pop_stack(M(2))
      value = value.and.BA_pop_stack(gion(2))
      value = value.and.BA_pop_stack(rion2(2))
      value = value.and.BA_pop_stack(rion1(2))
      if (.not. value) call errquit('psi_history_generate_coef',2,
     &       MA_ERR)

      end if

      call D3dB_Vector_SumAll(size-1,coef)
      
      return
      end

*     *******************************
*     *             				*
*     *     psi_history_write 		*
*     *             				*
*     *******************************

      subroutine psi_history_write(ispin,ne,psi)
      implicit none 
      integer ispin,ne(2)
      complex*16 psi(*)


*     **** local variables ****
      character*50 tag1,tag2
      integer      pointer

*     **** external functions ****
      logical  psi_history_check_last
      integer  psi_history,psi_history_size,psi_history_pointer_start
      external psi_history_check_last
      external psi_history,psi_history_size,psi_history_pointer_start
    

*     **** write out psi and ion histories ****
      if (psi_history().gt.0) then
      if (psi_history_check_last()) then
        call psi_history_add()
        pointer = psi_history_pointer_start() - 1
        if (pointer.le.0) pointer = psi_history_size()
        tag1 = 'nwpw_history_movecs'//CHAR(ICHAR('a')+pointer-1)
        tag2 = 'nwpw_history_ion'//CHAR(ICHAR('a')+pointer-1)
        call psi_write_filename(tag1,ispin,ne,psi) 
        call ion_write_filename(tag2)

      end if
      end if
      return
      end

*     *******************************
*     *                             *
*     *   psi_history_check_last    *
*     *                             *
*     *******************************

      logical function psi_history_check_last()
      implicit none 
#include "errquit.fh"

#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "bafdecls.fh"
#include "btdb.fh"

*     **** local variables ****
      integer MASTER
      parameter(MASTER=0)

      logical check,value
      integer ptr,nion,rion(2),ii,size
      integer icheck,msglen
      real*8  diff
      character*50 tag2

*     **** external functions ****
      integer  psi_history_pointer_start,psi_history_size
      integer  ion_nion
      real*8   ion_rion
      external psi_history_pointer_start,psi_history_size
      external ion_nion
      external ion_rion
     
      size = psi_history_size()
      if (size.le.0) then
        psi_history_check_last = .true.
        return
      end if

      nion = ion_nion()
      value = BA_push_get(mt_dbl,3*nion,'rion',rion(2),rion(1))
      if (.not. value) call errquit('psi_history_check_last',1, MA_ERR)

      ptr = psi_history_pointer_start() - 1
      if (ptr.le.0) ptr = size

      tag2 = 'nwpw_history_ion'//CHAR(ICHAR('a')+ptr-1)
      call ion_read_filename(tag2,dbl_mb(rion(1)))

      check = .false.
      do ii=1,nion
         diff = dbl_mb(rion(1)+3*(ii-1))-ion_rion(1,ii)
         if (dabs(diff).gt.1.0d-5) check=.true.

         diff = dbl_mb(rion(1)+3*(ii-1)+1)-ion_rion(2,ii)
         if (dabs(diff).gt.1.0d-5) check=.true.

         diff = dbl_mb(rion(1)+3*(ii-1)+2)-ion_rion(3,ii)
         if (dabs(diff).gt.1.0d-5) check=.true.
      end do

      value = BA_pop_stack(rion(2))
      if (.not. value) call errquit('psi_history_check_last',2, MA_ERR)

      icheck = 0
      if (check) icheck=1
      msglen = 1
      call BRDCST(9+MSGINT,icheck,mitob(msglen),MASTER)
      check = .false.
      if (icheck.eq.1) check=.true.
     

      psi_history_check_last = check
      return
      end

*     *******************************
*     *             				*
*     *           psi_history 		*
*     *             				*
*     *******************************

      integer function psi_history()
      implicit none 

#include "bafdecls.fh"
#include "btdb.fh"

*     **** local variables ****
      integer rtdb,history

*     **** external functions ****
      integer  control_rtdb
      external control_rtdb
     
      rtdb = control_rtdb()

      if (.not.btdb_get(rtdb,'nwpw:psi_history',mt_int,1,history)) then
         history = 0
      end if
      if (history.gt.0) history = history + 1

      psi_history = history
      return
      end


*     ***********************************
*     *             					*
*     *           psi_history_size 		*
*     *             					*
*     ***********************************

      integer function psi_history_size()
      implicit none 

#include "bafdecls.fh"
#include "btdb.fh"

*     **** local variables ****
      integer rtdb,size

*     **** external functions ****
      integer  control_rtdb
      external control_rtdb
     
      rtdb = control_rtdb()

      if (.not.btdb_get(rtdb,'nwpw:psi_history_size',mt_int,1,size)) 
     >   then
         size = 0
      end if

      psi_history_size = size
      return
      end

*     ***********************************
*     *             					*
*     *    psi_history_pointer_start 	*
*     *             					*
*     ***********************************

      integer function psi_history_pointer_start()
      implicit none 

#include "bafdecls.fh"
#include "btdb.fh"

*     **** local variables ****
      integer rtdb,start

*     **** external functions ****
      integer  control_rtdb
      external control_rtdb
     
      rtdb = control_rtdb()

      if (.not.btdb_get(rtdb,'nwpw:psi_history_pointer_start',
     >    mt_int,1,start)) then
         start = 1
      end if

      psi_history_pointer_start = start
      return
      end



*     ***********************************
*     *             					*
*     *        psi_history_shift	 	*
*     *             					*
*     ***********************************

      subroutine psi_history_shift()
      implicit none 
#include "errquit.fh"

#include "bafdecls.fh"
#include "btdb.fh"

*     **** local variables ****
      integer rtdb,start,size

*     **** external functions ****
      integer  control_rtdb,psi_history_size,psi_history_pointer_start
      external control_rtdb,psi_history_size,psi_history_pointer_start
     
      rtdb = control_rtdb()
      size  = psi_history_size()
      start = psi_history_pointer_start()


      start = start + 1
      if (start.gt.size) start = 1

      if (.not.btdb_put(rtdb,'nwpw:psi_history_pointer_start',
     >    mt_int,1,start)) then
        call errquit('psi_history_shift',0, RTDB_ERR)
      end if
      return
      end

*     ***********************************
*     *             					*
*     *        psi_history_adder	 	*
*     *             					*
*     ***********************************

      subroutine psi_history_adder()
      implicit none 
#include "errquit.fh"

#include "bafdecls.fh"
#include "btdb.fh"

*     **** local variables ****
      integer rtdb,size

*     **** external functions ****
      integer  control_rtdb,psi_history_size
      external control_rtdb,psi_history_size
     
      size  = psi_history_size()
      size  = size + 1

      rtdb = control_rtdb()
      if (.not.btdb_put(rtdb,'nwpw:psi_history_size',
     >    mt_int,1,size)) then
        call errquit('psi_history_adder',0, RTDB_ERR)
      end if

      return
      end

*     ***********************************
*     *             					*
*     *        psi_history_add		 	*
*     *             					*
*     ***********************************

      subroutine psi_history_add()
      implicit none 

*     **** local variables ****
      integer size,history

*     **** external functions ****
      integer  psi_history,psi_history_size
      external psi_history,psi_history_size
     
      size    = psi_history_size()
      history = psi_history()

      if (size.ge.history) then
        call psi_history_shift()
      else
        call psi_history_adder()
      end if

      return
      end


*     ***********************************
*     *                                 *
*     *           ion_write_filename    *
*     *                                 *
*     ***********************************

      subroutine ion_write_filename(filename)
      implicit none 
      character*(*) filename

#include "bafdecls.fh"
#include "errquit.fh"


*    *** local variables ***
      integer MASTER,taskid
      parameter(MASTER=0)

      logical value
      integer l,i
      integer nion,rion(2)
      character*255 full_filename
      double precision rx,ry,rz


*     ***** local functions ****
      double precision ion_rion
      external         ion_rion
      integer  ion_nion
      external ion_nion

      call nwpw_timing_start(50)
      call Parallel_taskid(taskid)

      nion = ion_nion()
      value = BA_push_get(mt_dbl,3*nion,'rion',rion(2),rion(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

      do i=1,nion
         rx = ion_rion(1,i)
         ry = ion_rion(2,i)
         rz = ion_rion(3,i)
         dbl_mb(rion(1)+3*(i-1))   = rx
         dbl_mb(rion(1)+3*(i-1)+1) = ry
         dbl_mb(rion(1)+3*(i-1)+2) = rz
      end do

*     **** open IONOUT binary file ****
      if (taskid.eq.MASTER) then
         call util_file_name_noprefix(filename,.false.,
     >                                .false.,
     >                       full_filename)
         l = index(full_filename,' ') -1
         call openfile(6,full_filename,l,'w',l)
         call iwrite(6,nion,1)
         call dwrite(6,dbl_mb(rion(1)),3*nion)
      end if

*     *** close ELCIN binary file ***
      if (taskid.eq.MASTER) then
        call closefile(6)
      end if

      value = BA_pop_stack(rion(2))
      call nwpw_timing_end(50)
      return
      end


*     ***********************************
*     *                                 *
*     *           ion_read_filename     *
*     *                                 *
*     ***********************************

      subroutine ion_read_filename(filename,rion)
      implicit none 
      character*(*) filename
      real*8 rion(*)

#include "bafdecls.fh"


*    *** local variables ***
      integer MASTER,taskid
      parameter(MASTER=0)

      integer l
      integer nion
      character*255 full_filename


*     ***** local functions ****
      integer  ion_nion
      external ion_nion

      call nwpw_timing_start(50)
      call Parallel_taskid(taskid)


*     **** open IONOUT binary file ****
      if (taskid.eq.MASTER) then
         call util_file_name_noprefix(filename,.false.,
     >                                .false.,
     >                       full_filename)
         l = index(full_filename,' ') -1
         call openfile(5,full_filename,l,'r',l)
         call iread(5,nion,1)
         call dread(5,rion,3*nion)
      end if

*     *** close ELCIN binary file ***
      if (taskid.eq.MASTER) then
        call closefile(5)
      end if

      call nwpw_timing_end(50)
      return
      end


