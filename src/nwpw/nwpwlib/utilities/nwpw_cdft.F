
*     *******************************************
*     *                                         *
*     *             nwpw_cdft_init              *
*     *                                         *
*     *******************************************
      subroutine nwpw_cdft_init(rtdb0)
      implicit none
      integer rtdb0

#include "bafdecls.fh"
#include "util.fh"
#include "stdio.fh"
#include "btdb.fh"
#include "errquit.fh"

*     **** local variables ****
      integer  MASTER
      parameter(MASTER=0)
      logical  omaster,oprint,value
      integer  taskid
      integer  nion,ii,jj,iat
      real*8   rr


*     **** nwpw_cdft common blocks ****
      logical cdft_on
      integer ucdft(2),qcdft(2),rtdb
      real*8  dielec
      common /nwpw_cdft_blk/ ucdft,qcdft,dielec,rtdb,cdft_on


*     ***** external functions *****
      logical  control_print
      external control_print
      integer  Parallel_threadid
      external Parallel_threadid
      integer  ion_nion,ion_rion_ptr,ion_katm
      external ion_nion,ion_rion_ptr,ion_katm
      real*8   ion_q
      external ion_q
      character*4 ion_atom
      external    ion_atom

c
      call Parallel_taskid(taskid)
      omaster =  (taskid.eq.MASTER).and.(Parallel_threadid().eq.MASTER)
      oprint = omaster.and.control_print(print_medium)


      if (.not.btdb_get(rtdb,'nwpw:cdft',mt_log,1,cdft_on)) 
     >   cdft_on = .false.

      if (cdft_on) then
         dielec = 1.0d0

         nion = ion_nion()
         value = BA_alloc_get(mt_dbl,nion,'ucdft',ucdft(2),ucdft(1))
         value = value.and.
     >           BA_alloc_get(mt_dbl,nion,'qcdft',qcdft(2),qcdft(1))
         if (.not. value) 
     >     call errquit('nwpw_cdft_init:out of heap memory',0,MA_ERR)

         if (.not.btdb_get(rtdb,'nwpw:qcdft',mt_dbl,
     >                     nion,dbl_mb(qcdft(1)))) then
            do ii=1,nion
               dbl_mb(qcdft(1)+ii-1) = 0.0d0
            end do
         end if

         if (.not.btdb_get(rtdb,'nwpw:ucdft',mt_dbl,
     >                     nion,dbl_mb(ucdft(1)))) then
            do ii=1,nion
               dbl_mb(ucdft(1)+ii-1) = 0.0d0
            end do
         end if



         if (oprint) then
            write(luout,*)
            write(luout,*) "initializing CDFT Coupling"
            write(luout,*) "---------------------------"
            write(luout,*)
            do ii=1,nion
               write(luout,101) ion_atom(ion_katm(ii)),ii,
     >                          dbl_mb(ucdft(1)+ii-1)
            end do
            write(luout,*)
         end if

      endif

  101 format(1x,"CDFT Potential ",A4," - u(",I1,") = ",F14.9,
     >       " Hartree")
         
      return
      end


*     *******************************************
*     *                                         *
*     *             nwpw_cdft_end               *
*     *                                         *
*     *******************************************
      subroutine nwpw_cdft_end()
      implicit none
      
#include "bafdecls.fh"
#include "btdb.fh"
#include "errquit.fh"

*     **** nwpw_cdft common blocks ****
      logical cdft_on
      integer ucdft(2),qcdft(2),rtdb
      real*8  dielec
      common /nwpw_cdft_blk/ ucdft,qcdft,dielec,rtdb,cdft_on

*     **** local variables ****
      logical value
      integer nion

*     **** external functions ****
      integer  ion_nion
      external ion_nion

c     *** Clean up cdft stuff, if necessary ***
      if (cdft_on) then
         nion = ion_nion()
         value = btdb_put(rtdb,'nwpw:ucdft',mt_dbl,
     >                           nion,dbl_mb(ucdft(1))) 
         value = value.and.btdb_put(rtdb,'nwpw:qcdft',mt_dbl,
     >                           nion,dbl_mb(qcdft(1))) 
         if (.not.value) 
     >      call errquit('nwpw_cdft_end:writing ucdft',0,MA_ERR)

         cdft_on = .false.
         value = BA_free_heap(qcdft(2))
         value = value.and.BA_free_heap(ucdft(2))
         if (.not.value) 
     >      call errquit('nwpw_cdft_end:heap memory',0,MA_ERR)
      endif
      return
      end

*     *******************************************
*     *                                         *
*     *              nwpw_cdft_on               *
*     *                                         *
*     *******************************************
      logical function nwpw_cdft_on()
      implicit none

*     **** nwpw_cdft common blocks ****
      logical cdft_on
      integer ucdft(2),qcdft(2),rtdb
      real*8  dielec
      common /nwpw_cdft_blk/ ucdft,qcdft,dielec,rtdb,cdft_on

      nwpw_cdft_on = cdft_on
      return
      end


*     *******************************************
*     *                                         *
*     *              nwpw_cdft_u_ptr            *
*     *                                         *
*     *******************************************
      integer function nwpw_cdft_u_ptr()
      implicit none

*     **** nwpw_cdft common blocks ****
      logical cdft_on
      integer ucdft(2),qcdft(2),rtdb
      real*8  dielec
      common /nwpw_cdft_blk/ ucdft,qcdft,dielec,rtdb,cdft_on

      nwpw_cdft_u_ptr = ucdft(1)
      return
      end 


*     *******************************************
*     *                                         *
*     *              nwpw_cdft_q_ptr            *
*     *                                         *
*     *******************************************
      integer function nwpw_cdft_q_ptr()
      implicit none

*     **** nwpw_cdft common blocks ****
      logical cdft_on
      integer ucdft(2),qcdft(2),rtdb
      real*8  dielec
      common /nwpw_cdft_blk/ ucdft,qcdft,dielec,rtdb,cdft_on

      nwpw_cdft_q_ptr = qcdft(1)
      return
      end 


*     *******************************************
*     *                                         *
*     *            nwpw_cdft_print              *
*     *                                         *
*     *******************************************
      subroutine nwpw_cdft_print(nga,nion_qm,qgaus)
      implicit none
      integer nga,nion_qm
      real*8 qgaus(nga,nion_qm)

#include "bafdecls.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"

*     **** nwpw_cdft common blocks ****
      logical cdft_on
      integer ucdft(2),qcdft(2),rtdb
      real*8  dielec
      common /nwpw_cdft_blk/ ucdft,qcdft,dielec,rtdb,cdft_on

*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)
      logical oprint,omaster
      integer ii,nion,j,indx
      real*8  Gsolv,sum

*     **** external functions ****
      character*4 ion_atom
      external    ion_atom
      integer     ion_nion,ion_katm,ion_rion_ptr
      external    ion_nion,ion_katm,ion_rion_ptr
      integer     ion_katm_qm
      external    ion_katm_qm
      real*8      nwpw_cdft_energy0,psp_zv
      external    nwpw_cdft_energy0,psp_zv
      integer     Parallel_threadid
      external    Parallel_threadid
      logical     control_print
      external    control_print

      call Parallel_taskid(taskid)
      omaster =  (taskid.eq.MASTER).and.(Parallel_threadid().eq.MASTER)
      oprint = omaster.and.control_print(print_medium)

      if (cdft_on) then
         nion  = ion_nion()
         do ii=1,nion_qm
            sum = 0.0d0
            do j=1,nga
              sum = sum + qgaus(j,ii)
            end do
            dbl_mb(qcdft(1)+ii-1) = psp_zv(ion_katm_qm(ii))-sum
         end do
         do ii=nion_qm+1,nion
            dbl_mb(qcdft(1)+ii-1) = psp_zv(ion_katm(ii))
         end do

         Gsolv = nwpw_cdft_energy0(nion,
     >                             dbl_mb(ucdft(1)),
     >                             dbl_mb(qcdft(1)))

         if (oprint) then
            write(luout,*)
            write(luout,*) "Generalized CDFT Solvation"
            write(luout,*) "--------------------------"
            write(luout,*)
            write(luout,'(2x,A,F8.2)') "Dielectric constant -eps- = ",
     >                              dielec
            write(luout,*)
            do ii=1,nion
               write(luout,101) ion_atom(ion_katm(ii)),ii,
     >                          dbl_mb(ucdft(1)+ii-1),
     >                          dbl_mb(qcdft(1)+ii-1)
            end do
            write(luout,*)
            write(luout,102) Gsolv,Gsolv*27.2116d0*23.06d0
            write(luout,*)
         end if

      end if

  101 format(2x,"CDFT Potential ",A4," - u(",I1,") = ",F8.3,
     >       " Hartrees - atom charge = ",F8.3)
  102 format(2x,"Gsolvation(Born) = ",E14.6," (",F8.3," kcal/mol)")
      return
      end


*     *******************************************
*     *                                         *
*     *            nwpw_cdft_energy0            *
*     *                                         *
*     *******************************************
      real*8 function nwpw_cdft_energy0(nion,u,q)
      implicit none
      integer nion
      real*8  u(nion),q(nion)

*     **** local variables ****
      integer MASTER,taskid,np
      parameter (MASTER=0)
      integer ii,jj,itask
      real*8  Gsolv,screen

!$OMP MASTER
      Gsolv = 0.0d0
      do ii=1,nion
         Gsolv = Gsolv + q(ii)*u(ii)
      end do
!$OMP END MASTER

      nwpw_cdft_energy0 = Gsolv
      return
      end

*     *******************************************
*     *                                         *
*     *            nwpw_cdft_energy             *
*     *                                         *
*     *******************************************
      real*8 function nwpw_cdft_energy()
      implicit none

#include "bafdecls.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"

*     **** nwpw_cdft common blocks ****
      logical cdft_on
      integer ucdft(2),qcdft(2),rtdb
      real*8  dielec
      common /nwpw_cdft_blk/ ucdft,qcdft,dielec,rtdb,cdft_on

*     **** external functions ****
      integer  ion_nion
      external ion_nion
      real*8   nwpw_cdft_energy0
      external nwpw_cdft_energy0
      
      nwpw_cdft_energy = nwpw_cdft_energy0(ion_nion(),
     >                             dbl_mb(ucdft(1)),
     >                             dbl_mb(qcdft(1)))
      return
      end

*     *******************************************
*     *                                         *
*     *            nwpw_cdft_dVdq               *
*     *                                         *
*     *******************************************
      subroutine nwpw_cdft_dVdq(nion,q,ucdft,u)
      implicit none
      integer nion
      real*8  q(nion),ucdft(nion)
      real*8  u(nion)


      call Parallel_shared_vector_copy(.true.,nion,ucdft,u)
      return
      end
