      logical function cpmd_qmmm_step(verlet,
     >                                rion,nion,
     >                                rcharge,qcharge,ncharge,
     >                                E,fion,fcharge)
      implicit none
      integer rtdb
      logical verlet
      real*8  rion(3,*)
      integer nion
      real*8  rcharge(3,*),qcharge(*)
      integer ncharge
      real*8  E(*)
      real*8  fion(3,*)
      real*8  fcharge(3,*)

#include "global.fh"
#include "bafdecls.fh"
#include "btdb.fh"
      

*     **** cpmd_qmmm_block ****
      integer nfft3d
      integer npack1
      integer ispin
      integer ne(2),nemax
      integer psi0(2),psi1(2),psi2(2)
      integer dn(2)
      integer Hpsi(2),psir(2)
      integer hml(2),lmd(2)
      integer icount,icount_shift
      real*8  emotion_time_shift,dt,eave,evar,eke
      logical mulliken
      common /cpmd_qmmm_block/ psi0,psi1,psi2,dn,hml,lmd,Hpsi,psir,
     >                         ispin,ne,npack1,nfft3d,nemax,
     >                         icount,icount_shift,
     >                         emotion_time_shift,dt,eave,evar,eke,
     >                         mulliken

      integer ms,ii,jj
      integer taskid,MASTER
      parameter (MASTER=0)

      logical  control_Nose
      external control_Nose

      call Parallel_taskid(taskid)


*     ***** first step using velocity ****
      call inner_loop_qmmm_step(verlet,ispin,ne,
     >             npack1,nfft3d,nemax,
     >             dcpl_mb(psi0(1)),
     >             dcpl_mb(psi1(1)),
     >             dcpl_mb(psi2(1)),
     >             dbl_mb(dn(1)),
     >             1,0,E,
     >             dbl_mb(hml(1)),dbl_mb(lmd(1)),
     >             dcpl_mb(Hpsi(1)),dcpl_mb(psir(1)),
     >             rion,nion,rcharge,qcharge,ncharge,fion,fcharge)

      if (verlet) then 


         eke = eke + E(3)
         icount = icount + 1

         if (taskid.eq.MASTER) then 

           write(6,1310) icount,E(1),E(2),E(3),E(4)
           call util_flush(6)

*          **** write out EMOTION data ****
           eave = E(25)/dble(icount+icount_shift)
           evar = E(26)/dble(icount+icount_shift)
           evar = evar - eave*eave
           if (control_Nose()) then
             write(31,1311) icount*dt + emotion_time_shift,
     >                    e(1),e(2),e(3),e(4),e(5),e(6),
     >                    e(7),e(8),e(9),e(10),
     >                    eave,evar
           else
             write(31,1311) icount*dt + emotion_time_shift,
     >                    e(1),e(2),e(3),e(4),e(5),e(6),
     >                    e(7),e(8),
     >                    eave,evar
           end if
           call util_flush(31)

*          **** write out EIGMOTION data -diagonal hml matrix ****
           if (mulliken) then
           write(32,1311) icount*dt,
     >       (( dbl_mb(hml(1)+ii-1+(ii-1)*ne(1)+(ms-1)*ne(1)*ne(1)), 
     >         ii=1,ne(ms)),ms=1,ispin)
           call util_flush(32)
           end if

*          **** write out HMOTION data - hml matrix ****
           if (mulliken) then
           write(34,1312) icount*dt,ispin
           do ms=1,ispin
             write(34,1313) ms,ne(ms),ne(ms)
             do ii=1,ne(ms)
               write(34,1311) 
     >         (dbl_mb(hml(1)+ii-1+(jj-1)*ne(1)+(ms-1)*ne(1)*ne(1)), 
     >          jj=1,ne(ms))
             end do
           end do
           call util_flush(34)
           end if

         end if


*        **** write xyz data ****
         call xyz_write()
         call MOTION_write((icount*dt))

*        **** write OMOTION data ****
         if (mulliken) call Orb_Write(dcpl_mb(psi1(1)))


      end if

      cpmd_qmmm_step = .true.
      return
 1309 FORMAT(I8,2E19.10,2E14.5,2F9.1)
 1310 FORMAT(I8,2E19.10,2E14.5,F14.2)
 1311 format(100e19.10)
 1312 format(e14.6,i3)
 1313 format(3i4)
      end



c $Id$
