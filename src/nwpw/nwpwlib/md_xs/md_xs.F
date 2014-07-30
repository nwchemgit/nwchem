c $Id: md_xs.F 25353 2014-03-17 14:12:25Z jhammond $


c     **********************************************
c     *                                            *
c     *                md_xs_analysis              *
c     *                                            *
c     **********************************************
      subroutine md_xs_analysis(rtdb,filename)
      implicit none
      integer rtdb
      character*(*) filename

#include "mafdecls.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "stdio.fh"

*     **** local variables ****
      integer taskid_i,taskid_j,npi,npj
      integer taskid,np,MASTER
      parameter (MASTER=0)

      logical oprint,value
      character*20 codeversion
      character*2 tag2
      character*4 tag4
      integer nframes,natoms,stride,natoms1
      integer iii,i,k,icount,icpu,pto,pfrom
      real*8  x,y,z,vx,vy,vz,time1,volume,unita(3,3)

      integer nkf
      integer kf(2),chif(2),rion(2),tag(2)

*     **** external functions ****
      integer  Parallel2d_convert_taskid_i
      external Parallel2d_convert_taskid_i
       

      call Parallel_np(np)
      call Parallel_taskid(taskid)
      oprint = (taskid.eq.MASTER)

      call Parallel2d_np_j(npi)
      call Parallel2d_np_j(npj)
      call Parallel2d_taskid_j(taskid_i)
      call Parallel2d_taskid_j(taskid_j)

#ifdef FEFF
      call feff_codeversion(codeversion)

      if (.not.rtdb_get(rtdb,'md_xs:stride',mt_int,1,stride))
     >   stride = 1

      if (oprint) then
         write(luout,*) 
     >   "==== ","MD XS analysis: ",codeversion," ===="
         write(luout,*)
         write(luout,*) 'MOTION filename:  ',trim(filename)
      end if

      if (taskid.eq.MASTER) then
         call MOTION_nframes(filename,nframes,natoms,volume,unita)
      end if
      call Parallel_Brdcst_ivalue(0,nframes)
      call Parallel_Brdcst_ivalue(0,natoms)
      call Parallel_Brdcst_value(0,volume)
      call Parallel_Brdcst_values(0,9,unita)

      if (oprint) then
         write(luout,*) 'number of frames: ',nframes
         write(luout,*) 'stride:           ',stride
         write(luout,*) 'natoms:           ',natoms
         write(luout,*) 'volume:           ',volume
      end if

*     **** allocate memory from stack ****
      nkf = 601
      value =  MA_push_get(mt_dbl,nkf,'kf',kf(2),kf(1))
      value = value.and.
     >        MA_push_get(mt_dbl,nkf,'chif',chif(2),chif(1))
      value = value.and.
     >        MA_push_get(mt_dbl,3*natoms,'rion',rion(2),rion(1)) 
      value = value.and.
     >        MA_push_get(mt_byte,2*natoms,'tag',tag(2),tag(1)) 
      if (.not.value) 
     >   call errquit("md_xs_analysis: out of stack",1,MA_ERR)


      if (taskid_i.eq.MASTER)
     >   open(unit=19,file=filename,form='formatted')

      call dcopy(nkf,0.0d0,0,dbl_mb(kf(1)),1)
      call dcopy(nkf,0.0d0,0,dbl_mb(chif(1)),1)
      icount = 0
      icpu   = 0
      do k=1,nframes
         if (taskid_i.eq.MASTER) then
            read(19,*) time1,natoms1,volume,unita
            do i=1,natoms1
               tag2 = '  '
               read(19,*) iii,tag2,tag4,x,y,z,vx,vy,vz
               byte_mb(tag(1)+2*(i-1))   = tag2(1:1)
               byte_mb(tag(1)+2*(i-1)+1) = tag2(2:2)
               dbl_mb(rion(1)+3*(i-1))   = x
               dbl_mb(rion(1)+3*(i-1)+1) = y
               dbl_mb(rion(1)+3*(i-1)+2) = z
            end do
         end if

         if (mod(k,stride).eq.0) then
            if (mod(icount,npj).eq.taskid_j) then
                 if (taskid_i.eq.MASTER) then
                    if (icpu.eq.MASTER) then
                       !*** compute md_xs with natoms1,tag,rion +--> kf,chif ***
                       call md_xs(rtdb,
     >                         natoms1,byte_mb(tag(1)),dbl_mb(rion(1)),
     >                         nkf,dbl_mb(kf(1)),dbl_mb(chif(1)))
                    else
                       !*** send natoms1,tag,rion to icpu ***
                       pto = Parallel2d_convert_taskid_i(icpu)
                       call Parallel_send_characters(pto,2*k,
     >                                        2*natoms1,byte_mb(tag(1)))
                       call Parallel_send_value(pto,2*k+1,
     >                                        3*natoms1,dbl_mb(rion(1)))
                    end if
                 else
                    if (icpu.eq.taskid_i) then
                       !*** receive natoms1,tag,rion from MASTER ***
                       pfrom = Parallel2d_convert_taskid_i(MASTER)
                       call Parallel_recv_characters(pfrom,2*k,
     >                                        2*natoms1,byte_mb(tag(1)))
                       call Parallel_recv_values(pfrom,2*k+1,
     >                                        3*natoms1,dbl_mb(rion(1)))

                       !*** compute md_xs with natoms1,tag,rion +--> kf,chif ***
                       call md_xs(rtdb,
     >                          natoms1,byte_mb(tag(1)),dbl_mb(rion(1)),
     >                          nkf,dbl_mb(kf(1)),dbl_mb(chif(1)))
                    end if
                 end if
                 icpu = mod(icpu+1,npi)
            end if
            icount = icount + 1
         end if
      end do
      if (taskid_i.eq.MASTER) close(19)

      x = 1.0d0/dble(icount)
      call Parallel_Vector_SumAll(nkf,dbl_mb(chif(1)))
      call dscal(nkf,x,dbl_mb(chif(1)),1)

c     **** deallocate stack ****
      value =           MA_pop_stack(tag(2))
      value = value.and.MA_pop_stack(rion(2))
      value = value.and.MA_pop_stack(chif(2))
      value = value.and.MA_pop_stack(kf(2))
      if (.not.value)
     >   call errquit("md_xs_analysis: popping stack",3,MA_ERR)

#else
      if (oprint) then
         write(luout,*) "==== ",
     >     "FEFF library not inlucded in compile: ",
     >     "MD XS analysis not performed"," ===="
      end if
#endif
      return
      end



#ifdef FEFF
c     **********************************************
c     *                                            *
c     *                md_xs                       *
c     *                                            *
c     **********************************************
      subroutine md_xs(rtdb,nion,symbols,rion,nkf,kf,chi)
      implicit none
      integer rtdb
      integer nion
      character*2 symbols(*)
      real*8      rion(3,*)
      integer nkf
      real*8 kf(*)
      real*8 chi(*)

#include "mafdecls.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "stdio.fh"

*     **** local variables ****
      integer taskid,MASTER,np
      parameter (MASTER=0)

      integer bigstr(2)

*     **** external functions ****


      call Parallel_taskid(taskid)
      call Parallel_np(np)


      return
      end


#endif

