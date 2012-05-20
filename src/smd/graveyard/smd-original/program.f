      PROGRAM moldyn

      implicit none

      include 'p_input.inc'
      include 'p_array.inc'

      integer iseed,istep
      integer ntype,nprnt,nstep,nequl
      integer ncons,ntcons,consatm
      integer nbond,ntbond,bondatm
      integer nshel,ntshel,shelatm

      real*8 ekin,ecoul,eshrt,ebond,eshel,ewald1
      real*8 tstep,ivv,etime1,etime2

      logical lveloc,lupdate

      dimension ivv(mxatms,3)
      dimension consatm(mxcons,2),bondatm(mxbond,2),shelatm(mxshel,2)

      etime1=0.0
      etime2=0.0
      call cpu_time(etime1)

      open(unit=output,file='OUTPUT')

      call md_inpt(iseed,tstep,nstep,nequl,nprnt,ntype,ncons,consatm,
     $    nbond,bondatm,nshel,shelatm,lveloc)

      call md_init(iseed,ntype,ncons,ntcons,consatm,nbond,ntbond,
     $    bondatm,nshel,ntshel,shelatm,lveloc,ewald1)

      call cpu_time(etime2)
      write(output,'(/,a,f20.3)')'Set-up CPU time : ',(etime2-etime1)

      call write_velo(70)
      do istep=1,nstep

       call verlt_test(tstep,ivv,lupdate)

       if(lupdate)call list_verlt()

       call md_frce(ntype,ecoul,eshrt,ebond,ntbond,eshel,ntshel,ewald1)

       if(ntcons.eq.0)call inte_leapf(tstep,ekin)
       if(ntcons.gt.0)call inte_shake(tstep,ntcons,ekin)
 
c       if(istep.le.nequl)call md_scle(ntshel)

       call print_output(istep,nprnt,ekin,ecoul,eshrt,ebond,eshel)

      enddo

      call print_final()

      call cpu_time(etime2)
      write(output,'(/,a,f20.3)')'Total CPU time : ',(etime2-etime1)

      END
c $Id$
