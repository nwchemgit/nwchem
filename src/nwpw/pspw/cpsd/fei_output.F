*
* $Id: fei_output.F,v 1.1 2002-10-09 21:08:00 bylaska Exp $
*
      subroutine fei_init()
      implicit none

*     **** local variables ****
      integer taskid,MASTER
      parameter(MASTER=0)

*     **** external functions ****
      logical  control_Fei
      external control_Fei


      if (control_Fei()) then

        call Parallel_taskid(taskid)
        if (taskid.eq.MASTER) open(unit=35,file='FEI',form='formatted')

      end if
      return
      end

      subroutine fei_end()
      implicit none

*     **** local variables ****
      integer taskid,MASTER
      parameter(MASTER=0)

*     **** external functions ****
      logical  control_Fei
      external control_Fei

      if (control_Fei()) then

        call Parallel_taskid(taskid)
        if (taskid.eq.MASTER) close(unit=35)

      end if
      return
      end


      subroutine fei_output(E,ion_fion)
      implicit none
      real*8 E
      integer nion
      real*8  ion_fion(3,*)

*     **** local variables ****
      integer taskid,MASTER
      parameter(MASTER=0)

      integer ii,nion

*     **** external functions ****
      logical  control_Fei
      integer  ion_nion
      real*8   lattice_unita,ion_rion
      external control_Fei
      external ion_nion
      external lattice_unita,ion_rion

      if (control_Fei()) then

        call Parallel_taskid(taskid)

        if (taskid.eq.MASTER) then
        nion = ion_nion()
        write(35,100) nion
        write(35,101) E
        write(35,102) lattice_unita(1,1),
     >                lattice_unita(2,2),
     >                lattice_unita(3,3)
        do ii=1,nion
          write(35,103) ii,
     >                  ion_rion(1,ii),ion_rion(2,ii),ion_rion(3,ii),
     >                  ion_fion(1,ii),ion_fion(2,ii),ion_fion(3,ii)
        end do
        call util_flush(35)
        end if

  100   format(I8)
  101   format(F20.6)
  102   format(3F10.5)
  103   format(I5,6F10.5)

      end if

      return
      end
