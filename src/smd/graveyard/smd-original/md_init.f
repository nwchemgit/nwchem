      SUBROUTINE md_init(iseed,ntype,ncons,ntcons,consatm,nbond,ntbond,
     $    bondatm,nshel,ntshel,shelatm,lveloc,ewald1)

      implicit none

      include 'p_array.inc'
      include 'p_const.inc'
      include 'p_input.inc'
      include 'cm_atom.inc'
      include 'cm_latt.inc'
      include 'cm_temp.inc'

      integer iseed,ntype
      integer ncons,ntcons,consatm
      integer nbond,ntbond,bondatm
      integer nshel,ntshel,shelatm

      real*8 x,det,ewald1

      logical lveloc

      dimension consatm(mxcons,2),bondatm(mxbond,2),shelatm(mxshel,2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call tool_randm(iseed,x)
      call tool_invrt(latt,rlatt,det)
      call tool_volme(latt,vol)
      call tool_rebox(natms,mxatms,latt,rlatt,ccc)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call ewald_setp()
c Self-interaction energy
      call ewald_self(ewald1)
      write(*,*) "ewald self",ewald1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Create exclude list
      call list_excld(ntype)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Create constraints list
      call list_const(ncons,ntcons,consatm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Create bonds list
      call list_bonds(nbond,ntbond,bondatm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Create shells list
      call list_shell(nshel,ntshel,shelatm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Temperature parameters
      degfree=dble(3*(natms-ntshel)-3-ntcons)
      targetke=degfree*temp*boltzmann*0.5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Initial random velocities
      if(.not.lveloc)call init_velo(iseed)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      return

      END
c $Id$
