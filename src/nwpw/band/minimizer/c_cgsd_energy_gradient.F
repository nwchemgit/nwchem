c     $Id: c_cgsd_energy_gradient.F,v 1.1 2002-01-20 23:49:22 bylaska Exp $                       

*     *******************************
*     *								*
*     *	  c_cgsd_energy_gradient	*
*     *								*
*     *******************************

      subroutine c_cgsd_energy_gradient(G1)
      implicit none
      real*8 G1(3,*)

      logical allow_translation
      integer MASTER
      parameter (MASTER=0)
      integer i,k,taskid,nion
      real*8  GG
      real*8  fmx,fmy,fmz
      real*8  fmx2,fmy2,fmz2

*     **** external functions ****
      logical     cpsp_semicore
      logical     control_allow_translation
      character*2 ion_aname
      integer     ion_katm,ion_nion,control_version
      real*8      ion_rion
      external cpsp_semicore
      external control_allow_translation
      external ion_aname,pspw_qmmm_aname
      external ion_katm,ion_nion,control_version
      external ion_rion


      allow_translation = control_allow_translation()
      nion = ion_nion()

      call dcopy(3*nion,0.0d0,0,G1,1)

      call cpsi_1force(G1)
c      if (cpsp_semicore(0)) call c_electron_semicoreforce(G1)

      call c_ewald_f(G1)


      if (.not.allow_translation) then
        call center_F_mass(G1,fmx,fmy,fmz)
        do i=1,nion
         G1(1,i) = G1(1,i) - fmx
         G1(2,i) = G1(2,i) - fmy
         G1(3,i) = G1(3,i) - fmz
        end do
      end if
      call center_F_mass(G1,fmx2,fmy2,fmz2)

      GG = 0.0d0
      do i=1,nion
         GG = GG + G1(1,i)**2 + G1(2,i)**2 + G1(3,i)**2
      end do

      call Parallel_taskid(taskid)
      if (taskid.eq.MASTER) then
        write(6,1300)
        write(6,1301)
        write(6,1302)
        write(6,1304)
        if (.not.allow_translation) write(6,1400) fmx,fmy,fmz
        write(6,1304)
        write(6,1410)
        write(6,1420)
        write(6,1190)(i,ion_aname(I),
     >                  (ion_rion(K,I),K=1,3),I=1,ion_nion())
        write(6,1421)
        write(6,1190)(i,ion_aname(I),
     >                  (G1(K,I),K=1,3),I=1,ion_nion())

        write(6,1210) fmx2,fmy2,fmz2  
        write(6,1425)
        write(6,1426) dsqrt(GG)
      end if

c     call dscal(3*nion,(-1.0d0),G1,1)

      return
 1190 FORMAT(5X, I4, A3  ,' (',3F11.5,' )')
 1210 FORMAT(5X,' C.O.M.',' (',3F11.5,' )')
 1300 FORMAT(//'========================')
 1301 FORMAT(  '= Gradient calculation =')
 1302 FORMAT(  '========================')
 1304 FORMAT(/)
 1400 FORMAT('Translation force removed: (',3F11.5,')')
 1410 FORMAT(10X,'=============  Ion Gradients =================')
 1425 FORMAT(10X,'===============================================')
 1426 FORMAT(10X,'|F|=',E21.10//)
 1420 FORMAT( ' Ion Positions:')
 1421 FORMAT( ' Ion Forces:')
      end

