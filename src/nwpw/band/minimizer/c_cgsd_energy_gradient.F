c
c     $Id$         c

*     *******************************
*     *				    *
*     *	  c_cgsd_energy_gradient    *
*     *				    *
*     *******************************

      subroutine c_cgsd_energy_gradient(G1)
      implicit none
      real*8 G1(3,*)

#include "stdio.fh"
#include "util.fh"

*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)

      logical allow_translation,oprint,hprint
      integer i,k,nion
      real*8  GG,fmax,fatom
      real*8  fmx,fmy,fmz
      real*8  fmx2,fmy2,fmz2

*     **** external functions ****
      logical     cpsp_semicore,control_print
      logical     control_allow_translation,ion_q_FixIon
      character*4 ion_aname
      integer     ion_katm,ion_nion,control_version
      real*8      ion_rion
      external cpsp_semicore,control_print
      external control_allow_translation,ion_q_FixIon
      external ion_aname
      external ion_katm,ion_nion,control_version
      external ion_rion
      logical  ion_disp_on
      external ion_disp_on


      allow_translation = control_allow_translation()
      nion = ion_nion()

      call dcopy(3*nion,0.0d0,0,G1,1)

      call cpsi_1force(G1)
      if (cpsp_semicore(0)) call c_electron_semicoreforce(G1)

      call ewald_f(G1)
      if (ion_disp_on()) call ion_disp_force(G1)

*     **** remove ion forces using ion_FixIon ****
      call ion_FixIon(G1)


      if (.not.allow_translation) then
        call center_F_mass(G1,fmx,fmy,fmz)
        do i=1,nion
         G1(1,i) = G1(1,i) - fmx
         G1(2,i) = G1(2,i) - fmy
         G1(3,i) = G1(3,i) - fmz
        end do
      end if
      call center_F_mass(G1,fmx2,fmy2,fmz2)

      GG   = 0.0d0
      fmax = 0.0d0
      do i=1,nion
         GG = GG + G1(1,i)**2 + G1(2,i)**2 + G1(3,i)**2
         fatom= dsqrt(G1(1,i)**2 + G1(2,i)**2 + G1(3,i)**2)
         if (fatom.gt.fmax) fmax = fatom
      end do

      call Parallel_taskid(taskid)
      hprint = ((taskid.eq.MASTER).and.control_print(print_high))
      oprint = ((taskid.eq.MASTER).and.control_print(print_medium))


      if (oprint) then
        write(luout,1300)
        write(luout,1301)
        write(luout,1302)
        write(luout,1304)
        if (.not.allow_translation) write(luout,1400) fmx,fmy,fmz
        write(luout,1304)
        write(luout,1410)
        if (hprint) then
           write(luout,1420)
           do I=1,ion_nion()
             if (ion_q_FixIon(I)) then
              write(6,1191) I,ion_aname(I),(ion_rion(K,I),K=1,3)
             else
              write(6,1190) I,ion_aname(I),(ion_rion(K,I),K=1,3)
             end if
           end do
        end if
        write(luout,1421)
        write(luout,1190)(i,ion_aname(I),
     >                  (G1(K,I),K=1,3),I=1,ion_nion())

        write(luout,1210) fmx2,fmy2,fmz2  
        write(luout,1425)
        write(luout,1426) dsqrt(GG),
     >                    dsqrt(GG)/dble(nion),
     >                    fmax,fmax*(27.2116d0/0.529177d0)
      end if

c     call dscal(3*nion,(-1.0d0),G1,1)

      return
 1190 FORMAT(5X, I4, A5  ,' (',3F11.5,' )')
 1191 FORMAT(5X, I4, A5  ,' (',3F11.5,' ) - fixed')
 1210 FORMAT(5X,'   C.O.M.',' (',3F11.5,' )')
 1300 FORMAT(//'========================')
 1301 FORMAT(  '= Gradient calculation =')
 1302 FORMAT(  '========================')
 1304 FORMAT(/)
 1400 FORMAT('Translation force removed: (',3F11.5,')')
 1410 FORMAT(10X,'=============  Ion Gradients =================')
 1420 FORMAT( ' Ion Positions:')
 1421 FORMAT( ' Ion Forces:')
 1425 FORMAT(10X,'===============================================')
 1426 FORMAT(10X,'|F|       =',E15.6,
     >      /10x,'|F|/nion  =',E15.6,
     >      /10x,'max|Fatom|=',E15.6,1x,'(',F8.3,'eV/Angstrom)'//)
      end

