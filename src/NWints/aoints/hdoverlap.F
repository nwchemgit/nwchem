      subroutine hdoverlap ( H, l_buf, l_scr, basis, geom, nproc,
     &                          nat, max1e, mem1, hdol, nmot)
C LB based loosely on grad1.  Program to calculate half-derivative overlap
C integrals, 1/2[<i|j^x>-<i^x|j>].  This capability already exists
C within intd_1eov, it just needed some teasing out.


      implicit none

#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
#include "bas.fh"
#include "rtdb.fh"
#include "sym.fh"
#include "bq_params.fh"
!#include "cdft.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "nwc_const.fh"
#include "basP.fh"

C-------------------------parameters--------------------------------
      INTEGER :: basis, lbus, max1e, mem1, nproc, nat, nbf
      INTEGER :: nmot
      INTEGER :: l_buf, l_scr, geom
      INTEGER :: hdol !handle for hdol GA

      double precision :: H ! integral derivatives
      DOUBLE PRECISION :: scr ! scratch
      DOUBLE PRECISION :: SQHDOL (3*nat,nmot,nmot) 

      dimension H ( l_buf ), scr(l_scr)

      logical oskel   ! symmetry?

C-------------------------local variables--------------------------

      integer ijatom, next, iat1, iat2, iat3, ish1, ish2,
     $     iab1f, iab1l, iab2f, iab2l, iac1f, iac1l, iac2f, iac2l,
     $     if1, il1, if2, il2,
     $     icart, ic, nint, ip1, ip2

      double precision crd1, crd2 ! atomic coordinates
      dimension crd1(3), crd2(3)

      integer idatom
      dimension idatom(2)

      double precision qfac

      logical status, pointforce,dobq

      character*16 name

      integer bq_ncent
      integer i_qbq,i_cbq
      double precision r12
c
      integer nxtask, task_size
      external nxtask
      LOGICAL :: hlfdrv !true = calculate half-derivative overlap integrals
      INTEGER :: atmidx1, atmidx2, nmss2, idx3  ! indices the SQHDOL matrices
      INTEGER :: contf, contl, shlf, shll ! more indices
      INTEGER :: orbstrt(nat) ! gives the orbital label where each atom
                              ! center starts.  There may be a better
                              ! built-in way to achieve this
      INTEGER :: i,j,k ! counters
      INTEGER :: itype, jtype ! angular momentum of contraction
      INTEGER :: nintsph      ! nint after conversion to spherical
      INTEGER :: dum1, dum2, dum3 ! dummy variables
      INTEGER :: dims(3)
      INTEGER, PARAMETER :: istrt(3) = (/1,1,1/)


c LB
!      WRITE(LUOUT,*)"LB Debug, now in hdoverlap"
      SQHDOL=0.0D0
!      WRITE(LUOUT,*)"LB Debug, iab1f=",iab1f,"iab2f=",iab2f
!      WRITE(LUOUT,*)"nmot=",nmot
      nbf=nmot ! convert to INTEGER*8
      orbstrt=0! construct a vector containing the orbital
               ! index at which each center starts.  There
               ! could be an easier way to do this.  MAKE SURE THIS IS
               ! COUNTING BASED ON SPHERICAL, NOT CARTESIAN!
      orbstrt(1)=1
      DO i = 1, nat-1
        orbstrt(i+1)=orbstrt(i)
        status = bas_ce2cnr(basis,i,contf,contl) 
        DO j = contf, contl
          status = bas_cn2bfr(basis,j,shlf,shll)
          orbstrt(i+1)=orbstrt(i+1)+shll-shlf+1
        ENDDO
      ENDDO
C LB
      task_size = 1
      status = rtdb_parallel(.true.) ! Broadcast reads to all processes

      pointforce = geom_include_bqbq(geom)
      dobq = geom_extbq_on()
      call hf_print_set(1)

      ijatom = -1
      next = nxtask(nproc,task_size)
      ! Loop Through Bra atom

      do 90, iat1 = 1, nat

        ! Loop through ket atom

        do 80, iat2 = 1, iat1

          ijatom = ijatom + 1
C LB
!          WRITE(LUOUT,*)"LB Debug, iat1=",iat1,"iat2=",iat2
!     &                 ,"ijatom=",ijatom
C LB
          if ( ijatom .eq. next ) then

            status = bas_ce2bfr(basis,iat1,iab1f,iab1l)
            status = bas_ce2bfr(basis,iat2,iab2f,iab2l)

            if (iab1f.le.0 .or. iab2f.le.0) then
c     
c     At least one center has no functions on it ... next atom
c
              goto 1010
            endif

            if (oskel) then
               if (.not. sym_atom_pair(geom, iat1, iat2, qfac))
     $              goto 1010
            else
               qfac = 1.0d0
            endif

            status = bas_ce2cnr(basis,iat1,iac1f,iac1l)
            status = bas_ce2cnr(basis,iat2,iac2f,iac2l)

!            call ga_get (g_dens, iab1f,iab1l,iab2f,iab2l,dens,max_at_bf)
!            call ga_get(g_wdens,iab1f,iab1l,iab2f,iab2l,wdens,max_at_bf)

            ! Loop through bra atom shells

            do 70, ish1 = iac1f, iac1l
              if ( iat1.eq.iat2 ) iac2l = ish1

              ! Loop through ket atom shells

              do 60, ish2 = iac2f, iac2l
C LB
!                WRITE(LUOUT,*)"LB Debug, ishl=",ish1,"ish2=",ish2
C LB

C               shell block in atomic (D/Dw)-matrix block
                status = bas_cn2bfr(basis,ish1,if1,il1)
C LB
!                WRITE(LUOUT,*)"LB hdoverlap, if1=",if1,"il1=",il1
C LB
                if1 = if1 - iab1f + 1
                il1 = il1 - iab1f + 1
                status = bas_cn2bfr(basis,ish2,if2,il2)
C LB
!                WRITE(LUOUT,*)"LB hdoverlap, if2=",if2,"il2=",il2
C LB
                if2 = if2 - iab2f + 1
                il2 = il2 - iab2f + 1

                nint = ( il1 - if1 + 1 ) * ( il2 - if2 + 1 )

C               half-derivative overlap derivatives
C LB
                hlfdrv=.true.
!                WRITE(LUOUT,*)"LB hd, iat1=",iat1,"iat2=",iat2
!                WRITE(LUOUT,*)"LB hd, ish1=",ish1,"ish2=",ish2
C LB
                call intd_1eovQ(basis,ish1,basis,ish2,l_scr,scr,
     &               l_buf,H,idatom,hlfdrv)
C LB
!                WRITE(LUOUT,*)"LB Debug, back from intd_1eovQ"

                status=bas_continfo(basis,ish1,itype,dum1,dum2,dum3)
                status=bas_continfo(basis,ish2,jtype,dum1,dum2,dum3)
                nintsph=(2*itype+1)*(2*jtype+1)
!                WRITE(LUOUT,*)"nintsph=",nintsph

!                WRITE(LUOUT,*)"Ova="
!                WRITE(LUOUT,*)H(1:nintsph*6)

C LB    We now have H (BO2i) which includes the integrals organized by
C       subshell and dx dy dz wrt both atom centers involved.  We now
C       need to get these into a matrix in the AO space, organized by
C       atom center derivative. I am not going to do this by symmetry
C       block right now since the initial SONADCT capability will only
C       be available in C1 symmetry.
        
!                WRITE(LUOUT,*)"LB hdoverlap, orbstrt=",orbstrt
!                WRITE(LUOUT,*)"Bra Atom=",iat1,"shell=",ish1
!                WRITE(LUOUT,*)"Ket Atom=",iat2,"shell=",ish2
                atmidx1=iat1*3-2 !index of atom center derivative
                atmidx2=iat2*3-2 !index of atom center derivative
!                WRITE(LUOUT,*)"if1=",if1,"il1=",il1
!                WRITE(LUOUT,*)"if2=",if2,"il2=",il2
!                WRITE(LUOUT,*)"iac1f=",iac1f,"iac1l=",iac1l
!                WRITE(LUOUT,*)"iac2f=",iac2f,"iac2l=",iac2l

                j=0
                IF (bas_spherical(basis)) THEN
                  il1=if1+2*itype
                  il2=if2+2*jtype
                  nint=nintsph
!                  WRITE(LUOUT,*)"in spherical region,
!     &                           il1=",il1,"il2=",il2
                ENDIF
                DO i = if1, il1
                  j=j+1
                  nmss2=il2-if2+1
                  idx3=(j-1)*nmss2
                  IF (iat1.EQ.iat2) THEN
                    !basis functions from same atom
                    SQHDOL(atmidx1,  orbstrt(iat1)+i-1,
     &                orbstrt(iat2)+if2-1:orbstrt(iat2)+il2-1)=
     &                -H(idx3+1:idx3+nmss2+1            )!x1
                    SQHDOL(atmidx1+1,orbstrt(iat1)+i-1,
     &                orbstrt(iat2)+if2-1:orbstrt(iat2)+il2-1)=
     &                -H(idx3+nint+1:idx3+nint+nmss2    )!y1
                    SQHDOL(atmidx1+2,orbstrt(iat1)+i-1,
     &                orbstrt(iat2)+if2-1:orbstrt(iat2)+il2-1)=
     &                -H(idx3+2*nint+1:idx3+2*nint+nmss2)!z1
                  ELSE
                    !basis functions from different atoms
                    SQHDOL(atmidx1,  orbstrt(iat1)+i-1,
     &                orbstrt(iat2)+if2-1:orbstrt(iat2)+il2-1)=
     &                H(idx3+1:idx3+nmss2+1            )!x1
                    SQHDOL(atmidx1+1,orbstrt(iat1)+i-1,
     &                orbstrt(iat2)+if2-1:orbstrt(iat2)+il2-1)=
     &                H(idx3+nint+1:idx3+nint+nmss2    )!y1
                    SQHDOL(atmidx1+2,orbstrt(iat1)+i-1,
     &                orbstrt(iat2)+if2-1:orbstrt(iat2)+il2-1)=
     &                H(idx3+2*nint+1:idx3+2*nint+nmss2)!z1

                    SQHDOL(atmidx2 ,orbstrt(iat1)+i-1,
     &                orbstrt(iat2)+if2-1:orbstrt(iat2)+il2-1)=
     &                H(idx3+1:idx3+nmss2+1)!x2
                    SQHDOL(atmidx2+1,orbstrt(iat1)+i-1,
     &                orbstrt(iat2)+if2-1:orbstrt(iat2)+il2-1)=
     &                H(idx3+nint+1:idx3+nint+nmss2)!y2
                    SQHDOL(atmidx2+2,orbstrt(iat1)+i-1,
     &                orbstrt(iat2)+if2-1:orbstrt(iat2)+il2-1)=
     &                H(idx3+2*nint+1:idx3+2*nint+nmss2)!z2
                  ENDIF
                ENDDO

!                DO i=1,3*nat
!                  WRITE(LUOUT,*)"SQHDOL",i
!                  CALL output(SQHDOL(i,:,:),1,nmot,1,nmot,nmot,nmot,1)
!                ENDDO


 60           continue
 70         continue
!
 1010       continue
           next = nxtask(nproc,task_size)
          endif
!
 80     continue
 90   continue
      next = nxtask(-nproc,task_size)
      !copy temp SQHDOL matrix too hdol global array
      dims(1) = 3*nat
      dims(2) = nmot
      dims(3) = nmot
      CALL nga_put(hdol,istrt,dims,SQHDOL,dims(1:2))

cgk debug
!      write(*,*)'gk: leaving hdoverlap'
cgk end

      return
      end
