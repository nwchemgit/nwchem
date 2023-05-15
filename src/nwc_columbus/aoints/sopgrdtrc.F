      subroutine sopgrdtrc
     $     ( H, lbuf, scr, lscr, wdens, 
     $     frc_sow, 
     $     g_dens, g_wdens, basis, geom, nproc, nat, 
     $     max_at_bf, oskel, 
     &     frc_sox, frc_soy, frc_soz, densx, densy, densz)
c based on grad1_so in nwdft/sodft -gsk
c keep the parallel stuff.  Sure why not?

C     spin orbit potential gradient using effective densities from 
c     columbus

      implicit none

#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
#include "bas.fh"
#include "rtdb.fh"
#include "sym.fh"

C-------------------------parameters--------------------------------
      integer lbuf, lscr,
     $     g_dens(3),        ! density matrix (summed if ROHF, UHF)
     $     g_wdens,       ! weighted density (Lagrangian)
     $     basis, geom, nproc, nat, max_at_bf

      double precision H, ! integral derivatives
     $     scr, 
     $     densx, densy, densz, ! local density block
     $     wdens,               ! local weighted density block
     $     frc_sow, frc_sox, frc_soy, frc_soz ! forces arrays

      dimension H ( lbuf ), 
     $     frc_sow(3, nat), frc_sox(3, nat), frc_soy(3,nat),
     $     frc_soz(3,nat), scr(lscr),
     $     wdens(max_at_bf,max_at_bf),
     $     densx(max_at_bf,max_at_bf), densy(max_at_bf,max_at_bf),
     $     densz(max_at_bf,max_at_bf) 

      logical oskel   ! symmetry?

C-------------------------local variables--------------------------

      integer ijatom, next, iat1, iat2, iat3, ish1, ish2,
     $     iab1f, iab1l, iab2f, iab2l, iac1f, iac1l, iac2f, iac2l,
     $     if1, il1, if2, il2,
     $     icart, ic, ip1, ip2
 
      integer i, j

      double precision crd1, crd2 ! atomic coordinates
      dimension crd1(3), crd2(3)

      integer idatom
      dimension idatom(2)

      double precision dE, dx, dy, dz, qfac, fact, q1, q2

      logical status, pointforce

      character*16 name
      integer nxtask, task_size
      external nxtask
cgk debug
      character*2 deriv(3)
      data deriv/'dx','dy','dz'/
cgk end

      task_size = 1
      status = rtdb_parallel(.true.) ! Broadcast reads to all processes

      pointforce = geom_include_bqbq(geom)

      call hf_print_set(1)

      do j=1, nat
       do i=1,3
        frc_sox(i,j)=0d0
        frc_soy(i,j)=0d0
        frc_soz(i,j)=0d0
        frc_sow(i,j)=0d0
       enddo
      enddo

cgk debug
*     write(*,*)'gk: printing so blocks'
cgk end
      ijatom = -1
      next = nxtask(nproc,task_size)
      do iat1=1, nat
        do iat2 = 1, nat
          ijatom = ijatom + 1
          if ( ijatom .eq. next ) then
            status = bas_ce2bfr(basis,iat1,iab1f,iab1l)
            status = bas_ce2bfr(basis,iat2,iab2f,iab2l)
            if (iab1f.le.0 .or. iab2f.le.0) goto 999
            status = bas_ce2cnr(basis,iat1,iac1f,iac1l)
            status = bas_ce2cnr(basis,iat2,iac2f,iac2l)
            call ga_get (g_dens(1), 
     &           iab1f,iab1l,iab2f,iab2l,densx,max_at_bf) 
            call ga_get (g_dens(2), 
     &           iab1f,iab1l,iab2f,iab2l,densy,max_at_bf)
            call ga_get (g_dens(3), 
     &           iab1f,iab1l,iab2f,iab2l,densz,max_at_bf)
            do ish1 = iac1f, iac1l
              do ish2 = iac2f, iac2l
                status = bas_cn2bfr(basis,ish1,if1,il1)
                status = bas_cn2bfr(basis,ish2,if2,il2)
cgk debug
*      write(*,'(a,2i4)')'gk: int block for shells, ', ish1, ish2
cgk end
                call intd_1eso(basis,ish1,basis,ish2,lscr,scr,
     &               lbuf,H)
                do iat3 = 1, nat
                   do icart = 1, 3
*                    call print_soblock(H,if1,il1,if2,il2,iat3,icart,1,
*    &                 nat)
                     call trc_soblock(H,if1,il1,if2,il2,iat3,icart,1,
     &                 nat,iab1f,iab2f,max_at_bf,frc_soz,densz)
*                    call print_soblock(H,if1,il1,if2,il2,iat3,icart,2,
*    &                 nat)
                     call trc_soblock(H,if1,il1,if2,il2,iat3,icart,2,
     &                 nat,iab1f,iab2f,max_at_bf,frc_soy,densy)
*                    call print_soblock(H,if1,il1,if2,il2,iat3,icart,3,
*    &                 nat)
                     call trc_soblock(H,if1,il1,if2,il2,iat3,icart,3,
     &                 nat,iab1f,iab2f,max_at_bf,frc_sox,densx)
                   enddo
                enddo
              enddo
            enddo
  999       continue
            next = nxtask(nproc,task_size)
          endif
        enddo
      enddo

*      write(*,'("forces sow",9e13.5)')((frc_sow(i,j),i=1,3),j=1,nat)
*      write(*,'("forces sox",9e13.5)')((frc_sox(i,j),i=1,3),j=1,nat)
*      write(*,'("forces soy",9e13.5)')((frc_soy(i,j),i=1,3),j=1,nat)
*      write(*,'("forces soz",9e13.5)')((frc_soz(i,j),i=1,3),j=1,nat)
      next = nxtask(-nproc, task_size)
      do j=1, nat
       do i=1,3
*       frc_sox(i,j)=0d0
*       frc_soy(i,j)=0d0
*       frc_soz(i,j)=0d0
        frc_sow(i,j)=0d0
       enddo
      enddo
   

      ijatom = -1
      next = nxtask(nproc,task_size)
      do 90, iat1 = 1, nat
        do 80, iat2 = 1, iat1

          ijatom = ijatom + 1
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

            call ga_get(g_wdens,
     &           iab1f,iab1l,iab2f,iab2l,wdens,max_at_bf)
*           call ga_get (g_dens(1), 
*    &           iab1f,iab1l,iab2f,iab2l,densx,max_at_bf) 
*           call ga_get (g_dens(2), 
*    &           iab1f,iab1l,iab2f,iab2l,densy,max_at_bf)
*           call ga_get (g_dens(3), 
*    &           iab1f,iab1l,iab2f,iab2l,densz,max_at_bf)
           do 70, ish1 = iac1f, iac1l
             if ( iat1.eq.iat2 ) iac2l = ish1
             do 60, ish2 = iac2f, iac2l

C               shell block in atomic (D/Dw)-matrix block
                status = bas_cn2bfr(basis,ish1,if1,il1)
                if1 = if1 - iab1f + 1
                il1 = il1 - iab1f + 1
                status = bas_cn2bfr(basis,ish2,if2,il2)
                if2 = if2 - iab2f + 1
                il2 = il2 - iab2f + 1


C               overlap derivatives
                call intd_1eov(basis,ish1,basis,ish2,lscr,scr,
     &               lbuf,H,idatom)

C     Dsow x S

                if ( idatom(1) .ge. 1 ) then
C               idatom(1).ge.0 <=> idatom(2).ge.0 (no check necessary)
                  ic = 1
                  do 28, icart = 1, 3
                    de = 0.D0
                    do 22, ip1 = if1, il1
                      do 20, ip2 = if2, il2
                        dE = dE + wdens(ip1,ip2) * H(ic)
                        ic = ic + 1
 20                   continue
 22                 continue
                    dE = dE * qfac
                    frc_sow(icart,idatom(1)) = frc_sow(icart,idatom(1)) 
     $                                      - dE - dE
                    frc_sow(icart,idatom(2)) = frc_sow(icart,idatom(2)) 
     $                                      + dE + dE
 28               continue
                endif

C     1el. so. derivatives
*               call intd_1eso(basis,ish1,basis,ish2,lscr,scr,
*    &               lbuf,H)
C     Dso x Hso
                
*               ic=1
*               do 150, iat3 = 1, nat
*                  do 140, icart = 1, 3
c     z componet 
cgk debug
*      write(*,*)'gk: z so trace'
*      write(*,2010)'deriv: ',iat3,deriv(icart)
 2010  format(a,i1,x,a2)
cgk end
*                     dE = 0.D0
*                     do 131, ip1 = if1, il1
*                        do 128, ip2 = if2, il2
cgk debug
*      write(*,2000)'i=',ip1+iab1f-1,' j=',ip2+iab2f-1,' den=',
*    & densz(ip1,ip2),' gint=',H(ic)
 2000 format(a,i3,a,i3,a,f20.10,a,f20.10)
cgk end
*                           dE = dE + densz(ip1,ip2)*H(ic)
*                           ic = ic + 1
*128                     continue
*131                  continue
*                     if ( iat1.ne.iat2 .or. ish1.ne.ish2 ) dE = dE + dE
*                     dE = dE * qfac
*                     frc_soz(icart,iat3) = frc_soz(icart,iat3) + dE
c     y componet 
cgk debug
*      write(*,*)'gk: y so trace'
*      write(*,2010)'deriv: ',iat3,deriv(icart)
cgk end
*                     dE = 0.D0
*                     do 230, ip1 = if1, il1
*                        do 231, ip2 = if2, il2
cgk debug
*      write(*,2000)'i=',ip1+iab1f-1,' j=',ip2+iab2f-1,' den=',
*    & densz(ip1,ip2),' gint=',H(ic)
cgk end
*                           dE = dE + densy(ip1,ip2)*H(ic)
*                           ic = ic + 1
*231                     continue
*230                  continue
*                     if ( iat1.ne.iat2 .or. ish1.ne.ish2 ) dE = dE + dE
*                     dE = dE * qfac
*                     frc_soy(icart,iat3) = frc_soy(icart,iat3) + dE
c     x component 
cgk debug
*      write(*,*)'gk: x so trace'
*      write(*,2010)'deriv: ',iat3,deriv(icart)
cgk end
*                     dE = 0.D0
*                     do 250, ip1 = if1, il1
*                        do 251, ip2 = if2, il2
cgk debug
*      write(*,2000)'i=',ip1+iab1f-1,' j=',ip2+iab2f-1,' den=',
*    & densz(ip1,ip2),' gint=',H(ic)
cgk end
*                           dE = dE + densx(ip1,ip2)*H(ic)
*                           ic = ic + 1
*251                     continue
*250                  continue
*                     if ( iat1.ne.iat2 .or. ish1.ne.ish2 ) dE = dE + dE
*                     dE = dE * qfac
*                     frc_sox(icart,iat3) = frc_sox(icart,iat3) + dE
*140               continue
*150            continue
                
 60          continue
 70        continue
          
 1010      continue

          next = nxtask(nproc,task_size)
          endif

 80     continue
 90   continue
      next = nxtask(-nproc,task_size)
       write(*,'("forces sow",9e13.5)')((frc_sow(i,j),i=1,3),j=1,nat)
       write(*,'("forces sox",9e13.5)')((frc_sox(i,j),i=1,3),j=1,nat)
       write(*,'("forces soy",9e13.5)')((frc_soy(i,j),i=1,3),j=1,nat)
       write(*,'("forces soz",9e13.5)')((frc_soz(i,j),i=1,3),j=1,nat)
 
      return
      end
