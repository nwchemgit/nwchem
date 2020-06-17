c
c     Compute <S**2> expectation value
c
      subroutine tddft_s2(tda,nfc,nocc,nmo,nfv,nov,
     &            x1,x2,y1,y2,corr,s2)
c
      implicit none        
c
      logical tda
      integer nfc(2),nocc(2),nmo(2),nfv(2),nov(2)
      double precision x1(*),x2(*),y1(*),y2(*),corr(*)
      double precision s2
c
c     local variables
c
      double precision na,nb
      integer ja1,ka1,ja2,ka2,jb2,kb2,jb1,jb3,ja3,kb1
      integer la1,la2,lb1,lb2
      double precision s2tol
      double precision rc1,rc2
      double precision xa1,xa2,ya1,ya2
      double precision xb1,xb2,yb1,yb2
#include "global.fh"
c
c     Preliminaries
c
      na=dble(nocc(1)-nfc(1))
      nb=dble(nocc(2)-nfc(2))
      s2tol = 1.0d-05
c
      s2=0d0
c
c     == alpha excitation ket ==
c
      la1=0
cold      do ja1=nfc(1)+1,nocc(1)
      do ja1=ga_nodeid()+nfc(1)+1,nocc(1),ga_nnodes()
       do ka1=nocc(1)+1,nmo(1)-nfv(1)
cold          la1=la1+1
          la1=
     I     (ja1-nfc(1)-1)*(nmo(1)-nfv(1)-nocc(1))+ka1-nocc(1)
          xa1= x1(la1)  !dbl_mb(k_x(1)+la1-1)
          ya1=0.0d0
          if (.not.tda) ya1 = y1(la1) !dbl_mb(k_y(1)+la1-1)
          if (dabs(xa1+ya1).gt.s2tol) then
c a alpha and j beta exchange
            do jb1=nfc(2)+1,nocc(2)
c alpha excitation bra
               la2=0
               do ja2=nfc(1)+1,nocc(1)
                  if (ja1.eq.ja2) then
                  do ka2=nocc(1)+1,nmo(1)-nfv(1)
cold                     la2=la2+1
                     la2=
     I               (ja2-nfc(1)-1)*(nmo(1)-nfv(1)-nocc(1))+ka2-nocc(1)
                     
                        xa2= x1(la2) !dbl_mb(k_x(1)+la2-1)
                        ya2=0.0d0
                        if (.not.tda) ya2= y1(la2) !dbl_mb(k_y(1)+la2-1)
                        rc1=corr((jb1-1)*nmo(1)+ka2)     !dbl_mb(k_corr+(jb1-1)*nmo(1)+ka2-1)
                        rc2=corr((jb1-1)*nmo(1)+ka1)     !dbl_mb(k_corr+(jb1-1)*nmo(1)+ka1-1)
                        s2=s2-(xa2-ya2)*(xa1+ya1)*rc1*rc2
                  enddo
                  endif
               enddo
c beta excitation bra
               lb2=0
               do jb2=nfc(2)+1,nocc(2)
                  if (jb1.eq.jb2) then
                  do kb2=nocc(2)+1,nmo(2)-nfv(2)
cold                    lb2=lb2+1
                     lb2=
     I              (jb2-nfc(2)-1)*(nmo(2)-nfv(2)-nocc(2))+kb2-nocc(2)
                        xb2=x2(lb2)  !dbl_mb(k_x(2)+lb2-1)
                        yb2=0.0d0
                        if (.not.tda) yb2=y2(lb2) !dbl_mb(k_y(2)+lb2-1)
                        rc1=corr((jb2-1)*nmo(1)+ja1)  !dbl_mb(k_corr+(jb2-1)*nmo(1)+ja1-1)
                        rc2=corr((kb2-1)*nmo(1)+ka1)  !dbl_mb(k_corr+(kb2-1)*nmo(1)+ka1-1)
                        s2=s2-(xb2-yb2)*(xa1+ya1)*rc1*rc2
                  enddo
                  endif
               enddo
            enddo
c i alpha and j beta exchange
            do ja3=nfc(1)+1,nocc(1)
               if (ja3.ne.ja1) then
                 do jb1=nfc(2)+1,nocc(2)
c same alpha excitation bra
                    rc1=corr((jb1-1)*nmo(1)+ja3) !dbl_mb(k_corr+(jb1-1)*nmo(1)+ja3-1)
                    s2=s2-(xa1-ya1)*(xa1+ya1)*rc1*rc1
                  enddo
               endif
            enddo
          endif
       enddo
      enddo
c
c     == beta excitation ket ==
c
      lb1=0
cold      do jb1=nfc(2)+1,nocc(2)
      do jb1=ga_nodeid()+nfc(2)+1,nocc(2),ga_nnodes()
         do kb1=nocc(2)+1,nmo(2)-nfv(2)
cold          lb1=lb1+1
            lb1=
     I         (jb1-nfc(2)-1)*(nmo(2)-nfv(2)-nocc(2))+kb1-nocc(2)
              xb1=x2(lb1) !dbl_mb(k_x(2)+lb1-1)
              yb1=0.0d0
              if (.not.tda) yb1=y2(lb1) !dbl_mb(k_y(2)+lb1-1)
              if (dabs(xb1+yb1).gt.s2tol) then
c i alpha and b beta exchange
                do ja1=nfc(1)+1,nocc(1)
c beta excitation bra
                  lb2=0
                  do jb2=nfc(2)+1,nocc(2)
                    if (jb1.eq.jb2) then
                    do kb2=nocc(2)+1,nmo(2)-nfv(2)
cold                       lb2=lb2+1
                       lb2=
     I            (jb2-nfc(2)-1)*(nmo(2)-nfv(2)-nocc(2))+kb2-nocc(2)
                        xb2=x2(lb2) !dbl_mb(k_x(2)+lb2-1)
                        yb2=0.0d0
                        if (.not.tda) yb2=y2(lb2) !dbl_mb(k_y(2)+lb2-1)
                        rc1=corr((kb1-1)*nmo(1)+ja1)  !dbl_mb(k_corr+(kb1-1)*nmo(1)+ja1-1)
                        rc2=corr((kb2-1)*nmo(1)+ja1)  !dbl_mb(k_corr+(kb2-1)*nmo(1)+ja1-1)
                        s2=s2-(xb2-yb2)*(xb1+yb1)*rc1*rc2
                    enddo
                    endif
                  enddo
c alpha excitation bra
                  la2=0
                  do ja2=nfc(1)+1,nocc(1)
                    if (ja1.eq.ja2) then
                    do ka2=nocc(1)+1,nmo(1)-nfv(1)
cold                      la2=la2+1
                       la2=
     I                (ja2-nfc(1)-1)*(nmo(1)-nfv(1)-nocc(1))+ka2-nocc(1)
                        xa2=x1(la2) !dbl_mb(k_x(1)+la2-1)
                        ya2=0.0d0
                        if (.not.tda) ya2=y1(la2)!dbl_mb(k_y(1)+la2-1)
                        rc1=corr((jb1-1)*nmo(1)+ja2)   !dbl_mb(k_corr+(jb1-1)*nmo(1)+ja2-1)
                        rc2=corr((kb1-1)*nmo(1)+ka2)   !dbl_mb(k_corr+(kb1-1)*nmo(1)+ka2-1)
                        s2=s2-(xa2-ya2)*(xb1+yb1)*rc1*rc2
                    enddo
                    endif
                  enddo
                enddo
c i alpha and j beta exchange
                do jb3=nfc(2)+1,nocc(2)
                  if (jb3.ne.jb1) then
                    do ja1=nfc(1)+1,nocc(1)
c same beta excitation bra
                      rc1=corr((jb3-1)*nmo(1)+ja1)     !dbl_mb(k_corr+(jb3-1)*nmo(1)+ja1-1)
                      s2=s2-(xb1-yb1)*(xb1+yb1)*rc1*rc1
                    enddo
                  endif
                enddo
              endif
         enddo
      enddo
      call ga_dgop(2018,s2,1,'+')
      
      s2=s2+na-(na-nb)/2.0d0+((na-nb)/2.0d0)**2
c
      return
      end
