      program testnat
* $Id$
      implicit none
      integer iat, jat, kat, lat, nat
      integer ibk, jbk, kbk, lbk
      integer iathi, iatlo
      integer jathi, jatlo
      integer kathi, katlo
      integer lathi, latlo
      integer ijunk, icount
c
      integer i,j,k,l
      integer isym2,isym4
      ISYM2(I,J)=MAX(I,J)*(MAX(I,J)-1)/2+MIN(I,J)
      ISYM4(I,J,K,L)=MAX(ISYM2(I,J),ISYM2(K,L))*
     &    (MAX(ISYM2(I,J),ISYM2(K,L))-1)/2+
     &    MIN(ISYM2(I,J),ISYM2(K,L))
c
      nat = 4
      iat = 0
      jat = 0
      kat = 0
      lat = 0
c
      icount = 0
      do kat = nat, 1, -1
        do iat = nat, kat, -1
          do jat = iat, 1, -1
            lathi = kat
            if (iat.eq.kat) lathi = jat
            do lat = 1,lathi
              icount = icount + 1
              ijunk = isym4(iat,jat,kat,lat)
              write(81,10000)ijunk,iat,jat,kat,lat
*              write(6,10000)ijunk,iat,jat,kat,lat
*              write(81,10000)kat,iat
            enddo
          enddo
        enddo
      enddo
      write(6,*)' original', icount
      write(6,*)' '
c
      ibk = 2
      jbk = 2
      kbk = 2
      lbk = 5
      iat = 0
      jat = 0
      kat = 0
      lat = 0
c
      icount = 0
      do kathi = nat, 1, -kbk
        katlo = max(1,kathi-kbk+1)
        do iathi = nat, 1, -ibk
          iatlo = max(1,iathi-ibk+1)
          do jathi = nat, 1, -jbk
            jatlo = max(1,jathi-jbk+1)
            do iat = iathi, iatlo, -1
              do jat = jathi,jatlo,-1
                do kat = kathi,katlo,-1
                  if (iat.ge.jat) then
                    if (iat.ge.kat) then
                      lathi = kat
                      if (iat.eq.kat) lathi = jat
                      do lat = 1,lathi
                        ijunk = isym4(iat,jat,kat,lat)
                        write(91,10000)ijunk,iat,jat,kat,lat
*                    write(91,10000)kat,iat
*                        write(6,10000)ijunk,iat,jat,kat,lat
                        icount = icount + 1
                      enddo
                    endif
                  endif
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      write(6,*)' blocked ', icount
c
      call exit
10000 format(i7, i5, i5, i5, i5)
      end

      
      
