C> \ingroup selci
C> @{
      subroutine selci_putinb(iflham,ibase,numf,hbuf,ninhb,jbuf,ninjb)
*
* $Id$
*
#include "implicit.fh"
#include "cbuff.fh"
      dimension hbuf(numf,*),jbuf(2,*)
c
c     put hbuf into the output buffer
c     the output buffer is divided into real and integer*2 parts
c     the integer part is (I,nj,J1,J2,...,Jnj),(I',...),...
c     the real part is (hij1,hij2,...,hijnj),(hi'j1'...)...
c
c     if I is input as -1 then set end of file flag
c
c     in the recored if I is iendr then end of record
c                    if I is iendf then end of last record and the file
c
      if (ibase.eq.-1) then
         ibuff(niinb+1) = iendf
         call selci_swrite(iflham,rbuff,lenbuf)
         nrec = nrec + 1
         nval = nval + nrinb
         return
      endif
c
      do 10 ii = 1,numf
         if (niinb+4.gt.lenibf) then
            ibuff(niinb+1) = iendr
            call selci_swrite(iflham,rbuff,lenbuf)
            nrec = nrec + 1
            nval = nval + nrinb
            niinb = 0
            nrinb = 0
         endif
c
c     put I in the buffer and save space for nj
c
         ibuff(niinb+1) = ibase+ii - imagic
         nj = 0
         iptnj = niinb+2
         niinb = niinb+2
c
         kh = 0
         do 20 jj = 1,ninjb
            jbase = jbuf(1,jj) - imagic
            numf2 = jbuf(2,jj)
            do 30 kk = 1,numf2
               kh = kh+1
               if (abs(hbuf(ii,kh)).gt.1.0d-10) then
                  nj = nj + 1
                  nrinb = nrinb+1
                  niinb = niinb+1
                  rbuff(nrinb) = hbuf(ii,kh)
                  ibuff(niinb) = jbase+kk
                  if (niinb+1.eq.lenibf .or. nrinb.eq.lenrbf) then
                     ibuff(iptnj) = nj
                     ibuff(niinb+1) = iendr
                     call selci_swrite(iflham,rbuff,lenbuf)
                     ibuff(1) = ibase+ii - imagic
                     nrec = nrec + 1
                     nval = nval + nrinb
                     iptnj = 2
                     niinb = 2
                     nrinb = 0
                     nj = 0
                  endif
               endif
 30         continue
 20      continue
         ibuff(iptnj) = nj
 10   continue
c
      ninhb = 0
      ninjb = 0
c
      end
C> @}
