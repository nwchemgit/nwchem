      subroutine kbppv3(psp_filename,formatted_filename,
     >                  nfft1,nfft2,nfft3,unita)
*
* $Id$
*
      implicit none

      character*50     psp_filename,formatted_filename
      integer          nfft1,nfft2,nfft3
      double precision unita(3,3)



      integer nrmax     
      double precision vp(NRMAX,0:3)
      double precision wp(NRMAX,0:3)
      double precision rho(nrmax)

      double precision vl(nfft1/2+1,nfft2,nfft3)
      double precision vnl(nfft1/2+1,nfft2,nfft3,9)
      double precision G(nfft1/2+1,nfft2,nfft3,3)




      
      double precision vnlnrm(9)
      double precision unitg(3,3)
      double precision rc(0:4)
  

      character*2 atom

  
*     pseudopotential data
      OPEN(UNIT=11,FILE=psp_filename,STATUS='OLD',FORM='FORMATTED')
      READ(11,'(A2)',ERR=9116,END=9117) ATOM
      READ(11,*,ERR=9116,END=9117) ZV,AMASS,LMAX
      READ(11,*,ERR=9116,END=9117) (RC(I),I=0,LMAX)
      READ(11,*,ERR=9116,END=9117) NRHO,DRHO
      IF(NRHO.GT.NRMAX) THEN
        IERR=3
        GO TO 9999
      ENDIF
      READ(11,*,ERR=9116,END=9117)(RHO(I),(VP(I,J),J=0,LMAX),I=1,NRHO)
      READ(11,*,ERR=9116,END=9117)(RHO(I),(WP(I,J),J=0,LMAX),I=1,NRHO)
      CLOSE(11)

      WRITE(*,'(A,I1)') 'Lmax = ',LMAX
      WRITE(*,'(a,$)') "enter the highest L desired =>"
      READ(*,*) LMAX0
      LMAX=MIN(LMAX,LMAX0)

*     preparation of constants
      CALL SETUP(IERR)
      IF(IERR.NE.0) THEN
        IERR=IERR+100
        GO TO 9999
      ENDIF

      CALL KBPP(IERR)
      IF(IERR.NE.0) THEN
        IERR=IERR+200
        GO TO 9999
      ENDIF

      l = indx(formatted_filename,' ') - 1
      call openfile(6,formatted_filename,'w',l)     
         call iwrite(6,nfft1,1)
         call iwrite(6,nfft2,1)
         call iwrite(6,nfft3,1)
         call dwrite(6,unita,9)
         call cwrite(6,atom,2)
         call dwrite(6,amass,1)
         call dwrite(6,zv,1)
         call iwrite(6,lmax,1)
         call dwrite(6,rc,lmax+1)
         call dwrite(6,vnlnrm,lmmax)
         call dwrite(6,vl,(nfft1/2+1)*nfft2*nfft3)
         call dwrite(6,vnl,(nfft1/2+1)*nfft2*nfft3*lmmax)
      call closefile(6)

      
      IERR=0
      GO TO 9999

 9110 IERR=10
      GO TO 9999
 9111 IERR=11
      GO TO 9999
 9116 IERR=16
      GO TO 9999
 9117 IERR=17
      GO TO 9999

 9999 IF(IERR.EQ.0) THEN
        WRITE(6,*) ' JOB HAS BEEN COMPLETED.  CODE=',IERR
      ELSE
        WRITE(6,*) ' JOB HAS BEEN TERMINATED DUE TO CODE=',IERR
      ENDIF
      STOP
      END


