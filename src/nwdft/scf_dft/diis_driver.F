      Subroutine diis_driver(ispin,ipol,icall,
     &     MXhist,nbf,geom,AO_bas_han,
     &     iga_fock,iga_dens,iga_tmp,iga_svecs,iga_diis,isvals)

C$Id: diis_driver.F,v 1.2 1995-10-17 08:54:59 d3g681 Exp $
      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)

      Logical LResult

      Integer nhist(2), icall(2)
      integer geom,AO_bas_han

      Save lpulay,npulay,ipulay,nhist

#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"

#if defined(DFT_TIMER)
#include "timer.fh"
#endif
c     
c     Use the DIIS scheme to accelerate convergence.
c     
c*******************************************************************************
      if( icall(ispin).lt.0 )then
        LResult = MA_Free_Heap(lpulay)
        return
      end if
      me=ga_nodeid()
      nproc=ga_nnodes()
c
#if defined(DFT_TIMER)
      call ga_sync
      call t_stat(10,-1)
#endif
c
      one = 1.d0
      zero = 0.d0
c
c     Allocate scratch space.
c
      if( icall(1).eq.0 )then
        npulay  = ipol*(MXhist+1)*(MXhist+2)/2
        LResult = MA_Alloc_Get(MT_Dbl,npulay,'Pulay',lpulay,ipulay)
      end if

      LResult = MA_Push_Get(MT_Dbl,MXhist+1,'Fock wghts',lfwght,ifwght)
c
c     Increment call counter.
c
      icall(ispin) = icall(ispin) + 1
c
#if defined(DFT_TIMER)
      call dft_second(tp1)
#endif
c       
c     Construct a temporary global array.
c       
      call gacrea(geom,AO_bas_han,
     &     nbf,nbf,'ga temp2',iga_tmp2,'atom')
c
C     Mat mult Fock * s-1/2
c
      iw=2
      call diis_builds12(isvals,iga_svecs,iga_tmp,iga_tmp2,nbf,iw)
      call ga_dgemm('N','N',nbf,nbf,nbf,
     *     one,iga_fock,iga_tmp,zero,iga_tmp2)
c
C     Mat mult S+1/2*P (transp)
c
      call gacrea(geom,AO_bas_han,
     &     nbf,nbf,'ga temp3',iga_tmp3,'atom')
      iw=3
      call diis_builds12(isvals,iga_svecs,iga_tmp,iga_tmp3,nbf,iw)
      call ga_dgemm('N','T',nbf,nbf,nbf,
     *     one,iga_tmp,iga_dens,zero,iga_tmp3)
c
C     Mat mult (S+*P(transp))(transp)*(F*S-)
c
      call ga_dgemm('N','N',nbf,nbf,nbf,
     *     one,iga_tmp3,iga_tmp2,zero,iga_tmp)
c
c     Destroy temp globals.
c
      call gadest(iga_tmp2)
      call gadest(iga_tmp3)

      LResult = MA_Push_Get(MT_Dbl,nbf,'row',lrow,irow)
      LResult = MA_Push_Get(MT_Dbl,nbf,'column',lcol,icol)

      call diis_refvec(nbf,iga_tmp,dbl_mb(irow),DBL_MB(icol))

      LResult = MA_Pop_Stack(lcol)
      LResult = MA_Pop_Stack(lrow)

      call ga_sync
c
#if defined(DFT_TIMER)
      call dft_second(tp2)
      subsect(1)='  errvec   '
      call tcllct(t_sub(1,1),tp2,tp1)
      call dft_second(tp1)
#endif
c
c     Compute Pulays B matrix.
c
      call diis_core(Dbl_MB(ipulay),MXhist,nhist,
     &     icall,nbf,ispin,ipol,
     &     iga_diis,iga_tmp)
      lenp = (mxhist+1)*(mxhist+2)/2
c
#if defined(DFT_TIMER)
      call dft_second(tp2)
      subsect(2)='  diis   '
      call tcllct(t_sub(1,2),tp2,tp1)
      call dft_second(tp1)
#endif
c
c     Solve DIIS equation for Hamiltonian weights.
c
      call diis_bsolve(Dbl_MB(ipulay+lenp*(ispin-1)),Dbl_MB(ifwght),
     &     MXhist,nhist(ispin),icall(ispin))
c
#if defined(DFT_TIMER)
      call dft_second(tp2)
      subsect(3)='  bsolve   '
      call tcllct(t_sub(1,3),tp2,tp1)
      call dft_second(tp1)
#endif
c
c     Compute an extrapolated Hamiltonian using the weights.
c
      call diis_hamwgt(Dbl_MB(ifwght),MXhist,nhist,nbf,icall,
     &     ispin,ipol,
     &     iga_fock,iga_diis,iga_tmp)
c
#if defined(DFT_TIMER)
      call dft_second(tp2)
      subsect(4)='  hamwght '
      call tcllct(t_sub(1,4),tp2,tp1)
#endif
c
      LResult = MA_Pop_Stack(lfwght)
c
#if defined(DFT_TIMER)
      ltime = 7
      if( me.eq.0 )then
        write(ltime,*)
        write(ltime,9996) '      Analysis of DIIS'
        write(ltime,*)
        write(ltime,9997) ' number of nodes:  ',nproc
        write(ltime,*)
      end if
      subname='DIIS  '
      call t_stat(4,1)
 9996 format(10x,a)
 9997 format(10x,a,i6)
#endif
c
      return
      end
C      SUBROUTINE MATOUT(A,NR,NC,M)                                      
C      IMPLICIT  REAL*8(A-H,O-Z), INTEGER(I-N)                           
CC                                                                       
CC     Rectangular array output.                                         
CC     input:                                                            
CC       A    Array to be written to for006.                             
CC       NR   The row order of A.                                        
CC       NC   The column order of A.                                     
CC       M    The row dimension of A in the calling routine.             
CC                                                                       
C      DIMENSION A(M,NC)                                                 
C      CHARACTER*15 LINE(5)                                              
C      DATA LINE /5*' _____________ '/                                   
C      MAXCOL=0                                                          
C      NSETS=(NC-1)/5+1                                                  
C      DO 100 NS=1,NSETS                                                 
C        MINCOL=MAXCOL                                                   
C        IF (NS.EQ.NSETS) THEN                                           
C          NUMCOL=NC-MINCOL                                              
C        ELSE                                                            
C          NUMCOL=5                                                      
C        END IF                                                          
C        MAXCOL=MINCOL+NUMCOL                                            
C        MINCOL=MINCOL+1                                                 
C        WRITE (6,1000) (I,I=MINCOL,MAXCOL)                              
C        WRITE (6,1010) (LINE(I),I=1,NUMCOL)                             
C        DO 90 I=1,NR                                                    
C          WRITE (6,1020) I,(A(I,J),J=MINCOL,MAXCOL)                     
C   90   CONTINUE                                                        
C  100 CONTINUE                                                          
C      WRITE (6,1030)                                                    
C 1000 FORMAT (/,5X,5(6X,I3,6X))                                         
C 1010 FORMAT (5X,5A15)                                                  
C 1020 FORMAT (1X,I3,1X,5(E14.7,1X))                                     
C 1030 FORMAT (/)                                                        
C      RETURN                                                            
C      END                                                               
C      SUBROUTINE TPRNTD(A,NR)                                           
C      IMPLICIT  REAL*8(A-H,O-Z), INTEGER(I-N)                           
C      CHARACTER*15 LINE(5)                                              
CC                                                                       
CC     This routine writes a triangular packed real array, A,            
CC     of order NR to for006.                                            
CC                                                                       
C      DIMENSION A(NR*(NR+1)/2)                                          
C      DATA LINE /5*' _____________ '/                                   
C      ISYM2(I,J)=I*(I-1)/2+J                                            
C      MAXCOL=0                                                          
C      NSETS=(NR-1)/5+1                                                  
C      DO 100 NS=1,NSETS                                                 
C        MINCOL=MAXCOL                                                   
C        IF (NS.EQ.NSETS) THEN                                           
C          NUMCOL=NR-MINCOL                                              
C        ELSE                                                            
C          NUMCOL=5                                                      
C        END IF                                                          
C        MAXCOL=MINCOL+NUMCOL                                            
C        MINCOL=MINCOL+1                                                 
C        WRITE(6,1000)(I,I=MINCOL,MAXCOL)                                
C        WRITE(6,1010)(LINE(I),I=1,NUMCOL)                               
C        DO 90 I=MINCOL,NR                                               
C          MXCOL = MIN0(MAXCOL,I)                                        
C          WRITE (6,1020) I,(A(ISYM2(I,J)),J=MINCOL,MXCOL)               
C   90   CONTINUE                                                        
C  100 CONTINUE                                                          
C      WRITE (6,1030)                                                    
C 1000 FORMAT (/,5X,5(6X,I3,6X))                                         
C 1010 FORMAT (5X,5A15)                                                  
C 1020 FORMAT (1X,I3,1X,5(D14.7,1X))                                     
C 1030 FORMAT (/)                                                        
C      RETURN                                                            
C      END                                                               
