      Subroutine dftg_cdfit( rtdb, geom, 
     &     AO_bas_han,CD_bas_han,
     &     nbf_ao,nbf_cd,nat,tol2e,
     $     scr, lscr, buf, lbuf,
     &     vec, PPP,max_sh_bf,
     &     iga_dens,force, CD_coef )
      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)
#include <mafdecls.h>
#include <rtdb.fh>
#include <bas.fh>
#include <global.fh>
#include <tcgmsg.fh>
c****
c**** handles
c****
      Integer rtdb,geom,AO_bas_han,CD_bas_han
      Integer iga_dens
C
C     scratch arrays for calls to integrals package
C
      dimension scr(lscr), buf(lbuf), idatom(4)
C
C     scratch arrays for dgem
C
      dimension vec(max_sh_bf*max_sh_bf),PPP(max_sh_bf,max_sh_bf)

c-->  Charge Density Expansion Coefficients
      
      Dimension CD_coef(nbf_cd)
c**** 
c**** force vector
c**** 
      Real*8 force(3,nat)
      logical doesit,isit
c--   >
      nproc  = ga_nnodes()
      me = ga_nodeid()
C
C     read CD_coeff vecor written by fitcd
C
      if(me.eq.0) then
      inquire(file='cdfit.dat',exist=doesit,opened=isit)
      if(.not.isit)then
        if(.not.doesit)then
          open(unit=79,file='cdfit.dat',access='sequential',
     &         form='unformatted',status='new')
        else
          open(unit=79,file='cdfit.dat',access='sequential',
     &         form='unformatted',status='old')
          endif
        endif
          read(79) (CD_coef(i),i=1,nbf_cd)
      rewind 79
      endif

      call ga_brdcst(24100,CD_coef,mdtob(nbf_cd),0)

c      
c     Determine the characteristics of the AO and CD Gaussian basis sets.
c      
      
      if( .not. bas_numcont(AO_bas_han,nshells_ao) )then
        call errquit('Exiting in grad_cdfit.',2)
      end if
      
      if( .not. bas_numcont(CD_bas_han,nshells_cd) )then
        call errquit('Exiting in grad_cdfit.',3)
      end if
      

C
C     2el. 3-c integral derivatives
C
      next=nxtval(nproc)+1
      do 205 ishp = 1,nshells_ao
        if(ishp.eq.next) then
          next=nxtval(nproc)+1
          if( .not. bas_cn2bfr(AO_bas_han,ishp,ifirstp,ilastp) )
     &         call errquit('Exiting in grad_cdfit.',6)
          nshbfp=ilastp-ifirstp+1
          
          do  ishq = 1,ishp
            if( .not. bas_cn2bfr(AO_bas_han,ishq,ifirstq,ilastq) )
     &           call errquit('Exiting in grad_cdfit.',7)
            nshbfq=ilastq-ifirstq+1
            
            Nao2 = nshbfp*nshbfq

c
c           get shell block of DM
c
            call ga_get (iga_dens, ifirstq, ilastq, ifirstp, ilastp,
     &                   PPP, max_sh_bf)

c
c           Schwarz screening on product MAX(ABS(DM(i,j)))*schwarz(i,j)
c
            PPP_max=0.d0
            do nu=1,nshbfp
              do mu=1,nshbfq
                PPP_max=max(PPP_max,abs(PPP(mu,nu)))
              enddo
            enddo
            if(PPP_max*schwarz_shell(ishp,ishq).gt.tol2e) then
              
              do ishd=1,nshells_cd

                if( .not. bas_cn2bfr(CD_bas_han,ishd,ifirstd,ilastd) )
     &               call errquit('Exiting in grad_cdfit.',8)
                nshbfd=ilastd-ifirstd+1
                Nint = nshbfp*nshbfq*nshbfd

                call intd_2e3c(CD_bas_han, ishd, 
     &               AO_bas_han, ishp, ishq,
     &               lscr, scr, lbuf, buf, idatom)

                icount=1
                do iat=1,4
                  if(idatom(iat).gt.0) then
                    do icart=1,3

                      call dgemv('N',Nao2,nshbfd,1.D0,buf(icount),Nao2,
     &                     CD_coef(ifirstd),1,0.D0,vec,1)
                      icount=icount+Nint
                      
                      dE=0.d0
                      icount2=1
                      do ifp=1,nshbfp
                        dE=dE+ddot(
     *                       nshbfq,PPP(1,ifp),1,vec(icount2),1)
                        icount2 = icount2 + nshbfq
                      enddo
                      if(ishp.ne.ishq) dE=dE+dE
                      
                      force(icart,idatom(iat)) = 
     =                     force(icart,idatom(iat)) + dE

                    enddo

                  else
                    icount = icount+3*Nint
                  endif
                enddo
                

              enddo

            endif
            
          enddo
          
        endif

  205 continue

      next=nxtval(-nproc)

C
C     2el. 2-c integral derivatives
C
      next=nxtval(nproc)+1
      do 305 ishp = 1,nshells_cd
        if(ishp.eq.next) then
          next=nxtval(nproc)+1
          if( .not. bas_cn2bfr(CD_bas_han,ishp,ifirstp,ilastp) )
     &         call errquit('Exiting in grad_cdfit.',9)
          nshbfp=ilastp-ifirstp+1
          
          do  ishq = 1,ishp
            if( .not. bas_cn2bfr(CD_bas_han,ishq,ifirstq,ilastq) )
     &           call errquit('Exiting in grad_cdfit.',10)
            nshbfq=ilastq-ifirstq+1
            Nint = nshbfp*nshbfq
            call intd_2e2c(CD_bas_han, ishp,  CD_bas_han, ishq, 
     &           lscr, scr, lbuf, buf, idatom)

            icount = 1
            do iat=1,2
              if(idatom(iat).gt.0) then
                ist=1
                do icart = 1,3
                  call dgemv('n',nshbfq,nshbfp,
     &                 1.d0,buf(icount),nshbfq,
     &                 CD_coef(ifirstp),1,
     &                 0.d0,vec(ist),1)
                  ist=ist+nshbfq
                  icount=icount+Nint
                enddo
                fac=-1.d0
                if(ishp.eq.ishq) fac=-.5d0
                call dgemv('t',nshbfq,3,
     &               fac,vec,nshbfq,
     &               CD_coef(ifirstq),1,
     &               1.d0,force(1,idatom(iat)),1)
              else
                icount = icount + 3 * Nint
              endif
            enddo
          
          enddo

        endif

  305 continue

      next=nxtval(-nproc)

      return 
      end 
c****
