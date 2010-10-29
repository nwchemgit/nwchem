c
c $Id$
c

      subroutine smd_leapf_shake(natms,
     >                           ntcons,
     >                           tstep,
     >                           ekin,
     >                           mass,
     >                           icon1,
     >                           icon2,
     >                           consdist,
     >                           ncc,
     >                           nvv,
     >                           dcc,
     >                           nrij,
     >                           orij,
     >                           fff,
     >                           vvv,
     >                           ccc)

      implicit none
c
      integer natms
      integer ntcons
      double precision tstep
      double precision ekin
      double precision mass(natms)
      integer          icon1(ntcons)
      integer          icon2(ntcons)
      double precision consdist(ntcons)
      double precision ncc(3,natms),nvv(3,natms),dcc(3,natms)
      double precision nrij(3,ntcons),orij(3,ntcons)
      double precision fff(3,natms)
      double precision vvv(3,natms)
      double precision ccc(3,natms)

      integer i,j,it,maxit,iatm1,iatm2

      double precision   force,rma,rmb
      double precision   tmpvx,tmpvy,tmpvz
      double precision   tstepsq,rtstep,tol,mxdiff
      double precision   nrijsq,dijsq,diffsq,dotprod
      character*(42) pname

      pname = "smd_leapf_shake"
      write(*,*) "in "//pname

      ekin=0.0
      tstepsq=tstep**2
      rtstep=1.0/tstep
      tol=1.0d-8
      maxit=100
      mxdiff = 0.0d0


      do i=1,ntcons

       iatm1=icon1(i)
       iatm2=icon2(i)

       orij(1,i)=ccc(1,iatm1)-ccc(1,iatm2)
       orij(2,i)=ccc(2,iatm1)-ccc(2,iatm2)
       orij(3,i)=ccc(3,iatm1)-ccc(3,iatm2)

      enddo

      call smd_lat_rebox(ntcons,orij)

      do i=1,natms

       ncc(1,i)=ccc(1,i)
       ncc(2,i)=ccc(2,i)
       ncc(3,i)=ccc(3,i)

       nvv(1,i)=vvv(1,i)+fff(1,i)*tstep/mass(i)
       nvv(2,i)=vvv(2,i)+fff(2,i)*tstep/mass(i)
       nvv(3,i)=vvv(3,i)+fff(3,i)*tstep/mass(i)

       ccc(1,i)=ncc(1,i)+tstep*nvv(1,i)
       ccc(2,i)=ncc(2,i)+tstep*nvv(2,i)
       ccc(3,i)=ncc(3,i)+tstep*nvv(3,i)

      enddo 


      do i=1,maxit

       do j=1,ntcons

        iatm1=icon1(j)
        iatm2=icon2(j)

        nrij(1,j)=ccc(1,iatm1)-ccc(1,iatm2)
        nrij(2,j)=ccc(2,iatm1)-ccc(2,iatm2)
        nrij(3,j)=ccc(3,iatm1)-ccc(3,iatm2)

       enddo

       call smd_lat_rebox(ntcons,nrij) 

       do j=1,natms
        dcc(1,j)=0.0
        dcc(2,j)=0.0
        dcc(3,j)=0.0
       enddo

       do j=1,ntcons

        iatm1=icon1(j)
        iatm2=icon2(j)

        nrijsq=nrij(1,j)**2+nrij(2,j)**2+nrij(3,j)**2
        dijsq=consdist(j)**2
        diffsq=dijsq-nrijsq
        mxdiff=max(mxdiff,abs(diffsq)/consdist(j))

        dotprod=orij(1,j)*nrij(1,j)
     $         +orij(2,j)*nrij(2,j)
     $         +orij(3,j)*nrij(3,j)

        rma= tstepsq/mass(iatm1)
        rmb=-tstepsq/mass(iatm2)
        force=diffsq/(-2.0*(rma-rmb)*dotprod)

        dcc(1,iatm1)=dcc(1,iatm1)-rma*orij(1,j)*force
        dcc(2,iatm1)=dcc(2,iatm1)-rma*orij(2,j)*force
        dcc(3,iatm1)=dcc(3,iatm1)-rma*orij(3,j)*force
        dcc(1,iatm2)=dcc(1,iatm2)-rmb*orij(1,j)*force
        dcc(2,iatm2)=dcc(2,iatm2)-rmb*orij(2,j)*force
        dcc(3,iatm2)=dcc(3,iatm2)-rmb*orij(3,j)*force


       enddo

       do j=1,ntcons

        iatm1=icon1(j)
        iatm2=icon2(j)

        ccc(1,iatm1)=ccc(1,iatm1)+0.5*dcc(1,iatm1)
        ccc(2,iatm1)=ccc(2,iatm1)+0.5*dcc(2,iatm1)
        ccc(3,iatm1)=ccc(3,iatm1)+0.5*dcc(3,iatm1)
        ccc(1,iatm2)=ccc(1,iatm2)+0.5*dcc(1,iatm2)
        ccc(2,iatm2)=ccc(2,iatm2)+0.5*dcc(2,iatm2)
        ccc(3,iatm2)=ccc(3,iatm2)+0.5*dcc(3,iatm2)

       enddo

       mxdiff=mxdiff*0.5

       if(mxdiff.lt.tol)goto 100

      enddo

100   continue

      do i=1,natms

       nvv(1,i)=(ccc(1,i)-ncc(1,i))*rtstep
       nvv(2,i)=(ccc(2,i)-ncc(2,i))*rtstep
       nvv(3,i)=(ccc(3,i)-ncc(3,i))*rtstep

       tmpvx=0.5*(nvv(1,i)+vvv(1,i))
       tmpvy=0.5*(nvv(2,i)+vvv(2,i))
       tmpvz=0.5*(nvv(3,i)+vvv(3,i))

       ekin=ekin+mass(i)*(tmpvx**2+tmpvy**2+tmpvz**2)

       fff(1,i)=(nvv(1,i)-vvv(1,i))*mass(i)*rtstep
       fff(2,i)=(nvv(2,i)-vvv(2,i))*mass(i)*rtstep
       fff(3,i)=(nvv(3,i)-vvv(3,i))*mass(i)*rtstep

       vvv(1,i)=nvv(1,i)
       vvv(2,i)=nvv(2,i)
       vvv(3,i)=nvv(3,i)

      enddo


      ekin=0.5*ekin

      write(*,*) "out "//pname
    
      return

      END
