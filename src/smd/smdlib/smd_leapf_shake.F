c
c $Id: smd_leapf_shake.F,v 1.1 2008/10/01 22:32:30 marat Exp $
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
      double precision ncc(natms,3),nvv(natms,3),dcc(natms,3)
      double precision nrij(ntcons,3),orij(ntcons,3)
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

       orij(i,1)=ccc(1,iatm1)-ccc(1,iatm2)
       orij(i,2)=ccc(2,iatm1)-ccc(2,iatm2)
       orij(i,3)=ccc(3,iatm1)-ccc(3,iatm2)

      enddo

      call smd_lat_rebox(ntcons,orij)

      do i=1,natms

       ncc(i,1)=ccc(1,i)
       ncc(i,2)=ccc(2,i)
       ncc(i,3)=ccc(3,i)

       nvv(i,1)=vvv(1,i)+fff(1,i)*tstep/mass(i)
       nvv(i,2)=vvv(2,i)+fff(2,i)*tstep/mass(i)
       nvv(i,3)=vvv(3,i)+fff(3,i)*tstep/mass(i)

       ccc(1,i)=ncc(i,1)+tstep*nvv(i,1)
       ccc(2,i)=ncc(i,2)+tstep*nvv(i,2)
       ccc(3,i)=ncc(i,3)+tstep*nvv(i,3)

      enddo 


      do i=1,maxit

       do j=1,ntcons

        iatm1=icon1(j)
        iatm2=icon2(j)

        nrij(j,1)=ccc(1,iatm1)-ccc(1,iatm2)
        nrij(j,2)=ccc(2,iatm1)-ccc(2,iatm2)
        nrij(j,3)=ccc(3,iatm1)-ccc(3,iatm2)

       enddo

       call smd_lat_rebox(ntcons,nrij) 

       do j=1,natms
        dcc(j,1)=0.0
        dcc(j,2)=0.0
        dcc(j,3)=0.0
       enddo

       do j=1,ntcons

        iatm1=icon1(j)
        iatm2=icon2(j)

        nrijsq=nrij(j,1)**2+nrij(j,2)**2+nrij(j,3)**2
        dijsq=consdist(j)**2
        diffsq=dijsq-nrijsq
        mxdiff=max(mxdiff,abs(diffsq)/consdist(j))

        dotprod=orij(j,1)*nrij(j,1)
     $         +orij(j,2)*nrij(j,2)
     $         +orij(j,3)*nrij(j,3)

        rma= tstepsq/mass(iatm1)
        rmb=-tstepsq/mass(iatm2)
        force=diffsq/(-2.0*(rma-rmb)*dotprod)

        dcc(iatm1,1)=dcc(iatm1,1)-rma*orij(j,1)*force
        dcc(iatm1,2)=dcc(iatm1,2)-rma*orij(j,2)*force
        dcc(iatm1,3)=dcc(iatm1,3)-rma*orij(j,3)*force
        dcc(iatm2,1)=dcc(iatm2,1)-rmb*orij(j,1)*force
        dcc(iatm2,2)=dcc(iatm2,2)-rmb*orij(j,2)*force
        dcc(iatm2,3)=dcc(iatm2,3)-rmb*orij(j,3)*force


       enddo

       do j=1,ntcons

        iatm1=icon1(j)
        iatm2=icon2(j)

        ccc(1,iatm1)=ccc(1,iatm1)+0.5*dcc(iatm1,1)
        ccc(2,iatm1)=ccc(2,iatm1)+0.5*dcc(iatm1,2)
        ccc(3,iatm1)=ccc(3,iatm1)+0.5*dcc(iatm1,3)
        ccc(1,iatm2)=ccc(1,iatm2)+0.5*dcc(iatm2,1)
        ccc(2,iatm2)=ccc(2,iatm2)+0.5*dcc(iatm2,2)
        ccc(3,iatm2)=ccc(3,iatm2)+0.5*dcc(iatm2,3)

       enddo

       mxdiff=mxdiff*0.5

       if(mxdiff.lt.tol)goto 100

      enddo

100   continue

      do i=1,natms

       nvv(i,1)=(ccc(1,i)-ncc(i,1))*rtstep
       nvv(i,2)=(ccc(2,i)-ncc(i,2))*rtstep
       nvv(i,3)=(ccc(3,i)-ncc(i,3))*rtstep

       tmpvx=0.5*(nvv(i,1)+vvv(1,i))
       tmpvy=0.5*(nvv(i,2)+vvv(2,i))
       tmpvz=0.5*(nvv(i,3)+vvv(3,i))

       ekin=ekin+mass(i)*(tmpvx**2+tmpvy**2+tmpvz**2)

       fff(1,i)=(nvv(i,1)-vvv(1,i))*mass(i)*rtstep
       fff(2,i)=(nvv(i,2)-vvv(2,i))*mass(i)*rtstep
       fff(3,i)=(nvv(i,3)-vvv(3,i))*mass(i)*rtstep

       vvv(1,i)=nvv(i,1)
       vvv(2,i)=nvv(i,2)
       vvv(3,i)=nvv(i,3)

      enddo


      ekin=0.5*ekin

      write(*,*) "out "//pname
    
      return

      END
