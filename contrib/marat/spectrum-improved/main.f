       program main
       implicit none
       logical ofile
       integer nmax
       parameter(nmax=20)
       double precision lmax(nmax),sigma(nmax),sum(nmax)
       double precision f(nmax),e(nmax)
       character*(180) buffer
       integer i,n
       double precision l1,l2,tmp
c
       inquire(file="spectrum.dat",exist=ofile)
       if(ofile) then
          open(10,file="spectrum.dat",
     $            form='formatted',status='old')
       else
         write(6,*) "Please provide spectrum.dat file in the format"
         return
       end if 
       n = 0
10     continue
       buffer = " "
       read(10,'(A180)',ERR=30,END=30) buffer
       if(buffer.ne." ".and.buffer(1:1).ne."#") then
         n = n+1
         if(n.gt.nmax) then
           write(6,*) "increase nmax"
           return
         end if
         read(buffer,*) lmax(n),f(n),sigma(n)
       end if
       goto 10
30     close(10)
       l1=lmax(1)
       l2=lmax(1)
       do i=1,n
         write(*,1000) lmax(i),f(i),sigma(i)
         call integrate(lmax(i),sigma(i),sum(i))
         e(i)=f(i)/(0.043*sum(i))
         write(*,1001) e(i)
         tmp = lmax(i)-5*sigma(i)
         if(tmp.lt.l1) l1=tmp
         tmp = lmax(i)+5*sigma(i)
         if(tmp.gt.l2) l2=tmp
       end do
       call plot(n,l1,l2,f,lmax,sigma,e)
1000   FORMAT("Root :", 3G12.6)
1001   FORMAT("Max. extinction :", G12.6)
       end
c
       subroutine integrate(lmax,sigma,sum)
       implicit none
       double precision lmax,sigma,sum
c       
       integer i,n
       double precision f,tmp
       double precision x,dx
       double precision dl

       dx = 0.1d0
       dl=5.0d0*sigma
       n = int(2*dl/dx)
       x = lmax-dl
       f = exp(-((lmax-x)/sigma)**2)/x**2
       sum = 0.5d0*f
       do i=1,n-1
         x=x+dx
         f = exp(-((lmax-x)/sigma)**2)/x**2
         sum = sum + f
       end do
       x=x+dx
       f = exp(-((lmax-x)/sigma)**2)/x**2
       sum = sum +0.5d0*f
       sum = sum*2*dl/dble(n)
       return
       
       end
c
       subroutine plot(n,l1,l2,f,lmax,sigma,e)
       implicit none
       integer n
       double precision l1,l2
       double precision f(n),lmax(n),sigma(n),e(n)
c       
       integer i
       integer k,kmax
       double precision s,tmp
       double precision x,dx
       double precision dl

       open(11,file="spectrum.txt",
     $            form='formatted',status='unknown')
c
       do i=1,n
         write(11,1000) lmax(i),f(i),sigma(i),e(i)
       end do
1000   FORMAT("# ",4G12.6)
c       
       dx = 2.0d0
       kmax = int((l2-l1)/dx)
       x = l1
       do k=1,kmax-1
         x=x+dx
         s = 0.0d0
         do i=1,n
           s =s+e(i)*exp(-((lmax(i)-x)/sigma(i))**2)
         end do
         write(11,*) x,s
       end do
       close(11)
       return
       
       end


c $Id$
