program tryrealimage
  implicit none
  !integer, parameter :: neq = 1001  !1001  ! neq must be odd
  !integer, parameter :: timepoint = 2401
  !integer, parameter :: mdim = neq !neq !/2 + 100
  integer :: neq
  integer :: timepoint
  real*8, allocatable ::  wfs(:,:) 
  !real*8 :: wfs(timepoint,neq)
  real,parameter :: pi=3.1415926
  character*16 :: outputfile
  character*20 :: inputfile
  integer :: i, j


  write(*,*) "Enter the dimension of the image array:"
  read(*,*) timepoint, neq
  allocate(wfs(timepoint,neq))
    
  write(*,*) "Input file name (bin):"
  read(*,*) inputfile

  write(*,*) "Output file name (ppm):"
  read(*,*) outputfile
  !fname ="flag3.ppm"
  open(160,file=inputfile,form='unformatted')
  read(160) wfs
  close(160)
!  do i=1, timepoint
!    do j=1, mdim
!      write(99,*) i, j, wfs(i,j)
!    end do
!  end do
  write(*,*) "=================="
  call realimage(timepoint,neq,wfs,outputfile)

  stop
end





 !***********************************************

subroutine realimage(ndim,mdim,r1,fname)
implicit none
integer ndim,mdim
integer m,n,i
!real*8 r(ndim,mdim),r1(ndim,mdim),rmax,rmin,rall,rzero,rmi,rma
double precision :: r(ndim,mdim),r1(ndim,mdim),rmax,rmin,rall,rzero,rmi,rma
integer  colormax(3),colormin(3),colorzero(3),colorall(3)
integer pixel(ndim,mdim,3),dc(3)
integer dcminzero(3), dcmaxmin(3),dcallmax(3)
     !**the 3 ##s in the byte array stand for r,g,b
character*16 fname !**filename; set to trypix.ppm


!open(15,file='controls.txt')

r=sqrt(r1)
write(*,*) 'inside the subroutine'

!**two-color mix. To be changed for a more complex one later****
colormax =(/150,150,30/)   !**gold? yellow?**
colormin =(/10,10,255/)
colorall = (/255,255,250/)
colorzero=(/5,0,10/)

dcminzero=float(colormin-colorzero)
dcmaxmin=float(colormax-colormin)
dcallmax=float(colorall-colormax)


rmi=minval(r)
rma=maxval(r)
rzero=rmi+0.99*2*(rma-rmi)/10.
rmin=rmi+0.99*3.5*(rma-rmi)/10.
rmax=rmi+0.99*7.*(rma-rmi)/10.
rall=rmi+0.99*9.5*(rma-rmi)/10.

do n=1,ndim
 do m=1,mdim
  if(r(n,m)<rzero)then  ! simplify later
    pixel(n,m,:)=colorzero
   else
   if(r(n,m)<rmin)then  ! simplify later
     dc=int(dcminzero*(r(n,m)-rzero)/(rmin-rzero))
     pixel(n,m,:)=colorzero+dc
        else
    if(r(n,m)<rmax)then         ! simplify later
      dc=int(dcmaxmin*(r(n,m)-rmin)/(rmax-rmin))
      pixel(n,m,:)=colormin+dc
         else
     if(r(n,m)<rall)then        ! simplify later
       dc=int(dcallmax*(r(n,m)-rmax)/(rall-rmax))
       pixel(n,m,:)=colormax+dc
          else
        pixel(n,m,:)=colorall
          endif
         endif
        endif
       endif
 enddo
 enddo

 write(*,*) 'pixs set'
!*****two-color mix finished***************************

write(*,*) 'array done'

!**now writing to a file***
!**header***
 open(11,file=fname)
 write(11,1002)'P3'     !*'P6' for binary output coding later. 'P3' for integer output
 write(11,*) ndim,mdim
 write(11,1001)255  !max color ->1byte
 write(11,*)  '#'  !sign of comment
 close(11)

 write(*,*) 'header done'
!**contents - binary*****
 open(11,file=fname,ACCESS= 'APPEND')  !**for integer coding
 write(*,*)  'second time opened'

!**integer:
 do m=1,mdim
!  do n=1,ndim
    write(11,1001)  ((pixel(n,m,i),i=1,3),n=1,ndim)
 !  enddo
 enddo
 close(11)

! close(15)

1001    FORMAT ( 3000 (I4,1X) )
1002    FORMAT ( (A2) )
1003    FORMAT ( B3,\ )

return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!SECOND


! -------------------------------------------------------------------------------
!
!subroutine realimageM(r,fname,rmin,rmax)
!!use stack_vars
!use dimensions
!implicit none
!integer mdim
!integer m,n,i
!real*8  r(ndim,ndim),rmax,rmin,deltacolor(3) !,rmaxM !, ,deltacolorM(3)
!integer  colormax(3),colormin(3),pixel(ndim,ndim,3),dc(3) !,colorall(3)
!     !**the 3 ##s in the byte array stand for r,g,b
!character*14 fname !**filename; set to trypix.ppm

! ! write(*,*) '        real image started'

!mdim = ndim

!! - outputt data into a text file -
!!       open(10,file=fname//'RealImage.txt')
!!       do n = 1,ndim
!!       do m = 1,ndim
!!         write(10,*)   n,m, r(n,m)
!!        enddo
!!        enddo
!!       close(10)

!!**two-color mix. To be changed for a more complex one later****
!colormax =(/255,255,255/)   !**gold? yellow?**
!colormin =(/0,0,0/)

!deltacolor=dble(colormax-colormin)

!!rmax=maxval(abs(r))

!do n=1,ndim
! do m=1,mdim
!    if(r(n,m)<rmin)then
!      pixel(n,m,:)=colormin
!     else
!       if(r(n,m)>rmax)then
!         pixel(n,m,:)=colormax
!        else
!         dc=int(deltacolor*(r(n,m)-rmin)/(rmax-rmin))
!         pixel(n,m,:)=colormin+dc
!        endif
!     endif
! enddo
! enddo

!!*****two-color mix finished***************************

!!**now writing to a file***
!!**header***
! open(11,file=fname)
! write(11,1002)'P3'     !*'P6' for binary output coding later. 'P3' for
!integer output
! write(11,1001) ndim,mdim
! write(11,1001)255  !max color ->1byte
! write(11,*)  '#'  !sign of comment
! close(11)


!!**contents - binary*****
!! open(11,file=fname,form='unformatted',ACCESS= 'APPEND')   !*** for binary coding
! open(11,file=fname,ACCESS= 'APPEND')  !**for integer coding

!!**integer:
! do m=1,mdim
!!  do n=1,ndim
!    write(11,1001)  ((pixel(n,m,i),i=1,3),n=1,ndim)
! !  enddo
! enddo
! close(11)

! close(15)

!1001    FORMAT ( 4000 (I4,1X) )
!1002    FORMAT ( (A2) )
!1003    FORMAT ( B3,\ )

!!write(*,*) ' real image finished'

!return
!end
!!
!!
!! ------------------------------------------------------
