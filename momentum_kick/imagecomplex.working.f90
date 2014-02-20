program trycompleximage
  implicit none
  integer, parameter :: neq = 201  !1001  ! neq must be odd
  integer, parameter :: timepoint = 101
  integer, parameter :: mdim = neq !neq/2 + 100
  double complex wfs(timepoint,mdim) 
  double complex :: w(timepoint, 91)
  integer n
  real,parameter :: pi=3.1415926

  open(160,file="x_complex_image.bin",form='unformatted')
  read(160) wfs
  close(160)
  w = wfs(1:timepoint, 30:120)
  
  !call compleximage(timepoint,91,w,'momentum_kick_k.ppm')
  call compleximage(timepoint,mdim,wfs,'momentum_kick_x.ppm')
  
  stop
end


!***************************************************
subroutine compleximage(ndim,mdim,wf,fname)
implicit none	
real*8,parameter :: pi=3.1415926
integer ndim,mdim
integer m,n,i
double complex wf(ndim,mdim)
real*8 rmax,rmin,rmiddle, rr,deltacolor(3),arg
integer  colormax(3),colormin(3),pixel(ndim,mdim,3),dc(3),pixh(3)
      !**the 3 ##s in the byte array stand for r,g,b
character*20 fname !**filename; set to trypix.ppm



rr=maxval(abs(wf))**2 *1.3
rmin=minval(abs(wf))**2.
rmiddle=rmin+(rr-rmin)/3.
rmax=rr !make simpler! maybe,make integer


!**two-color mix. To be changed for a more complex one later****
colormax =(/255,255,255/)   !**gold? yellow?**
colormin =(/0,0,0/)

do n=1,ndim
 do m=1,mdim
!***Michael's color coding - strong complementary colors *****
 	 i=int(mod(arg(wf(n,m))+0.*pi/3.,2.*pi)*6.*255./2./pi) !integer ariphmetics is faster
 	 !write(30,*) n, m, wf(n,m)!mod(arg(wf(n,m))+0.*pi/3.,2.*pi) 
	 if(i<255)then
	   pixh=(/255,i,0/)
	  else if (i<2*255) then
	   pixh=(/510-i,255,0/)
	  else if (i<3*255) then
	   pixh=(/0,255,i-510/)
	  else if (i<1020) then
	   pixh=(/0,1020-i,255/)
   	  else if (i<5*255) then
	   pixh=(/i-1020,0,255/)
	  else 
	   pixh=(/255,0,6*255-i/)
	 endif

!!**smoother cos-like coding; strong primary colors***
!  pixh=0
! rr = arg(wf(n,m))
!    if(rr<(2.*pi/3.)) then
!	  pixh(1)=int(255.*cos(rr*0.75)**2.)
!	 else
!	  pixh(3)=int(255.*cos((rr-4.*pi/3.)*0.75)**2. )
!      endif
!	if(rr<(4.*pi/3.)) then
!	  pixh(2)=int(255.*cos((rr-2.*pi/3.)*0.75)**2.)
!	 else
!	  pixh(1)=int(255.*cos((rr-2.*pi)*0.75)**2.)
!	  endif

   pixel(n,m,:)=pixh
    rr=abs(wf(n,m))**2.
	
!*******now abs(wf) shifts the color to colormin or colormax
	 if(rr<rmiddle)then	!**** make integer later to make faster
	   deltacolor=float(pixh-colormin)
       dc=int(deltacolor*(rr-rmiddle)/(rmiddle-rmin))
       pixel(n,m,:)=pixh+dc
	  else
   	   deltacolor=float(colormax-pixh)
       dc=int(deltacolor*(rr-rmiddle)/(rmax-rmiddle))
       pixel(n,m,:)=pixh+dc
	  endif
  enddo
 enddo

 write(*,*) 'pixs set'
!*****two-color mix finished***************************

write(*,*) 'array done'

!**now writing to a file***
!**header***
 open(11,file=fname)
 write(11,1002)'P3'	!*'P6' for binary output coding later. 'P3' for integer output
 write(11,*) ndim,mdim
 write(11,1001)255  !max color ->1byte
 write(11,*)  '#  CREATED BY ZHENIA THE COLOR MASTER'
 write(11,*)  '#  TOUGHT BY MICHAEL THE KING OF CODES'  !sign of comment
 close(11)

! open(11,file=fname,form='unformatted',ACCESS= 'APPEND')   !*** for binary coding
 open(11,file=fname,ACCESS= 'APPEND')  !**for integer coding
! write(11) pixel

!**integer:
 do m=1,mdim
     write(11,1001)  ((pixel(n,m,i),i=1,3),n=1,ndim)
  enddo
 close(11)

1001    FORMAT ( 2001 (I3,1X) ) 
1002    FORMAT ( (A2) ) 
1003    FORMAT ( B3,\ )

return
end

!************************************************
real function arg(z)
complex z
real,parameter :: pi=3.1415926
real XX,YY,FI

XX=real(z)
YY=imag(z)

FI=ABS(ATAN(YY/XX))
IF((XX.LT.0.).AND.(YY.GT.0.)) FI=Pi-FI
IF((XX.LT.0.).AND.(YY.LT.0.)) FI=Pi+FI
IF((XX.GT.0.).AND.(YY.LT.0.)) FI=2*Pi-FI 

arg=FI

return
end
