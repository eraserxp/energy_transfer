program plots
implicit none
real,parameter :: pi=3.1415926
integer, parameter :: ndim = 101
integer, parameter :: mdim = 201
integer    n,m
double precision       probs(ndim,mdim)
double precision  rr(ndim,mdim),ii(ndim,mdim)
double precision    r(ndim,mdim)
double complex wf(ndim,mdim)

!double precision       probs(ndim,mdim)
!double precision  rr(ndim,mdim),ii(ndim,mdim)
!double precision    r(ndim,mdim)

character*13 fname,fname1


! -----  X - part, 3 versions ---
open(160,file='x_real_image.bin',form='unformatted')
read(160) rr  ! or read(160) potentialC!
close(160)

open(160,file='x_complex_image.bin',form='unformatted')
read(160) ii  ! or read(160) potentialC!
close(160)

!do n=1,ndim
!do m=1,mdim
! probs(n,m)=exp(-dble((n-20)**2)/dble((20)**2))
!enddo 
!enddo 

wf = rr+ (0.,1.)*ii
r=abs(wf)
probs = r**2

fname='X_CPL.ppm'
call compleximageTruncated(ndim,mdim,wf,fname)
fname='X_RP1.ppm'
call realimage(ndim,mdim,r,fname)
fname='X_RP2.ppm'
call realimageM(ndim,mdim,probs,fname,0.d0)

! -----  K - part, 3 versions ---
open(160,file='k_real_image.bin',form='unformatted')
read(160) rr  ! or read(160) potentialC!
close(160)

open(160,file='k_complex_image.bin',form='unformatted')
read(160) ii  ! or read(160) potentialC!
close(160)

!do n=1,ndim
!do m=1,mdim
! probs(n,m)=exp(-dble((n-20)**2)/dble((20)**2))
!enddo 
!enddo 

wf = rr+ (0.,1.)*ii
r=abs(wf)
probs = r**2
fname='K_CPL.ppm'
call compleximageTruncated(ndim,mdim,wf,fname)
fname='K_RP1.ppm'
call realimage(ndim,mdim,r,fname)
fname='K_RP2.ppm'
call realimageM(ndim,mdim,probs,fname,0.d0)

stop
end
! --------------------------------------------------------------
subroutine compleximageTruncated(ndim,mdim,wf,fname)
implicit none	
double precision,parameter :: pi=3.1415926
integer ndim,mdim
integer m,n,i
double complex wf(ndim,mdim)
double precision rmax,rmin,rmiddle, rr,deltacolor(3),arg
integer  colormax(3),colormin(3),pixel(ndim,mdim,3),dc(3),pixh(3)
      !**the 3 ##s in the byte array stand for r,g,b
character*30 fname !**filename; set to trypix.ppm



rr=maxval(abs(wf))**2.
rmin=minval(abs(wf))**2.
rmiddle=rmin+(rr-rmin)/12.
rmax=rr !make simpler! maybe,make integer


!**two-color mix. To be changed for a more complex one later****
colormax =(/255,255,255/)   !**gold? yellow?**
colormin =(/0,0,0/)

do n=1,ndim
 do m=1,mdim
!***Michael's color coding - strong complementary colors *****
 	 i=int(mod(arg(wf(n,m))+0.*pi/3.,2.*pi)*6.*255./2./pi) !integer ariphmetics is faster
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
 write(11,1001) ndim,mdim
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

1001    FORMAT ( 1005 (I4,1X) ) 
1002    FORMAT ( (A2) ) 
1003    FORMAT ( B3,\ )

return
end

!************************************************
double precision function arg(z)
double complex z
double precision,parameter :: pi=3.1415926
double precision XX,YY,FI

XX=real(z)
YY=imag(z)

FI=ABS(ATAN(YY/XX))
IF((XX.LT.0.).AND.(YY.GT.0.)) FI=Pi-FI
IF((XX.LT.0.).AND.(YY.LT.0.)) FI=Pi+FI
IF((XX.GT.0.).AND.(YY.LT.0.)) FI=2*Pi-FI 

arg=FI

return
end
!
! ----------------------------------------------------------------------------------------
!
!
! ---------------------------------------------------------------------
!
subroutine realimage(ndim,mdim,r,fname)
implicit none
integer ndim,mdim
integer m,n,i
double precision r(ndim,mdim),rmax,rmin,rall,rzero,rmi,rma 
integer  colormax(3),colormin(3),colorzero(3),colorall(3)
integer pixel(ndim,mdim,3),dc(3)
integer dcminzero(3), dcmaxmin(3),dcallmax(3)
      !**the 3 ##s in the byte array stand for r,g,b
character*13 fname !**filename; set to trypix.ppm


write(*,*) 'inside the subroutine'

open(15,file='controls.txt')
! goto 200

!**two-color mix. To be changed for a more complex one later****
! - strait -
!colormax =(/220,220,40/)   !**gold? yellow?**
!colormin =(/150,0,180/)
!colorall = (/255,255,255/)
!colorzero=(/30,0,30/)

! - inverted -
colormin =(/220,220,40/)   !**gold? yellow?**
colormax =(/150,0,180/)
colorzero = (/255,255,255/)
colorall  =(/30,0,30/)

dcminzero=float(colormin-colorzero)
dcmaxmin=float(colormax-colormin)
dcallmax=float(colorall-colormax)

r=sqrt(r)

rmi=minval(r)
rma=maxval(r)
rzero=rmi+1.*(rma-rmi)/10.
rmin=rmi+5.*(rma-rmi)/10.
rmax=rmi+8*(rma-rmi)/10.
rall=rmi+10.*(rma-rmi)/10.



do n=1,ndim
 do m=1,mdim
   if(r(n,m)<rzero)then	 ! simplify later
     pixel(n,m,:)=colorzero
    else
    if(r(n,m)<rmin)then	 ! simplify later
      dc=int(dcminzero*(r(n,m)-rzero)/(rmin-rzero))
      pixel(n,m,:)=colorzero+dc
 	 else
     if(r(n,m)<rmax)then	 ! simplify later
       dc=int(dcmaxmin*(r(n,m)-rmin)/(rmax-rmin))
       pixel(n,m,:)=colormin+dc
	  else
      if(r(n,m)<rall)then	 ! simplify later
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
 write(11,1002)'P3'	!*'P6' for binary output coding later. 'P3' for integer output
 write(11,1001) ndim,mdim
 write(11,1001)255  !max color ->1byte
 write(11,*)  '#'  !sign of comment
 close(11)

 write(*,*) 'header done'
!**contents - binary*****
! open(11,file=fname,form='unformatted',ACCESS= 'APPEND')   !*** for binary coding
 open(11,file=fname,ACCESS= 'APPEND')  !**for integer coding
 write(*,*)  'second time opened'
! write(11) pixel

!write(*,*) pixel(1,100,:)

!**integer:
 do m=1,mdim
!  do n=1,ndim
     write(11,1001)  ((pixel(n,m,i),i=1,3),n=1,ndim)
 !  enddo
  enddo
 close(11)

 close(15)

1001    FORMAT ( 2000 (I4,1X) ) 
1002    FORMAT ( (A2) ) 
1003    FORMAT ( B3,\ )

200 continue

return
end
!
! -------------------------------------------------------------------------------
!
subroutine realimageM(ndim,mdim,r,fname,rmin)
!use stack_vars
implicit none
integer ndim,mdim
integer m,n,i
double precision  r(ndim,mdim),rmax,rmin,deltacolor(3) !,rmaxM !, ,deltacolorM(3)
integer  colormax(3),colormin(3),pixel(ndim,mdim,3),dc(3) !,colorall(3)
      !**the 3 ##s in the byte array stand for r,g,b
character*13 fname !**filename; set to trypix.ppm

write(*,*) ' real image started'

!**two-color mix. To be changed for a more complex one later****
colormin =(/240,240,180/)   !**gold? yellow?**
colormax =(/0,0,140/)

deltacolor=dble(colormax-colormin)

!write(*,*) 'real 197'

!rmin = 0.
!rmin form outside =minval(r)*0.d0
write(*,*) 'real 201'
rmax=maxval(abs(r)) !3.

write(*,*) minval(r),maxval(r)
write(*,*) rmin,rmax

do n=1,ndim
 do m=1,mdim
     if(r(n,m)>rmax)then
        pixel(n,m,:) = colormax
      else
       dc=int(deltacolor*(r(n,m)-rmin)/(rmax-rmin))
       pixel(n,m,:)=colormin+dc
      endif
  enddo
 enddo

!*****two-color mix finished***************************



!**now writing to a file***
!**header***
 open(11,file=fname)
 write(11,1002)'P3'	!*'P6' for binary output coding later. 'P3' for integer output
 write(11,1001) ndim,mdim
 write(11,1001)255  !max color ->1byte
 write(11,*)  '#'  !sign of comment
 close(11)


!**contents - binary*****
! open(11,file=fname,form='unformatted',ACCESS= 'APPEND')   !*** for binary coding
 open(11,file=fname,ACCESS= 'APPEND')  !**for integer coding
 write(*,*)  'second time opened'
! write(11) pixel

!write(*,*) pixel(1,100,:)

!write(11) pixel
!**binary:
!do m=1,mdim
! write(11,1003) (pixel(n,m,1),pixel(n,m,2),pixel(n,m,3),n=1,ndim)
! enddo


!**integer:
 do m=1,mdim
!  do n=1,ndim
     write(11,1001)  ((pixel(n,m,i),i=1,3),n=1,ndim)
 !  enddo
  enddo
 close(11)

 close(15)

1001    FORMAT ( 4000 (I4,1X) ) 
1002    FORMAT ( (A2) ) 
1003    FORMAT ( B3,\ )

!write(*,*) ' real image finished'

return
end
!



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
