!!!!!!calculation overlap area between two atoms from different layer!!!!!
Module ran_mod
	Implicit None
contains
	function reg (x1 , y1 , x2 , y2)
	implicit none
	double precision,parameter :: pi = 3.141592653589793239
	double precision :: x1,y1,x2,y2,r,d,reg
	d=dsqrt((x1−x2)**2+(y1−y2)**2)
	r=0.5d0*1.42d0
	if(d.le.2d0*r)then
		reg=((pi*r**2)*2d0*dacos(d/2d0/r)/(2d0*pi)−dsqrt(r**2−(d/2d0)**2)*(d/2d0))*2d0/(pi*r**2d0)
	else
		reg=0d0
	endif
	end function reg
End Module ran_mod

program main
use ran_mod
implicit none
double precision :: xx,yy,pi = 3.141592653589793239,xmin1,xmin2,xmax1,xmax2
integer::i,j,k,nline,natom,step,frame
double precision,DIMENSION(:,:)::atom(1000000,12),box(3,2)
double precision,DIMENSION(:)::tmp(1000000)
character*200::up,down,filenum,in

do k=1,500
	write(filenum,*) k
	in='../'//trim(adjustl(filenum))//’.dump’
	open(12,file=in,status='unknown')
	write (*,*) in
	read(12,*)
	read(12,*) step
	read(12,*)
	read(12,*) natom
	read(12,*)
	read(12,*) box(1,1),box(1,2)
	read(12,*) box(2,1),box(2,2)
	read(12,*) box(3,1),box(3,2)
	read(12,*)
	do i=1,natom
	read(12,*) atom(i,:)
	enddo
	write(*,*) k,'read done'
	tmp=0d0
	do i =1,natom
		do j =1,natom
			xx=dabs(atom(i,3)−atom(j,3))
			if(dabs(atom(i,4)−atom(j,4)) .gt. dabs(box(2,2)−box(2,1))/2d0)then
				yy=abs(box(2,2)−box(2,1))−abs(atom(i,4)−atom(j,4))
			else
				yy=abs(atom(i,4)−atom(j,4))
			endif
			if(atom(i,2) .ne. atom(j,2))then
				tmp(i)=tmp(i)+reg(0d0,0d0,xx,yy) !sum of overlap area
			endif
		enddo
	enddo
!!!!!!!!!!!!!!!output!!!!!!!!!!!!!!!
	up='up'//trim(adjustl(filenum))//’.dat’
	down='down'//trim(adjustl(filenum))//’.dat’
	xmin1=minval(atom(1:natom/2,3))/10
	xmin2=minval(atom(natom/2+1:natom,3))/10
	xmax1=maxval(atom(1:natom/2,3))/10
	xmax2=maxval(atom(natom/2+1:natom,3))/10
	write (*,*) xmin1,xmin2

!!!!!!!!!!!!!!!LRI along x direction for bottom layer!!!!!!!!!!!!!!!
	open(55,file=down, status='unknown')
	do i=1,natom/4
		if((atom(2*i,3)+atom(2*i−1,3))/20d0 > xmin1 . and . (atom(2*i,3 )+atom(2*i−1,3))/20d0 <xmax2) then
		write(55,'(5F12.5)') (atom(2*i,3 )+atom(2*i−1,3))/20d0−xmin1,(atom(2*i,6)+atom(2*i−1,6))/2d0,(atom(2*i,7)+atom(2*i−1,7))/2d0,(atom(2*i,9)+atom(2*i−1,9))/2d0,(tmp(2*i)+tmp(2*i−1))/2d0 ! ! ! x , stress_xx , stress_yy , stress_xy,lri
		endif
	enddo
	close(55)
	
!!!!!!!!!!!!!!!LRI along x direction for top layer!!!!!!!!!!!!!!!
	open(66,fil e=up, status='unknown')
	do i=natom/4+1,natom/2
		if((atom(2*i,3)+atom(2*i−1,3))/20d0>xmin2 .and. (atom(2*i,3)+atom(2*i−1,3))/20d0 <xmax1)then
			write(66,'(5F12.5)') (atom(2*i,3)+atom(2*i−1,3))/20d0−xmin2,(atom(2*i,6)+atom(2*i−1,6))/2d0,(atom(2*i,7)+atom(2*i−1,7))/2d0,(atom(2*i,9)+atom(2*i−1,9))/2d0,(tmp(2*i)+tmp(2*i−1))/2d0
		endif
	enddo
	close(66)
	close(12)
	enddo
end program main