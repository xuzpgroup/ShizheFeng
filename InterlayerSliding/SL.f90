Module ran_mod
 Implicit None
contains
	function reg(x1,y1,x2,y2)
	implicit none
	integer :: flag
	double precision,parameter :: pi = 3.141592653589793239
	double precision :: x1,y1,x2,y2,r,d,reg
	d=dsqrt((x1−x2)**2+(y1−y2)**2)
	r =0.5d0*1.42d0
	if(d.le.2d0*r)then
		reg=((pi*r**2)*2d0*dacos(d/2d0/r)/(2d0*pi)−dsqrt(r**2−(d/2d0)**2)*(d/2d0))*2d0/(pi*r**2d0)
	else
		reg = 0d0
	endif
	end function reg
End Module ran_mod

program main
 use ran_mod
 implicit none
 double precision :: xx,yy,zz,rr,pi=3.141592653589793239,r0,k1,k2,k3,rr0,rr1,rr2,rr3,la
 double precision :: xmax0,xmax,xmin0,xmin,m1,m2,as,theta,l0,tp,e,xl,xmin2,xmax2
 integer :: i,j,jj,jjj,m,natom,natom1,natom2,step,n_min,n_min2,nn,n1,n2,k
 double precision,DIMENSION(:,:) :: box0(3 ,2),atom0(300000,7),box(3,2),atom(300000,5),ini(300000,5),r2(6,2)
 double precision,DIMENSION(:,:) :: min0(300000,4),min(10,10),near(300000,3),near2(300000,6)
 double precision,DIMENSION(:,:) :: t1(4,3),final(300000,10)
 double precision,DIMENSION(:) :: temp(4),min1(3),min2(3),x_temp(2),sl(1000),temp2(7),temp0(5),near3(300000)
 character*200::up,filenum,du,down,in
 double precision,allocatable :: tmp(:),tmp1(:)
 la =1.42
 natom1=0
 natom2=0
 n1=0
 n2=0
 atom0(:,:)=0
 sl(:)=0
!!!!!!!!!initial sate!!!!!!!!!!!
 open(1,file='../l300.data')
	read(1,*)
	read(1,*) natom
	do i =1,3
		read(1,*)
	enddo
	read(1,*) box0(1,1),box0(1 ,2)
	read(1,*) box0(2,1),box0(2 ,2)
	read(1,*) box0(3,1),box0(3 ,2)
	do i=1,8
		read(1,*)
	enddo
	do i=1,natom
	read(1,*) atom0(i,1:7)
		if(atom0(i,2)==1) then
			natom1=natom1+1
		elseif(atom0(i,2)==2) then
			natom2=natom2+1
		endif
	enddo
	write (*,*) natom , natom1 , natom2 ! ! ! atom number: all , bottom layer , top layer
!!!!!!!!!sort atom id!!!!!!!!!
	do jj=natom−1,1,−1
		do j=1,jj
		if(atom0(j,1)>atom0(j+1,1))then
			temp2(:)=atom0(j,:)
			atom0(j,:)=atom0(j+1,:)
			atom0(j+1,:)=temp2(:)
		endif
		enddo
	enddo
	xmax0=maxval(atom0(1:natom1,5))
	xmin0=minval(atom0(1:natom1,5))
!!!!!!!!!find three bonded neighbors from same layer
 near(:,:)=0
 near2(:,:)=0
 do i=1,natom
	n_min=0
	n_min2=0
	do j=1,natom
		if(atom0(i,2)==atom0(j,2))then
			xx=dabs(atom0(i,5)−atom0(j,5))
			if(dabs(atom0(i,6)−atom0(j,6)) .gt. dabs(box0(2,2)−box0(2,1))/2d0)then
				yy=abs(box0(2,2)−box0(2,1))−abs(atom0(i,6)−atom0(j,6))
			else
				yy=abs(atom0(i,6)−atom0(j,6))
			endif
			rr=xx**2+yy**2
			if(rr <1.5**2 .and. rr>0.5**2)then
				n_min=n_min+1
				near(i,n_min)=j
				if(yy/xx<0.2)then
					near3(i)=j
				endif
			elseif(rr>1.6**2.and.rr<2.6**2)then
				n_min2=n_min2+1
				near2(i,n_min2)=j
			endif
		endif
	enddo
!!!!!!!!!!!!!sort by intralayer distance
	allocate(tmp(n_min),tmp1(n_min2))
	do jj=n_min−1,1,−1
		do j=1,jj
		if(near(j,1)>near(j+1,1))then
			tmp(:)=near(j,:)
			near(j,:)=near(j+1,:)
			near(j+1,:)=tmp(:)
		endif
		enddo
	enddo
	do jj=n_min2−1,1,−1
		do j=1,jj
			if(near2(j,1)>near2(j+1,1))then
				tmp1(:)=near2(j,:)
				near2(j,:)=near2(j+1,:)
				near2(j+1,:)=tmp1(:)
			endif
		enddo
	enddo
	DEALLOCATE(tmp,tmp1)
	enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	do k=1,400
!!!!!!!!!!!!!load trajectory
		write(filenum,*)k
		in='../'//trim(adjustl(filenum))//’.dump’
		open(12,file=in,status='unknown')
		write(*,*) in
		final(:,:)=0
		atom(:,:)=0
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
			read(12,*) atom(i,1:5)
		enddo
		close(12)
!!!!!!!!!!!!!sort atom id!!!!!!!!!!!!!
	do jj=natom−1,1,−1
		do j=1,jj
			if(atom(j,1)>atom(j+1,1))then
			temp0(:)=atom(j,:)
			atom(j,:)=atom(j+1,:)
			atom(j+1,:)=temp0(:)
			endif
		enddo
	enddo
	
	xmax=maxval(atom(1:natom1,3))
	xmin=minval(atom(1:natom1,3))
	xmin2=minval(atom(natom1+1:natom,3))
	xmax2=maxval(atom(natom1+1:natom,3))
	
	final(:,:)=0
	
!!!!!!!!!!!!!find three neighbor atoms from different layer
	do i=1,natom
		m=0
		n_min=0
		min0(:,:)=0
		min(:,:)=0
		do j=1,natom
			if(atom(i,2) .ne. atom(j,2))then
				if(dabs(atom(i,3)−atom(j,3)) .gt. dabs(xmax−xmin+la)/2d0)then
				xx=abs(xmax−xmin+la)−dabs(atom(i,3)−atom(j,3))
				else
				xx=dabs(atom(i,3)−atom(j,3))
				endif
				if(dabs(atom(i,4)−atom(j,4)) .gt. dabs(box(2,2)−box(2,1))/2d0)then
				yy=abs(box(2,2)−box(2,1))−abs(atom(i,4)−atom(j,4))
				else
				yy=abs(atom(i,4)−atom(j,4))
				endif
				zz=abs(atom(i,5)−atom(j,5))
				rr=xx**2+yy**2+zz**2
				if(rr<4**2.and.rr>0.2)then
					n_min=n_min+1
					min0(n_min,1)=j
					min0(n_min,2)=rr
				endif
!!!!!!!!!!!!ri
			final(i,8)=final(i,8)+reg(0d0,0d0,xx,yy)
			endif
		enddo
!!!!!!!!!!!!sort by distance
	do jj=n_min−1,1,−1
		do j=1,jj
			if(min0(j,2)>min0(j+1,2))then
				temp(:)=min0(j,:)
				min0(j,:)=min0(j+1,:)
				min0(j+1,:)=temp(:)
			endif
		enddo
	enddo
!!!!!!!!!!!!!!chose nonbonded three atoms (# j ) in on hexagon
	min(1,1)=min0(1,1)
	m=1
	do jjj=1,6
		xx=dabs(atom(near2(min0(1,1),jjj),3)−atom(i,3))
		yy=dabs(atom(near2(min0(1,1),jjj),4)−atom(i,4))
		r2(jjj,1)=xx**2+yy**2
		r2(jjj,2)=near2(min0(1,1),jjj)
	enddo
	
	do jj=6−1,1,−1
		do j=1,jj
			if(r2(j,1)>r2(j+1,1))then
				x_temp(:)=r2(j,:)
				r2(j,:)=r2(j+1,:)
				r2(j+1,:)=x_temp(:)
			endif
		enddo
	enddo
	
!!!!!!!!!!!!!!!!!!!
	min(2,1)=r2(1,2)
	min(3,1)=r2(2,2)
	min(1,2:4)=atom(min(1,1),3:5)
	min(1,5:7)=atom0(min(1,1),5:7)
	min(2,2:4)=atom(min(2,1),3:5)
	min(2,5:7)=atom0(min(2,1),5:7)
	min(3,2:4)=atom(min(3,1),3:5)
	min(3,5:7)=atom0(min(3,1),5:7)
!!!!!!!!!!!!!!find kij fo rx(i)=kij*x(j)
	t1(:,:)=0
	k1=0
	k2=0
	k3=0
	t1(1,1)=atom(i,3)−min(1,2)
	t1(1,2)=atom(i,4)−min(1,3)
	t1(1,3)=atom(i,5)−min(1,4)
	t1(2,1)=min(2,2)−min(1,2)
	t1(2,2)=min(2,3)−min(1,3)
	t1(2,3)=min(2,4)−min(1,4)
	t1(3,1)=min(3,2)−min(1,2)
	t1(3,2)=min(3,3)−min(1,3)
	t1 (3 ,3)=min(3 ,4 )−min(1 ,4 )
	k1=(t1(1,1)*t1 (2,1)+t1(1,2)*t1(2,2)+t1(1,3)*t1(2,3))/(t1(2,1)**2+t1(2,2)**2+t1(2,3)**2)
	k2=(t1(1,1)*t1(3,1)+t1(1,2)*t1(3,2)+t1(1,3)*t1(3,3))/(t1(3,1)**2+t1(3,2)**2+t1(3,3)**2)
	
	if(k1==0 .or. k2==0)then
		final(i,1)=min(1,5)
		final(i,2)=min(1,6)
		final(i,3)=min(1,7)
		final(i,4:8)=0
		else
		t1(4,1)=k1*t1(2,1)+k2*t1(3,1)
		t1(4,2)=k1*t1(2,2)+k2*t1(3,2)
		t1(4,3)=k1*t1(2,3)+k2*t1(3,3)
		k3=(t1(1,1)*t1(4,1)+t1(1,2)*t1(4,2)+t1(1,3)*t1(4,3))/(t1(4,1)**2+t1(4,2)**2+t1(4,3)**2)
		final(i,1)=(1−k1*k3−k2*k3)*min(1,5)+k1*k3*min(2,5)+k2*k3*min(3,5)
		final(i,2)=(1−k1*k3−k2*k3)*min(1,6)+k1*k3*min(2,6)+k2*k3*min(3 ,6)
		final(i,3)=(1−k1*k3−k2*k3)*min(1,7)+k1*k3*min(2,7)+k2*k3*min(3,7)
		m1=(min(1,5)+min(2,5)+min(3,5))/3
		m2=(min(1,6)+min(2,6)+min(3,6))/3
		jj=1
		if(((final(i,1)−m1)**2+(final(i,2)−m2)**2)<((final(i,1)−min(jj,5))**2+(final(i,2)−min(jj,6))**2))then
			final(i,4)=abs(final(i,1)−m1)
			final(i,5)=abs(final(i,2)−m2)
		else
			final(i,4)=abs(final(i,1)−min(jj,5)) !!!slip length in x direction
			final(i,5)=abs(final(i,2)−min(jj,6)) !!!slip length in y direction
		endif
			final(i,6)=sqrt(final(i,4 )**2+final(i,5 )**2) !!!! sli p length
			final(i,7)=final(i,6)/la !!!!!!!!SL~[0 ,1]
	endif
enddo
	
	do i=1,natom
		if((final(i,8)−final(near3(i),8))*(atom(i,3)−atom(near3(i),3)).le.0)then !!!!AB
			final(i,10)=0
			final(i,9)=final(i,7)
		elseif((final(i,8)−final(near3(i),8))*(atom(i,3)−atom(near3(i),3))>0)then !!!!BA
			final(i,10)=1
			final(i,9)=1−final(i,7)
		endif
enddo

!!!!!!!!!!!!!!!output
 write(*,*)"output"
 do i=1,natom
 	if(i<=natom1)then
		final(i,10)=1−final(i,9)
	endif
	
	write(filenum,*) k
	down='down'//trim(adjustl(filenum))//’.dat’
	up='up'//trim(adjustl(filenum))//’.dat’
	open(14,file=down,status='unknown')
	open(15,file=up,status='unknown')
	nn=0
	l0=10.0
	do i=1,natom1
		nn=nn+1
		if(atom0(i,6)<(box0(2,2)−l0) .and. atom0(i,6)>(box0(2,1)+l0)) then
			if(atom0(i,5)>xmin2 .and. atom0(i,5)<xmax) then
				if(isnatomn(final(i,7))) then
				write(14,"(10F16.5)") atom0(i,5)/10,atom0(i,6)/10,0.0,0.0
				else
				write(14,"(10F16.5)") atom0(i,5)/10,atom0(i,6)/10,final(i,7),final(i,9)
				endif
			endif
		endif
	enddo
	
	do i=natom1+1,natom
		nn=nn+1
		if(atom0(i,6)<(box0(2,2)−l0) .and. atom0(i,6)>(box0(2,1)+l0)) then
			if(atom0(i,5)>xmin2 .and. atom0(i,5)<xmax) then
				if(isnatomn(final(i,7))) then
				write(15,"(4F16.5)") atom0(i,5)/10, atom0(i,6)/10 ,0.0 ,0.0
				else
				write(15,"(4F16.5)") atom0(i,5)/10, atom0(i,6)/10, final(i,7) , fi n al ( i , 9 )
				endif
			endif
		endif
	enddo
	close (14)
	close (15)
 enddo
end program main