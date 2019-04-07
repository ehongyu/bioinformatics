	character bb(10000)*4,cc(10000)*62
	character bb1(10000)*4,xx(10000)*66,xx1(10000)*66
	character infile*10, oufile*10
	integer nx(10000),ny(10000)
 	nn=0
	do i=1,10000
	read(5,3,end=10)bb(i),cc(i)
3	format(12x,a4,a62)
	if(bb(i)(1:1).eq.'H'.or.bb(i)(2:2).eq.'H')then
	else
	nn=nn+1
	xx(nn)=bb(i)//cc(i)
	bb1(nn)=bb(i)
	endif
	end do
10	no=nn
	m1=0
	m2=0
	do i=1,nn
	if(bb1(i).eq.' N  ')then
	m1=m1+1
	nx(m1)=i
	else if(bb1(i).eq.' C  ')then
	m2=m2+1
	ny(m2)=i
	end if
	end do
	do i=1,m1-1
	do j=nx(i),nx(i)+1
	xx1(j)=xx(j)
	end do
	do k=nx(i)+2,nx(i)+3
	xx1(k)=xx(k-nx(i)+ny(i)-2)
	end do	
	do l=nx(i)+4,nx(i+1)-1
	xx1(l)=xx(l-2)
	end do
	end do
	do j=nx(m1),nn
	xx1(j)=xx(j)
	end do
	do i=1,nn-1
	write(6,4)i,xx1(i)
4	format('ATOM',3x,i4,1x,a66)
	end do
	stop
	end
