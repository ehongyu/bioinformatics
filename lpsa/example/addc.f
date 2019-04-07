	character bb(1000)*4,cc(1000)*62
	character xx(1000)*66
	nn=0
	do i=1,1000
	read(5,3,end=10)bb(i),cc(i)
3	format(12x,a4,a62)
	if(bb(i)(3:3).eq.bb(i-1)(3:3).and.bb(i)(3:3).ne.' ')then
	bb(i-1)(4:4)='1'
	bb(i)(4:4)='2'
	else if (bb(i)(3:3).eq.bb(i-2)(3:3).and.bb(i)(3:3).ne.' ') then
	bb(i-2)(4:4)='1'
	bb(i)(4:4)='2'
	endif
	cc(i)(57:60)='4PTI'
	nn=nn+1
	end do
10	do i=1,nn
	xx(i)=bb(i)//cc(i)
	write(6,4)i,xx(i)
4	format('ATOM',3x,i4,1x,a66)
	end do
	end
