	character*60 type,t
	open(7,file='4pti.pdb',status='old',readonly)
	mm=0
	xx=0.
	do i=1,1000
	read(5,15,end=18)type,a,b,c
	read(7,15,end=18)t,a1,b1,c1
15	format(13x,a3,14x,3f8.3)
	if (i.ge.94.and.i.le.124) then
	if (t.eq.'N  '.or.t.eq.'C  '.or.t.eq.'O  '.or.t.eq.'CA '.or.
     1  t.eq.' N '.or.t.eq.' C '.or.t.eq.' O '.or.t.eq.' CA') then
	mm=mm+1
	xx=xx+(a1-a)**2.+(b1-b)**2.+(c1-c)**2.
	endif
	endif
	end do
18	xx=(xx/mm)**.5
	close(5)
	close(7)
	print*,"rmsd=",xx
	print*,mm
	stop
	end
