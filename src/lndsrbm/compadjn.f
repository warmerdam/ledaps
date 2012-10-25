        Program compadjn
	real dy,dx
	character*20 ipa1
	real ang,fac
	call getarg(1,ipa1)
	read(ipa1,'(F10.5)') dx
	call getarg(2,ipa1)
	read(ipa1,'(F10.5)') dy
	fac=atan(1.)/45.
	ang=atan(dx/dy)/fac
	write(6,*) ang
	stop
	end
	
