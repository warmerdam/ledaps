         program comptemp
	 real temp(5)
	 real time(5)
	 real sctime
	 character*20 ipa1
	 
	 call getarg(1,ipa1)
	 read(ipa1,'(F10.4)') sctime
	 
	 do i=1,4
	 read(5,*) temp(i)
	 time(i)=(i-1)*6
	 enddo
	 temp(5)=temp(1)
	 time(5)=24.0
	 if (sctime.lt.0.01) sctime=0.01
	 
	 i=0
    5    i=i+1
	 if (sctime.gt.time(i)) goto 5
	 slp=(temp(i)-temp(i-1))/6.
	 sctemp=temp(i-1)+slp*(sctime-time(i-1))
	 write(6,*) sctemp
	 
	 
	 stop
	 end
	 
	 
	 
	 
	 
	 
