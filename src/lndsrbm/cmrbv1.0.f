         Program readl7sr
         
         character*200 filein,fileth
         integer ii
         integer sfstart, sfselect, sfrdata, sfendacc, sfend
	 integer sfscatt
         integer sd_id, sds_id, sds_index, status
         integer start(5), edges(5), stride(5)
         integer nstart(2), nedges(2), nstride(2)
         integer DFACC_READ,DFACC_RDWR,DFNT_CHAR8
         parameter (DFACC_READ = 1)
         parameter (DFACC_RDWR = 3)
	 integer nc,nr,mnc,mnr
         parameter (DFNT_CHAR8 = 4)
	 parameter (mnc=9000)
	 parameter (mnr=9000)

         integer(2), allocatable :: qa(:,:)
         integer(2), allocatable :: temp(:,:)
         integer(2), allocatable :: band1(:,:)
         integer(2), allocatable :: band3(:,:)
         integer(2), allocatable :: band2(:,:)
         integer(2), allocatable :: band5(:,:)
	 integer*2 cloud,cloudadja,cloudshad,tmpbit,snow
	 real tclear,ts,tv,fs,tcloud,cfac,cldh,pixsize,tna
	 real fack,facl
	 character*80 sds_name
	 integer rank,data_type
	 integer n_attrs
	 integer dim_sizes(5)
	 integer mband5,mband5k,mband5l
	 double precision mclear
	 real pclear
	 integer*8 nbclear,nbval,nbcloud
	 real t6
	 integer cldhmin,cldhmax
	 
	 cloud=8
	 cloudadja=12
	 cloudshad=9
	 snow=10
	 tmpbit=13
	 cfac=6.
	 dtr=atan(1.)/45.
	 pixsize=28.5 
	 
	 
         
         call getarg(1,filein)
         ii=index(filein," ")-1
c         call getarg(2,fileth)
c         jj=index(fileth," ")-1
	 write(6,*) "temp clear [k] , ts,tv,fs,truenorthadjustement "
	 read(5,*)  tclear,ts,tv,fs,tna
	 tclear=tclear-273.0
	 
        
                 
cc Read from the HDF file 
         sd_id= sfstart(filein(1:ii),DFACC_RDWR)
         
         write(6,*) "sd_id", sd_id
         if (sd_id.lt.0) then
         write(6,'(A25,1X,A80)')
     &  "WW can not open input ",filein
          stop
	  endif
c reading qa       
c allocate memory for qa
       sds_index = 7
       sds_id    = sfselect(sd_id, sds_index)
       write(6,*) "sds_id", sds_id
       status= sfginfo(sds_id, sds_name, rank, dim_sizes, data_type,
     s  n_attrs)
       write(6,*) "sdsname ",sds_name
       write(6,*) dim_sizes(1),dim_sizes(2)
       nc= dim_sizes(1)
       nr=dim_sizes(2)
	allocate (qa(nc,nr),stat=ierr)
	allocate (temp(nc,nr),stat=ierr)
	allocate (band1(nc,nr),stat=ierr)
	allocate (band2(nc,nr),stat=ierr)
	allocate (band3(nc,nr),stat=ierr)
	allocate (band5(nc,nr),stat=ierr)

cc read qa  data 
       start(1)=0
       start(2) = 0
       edges(1) = nc
       edges(2) = nr
       stride(1) = 1
       stride(2) = 1
       status = sfrdata(sds_id, start, stride, edges,qa)
       write(6,*) "status", status
c reading temperature       
       sds_index = 8
       sds_id    = sfselect(sd_id, sds_index)
       write(6,*) "sds_id", sds_id
       status = sfrdata(sds_id, start, stride, edges,temp)
       write(6,*) "status", status
c reading band1 (blue)
       sds_index =0
       sds_id    = sfselect(sd_id, sds_index)
       write(6,*) "sds_id", sds_id
       status = sfrdata(sds_id, start, stride, edges,band1)
       write(6,*) "status", status
c reading band3 (red)
       sds_index =2
       sds_id    = sfselect(sd_id, sds_index)
       write(6,*) "sds_id", sds_id
       status = sfrdata(sds_id, start, stride, edges,band3)
       write(6,*) "status", status
c reading band5 (swir)
       sds_index =4
       sds_id    = sfselect(sd_id, sds_index)
       write(6,*) "sds_id", sds_id
       status = sfrdata(sds_id, start, stride, edges,band5)
       write(6,*) "status", status
c reading band2 (green)
       sds_index =1
       sds_id    = sfselect(sd_id, sds_index)
       write(6,*) "sds_id", sds_id
       status = sfrdata(sds_id, start, stride, edges,band2)
       write(6,*) "status", status


c printing some values for checking out
       i=3953
       j=3153
       write(6,*) " i,j ", i,j
       write(6,*) "qa(j,i) ",qa(j,i)
       write(6,*) "btest ", btest(qa(j,i),cloud),
     s   btest(qa(j,i),cloudshad)
       write(6,*) "temp(j,i) ",temp(j,i)
       write(6,*) "band1(j,i) ",band1(j,i)
       write(6,*) "band3(j,i) ",band3(j,i)
       write(6,*) "band5(j,i) ",band5(j,i)

       
       
c note i is the line and j is the collumn  
       write(6,*) "resetting cloud,cloudadja,cloudshad bits,tmpbit"
       do i=1,nr
       do j=1,nc
c reset cloud,cloud shadow,adjacent cloud bit
       if (btest(qa(j,i),cloud)) qa(j,i)=qa(j,i)-(2**cloud)
       if (btest(qa(j,i),cloudadja)) qa(j,i)=qa(j,i)-(2**cloudadja)
       if (btest(qa(j,i),cloudshad)) qa(j,i)=qa(j,i)-(2**cloudshad)
       if (btest(qa(j,i),tmpbit)) qa(j,i)=qa(j,i)-(2**tmpbit)
       enddo
       enddo
       
       
c update the cloud mask
       write(6,*) "updating cloud mask"
       nbclear=0
       mclear=0.
       nbval=0
       nbcloud=0
      
c compute the average temperature of the clear data       
       do i=1,nr
       do j=1,nc
       if (band1(j,i).ne.-9999) then
       nbval=nbval+1
       anom=band1(j,i)-band3(j,i)/2.
       t6=temp(j,i)/100.
       if (btest(qa(j,i),snow)) then
       continue
       else
       if ((anom.gt.300.).and.(band5(j,i).gt.300)
     s   .and.(t6.lt.tclear)) then
       qa(j,i)=qa(j,i)+(2**cloud)
       nbcloud=nbcloud+1
       else
       if ((band1(j,i).gt.3000).and.(t6.lt.tclear)) then
       qa(j,i)=qa(j,i)+(2**cloud)
       nbcloud=nbcloud+1
       else
       mclear=mclear+t6/10000.
       nbclear=nbclear+1
       endif
       endif
       endif
       endif
       enddo
       enddo
c 
        if (nbclear.gt.0) then
        mclear=mclear*10000./nbclear
	endif
	pclear=(nbclear*100./nbval)
	write(6,*) nbclear-nbval,nbcloud
	write(6,*) "average clear temperature  %clear", mclear,pclear
	if (pclear.gt.5.) then
	tclear=mclear       
        endif
	     
c update the adjacent cloud bit
       write(6,*) "updating adjacent cloud bit"
       do i=1,nr
       do j=1,nc
       if (btest(qa(j,i),cloud)) then
c       write(6,*) "i found a cloud at ",i,j,"with t=",temp(j,i)/100.
       do k=i-5,i+5
       do l=j-5,j+5
       if ((k.ge.1).and.(k.le.nr).and.(l.ge.1).and.(l.le.nc)) then
       if (btest(qa(l,k),cloudadja).OR.btest(qa(l,k),cloud)) then
       continue
       else
       qa(l,k)=qa(l,k)+(2**cloudadja)
       endif
       endif
       enddo
       enddo
       endif
       enddo
       enddo
       
C compute cloud shadow  
       write(6,*) "looking for cloud shadow"
       fack=cos(fs*dtr)*tan(ts*dtr)/pixsize
       facl=sin(fs*dtr)*tan(ts*dtr)/pixsize
       do i=1,nr
       do j=1,nc
       if ((i.eq.5948).and.(j.eq.5291)) then
       write(6,*) "we are at 5948,5291"
       endif
       if (btest(qa(j,i),cloud)) then
           tcloud=temp(j,i)/100.
           cldh=(tclear-tcloud)*1000./cfac
	   if (cldh.lt.0.) cldh=0.
	   cldhmin=cldh-1000.
	   cldhmax=cldh+1000.
	   mband5=9999
           do icldh=cldhmin/10,cldhmax/10
	    cldh=icldh*10.
            k=i+fack*cldh
            l=j-facl*cldh
            if ((k.ge.1).and.(k.le.nr).and.(l.ge.1)
     &	    .and.(l.le.nc)) then
            if ((band5(l,k).lt.800).and.
     &	    ((band2(l,k)-band3(l,k)).lt.100)) then
            if (btest(qa(l,k),cloudadja).OR.btest(qa(l,k),cloud)
     &       .or.btest(qa(l,k),cloudshad)) then
             continue
             else
c store the value of band5 as well as the l and k value
             if (band5(l,k).lt.mband5) then
	     mband5=band5(l,k)
	     mband5k=k
	     mband5l=l
	     endif    
             endif
             endif
             endif
	     enddo
	     if (mband5.lt.9999) then
	     l=mband5l
	     k=mband5k
             qa(l,k)=qa(l,k)+(2**cloudshad)
	     endif
         endif
	 enddo
	 enddo
c dilate the cloud shadowc
c       goto 99
       write(6,*) "dilating cloud shadow"
       do i=1,nr
       do j=1,nc
       if (btest(qa(j,i),cloudshad)) then
       do k=i-3,i+3
       do l=j-3,j+3
       if ((k.ge.1).and.(k.le.nr).and.(l.ge.1).and.(l.le.nc)) then
       if (btest(qa(l,k),cloudadja).OR.btest(qa(l,k),cloud)
     &   .OR.btest(qa(l,k),cloudshad)) then
       continue
       else
       if (btest(qa(l,k),tmpbit)) then
       continue
       else
       qa(l,k)=qa(l,k)+(2**tmpbit)
       endif
       endif
       endif
       enddo
       enddo
       endif
       enddo
       enddo
c update the cloud shadow       
       write(6,*) "updating cloud shadow"
       do i=1,nr
       do j=1,nc
       if (btest(qa(j,i),tmpbit)) then
       qa(j,i)=qa(j,i)+(2**cloudshad)
       qa(j,i)=qa(j,i)-(2**tmpbit)
       endif
       enddo
       enddo
       
       
c updating the qa and closing
 99     sds_index=7
        sds_id= sfselect(sd_id,sds_index)
        status=sfwdata(sds_id, start, stride,edges,qa)
        write(6,*) "status", status
	status = sfendacc(sds_id)
        write(6,*) "status", status
c update the global attribute to reflect cloud processing
         status=sfscatt(sd_id, "Cloud Mask Algo Version ",DFNT_CHAR8,
     s  22, "CMReflectanceBasedv1.0")
         write(6,*) "status write attribute", status
c close HDF file
       status = sfend(sd_id)
       deallocate(qa,temp,band1,band3,band2,band5)
       stop
       end
