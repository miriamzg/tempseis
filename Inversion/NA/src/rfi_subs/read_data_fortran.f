




	Subroutine readdata_fortran(fname,obs)

	include 'rfi_param.inc'

	character*40	fname
	real*8		obs(maxdata)
	integer		npoints
	



	open(666, file=fname, status="old")
	read(666,*) npoints
	do i=1,npoints
		read(666,*) time, obs(i)
	end do
	close(666)


	return
	end
