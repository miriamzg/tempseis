
	program readsac

	parameter(max=40000)
	character nin*100, nout*100
	real*4 x(max), y(max)

	nin='../../data/rfi_files/ORF/NA.COCO..LHZ.sac.disp_S2'
	call rsac1(nin, x, npts, beg, dt, max, nerr)

* In the subroutine above, nin, x, and max are input, and
*	the rest are output.

	write(6,61) dt
61	format('The sampling interval is ', f6.3, ' seconds')
	print*,'The number of points in the file is ', npts
	print*,'The 10th data point is ', x(10)

** The array x is read. Now let us do something about it and
*	write the result back as a new SAC  file

	do i=1,npts
	y(i)=abs(x(i))
	enddo

	nout='output.sac'
	beg=0.0
	call wsac1(nout, y, npts, beg, dt, nerr)

	stop
	end

***********************************************
** To compile the program, copy and paste the
*following 3 lines into a file called Makefile and type make
* Of course, you need to remove the * and make sure that
* the 2nd line is a single Tab

*readsac.exe:readsac.f
*	g77 -g -m32 readsac.f -o readsac.exe \
*/home/sgao/lib/sacio.a 
