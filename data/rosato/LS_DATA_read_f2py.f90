!07-01-2016: this code reads the line shape database
!Each line shape has 1000 points
!----------------------------------------------------------------
module param_file
	real::density_val(10)=(/1.e13,2.15e13,4.64e13,1.e14,2.15e14,4.64e14,1.e15,2.15e15,4.64e15,1.e16/)
	real::temperature_val(5)=(/.316,1.,3.16,10.,31.6/)
	real::B_val(6)=(/0.,1.,2.,2.5,3.,5./)
end module
!****************************************************************
!program LS_DATA_read
!!Main program
!!- input: 'in.txt' (input parameters)
!!         '.\database\' (database)
!!- output: 'ls.txt' (line shape)
!	implicit none
!	real,allocatable::ls_perp(:),ls_par(:)
!	real,dimension(0:1,0:1,0:1,1000)::w_arr,ls_arr
!	real::Ne,Te,B,theta,wmax
!	integer::nq,npts,iN,iT,iB
!	character(21)::dir
!	character(66)::name
!	call in_param(nq,Ne,Te,B,theta,wmax,npts)
!	allocate(ls_perp(npts),ls_par(npts))
!	call set_bounds(Ne,Te,B,iN,iT,iB)
!	call set_name_file(nq,iN,iT,iB,0,dir,name)	!perpendicular observation
!	call read_file(dir,name,w_arr,ls_arr)
!	call ls_interpol(nq,Ne,Te,B,wmax,npts,w_arr,ls_arr,iN,iT,iB,ls_perp)
!	call set_name_file(nq,iN,iT,iB,1,dir,name)	!parallel observation
!	call read_file(dir,name,w_arr,ls_arr)
!	call ls_interpol(nq,Ne,Te,B,wmax,npts,w_arr,ls_arr,iN,iT,iB,ls_par)
!	call write_ls(npts,wmax,theta,ls_perp,ls_par)
!	deallocate(ls_perp,ls_par)
!end program
!****************************************************************
subroutine in_param(dir_in_len,dir_in,nq,Ne,Te,B,theta,wmax,npts)
!Reads or generates the input parameters in 'in.txt'. MODIFIED TO ONLY WRITE TO A NEW IN.TXT FILE
	! added dir input for absolute paths
	implicit none
	integer, intent(in)::dir_in_len
	! this part is a hack to allow absolute paths of varying lengths to be inputted, it could be bad practice, but seems
	! to work. For strange errors, this is probably a good place to start.
	character(dir_in_len), intent(in) ::dir_in
	integer, intent(in) ::nq,npts
	real, intent(in) ::Ne,Te,B,theta,wmax
	logical::error
	error=				(nq.le.2).or.(nq.ge.8)			&
				.or.	(Ne.lt.1.e13).or.(Ne.gt.1.e16)	&
				.or.	(Te.lt..316).or.(Te.gt.31.6)	&
				.or.	(B.lt.0.).or.(B.gt.5.)
	if(error) then
		write(*,*)'Error!'
		write(*,*)
		stop
	end if

	open(unit=20,file=dir_in//'in.txt',status='old')
	write(20,"('-----------------------------------------------------')")
	write(20,"('Initial principal quantum number         : ',i8)")nq
	write(20,"('Density (cm-3)                           : ',es10.2)")Ne
	write(20,"('Temperature (eV)                         : ',es10.2)")Te
	write(20,"('Magnetic field (T)                       : ',es10.2)")B
	write(20,"('Angle (degrees)                          : ',es10.2)")theta
	write(20,"('-----------------------------------------------------')")
	write(20,"('Delta_omega max (eV)                     : ',es10.2)")wmax
	write(20,"('Number of points                         : ',i8)")npts
	write(20,"('-----------------------------------------------------')")
	close(20)
end subroutine
!****************************************************************
subroutine set_bounds(Ne,Te,B,iN,iT,iB)
!Sets the upper bound of each interval:
!iN corresponds to the value of N2 if Ne is in [N1,N2[ etc.
	use param_file
	implicit none
	real, intent(in) ::Ne,Te,B
	integer, intent(out) ::iN,iT,iB
	logical::error
	do iN=1,10
		if(Ne.lt.density_val(iN)) exit
	end do
	do iT=1,5
		if(Te.lt.temperature_val(iT)) exit
	end do
	do iB=1,6
		if(B.lt.B_val(iB)) exit
	end do
	error=(iN.eq.11).or.(iT.eq.6).or.(iB.eq.7)
	if(error) then
		write(*,*)'Error!'
		write(*,*)
		stop
	end if
end subroutine
!****************************************************************
subroutine set_name_file(nq,iN,iT,iB,iangle,name)
!Sets the name of files that must be used in the database
	implicit none
	integer, intent(in) ::nq,iN,iT,iB,iangle
	character(66), intent(out)::name
	integer::i,iN1,iT1,iB1
	! JSALLCOCK: I have changed the back-slashes to forward slashes here to work on MAC
	i=1
	do iN1=0,1
		do iT1=0,1
			if(iN1*iT1.eq.1) exit
			do iB1=0,1
				name(i:i+10)='ls0'//trim(char(ichar('0')+iN-1+iN1))//trim(char(ichar('0')+ &
						iT-1+iT1))//trim(char(ichar('0')+iB-1+iB1))//trim(char(ichar('0')+ &
						iangle))//'.txt'
				if(iN.eq.10) name(i:i+10)='ls10'//trim(char(ichar('0')+iT-1+iT1))//trim(char(ichar('0')+&
						iB-1+iB1))//trim(char(ichar('0')+iangle))//'.txt'
				i=i+11
			end do
		end do
	end do
end subroutine
!****************************************************************
subroutine read_file(dir_len,dir,name,w_arr,ls_arr)
!Reads the files from database
	implicit none
	integer, intent(in)::dir_len
	! this part is a hack to allow absolute paths of varying lengths to be inputted, it could be bad practice, but seems
	! to work. For strange errors, this is probably a good place to start.
	character(dir_len), intent(in) ::dir
	character(66), intent(in) ::name
	real,dimension(0:1,0:1,0:1,1000), intent(out) ::w_arr,ls_arr
	integer::i,iN1,iT1,iB1,iw
	w_arr=0.
	ls_arr=0.
	i=1
	do iN1=0,1
		do iT1=0,1
			if(iN1*iT1.eq.1) exit
			do iB1=0,1
				open(unit=10,file=trim(dir)//trim(name(i:i+10)),status='old')
				do iw=1,1000
					read(10,*)w_arr(iN1,iT1,iB1,iw),ls_arr(iN1,iT1,iB1,iw)
				end do
				close(10)
				i=i+11	
			end do
		end do
	end do
end subroutine
!****************************************************************
subroutine ls_interpol(nq,Ne,Te,B,wmax,npts,w_arr,ls_arr,iN,iT,iB,ls)
!Calculates the line shape by interpolation
	use param_file
	implicit none
	integer, intent(in) ::nq,npts,iN,iT,iB
	real, intent(in) ::Ne,Te,B,wmax
	real,dimension(0:1,0:1,0:1,1000), intent(in)::w_arr,ls_arr
	real,dimension(npts), intent(out)::ls
	real,dimension(0:1,0:1,0:1,npts)::ls_arr2
	real,dimension(0:1,0:1,npts)::ls_arr3
	real::w
	integer::i,iN1,iT1,iB1,iw
	ls_arr2=0.
	ls_arr3=0.
	if(iB.ne.2) then		!B >= B_val(2)
		do i=1,npts
			w=-wmax+real(i-1)*2.*wmax/real(npts-1)
			do iN1=0,1
				do iT1=0,1
					if(iN1*iT1.eq.1) exit
					do iB1=0,1
						do iw=1,1000
							if(w.lt.(w_arr(iN1,iT1,iB1,iw)*B/B_val(iB-1+iB1))) exit ! JA: this part ensures that the lineshape is zero outside of the tabulated detuning values
						end do
						if((iw.ne.1).and.(iw.lt.1000)) ls_arr2(iN1,iT1,iB1,i)	=			ls_arr(iN1,iT1,iB1,iw-1)						&
																					+	(w*B_val(iB-1+iB1)/B-w_arr(iN1,iT1,iB1,iw-1))		&
																					*	(ls_arr(iN1,iT1,iB1,iw)-ls_arr(iN1,iT1,iB1,iw-1))	&
																					/	(w_arr(iN1,iT1,iB1,iw)-w_arr(iN1,iT1,iB1,iw-1))
					end do
					ls_arr3(iN1,iT1,i)	=		((B-B_val(iB-1))/(B_val(iB)-B_val(iB-1)))	*	(B_val(iB)  /B)	*	ls_arr2(iN1,iT1,1,i)	&
											+	(  (B_val(iB)-B)/(B_val(iB)-B_val(iB-1)))	*	(B_val(iB-1)/B)	*	ls_arr2(iN1,iT1,0,i)
				end do
			end do
		end do
	else if(B.eq.0.) then	!B = 0.
		do i=1,npts
			w=-wmax+real(i-1)*2.*wmax/real(npts-1)
			do iN1=0,1
				do iT1=0,1
					if(iN1*iT1.eq.1) exit
						do iw=1,1000
							if(w.lt.w_arr(iN1,iT1,0,iw)) exit
						end do
						if((iw.ne.1).and.(iw.lt.1000)) ls_arr2(iN1,iT1,0,i)	=			ls_arr(iN1,iT1,0,iw-1)						&
																					+	(w-w_arr(iN1,iT1,0,iw-1))		&
																					*	(ls_arr(iN1,iT1,0,iw)-ls_arr(iN1,iT1,0,iw-1))	&
																					/	(w_arr(iN1,iT1,0,iw)-w_arr(iN1,iT1,0,iw-1))
						ls_arr3(iN1,iT1,i)	=	ls_arr2(iN1,iT1,0,i)
				end do
			end do
		end do
	else					!0. < B < B_val(2)
		do i=1,npts
			w=-wmax+real(i-1)*2.*wmax/real(npts-1)
			do iN1=0,1
				do iT1=0,1
					if(iN1*iT1.eq.1) exit
						do iw=1,1000
							if(w.lt.w_arr(iN1,iT1,0,iw)) exit
						end do
						if((iw.ne.1).and.(iw.lt.1000)) ls_arr2(iN1,iT1,0,i)	=			ls_arr(iN1,iT1,0,iw-1)						&
																					+	(w-w_arr(iN1,iT1,0,iw-1))		&
																					*	(ls_arr(iN1,iT1,0,iw)-ls_arr(iN1,iT1,0,iw-1))	&
																					/	(w_arr(iN1,iT1,0,iw)-w_arr(iN1,iT1,0,iw-1))
						do iw=1,1000
							if(w.lt.(w_arr(iN1,iT1,1,iw)*B/B_val(2))) exit
						end do
						if((iw.ne.1).and.(iw.lt.1000)) ls_arr2(iN1,iT1,1,i)	=			ls_arr(iN1,iT1,1,iw-1)						&
																					+	(w*B_val(2)/B-w_arr(iN1,iT1,1,iw-1))		&
																					*	(ls_arr(iN1,iT1,1,iw)-ls_arr(iN1,iT1,1,iw-1))	&
																					/	(w_arr(iN1,iT1,1,iw)-w_arr(iN1,iT1,1,iw-1))
						ls_arr3(iN1,iT1,i)	=		((B-B_val(iB-1))/(B_val(iB)-B_val(iB-1)))	*	(B_val(iB)  /B)	*	ls_arr2(iN1,iT1,1,i)	&
												+	(  (B_val(iB)-B)/(B_val(iB)-B_val(iB-1)))						*	ls_arr2(iN1,iT1,0,i)
				end do
			end do
		end do	
	end if
	do i=1,npts
		w=-wmax+real(i-1)*2.*wmax/real(npts-1)
		ls(i)=	 3.*log10(Ne/density_val(iN-1))											*	ls_arr3(1,0,i)	&
				+2.*log10(Te/temperature_val(iT-1))										*	ls_arr3(0,1,i)	&
				+(1.-3.*log10(Ne/density_val(iN-1))-2.*log10(Te/temperature_val(iT-1)))	*	ls_arr3(0,0,i)
	end do
end subroutine
!****************************************************************
subroutine write_ls(npts,wmax,theta,ls_perp,ls_par)
!Writes the line shape in 'ls.txt'
	implicit none
	integer::npts
	real::wmax,theta
	real,dimension(npts)::ls_perp,ls_par
	real::w
	integer::i
	theta=theta*3.141593/180.
	open(unit=20,file='ls.txt',status='replace')
	do i=1,npts
		w=-wmax+real(i-1)*2.*wmax/real(npts-1)
		write(20,*)w,ls_perp(i)*sin(theta)*sin(theta)+ls_par(i)*cos(theta)*cos(theta)
	end do
	close(20)
end subroutine
!***************************** END ******************************
