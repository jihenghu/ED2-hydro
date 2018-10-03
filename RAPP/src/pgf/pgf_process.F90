!==========================================================================================!
!==========================================================================================!
!  Subroutine pgf_process()                                                    !
!                                                                                          !
!     This subroutine will read the pgf dataset in netcdf format, process the data     	   !
! , and write them into .h5 file that ED2 model can read								   !
!------------------------------------------------------------------------------------------!
subroutine pgf_process(year,pgf_res)
	use mod_maxdims		, only : maxstr
	use mod_ioopts		, only : inpath,	&
								 outpref,  	&
								 lonw,lone, &
								 lats,latn
    use mod_netcdf     , only : ncid             & ! intent(out)
                             , ndimensions      & ! intent(out)
                             , nvariables       & ! intent(out)
                             , nglobals         & ! intent(out)
                             , unlimiteddimid   & ! intent(out)
                             , xtype            & ! intent(out)
                             , timeid           & ! intent(out)
                             , dummy_vname      & ! intent(out)
                             , ndims            & ! intent(out)
                             , dimids           & ! intent(out)
                             , natts            & ! intent(out)
                             , dimglobal        & ! intent(out)
                             , globalid         ! ! intent(out)
   use netcdf
   use mod_ncdf_globio, only : ncdf_load_err
   use hdf5
   use hdf5_utils

	implicit none
	!----------------- Argument
	integer, 		intent(in)		:: 		year
	real,   		intent(in)		:: 		pgf_res  ! spatial resolution
	!----------------- local parameters
	integer,		parameter		:: 		pgf_var_num = 7
	real,		parameter		:: 		pgf_t_res = 3600. * 3. ! temporal resolution in sec
	integer,	dimension(13), parameter :: 	month_first_day = &
		(/1,32,60,91,121,152,182,213,244,274,305,335,366/)
	integer,	dimension(13), parameter :: 	month_first_day_leap = &
		(/1,32,61,92,122,153,183,214,245,275,306,336,367/)
	! dlwrf, dswrf, prcp, press, shum, tas, wind
	!----------------- local variables
	character(len=10),	dimension(12)	:: month_str
	integer							::    var_id
	integer							::    i_var
	integer							::    i_lon, i_lat
	integer							::    iatt
	character(len=maxstr)			::    var_name
	character(len=maxstr)			::    output_var_name
	character(len=maxstr)			::    file_name, output_name,dim_name
	character(len=maxstr)			::    att_value
	integer							::	  i_month
	integer							::	  x_num,y_num,t_num,x_num_w,x_num_e
	integer							::    x_start,y_start,t_start,x_start_w,x_start_e
	real							::    lon_start, lat_start
    logical                                     :: file_is_there
    integer                                     :: ierr
	integer										:: dim_len
    logical                         ::    data_buff_flag
	real, pointer, dimension(:,:,:)	::    data_buff
	real, pointer, dimension(:,:,:)	::    data_buff_w ! for west hemisphere
	real, pointer, dimension(:,:,:)	::    data_buff_e ! for east hemisphere
	real, pointer, dimension(:,:,:)	::    vbdsf_buff,nbdsf_buff,vddsf_buff,nddsf_buff
	real, pointer, dimension(:,:,:)	::    out_buff
	real, pointer, dimension(:,:)	::	  lon_matrix, lat_matrix
	real, pointer, dimension(:)		:: 	  lat_buff
   integer, dimension(2)                  :: grid_dims
   integer, dimension(3)                  :: var_dims

	!--------------------------------------------------------

    if (pgf_res /= 1.0 .and. pgf_res /= 0.5) then
        print*,'pgf_res',pgf_res
        call fatal_error('pgf_res is wrong, should be either 0.5 or 1.0 deg'    &
                      ,'pgf_process','pgf_process.f90')
    endif



	! some constants
	month_str(1) = "JAN"
	month_str(2) = "FEB"
	month_str(3) = "MAR"
	month_str(4) = "APR"
	month_str(5) = "MAY"
	month_str(6) = "JUN"
	month_str(7) = "JUL"
	month_str(8) = "AUG"
	month_str(9) = "SEP"
	month_str(10) = "OCT"
	month_str(11) = "NOV"
	month_str(12) = "DEC"

	! calculate dimensions we need to extract and write the data
	!-------------------------------------------------------- 
	! lon and lat
	! pgf data starts from -89.75 and 0.25
    ! however, we will write from -180 to 180 in longitude

    if (lone > 180.) then
        call fatal_error('longitude should not be larger than 180.'    &
                      ,'pgf_process','pgf_process.f90')
    endif

	x_num = (nint(lone / pgf_res) - nint(lonw / pgf_res))  ! longitude
    x_num_w = min(0,nint(lone / pgf_res)) - min(0, nint(lonw / pgf_res))
    x_num_e = max(0,nint(lone / pgf_res)) - max(0, nint(lonw / pgf_res))
    print*,'x_num_w',x_num_w,'x_num_e',x_num_e

    if (lonw < 0. .and. lone > 0.) then
        ! separated in two hemisphere
		x_start = nint((360. + lonw) / pgf_res) + 1
		lon_start = (x_start - 1) * pgf_res + pgf_res/2. - 360.
        x_start_w = nint((360. + lonw) / pgf_res) + 1
        x_start_e = 1
    elseif (lonw < 0. .and. lone < 0.) then
		x_start = nint((360. + lonw) / pgf_res) + 1
		lon_start = (x_start - 1) * pgf_res + pgf_res/2. - 360.
        x_start_w = 0
        x_start_e = 0

    else
		x_start = nint(lonw / pgf_res + 1)
		lon_start = (x_start - 1) * pgf_res + pgf_res/2.
    endif

    print*,'x_start',x_start,'x_start_w',x_start_w,'x_start_e',x_start_e

	y_num = (nint(latn / pgf_res) - nint(lats / pgf_res))
	y_start = nint((lats + 90) / pgf_res) + 1
	lat_start = (y_start - 1) * pgf_res - 90 + pgf_res/2.

	allocate(lon_matrix(x_num,y_num))
	allocate(lat_matrix(x_num,y_num))

	do i_lon = 1 , x_num
		lon_matrix(i_lon,:) = lon_start + (i_lon - 1) * pgf_res
	enddo

	do i_lat = 1 , y_num
		lat_matrix(:,i_lat) = lat_start + (y_num  - i_lat) * pgf_res
	enddo


	! loop over all the months to save memory
	month_loop:do i_month = 1, 12
		! calculate time dimensions
		if (isLeap(year)) then
			t_start = (month_first_day_leap(i_month) - 1) * 8 + 1   ! every day has 8 data points
			t_num = (month_first_day_leap(i_month+1) - &
					 month_first_day_leap(i_month)) * 8   ! every day has 8 data points
		else
			t_start = (month_first_day(i_month) - 1) * 8 + 1   ! every day has 8 data points
			t_num = (month_first_day(i_month+1) - &
					 month_first_day(i_month)) * 8   ! every day has 8 data points
		endif
		
		allocate(data_buff(x_num,y_num,t_num))
        data_buff_flag = .false.
        ! to see whether it is necessary to allocate data_buff_w and data_buff_e
        if (x_num_w > 0 .and. x_num_e > 0) then
            allocate(data_buff_w(x_num_w,y_num,t_num))
            allocate(data_buff_e(x_num_e,y_num,t_num))
            data_buff_flag = .true.
        endif

        
		allocate(out_buff(t_num,x_num,y_num))

		! create output file
   		write(output_name,fmt='(a,i4.4,a,a)')                                                   &
                      trim(outpref),year,trim(month_str(i_month)),'.h5'
		
   		write (unit=*,fmt='(a,1x,2a)') '         [|] Writing file:',trim(output_name),'...'

  	    call shdf5_open_f(trim(output_name),'W',1)

   		! write lat and lon
		grid_dims = (/x_num, y_num/)
		call shdf5_orec_f( 2, grid_dims, trim('lon'),rvara=lon_matrix)
		call shdf5_orec_f( 2, grid_dims, trim('lat'),rvara=lat_matrix)



		! get lat
			! get file info

 !     		write(file_name,fmt='(4a,i4.4,a,i4.4,a)')                                                   &
 !                     trim(inpath),'/',trim('dlwrf'),'_3hourly_',year,'-',year,'.nc'
 !     		inquire(file=file_name,exist=file_is_there)
!			
!			if(.not. file_is_there) then
!				print*,'file does not exit. Program stops'
!				stop;
!			endif
			
			! open the file

 !     		ierr = nf90_open(file_name,NF90_NOWRITE,ncid)
 !     		if (ierr /= NF90_NOERR) then
!		         call ncdf_load_err(ierr)
!				 print*,'Error when opening the file:'
!				 print*,file_name
!		    end if

!      		ierr = nf90_inquire(ncid,ndimensions,nvariables,nglobals,unlimiteddimid)

			! get variable id
!      		ierr = nf90_inq_varid(ncid,'lat',var_id)
!		    if (ierr /= NF90_NOERR) then
!        	 call ncdf_load_err(ierr)
!				 print*,'Error when inquiring the variable ',var_name,' in the file:'
!				 print*,file_name
!    		end if
			
!         	ierr = nf90_inquire_variable(ncid,var_id,dummy_vname,xtype,ndims,dimids,natts)
!			print*,'ndims',ndims

			! get variable
!			allocate(lat_buff(y_num))
!			ierr = nf90_get_var(ncid,var_id,lat_buff,start = (/y_start/),count= (/y_num/))
!			print*,'lat_buff',lat_buff
!			deallocate(lat_buff)
		! loop over the seven variableps
			
		do i_var = 1, 7
			select case (i_var)
				case (1)   ! dlwrf
					var_name = 'dlwrf'
					output_var_name = 'dlwrf'
				case (2)   ! dswrf
					var_name = 'dswrf'
					output_var_name = 'dswrf'
				case (3)   ! prcp
					var_name = 'prcp'
					output_var_name = 'prate'
				case (4)   ! pres
					var_name = 'pres'
					output_var_name = 'pres'
				case (5)   ! shum
					var_name = 'shum'
					output_var_name = 'sh'
				case (6)   ! tas
					var_name = 'tas'
					output_var_name = 'tmp'
				case (7)   ! wind
					var_name = 'wind'
					output_var_name = 'ugrd'
			end select



			! get file info

      		write(file_name,fmt='(4a,i4.4,a,i4.4,a)')                                                   &
                      trim(inpath),'/',trim(var_name),'_3hourly_',year,'-',year,'.nc'
      		inquire(file=file_name,exist=file_is_there)
			
			if(.not. file_is_there) then
				print*,'file does not exit. Program stops'
				print*,file_name
				stop;
			endif
			
			! open the file

      		ierr = nf90_open(file_name,NF90_NOWRITE,ncid)
      		if (ierr /= NF90_NOERR) then
		         call ncdf_load_err(ierr)
				 print*,'Error when opening the file:'
				 print*,file_name
		    end if

      		ierr = nf90_inquire(ncid,ndimensions,nvariables,nglobals,unlimiteddimid)

			! get variable id
      		ierr = nf90_inq_varid(ncid,var_name,var_id)
		    if (ierr /= NF90_NOERR) then
        	 call ncdf_load_err(ierr)
				 print*,'Error when inquiring the variable ',var_name,' in the file:'
				 print*,file_name
    		end if
			
!         	ierr = nf90_inquire_variable(ncid,var_id,dummy_vname,xtype,ndims,dimids,natts)

			! get variable
            ! check whether necessary to read two databuff
            if (x_num_w > 0 .and. x_num_e > 0) then
                ! read two data buff
			    ierr = nf90_get_var(ncid,var_id,data_buff_w,&
                                    start = (/x_start_w,y_start,t_start/),&
                                    count= (/x_num_w,y_num,t_num/))
			    ierr = nf90_get_var(ncid,var_id,data_buff_e,&
                                    start = (/x_start_e,y_start,t_start/),&
                                    count= (/x_num_e,y_num,t_num/))
                ! save into data_buff
                data_buff(1:x_num_w,:,:) = data_buff_w(:,:,:)
                data_buff(x_num_w+1:x_num,:,:) = data_buff_e(:,:,:)
            else

    			ierr = nf90_get_var(ncid,var_id,data_buff,&
                                    start = (/x_start,y_start,t_start/),&
                                    count= (/x_num,y_num,t_num/))
            endif

		    ! close the file
			ierr = nf90_close(ncid)
		    if (ierr /= NF90_NOERR) then
        	 call ncdf_load_err(ierr)
				 print*,'Error when closing the file:'
				 print*,file_name
    		end if
			!print*,var_name,data_buff(1,:,1)
			!------------------------  data process and write  ------------!

				

			if (i_var == 2) then
				! short wave, we have to process them into different bands...
         		call swap_dim(x_num,y_num,t_num,x_num,y_num,0.,data_buff,out_buff)
				! allocate the buffer for four different components of dswrf
				allocate(vbdsf_buff(t_num,x_num,y_num))
				allocate(vddsf_buff(t_num,x_num,y_num))
				allocate(nbdsf_buff(t_num,x_num,y_num))
				allocate(nddsf_buff(t_num,x_num,y_num))

				! attibute values to the four buffer
				call split_radiation(x_start,x_num,y_start,y_num,t_start,t_num,&
									 out_buff,lat_matrix,lon_matrix, &
									 vbdsf_buff,vddsf_buff,nbdsf_buff,nddsf_buff)
!				print*,'quotient',(vbdsf_buff(1,1,1) + vddsf_buff(1,1,1) + &
!						nbdsf_buff(1,1,1) + nddsf_buff(1,1,1)) / out_buff(1,1,1), &
!						'total radiation', out_buff(1,1,1)
				! write them into h5 file

				var_dims = (/t_num, x_num, y_num/)
				call shdf5_orec_f( 3, var_dims, trim('vbdsf'),rvara=vbdsf_buff)
				call shdf5_orec_f( 3, var_dims, trim('vddsf'),rvara=vddsf_buff)
				call shdf5_orec_f( 3, var_dims, trim('nbdsf'),rvara=nbdsf_buff)
				call shdf5_orec_f( 3, var_dims, trim('nddsf'),rvara=nddsf_buff)

				deallocate(nddsf_buff)
				deallocate(nbdsf_buff)
				deallocate(vddsf_buff)
				deallocate(vbdsf_buff)
			else
			! rewrite the data_buff to output_buff because in h5 file time
			! dimension comes first and lat dimension starts from the largest
         		call swap_dim(x_num,y_num,t_num,x_num,y_num,0.,data_buff,out_buff)
				var_dims = (/t_num, x_num, y_num/)
				call shdf5_orec_f( 3, var_dims, trim(output_var_name),rvara=out_buff)

			endif


		
		enddo

		! write vgrd as dummy var
		out_buff = 0.
		var_dims = (/t_num, x_num, y_num/)
		call shdf5_orec_f( 3, var_dims, trim('vgrd'),rvara=out_buff)


		! write hgt
!		out_buff = 150.
!		var_dims = (/t_num, x_num, y_num/)
!		call shdf5_orec_f( 3, var_dims, trim('hgt'),rvara=out_buff)

        if (data_buff_flag) then
            deallocate(data_buff_w)
            deallocate(data_buff_e)
        endif

		deallocate(out_buff)
		deallocate(data_buff)
   call shdf5_close_f()
	enddo month_loop

		deallocate(lat_matrix)
		deallocate(lon_matrix)

contains
logical function isLeap(year)
	implicit none
	integer,   intent(in)	:: year
	isLeap = (mod(year,4) == 0 .and. mod(year,100) > 0) .or. &
			 (mod(year,400) == 0 .and. mod(year,100) == 0)
end function isLeap

end subroutine


subroutine swap_dim(mxi,myi,mtp,mxo,myo,eoff,arin,arout)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                        , intent(in)  :: mxi   ! # of X points, input array
   integer                        , intent(in)  :: myi   ! # of Y points, input array
   integer                        , intent(in)  :: mtp   ! # of T points, both arrays
   integer                        , intent(in)  :: mxo   ! # of X points, output array
   integer                        , intent(in)  :: myo   ! # of Y points, output array
   integer                        , intent(in)  :: eoff  ! Edge offset of output array
   real   , dimension(mxi,myi,mtp), intent(in)  :: arin  ! Input array  (work)
   real   , dimension(mtp,mxo,myo), intent(out) :: arout ! Output array (buffer)
   !----- Local variables. ----------------------------------------------------------------!
   integer                                      :: xi    ! Index of input array, X
   integer                                      :: yi    ! Index of input array, Y
   integer                                      :: tt    ! Time index, both arrays
   integer                                      :: xo    ! Index of output array, X
   integer                                      :: yo    ! Index of output array, Y
   logical                                      :: okx   ! Flag for correct X dimensions
   logical                                      :: oky   ! Flag for correct Y dimensions
   !---------------------------------------------------------------------------------------!

   !----- Quick dimension check. ----------------------------------------------------------!
   okx = mxo == mxi - 2 * eoff
   oky = myo == myi - 2 * eoff

   !----- Run only when the input vs. output dimensions are consistent. -------------------!
   if (okx .and. oky) then
      yoloop: do yo=1,myo
         yi = myi - (yo + eoff) + 1

         xoloop: do xo=1,mxo
            xi = xo + eoff

            ttloop: do tt=1,mtp
               arout(tt,xo,yo) = arin(xi,yi,tt) 
            end do ttloop

         end do xoloop
      end do yoloop
   else
      write (unit=*,fmt='(a,1x,i5)') 'MXI  =',mxi
      write (unit=*,fmt='(a,1x,i5)') 'MYI  =',myi
      write (unit=*,fmt='(a,1x,i5)') 'MXO  =',mxo
      write (unit=*,fmt='(a,1x,i5)') 'MYO  =',myo
      write (unit=*,fmt='(a,1x,i5)') 'EOFF =',eoff
      write (unit=*,fmt='(a,1x,l1)') 'MXO = MXI-2EOFF?',okx
      write (unit=*,fmt='(a,1x,l1)') 'MYO = MYI-2EOFF?',oky
      call fatal_error('Mismatch between input and output grid (check values above...)'    &
                      ,'swap_dim','pgf_process.f90')
   end if

   return
end subroutine swap_dim


subroutine split_radiation(x_start,x_num,y_start,y_num,t_start,t_num,&
									 out_buff,lat_matrix,lon_matrix, &
									 vbdsf_buff,vddsf_buff,nbdsf_buff,nddsf_buff)
implicit none
!----------------------- input------------------
integer,			intent(in)				::	x_start
integer,			intent(in)				::	x_num
integer,			intent(in)				::	y_start
integer,			intent(in)				::	y_num
integer,			intent(in)				::	t_start
integer,			intent(in)				::	t_num
real,		dimension(t_num,x_num,y_num),	intent(in)	:: out_buff	
real,		dimension(x_num,y_num),	intent(in)	:: lon_matrix	
real,		dimension(x_num,y_num),	intent(in)	:: lat_matrix	
real,		dimension(t_num,x_num,y_num),	intent(out)	:: vbdsf_buff	
real,		dimension(t_num,x_num,y_num),	intent(out)	:: vddsf_buff	
real,		dimension(t_num,x_num,y_num),	intent(out)	:: nbdsf_buff	
real,		dimension(t_num,x_num,y_num),	intent(out)	:: nddsf_buff	
!--------------------  constants--------------------
real,				parameter				:: solar_cons = 1.3533e3   ! W/m2
real,				parameter				:: pi = 3.141592653
real,				parameter				:: pio180 = pi / 180.
!-------------------   local variable
integer										:: i_time, i_x, i_y
integer										:: doy
real										:: utc_sec
real										:: d0, d02
real										:: solfac
real										:: t1, declin, t2, eqn_of_time
real										:: sun_lon,sunx,suny,sunz
real										:: direct_frac
! vars depdent on lat and lon
real,	pointer, dimension(:,:)				:: wnx,wny,wnz,cosz,rshort_clear

allocate(wnx(x_num,y_num))
allocate(wny(x_num,y_num))
allocate(wnz(x_num,y_num))
allocate(cosz(x_num,y_num))
allocate(rshort_clear(x_num,y_num))

! loop over all the time frame
do i_time = 1, t_num
	! find out doy and utc for this time
	doy = (t_start - 1) / 8 + 1 + floor(real(i_time) / 8.)
	utc_sec = mod(i_time-1,8) * 3 * 3600    ! pgf is 3hourly data
	
	d0 = 2 * pi * (doy - 1) / 365
	d02 = d0 * 2
	solfac = 1.00011 + 0.034221 * cos(d0) + 0.00128 * sin(d0) + &
		    0.000719 * cos(d02) + 0.000077 * sin(d02)
	t1 = pi * 2 * doy / 366
	declin = 0.322003 - 22.971 * cos(t1)- 0.357898 * cos(t1 * 2) - &
		    0.14398 * cos(t1 * 3) + &
		    3.94638 * sin(t1) + &
		    0.019334 * sin(t1 * 2) + &
		    0.05928 * sin(t1 * 3)

	t2 = (279.134 + 0.985647 * doy) * pio180
	eqn_of_time = 5.0323 - 100.976 * sin(t2) + &
			    595.275 * sin(t2 * 2) + 3.6858 * sin(t2 * 3) - &
				12.47 * sin(t2 * 4) - &
			    430.847 * cos(t2) + 12.5024 * cos(t2 * 2) + &
				18.25 * 	cos(t2 * 3)
	sun_lon = 180. - 360. * (utc_sec + eqn_of_time) / 86400.
	sunx = cos(declin * pio180) * cos(sun_lon * pio180)
	suny = cos(declin * pio180) * sin(sun_lon * pio180)
	sunz = sin(declin * pio180)

	wnx = cos(lat_matrix * pio180) * cos(lon_matrix * pio180)
	wny = cos(lat_matrix * pio180) * sin(lon_matrix * pio180)
	wnz = sin(lat_matrix * pio180)

	cosz = wnx * sunx + wny * suny + wnz * sunz
	rshort_clear = max(0., solar_cons * solfac * cosz)

	! circulate over the map
	direct_frac = 0.
	do i_x = 1, x_num
		do i_y = 1, y_num
			if (rshort_clear(i_x,i_y) > 0.) then
				direct_frac = min(0.9, &
				max(0.1, out_buff(i_time,i_x,i_y) / rshort_clear(i_x,i_y)))
			else
				direct_frac = 0.9
			endif

            ! for debug purposes
!            if (out_buff(i_time,i_x,i_y) > solar_cons) then
!                print*,'i_x,i_y',i_x,i_y
!                print*,'i_time',i_time
!                print*,'doy',doy
!            endif

			nbdsf_buff(i_time,i_x,i_y) = out_buff(i_time,i_x,i_y) * direct_frac  * 0.55
			nddsf_buff(i_time,i_x,i_y) = out_buff(i_time,i_x,i_y) * (1 - direct_frac)  * 0.55
			vbdsf_buff(i_time,i_x,i_y) = out_buff(i_time,i_x,i_y) * direct_frac  * 0.45
			vddsf_buff(i_time,i_x,i_y) = out_buff(i_time,i_x,i_y) * (1 - direct_frac)  * 0.45
		enddo
	enddo



enddo

deallocate(rshort_clear)
deallocate(cosz)
deallocate(wnz)
deallocate(wny)
deallocate(wnx)

end subroutine split_radiation
