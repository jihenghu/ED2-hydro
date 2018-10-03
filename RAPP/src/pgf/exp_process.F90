!==========================================================================================!
!==========================================================================================!
!  Subroutine exp_process()                                                    !
!                                                                                          !
!     Currently This subroutine will read the pgf dataset in netcdf format, process the data!
! calculate VPD, and write them into .h5 file that ED2 model can read					   !
!------------------------------------------------------------------------------------------!
subroutine exp_process(year)
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
	!----------------- local parameters
	real,		parameter		:: 		pgf_res = 0.5 ! deg resolution
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
	integer							::	  x_num,y_num,t_num, i_x, i_y, i_t
	integer							::	  x_start,y_start,t_start
	integer							::	  x_start_1,x_num_1,x_num_2, x_start_2 		!in case the region of interest spans 0 or 360
	real							::    lon_start, lat_start
    logical                                     :: file_is_there
    integer                                     :: ierr
	integer										:: dim_len
	real, pointer, dimension(:,:,:)	::    tas_buff, shum_buff, pres_buff, vpd_buff
	real, pointer, dimension(:,:,:)	::    temp_buff					! in case the region of interest spans 0 or 360
	real, pointer, dimension(:,:,:)	::    out_buff
	real, pointer, dimension(:,:)	::	  lon_matrix, lat_matrix
	real, pointer, dimension(:)		:: 	  lat_buff
   integer, dimension(2)                  :: grid_dims
   integer, dimension(3)                  :: var_dims

	!--------------------------------------------------------
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
	x_num = (nint((lone - lonw) / pgf_res) + 1)  ! longitude
	x_num_1 = -1
	x_start_1 = -1
	x_num_2 = -1
	x_start_2 = -1


	y_num = nint((latn  - lats) / pgf_res) + 1
	y_start = nint((lats + 89.75) / pgf_res) + 1
	lat_start = (y_start - 1) * pgf_res - 89.75


	! consider different conditions
	if (lone >= 0. .and. lonw >= 0.) then
	! 1. the user used 0-360
	!    same as the pgf, use directly
		x_start = nint((lonw - 0.25) / pgf_res) + 1
		lon_start = (x_start - 1) * pgf_res + 0.25
		
	else	
	! 2. the user used -180 ~ 180
	!    at least lonw is < 0.
		x_start = nint((359.75 + lonw) / pgf_res) + 1
		lon_start = (x_start - 1) * pgf_res + 0.25 - 360.
		
		if (lone < 0.) then
		! 2.1 lone is also < 0.
		! x_start and x_num are enough to extract data
		! no need to do anything here
		else
		! 2.2 lone >=0, corssed the 360 line
			x_start_2 = 1
			x_num_2 = (nint((lone - 0.25) / pgf_res) + 1)
			x_num_1 = (nint((0.25 - lonw) / pgf_res))
			x_start_1 = x_start

		endif
			

	endif

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
		

		allocate(tas_buff(x_num,y_num,t_num))
		allocate(shum_buff(x_num,y_num,t_num))
		allocate(pres_buff(x_num,y_num,t_num))
		allocate(vpd_buff(x_num,y_num,t_num))

		! create output file
   		write(output_name,fmt='(a,i4.4,a,a)')                                                   &
                      trim(outpref),year,trim(month_str(i_month)),'.h5'
		
   		write (unit=*,fmt='(a,1x,2a)') '         [|] Writing file:',trim(output_name),'...'

	   call shdf5_open_f(trim(output_name),'W',1)

   		! write lat and lon
		grid_dims = (/x_num, y_num/)
		call shdf5_orec_f( 2, grid_dims, trim('lon'),rvara=lon_matrix)
		call shdf5_orec_f( 2, grid_dims, trim('lat'),rvara=lat_matrix)



		! loop over pressure, humidity, temperature
			
		do i_var = 1, 3
			select case (i_var)
				case (1)   ! pres
					var_name = 'pres'
					output_var_name = 'pres'
				case (2)   ! shum
					var_name = 'shum'
					output_var_name = 'sh'
				case (3)   ! tas
					var_name = 'tas'
					output_var_name = 'tmp'
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
			if (x_start_1 < 0) then
			! the simple case
				select case (i_var)
					case (1)   ! pres
				        ierr = nf90_get_var(ncid,var_id,pres_buff,start = (/x_start,y_start,t_start/),count= (/x_num,y_num,t_num/))
					case (2)   ! shum
			    	    ierr = nf90_get_var(ncid,var_id,shum_buff,start = (/x_start,y_start,t_start/),count= (/x_num,y_num,t_num/))
					case (3)   ! tas
				        ierr = nf90_get_var(ncid,var_id,tas_buff,start = (/x_start,y_start,t_start/),count= (/x_num,y_num,t_num/))
				end select
			else
			! the more complext case
			! we need to read twice
				select case (i_var)
					case (1)   ! pres
!						print*,'x_start_1',x_start_1,'x_start_2',x_start_2,'x_num_1',x_num_1,'x_num_2',x_num_2,&
!								'x_start',x_start,'x_num',x_num
						allocate(temp_buff(x_num_1,y_num,t_num))
				        ierr = nf90_get_var(ncid,var_id,temp_buff,start = (/x_start_1,y_start,t_start/),count= (/x_num_1,y_num,t_num/))
						pres_buff(1:x_num_1,1:y_num,1:t_num) = temp_buff
						deallocate(temp_buff)

						allocate(temp_buff(x_num_2,y_num,t_num))
				        ierr = nf90_get_var(ncid,var_id,temp_buff,start = (/x_start_2,y_start,t_start/),count= (/x_num_2,y_num,t_num/))
						pres_buff((x_num_1+1):x_num,1:y_num,1:t_num) = temp_buff
!						print*,temp_buff
						deallocate(temp_buff)
					case (2)   ! shum
						allocate(temp_buff(x_num_1,y_num,t_num))
				        ierr = nf90_get_var(ncid,var_id,temp_buff,start = (/x_start_1,y_start,t_start/),count= (/x_num_1,y_num,t_num/))
						shum_buff(1:x_num_1,1:y_num,1:t_num) = temp_buff
						deallocate(temp_buff)

						allocate(temp_buff(x_num_2,y_num,t_num))
				        ierr = nf90_get_var(ncid,var_id,temp_buff,start = (/x_start_2,y_start,t_start/),count= (/x_num_2,y_num,t_num/))
						shum_buff((x_num_1+1):x_num,1:y_num,1:t_num) = temp_buff
!						print*,temp_buff
						deallocate(temp_buff)
					case (3)   ! tas
						allocate(temp_buff(x_num_1,y_num,t_num))
				        ierr = nf90_get_var(ncid,var_id,temp_buff,start = (/x_start_1,y_start,t_start/),count= (/x_num_1,y_num,t_num/))
						tas_buff(1:x_num_1,1:y_num,1:t_num) = temp_buff
						deallocate(temp_buff)

						allocate(temp_buff(x_num_2,y_num,t_num))
				        ierr = nf90_get_var(ncid,var_id,temp_buff,start = (/x_start_2,y_start,t_start/),count= (/x_num_2,y_num,t_num/))
						tas_buff((x_num_1+1):x_num,1:y_num,1:t_num) = temp_buff
!						print*,temp_buff
						deallocate(temp_buff)
				end select

			endif

		    ! close the file
			ierr = nf90_close(ncid)
		    if (ierr /= NF90_NOERR) then
        	 call ncdf_load_err(ierr)
				 print*,'Error when closing the file:'
				 print*,file_name
    		end if
			!print*,var_name,data_buff(1,:,1)

		
		enddo

			!------------------------  data process and write  ------------!
			!calculate VPD first
			do i_x = 1, x_num
				do i_y = 1, y_num
					do i_t = 1, t_num
						! Here we use August-Roche-Magnus approximation for
						! saturated water vapor pressure
						vpd_buff(i_x,i_y,i_t) = &
							610.94 * exp(17.625 * (tas_buff(i_x,i_y,i_t) - 273.15) / &
							(tas_buff(i_x,i_y,i_t) - 273.15 + 243.04)) - &  ! in Pa
							shum_buff(i_x,i_y,i_t) / (0.622 + 0.378 * shum_buff(i_x,i_y,i_t)) * &
							pres_buff(i_x,i_y,i_t)
					enddo
				enddo
			enddo

		! write the output
		allocate(out_buff(t_num,x_num,y_num))
		var_dims = (/t_num, x_num, y_num/)

				

		! rewrite the data_buff to output_buff because in h5 file time
		! dimension comes first and lat dimension starts from the largest
!   		call swap_dim(x_num,y_num,t_num,x_num,y_num,0.,tas_buff,out_buff)
!		call shdf5_orec_f( 3, var_dims, trim('tas'),rvara=out_buff)

 !  		call swap_dim(x_num,y_num,t_num,x_num,y_num,0.,pres_buff,out_buff)
!		call shdf5_orec_f( 3, var_dims, trim('pres'),rvara=out_buff)

 !  		call swap_dim(x_num,y_num,t_num,x_num,y_num,0.,shum_buff,out_buff)
!		call shdf5_orec_f( 3, var_dims, trim('shum'),rvara=out_buff)

   		call swap_dim(x_num,y_num,t_num,x_num,y_num,0.,vpd_buff,out_buff)
		call shdf5_orec_f( 3, var_dims, trim('vpd'),rvara=out_buff)

		deallocate(tas_buff)
		deallocate(pres_buff)
		deallocate(shum_buff)
		deallocate(vpd_buff)
		deallocate(out_buff)

		call shdf5_close_f()

	enddo month_loop

		deallocate(lon_matrix)
		deallocate(lat_matrix)

contains
logical function isLeap(year)
	implicit none
	integer,   intent(in)	:: year
	isLeap = (mod(year,4) == 0 .and. mod(year,100) > 0) .or. &
			 (mod(year,400) == 0 .and. mod(year,100) == 0)
end function isLeap

end subroutine exp_process



