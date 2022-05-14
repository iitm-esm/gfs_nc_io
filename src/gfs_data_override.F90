module gfs_data_override_mod
    use fms, only: fms_init 
    use fms, only: data_override_init, data_override
    use fms, only: set_calendar_type, GREGORIAN, set_date, time_type
    use mpi_def, only: icolor, mc_comp, mc_io
    USE ESMF_Mod, only: ESMF_Time, ESMF_TimeGet, ESMF_Clock, ESMF_ClockGet

    implicit none

    type(time_type) :: data_override_time

    interface gfs_data_override
        module procedure gfs_data_override2d
        module procedure gfs_data_override3d
    end interface gfs_data_override

    contains

    subroutine gfs_data_override_init(lon_in, lat_in, clock)
        real, intent(in), dimension(:,:) :: lon_in, lat_in !> in radians
        TYPE(ESMF_Clock), intent(in) :: clock

        if (icolor/=2) then 
          call fms_init(localcomm=mc_comp, alt_input_nml_path='gfs_input.nml')
        else
          call fms_init(localcomm=mc_io, alt_input_nml_path='gfs_input.nml')
        endif

        call set_calendar_type(GREGORIAN)
        call data_override_init(gfs_lon_in=lon_in, gfs_lat_in=lat_in)
        call set_gfs_data_override_time(clock) 

    end subroutine gfs_data_override_init
    

    subroutine set_gfs_data_override_time(clock)
        TYPE(ESMF_Clock), intent(in) :: clock
        TYPE(ESMF_Time) :: Time
        integer :: itm(6)
        CALL ESMF_ClockGet(clock, currTime = Time)
        CALL ESMF_TimeGet (Time, yy=itm(1), mm=itm(2), dd=itm(3), h=itm(4), m=itm(5), s=itm(6))
        data_override_time = set_date(itm(1),itm(2),itm(3),itm(4),itm(5),itm(6))
    end subroutine set_gfs_data_override_time

    subroutine gfs_data_override2d(field_name,field,override)
        character(len=*), intent(in) :: field_name
        real, intent(out), dimension(:,:) :: field
        logical, intent(out), optional :: override

        call data_override('GFS',field_name,field,time=data_override_time,override=override)

    end subroutine gfs_data_override2d

    subroutine gfs_data_override3d(field_name,field,override)
        character(len=*), intent(in) :: field_name
        real, intent(out), dimension(:,:,:) :: field
        logical, intent(out), optional :: override

        call data_override('GFS',field_name,field,time=data_override_time,override=override)
    end subroutine gfs_data_override3d

end module gfs_data_override_mod