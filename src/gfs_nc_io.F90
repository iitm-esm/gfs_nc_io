module gfs_nc_io_mod


    use fms, only: handle_error=>mpp_error, FATAL, WARNING, NOTE, &
            mpp_pe, mpp_root_pe, fms_init, mpp_gather, mpp_alltoall, &
            mpp_get_compute_domain, diag_axis_add_attribute, mpp_max, &
            time_type, print_date, set_date, set_time, print_time, &
            diag_manager_init, diag_manager_end, diag_send_complete, &
            diag_axis_init, register_static_field, diag_manager_set_time_end, &
            register_diag_field, send_data, assignment(=), data_override_init, &
            set_calendar_type, GREGORIAN, set_date, time_type, date_to_string, &
            domain2d, mpp_define_domains, mpp_npes, mpp_sync, mpp_broadcast, &
            data_override
    
    use mpp_io_mod, only: mpp_open, mpp_close, mpp_write_meta, MPP_OVERWR, MPP_SINGLE, &
            axistype, MPP_NETCDF, mpp_write, MPP_ASCII

    use constants_mod,    only: RAD_TO_DEG, PI

    USE ESMF_Mod, only: ESMF_Time, ESMF_TimeGet, ESMF_Clock, ESMF_ClockGet

    use mpi_def, only: icolor, mc_comp, mc_io
    use layout1, only : ipt_lats_node_r, me 
    use coordinate_def, only: ak5, bk5  

    implicit none

    private

    integer, parameter :: max_fields = 100

    type(domain2d) :: domain

    type field_id_type
        integer :: id = 0
        character(len=32) :: name = ''
    end type field_id_type

    type(field_id_type), dimension(max_fields) :: field_ids

    integer :: max_nlevs=100, n_nlevs=0, total_eles=0, n_fields=0

    integer :: nlonGlobal, nlatLocal, nlatGlobal
    integer :: lon_id, lat_id

    integer, allocatable :: nlevs(:,:)
    type(time_type) :: currtime, time_step, starttime

    logical :: diag_active = .false., gfs_io_pe=.true.
  
    public :: gfs_nc_io_init, gfs_nc_io_end, gfs_register_diag_field, gfs_send_data, &
                gfs_register_static_field, gfs_diag_send_complete, set_gfs_nc_io_time, &
                gfs_data_override

    interface gfs_send_data
       module procedure update_opdata_2d_o, update_opdata_3d_o, update_opdata_2d, &
                    update_opdata_3d
    end interface gfs_send_data

    interface gfs_data_override
        module procedure gfs_data_override0d, gfs_data_override2d, gfs_data_override3d, &
                    gfs_data_override2d_o, gfs_data_override3d_o
    end interface gfs_data_override
  
  
    contains

    subroutine set_domain()
        integer, allocatable :: nlat_in_pes(:)

        allocate(nlat_in_pes(mpp_npes())) 

        call mpp_gather([nlatLocal],nlat_in_pes)
        call mpp_broadcast(nlat_in_pes,size(nlat_in_pes,1),mpp_root_pe())

        call mpp_define_domains([1,nlonGlobal,1,nlatGlobal], [1, mpp_npes()], domain, yextent=nlat_in_pes)

        deallocate(nlat_in_pes)
    end subroutine set_domain


    subroutine write_diag_post_nml(startdate, enddate, deltim, calendar_type)
      integer, intent(in) :: startdate(6), enddate(6), deltim, calendar_type

      integer :: ounit

      namelist/diag_post_nml/ startdate, calendar_type, deltim, enddate

      if (mpp_pe()==mpp_root_pe()) then
         call mpp_open(ounit,'diag_post.nml',action=MPP_OVERWR, &
               form=MPP_ASCII,threading=MPP_SINGLE)
         write(ounit,nml=diag_post_nml)
         call mpp_close(ounit)
      endif

    end subroutine write_diag_post_nml


    subroutine set_gfs_nc_io_time(clock)
        TYPE(ESMF_Clock), intent(in) :: clock

        TYPE(ESMF_Time) :: Time
        integer :: itm(6)

        if (gfs_io_pe) return 
        
        CALL ESMF_ClockGet(clock, currTime = Time)
        CALL ESMF_TimeGet (Time, yy=itm(1), mm=itm(2), dd=itm(3), h=itm(4), m=itm(5), s=itm(6))
        currtime = set_date(itm(1),itm(2),itm(3),itm(4),itm(5),itm(6))

        diag_active = .true. 
    end subroutine set_gfs_nc_io_time


    subroutine gfs_nc_io_init(xlat, xlon, global_lats_r, lonsperlar, clock, dt_sec)
        real, intent(in) :: xlat(:,:), xlon(:,:)
        integer, intent(in) :: global_lats_r(:), lonsperlar(:)
        integer, intent(in) :: dt_sec
        TYPE(ESMF_Clock), intent(in) :: clock

        real :: xlonf(size(xlon,1)), dlon

        integer :: i, lev_id, id, k, itm(6), j
        integer :: latInd(size(xlat,2)), nk, ounit, lonsperlat(size(xlat,2))
        Character (len=32) :: tmpc
        logical :: used
        type(axistype) :: ak_axis, bk_axis
        real :: rtmp(size(ak5,1))
        TYPE(ESMF_Time) :: Time
        integer :: js, je, is, ie, max_nlat
        real, allocatable :: xlatf(:)
        integer, dimension(size(global_lats_r,1)) :: latIndGlobal, lonsperlatGlobal
        integer, allocatable :: itmpLocal(:), itmp(:)

        if (icolor/=2) then 
          call fms_init(localcomm=mc_comp, alt_input_nml_path='gfs_input.nml')
          gfs_io_pe = .false.
        else
          call fms_init(localcomm=mc_io, alt_input_nml_path='gfs_input.nml')
          return
        endif

        call set_calendar_type(GREGORIAN)

        nlonGlobal = size(xlon,1)
        nlatLocal = size(xlat,2)
        nlatGlobal = size(global_lats_r,1)

        call set_domain

        nk = size(ak5,1)
        if (mpp_pe()==mpp_root_pe()) then
           call mpp_open(ounit,'ak_bk_out.nc',action=MPP_OVERWR, &
                 form=MPP_NETCDF,threading=MPP_SINGLE)
           rtmp = [(ak5(k),k=nk,1,-1)] * 1000. ! cb -> Pascal
           call mpp_write_meta(ounit, ak_axis, 'ak', 'Pascal', 'coef_a', data=rtmp)
           rtmp = [(bk5(k),k=nk,1,-1)]
           call mpp_write_meta(ounit, bk_axis, 'bk', '1', 'coef_b', data=rtmp)
           call mpp_write(ounit,ak_axis)
           call mpp_write(ounit,bk_axis)
           call mpp_close(ounit)
        endif

        dlon = 2. * PI / nlonGlobal
        xlonf = 0.

        do i = 2, nlonGlobal
           xlonf(i) = xlonf(i-1) + dlon
        end do

        allocate(xlatf(nlatGlobal))
        xlatf = -1000.0

        call mpp_get_compute_domain(domain, is, ie, js, je)

        xlatf(js:je) = xlat(1,:)

        call diag_manager_init()

        total_eles = 0
        do i=1,nlatLocal
           latInd(i)=global_lats_r(ipt_lats_node_r-1+i)
           total_eles = total_eles + lonsperlar(latInd(i))
           lonsperlat(i) = lonsperlar(latInd(i))
        enddo

        lon_id = diag_axis_init(name='lon', data=xlonf(:)*RAD_TO_DEG, &
                   units='degrees_east' , cart_name='X', long_name='longitude', domain2=domain)
        call diag_axis_add_attribute(lon_id, 'lonsperlat', lonsperlat)

        lat_id = diag_axis_init(name='lat', data=xlatf(:)*RAD_TO_DEG, &
                   units='degrees_north' , cart_name='Y', long_name='latitude', domain2=domain)
        call diag_axis_add_attribute(lat_id, 'decomp_gfs', latInd)


        allocate(nlevs(2,max_nlevs))

        call set_gfs_nc_io_time(clock)
        diag_active = .false. ! this is needed because of tldfi in gfs
        starttime = currtime

        CALL ESMF_ClockGet(clock, stopTime = Time)
        CALL ESMF_TimeGet (Time, yy=itm(1), mm=itm(2), dd=itm(3), h=itm(4), m=itm(5), s=itm(6))
        call diag_manager_set_time_end(set_date(itm(1), itm(2), itm(3), itm(4), itm(5), itm(6)))

        time_step = set_time(dt_sec)

        deallocate(xlatf)

        call data_override_init(gfs_domain_in=domain, gfs_lon_in=xlon, gfs_lat_in=xlat)

    end subroutine gfs_nc_io_init


    subroutine gfs_diag_send_complete()
        if (gfs_io_pe) return
       call diag_send_complete(time_step)
    end subroutine gfs_diag_send_complete

    integer function gfs_register_diag_field(name, long_name, units, range, standard_name, levs)
        character (len=*), intent(in) :: name
        character (len=*), intent(in), optional :: long_name, units
        real, intent(in), optional :: range(2)
        character (len=*), intent(in), optional :: standard_name
        integer, intent(in), optional :: levs

        character (len=32) :: tmp
        integer :: n, i, lev_id

        if (gfs_io_pe) return
        if (.not.present(levs)) then
            gfs_register_diag_field = register_diag_field('gfs', trim(name), &
                (/lon_id, lat_id/), starttime, long_name=long_name, &
                units=units, range=range, standard_name=standard_name, total_elements=total_eles)
        else
            lev_id = get_level_id(levs)
            gfs_register_diag_field = register_diag_field('gfs', trim(name), &
                (/lon_id, lat_id, lev_id/), starttime, &
                long_name=long_name, units=units, range=range, &
                standard_name=standard_name, total_elements=total_eles*levs)
        endif

    end function gfs_register_diag_field

    integer function gfs_register_static_field(name, long_name, units, range, standard_name, levs)
      character (len=*), intent(in) :: name
      character (len=*), intent(in), optional :: long_name, units
      real, intent(in), optional :: range(2)
      character (len=*), intent(in), optional :: standard_name
      integer, intent(in), optional :: levs

      integer :: n, i, lev_id

      if (gfs_io_pe) return
      if (.not.present(levs)) then
         gfs_register_static_field = register_static_field('gfs', trim(name), (/lon_id, lat_id/), &
              long_name=long_name, units=units, range=range, standard_name=standard_name, total_elements=total_eles)
      else
         lev_id = get_level_id(levs)
         gfs_register_static_field = register_static_field('gfs', trim(name), (/lon_id, lat_id, lev_id/), &
              long_name=long_name, units=units, range=range, standard_name=standard_name, total_elements=total_eles*levs)
      endif

    end function gfs_register_static_field

#include <update_opdata.inc>
#include <data_override.inc>

    function get_field_id(name, static, long_name, units, range, standard_name, levs)
        character (len=*), intent(in) :: name
        character (len=*), intent(in), optional :: long_name, units
        real, intent(in), optional :: range(2)
        character (len=*), intent(in), optional :: standard_name
        integer, intent(in), optional :: levs
        logical, intent(in), optional :: static
        integer :: get_field_id

        integer :: i  

        get_field_id = 0

        do i = 1, n_fields
            if (trim(name)==trim(field_ids(i)%name)) then
                get_field_id = field_ids(i)%id
                return
            endif
        enddo

!$omp single
        n_fields=n_fields+1
        if (n_fields > max_fields) call handle_error(fatal, 'gfs_diag_manager_mod: n_fields>max_fields')
        if (present(static) .and. static) then
            field_ids(n_fields)%id = gfs_register_static_field(name, long_name=long_name, units=units, &
                        range=range, standard_name=standard_name, levs=levs)
        else
            field_ids(n_fields)%id = gfs_register_diag_field(name, long_name=long_name, units=units, &
                        range=range, standard_name=standard_name, levs=levs)
        endif
        field_ids(n_fields)%name = name
!$omp end single
        get_field_id = field_ids(n_fields)%id

    end function get_field_id


    function get_level_id(levs)
        integer, intent(in) :: levs
        integer :: get_level_id

        integer :: i  
        character (len=32) :: tmp

          do i = 1, n_nlevs
             if (levs==nlevs(1,i)) then
                get_level_id=nlevs(2,i)
                return
             endif
          enddo

          n_nlevs=n_nlevs+1
          if (n_nlevs > max_nlevs) call handle_error(fatal, 'gfs_register_diag_field: n_nlevs>max_nlevs, increase max_nlevs!')
          nlevs(1,n_nlevs) = levs
          write(tmp,*) levs
          tmp=trim(adjustl(tmp))
          nlevs(2,n_nlevs) = diag_axis_init(name='lev'//trim(tmp), data=(/(real(i),i=1,levs)/), units='', &
                                      cart_name='Z', long_name='levels')
          get_level_id = nlevs(2,n_nlevs)

    end function get_level_id

    subroutine gfs_nc_io_end()
      if (gfs_io_pe) return
      call diag_manager_end(currtime)
    end subroutine gfs_nc_io_end

end module gfs_nc_io_mod