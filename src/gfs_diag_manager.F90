module gfs_diag_manager_mod

    use fms, only: handle_error=>mpp_error, FATAL, WARNING, NOTE, &
                    mpp_pe, mpp_root_pe, fms_init, mpp_gather, mpp_alltoall, &
                    mpp_get_compute_domain, diag_axis_add_attribute
    use fms, only: time_type, print_date, set_date, set_time, print_time
    use fms, only: diag_manager_init, diag_manager_end, diag_send_complete
    use fms, only: diag_axis_init, register_static_field, diag_manager_set_time_end
    use fms, only: register_diag_field, send_data, assignment(=)
    use fms, only: set_calendar_type, GREGORIAN, set_date, time_type, date_to_string
    use fms, only: domain2d, mpp_define_domains, mpp_npes, mpp_sync, mpp_broadcast
    
    use mpp_io_mod, only: mpp_open, mpp_close, mpp_write_meta, MPP_OVERWR, MPP_SINGLE, &
                          axistype, MPP_NETCDF, mpp_write, MPP_ASCII

    use constants_mod,    only: RAD_TO_DEG, PI

    USE ESMF_Mod, only: ESMF_Time, ESMF_TimeGet, ESMF_Clock, ESMF_ClockGet

    use mpi_def, only: icolor, mc_comp, mc_io
    use layout1, only : ipt_lats_node_r
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

    integer :: lonr, lats_node_r, latr
    integer :: lon_id, lat_id

    integer, allocatable :: nlevs(:,:)
    type(time_type) :: currtime, time_step, starttime

    logical :: diag_active = .false., gfs_io_pe=.true.
  
    public :: gfs_diag_manager_init, gfs_diag_manager_end, register_var, &
             gfs_send_data, register_static, gfs_diag_send_complete, set_gfs_diag_manager_time

   interface gfs_send_data
      module procedure update_opdata_2d_o
      module procedure update_opdata_3d_o
      module procedure update_opdata_2d
      module procedure update_opdata_3d
   end interface gfs_send_data
  
  
   contains

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


    subroutine set_gfs_diag_manager_time(clock)
        TYPE(ESMF_Clock), intent(in) :: clock

        TYPE(ESMF_Time) :: Time
        integer :: itm(6)

        CALL ESMF_ClockGet(clock, currTime = Time)
        CALL ESMF_TimeGet (Time, yy=itm(1), mm=itm(2), dd=itm(3), h=itm(4), m=itm(5), s=itm(6))
        currtime = set_date(itm(1),itm(2),itm(3),itm(4),itm(5),itm(6))

        diag_active = .true. 
    end subroutine set_gfs_diag_manager_time

    subroutine gfs_diag_manager_init(xlat, xlon, global_lats_r, lonsperlar, clock, dt_sec)
        real, intent(in) :: xlat(:,:), xlon(:,:)
        integer, intent(in) :: global_lats_r(:), lonsperlar(:)
        integer, intent(in) :: dt_sec
        TYPE(ESMF_Clock), intent(in) :: clock

        real :: xlonf(size(xlon,1)), dlon

        integer :: i, lev_id, id, k, itm(6)
        integer :: tmp(size(xlat,2)+1), nk, ounit, lonsperlat(size(xlat,2))
        Character (len=32) :: tmpc
        logical :: used
        type(axistype) :: ak_axis, bk_axis
        real :: rtmp(size(ak5,1))
        TYPE(ESMF_Time) :: Time
        integer, allocatable :: nlat_in_pes(:) 
        integer :: nlat_this_pe(1), js, je, is, ie
        real, allocatable :: xlatf(:)

        if (icolor/=2) then 
          call fms_init(localcomm=mc_comp, alt_input_nml_path='gfs_input.nml')
          gfs_io_pe = .false.
        else
          call fms_init(localcomm=mc_io, alt_input_nml_path='gfs_input.nml')
          return
        endif

        call set_calendar_type(GREGORIAN)

        lonr = size(xlon,1)
        lats_node_r = size(xlat,2)
        latr = size(global_lats_r,1)

        nlat_this_pe = lats_node_r

        allocate(nlat_in_pes(mpp_npes())) 

        call mpp_gather(nlat_this_pe,nlat_in_pes)
        call mpp_broadcast(nlat_in_pes,size(nlat_in_pes,1),mpp_root_pe())

        call mpp_define_domains([1,lonr,1,latr], [1, mpp_npes()], domain, yextent=nlat_in_pes)

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

        dlon = 2. * PI / lonr
        xlonf = 0.

        do i = 2, lonr
           xlonf(i) = xlonf(i-1) + dlon
        end do

        allocate(xlatf(latr))
        xlatf = -1000.0

        call mpp_get_compute_domain(domain, is, ie, js, je)

        xlatf(js:je) = xlat(1,:)

        call diag_manager_init()

        tmp(1)=latr
        total_eles = 0
        do i=1,lats_node_r
           tmp(i+1)=global_lats_r(ipt_lats_node_r-1+i)
           total_eles = total_eles + lonsperlar(tmp(i+1))
           lonsperlat(i) = lonsperlar(tmp(i+1))
        enddo

        lon_id = diag_axis_init(name='lon', data=xlonf(:)*RAD_TO_DEG, &
                   units='degrees_east' , cart_name='X', long_name='longitude', domain2=domain)
        call diag_axis_add_attribute(lon_id, 'lonsperlat', lonsperlat)
        lat_id = diag_axis_init(name='lat', data=xlatf(:)*RAD_TO_DEG, &
                   units='degrees_north' , cart_name='Y', long_name='latitude', domain2=domain)
        call diag_axis_add_attribute(lat_id, 'decomp_gfs', tmp)

        allocate(nlevs(2,max_nlevs))


        call set_gfs_diag_manager_time(clock)
        diag_active = .false. ! this is needed because of initial filtering in gfs
        starttime = currtime

        CALL ESMF_ClockGet(clock, stopTime = Time)
        CALL ESMF_TimeGet (Time, yy=itm(1), mm=itm(2), dd=itm(3), h=itm(4), m=itm(5), s=itm(6))
        call diag_manager_set_time_end(set_date(itm(1), itm(2), itm(3), itm(4), itm(5), itm(6)))

        time_step = set_time(dt_sec)

        deallocate(xlatf, nlat_in_pes)

    end subroutine gfs_diag_manager_init


    subroutine gfs_diag_send_complete()
        if (gfs_io_pe) return
       call diag_send_complete(time_step)
    end subroutine gfs_diag_send_complete

    integer function register_var(name, long_name, units, range, standard_name, levs)
        character (len=*), intent(in) :: name
        character (len=*), intent(in), optional :: long_name, units
        real, intent(in), optional :: range(2)
        character (len=*), intent(in), optional :: standard_name
        integer, intent(in), optional :: levs

        character (len=32) :: tmp
        integer :: n, i, lev_id

        if (gfs_io_pe) return
        if (.not.present(levs)) then
            register_var = register_diag_field('gfs', trim(name), (/lon_id, lat_id/), starttime, &
                long_name=long_name, units=units, range=range, standard_name=standard_name)
        else
            lev_id = get_level_id(levs)
            register_var = register_diag_field('gfs', trim(name), (/lon_id, lat_id, lev_id/), starttime, &
                long_name=long_name, units=units, range=range, standard_name=standard_name)
        endif

    end function register_var

    integer function register_static(name, long_name, units, range, standard_name, levs)
      character (len=*), intent(in) :: name
      character (len=*), intent(in), optional :: long_name, units
      real, intent(in), optional :: range(2)
      character (len=*), intent(in), optional :: standard_name
      integer, intent(in), optional :: levs

      integer :: n, i, lev_id

      if (gfs_io_pe) return
      if (.not.present(levs)) then
         register_static = register_static_field('gfs', trim(name), (/lon_id, lat_id/), &
              long_name=long_name, units=units, range=range, standard_name=standard_name)
      else
         lev_id = get_level_id(levs)
         register_static = register_static_field('gfs', trim(name), (/lon_id, lat_id, lev_id/), &
              long_name=long_name, units=units, range=range, standard_name=standard_name)
      endif

    end function register_static


    subroutine update_opdata_2d_o(field_id, field, istrt, im, lan)
      integer, intent(in) :: field_id
      real, intent(in) :: field(:)
      integer, intent(in) :: lan, istrt, im

      !local
      real :: field1(im,1)
      real :: wgt1(im,1)
      logical :: mask(im,1)
      logical :: used
      integer :: is, ie 

        if (gfs_io_pe) return
        if (.not. diag_active) return
      if (field_id<=0) return

      is = istrt
      ie = is+im-1

      field1(1:im,1)=field(:)

      used = send_data(field_id, field1(:,:), currtime, js_in=lan, je_in=lan, is_in=is, ie_in=ie)
    end subroutine update_opdata_2d_o

    subroutine update_opdata_3d_o(field_id, field, istrt, im, lan)
      integer, intent(in) :: field_id
      real, intent(in) :: field(:,:)
      integer, intent(in) :: lan, im, istrt
      !local
      real :: field2(im,1,size(field,2))
      real :: wgt1(im,1,size(field,2))
      logical :: mask(im,1,size(field,2))
      logical :: used
      integer :: k, is, ie

        if (gfs_io_pe) return
        if (.not. diag_active) return
      if (field_id<=0) return

      is = istrt
      ie = is+im-1

      field2(1:im,1,:) = field(:,:)

      used = send_data(field_id, field2(:,:,:), currtime, js_in=lan, je_in=lan, is_in=is, ie_in=ie)
    end subroutine update_opdata_3d_o

    subroutine update_opdata_2d(field, id, name, static, long_name, units, range, standard_name)
        real, intent(in) :: field(:,:)
        integer, intent(in), optional :: id
        character (len=*), intent(in), optional :: name
        character (len=*), intent(in), optional :: long_name, units
        real, intent(in), optional :: range(2)
        character (len=*), intent(in), optional :: standard_name
        logical, intent(in), optional :: static

        !local
        logical :: used
        integer :: field_id

        if (gfs_io_pe) return
        if (.not. diag_active) return

        if (present(id).and.(id<=0)) return
        if (.not. present(id) .and. .not. present(name) ) call handle_error(fatal, 'gfs_diag_manager: both id and name not present for field')

        if (present(id).and.(id>0)) then
            used = send_data(id, field(:,:), currtime)
            return
        endif

        field_id = get_field_id(name, static=static, long_name=long_name, units=units, &
                        range=range, standard_name=standard_name)

        used = send_data(field_id, field(:,:), currtime)

    end subroutine update_opdata_2d


    subroutine update_opdata_3d(field, id, name, static, long_name, units, range, standard_name)
        real, intent(in) :: field(:,:,:)
        integer, intent(in), optional :: id
        character (len=*), intent(in), optional:: name
        character (len=*), intent(in), optional :: long_name, units
        real, intent(in), optional :: range(2)
        character (len=*), intent(in), optional :: standard_name
        logical, intent(in), optional :: static

        !local
        logical :: used
        integer :: field_id
        integer :: levs

        if (gfs_io_pe) return
        if (.not. diag_active) return

        if (present(id).and.(id<=0)) return
        if (.not. present(id) .and. .not. present(name) ) call handle_error(fatal, 'gfs_diag_manager: both id and name not present for field')

        if (present(id).and.(id>0)) then
            used = send_data(id, field(:,:,:), currtime)
            return
        endif

        levs = size(field,3)
        field_id = get_field_id(name, static=static, long_name=long_name, units=units, &
                        range=range, standard_name=standard_name, levs=levs)
        used = send_data(id, field(:,:,:), currtime)

    end subroutine update_opdata_3d

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

        n_fields=n_fields+1
        if (n_fields > max_fields) call handle_error(fatal, 'gfs_diag_manager_mod: n_fields>max_fields')
        if (present(static) .and. static) then
            field_ids(n_fields)%id = register_static(name, long_name=long_name, units=units, &
                        range=range, standard_name=standard_name, levs=levs)
        else
            field_ids(n_fields)%id = register_var(name, long_name=long_name, units=units, &
                        range=range, standard_name=standard_name, levs=levs)
        endif
        field_ids(n_fields)%name = name
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
          if (n_nlevs > max_nlevs) call handle_error(fatal, 'register_var: n_nlevs>max_nlevs, increase max_nlevs!')
          nlevs(1,n_nlevs) = levs
          write(tmp,*) levs
          tmp=trim(adjustl(tmp))
          nlevs(2,n_nlevs) = diag_axis_init(name='lev'//trim(tmp), data=(/(real(i),i=1,levs)/), units='', &
                                      cart_name='Z', long_name='levels')
          get_level_id = nlevs(2,n_nlevs)

    end function get_level_id

    subroutine gfs_diag_manager_end()
      if (gfs_io_pe) return
      call diag_manager_end(currtime)
    end subroutine gfs_diag_manager_end

end module gfs_diag_manager_mod
