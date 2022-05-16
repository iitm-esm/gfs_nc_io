

    subroutine UPDATE_OPDATA_ (field, id, name, static, long_name, units, range, standard_name)
        real, intent(in), DIMENSION_ :: field
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
            used = send_data(id, field, currtime)
            return
        endif

#ifdef LEVS_
            levs = size(field, LEVS_ )
            field_id = get_field_id(name, static=static, long_name=long_name, units=units, &
                            range=range, standard_name=standard_name, levs=levs)
#else
            field_id = get_field_id(name, static=static, long_name=long_name, units=units, &
                            range=range, standard_name=standard_name)
#endif

        used = send_data(field_id, field, currtime)

    end subroutine UPDATE_OPDATA_



    subroutine O_UPDATE_OPDATA_ (field, istrt, im, lan, id, name, static, long_name, units, range, standard_name)
        real, intent(in), O_DIMENSION_ :: field
        integer, intent(in) :: istrt, im, lan
        integer, intent(in), optional :: id
        character (len=*), intent(in), optional :: name
        character (len=*), intent(in), optional :: long_name, units
        real, intent(in), optional :: range(2)
        character (len=*), intent(in), optional :: standard_name
        logical, intent(in), optional :: static

        !local
        logical :: used
        integer :: field_id, is, ie, js, je, levs

#ifdef LEVS_
        real :: field1(1:im,1,size(field,2))
#else
        real :: field1(1:im,1)
#endif

        if (gfs_io_pe) return
        if (.not. diag_active) return


        is = istrt; ie = istrt + im -1
        js = lan; je = lan 

        if (present(id).and.(id<=0)) return
        if (.not. present(id) .and. .not. present(name) ) call handle_error(fatal, 'gfs_diag_manager: both id and name not present for field')

        if (present(id).and.(id>0)) then
            used = send_data(field_id, field1, currtime, is_in=is, js_in=js)
            return
        endif


#ifdef LEVS_
            field1(:,1,:) = field(:,:)
            levs = size(field1, LEVS_ )
            field_id = get_field_id(name, static=static, long_name=long_name, units=units, &
                            range=range, standard_name=standard_name, levs=levs)
#else
            field1(:,1) = field(:)
            field_id = get_field_id(name, static=static, long_name=long_name, units=units, &
                            range=range, standard_name=standard_name)
#endif
        used = send_data(field_id, field1, currtime, is_in=is, js_in=js)

    end subroutine O_UPDATE_OPDATA_