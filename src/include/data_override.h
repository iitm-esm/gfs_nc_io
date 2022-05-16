
    subroutine DATA_OVERRIDE_(field_name, field, override)
        character(len=*), intent(in) :: field_name
        real, DIMENSION_ :: field
        logical, intent(out), optional :: override

        if (gfs_io_pe) return
        call data_override('GFS', field_name, field, time=currtime, override=override)

    end subroutine DATA_OVERRIDE_

