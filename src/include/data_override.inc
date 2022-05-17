
#undef DATA_OVERRIDE_
#define DATA_OVERRIDE_ gfs_data_override2d
#undef DIMENSION_
#define DIMENSION_ intent(inout), dimension(:,:)
#include <data_override.h>

#undef DATA_OVERRIDE_
#define DATA_OVERRIDE_ gfs_data_override3d
#undef DIMENSION_
#define DIMENSION_ intent(inout), dimension(:,:,:)
#include <data_override.h>

#undef DATA_OVERRIDE_
#define DATA_OVERRIDE_ gfs_data_override0d
#undef DIMENSION_
#define DIMENSION_ intent(out)
#include <data_override.h>

#undef O_DATA_OVERRIDE_
#define O_DATA_OVERRIDE_ gfs_data_override2d_o
#undef O_DIMENSION_
#define O_DIMENSION_ dimension(:)
#undef LEVS_
#include <data_override_o.h>

#undef O_DATA_OVERRIDE_
#define O_DATA_OVERRIDE_ gfs_data_override3d_o
#undef O_DIMENSION_
#define O_DIMENSION_ dimension(:,:)
#undef LEVS_
#define LEVS_ 3
#include <data_override_o.h>