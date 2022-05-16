

#undef UPDATE_OPDATA_ 
#define UPDATE_OPDATA_ update_opdata_2d
#undef O_UPDATE_OPDATA_
#define O_UPDATE_OPDATA_ update_opdata_2d_o
#undef DIMENSION_
#define DIMENSION_ dimension(:,:)
#undef O_DIMENSION_
#define O_DIMENSION_ dimension(:)
#undef LEVS_
#include <update_opdata.h>

#undef UPDATE_OPDATA_ 
#define UPDATE_OPDATA_ update_opdata_3d
#undef DIMENSION_
#define DIMENSION_ dimension(:,:,:)
#undef O_UPDATE_OPDATA_
#define O_UPDATE_OPDATA_ update_opdata_3d_o
#undef O_DIMENSION_
#define O_DIMENSION_ dimension(:,:)
#undef LEVS_
#define LEVS_ 3
#include <update_opdata.h>