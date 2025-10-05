! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
!|  _ _ _   _   _   _      |                                   |
!| |_|_|_| |_| |_| |_|     |  The Great Finite Volume Library  |
!| |_|_    |_|_|_| |_|     |  Author: Ricardo Costa            |
!| |_|_|    \_\_/  |_|_ _  |  Version: 1.0                     |
!| |_|       \_/   |_|_|_| |  Release: January, 2022           |
!|_ _ _ _ _ _ _ _ _ _ _ _ _|_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _|
! About: scalar allocator
! Modification: February, 2025

! ============================================================================
! CREATE MODELS
! ============================================================================

! ----------------------------------------------------------------------------
! physical unknown models
! ----------------------------------------------------------------------------

#undef __fvl_model__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_model__ t
#define __fvl_mesh__ mesh
#define __fvl_type__ fvl_field2d_cellmeansfield
#include "fvl_createscalarmodel2d_src.f"

#undef __fvl_model__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_model__ td
#define __fvl_mesh__ mesh
#define __fvl_type__ fvl_field2d_cellmeansfield
#include "fvl_createscalarmodel2d_src.f"

#undef __fvl_model__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_model__ tdd
#define __fvl_mesh__ mesh
#define __fvl_type__ fvl_field2d_cellmeansfield
#include "fvl_createscalarmodel2d_src.f"

#undef __fvl_model__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_model__ tddd
#define __fvl_mesh__ mesh
#define __fvl_type__ fvl_field2d_cellmeansfield
#include "fvl_createscalarmodel2d_src.f"

#undef __fvl_model__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_model__ p
#define __fvl_mesh__ mesh
#define __fvl_type__ fvl_field2d_cellmeansfield
#include "fvl_createscalarmodel2d_src.f"

#undef __fvl_model__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_model__ u
#define __fvl_mesh__ diamondmesh
#define __fvl_type__ fvl_field2d_cellmeansfield
#include "fvl_createvectormodel2d_src.f"

#undef __fvl_model__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_model__ xi
#define __fvl_mesh__ mesh
#define __fvl_type__ fvl_field2d_cellmeansfield
#include "fvl_createscalarmodel2d_src.f"

! ----------------------------------------------------------------------------
! physical property models
! ----------------------------------------------------------------------------

#undef __fvl_model__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_model__ k
#define __fvl_mesh__ mesh
#define __fvl_type__ fvl_field2d_edgequadratsfield
#include "fvl_createscalarmodel2d_src.f"

#undef __fvl_model__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_model__ cp
#define __fvl_mesh__ mesh
#define __fvl_type__ fvl_field2d_edgequadratsfield
#include "fvl_createscalarmodel2d_src.f"

#undef __fvl_model__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_model__ mu
#define __fvl_mesh__ diamondmesh
#define __fvl_type__ fvl_field2d_edgequadratsfield
#include "fvl_createscalarmodel2d_src.f"

#undef __fvl_model__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_model__ rho
#define __fvl_mesh__ mesh
#define __fvl_type__ fvl_field2d_edgequadratsfield
#include "fvl_createscalarmodel2d_src.f"

#undef __fvl_model__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_model__ rhod
#define __fvl_mesh__ diamondmesh
#define __fvl_type__ fvl_field2d_edgequadratsfield
#include "fvl_createscalarmodel2d_src.f"

#undef __fvl_model__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_model__ alpha
#define __fvl_mesh__ mesh
#define __fvl_type__ fvl_field2d_edgequadratsfield
#include "fvl_createscalarmodel2d_src.f"

#undef __fvl_model__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_model__ nu
#define __fvl_mesh__ diamondmesh
#define __fvl_type__ fvl_field2d_edgequadratsfield
#include "fvl_createscalarmodel2d_src.f"

#undef __fvl_model__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_model__ vp
#define __fvl_mesh__ diamondmesh
#define __fvl_type__ fvl_field2d_edgequadratsfield
#include "fvl_createscalarmodel2d_src.f"

#undef __fvl_model__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_model__ uc
#define __fvl_mesh__ mesh
#define __fvl_type__ fvl_field2d_edgequadratsfield
#include "fvl_createvectormodel2d_src.f"

#undef __fvl_model__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_model__ vc
#define __fvl_mesh__ diamondmesh
#define __fvl_type__ fvl_field2d_edgequadratsfield
#include "fvl_createvectormodel2d_src.f"

! ----------------------------------------------------------------------------
! physical source models
! ----------------------------------------------------------------------------

#undef __fvl_model__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_model__ tc
#define __fvl_mesh__ diamondmesh
#define __fvl_type__ fvl_field2d_cellquadratsfield
#include "fvl_createscalarmodel2d_src.f"

#undef __fvl_model__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_model__ ft
#define __fvl_mesh__ mesh
#define __fvl_type__ fvl_field2d_cellquadratsfield
#include "fvl_createscalarmodel2d_src.f"

#undef __fvl_model__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_model__ ftd
#define __fvl_mesh__ mesh
#define __fvl_type__ fvl_field2d_cellquadratsfield
#include "fvl_createscalarmodel2d_src.f"

#undef __fvl_model__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_model__ ftdd
#define __fvl_mesh__ mesh
#define __fvl_type__ fvl_field2d_cellquadratsfield
#include "fvl_createscalarmodel2d_src.f"

#undef __fvl_model__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_model__ fpu
#define __fvl_mesh__ diamondmesh
#define __fvl_type__ fvl_field2d_cellquadratsfield
#include "fvl_createvectormodel2d_src.f"

#undef __fvl_model__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_model__ gpu
#define __fvl_mesh__ mesh
#define __fvl_type__ fvl_field2d_cellquadratsfield
#include "fvl_createscalarmodel2d_src.f"

#undef __fvl_model__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_model__ g
#define __fvl_mesh__ diamondmesh
#define __fvl_type__ fvl_field2d_cellquadratsfield
#include "fvl_createvectormodel2d_src.f"

#undef __fvl_model__
#undef __fvl_mesh__
#undef __fvl_type__
#define __fvl_model__ beta
#define __fvl_mesh__ diamondmesh
#define __fvl_type__ fvl_field2d_cellquadratsfield
#include "fvl_createscalarmodel2d_src.f"

! end of file
