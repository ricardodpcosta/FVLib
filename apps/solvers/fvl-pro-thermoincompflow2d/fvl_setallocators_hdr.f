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
! SET ALLOCATORS
! ============================================================================

! ----------------------------------------------------------------------------
! function lists
! ----------------------------------------------------------------------------

fvl_scalarfield2d_functionlist=fvl_scalarfunctionlist2d()
fvl_vectorfield2d_functionlist=fvl_vectorfunctionlist2d()

! ----------------------------------------------------------------------------
! physical unknown models
! ----------------------------------------------------------------------------

#undef __fvl_model__
#define __fvl_model__ t
#include "fvl_setscalarallocator2d_hdr.f"

#undef __fvl_model__
#define __fvl_model__ p
#include "fvl_setscalarallocator2d_hdr.f"

#undef __fvl_model__
#define __fvl_model__ u
#include "fvl_setvectorallocator2d_hdr.f"

! ----------------------------------------------------------------------------
! physical property models
! ----------------------------------------------------------------------------

#undef __fvl_model__
#define __fvl_model__ k
#include "fvl_setscalarallocator2d_hdr.f"

#undef __fvl_model__
#define __fvl_model__ cp
#include "fvl_setscalarallocator2d_hdr.f"

#undef __fvl_model__
#define __fvl_model__ mu
#include "fvl_setscalarallocator2d_hdr.f"

#undef __fvl_model__
#define __fvl_model__ rho
#include "fvl_setscalarallocator2d_hdr.f"

#undef __fvl_model__
#define __fvl_model__ alpha
#include "fvl_setscalarallocator2d_hdr.f"

#undef __fvl_model__
#define __fvl_model__ nu
#include "fvl_setscalarallocator2d_hdr.f"

#undef __fvl_model__
#define __fvl_model__ vp
#include "fvl_setscalarallocator2d_hdr.f"

! ----------------------------------------------------------------------------
! physical source models
! ----------------------------------------------------------------------------

#undef __fvl_model__
#define __fvl_model__ ft
#include "fvl_setscalarallocator2d_hdr.f"

#undef __fvl_model__
#define __fvl_model__ fpu
#include "fvl_setvectorallocator2d_hdr.f"

#undef __fvl_model__
#define __fvl_model__ gpu
#include "fvl_setscalarallocator2d_hdr.f"

#undef __fvl_model__
#define __fvl_model__ g
#include "fvl_setvectorallocator2d_hdr.f"

#undef __fvl_model__
#define __fvl_model__ beta
#include "fvl_setscalarallocator2d_hdr.f"

! ----------------------------------------------------------------------------
! specific fields
! ----------------------------------------------------------------------------

call fpu_fieldselector%setallocator("boussinesq",fvl_boussinesq_forcesource2d_allocate)

! end of file
