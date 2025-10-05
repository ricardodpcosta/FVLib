! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
!|  _ _ _   _   _   _      |                                   |
!| |_|_|_| |_| |_| |_|     |  The Great Finite Volume Library  |
!| |_|_    |_|_|_| |_|     |  Author: Ricardo Costa            |
!| |_|_|    \_\_/  |_|_ _  |  Version: 1.0                     |
!| |_|       \_/   |_|_|_| |  Release: January, 2022           |
!|_ _ _ _ _ _ _ _ _ _ _ _ _|_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _|
! About: create meshes
! Modification: February, 2025

! ============================================================================
! VARIABLES
! ============================================================================

type(fvl_ptr_runtime)::runtimes(2)
type(fvl_mesh3d),target::mesh
type(fvl_mesh3d),pointer::diamondmesh
type(fvl_patch3d),target::cellspatch,edgespatch,boundarypatch
type(fvl_patch3d),pointer::diamond_cellspatch,diamond_edgespatch,diamond_boundarypatch

! end of file
